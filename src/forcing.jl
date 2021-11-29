##############
# Trait system
##############

abstract type ForcingTrait end
abstract type ForcingBehaviourTrait end
abstract type ForcingFactoryTrait end

# behaviour traits
# ------------------
abstract type DefaultForcing <: ForcingBehaviourTrait end

# factory traits
# ---------------
struct FromSerializedForcing <: ForcingFactoryTrait end
struct FromSimstratCalibration <: ForcingFactoryTrait
	ν::Float64
	η::Float64

	T0::Float64
	p0::Float64
	R::Float64
end

FromSimstratCalibration() = begin
		ν = 1.5e-6  # Kinematic Viscosity of Water
		η = 1/(2*pi)  # according to Lorke et al. #0.29 # calibration constant

		T0 = 273.15  # Standard Temperature [K]
		p0 = 1e5	# Standard Pressure [Pa]
		R = 8.314  # universal gas constant

		FromSimstratCalibration(
			ν,
			η,

			T0,
			p0,
			R
			)
end

#############
# Type System
#############
abstract type LakeForcing end

####################################
# MOBDiversityForcing implementation
####################################
struct MOBDiversityForcing{FactoryTrait, BehaviourTrait} <: LakeForcing
	water_temperature::Array{Float64,2}
	diffusivity::Array{Float64,2}
	dissipation::Array{Float64,2}
	pistonvelocity::Array{Float64,2}
end

# read a simstrat output file into a DataFrame
function read_simstrat_output_file(simstrat_output_path::AbstractString, variable_name::Symbol)
	#df = CSV.read(simstrat_output_path, delim=' ', ignorerepeated=true, header=true)
	#deletecols!(df, names(df)[end-1:end])
	df = CSV.read(simstrat_output_path, DataFrame; delim=',', header=true)
	df=stack(df, names(df)[2:end])
	names!(df, [:Depth, variable_name, :Timestamp])
	df[!,:Depth] = -[parse(Float64, "$d") for d in df[!,:Depth]]
	sort!(df, [:Timestamp, :Depth])
	return df
end

function simstrat_output_interpolant(df::DataFrame, tgrid::AbstractArray)
	df = df[(df[!,:Timestamp].>=minimum(tgrid)).&(df[!,:Timestamp].<=maximum(tgrid)), :]
	x = sort(unique(df[!,:Timestamp]))
	y = sort(unique(df[!,:Depth]))
	z = reshape(df[!,names(df)[2]], length(y), length(x))
	interpolate((y,x), z, Gridded(Linear()))
end

function sample_simstrat_output(df::DataFrame, tgrid::AbstractArray, zgrid::AbstractArray)
	reshape(simstrat_output_interpolant(df, tgrid)(zgrid[:],tgrid), length(zgrid), length(tgrid))'
end

(::Type{from_file{MOBDiversityForcing{FromSimstratCalibration, BehaviourTrait}}})(simstrat_calibration_path::AbstractString, tgrid::AbstractArray, grid::StaggeredGrid) where BehaviourTrait<:ForcingBehaviourTrait = begin
	# load rawdata
	rawdata_paths = [	("T_out.dat", :Temperature),
						("nuh_out.dat", :Diffusivity),
						("eps_out.dat", :Dissipation)]
	temperature_df, diffusivity_df, dissipation_df = [read_simstrat_output_file(joinpath(simstrat_calibration_path, fname), label) for (fname, label) in rawdata_paths]

	# interpolate to simulation grid
	water_temperature = sample_simstrat_output(temperature_df, tgrid, grid.z_centres)
	diffusivity = sample_simstrat_output(diffusivity_df, tgrid, grid.z_faces)
	diffusivity = diffusivity.*(3600.0*24.0)
	dissipation = sample_simstrat_output(dissipation_df, tgrid, [grid.z_faces[1]])

	factorytype = FromSimstratCalibration()
	pistonvelocity = k(factorytype, dissipation, water_temperature.+273.15)

	# return datastructure
	MOBDiversityForcing{FromSimstratCalibration, BehaviourTrait}(water_temperature, diffusivity, dissipation,pistonvelocity)
end

Sc_ch4_Wanninkhof(T) = 1909.4 - 120.78*(T-273.15) + 4.1555*((T-273.15)^2) - 0.080578*((T-273.15)^3) + 0.00065777*((T-273.15)^4)
k(factory::FromSimstratCalibration, dissipation::AbstractArray{Float64,2}, water_temperature::AbstractArray{Float64,2}) = factory.η.*(dissipation.*factory.ν).^0.25.*[Sc_ch4_Wanninkhof(T)^(-0.5) for T in water_temperature[:, 1]]

save(f::MOBDiversityForcing, path::AbstractString) = begin
    open(joinpath(path, "diffusivity.dat"), "w") do io
		writedlm(io, f.diffusivity, ',')
	end
	open(joinpath(path, "water_temperature.dat"), "w") do io
		writedlm(io, f.water_temperature, ',')
	end
	open(joinpath(path, "dissipation.dat"), "w") do io
		writedlm(io, f.dissipation, ',')
	end
	open(joinpath(path, "pistonvelocity.dat"), "w") do io
		writedlm(io, f.pistonvelocity, ',')
	end
end

function save_mmaped(f::MOBDiversityForcing, path::AbstractString)
	for (name, data) in [("water_temperature.bin", f.water_temperature), ("diffusivity.bin", f.diffusivity), ("dissipation.bin", f.dissipation), ("pistonvelocity.bin", f.pistonvelocity)]
		open(joinpath(path, name), "w+") do io
			write(io, size(data,1))
			write(io, size(data,2))
			write(io, data)
		end
	end
end

(::Type{from_file{MOBDiversityForcing{FromSerializedForcing, BehaviourTrait}}})(modeldefinition_base::AbstractString) where BehaviourTrait <: ForcingBehaviourTrait = begin
	ios=[]
	datasets=[]
	for name in ["water_temperature.bin", "diffusivity.bin", "dissipation.bin", "pistonvelocity.bin"]
		io = open(joinpath(modeldefinition_base, name))
		push!(ios, io)
		m = read(io, Int)
		n = read(io, Int)
		push!(datasets,Mmap.mmap(io, Array{Float64,2}, (m,n)))
	end
	(ios, MOBDiversityForcing{FromSerializedForcing, BehaviourTrait}(datasets[1], datasets[2], datasets[3], datasets[4]))
end


##################################################
struct MeteorologicalForcing
	air_temperature::Interpolation
	cloud_cover::Interpolation
	vapour_pressure::Interpolation
	global_radiation::Interpolation
	wind_speed::Interpolation
end

(::Type{from_file{MeteorologicalForcing}})(path::AbstractString) = begin
	# load simstrat forcing dataset
	# ==============================
	forcing_df = CSV.read(path, DataFrame; delim="\t")  # read csv
	forcing_df = disallowmissing!(forcing_df[completecases(forcing_df),:])  # remove missing values

	# interpolators for individual data series
	# =========================================
	timestamps = forcing_df[:t]
	air_temperature = LinearInterpolation(timestamps, forcing_df[Symbol("Tair (°C)")].+273.15)
	cloud_cover = LinearInterpolation(timestamps, forcing_df[Symbol("cloud coverage")])
	vapour_pressure = LinearInterpolation(timestamps, forcing_df[Symbol("vap (mbar)")])

	# combine wind speed components
	u10 = forcing_df[Symbol("u (m/s)")]
	v10 = forcing_df[Symbol("v (m/s)")]
	wind_speed = LinearInterpolation(timestamps, sqrt.(u10.^2+v10.^2))

	# global shortwave irradiation
	G = forcing_df[Symbol("Fsol (W/m2)")]
	C = forcing_df[Symbol("cloud coverage")]

	## Wuest et al.
	#Fdir = (1.-C)./((1.-C)+0.5*C)
	#Fdiff = (0.5*C)./((1.-C)+0.5*C)
	#Alb_dir = 0.2
	#Alb_diff = 0.066
	#Hs = G.*Fdir*(1-Alb_dir)+G.*Fdiff*(1-Alb_diff)

	## Simstrat
	Hs = (1-0.08)*G
	global_radiation = LinearInterpolation(timestamps, Hs)

	# return data object
	MeteorologicalForcing(air_temperature, cloud_cover, vapour_pressure, global_radiation, wind_speed)
end