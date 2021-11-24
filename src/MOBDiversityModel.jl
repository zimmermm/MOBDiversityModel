module MOBDiversityModel
using Reexport

@reexport using FiniteVolumeRDS
using Dates
using Interpolations
using CSV, DelimitedFiles, DataFrames, JSON
using TextParse
using TextParse: Record, Field, Numeric, tryparsenext
using Mmap
using LambertW

export	Interpolation, LinearInterpolation,
		top, bottom
export	datetime_to_days, days_to_datetime
export	from_file
export	MOBDiversityLakeRDS, update!
export	MOBDiversityForcing, FromSerializedForcing, FromSimstratCalibration, DefaultForcing

struct from_file{T} end

include("interpolations.jl")
include("dateconversion.jl")
include("profileshapes.jl")
include("forcing.jl")
include("growthmodel.jl")
include("mobdiversity_rds.jl")

function read_modeldefinition(path::String)
	modeldefinition = Dict()
	open(path, "r") do jsonfile
		jsontxt = read(jsonfile, String)  # file information to string
		modeldefinition=JSON.parse(jsontxt)  # parse and transform data
	end
	return modeldefinition
end

function decode_growthmodel(trait::Array{Float64,1}, baseline_mortality::Float64)
	#eoc
	Tmin = 7.0
	Tmax = 33.0
	Topt = 15.0
	if trait[1] >= 0.5
		# warm adapted
		trait_scale = (trait[1]-0.5)*2
		Tmin = Tmin+5.0*trait_scale
		Tmax = Tmax+22.0*trait_scale
		Topt = Topt+10.0*trait_scale
	else
		# cold adapted
		trait_scale = 1.0-trait[1]*2
		Tmax = Tmax-22.0*trait_scale
		Tmin = Tmin-6.0*trait_scale
		Topt = Topt-10.0*trait_scale
	end
	# translate kinetics
	Km_CH4 = 10.0^(2.5-3.5*trait[2])*1000.0
	Vmax = (0.1+1.9*trait[3])*1e-5
	
	# narrow
	#Tmin = 7.0
	#Tmax = 50.0
	#Topt = 20.0
	#if trait[1] >= 0.5
	#	# warm adapted
	#	trait_scale = (trait[1]-0.5)*2
	#	Tmin = Tmin+5.0*trait_scale
	#	Tmax = Tmax+20.0*trait_scale
	#	Topt = Topt+10.0*trait_scale
	#else
	#	# cold adapted
	#	trait_scale = 1.0-trait[1]*2
	#	Tmax = Tmax-20.0*trait_scale
	#	Tmin = Tmin-5.0*trait_scale
	#	Topt = Topt-15.0*trait_scale
	#end
	# translate kinetics
	#Km_CH4 = 10.0^(2.5-3.5*trait[2])*1000.0
	#Vmax = (0.1+1.9*trait[3])*1e-5
	## investments[2:end] have to sum to 1
	#Tmin = 7.0
	#Tmax = 50.0
	#Topt = 20.0
	#if investments[1] > 0.5
	#	Tmin = Tmin+5.0*investments[2]
	#	Tmax = Tmax+20.0*investments[2]
	#	Topt = Topt+10.0*investments[2]
	#else
	#	Tmax = Tmax-20.0*investments[2]
	#	Tmin = Tmin-5.0*investments[2]
	#	Topt = Topt-15.0*investments[2]
	#end
	#Km_CH4 = 10.0^(1.5-2.5*investments[3])*1000.0
	#Vmax = (1.2+2.0*investments[4])*1.0e-5
	println("Tmin: ", Tmin, ", Topt: ", Topt, ", Tmax: ", Tmax, ", Km_CH4: ", Km_CH4, "Vmax: ", Vmax)
	return MOBMonodRatkowskiModel(Tmin, Tmax, Topt, Vmax, Km_CH4, 300.0, 0.3, baseline_mortality, :infer)
end

(::Type{from_file{:MOBDiversityModel}})(modeldefinition_path::AbstractString; verbose=true, use_precomputed_forcing=false) = begin
	modeldefinition = read_modeldefinition(modeldefinition_path)
	# simulation
	##############
	simulation_start = datetime_to_days(Dates.DateTime(modeldefinition["Simulation"]["Start"], Dates.DateFormat("d.m.y")))
	simulation_end = datetime_to_days(Dates.DateTime(modeldefinition["Simulation"]["End"], Dates.DateFormat("d.m.y")))
	dt = modeldefinition["Simulation"]["dt"]
	t = collect(simulation_start:dt:simulation_end)
	# morphology
	##############
	morphology = modeldefinition["Morphology"]
	bathymetry = morphology["Bathymetry"]
	nz = morphology["Gridcells"]
	grid = StaggeredGrid(	0.0,
							morphology["Maxdepth"],
							nz,
							convert(Array{Float64,1}, bathymetry["zareas"]),
							convert(Array{Float64,1}, bathymetry["areas"]))
	# forcing
	##########
	modeldefinition_base = Base.Filesystem.dirname(modeldefinition_path)
	if !use_precomputed_forcing
		forcing = from_file{MOBDiversityForcing{FromSimstratCalibration, DefaultForcing}}(modeldefinition["Forcing"]["Simstrat calibration"], t, grid)
		#save(forcing, modeldefinition_base)
		save_mmaped(forcing, modeldefinition_base)
		forcing = nothing
	end

	# Parameter
	###########
	sediment_exchange_k = modeldefinition["Parameter"]["sediment_exchange_k"]
	sediment_diffusivity = modeldefinition["Parameter"]["sediment_diffusivity"]
	sediment_vmax = modeldefinition["Parameter"]["sediment_vmax"]
	OM_sedimentation_rate = modeldefinition["Parameter"]["OM_sedimentation_rate"]
	sediment_kI_o2 = modeldefinition["Parameter"]["sediment_kI_o2"]
	o2_consumption_rate = modeldefinition["Parameter"]["o2_consumption_rate"]
	background_mob_rate = modeldefinition["Parameter"]["background_mob_rate"]
	baseline_mortality = modeldefinition["Parameter"]["baseline_mortality"]

	growthmodel_definitions = modeldefinition["Growthmodel"]
	growthmodels::Array{MOBGrowthModel,1} = repeat([NoGrowth()], outer=[length(growthmodel_definitions)])
	for (i, investments) in enumerate(growthmodel_definitions)
		growthmodels[i] = decode_growthmodel(convert(Array{Float64,1}, investments), baseline_mortality)
	end

	(ios, forcing) = from_file{MOBDiversityForcing{FromSerializedForcing, DefaultForcing}}(modeldefinition_base)
	#mob1 = MOBMonodModel(1.8e-05, 18000.0, 1000.0, 0.13)
	#mob2 = MOBMonodModel(0.5e-05, 500.0, 1000.0, 0.13)
	#thottathil_mox = ThottathilModel()
	#mob1 = MOBMonodRatkowskiModel(0.0, 50.0, 20.0, 1.8e-5, 500.0, 300.0, 0.13, :infer)
	#mob2 = MOBMonodRatkowskiModel(0.0, 50.0, 18.0, 1.5e-5, 700.0, 300.0, 0.13, :infer)
	#mob3 = MOBMonodRatkowskiModel(0.0, 50.0, 21.0, 1.8e-5, 600.0, 300.0, 0.13, :infer)
	#mob4 = MOBMonodRatkowskiModel(0.0, 50.0, 7.0, 0.9e-5, 17000.0, 300.0, 0.13, :infer)
	#mob5 = MOBMonodRatkowskiModel(0.0, 50.0, 10.0, 0.8e-5, 15000.0, 300.0, 0.13, :infer)
	#mob6 = MOBMonodRatkowskiModel(0.0, 50.0, 5.0, 0.8e-5, 20000.0, 300.0, 0.13, :infer)
	#growthmodels::Array{MOBGrowthModel,1} = [mob1, mob2, mob3, mob4, mob5, mob6]
	rds = MOBDiversityLakeRDS(forcing, growthmodels, sediment_exchange_k, sediment_diffusivity, sediment_vmax, OM_sedimentation_rate, sediment_kI_o2, o2_consumption_rate, background_mob_rate, baseline_mortality)
	nvars = 5+length(rds.growthmodels) #11^3
	fvd = FiniteVolumeDiscretization(grid,nvars,modeldefinition["Simulation"]["dt"])
	u0 = repeat([1.0e-9], outer=[nz,nvars])
	u0[:,1].=1000
	u0[:,2].=320000
	u0[:,6:end] .= 3500
	return (fvd, u0, rds, simulation_start, simulation_end, ios, modeldefinition["Simulation"]["repetitions"], modeldefinition["Simulation"]["thinning"])
end

end # module
