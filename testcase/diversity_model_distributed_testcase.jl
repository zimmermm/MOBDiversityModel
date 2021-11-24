using Distributed
addprocs(6)

@everywhere using FiniteVolumeRDS
@everywhere using MOBDiversityModel
@everywhere using Serialization
@everywhere using DataFrames, CSV, Printf, DelimitedFiles, JSON
@everywhere using Interpolations
@everywhere using Mmap
@everywhere using LinearAlgebra
@everywhere using SharedArrays
@everywhere using Combinatorics
@everywhere using Dates
@everywhere using PyPlot, Colors

@everywhere const results_base =".\\results\\"
@everywhere const model_workspace = ".\\modeldefinition\\diversity\\"
@everywhere const modeldefinition_path = joinpath(model_workspace, "modeldefinition.json")

########################################
# Basic Functionality
########################################

function load_modeldefinition()
	from_file{:MOBDiversityModel}(modeldefinition_path, use_precomputed_forcing=false)
end

function load_modeldefinition(year, seasonal, tn, background_mob_rate, baseline_mortality)
	update_modeldefinition(year, seasonal, tn, background_mob_rate, baseline_mortality)
	from_file{:MOBDiversityModel}(modeldefinition_path, use_precomputed_forcing=false)
end

function run()
	run(load_modeldefinition())
end

function run(modeldefinition)
 	(fvd, u0, rds, simulation_start, simulation_end, ios, repetitions, thinning) = modeldefinition
 	sol = solve(fvd, u0, simulation_start, simulation_end, rds; repetitions=repetitions, solutiontype=ThinnedFiniteVolumeRDSSolution(thinning, joinpath(results_base, "solution.dat")))
 	FiniteVolumeRDS.save(sol)
 	finalize(sol.sol)
 	sol = nothing
 	GC.gc()
 
 	Serialization.serialize(joinpath(model_workspace, "discretization.jl"), fvd)
 	for io in ios
 		close(io)
 	end
end

function plot_solution()
	modeldefinition_base = Base.Filesystem.dirname(modeldefinition_path)
	sol = ThinnedFiniteVolumeRDSSolution(joinpath(results_base, "solution.dat"))
	fvd = Serialization.deserialize(joinpath(model_workspace, "discretization.jl"))
	(ios, forcing) = from_file{MOBDiversityForcing{FromSerializedForcing, DefaultForcing}}(modeldefinition_base)
	npar = size(sol)[2]
	nmox = npar-5
	modeldefinition = MOBDiversityModel.read_modeldefinition(joinpath(model_workspace, "modeldefinition.json"))
	simulation_start = datetime_to_days(DateTime(modeldefinition["Simulation"]["Start"], DateFormat("d.m.y")))
	simulation_end = datetime_to_days(DateTime(modeldefinition["Simulation"]["End"], DateFormat("d.m.y")))
	dt = modeldefinition["Simulation"]["dt"]
	repetitions = modeldefinition["Simulation"]["repetitions"]

	fig, ax = plt.subplots(1,1,figsize=(10,10))
	mob_sums = zeros(length(sol), nmox)
	for t in 1:length(sol)
		for mob_i in 1:nmox
			mob_sums[t, mob_i] = sum(sol[t,mob_i])
		end
	end
	plot_sol = log10.(mob_sums[100:end,:]'./1000)
	plot_sol[plot_sol .< -3] .= -3
	println(maximum(plot_sol))
	println(minimum(plot_sol))
	p = ax.contourf(plot_sol)
	fig.colorbar(p, ax=ax)
end

function plot_relative_abundance()
	modeldefinition_base = Base.Filesystem.dirname(modeldefinition_path)
	sol = ThinnedFiniteVolumeRDSSolution(joinpath(results_base, "solution.dat"))
	fvd = Serialization.deserialize(joinpath(model_workspace, "discretization.jl"))
	(ios, forcing) = from_file{MOBDiversityForcing{FromSerializedForcing, DefaultForcing}}(modeldefinition_base)
	npar = size(sol)[2]
	nmox = npar-5

	fig, ax = plt.subplots(nmox,1,figsize=(10,10))
	burnin = 2*365#convert(Int64, trunc(length(sol)/5))
	print(size(sol))
	total_mob = zeros(size(sol)[1],size(sol)[3]-burnin+1)
	for i in 1:nmox
		total_mob .= total_mob .+ sol[burnin:end,5+i].*12e-6./0.42e-12./1.0e6
	end
	for i in 1:nmox
		p = ax[i].contourf(sol[burnin:end,5+i].*12e-6./0.42e-12./1.0e6./total_mob)
		year1 = sum(sol[end-3*365:end-2*365,i])
		year2 = sum(sol[end-1*365:end,i])
		println(year1)
		println(year2)
		println(abs(year2-year1)/(year1+year2))
		ax[i].invert_yaxis()
		fig.colorbar(p, ax=ax[i])
	end

	savefig(joinpath(results_base, "relative_abundance.pdf"), format="pdf", dpi=600)
	finalize(sol.sol)
end


########################################
# Performance Score
########################################

@everywhere function performance!(refs,samples, sample_nr, modeldefinition,year,narrow)
	println("Sample ", sample_nr, " of ", size(samples)[1], " (", sample_nr/size(samples)[1],"): ", samples[sample_nr,1:3])
	(fvd, u0, rds, simulation_start, simulation_end, ios, repetitions, thinning) = modeldefinition
	# translate trait sample to growthmodel
	if narrow
		rds.growthmodels[:] .= vcat(refs,growthmodel_narrow(samples[sample_nr, 1:3], rds.baseline_mortality))
	else
		rds.growthmodels[:] .= vcat(refs,growthmodel(samples[sample_nr, 1:3], rds.baseline_mortality))
	end
		
	# run model
	sol = solve(fvd, u0, simulation_start, simulation_end, rds; repetitions=repetitions, solutiontype=ThinnedMemoryBoundFiniteVolumeRDSSolution(thinning))

	for growthmodel_idx in eachindex(rds.growthmodels)
		start_idx = 4+(growthmodel_idx-1)*repetitions
		end_idx = start_idx+repetitions-1
		# evaluate solution
		if year != 2016
			samples[sample_nr, start_idx:end_idx] .= [sum(sol[end-rep*364:end-(rep-1)*364,5+growthmodel_idx]) for rep in 1:repetitions]
		else
			samples[sample_nr, start_idx:end_idx] .= [sum(sol[end-rep*365:end-(rep-1)*365,5+growthmodel_idx]) for rep in 1:repetitions]
		end
	end

	sol = nothing
	return nothing
end

# translate trait to growthmodel
@everywhere function growthmodel(trait::Union{AbstractArray,SubArray},baseline_mortality::Float64)
	# normalized trait space consists of 3 dimensions [T, Km, Vmax] from 0.0 to 1.0
	# translate temperature range
	println("# Initialize growthmodel")
	Tmin = 7.0
	Tmax = 50.0
	Topt = 20.0
	if trait[1] >= 0.5
		# warm adapted
		trait_scale = (trait[1]-0.5)*2
		Tmin = Tmin+5.0*trait_scale
		Tmax = Tmax+20.0*trait_scale
		Topt = Topt+10.0*trait_scale
	else
		# cold adapted
		trait_scale = 1.0-trait[1]*2
		Tmax = Tmax-20.0*trait_scale
		Tmin = Tmin-5.0*trait_scale
		Topt = Topt-15.0*trait_scale
	end
	# translate kinetics
	Km_CH4 = 10.0^(2.5-3.5*trait[2])*1000.0
	Vmax = (0.1+1.9*trait[3])*1e-5
	println("Tmin: ", Tmin, " Topt: ", Topt, " Tmax: ", Tmax, " Km: ", Km_CH4, " Vmax: ", Vmax)
	return MOBDiversityModel.MOBMonodRatkowskiModel(Tmin, Tmax, Topt, Vmax, Km_CH4, 300.0, 0.3, baseline_mortality, :infer)
end

@everywhere function growthmodel_narrow(trait::Union{AbstractArray,SubArray},baseline_mortality::Float64)
	# normalized trait space consists of 3 dimensions [T, Km, Vmax] from 0.0 to 1.0
	# translate temperature range
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
	
	println("Tmin: ", Tmin, " Topt: ", Topt, " Tmax: ", Tmax, " Km: ", Km_CH4, " Vmax: ", Vmax)
	return MOBDiversityModel.MOBMonodRatkowskiModel(Tmin, Tmax, Topt, Vmax, Km_CH4, 300.0, 0.3, baseline_mortality, :infer)
end

function update_modeldefinition(year, seasonal, tn, background_mob_rate, baseline_mortality)
	modeldefinition = MOBDiversityModel.read_modeldefinition(modeldefinition_path)
	modeldefinition["Simulation"]["Start"] = Dates.format(Date(year, 1, 1), "dd.mm.yyyy")
	modeldefinition["Simulation"]["End"] = Dates.format(Date(year, 12, 31), "dd.mm.yyyy")
	
	if seasonal
		modeldefinition["Forcing"]["Simstrat calibration"] = "D:/Matthias/results/physicalmodel/seasonal/3_methane_iteration_revised/Rotsee_Output/"
		modeldefinition["Forcing"]["Weather"] = "D:/Matthias/results/physicalmodel/seasonal/3_methane_iteration_revised/Rotsee_Input/rotsee_forcing.dat"
	else
		modeldefinition["Forcing"]["Simstrat calibration"] = "D:/Matthias/results/diversity/physicalmodel/Rotsee_Output/"
		modeldefinition["Forcing"]["Weather"] = "D:/Matthias/results/diversity/physicalmodel/Rotsee_Input/rotsee_forcing.dat"
	end
	
	modeldefinition["Growthmodel"] = zeros(3, tn+1)
	
	modeldefinition["Parameter"]["background_mob_rate"] = background_mob_rate
	modeldefinition["Parameter"]["baseline_mortality"] = baseline_mortality
	
	open(modeldefinition_path, "w") do jsonfile
		JSON.print(jsonfile, modeldefinition,4)
	end
end



function map_Tn(n,highres,ref, year,seasonal,narrow,background_mob_rate,baseline_mortality,name)
	n_samples=11^3
	dtrait = 0.1
	nparallel = 11
	if highres
		n_samples=21^3
		dtrait = 0.05
		nparallel = 21
	end	
	
	modeldefinition = load_modeldefinition(year, seasonal, n,background_mob_rate,baseline_mortality)
	repetitions = modeldefinition[7]
	# pre allocate sample array
	samples::SharedArray{Float64,2} = zeros(n_samples, 3+(n+1)*repetitions)
	horizon::Int64 = 0
	# initialize start samples
	for x in 0.0:dtrait:1.0
		for y in 0.0:dtrait:1.0
			@sync @distributed for i in 1:nparallel
				samples[horizon+i, 1:3] .= [x,y,0.0+(i-1)*dtrait]
				performance!(ref,samples,horizon+i,modeldefinition,year,narrow)
			end
			horizon = horizon+nparallel
		end
	end

	disk_samples::Array{Float64,2} = zeros(n_samples, 3+(n+1)*repetitions)
	disk_samples[:] .= samples[:]
	Serialization.serialize(joinpath(model_workspace, string("T", n, "_grid_", n_samples, "_", year, name, ".jl")), disk_samples)
end


function map_traitspace_manuscript()
	baseline_mortality = 0.0024
	background_mob_rate = 0.0
	
	eco_traits_community1 = [
		[0.9, 0.9, 0.8], 
		[0.3, 0.5, 0.54], 
		[0.05, 0.1, 0.34],
		[0.47, 0.99, 0.36]
	]
	eco_traits_community2 = [
		[0.1,0.9,0.1],
        [0.35,0.9,0.9],
		[0.2,0.35,0.65]
	]
	eco_traits_community3 = [
		[0.12, 0.75, 0.85],
		[0.47, 0.99, 0.36],
		[0.75, 0.7, 0.73]
	]
	# check T4 and relaxation
	baseline_mortality = 0.0024
	background_mob_rate = 0.0
	
	#Community1-T3
	ref = [growthmodel_narrow(eco_traits_community1[i], baseline_mortality) for i in 1:3]
	map_Tn(3, true,ref, 2005,false,true,background_mob_rate,baseline_mortality,"_longterm_eco_community1")
	
	#Community2-T3
	ref = [growthmodel_narrow(eco_traits_community2[i], baseline_mortality) for i in 1:3]
	map_Tn(3, true,ref, 2005,false,true,background_mob_rate,baseline_mortality,"_longterm_eco_community2")
	
	#Community3-T2
	#ref = [growthmodel_narrow(eco_traits_community3[i], baseline_mortality) for i in 1:2]
	#map_Tn(2,false,ref, 2016,true,true,background_mob_rate,baseline_mortality,"_seasonal_eco_community3_trait12")
	
	#Community1-T3
	#ref = [growthmodel_narrow(eco_traits_community1[i], baseline_mortality) for i in 1:3]
	#map_Tn(3, true,ref, 2016,true,true,background_mob_rate,baseline_mortality,"_seasonal_eco_community1")
	
	#Community2-T3
	#ref = [growthmodel_narrow(eco_traits_community2[i], baseline_mortality) for i in 1:3]
	#map_Tn(3, true,ref, 2016,true,true,background_mob_rate,baseline_mortality,"_seasonal_eco_community2")
	
	#Community3-T3
	#ref = [growthmodel_narrow(eco_traits_community3[i], baseline_mortality) for i in 1:3]
	#map_Tn(3, true,ref, 2016,true,true,background_mob_rate,baseline_mortality,"_seasonal_eco_community3")
	
	#community 1
	
	#Community1-T1
	#map_Tn(1,false,[growthmodel_narrow(eco_traits_community1[1], baseline_mortality)], 2016,true,true,background_mob_rate,baseline_mortality,"_seasonal_eco_community1_trait1")
	#map_Tn(1,true,[growthmodel_narrow(eco_traits_community1[2], baseline_mortality)], 2016,true,true,background_mob_rate,baseline_mortality,"_seasonal_eco_community1_trait2")
	#map_Tn(1,false,[growthmodel_narrow(eco_traits_community1[3], baseline_mortality)], 2016,true,true,background_mob_rate,baseline_mortality,"_seasonal_eco_community1_trait3")
	#map_Tn(1,false,[growthmodel_narrow(eco_traits_community1[4], baseline_mortality)], 2016,true,true,background_mob_rate,baseline_mortality,"_seasonal_eco_community1_trait4")
	
	#T2
	#ref = [growthmodel_narrow(eco_traits_community1[i], baseline_mortality) for i in 1:2]
	#map_Tn(2,true,ref, 2016,true,true,background_mob_rate,baseline_mortality,"_seasonal_eco_community_2_trait12")
	
	#T4
	#ref = [growthmodel_narrow(eco_traits_community1[i], baseline_mortality) for i in 1:4]
	#map_Tn(4, false,ref, 2016,true,true,background_mob_rate,baseline_mortality,"_seasonal_eco")
	
	#community 2
	#T2
	#ref = [growthmodel_narrow(eco_traits_community2[i], baseline_mortality) for i in 1:2]
	#map_Tn(2,false,ref, 2016,true,true,background_mob_rate,baseline_mortality,"_seasonal_eco_trait12")
	
	
	#community 3
	#map_Tn(1,false,[growthmodel_narrow(eco_traits_community3[1], baseline_mortality)], 2016,true,true,background_mob_rate,baseline_mortality,"_seasonal_eco_community3_trait1")
	
	#T0
	#map_Tn(0,true,[], 2016,true,true,background_mob_rate,baseline_mortality,"_seasonal_eco")
end





