# macro to define lake physics functions in short notation
macro RDSfn(fnexpr)
	fnrepr = repr(fnexpr)
	fnsignature = split(fnrepr[3:end-1], "=")[1]
	# first argument of function signature
	fnfirstarg = split(split(match(r"\((.*?)\)", fnsignature)[1], ",")[1],"::")
	fn_classname = fnfirstarg[1]
	if length(fnfirstarg) > 1
		# look for field names and prefix them with the class
		fn_class_fieldnames = fieldnames(eval(Meta.parse(fnfirstarg[2])))
		for field in fn_class_fieldnames
			fnrepr = replace(fnrepr, Regex(string("(?<![\\w_])(|\\d+)(", field, ")(?![\\w_])"))=>SubstitutionString(string("\\g<1>",fn_classname, ".", field)))
		end
	else
		error(string("You have to define the type of ", fn_classname ," in the function signature"))
	end
	eval(Meta.parse(fnrepr))
end

struct MOBDiversityLakeRDS <: ReactionDiffusionSystem
	forcing::MOBDiversityForcing
	growthmodels::Array{MOBGrowthModel,1}
	sediment_exchange_k::Float64
	sediment_diffusivity::Float64
	sediment_vmax::Float64
	OM_sedimentation_rate::Float64
	sediment_kI_o2::Float64
	o2_consumption_rate::Float64
	background_mob_rate::Float64
	baseline_mortality::Float64
	T0::Float64
	p0::Float64
	R::Float64
end
MOBDiversityLakeRDS(forcing::MOBDiversityForcing, growthmodels::Array{<:MOBGrowthModel,1}, sediment_exchange_k::Float64, sediment_diffusivity::Float64, sediment_vmax::Float64, OM_sedimentation_rate::Float64, sediment_kI_o2::Float64, o2_consumption_rate::Float64, background_mob_rate::Float64, baseline_mortality::Float64) = MOBDiversityLakeRDS(forcing, growthmodels, sediment_exchange_k, sediment_diffusivity, sediment_vmax, OM_sedimentation_rate, sediment_kI_o2, o2_consumption_rate, background_mob_rate, baseline_mortality, 273.15, 1e5, 8.314)

FiniteVolumeRDS.update!(rds::MOBDiversityLakeRDS, u::AbstractArray{Float64,2}, t::Float64, tidx::Int64, trep::Float64, tidx_rep::Int64, diffusivities::Array{Float64,2}, fluxes::Array{Float64,2}, sources::Array{Float64,2}) = begin
	# state vars: 1 ch4, 2 o2, 3 sed_ch4, 4 sed_o2, 5 sed_org, mob...
	nvars = size(diffusivities)[2]
	# diffusion
	diffusivities[:] .= repeat(rds.forcing.diffusivity[tidx_rep,:], outer=[1,nvars])[:]
	diffusivities[:,3] .= rds.sediment_diffusivity
	diffusivities[:,4] .= rds.sediment_diffusivity
	diffusivities[:,5] .= rds.sediment_diffusivity

	# reset sources
	sources[:,:] .= 0.0

	# exchange with atmosphere
	######################################################
	fluxes[1,1] = rds.forcing.pistonvelocity[tidx_rep]*3600.0*24.0*(Ceq_ch4(rds,rds.forcing.water_temperature[tidx_rep, 1]+273.15, 0.0)*1.0e6-u[1,1])
	fluxes[1,2] = rds.forcing.pistonvelocity[tidx_rep]*3600.0*24.0*(Ceq_o2(rds,rds.forcing.water_temperature[tidx_rep, 1]+273.15, 0.0)*1.0e6-u[1,2])

	# sedimentation of organic matter
	###################################
	fluxes[1,5] = rds.OM_sedimentation_rate*3600.0*24.0

	# sediment methane production and exchange with water column
	######################################################
	# methane
	sediment_exchange_k = rds.sediment_exchange_k
	fluxes[end,1] = sediment_exchange_k*3600.0*24.0*(u[end,1]-u[1,3])  # into water column
	fluxes[1,3] = sediment_exchange_k*3600.0*24.0*(u[end,1]-u[1,3])  # out of sediment
	# oxygen
	fluxes[end,2] = sediment_exchange_k*3600.0*24.0*(u[end,2]-u[1,4])  # into water column
	fluxes[1,4] = sediment_exchange_k*3600.0*24.0*(u[end,2]-u[1,4])  # out of sediment
	# sediment methane production (assume monod type inhibition by oxygen)
	sources[:,3] .= (rds.sediment_vmax*3600.0*24.0).*u[:,5].*(rds.sediment_kI_o2./(rds.sediment_kI_o2.+u[:,4]))
	sources[:,5] .= .-sources[:,3]


	# water column methane oxidation and oxygen consuption
	######################################################

	# photosynthesis
	#chla = 50000.0.*rds.forcing.water_temperature[tidx_rep, :]./(rds.forcing.water_temperature[tidx_rep, :].+273.15.+10.0)

	# background oxygen consuption
	#bod = 1.13.^(rds.forcing.water_temperature[tidx_rep, :].-273.15.-20).*rds.o2_consumption_rate
	#sources[:,2] .= -rds.o2_consumption_rate.*30000#u[:,2]
	#sources[:,4] .= -rds.o2_consumption_rate.*u[:,5]
	sources[:,2] .= -rds.o2_consumption_rate#u[:,2]
	#sources[:,4] .= -rds.o2_consumption_rate.*64.0

	# growth
	rates = [mox(growthmodel, u[:,1], u[:,2], rds.forcing.water_temperature[tidx_rep, :]) for growthmodel in rds.growthmodels]
	for (i, rate) in enumerate(rates)
		#growth
		sources[:,5+i] .= rate[1].*u[:,5+i].+rds.background_mob_rate
		# ch4 consuption
		sources[:,1] .= sources[:,1] .+ rate[2].*u[:,5+i]
		# o2 consumption
		sources[:,2] .= sources[:,2] .+ rate[3].*u[:,5+i]
	end
	nothing
end

@RDSfn Ceq_ch4(rds::MOBDiversityLakeRDS,T,S) = exp(-68.8862 + 101.4956*(100/T) + 28.7314*log(T/100) + S*(-0.076146 + 0.043970*(T/100) - 0.006872*(T/100)^2))*(p0/(R*T0))*1.8e-6
@RDSfn Ceq_o2(rds::MOBDiversityLakeRDS,T,S) = exp(-58.3877 + 85.8079*(100/T) + 23.8439*log(T/100) + S*(-0.034892 + 0.015568*(T/100) - 0.0019387*(T/100)^2))*(p0/(R*T0))*0.2095


#
## Air/Water transfer velocity
##############################
#
## Schmitt-Number
#Sc_ch4_Wanninkhof(T) = 1909.4 - 120.78*(T-273.15) + 4.1555*((T-273.15)^2) - 0.080578*((T-273.15)^3) + 0.00065777*((T-273.15)^4)
#@RDSfn Ceq_ch4(p::MOBDiversityLakeRDS{<:DefaultRDS},T,S) = exp(-68.8862 + 101.4956*(100/T) + 28.7314*log(T/100) + S*(-0.076146 + 0.043970*(T/100) - 0.006872*(T/100)^2))*(p0/(R*T0))*1.8e-6
## Gas Transfer Velocity
#
#@RDSfn k(p::MOBDiversityLakeRDS{<:DefaultRDS}, u,t) = η*(ϵ(p,u,t)*ν)^0.25*Sc_ch4_Wanninkhof(u[5])^(-0.5)
#@RDSfn k_u(p::MOBDiversityLakeRDS{<:DefaultRDS}, u,t) = η*(ϵ_u(p,u,t)*ν)^0.25*Sc_ch4_Wanninkhof(u[5])^(-0.5)
#@RDSfn k_B(p::MOBDiversityLakeRDS{<:DefaultRDS}, u,t) = η*(ϵ_B(p,u,t)*ν)^0.25*Sc_ch4_Wanninkhof(u[5])^(-0.5)
#
#@RDSfn Fatm(p::MOBDiversityLakeRDS{<:DefaultRDS}, u,t) = -k(p,u,t)*(u[3]-Ceq_ch4(p,u[5],0)*1e6)
#

#@enum STATEVAR CH4=1 O2=2 MOB1=3
#Base.to_index(s::STATEVAR) = Int(s)
#
#function mobdiversity_ode(p::MOBDiversityLakeRDS{<:DefaultRDS}, u,t,tidx,sources,fluxes,diffusivities)
#	# ch4, o2, MOB1, MOB2, MOB3, ...
#	# u[:,CH4]
#	nothing
#end