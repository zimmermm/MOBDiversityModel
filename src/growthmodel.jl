#===================================
Growth Model Implementations
===================================#

# abstrat type
abstract type MOBGrowthModel end


# No growth
# =============================
struct NoGrowth <: MOBGrowthModel end

mox(model::NoGrowth, ch4::Array{Float64,1}, o2::Array{Float64,1}, temperature::Array{Float64,1}) = zeros(length(ch4))
μ(model::NoGrowth, ch4::Array{Float64,1}, o2::Array{Float64,1}, temperature::Array{Float64,1}) = zeros(length(ch4))


# Monod kinetic implementation
# =============================
struct MOBMonodModel <: MOBGrowthModel
	Vmax::Float64
	Km_ch4::Float64
	Km_o2::Float64
	y::Float64
end

mortality(model::MOBMonodModel, ch4::Array{Float64,1}, o2::Array{Float64,1}, temperature::Array{Float64,1}) = 0.0024.+0.022.*o2./(model.Km_o2.+o2)

mox(model::MOBMonodModel, ch4::Array{Float64,1}, o2::Array{Float64,1}, temperature::Array{Float64,1}) = begin
	mox_rates=(model.Vmax*3600.0*24.0).*(ch4./(model.Km_ch4.+ch4)).*(o2./(model.Km_o2.+o2))
	mox_rates[mox_rates.<0.0] .= 0.0
	μ = mox_rates.*model.y.-mortality(model, ch4, o2, temperature)
	o2_rates = .-(2.0-model.y).*mox_rates
	ch4_rates = .-mox_rates
	return (μ, ch4_rates, o2_rates)
end


# Monod kinetic implementation
# =============================
struct MOBMonodRatkowskiModel <: MOBGrowthModel
	# ratkowski
	Tmin::Float64
	Tmax::Float64
	b::Float64
	d::Float64
	# monod
	Km_ch4::Float64
	Km_o2::Float64
	KI_o2::Float64
	y::Float64
	baseline_mortality::Float64
	# mortality
end

@inline ratkowski_d(Tmin, Tmax, Topt, branch) = (
									Tmin*lambertw(-((exp(Tmax/(Tmin-Topt)-Topt/(Tmin-Topt))*(Tmax-Topt))/(Topt-Tmin)), branch)
									-Topt*lambertw(-((exp(Tmax/(Tmin-Topt)-Topt/(Tmin-Topt))*(Tmax-Topt))/(Topt-Tmin)), branch)
										-Tmax+Topt)/((Tmin-Topt)*(Topt-Tmax))
@inline ratkowski_b(Tmin, Tmax, d, Topt, muopt) = sqrt(muopt/((Tmin-Topt)^2*(exp(d*(Topt-Tmax))-1)^2))
MOBMonodRatkowskiModel(Tmin::Float64, Tmax::Float64, Topt::Float64, Vmax::Float64, Km_ch4::Float64, Km_o2::Float64, y::Float64, baseline_mortality::Float64, s::Symbol) = MOBMonodRatkowskiModel(Tmin, Tmax, Topt, Vmax, Km_ch4, Km_o2, y, baseline_mortality, Val{s})
function MOBMonodRatkowskiModel(Tmin::Float64, Tmax::Float64, Topt::Float64, Vmax::Float64, Km_ch4::Float64, Km_o2::Float64, y::Float64, baseline_mortality::Float64, s::Type{Val{:infer}})
	d = ratkowski_d(Tmin, Tmax, Topt, 0)
	if d == 0
		d = ratkowski_d(Tmin, Tmax, Topt, -1)
	end
	b = ratkowski_b(Tmin, Tmax, d, Topt, Vmax)
	MOBMonodRatkowskiModel(Tmin, Tmax, b, d, Km_ch4, Km_o2, 200000.0, y, baseline_mortality)
end


@inline function ch4_limitation(model::MOBMonodRatkowskiModel, ch4::Array{Float64,1})
	ch4./(model.Km_ch4.+ch4)
end

@inline function o2_limitation(model::MOBMonodRatkowskiModel, o2::Array{Float64,1})
	o2./(model.Km_o2.+o2).*model.KI_o2./(model.KI_o2.+o2)
	#exp.(-0.01.*o2./1000).-exp.(-(0.01+0.18).*o2./1000)
end

@inline function ratkowski(model::MOBMonodRatkowskiModel, T::Array{Float64,1})
	rate = (model.b.*(T.-model.Tmin).*(1.0.-exp.(model.d.*(T.-model.Tmax)))).^2
	rate[T.<model.Tmin] .= 0.0
	rate[T.>model.Tmax] .= 0.0
	return rate
end

@inline mortality(model::MOBMonodRatkowskiModel, ch4::Array{Float64,1}, o2::Array{Float64,1}, temperature::Array{Float64,1}) = model.baseline_mortality.+0.022.*o2./(100000.0.+o2)

@inline mox(model::MOBMonodRatkowskiModel, ch4::Array{Float64,1}, o2::Array{Float64,1}, temperature::Array{Float64,1}) = begin
	mox_rates = (3600.0*24.0).*ratkowski(model, temperature).*ch4_limitation(model, ch4).*o2_limitation(model, o2)
	
	μ = mox_rates.*model.y.-mortality(model, ch4, o2, temperature)
	#μ[:] .= 0.13*1.5e-5*(3600.0*24.0)

	#μ[:] .= 0.0
	o2_rates = -(2.0-model.y).*mox_rates
	#o2_rates[μ.<0.0].=0.0
	ch4_rates = -mox_rates
	#ch4_rates[μ.<0.0].=0.0
	return (μ, ch4_rates, o2_rates)
end

# Thottathil Kinetics
# =============================
struct ThottathilModel <: MOBGrowthModel
	β0::Float64
	β1::Float64
	β2::Float64
	α::Float64
	γ::Float64
end

ThottathilModel() = ThottathilModel(20.08, 0.79, -5669.61, 0.01, 0.18)

mox(model::ThottathilModel, ch4::Array{Float64,1}, o2::Array{Float64,1}, temperature::Array{Float64,1}) = begin
	ln_mox = @. model.β0+model.β1*log(ch4)+model.β2*(1.0/temperature)+log(exp(-model.γ*o2)-exp(-(model.α+model.γ)*o2))
	μ = zeros(length(ch4))
	ch4_rates = -exp.(ln_mox)
	o2_rates = -exp.(ln_mox)
	return (μ, ch4_rates, o2_rates)
end