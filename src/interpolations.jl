#========================
Linear interpolation
=========================#

# generic interpolator datastructure
struct Interpolation{InterpolationMethod}
	x::Array{<:Real,1}  # 1d domain
	y::Array{<:Real,1}  # 1d data
	at::Function  # interpolator
end

top(i::Interpolation{<:Any}) = i.x[1]
bottom(i::Interpolation{<:Any}) = i.x[end]

# linear interpolator
interp1d(x::Array{T,1}, y::Array{T,1}) where T<:Number = begin
							itp = interpolate((x,), y, Gridded(Linear()))
							(at) -> itp(at)
						end
precompile(interp1d, (Array{Float64,1}, Array{Float64,1}))
interp1d!(x,y) = interp1d(convert(Array{Float64}, x), convert(Array{Float64}, y)) 

# implementation for linear interpolation
const LinearInterpolation = Interpolation{:Linear}
(::Type{Interpolation{:Linear}})(x::Array{<:Real,1}, y::Array{<:Real,1}) = Interpolation{:Linear}(x, y, interp1d(x, y))
(::Type{Interpolation{:Linear}})(x::Array{<:Float64,1}, y::Array{<:Int64,1}) = Interpolation{:Linear}(x, y, interp1d!(x, y))

# linear interpolator
interp2d(x::Array{T,1}, y::Array{T,1}, z::Array{T,1}) where T<:Number = begin
							itp = interpolate((x,y), z, Gridded(Linear()))
							(atx, aty) -> itp(atx, aty)
						end