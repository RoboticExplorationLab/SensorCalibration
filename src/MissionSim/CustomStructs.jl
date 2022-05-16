module CustomStructs

# UPDATE CONSTRUCTORS (with random? Zeros? It is a huge pain rn to simplify)
# Add in STATE constructor (and for all?) that do and dont use SVector


using EarthAlbedo
# include("/home/benjj/.julia/dev/EarthAlbedo.jl/src/EarthAlbedo.jl");  using .EarthAlbedo 
using SatelliteDynamics
using StaticArrays

export MAGNETOMETER, DIODES, SATELLITE, SENSORS, ALBEDO, GROUND_TRUTH, ESTIMATES, TRIVIAL, FLAGS, NOISE, STATE


mutable struct MAGNETOMETER
    scale_factors::Array{<:Real, 1}          # Linear scale factors for soft iron materials (a, b, c)    |   [3,]
    non_ortho_angles::Array{<:Real, 1}       # Offset from purely orthogonal vectors (ρ, λ, ϕ)           |   [3,]
    bias::Array{<:Real, 1}                   # Constant bias on each axis (DC offset) (x₀, y₀, z₀)       |   [3,] 
end
Base.deepcopy(s::MAGNETOMETER) = MAGNETOMETER(deepcopy(s.scale_factors), deepcopy(s.non_ortho_angles), deepcopy(s.bias)) #, deepcopy(s.induced_scale_factors))

mutable struct DIODES
    calib_values::Array{Float64, 1}      #   [number of diodes, ]
    azi_angles::Array{Float64, 1}        #   [number of diodes, ]
    elev_angles::Array{Float64, 1}       #   [number of diodes, ]
end
Base.deepcopy(s::DIODES) = DIODES(deepcopy(s.calib_values), deepcopy(s.azi_angles), deepcopy(s.elev_angles))


# WHY is state necessary here...?
mutable struct SATELLITE
    J::Array{Float64, 2}                   # Inertia Matrix of satellite  (array of 1 instead of 2 because diagonal)    |   [3 x 3]
    magnetometer::MAGNETOMETER            # Mag calibration values      |   [3, 3, 3, 3 x number of diodes]
    diodes::DIODES                        # Diode calibration values 
    state#::Array{Float64, 1}               # Satellite State             | [q⃗ q₀ β⃗]
    covariance #::Array{<:Real, 2}
end
Base.deepcopy(s::SATELLITE) = SATELLITE(deepcopy(s.J), deepcopy(s.magnetometer), deepcopy(s.diodes), deepcopy(s.state), deepcopy(s.covariance))

mutable struct SENSORS 
    magnetometer::Array{Float64, 1}  
    diodes::Array{Float64, 1}
    gyro::Array{Float64, 1}
    gps::Array{Float64, 1}           
end

struct NOISE
    diodes::Array{Float64, 1}
    gyro::Array{Float64, 1}
    gps::Array{Float64, 1}
    sᴮ_rot::Array{Float64, 2}
    Bᴮ_rot::Array{Float64, 2}
end
    
struct ALBEDO
    refl#::refl_struct
    cell_centers_ecef::Array{<:Real, 3}
end

mutable struct GROUND_TRUTH # Same as sim_results?
    t_hist::Epoch
    Bᴵ_hist::Array{Float64, 1}
    sᴵ_hist::Array{Float64, 1}
    ŝᴮ_hist::Array{Float64, 1}
    Bᴮ_hist::Array{Float64, 1}
end


struct TRIVIAL
    junk
end


mutable struct FLAGS
    in_sun::Bool 
    magnetometer_calibrated::Bool
    diodes_calibrated::Bool
    detumbling::Bool
    calibrating::Bool
end

struct STATE{S, T}
    r::SVector{3, T}   # Position (Cartesian)
    v::SVector{3, T}   # Velocity
    q::SVector{4, T}   # Scalar-first quaternion
    ω::SVector{3, T}   # Angular velocity
    β::SVector{3, T}   # Gyroscope bias
    x::SVector{S, T}

    function STATE(r::SVector{3, T}, v::SVector{3, T}, q::SVector{4, T}, ω::SVector{3, T}, β::SVector{3, T}) where {T} 
        S = 16 #length(r) + length(v) + length(q) + length(ω) + length(β)
        x = SVector{S, T}([r; v; q; ω; β])

        new{S, T}(r, v, q, ω, β, x)
    end

    function STATE(r::Vector{T}, v::Vector{T}, q::Vector{T}, ω::Vector{T}, β::Vector{T}) where {T}
        @warn "Initializing STATE without static arrays! Automatically converting..."
        _r = SVector{3, T}(r)
        _v = SVector{3, T}(v)
        _q = SVector{4, T}(q)
        _ω = SVector{3, T}(ω)
        _β = SVector{3, T}(β)

        STATE(_r, _v, _q, _ω, _β)  # Call the default constructor
    end
end
    # Define addition, multiplication, and subtraction?  (NEEDED for RK4)
function Base.:+(x₁::STATE, x₂::STATE) 
    r = x₁.r + x₂.r 
    v = x₁.v + x₂.v 
    q = x₁.q + x₂.q 
    ω = x₁.ω + x₂.ω 
    β = x₁.β + x₂.β 

    return STATE(r, v, q, ω, β)
end

function Base.:*(k::Real, x::STATE)
    r = k * x.r 
    v = k * x.v  
    q = k * x.q  
    ω = k * x.ω  
    β = k * x.β  

    return STATE(r, v, q, ω, β)
end

function Base.:/(x::STATE, k::Real)
    r = x.r / k
    v = x.v / k
    q = x.q / k  
    ω = x.ω / k
    β = x.β / k 

    return STATE(r, v, q, ω, β)
end


end