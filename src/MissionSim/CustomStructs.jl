module CustomStructs

# UPDATE CONSTRUCTORS (with random? Zeros? It is a huge pain rn to simplify)
# Add in STATE constructor (and for all?) that do and dont use SVector


using EarthAlbedo
# include("/home/benjj/.julia/dev/EarthAlbedo.jl/src/EarthAlbedo.jl");  using .EarthAlbedo 
using SatelliteDynamics
using StaticArrays, Distributions, LinearAlgebra

export MAGNETOMETER, DIODES, SATELLITE, SENSORS, ALBEDO, GROUND_TRUTH, ESTIMATES, TRIVIAL, FLAGS, NOISE, STATE

# Add in random constructors? Make immutable?

struct MAGNETOMETER{T}
    scale_factors::SVector{3, T}         # Linear scale factors for soft iron materials (a, b, c)    |   [3,]
    non_ortho_angles::SVector{3, T}      # Offset from purely orthogonal vectors (ρ, λ, ϕ)           |   [3,]
    bias::SVector{3, T}                  # Constant bias on each axis (DC offset) (x₀, y₀, z₀)       |   [3,] 

    function MAGNETOMETER(s::SVector{3, T}, noa::SVector{3, T}, b::SVector{3, T}) where {T}
        """ Primary Constructor """
        new{T}(s, noa, b)
    end

    function MAGNETOMETER(s, noa, b) 
        @warn "Not using static vectors - converting ..."
        T = typeof(s[1])
        MAGNETOMETER(SVector{3, T}(s), SVector{3, T}(noa), SVector{3, T}(b))
    end

    function MAGNETOMETER() 
        """ Randomly generate """
        sf  = 1.0 .+ 0.2 * rand(3) 
        noa = 0.2 * rand(3) 
        b   = 0.25 * rand(3)
        MAGNETOMETER(SVector{3, Float64}(sf), SVector{3, Float64}(noa), SVector{3, Float64}(b))
    end
end

# mutable struct MAGNETOMETER
#     scale_factors::Array{<:Real, 1}         
#     non_ortho_angles::Array{<:Real, 1}       # Offset from purely orthogonal vectors (ρ, λ, ϕ)           |   [3,]
#     bias::Array{<:Real, 1}                   # Constant bias on each axis (DC offset) (x₀, y₀, z₀)       |   [3,] 

#     function MAGNETOMETER(s, a, b) 
#         new(s, a, b)
#     end

#     function MAGNETOMETER() 
#         """ Randomly generate """
#         sf  = 1.0 .+ 0.2 * rand(3) 
#         noa = 0.2 * rand(3) 
#         b   = 0.25 * rand(3)
#         MAGNETOMETER(sf, noa, b)
#     end
# end
Base.deepcopy(s::MAGNETOMETER) = MAGNETOMETER(deepcopy(s.scale_factors), deepcopy(s.non_ortho_angles), deepcopy(s.bias)) #, deepcopy(s.induced_scale_factors))

struct DIODES{S, T}
    calib_values::SVector{S, T}       #   [number of diodes, ]
    azi_angles::SVector{S, T}         #   [number of diodes, ]
    elev_angles::SVector{S, T}        #   [number of diodes, ]

    function DIODES(cv::SVector{S, T}, aa::SVector{S, T}, ea::SVector{S, T}) where {S, T}
        new{S, T}(cv, aa, ea)
    end

    function DIODES(cv, aa, ea)
        @warn "Improper declaration of DIODES; converting to static..."
        S = length(cv)
        T = typeof(cv[1])
        DIODES(SVector{S, T}(cv), SVector{S, T}(aa), SVector{S, T}(ea))
    end

    function DIODES(; N = 6)
        """ Randomly Generate """
        cv = 1.0 .+ 0.15 * rand(N) 

        if N == 6
            ea = [(-pi/4); (pi/4);  0.0;    0.0;    (pi/4); (-pi/4)] 
            aa = [0.0;      pi;     (pi/2); (-pi/2); 0.0;    pi]  
            ea = ea + rand(Normal(0.0, deg2rad(5.0)), N)
            aa = aa + rand(Normal(0.0, deg2rad(5.0)), N)
        else
            ea = deg2rad.(rand(-90:90, N))
            aa = deg2rad.(rand(-179:180, N))
        end

        DIODES(SVector{N, Float64}(cv), SVector{N, Float64}(ea), SVector{N, Float64}(aa))
    end
end

# mutable struct DIODES
#     calib_values::Array{Float64, 1}      #   [number of diodes, ]
#     azi_angles::Array{Float64, 1}        #   [number of diodes, ]
#     elev_angles::Array{Float64, 1}       #   [number of diodes, ]

#     function DIODES(cv, aa, ea)
#         new(cv, aa, ea)
#     end

#     function DIODES(; N = 6)
#         """ Randomly Generate """
#         cv = 1.0 .+ 0.15 * rand(N) 
#         ea = [(-pi/4); (pi/4);  0.0;    0.0;    (pi/4); (-pi/4)] 
#         aa = [0.0;      pi;     (pi/2); (-pi/2); 0.0;    pi]  

#         ea = ea + rand(Normal(0.0, deg2rad(5.0)), N)
#         aa = aa + rand(Normal(0.0, deg2rad(5.0)), N)

#         DIODES(cv, ea, aa)
#     end
# end
Base.deepcopy(s::DIODES) = DIODES(deepcopy(s.calib_values), deepcopy(s.azi_angles), deepcopy(s.elev_angles))


# WHY is state necessary here...?
mutable struct SATELLITE{S, T}
    J::SMatrix{3, 3, T, 9}              # Inertia Matrix of satellite  (array of 1 instead of 2 because diagonal)    |   [3 x 3]
    magnetometer::MAGNETOMETER{T}       # Mag calibration values      |   [3, 3, 3, 3 x number of diodes]
    diodes::DIODES{S, T}                # Diode calibration values 
    # state#::Array{Float64, 1}               # Satellite State             | [q⃗ q₀ β⃗]
    # covariance #::Array{<:Real, 2}

    function SATELLITE(J::SMatrix{3, 3, T, 9}, mag::MAGNETOMETER{T}, dio::DIODES{S, T}) where {S, T}
        new{S, T}(J, mag, dio)
    end

    function SATELLITE(J, mag::MAGNETOMETER, dio::DIODES, state, cov)
        @warn "USING OLD SATELLITE STRUCT!!"
        SATELLITE(J, mag, dio)
    end

    function SATELLITE(; _J = nothing, mag::MAGNETOMETER{T} = MAGNETOMETER(), dio::DIODES{S, T} = DIODES()) where {S, T}
        """ Random SATELLITE """
        m = 1.5 + 0.5 * rand()
        l, w, h = 1 .+0.5 * rand(3)

        if !isnothing(_J)
            return SATELLITE(_J, mag, dio)
        else
            J = zeros(T, 3, 3)
            J[1, 1] = (m / 12.0) * (l^2 + w^2) 
            J[2, 2] = (m / 12.0) * (l^2 + h^2)
            J[3, 3] = (m / 12.0) * (h^2 + w^2)
            
            SATELLITE(SMatrix{3, 3, T, 9}(J), mag, dio)
        end
    end

end
Base.deepcopy(s::SATELLITE) = SATELLITE(deepcopy(s.J), deepcopy(s.magnetometer), deepcopy(s.diodes))#, deepcopy(s.state), deepcopy(s.covariance))

mutable struct SENSORS 
    magnetometer::Array{Float64, 1}  
    diodes::Array{Float64, 1}
    gyro::Array{Float64, 1}
    pos::Array{Float64, 1}           
end

struct NOISE
    diodes::Array{Float64, 1}
    gyro::Array{Float64, 1}
    pos::Array{Float64, 1}
    # sᴮ_rot::Array{Float64, 2}
    # Bᴮ_rot::Array{Float64, 2}
end
    
struct ALBEDO
    refl::REFL
    cell_centers_ecef::Array{<:Real, 3}
end

mutable struct GROUND_TRUTH # Same as sim_results?
    t::Epoch
    Bᴵ::Array{Float64, 1}
    sᴵ::Array{Float64, 1}
    ŝᴮ::Array{Float64, 1}
    Bᴮ::Array{Float64, 1}
    I::Array{Float64, 1}
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

    function STATE(; r = nothing, v = nothing, q = nothing, ω = nothing, β = nothing, a = nothing, Rₑ = 6378.1363e3, μ = 3.9860044188e14)
        """ Random Constructor - Primarily used for testing and convenience """
        
        function rand_quat()
            q = randn(4)
            return q / norm(q)
        end

        _r = isnothing(r) ? (Rₑ + (rand(20:100) * 10^4)) * [1.0, 0.0, 0.0] : r  
        _a = isnothing(a) ? norm(_r) * (1.0 + 0.5 * rand()) : a
        _v = isnothing(v) ? sqrt(μ*( (2/norm(_r)) - (1/_a))) * [0.0, 1.0, 0.0] : v 
        _q = isnothing(q) ? rand_quat() : q 
        _ω = isnothing(ω) ? randn(3) : ω
        _β = isnothing(β) ? randn(3) : β

        STATE(SVector{3, Float64}(_r), SVector{3, Float64}(_v), SVector{4, Float64}(_q), 
                SVector{3, Float64}(_ω), SVector{3, Float64}(_β))
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