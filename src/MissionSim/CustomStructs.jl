module CustomStructs

""" TO Do
 - Decide on default photodiode configuration (in diodes and sat state)
 - Add in plots, magnetometer noise to the 'NOISE' struct...?
 - Alphabetize? Organize somehow (in separate files?)
 - flags
"""

using EarthAlbedo
# include("/home/benjj/.julia/dev/EarthAlbedo.jl/src/EarthAlbedo.jl");  using .EarthAlbedo 
using SatelliteDynamics
using StaticArrays, Distributions, LinearAlgebra
using Plots

export MAGNETOMETER, DIODES, SATELLITE, SENSORS, ALBEDO, GROUND_TRUTH, ESTIMATES, TRIVIAL, FLAGS, NOISE, STATE
export SAT_STATE, SAT_COVARIANCE
export update_state

                                ##########################################################################
                                #                              MAGNETOMETER                              # 
                                ##########################################################################
"""
    MAGNETOMETER{T} -> scale_factors, non_ortho_angles, bias

      Magnetometer sensor with that includes scale factors along each axis, non-orthogonality angles
    (angular distance from perfectly orthogonal axes), and a constant bias. 

      Includes a primary constructor that allows for static and non-static inputs, as well as a  
    generator that can make an ideal magnetometer or a random noisy one. 
"""
struct MAGNETOMETER{T}
    scale_factors::SVector{3, T}         # Linear scale factors for soft iron materials (a, b, c)    |   [3,]
    non_ortho_angles::SVector{3, T}      # Offset from purely orthogonal vectors (ρ, λ, ϕ)           |   [3,]
    bias::SVector{3, T}                  # Constant bias on each axis (DC offset) (x₀, y₀, z₀)       |   [3,] 

    function MAGNETOMETER(s::SVector{3, T}, noa::SVector{3, T}, b::SVector{3, T}) where {T}
        """ Primary Constructor """
        new{T}(s, noa, b)
    end

    function MAGNETOMETER(s, noa, b) 
        """ Converts to static """
        T = typeof(s[1])
        MAGNETOMETER(SVector{3, T}(s), SVector{3, T}(noa), SVector{3, T}(b))
    end

    function MAGNETOMETER(; ideal = false) 
        """ Randomly generate """

        if ideal # Generate a perfect magnetometer
            sf  = SVector{3, Float64}(ones(3))   # No scale 
            noa = SVector{3, Float64}(zeros(3))  # No offset 
            b   = SVector{3, Float64}(zeros(3))  # No Bias
            return MAGNETOMETER(sf, noa, b);
        else # Make a noisy one
            sf  = 1.0 .+ 0.2 * rand(3) 
            noa = 0.2 * rand(3) 
            b   = 0.25 * rand(3)
            return MAGNETOMETER(SVector{3, Float64}(sf), SVector{3, Float64}(noa), SVector{3, Float64}(b))
        end
    end
end
Base.deepcopy(s::MAGNETOMETER) = MAGNETOMETER(deepcopy(s.scale_factors), deepcopy(s.non_ortho_angles), deepcopy(s.bias)) #, deepcopy(s.induced_scale_factors))







                                ##########################################################################
                                #                                 DIODES                                 # 
                                ##########################################################################
"""
    DIODES{S, T} -> calib_values, azi_angles, elev_angles

      Calibration values (~ scale factor) and surface normal (in azimuth and elevation angle) 
    for a set of sun-sensing photodiodes on the CubeSat. Allows for any number of diodes, but 
    the default is six, with one located on each side of the CubeSat. 
"""
struct DIODES{S, T}
    calib_values::SVector{S, T}       #   Calibration value (How much 1 unit of current is scaled by)
    azi_angles::SVector{S, T}         #   Azimuth angle for each photodiode surface normal
    elev_angles::SVector{S, T}        #   Elevation angle for each photodiode surface normal

    function DIODES(cv::SVector{S, T}, aa::SVector{S, T}, ea::SVector{S, T}) where {S, T}
        """ Primary Constructor """
        new{S, T}(cv, aa, ea)
    end

    function DIODES(cv, aa, ea)
        """ Constructor for non-static types """
        @warn "Improper declaration of DIODES; converting to static..."

        S = length(cv)
        T = typeof(cv[1])
        DIODES(SVector{S, T}(cv), SVector{S, T}(aa), SVector{S, T}(ea))
    end

    function DIODES(; ideal = false, N = 6)
        """ 
              Randomly generates a set of diodes. If using six diodes, the default orientation is
            set up and perturbed slightly. If some other number of diodes is used, each is assigned
            a completely random surface normal. Also allows for an 'ideal' set of diodes with perfect 
            scale and surface normals exactly as we would want.
        """

        if ideal 
            cv = ones(6)
            # ea = [(-pi/4); (pi/4);  0.0;    0.0;    (pi/4); (-pi/4)] 
            # aa = [0.0;      pi;     (pi/2); (-pi/2); 0.0;    pi]  
            ea = []

            aa = [0.0,  pi, pi/2, -pi/2,  0.0,   0.0]
            ea = [0.0, 0.0,  0.0,   0.0, pi/2, -pi/2]

            return DIODES(SVector{N, Float64}(cv), SVector{N, Float64}(ea), SVector{N, Float64}(aa))
        end

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

        return DIODES(SVector{N, Float64}(cv), SVector{N, Float64}(ea), SVector{N, Float64}(aa))
    end
end
Base.deepcopy(s::DIODES) = DIODES(deepcopy(s.calib_values), deepcopy(s.azi_angles), deepcopy(s.elev_angles))







                                ##########################################################################
                                #                                SAT STATE                               # 
                                ##########################################################################
"""
    SAT_STATE{N, T} -> q, β, C, α, ϵ

      State of the SATELLITE (not to be confused with the environment state 'state'). Comprised of 
    the attitude 'q' (as a scalar-first unit quaternion) and the gyroscope bias 'β', as well as the 
    diode calibration parameters: Calibration value 'C', and the azimuth 'α' and elevation 'ϵ' angles 
    for the surface normals of each diode.
"""
struct SAT_STATE{N, T}
    q::SVector{4, T}     # Scalar-first unit quaternion representing orientation of the satellite 
    β::SVector{3, T}     # Gyroscope bias 
    C::SVector{N, T}     # Current scale factor for each diode 
    α::SVector{N, T}     # Azimuth angle corresponding to each diode 
    ϵ::SVector{N, T}     # Elevation angle corresponding to each diode 

    function SAT_STATE(q::SVector{4, T}, β::SVector{3, T}, C::SVector{N, T}, α::SVector{N, T}, ϵ::SVector{N, T}) where {N, T}
        """ Primary Constructor """
        new{N, T}(q, β, C, α, ϵ)
    end

    function SAT_STATE(; q = nothing, β = nothing, C = nothing, α = nothing, ϵ = nothing, N = 6, ideal = false)
        """ Randomly generate either a noisy or ideal satellite state """

        function rand_quat()
            t = randn(4)
            return t / norm(t)
        end

        if ideal  # No rotation, no bias, unit scale factors, and perfect 

            if N != 6 
                @warn "Don't know how to make an ideal state if N != 6!"
            end

            _q = SVector{4, Float64}(1.0, 0, 0, 0)
            _β = SVector{3, Float64}(0.0, 0.0, 0.0)
            _C = SVector{6, Float64}(ones(N))
            _α = SVector{6, Float64}(0.0, deg2rad(180), deg2rad(90), deg2rad(-90), 0.0, 0.0)
            _ϵ = SVector{6, Float64}(0.0, 0.0, 0.0, 0.0, deg2rad(90), deg2rad(-90))

            return SAT_STATE(_q, _β, _C, _α, _ϵ)
        else

            _q = isnothing(q)  ?  rand_quat()  :  q 
            _β = isnothing(β)  ?  randn(3)     :  β 
            _C = isnothing(C)  ?  1.0 .+ 0.1 * randn(N)  :  C 
            _α = isnothing(α)  ?  deg2rad.(rand(0:359, N))  :  α 
            _ϵ = isnothing(ϵ)  ?  deg2rad.(rand(0:179, N))  : ϵ
    
            return SAT_STATE(SVector{4, Float64}(_q),  SVector{3, Float64}(_β),  SVector{N, Float64}(_C),
                                SVector{N, Float64}(_α),  SVector{N, Float64}(_ϵ))
        end
    end
end
Base.deepcopy(ss::SAT_STATE) = SAT_STATE(; q = deepcopy(ss.q), β = deepcopy(ss.β), C = deepcopy(ss.C), α = deepcopy(ss.α), ϵ = deepcopy(ss.ϵ))

"""
    update_state(x; q, β, C, α, ϵ)

      'Updates' portions of a SAT_STATE struct by creating a new struct 
    and copying relevant stuff over. 
"""
function update_state(x::SAT_STATE{N, T}; q = nothing, β = nothing, C = nothing, α = nothing, ϵ = nothing) where {N, T}

    _q = isnothing(q) ? copy(x.q) : q 
    _β = isnothing(β) ? copy(x.β) : β
    _C = isnothing(C) ? copy(x.C) : C
    _α = isnothing(α) ? copy(x.α) : α
    _ϵ = isnothing(ϵ) ? copy(x.ϵ) : ϵ

    return SAT_STATE(_q, _β, _C, _α, _ϵ)
end







                                ##########################################################################
                                #                             SAT COVARIANCE                             # 
                                ##########################################################################
""" Originally made so I could define functions on it, but I am not sure it is worht it..."""
# Unnecessary?
struct SAT_COVARIANCE{T}
    """ Kept together b/c may have cross terms """
    """ Probably lower triangular, but not necessarily forced to be so """
    # Too large for Static, (> 100)
    Σ::Matrix{T}   # Covariance matrix for satellite state
    N::Int         # Number of diodes

    function SAT_COVARIANCE(ϕ, β, C, α, ϵ)
        N = size(C, 1)
        T = typeof(C[1])

        ℓ = 6 + 3 * N 
        Σ = zeros(T, ℓ, ℓ)
        Σ[1:3, 1:3] .= ϕ
        Σ[4:6, 4:6] .= β 

        i₀ = 7; i = i₀ - 1 + N
        Σ[i₀:i, i₀:i]   .= C 

        i₀ = i + 1; i = i₀ - 1 + N
        Σ[i₀:i, i₀:i] .= α 

        i₀ = i + 1; i = i₀ - 1 + N
        Σ[i₀:i, i₀:i] .= ϵ

        new{T}(Σ, N)
    end

    ## NOISE RANDOM!
    function SAT_COVARIANCE(; σϕ = deg2rad(10), σβ = deg2rad(10), σC = 0.1, σα = deg2rad(3), σϵ = deg2rad(3), N = 6)
        """ random """
        Σϕ = diagm( (σϕ^2) * ones(3) )
        Σβ = diagm( (σβ^2) * ones(3) )
        ΣC = diagm( (σC^2) * ones(N) )
        Σα = diagm( (σα^2) * ones(N) )
        Σϵ = diagm( (σϵ^2) * ones(N) )
        
        SAT_COVARIANCE(Σϕ, Σβ, ΣC, Σα, Σϵ)
    end
end
function ϕ(cov::SAT_COVARIANCE{T}) where {T}
    return SMatrix{4, 4, T, 16}(cov.Σ[1:3, 1:3])
end
function β(cov::SAT_COVARIANCE{T}) where {T}
    return SMatrix{3, 3, T, 9}(cov.Σ[4:6, 4:6])
end
function C(cov::SAT_COVARIANCE{T}) where {T}
    N = cov.N
    i₀ = 4
    i  = i₀ + N - 1
    return SMatrix{N, N, T, N * N}(cov.Σ[i₀:i, i₀:i])
end
function α(cov::SAT_COVARIANCE{T}) where {T}
    N = cov.N
    i₀ = 4 + N
    i  = i₀ + N - 1
    return SMatrix{N, N, T, N * N}(cov.Σ[i₀:i, i₀:i])
end
function ϵ(cov::SAT_COVARIANCE{T}) where {T}
    N = cov.N
    i₀ = 4 + 2 * N
    i  = i₀ + N - 1
    return SMatrix{N, N, T, N * N}(cov.Σ[i₀:i, i₀:i])
end







                                ##########################################################################
                                #                               SATELLITE                                # 
                                ##########################################################################
"""
    SATELLITE{S, T} -> J, magnetometer, diodes, state, covariance

      Satellite struct that contains the inertia matrix for the CubeSat, as well
    as magnetometers and diodes. Also tracks the SATELLITE state and covariance 
    (not the environment state state)
"""
struct SATELLITE{S, T}
    J::SMatrix{3, 3, T, 9}               # Inertia Matrix of satellite  
    magnetometer::MAGNETOMETER{T}        # Mag calibration values      
    diodes::DIODES{S, T}                 # Diode calibration values 
    state::SAT_STATE{S, T}               # Satellite state 
    covariance::Matrix{T}                # Satellite state covariance 

    function SATELLITE(J::SMatrix{3, 3, T, 9}, mag::MAGNETOMETER{T}, dio::DIODES{S, T}, 
                        sta::SAT_STATE{S, T}, cov::Matrix{T})  where {S, T}
        """ Primary Constructor """
        new{S, T}(J, mag, dio, sta, cov)
    end

    function SATELLITE(; J = nothing, mag::MAGNETOMETER{T} = MAGNETOMETER(), dio::DIODES = DIODES(), 
                            sta::SAT_STATE = SAT_STATE(; N = 6), cov::Matrix{T} = SAT_COVARIANCE().Σ, ideal::Bool = false) where {T}
        """ Random SATELLITE """

        if ideal 
            _J = SMatrix{3, 3, Float64, 9}(0.2 * I(3))
            _mag = MAGNETOMETER(; ideal = true)
            _dio = DIODES(; ideal = true)
            _sta = SAT_STATE(; ideal = true)
            _cov = I(24)  # 6 + 3 * (N = 6)

            return SATELLITE(_J, _mag, _dio, _sta, _cov)
        end


        if !isnothing(J)
            return SATELLITE(J, mag, dio, sta, cov)
        else
            m = 1.5 + 0.5 * rand()
            l, w, h = 1 .+ 0.5 * rand(3)
    
            J = zeros(T, 3, 3)
            J[1, 1] = (m / 12.0) * (l^2 + w^2) 
            J[2, 2] = (m / 12.0) * (l^2 + h^2)
            J[3, 3] = (m / 12.0) * (h^2 + w^2)
            
            return SATELLITE(SMatrix{3, 3, T, 9}(J), mag, dio, sta, cov)
        end
    end

end
Base.deepcopy(s::SATELLITE) = SATELLITE(deepcopy(s.J), deepcopy(s.magnetometer), deepcopy(s.diodes), deepcopy(s.state), deepcopy(s.covariance))







                                ##########################################################################
                                #                                SENSORS                                 # 
                                ##########################################################################
"""
    SENSORS{N, T} -> magnetometer, diodes, gyro, pos

      Struct containing the measurements for each of the satellite sensors 
    (magnetometer, diodes, gyroscope, and something for position). Also 
    comes with a custom plotting recipe. 
"""
struct SENSORS{N, T}
    magnetometer::SVector{3, T}      # Measured magnetometer reading
    diodes::SVector{N, T}            # Measured photodiode currents (per diode)
    gyro::SVector{3, T}              # Measured angular velocity
    pos::SVector{3, T}               # Measured position 
end

function RecipesBase.plot(s::Vector{SENSORS}; split = false, kwargs...)
    N = size(s, 1)

    # Split apart and format as matrices for plotting
    mags = [s[i].magnetometer for i = 1:N]; mags = hcat(mags...)';
    dios = [s[i].diodes       for i = 1:N]; dios = hcat(dios...)';
    gyrs = [s[i].gyro         for i = 1:N]; gyrs = hcat(gyrs...)';
    poss = [s[i].pos          for i = 1:N]; poss = hcat(poss...)';

    # Make the plots
    pM = plot(mags, title = "Magnetometers", xlabel = "Index", ylabel = "Mag Field (μT)",   label = ["x" "y" "z"]; kwargs...)
    pD = plot(dios, title = "Diodes",        xlabel = "Index", ylabel = "Current (A?)",                        ; kwargs...)
    pG = plot(gyrs, title = "Gyroscope",     xlabel = "Index", ylabel = "Ang Vel (rad/s)", label = ["x" "y" "z"]; kwargs...)
    pP = plot(poss, title = "Position",      xlabel = "Index", ylabel = "Position (m)",    label = ["x" "y" "z"]; kwargs...)

    if split  # Return each plot separately
        return pM, pD, pG, pP
    else      # Return one plot containing subplots
        return plot(pM, pD, pG, pP, plot_title = "Sensors", layout = (2, 2))
    end
end




                                ##########################################################################
                                #                                 NOISE                                  # 
                                ##########################################################################
"""
    NOISE{T} -> diodes, gyro, pos 

      Struct containing the noise in each sensor measurement, except for 
    some resaon I don't have magnetometer noise. This is really only used 
    for debugging anyway, so...

    I may add in plotting for this at some point...
"""
struct NOISE{N, T}
    diodes::SVector{N, T}       # Noise in each measured diode current
    gyro::SVector{3, T}         # Noise in the gyroscope measurement
    pos::SVector{3, T}          # Noise in the position estimate 
end




                                ##########################################################################
                                #                                ALBEDO                                  # 
                                ##########################################################################
"""
    ALBEDO -> refl, cell_centers_ecef 

      Simple struct containing an EarthAlbedo.REFL struct (which contains 
    data, type, start_time, and stop_time), as well as the location of the
    center of each cell in Earth-Centered Earth-Fixed (ECEF) frame
"""
struct ALBEDO
    refl::REFL
    cell_centers_ecef::Array{<:Real, 3}
end







                                ##########################################################################
                                #                             GROUND TRUTH                               # 
                                ##########################################################################
"""
    GROUND_TRUTH{N, T} -> t, Bᴵ, sᴵ, ŝᴮ, Bᴮ, i

      Struct containing the ground truth of various things at a 
    given time stamp. Used to validate the performance of the algorithm, 
    and for aid in debugging. 

    Comes with a custom plotting function.
"""
struct GROUND_TRUTH{N, T}
    t::Epoch                        # Time stamp at which this data holds
    Bᴵ::SVector{3, T}               # Magnetic field vector in inertial frame
    sᴵ::SVector{3, T}               # Sat-Sun vector in inertial frame
    ŝᴮ::SVector{3, T}               # Unit Sat-Sun vector in body frame
    Bᴮ::SVector{3, T}               # Magnetic field vector in body frame
    I::SVector{N, T}                # Diode currents 
end
function RecipesBase.plot(gt::Vector{GROUND_TRUTH}; split = false, kwargs...)
    N = size(gt, 1)
    ts  = [gt[i].t - gt[1].t for i = 1:N];  # Convert to seconds 
    
    # Split apart and format as matrices for plotting
    Bᴵs = [gt[i].Bᴵ for i = 1:N]; Bᴵs = hcat(Bᴵs...)';
    sᴵs = [gt[i].sᴵ for i = 1:N]; sᴵs = hcat(sᴵs...)';
    ŝᴮs = [gt[i].ŝᴮ for i = 1:N]; ŝᴮs = hcat(ŝᴮs...)';
    Bᴮs = [gt[i].Bᴮ for i = 1:N]; Bᴮs = hcat(Bᴮs...)';
    Is  = [gt[i].I  for i = 1:N]; Is  = hcat(Is... )';

    # Generate plots
    pBᴵ = plot(ts, Bᴵs, title = "Inertial Mag Vec (Bᴵ)", xlabel = "Time (s)", ylabel = "Strength (μT)",  label = ["x" "y" "z"]; kwargs...)
    pBᴮ = plot(ts, Bᴮs, title = "Body Mag Vec (Bᴮ)",     xlabel = "Time (s)", ylabel = "Strength (μT)",  label = ["x" "y" "z"]; kwargs...)
    psᴵ = plot(ts, sᴵs, title = "Inertial Sun Vec (sᴵ)", xlabel = "Time (s)", ylabel = "Distance (m)",   label = ["x" "y" "z"]; kwargs...)
    psᴮ = plot(ts, ŝᴮs, title = "Body Sun Vec (ŝᴮ)",     xlabel = "Time (s)", ylabel = "Relative Scale", label = ["x" "y" "z"]; kwargs...)
    pI  = plot(ts, Is,  title = "Currents (I)",          xlabel = "Time (s)", ylabel = "Current (A)"; kwargs...)
    p   = plot()

    if split   # Return each plot separately
        return pBᴵ, pBᴮ, psᴵ, psᴮ, pI
    else       # Return as one big plot
        return plot(pBᴵ, pBᴮ, psᴵ, psᴮ, pI, p, plot_title = "Ground Truth", layout = (3, 2))
    end
end







                                ##########################################################################
                                #                                 STATE                                  # 
                                ##########################################################################
"""
    STATE{S, T} -> r, v, q, ω, β, x 

      Current state of the environment, used for propagating dynamics 
    and generating measurements. Not to be confused with the SATELLITE state
    SAT_STATE. 

    Comes with a custom plotting recipe, as well as a few basic functions 
    for addition, multiplication, and division of states.
"""
struct STATE{T}
    r::SVector{3, T}   # Position (Cartesian)
    v::SVector{3, T}   # Velocity
    q::SVector{4, T}   # Scalar-first quaternion
    ω::SVector{3, T}   # Angular velocity
    β::SVector{3, T}   # Gyroscope bias


    function STATE(r::SVector{3, T}, v::SVector{3, T}, q::SVector{4, T}, ω::SVector{3, T}, β::SVector{3, T}) where {T} 
        """ Primary Constructor """
        new{T}(r, v, q, ω, β)
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

""" Converts a struct to a vector """
function x(s::STATE{T}) where {T}
    vec = [s.r; s.v; s.q; s.ω; sβ]
    return SVector{16, T}(vec)
end

#NOT efficient, but whatever
function RecipesBase.plot(s::Vector{STATE{T}}; split = false, kwargs...) where {T}
    N = size(s, 1)

    # Split apart and format as matrices for plotting
    rs = [vcat([s[i].r;]...) for i = 1:N]; rs = hcat(rs...)';
    vs = [vcat([s[i].v;]...) for i = 1:N]; vs = hcat(vs...)';
    qs = [vcat([s[i].q;]...) for i = 1:N]; qs = hcat(qs...)';  # THIS LINE CAUSES WEIRD ERRORS
    ωs = [vcat([s[i].ω;]...) for i = 1:N]; ωs = hcat(ωs...)';
    βs = [vcat([s[i].β;]...) for i = 1:N]; βs = hcat(βs...)';

    pr = plot(rs, title = "Position (Cart)", xlabel = "Index", ylabel = "Position (m)",    label = ["x" "y" "z"]; kwargs...)
    pv = plot(vs, title = "Velocity",        xlabel = "Index", ylabel = "Velocity (m/s)",  label = ["x" "y" "z"]; kwargs...)
    pq = plot(qs, title = "Attitude (quat)", xlabel = "Index",                             label = ["x" "y" "z"]; kwargs...)
    pω = plot(ωs, title = "Ang Velocity",    xlabel = "Index", ylabel = "Ang Vel (rad/s)", label = ["x" "y" "z"]; kwargs...)
    pβ = plot(βs, title = "Gyro Bias",       xlabel = "Index", ylabel = "Bias (rad/s)",    label = ["x" "y" "z"]; kwargs...)
    p   = plot()

    if split 
        return pr, pv, pq, pω, pβ
    else
        return plot(pr, pv, pq, pω, pβ, p, plot_title = "State", layout = (3, 2))
    end
end

""" Normalizes a quaternion """
function normalize_quat(x::STATE)
    q = x.q / norm(x.q)
    return STATE(x.r, x.v, q, x.ω, x.β)
end

# Define addition, multiplication, and subtraction?  (NEEDED for RK4)
function Base.:+(x₁::STATE, x₂::STATE) 
    r = x₁.r + x₂.r 
    v = x₁.v + x₂.v 
    q = x₁.q + x₂.q 
    ω = x₁.ω + x₂.ω 
    β = x₁.β + x₂.β 

    # q /= norm(q)
    return STATE(r, v, q, ω, β)
end

function Base.:*(k::Real, x::STATE)
    r = k * x.r 
    v = k * x.v  
    q = k * x.q  
    ω = k * x.ω  
    β = k * x.β  

    # q /= norm(q)
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







                                ##########################################################################
                                #                                 FLAGS                                  # 
                                ##########################################################################
mutable struct FLAGS
    in_sun::Bool 
    magnetometer_calibrated::Bool
    diodes_calibrated::Bool
    detumbling::Bool
    calibrating::Bool
end

end # Module