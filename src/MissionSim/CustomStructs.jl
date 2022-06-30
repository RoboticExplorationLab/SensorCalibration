module CustomStructs
 

""" TO Do
 - Decide on default photodiode configuration (in diodes and sat state)
 - Add in plots, magnetometer noise to the 'NOISE' struct...?
 - Alphabetize? Organize somehow (in separate files?)
 - flags

 - SATELLITE has diodes in SAT_STATE && in DIODES and that can be way problematic (make 'diodes' a function call that pulls from state?)
 - SAT_STATE and DIODES both have the C, α, ϵ -> remove from sat_state


 try reduce(hcat, v) instead of (hcat([v;]...)) in plotting!

    Make a 'Structs/types' folder and have each one in its own, so we can add more functions cleanly
"""

using EarthAlbedo
using SatelliteDynamics
using StaticArrays, Distributions, LinearAlgebra
using Plots

# include("/home/benjj/.julia/dev/EarthAlbedo.jl/src/EarthAlbedo.jl");  using .EarthAlbedo 
using Infiltrator

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
            sf  = rand(Normal(1.0, 0.1), 3) 
            noa = rand(Normal(0.0, deg2rad(3.0)), 3) 
            b   = rand(Normal(0.0, deg2rad(4.0)), 3)
            return MAGNETOMETER(SVector{3, Float64}(sf), SVector{3, Float64}(noa), SVector{3, Float64}(b))
        end
    end
end
# Base.deepcopy(s::MAGNETOMETER) = MAGNETOMETER(deepcopy(s.scale_factors), deepcopy(s.non_ortho_angles), deepcopy(s.bias))







                                ##########################################################################
                                #                                 DIODES                                 # 
                                ##########################################################################
"""
    DIODES{S, T} -> calib_values, azi_angles, elev_angles

      Calibration values (~ scale factor) and surface normal (in azimuth and elevation angle) 
    for a set of sun-sensing photodiodes on the CubeSat. Allows for any number of diodes, but 
    the default is six, with one located on each side of the CubeSat. 
"""

@info "Add a plot for diodes (and for each struct)"

struct DIODES{S, T}
    calib_values::SVector{S, T}       #   Calibration value (How much 1 unit of current is scaled by)
    azi_angles::SVector{S, T}         #   Azimuth angle for each photodiode surface normal
    elev_angles::SVector{S, T}        #   Elevation angle for each photodiode surface normal

    function DIODES(cv::SVector{S, T}, aa::SVector{S, T}, ea::SVector{S, T}) where {S, T}

        # Wrap α ∈ (-pi, pi] and ϵ ∈ (-pi/2, pi/2]
        if any(ea .≤ -pi/2) || any(ea .> pi/2)
            # Convert to cartesian, and convert back to stay in range

            x, y, z = sph2cart(aa, ea)
            aa, ea = cart2sph(x, y, z)
            aa = SVector{S, T}(aa)
            ea = SVector{S, T}(ea)
        elseif any(aa .≤ -pi) || any(aa .> pi)
            aa = wrap(aa)
        end
        
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
            ea = [(-pi/4); (pi/4);  0.0;    0.0;    (pi/4); (-pi/4)] 
            aa = [0.0;      pi;     (pi/2); (-pi/2); 0.0;    pi]  

            # 
            # aa = [0.0,  pi, pi/2, -pi/2,  0.0,   0.0]
            # ea = [0.0, 0.0,  0.0,   0.0, pi/2, -pi/2]

            return DIODES(SVector{N, Float64}(cv), SVector{N, Float64}(aa), SVector{N, Float64}(ea))
        end

        # cv = rand(Normal(2.0, 0.3), N)   # On-Orbit... paper has 0.5 or 0.2 
        cv = rand(Normal(1.0, 0.2), N)

        if N == 6
            aa = [0.0;        pi;   (pi/2); (-pi/2);      0.0;      pi]  
            ea = [(-pi/4); (pi/4);     0.0;     0.0;   (pi/4); (-pi/4)] 


            # aa = [0.0,  pi, pi/2, -pi/2,  0.0,   0.0]
            # ea = [0.0, 0.0,  0.0,   0.0, pi/2, -pi/2]

            ea = ea + rand(Normal(0.0, deg2rad(5.0)), N)   # On-Orbit paper has 10 or 2 deg
            aa = aa + rand(Normal(0.0, deg2rad(5.0)), N)
        else
            ea = deg2rad.(rand(-90:90, N))
            aa = deg2rad.(rand(-179:180, N))
        end

        return DIODES(SVector{N, Float64}(cv), SVector{N, Float64}(aa), SVector{N, Float64}(ea))
    end
end
# Base.deepcopy(s::DIODES) = DIODES(deepcopy(s.calib_values), deepcopy(s.azi_angles), deepcopy(s.elev_angles))
# Base.deepcopy_internal(s::DIODES, dict::IdDict) = DIODES(deepcopy(s.calib_values), deepcopy(s.azi_angles), deepcopy(s.elev_angles))


function wrap(a::SVector{N, T}) where {N, T}
    v = atan.(sin.(a), cos.(a))

    # v = zeros(T, N)
    # for i = 1:N
    #     aᵢ = a[i]
    #     while aᵢ > pi
    #         aᵢ -= (2 * pi)
    #     end 
    #     while aᵢ ≤ -pi 
    #         aᵢ += (2 * pi)
    #     end
    #     v[i] = aᵢ
    # end

    return SVector{N, T}(v)
end

function sph2cart(α, ϵ)
    θ = (pi/2) .- ϵ  # Convert to from elevation to inclination or whatever
    x = cos.(α) .* sin.(θ)
    y = sin.(α) .* sin.(θ)
    z = cos.(θ)
    return x, y, z
end

function cart2sph(x, y, z)
    α = atan.(y, x)
    θ = atan.(sqrt.( (x.^2) .+ (y.^2)), z)  # Assumes unit 
    ϵ = pi/2 .- θ
    return α, ϵ
end





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
struct SAT_STATE{T}
    q::SVector{4, T}     # Scalar-first unit quaternion representing orientation of the satellite 
    β::SVector{3, T}     # Gyroscope bias 
    # C::SVector{N, T}     # Current scale factor for each diode 
    # α::SVector{N, T}     # Azimuth angle corresponding to each diode 
    # ϵ::SVector{N, T}     # Elevation angle corresponding to each diode 

    function SAT_STATE(q::SVector{4, T}, β::SVector{3, T}) where {T}
        """ Primary Constructor """

        q = (norm(q) == 1) ? q  : (q / norm(q))
        q = (q[1] < 0)  ?  -q : q                   # Because quaternions double cover, we force the scalar to be positive
        new{T}(q, β)
    end

    function SAT_STATE(; q = nothing, β = nothing, ideal = false)
        """ Randomly generate either a noisy or ideal satellite state """

        function rand_quat()
            t = randn(4)
            return t / norm(t)
        end

        if ideal  # No rotation, no bias, unit scale factors, and perfect 

            _q = SVector{4, Float64}(1.0, 0, 0, 0)
            _β = SVector{3, Float64}(0.0, 0.0, 0.0)

            return SAT_STATE(_q, _β)
        else

            _q = isnothing(q)  ?  rand_quat()  :  q 
            _β = isnothing(β)  ?  randn(3)     :  β 
    
            return SAT_STATE(SVector{4, Float64}(_q),  SVector{3, Float64}(_β))
        end
    end
end
# Base.deepcopy(ss::SAT_STATE) = SAT_STATE(; q = deepcopy(ss.q), β = deepcopy(ss.β))

"""
    update_state(x; q, β, C, α, ϵ)

      'Updates' portions of a SAT_STATE struct by creating a new struct 
    and copying relevant stuff over. 
"""
function update_state(x::SAT_STATE{T}; q = nothing, β = nothing) where {N, T}

    _q = isnothing(q) ? Array(x.q)  : Array(q) 
    _β = isnothing(β) ? Array(x.β)  : Array(β)
    # _q = isnothing(q) ? deepcopy( SVector{4, T}(vcat([x.q;]...))) : q 
    # _β = isnothing(β) ? SVector{3, T}(deepcopy(vcat([x.β;]...))) : SVector{3, T}(deepcopy(vcat([β;]...)))

    return SAT_STATE(SVector{4, T}(_q), SVector{3, T}(_β))
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
        Tt = typeof(ϕ[1])

        ℓ = 6 + 3 * N 
        Σ = zeros(Tt, ℓ, ℓ)
        Σ[1:3, 1:3] .= ϕ
        Σ[4:6, 4:6] .= β 

        i₀ = 7; i = i₀ - 1 + N
        Σ[i₀:i, i₀:i]   .= C 

        i₀ = i + 1; i = i₀ - 1 + N
        Σ[i₀:i, i₀:i] .= α 

        i₀ = i + 1; i = i₀ - 1 + N
        Σ[i₀:i, i₀:i] .= ϵ

        Σ = cholesky(Hermitian(Matrix(Σ))).U

        new{Tt}(Σ, N)
    end

    ## Not sure if noise is right!
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
    state::SAT_STATE{T}                  # Satellite state 
    covariance::Matrix{T}                # Satellite state covariance 

    function SATELLITE(J::SMatrix{3, 3, T, 9}, mag::MAGNETOMETER{T}, dio::DIODES{S, T}, 
                        sta::SAT_STATE{T}, cov::Matrix{T})  where {S, T}
        """ Primary Constructor  """
        new{S, T}(SMatrix{3, 3, T, 9}(J), mag, dio, sta, cov)
    end

    function SATELLITE(; J = nothing, mag::MAGNETOMETER{T} = MAGNETOMETER(), dio::DIODES = DIODES(), 
                            sta::SAT_STATE = SAT_STATE(), cov::Matrix{T} = SAT_COVARIANCE().Σ, ideal::Bool = false) where {T}
        """ Random SATELLITE """

        if ideal 
            _J = isnothing(J) ? SMatrix{3, 3, Float64, 9}(0.2 * I(3)) : J
            _mag = MAGNETOMETER(; ideal = true)
            _dio = DIODES(; ideal = true)
            _sta = SAT_STATE(; ideal = true)
            _cov = cholesky(Hermitian(Matrix(cov))).U #Matrix(1.0 * I(24))  # 6 + 3 * (N = 6)

            return SATELLITE(_J, _mag, _dio, _sta, Matrix(_cov))
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
# Base.deepcopy(s::SATELLITE) = SATELLITE(deepcopy(s.J), deepcopy(s.magnetometer), deepcopy(s.diodes), deepcopy(s.state), deepcopy(s.covariance))







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

""" Add info! """
function RecipesBase.plot(s::Vector{SENSORS{6, T}}, sensor::Symbol = :a; kwargs...) where {T}
    N = size(s, 1)

    # Split apart and format as matrices for plotting
    if sensor == :a
        mags = [s[i].magnetometer for i = 1:N]; mags = reduce(hcat, mags)'; # mags = hcat(mags...)';
        dios = [s[i].diodes       for i = 1:N]; dios = reduce(hcat, dios)'; # dios = hcat(dios...)';
        gyrs = [s[i].gyro         for i = 1:N]; gyrs = reduce(hcat, gyrs)'; # gyrs = hcat(gyrs...)';
        poss = [s[i].pos          for i = 1:N]; poss = reduce(hcat, poss)'; # poss = hcat(poss...)';

        # Make the plots
        pM = plot(mags, title = "Magnetometers", xlabel = "Index", ylabel = "Mag Field (μT)",   label = ["x" "y" "z"]; kwargs...)
        pD = plot(dios, title = "Diodes",        xlabel = "Index", ylabel = "Current (A?)",                         ; kwargs...)
        pG = plot(gyrs, title = "Gyroscope",     xlabel = "Index", ylabel = "Ang Vel (rad/s)", label = ["x" "y" "z"]; kwargs...)
        pP = plot(poss, title = "Position",      xlabel = "Index", ylabel = "Position (m)",    label = ["x" "y" "z"]; kwargs...)
        return plot(pM, pD, pG, pP, plot_title = "Sensors", layout = (2, 2))
    
    elseif sensor == :m 
        mags = [s[i].magnetometer for i = 1:N]; mags = reduce(hcat, mags)';
        ms = []
        for i = 1:3 
            m = plot(mags[:, i])
            push!(ms, m)
        end

        return plot(ms..., plot_title = "Magnetometer Measurements", layout = 3)

    elseif sensor == :d
        dios = [s[i].diodes for i = 1:N]; dios = reduce(hcat, dios)'; 
        ds = []
        nd = size(s[1].diodes, 1)
        for i = 1:nd
            d = plot(dios[:, i])
            push!(ds, d)
        end
        return plot(ds..., plot_title = "Diode Measurements", layout = nd) #(3, 2))
        
    elseif sensor == :g 
        gyrs = [s[i].gyro         for i = 1:N]; gyrs = reduce(hcat, gyrs)';
        gs = []
        labels = ["x", "y", "z"]
        for i = 1:3 
            g = plot(gyrs[:, i], label = labels[i])
            push!(gs, g)
        end
        return plot(gs..., plot_title = "Gyro Measurements", layout = 3)

    elseif sensor == :p
        poss = [s[i].pos          for i = 1:N]; poss = reduce(hcat, poss)';
        ps = []
        labels = ["x", "y", "z"]
        for i = 1:3 
            p = plot(poss[:, i], label = labels[i])
            push!(ps, p)
        end
        return plot(ps..., plot_title = "Position Measurements", layout = 3)

    else 
        println("Warning! Symbol $sensor not supported yet!")
        println("\tAvailable symbols: a, m, d, g, p")
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
function RecipesBase.plot(gt::Vector{GROUND_TRUTH{6, T}}; ds = 1, split = false, kwargs...) where {T}
    N = size(gt, 1)
    ts  = [gt[i].t - gt[1].t for i = 1:ds:N];  # Convert to seconds 
    
    # Split apart and format as matrices for plotting
    Bᴵs = [gt[i].Bᴵ for i = 1:ds:N]; Bᴵs = hcat(Bᴵs...)';
    sᴵs = [gt[i].sᴵ for i = 1:ds:N]; sᴵs = hcat(sᴵs...)';
    ŝᴮs = [gt[i].ŝᴮ for i = 1:ds:N]; ŝᴮs = hcat(ŝᴮs...)';
    Bᴮs = [gt[i].Bᴮ for i = 1:ds:N]; Bᴮs = hcat(Bᴮs...)';
    Is  = [gt[i].I  for i = 1:ds:N]; Is  = hcat(Is... )';

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

        # q = (norm(q) == 1) ? q  : (q / norm(q))  # Enforce unit norm constraint
        # q = (q[1] < 0) ? -q : q
        new{T}(r, v, q, ω, β)
    end

    function STATE(r::Vector{T}, v::Vector{T}, q::Vector{T}, ω::Vector{T}, β::Vector{T}) where {T}
        # @warn "Initializing STATE without static arrays! Automatically converting..."
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
function RecipesBase.plot(s::Vector{STATE{T}}, state = :a, ds = 1, kwargs...) where {T}
    N = size(s, 1)

    if state == :a
        # Split apart and format as matrices for plotting
        rs = [vcat([s[i].r;]...) for i = 1:ds:N]; rs = hcat(rs...)';
        vs = [vcat([s[i].v;]...) for i = 1:ds:N]; vs = hcat(vs...)';
        qs = [vcat([s[i].q;]...) for i = 1:ds:N]; qs = hcat(qs...)';  
        ωs = [vcat([s[i].ω;]...) for i = 1:ds:N]; ωs = hcat(ωs...)';
        mag_ω = [norm(s[i].ω) for i = 1:ds:N]
        βs = [vcat([s[i].β;]...) for i = 1:ds:N]; βs = hcat(βs...)';

        pr = plot(rs, title = "Position (Cart)", xlabel = "Index", ylabel = "Position (m)",    label = ["x" "y" "z"]; kwargs...)
        pv = plot(vs, title = "Velocity",        xlabel = "Index", ylabel = "Velocity (m/s)",  label = ["x" "y" "z"]; kwargs...)
        pq = plot(qs, title = "Attitude (quat)", xlabel = "Index",                             label = ["x" "y" "z"]; kwargs...)
        pω = plot(ωs, title = "Ang Velocity",    xlabel = "Index", ylabel = "Ang Vel (rad/s)", label = ["x" "y" "z"]; kwargs...);
            pω = plot!(mag_ω, label = "||ω||", c = :black, ls = :dash)
        pβ = plot(βs, title = "Gyro Bias",       xlabel = "Index", ylabel = "Bias (rad/s)",    label = ["x" "y" "z"]; kwargs...)
        p   = plot()

        return plot(pr, pv, pq, pω, pβ, p, plot_title = "State", layout = (3, 2))
    elseif state == :r
        pos = [vcat([s[i].r;]...) for i = 1:ds:N]; pos = hcat(pos...)';
        rs = []
        labels = ["x" "y" "z"]
        for i = 1:3
            r = plot(pos[:, i], label = labels[i])
            push!(rs, r)
        end
        return plot(rs..., layout = 3, plot_title = "Position (Cart)", xlabel = "Index", ylabel = "Position (m)"; kwargs...)

    elseif state == :v  
        vels = [vcat([s[i].v;]...) for i = 1:ds:N]; vels = hcat(vels...)';
        vs = []
        labels = ["x" "y" "z"]
        for i = 1:3
            v = plot(vels[:, i], label = labels[i])
            push!(vs, v)
        end
        return plot(vs..., plot_title = "Velocity", xlabel = "Index", ylabel = "Velocity (m/s)"; kwargs...)

    elseif state == :q 
        quats = [vcat([s[i].q;]...) for i = 1:ds:N]; quats = hcat(quats...)';
        qs = []
        for i = 1:4
            q = plot(quats[:, i])
            push!(qs, q)
        end
        return plot(qs..., plot_title = "Attitude (quat)", xlabel = "Index", ylabel = "Magnitutde"; kwargs...)


    elseif state == :ω 
        angs = [vcat([s[i].ω;]...) for i = 1:ds:N]; angs = hcat(angs...)';
        ωs = []
        labels = ["x" "y" "z"]
        for i = 1:3
            ω = plot(angs[:, i], label = labels[i])
            push!(ωs, ω)
        end
        return plot(ωs..., plot_title = "Ang Vel", xlabel = "Index", ylabel = "Vel (rad/s)"; kwargs...)


    elseif state == :β
        bias = [vcat([s[i].β;]...) for i = 1:ds:N]; bias = hcat(bias...)';
        βs = []
        labels = ["x" "y" "z"]
        for i = 1:3
            β = plot(bias[:, i], label = "labels[i]")
            push!(βs, β)
        end
        return plot(βs..., plot_title = "Gyro Bias", xlabel = "Index", ylabel = "Magnitutde"; kwargs...)


    else 
        println("Warning: State $state is not ready for plotting!")
        println("\tViable symbols: r, v, q, ω, β, a")
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
    init_detumble::Bool 
    magnetometer_calibrated::Bool
    diodes_calibrated::Bool
    final_detumble::Bool 
    in_sun::Bool 

    function FLAGS(; in_sun = false, mag_cal = false, dio_cal = false, init_detumble = false, final_detumble = false)
        new(init_detumble, mag_cal, dio_cal, final_detumble, in_sun)
    end
end

end # Module
