# System parameters 
using LinearAlgebra, SatelliteDynamics, Distributions, EarthAlbedo

# STRUCTS ----------------------------------------------------------------------#


# CONSTANTS --------------------------------------------------------------------#
    _Re = 6378136.3                 # Radius of the earth (m)) 
    _mu = (3.9860044188)*(10^(14))  # Standard Gravitational Parameter (m^3/s^2) 

    sat_mass = 1;  # (kg)
    sat_len = 0.1; # (m)

    _J = ((1/12) * (sat_mass) * ((sat_len^2) + (sat_len^2))) * I(3); # Inertia Matrix (units...?)
    _J = _J .* [0.85 1.0 1.18] #rand(Uniform(0.75, 1.25), 3)
    println("Inducing a non-uniform inertia matrix")

    _E_am0 = 1366.9 # Irradiance of sunlight (TSI - visible & infrared), W/m^2 

    CONSTANTS = (_Re = _Re, _mu = _mu, _E_am0 = _E_am0)

# STATE ------------------------------------------------------------------------#
    # Initial position and velocity (in terms of osculating orbital elements)
    sma = _Re + 550e3 + 5000 * randn()  # Semi-major axis       (m)
    ecc = 0.01 + 0.001 * randn()       # Eccentricity          ()
    inc = 96.0 + 5.0 * randn()         # Inclination           (deg)
    Ω   = 75.0 + 5.0 * randn()         # RAAN                  (deg)
    ω   = 30.0 + 3.0 * randn()         # Argument of perigee   (deg)
    M   = 0.0  + 2.0 * randn()         # Mean anomaly          (deg)
    
    oe0 = [sma, ecc, inc, Ω, ω, M]     # Initial state 
    eci0 = sOSCtoCART(oe0, use_degrees = true) # Convert to Cartesean

    # Initial attitude and angular velocity 
    r = randn(3)
    r = r / norm(r);       
    θ = rand(Uniform(0, pi))                         
    q0 =  [r * sin(θ/2); cos(θ/2)];
    q0 = q0 / norm(q0)
    w0 = 1.0 .* deg2rad.([10, 2, -7.2]) # HOW DO I SELECT REASONABLE VALUES HERE...?    
    # w0 = 0.05 .* w0

    # Initial Bias
    β0 = deg2rad(1.0) * randn(3)    

    x0 = [eci0[:]; q0[:]; w0[:]; β0[:]]   # Initial state   | [16,]


# DIODES -----------------------------------------------------------------------#
    _num_diodes = 6
    _sensor_scale_factors = rand(Normal(1.0, 0.15), _num_diodes)  # Calibration values for each photodiode

    # Installation angles for each photodiode
    _elev_angles = [0.0; 0.0;  0.0;    0.0;    (pi/4); (-pi/4)] #.+ (pi/6)
    _azi_angles  = [0.0; pi;   (pi/2); (-pi/2); 0.0;    (-pi/4)]   #.+ (pi/6)

    println("Diodes rotated 45° around ŷ")

    _elev_angles += rand(Normal(0.0, deg2rad(2.0)), _num_diodes)
    _azi_angles  += rand(Normal(0.0, deg2rad(2.0)), _num_diodes)

    # _elev_angles = rand(Normal( (pi/4), 0.4), _num_diodes) 
    # _elev_angles[_elev_angles .> (pi/2)] .= (pi/2);  #   Elevation ∈ [0, π/2]
    # _elev_angles[_elev_angles .< 0] .= 0

    # _azi_angles = rand(Normal( (pi), 1.0), _num_diodes)
    # _azi_angles[_azi_angles .>= (2*pi)] .= (2*pi - 0.1); # Azimuth ∈ [0, 2π)
    # _azi_angles[_azi_angles .< 0] .= 0


# MAGNETOMETER -----------------------------------------------------------------#
    _a, _b, _c = rand(Normal(1.0, 0.1), 3) # Linear scale factors 
    _βx₀, _βy₀, _βz₀ = rand(Normal(0.0, 2), 3) # Bias (μTeslas) (?)
    # println("bias should be 5%, its currently 0 in init and dynamics")
    _ρ, _ϕ, _λ = rand(Normal(0.0, deg2rad(3.0)), 3) # 0.0 * randn(3)        # Non-orthogonality non_ortho_angles



# GENERAL ----------------------------------------------------------------------#
    _dt = 1.0 #0.2 # 1.0  # (s)    
    # println("Running at 5 Hz!")
    _T = round(5.0 * orbit_period(oe0[1]) / _dt)  # Run for 3 orbits    
    _epc = Epoch(2021, 9, 1, 12, 0, 0, 0.0); # Initial time for sim
    _max_sim_length = Int(_T)

    SYSTEM = (_dt = _dt, _T = _T, _epc = _epc, _max_sim_length = _max_sim_length, _num_diodes = _num_diodes)

println("No real noise yet")