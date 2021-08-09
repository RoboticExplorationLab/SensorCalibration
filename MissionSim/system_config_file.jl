using LinearAlgebra, SatelliteDynamics, Distributions, EarthAlbedo

# CONSTANTS --------------------------------------------------------------------#
    const _Re = 6378136.3                 # Radius of the earth (m)) 
    const _mu = (3.9860044188)*(10^(14))  # Standard Gravitational Parameter (m^3/s^2) 

    const sat_mass = 1;  # (kg)
    const sat_len = 0.1; # (m)

    const J_temp = ((1/12) * (sat_mass) * ((sat_len^2) + (sat_len^2))) * I(3); # Inertia Matrix (units...?)
    const _J = J_temp .* [0.86 1.0 1.2] 

    const _E_am0 = 1366.9 # Irradiance of sunlight (TSI - visible & infrared), W/m^2 

    const CONSTANTS = (_Re = _Re, _mu = _mu, _E_am0 = _E_am0)

# STATE ------------------------------------------------------------------------#

    # # Initial position and velocity (in terms of osculating orbital elements)
    # orbit = "Generic"
    # sma = _Re + 550e3 + 5000 * randn()  # Semi-major axis       (m)
    # ecc = 0.01 + 0.001 * randn()       # Eccentricity          ()
    # inc = 96.0 + 5.0 * randn()         # Inclination           (deg)
    # Ω   = 90 + 75.0 + 5.0 * randn()         # RAAN                  (deg)
    # ω   = 30.0 + 5.0 * randn()         # Argument of perigee   (deg)
    # M   = 0.0  + 5.0 * randn()         # Mean anomaly          (deg)
    

    # """   ISS   """
    const orbit = "ISS"
    const ecc = 0.0001717 + 0.000001 * randn()
    const inc = 51.6426 + randn()
    const Ω   = 178.1369 + randn()
    const ω   = 174.7410 + randn()
    const M   = 330.7918 + 180  + 110 #+ randn() # +94/95 is just before sun
    const sma = (_Re + 421e3) / (1 + ecc)  # Apogee = semi_major * (1 + ecc)


    const oe0 = [sma, ecc, inc, Ω, ω, M]     # Initial state 
    const eci0 = sOSCtoCART(oe0, use_degrees = true) # Convert to Cartesean

    # Initial attitude and angular velocity 
    const r_temp = randn(3)
    const r = r_temp / norm(r_temp);       
    const θ = rand(Uniform(0, pi))                          
    const q0_temp =  [r * sin(θ/2); cos(θ/2)];
    const q0 = q0_temp / norm(q0_temp)
    const w0 = deg2rad.(rand(Normal(0.0, 5.0),3))
    # w0 = 1.25 .* deg2rad.([10, 2, -7.2]) # HOW DO I SELECT REASONABLE VALUES HERE...?    

    # Initial Bias
    const β0 = deg2rad(2.0) * randn(3)    

    const x0 = [eci0[:]; q0[:]; w0[:]; β0[:]]   # Initial state   | [16,]


# DIODES -----------------------------------------------------------------------#
    const _num_diodes = 6
    const _sensor_scale_factors = rand(Normal(1.0, 0.15), _num_diodes)  # Calibration values for each photodiode

    # Installation angles for each photodiode
    const elev_angles = [(-pi/4); (pi/4);  0.0;    0.0;    (pi/4); (-pi/4)] 
    const azi_angles  = [0.0;      pi;     (pi/2); (-pi/2); 0.0;    pi]  
    # _elev_angles = [0.0;     0.0;        0.0;    0.0;  (pi/2); (-pi/2)] 
    # _azi_angles  = [0.0;      pi;     (pi/2); (-pi/2); 0.0;    0.0]  

    const _elev_angles = elev_angles + rand(Normal(0.0, deg2rad(2.0)), _num_diodes)
    const _azi_angles  = azi_angles + rand(Normal(0.0, deg2rad(2.0)), _num_diodes)


# MAGNETOMETER -----------------------------------------------------------------#
    const _a, _b, _c = rand(Normal(1.0, 0.1), 3)           # Linear scale factors 
    const _βx₀, _βy₀, _βz₀ = rand(Normal(0.0, 2), 3)       # Bias (μTeslas) (?)  ################################
    const _ρ, _ϕ, _λ = rand(Normal(0.0, deg2rad(3.0)), 3)  # Non-orthogonality non_ortho_angles

# GENERAL ----------------------------------------------------------------------#
    const _dt = 0.2 #1.0  # (s)    
    const _run_freq = 1 / _dt
    const _T = round(1.25 * orbit_period(oe0[1]) / _dt)  # Run for 2.25 orbits
    const _epc = Epoch(2021, 9, 1, 11, 0, 0, 0.0); # Initial time for sim  ############## add in randomness to time?
    const _max_sim_length = Int(_T)
    const _ds_rate = Int(round(120/_dt))

    const SYSTEM = (_dt = _dt, _T = _T, _epc = _epc, _max_sim_length = _max_sim_length, _num_diodes = _num_diodes)
    
    const eclipse_threshold = 0.02

@info "Generating $orbit orbit and running at $_run_freq Hz"