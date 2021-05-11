using Random, Distributions # Used to generate random numbers for noise, photodiode values

# OTHER ------------------------------------------------------------------------
dscale = 1e6
tscale = 3600/5 
num_states = 16
num_meas = 19

# CONSTANTS --------------------------------------------------------------------
_Re = 6378136.3 / dscale;            # Radius of the earth (adjusted units) 
_mu = (tscale^2)*(3.9860044188)*(10^(14)) / (dscale ^ 3);  # (m^3/s^2) Standard Gravitational Parameter
_r0 = (550000 / dscale) + _Re;       # Distance between Earth and Satellite centers of masses
_a = _r0;                            # Semi-major axis of (elliptical) orbit (circular orbit if a = r0)
_v0 = sqrt(_mu*( (2/_r0) - (1/_a))); # Velocity for a stable orbit
sat_mass = 1;  # kg
sat_len = 0.1; # m
_J = ((1/12) * (sat_mass) * ((sat_len^2) + (sat_len^2)) * I(3)) / (dscale^2); # Should this be divided by dscale...? It *is* in meters
dynamics_params = [_mu, _Re, _r0, _J]; # Used to propagate dynamics

# SIMULATOR -------------------------------------------------------------------
epc = Epoch(2019, 1, 1, 12, 0, 0, 0.0); # Initial time for sim
dt = 1.0 / tscale;
sim_length = (120.0 * 60.0) / tscale;   # Seconds (adjusted)
T = Int(sim_length / dt)                # Knot points (?)

# INITIAL CONDITIONS ----------------------------------------------------------
r = [1; 1; 1/2]; r = r / norm(r);       
θ = pi/4;                               
q0 =  [r * sin(θ/2); cos(θ/2)];         # Initial satellite orientation
p0 = [0, _r0, 0];                       # Initial satellite position, m (adjusted)
v0 = [_v0, 0, 0];                       # Initial satellite velocity, m/s (adjusted)
w0 = [0.012, 0.02, 0.031] * tscale;     # Initial satellite angular velocity rad/sec (adjusted)
β0 = deg2rad(2.0) * randn(3)            # Initial gyroscope bias, rad

x0 = [p0; v0; q0; w0; β0];              # Initial state   | [16,]

# NOISE ------------------------------------------------------------------------
σ_β = deg2rad(0.22) * tscale  # Gyro noise (adjusted to match w)

σ_η_sun = deg2rad(5.0)        # Sensor noise for sun vector (radians)
σ_offset_sun = deg2rad(0.0)   # Initial offset for sun vector (radians)

σ_η_mag = deg2rad(5.0)        # Sensor noise for mag vector (radians)
σ_offset_mag = deg2rad(0.0)   # Initial offset for mag vector (radians)

σ_η_cur = 0.008               # Noise in photodiode current (A)

# PHOTODIODES ------------------------------------------------------------------
num_diodes = 6
calib_vals = rand(Normal(1, 0.2), num_diodes)  # Calibration values for each photodiode

ϵs = rand(Normal((pi/4), 0.3), num_diodes) # Elevation angles for each photodiode
ϵs[ϵs .> (pi/2)] .= (pi/2);  #   Elevation ∈ [0, π/2]
ϵs[ϵs .< 0] .= 0

αs = rand(Normal((pi), 1), num_diodes)       # Azimuth angles for each photodiode
αs[αs .>= (2*pi)] .= (2*pi - 0.1); # Azimuth ∈ [0, 2π)
αs[αs .< 0] .= 0