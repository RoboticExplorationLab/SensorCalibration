module Simulator

using ..CustomStructs
using LinearAlgebra, SatelliteDynamics, EarthAlbedo, Distributions

using Infiltrator

include("../rotationFunctions.jl")
include("../mag_field.jl")

# Primary functions
export rk4
export generate_measurements

# Simulator options:
export SIM


####################################################################
#                         SOFTWARE                                 #
####################################################################

struct SIM # For now, just used for multiple dispatch 
    junk
end

function dynamics(sat::SATELLITE, x, u, t)
    """ 
        Propagates the state dynamics ([r⃗, v⃗, (q⃗, q₀), ω, β])
    """

    if norm(x[1:3]) < 6378136.3 # Satellite has crashed!
        error("Impact at time $t")
    end

    ṙ = @view x[4:6]
    v̇ = accel_perturbations(t, x[1:6]);

    w = @view x[11:13] # angular velocity 
    q̇ = 0.5 * qmult(x[7:10], [w; 0]) # Scalar LAST quaternion

    ẇ = (sat.J^(-1)) * (u - cross(w, (sat.J*w)))

    σβ = 0.1 * deg2rad(0.22)   # NOT SURE IF THIS IS ACTUALLY THE RIGHT NUMBER
    δβ = σβ * randn(3)  #β̇  looks funny so I am using δ

    return [ṙ[:]; v̇[:]; q̇[:]; ẇ[:]; δβ[:]]
end

function accel_perturbations(epc::Epoch, x::Array{<:Real} ;
    mass::Real=1.0, area_drag::Real=0.01, coef_drag::Real=2.3,
    area_srp::Real=1.0, coef_srp::Real=1.8,
    n_grav::Integer=10, m_grav::Integer=10)
    """Accelerations for spacecraft in LEO, ForwardDiff friendly (From Kevin)"""

    # Extract position and velocity
    r = x[1:3]
    v = x[4:6]

    # Compute ECI to ECEF Transformation -> IAU2010 Theory
    PN = bias_precession_nutation(epc)
    E  = earth_rotation(epc)
    W  = polar_motion(epc)
    R  = W * E * PN

    # Compute sun and moon position
    r_sun  = sun_position(epc)
    r_moon = moon_position(epc)

    # Compute acceleration (eltype(x) makes this forward diff friendly)
    a = zeros(eltype(x), 3)

    # spherical harmonic gravity
    a += accel_gravity(x, R, n_grav, m_grav)

    # atmospheric drag
    ρ = density_harris_priester(epc,r)
    a += accel_drag([r;v],ρ,mass, area_drag, coef_drag, Array{Real, 2}(PN))

    # SRP
    nu = eclipse_conical(x, r_sun)
    a += nu*accel_srp(x, r_sun, mass, area_srp, coef_srp)

    # third body sun
    a += accel_thirdbody_sun(x, r_sun)

    # third body moon
    a += accel_thirdbody_moon(x, r_moon)

    return a
end


function rk4(model, x, u, t, h)
    """ Modified rk4 function """

    k1 = h * dynamics(model, x, u, t)
    k2 = h * dynamics(model, x + k1/2, u, t + h/2)
    k3 = h * dynamics(model, x + k2/2, u, t + h/2)
    k4 = h * dynamics(model, x + k3,   u, t + h)

    return (x + (1/6)*(k1+2*k2+2*k3 + k4))
end

# Still relies on old faulty eclipse
function generate_measurements(sim::SIM, sat::SATELLITE, alb::ALBEDO, x, t, CONSTANTS, dt)
    """ 
        Generates sensor measurements, including noise.

        Arguments:
        - sim: Used to determine which simulator to use                                 | SIM 
        - sat: TRUE satellite data, which is used to generate sensor measurements       | SATELLITE
                    J (inertia matrix, [3 x 3])          //  magnetometer (calibration values, [3,]),
                    diodes (photodiode parameters [3,])  //  state (attitude and bias, [7,])
        - alb: Albedo struct containing REFL data and cell centers                      | ALBEDO
                    REFL    //   cell_centers_ecef
        - x: Environmental state (pos, vel, att, ang vel, bias)                         | [16,]
        - t: Current time, as an epoch                                                  | Scalar{Epoch}
        - CONSTANTS: Various constants needed (Earth radius, μ , E_am₀)                 | [3,]       
        - dt: Time step                                                                 | Scalar
        
        Returns:
        - truth: TRUTH struct containing true values for current time and vectors in body/inertial frame    | TRUTH
        - sensors: SENSORS struct containing noisy values for each satellite sensor being simulated         | SENSOR
        - ecl: Scale factor used to track how much the satellite is eclipsed by the Earth (∈ [0, 1])        | Scalar
        - noise: NOISE struct containing the amount of noise added to the sensors                           | NOISE
    """  
    pos = @view x[1:3]
    q⃗   = @view x[7:9]; q₀ = x[10] # Vector, Scalar portions

    ᴮRᴵ = dcm_from_q(x[7:10])'     # ᴺRᴮ = I(3) + 2 * hat(q⃗) * (q₀ * I(3) + hat(q⃗)); # Equation from Kevin, quaternion -> DCM (eq 39)

    sᴵₑ = sun_position(t)  # Sun-Earth vector
    ecl = eclipse_conical(-pos, sᴵₑ) # > 0.98 ? 1.0 : 0.0
    sᴵ = sᴵₑ - pos         # Sun-Sat vector

    albedo_matrix, junk = albedo(pos, sᴵ, alb.refl)  #  NOT GREAT with @code_warntype or @time

    sᴵ = ecl * (sᴵ) 
    Bᴵ = IGRF13(pos, t)

    η_sun = generate_noise_matrix(deg2rad(3.0), dt)  # Not really used...
    η_mag = generate_noise_matrix(deg2rad(3.0), dt)

    sᴮ = η_sun * (ᴮRᴵ * (sᴵ)) # (noisy) Sun vector in body frame     Also not really used...
    Bᴮ = η_mag * (ᴮRᴵ * (Bᴵ)) # (noisy) Mag vector in body frame

    𝐬ᴮ = ecl > 0.01 ? (η_sun * ecl * (ᴮRᴵ * (sᴵ / norm(sᴵ)))) :  [0.0; 0.0; 0.0]  # Avoid the divide-by-zero problem

    𝐁ᴮ = η_mag * (ᴮRᴵ * (Bᴵ / norm(Bᴵ)))   # unit

    num_diodes =  size(sat.diodes.calib_values, 1)
    current_vals = zeros(num_diodes) # Currents being generated by each photodiode
    current_noises = zeros(num_diodes)

    C, α, ϵ = sat.diodes.calib_values, sat.diodes.azi_angles, sat.diodes.elev_angles 

    for i = 1:num_diodes 
        surface_normal = [(cos(ϵ[i])*cos(α[i])) (cos(ϵ[i])*sin(α[i])) sin(ϵ[i])]   # Photodiode surface normal 
        diode_albedo = compute_diode_albedo(albedo_matrix, alb.cell_centers_ecef, surface_normal, pos)

        current = (C[i] * surface_normal * 𝐬ᴮ) .+ (C[i] * diode_albedo / CONSTANTS._E_am0) # Calculate current, including noise and Earth's albedo 
        current_noise = rand(Normal(0.0, abs((0.05 * current[1]))))
        current = current .+ current_noise # Generate noise that scales with magnitude of the current

        current_vals[i] = current[1] * ecl 
        current_noises[i] = current_noise
    end

    current_vals[current_vals .< 0] .= 0 # Photodiodes don't generate negative current

    mag_calib_matrix = generate_mag_calib_matrix(sat)

    Bᴮ = (mag_calib_matrix * Bᴮ) + sat.magnetometer.bias # Noise already added!

    w̃ = x[11:13] .+ x[14:16] 
    gyro_noise = 0.1 * deg2rad(1.0) * randn(3)  
    w̃ += gyro_noise

    pos = x[1:3] 
    gps_noise = (5e4) * randn(3) 
    pos += gps_noise

    Bᴮ_unadjusted = (ᴮRᴵ * (Bᴵ))

    sensors = SENSORS(Bᴮ, current_vals, w̃, pos)
    truth = GROUND_TRUTH(t, Bᴵ, sᴵ, 𝐬ᴮ, Bᴮ_unadjusted)
    noise = NOISE(current_noises,
                    gyro_noise, gps_noise, 
                    η_sun, η_mag)

    return truth, sensors, ecl, noise
end


function generate_noise_matrix(σ, dt)
    """
        Generates a [3 x 3] noise rotation matrix given a standard deviation 
            First generates a noise vector and then converts that into a rotation matrix
            (Note that if the standard deviation provided is 0, the identity matrix is returned)
    """
    if σ != 0.0
        η_vec = σ * randn(3) # Generate a vector 
        skew = hat(η_vec)
        norm_η = norm(η_vec)

        R = (I(3) + (skew/norm_η)*sin(norm_η*dt) + ((skew/norm_η)^2)*(1 - cos(norm_η*dt))); # Rodrigues for matrix exponential (?)
    else
        R = I(3) 
    end
    
    return R
end

# ALSO exists in ESTIMATOR...
function compute_diode_albedo(albedo_matrix, cell_centers_ecef, surface_normal, sat_pos)
    """ 
        Estimates the effect of Earth's albedo on a specific photodiode (by using the surface normal of that diode)
            = cell_albedo * surface_normal^T * r_g, with r_g as a unit vector in the direction of the grid point on Earth

        Arguments:
        - albedo_matrix: Albedo values for each cell on the Earth's surface         | [num_lat x num_lon] 
        - surface_normal: Photodiode surface normal                                 | [3,]
        - sat_pos: Cartesian position of satellite                                  | [3,]

        Returns:
        - diode_albedo: Total effect of albedo on specified photodiode              | Scalar
    """    
    cell_albedos = zeros(size(albedo_matrix))

    rows, cols = size(albedo_matrix)
    diode_albedo = 0.0
    for r = 1:1:rows
        for c = 1:1:cols
            if albedo_matrix[r,c] != 0
                r_g = cell_centers_ecef[r,c,:] - sat_pos # Distance from satellite to cell center
                r_g = r_g / norm(r_g)  # Make unit

                cell_albedo = (albedo_matrix[r,c] * (surface_normal * r_g))[1]

                if cell_albedo > 0.0    # Can't be negative
                    diode_albedo = diode_albedo + cell_albedo 
                end
            end
        end
    end
    
    return diode_albedo
end

function generate_mag_calib_matrix(sat::SATELLITE)
    """ Generates the calibration matrix that alters the measured magnetic field vector in body frame """
    a, b, c = sat.magnetometer.scale_factors
    ρ, λ, ϕ = sat.magnetometer.non_ortho_angles

    T = [a        0.0              0.0;
        b*sin(ρ)  b*cos(ρ)         0.0;
        c*sin(λ)  c*sin(ϕ)*cos(λ)  c*cos(ϕ)*cos(λ)]

    return T
end


####################################################################
#                         HARDWARE                                 #
####################################################################
# function propagate_dynamics()
# end

# function generate_measurements()
# end


end