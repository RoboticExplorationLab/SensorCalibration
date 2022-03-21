""" Masurements.jl:

    Contains all relevant functions to generate state measurements, including the
        simulation of noise for sensors.

    TODO:
    - (When eclipse_conical is updated) Adjust -pos to pos
    - Adjust noise to not include η_sun/_mag? Or do include them as vectors?
    - Once EarthAlbedo is updated, use its get_diode_albedo
    (Speed up diode/Earth Albedo stuff)

"""

# Update noise
function generate_measurements(sat::SATELLITE, alb::ALBEDO, x, t, CONSTANTS, dt)
    """ 
        Generates sensor measurements, including noise.

        Arguments:
        - sat: TRUE satellite data, which is used to generate sensor measurements       | SATELLITE
                    J (inertia matrix, [3 x 3])          //  magnetometer (calibration values, [3,]),
                    diodes (photodiode parameters [3,])  //  state (attitude and bias, [7,])
        - alb: Albedo struct containing REFL data and cell centers                      | ALBEDO
                    REFL    //   cell_centers_ecef
        - x: Environmental state (pos, vel, att, ang vel, bias)                         | [16,]
        - t: Current time, as an epoch                                                  | Epoch
        - CONSTANTS: Various constants needed (Earth radius, μ , E_am₀)                 | [3,]       
        - dt: Time step                                                                 | Scalar
        
        Returns:
        - truth: TRUTH struct containing true values for current time and vectors in body/inertial frame    | TRUTH
        - sensors: SENSORS struct containing noisy values for each satellite sensor being simulated         | SENSOR
        - ecl: Scale factor used to track how much the satellite is eclipsed by the Earth (∈ [0, 1])        | Scalar
        - noise: NOISE struct containing the amount of noise added to the sensors                           | NOISE
    """  
    ᴮRᴵ = dcm_from_q(view(x, 7:10))'  

    # Bᴵ, Bᴮ, B̃ᴮ = generate_magnetic_field(view(x, 1:3), t, sat, ᴮRᴵ, dt)
    Bᴵ, Bᴮ, B̃ᴮ = generate_magnetic_field(view(x, 1:3), t, sat, ᴮRᴵ, dt, view(x, 17:19))
    w, w̃, ηω   = generate_gyro_measurement(x) 
    r, r̃, ηr   = generate_gps_measurement( x) 
    sᴵ, 𝐬ᴮ, ecl = update_sun_vectors(view(x, 1:3), t, ᴮRᴵ, dt)

    I, Ĩ, ηI = generate_diode_currents(sat, view(x, 1:3), alb, sᴵ, 𝐬ᴮ, ecl, CONSTANTS)

    sensors = SENSORS(B̃ᴮ, Ĩ, w̃, r̃)
    truth = GROUND_TRUTH(t, Bᴵ, sᴵ, 𝐬ᴮ, Bᴮ)
    junk_noise = zeros(3,3)
    noise = NOISE(ηI, ηω, ηr, junk_noise, junk_noise)

    return truth, sensors, ecl, noise
end

function generate_magnetic_field(pos, t, sat, ᴮRᴵ, dt, mag_bias)
 # function generate_magnetic_field(pos, t, sat, ᴮRᴵ, dt)
    """
        Generates the true magnetic field vector in inertial and body frames, 
            as well as the noisy sensor measurement.

        Arguments:
        - pos: Current (x, y, z) position of the satellite in ECEF                                          | [3,]
        - t: Current time, as an epoch                                                                      | Epoch
        - sat: TRUE satellite data, which is used to generate sensor measurements                           | SATELLITE
                    J (inertia matrix, [3 x 3])          //  magnetometer (calibration values, [3,]),
                    diodes (photodiode parameters [3,])  //  state (attitude and bias, [7,])
        - ᴮRᴵ: Rotation matrix from inertial frame to body frame                                            | [3 × 3]
        - dt: Time step                                                                                     | Scalar
        
        Returns:
        - Bᴵ: Magnetic field vector in inertial frame (noiseless)                                           | [3,]
        - Bᴮ: Magnetic field vector in body frame (noiseless)                                               | [3,]
        - B̃ᴮ: Measured magnetic field vector in body frame (noisy, biased, and scaled)                      | [3,]
    """

    Bᴵ = IGRF13(pos, t)     # Mag vector in inertial frame
    η_mag = generate_noise_matrix(σ_mag, dt)
    Bᴮ = (ᴮRᴵ * (Bᴵ))       # Mag vector in body frame

    mag_calib_matrix = generate_mag_calib_matrix(sat)
    # B̃ᴮ = (mag_calib_matrix * η_mag * Bᴮ) + sat.magnetometer.bias   # Measured mag vector in body frame
    B̃ᴮ = (mag_calib_matrix * η_mag * Bᴮ) + mag_bias

    return Bᴵ, Bᴮ, B̃ᴮ
end

function generate_gyro_measurement(state)
    """ 
        Generates the gyroscope measurement by adding noise into the true values.

        Arguments: 
        - state: Environmental state (pos, vel, att, ang vel, bias)                          | [16,]
        
        Returns:
        - w: Angular velocity                                                            | [3,]
            (NOTE that this is just returned for convenience, as it is already in state)
        - w̃: Measured angular velocity (including noise and bias)                        | [3,]
        - ηw: Gyroscope noise                                                            | [3,]
            (Note that this is really just tracked for debugging purposes)
    """
    w = @view state[11:13] 
    β = @view state[14:16]
    ηw = rand(Normal(μ_gyro_scale * norm(w), σ_gyro_scale * norm(w)), 3)  # Noise scales with actual value
    
    w̃ = w .+ β .+ ηw 

    return w, w̃, ηw
end

function generate_gps_measurement(state)
    """ 
        Generates the GPS measurement by adding noise into the true values.

        Arguments: 
        - state: Environmental state (pos, vel, att, ang vel, bias)                          | [16,]
        
        Returns:
        - r: Global position                                                             | [3,]
            (NOTE that this is just returned for convenience, as it is already in state)
        - r̃: Measured position (including noise)                                         | [3,]
        - ηr: Position noise                                                             | [3,]
            (Note that this is really just tracked for debugging purposes)
    """
    r = @view state[1:3]
    ηr = rand(Normal(μ_gps, σ_gps), 3)
    r̃ = r + ηr
    
    return r, r̃, ηr
end

# (Correct eclipse_conical)
function update_sun_vectors(pos, t, ᴮRᴵ, dt)
    """
        Generates the sun vector in the body and inertial frames, as well as whether or not the 
            sun is being eclipsed by the Earth (from the satellites perspective).
            
        Arguments:
        - pos: [x, y, z] position of the satellite in ECEF frame                                     | [3,]
        - t: Current time, as an epoch                                                               | Epoch
        - ᴮRᴵ: Rotation matrix from inertial frame to body frame                                     | [3 × 3]
        - dt: Time step                                                                              | Scalar
                                      
        Returns:
        - sᴵ: Satellite-sun vector in the inertial frame                                             | [3,]
        - 𝐬ᴮ: Unit satellite-sun vector in the body frame                                            | [3,]
        - ecl: Eclipse factor                                                                        | Scalar
        
    """
    sᴵₑ = sun_position(t) 
    sᴵ = sᴵₑ - pos 

    ecl = eclipse_conical(-pos, sᴵₑ)
    sᴵ *= ecl 

    η_sun = generate_noise_matrix(σ_sun, dt)   #I(3)
    𝐬ᴮ = ecl > 0.01 ? (η_sun * (ᴮRᴵ * (sᴵ / norm(sᴵ)))) :  SVector(0.0, 0.0, 0.0)  # Avoid the divide-by-zero problem

    return sᴵ, 𝐬ᴮ, ecl
end

function generate_diode_currents(sat, pos, alb, sᴵ, 𝐬ᴮ, ecl, CONSTANTS)
    """ 
        Generates photodiode currents and current measurements. 
            (NOTE that current generated cannot be negative so these are clipped at zero)

        Arguments:
        - sat: TRUE satellite data, which is used to generate sensor measurements                       | SATELLITE
                    J (inertia matrix, [3 x 3])          //  magnetometer (calibration values, [3,]),
                    diodes (photodiode parameters [3,])  //  state (attitude and bias, [7,])
        - pos: Current (x, y, z) position of the satellite in ECEF                                      | [3,]
        - alb: Albedo struct containing REFL data and cell centers                                      | ALBEDO
                    REFL    //   cell_centers_ecef
        sᴵ: Satellite-Sun vector in the inertial frame (noiseless)                                      | [3,]
        𝐬ᴮ: Unit satellite-sun vector in the body frame (noiseless)                                     | [3,]
        - CONSTANTS: Various constants needed (Earth radius, μ , E_am₀)                                 | [3,]       
        
        Returns:
        - current_vals: Currents generated by each photodiode                                           | [num_diodes,]
        - current_meas: Measured (noisy) currents for each photodiode                                   | [num_diodes,]
        - current_noises: Current noise added to the measurement                                        | [num_diodes,]
            (Note that this is really just tracked for debugging purposes)
    """  
    C, α, ϵ = sat.diodes.calib_values, sat.diodes.azi_angles, sat.diodes.elev_angles 

    num_diodes =  size(sat.diodes.calib_values, 1)
    current_vals = zeros(num_diodes) # Currents being generated by each photodiode
    current_noises = zeros(num_diodes)
    current_meas = zeros(num_diodes)

    albedo_matrix, _ = albedo(pos, sᴵ, alb.refl) 

    for i = 1:num_diodes 
        surface_normal = [(cos(ϵ[i])*cos(α[i])) (cos(ϵ[i])*sin(α[i])) sin(ϵ[i])]   # Photodiode surface normal 

        diode_albedo = compute_diode_albedo(albedo_matrix, alb.cell_centers_ecef, surface_normal, pos)

        current = ((C[i] * surface_normal * 𝐬ᴮ) .+ (C[i] * diode_albedo / CONSTANTS._E_am0))[1] # Calculate current, including noise and Earth's albedo 
        current_noise = rand(Normal(0.0, abs((σ_current_scale * current))))

        current_meas[i] = (current .+ current_noise) * ecl
        current_vals[i] = current * ecl 
        current_noises[i] = current_noise
    end

    current_meas[current_meas .< 0] .= 0
    current_vals[current_vals .< 0] .= 0 # Photodiodes don't generate negative current
    
    return current_vals, current_meas, current_noises
end 



####################
# HELPER FUNCTIONS #
####################

function generate_noise_matrix(σ, dt)
    """
        Generates a [3 x 3] noise rotation matrix given a standard deviation 
            First generates a noise vector and then converts that into a rotation matrix
            (Note that if the standard deviation provided is 0, the identity matrix is returned)
    """
    if σ != 0.0
        η_vec = rand(Normal(0.0, σ), 3)  # Generate a vector 
        skew = hat(η_vec)
        norm_η = norm(η_vec)

        R = (I(3) + (skew/norm_η)*sin(norm_η*dt) + ((skew/norm_η)^2)*(1 - cos(norm_η*dt))); # Rodrigues for matrix exponential (?)
    else
        R = I(3) 
    end
    
    return R
end


#TODO remove this once EarthAlbedo is updated
# ALSO exists in ESTIMATOR, and Earth_albedo (?)... -> make a bit more general and stick with EarthAlbedo's
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
    diode_albedo = 0.0
    r_g = zeros(Float64, 3)
    for r = 1:1:size(albedo_matrix, 1)
        for c = 1:1:size(albedo_matrix, 2)
            if albedo_matrix[r,c] != 0
                r_g .= view(cell_centers_ecef, r, c, :) .- sat_pos
                r_g .= r_g ./ norm(r_g)  # Make unit

                cell_albedo = (albedo_matrix[r,c] * dot(surface_normal, r_g))

                if cell_albedo > 0.0    # Can't be negative
                    diode_albedo += cell_albedo 
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


