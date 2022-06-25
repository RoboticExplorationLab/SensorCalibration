# [src/MissionSim/Simulator/measurements.jl]

"""
      Contains the code used in the Simulator to generate 
    the measurements for each of the sensors.

    To Do:
      Add:

      Fix/Update:
      - These noise values are almost entirely arbitrary; get better ones
      - diode_meas is reaaaaal slow (2.8ms / 19 / 413 KiB)  -> earth_albedo takes up 2.241 ms of that, but if we downsample albedo before...
      - IGRF in mag_meas is also slow ∼ μs
      - (Validate - σ_scales used when additive noise, σ for rotation noise, ya?)
      - (Verify - no need for μ in gyro, pos? No need for noise in sun?)
      - Should get_calib_matrix be static?
"""

# σ_current, gyro  no longer scales, but matches the values in On-Orbit. I do noise on mag with rotation, not noise, but I ran some 
#   tests to determine what the angular distance corresponded to 0.25deg
"""
    generate_measurements(sat, alb, x, t, dt; E_am₀, σB, σ_gyro_scale, σr, σ_current_scale)

      Uses the current environment state and satellite parameters to generate both noisy measurements and what the 
    noise-free measurements would be. Sensors include magnetometers, gyroscopes, and photodiodes. Additionally, the 
    sun vector is estimated, and position data is provided (based off of true position in STATE x).
    
      Optional arguments provide the ability to control the noise values for each sensor. 

    Arguments:
      - `sat`:  (True) SATELLITE struct, which is used to generate sensor measurements     |  SATELLITE
                    (J  //  magnetometer  //  diodes)
      - `alb`:  struct containing REFL data and cell centers                               |  ALBEDO
                    (REFL  //  cell_centers_ecef)
      - `x`:    struct containing current environment state                                |  STATE 
                    (r  //  v  //  q  //  ω  //  β)
      - `t`:    Current time, as an epoch                                                  |  Epoch
      - `dt`:   Time step                                                                  |  Scalar 

      - `E_am₀`:  (Opt) Irradiance of sunlight (TSI - visible & infrared). Default is 1366.0 W/m²     |  Scalar
      - `σB`:     (Opt) stdev of noise used when generating the mag vector in body frame, in radians  |  Scalar
      - `σ_gyro_scale`: (Opt) Scale that modifies ||gyro|| to get the stdev of the gyro noise         |  Scalar
      - `σr`:     (Opt) stdev of noise used when generating the position measurement, in m            |  Scalar 
      - `σ_current_scale`: (Opt) Scale that modifies ||diode current|| to get stdev of diode noise    |  Scalar 

    Returns:
      - `truth`:    Struct containing noise-free values, used for evaluation and debugging   |  GROUND_TRUTH 
      - `sensors`:  Struct containing noisy sensor "measurements" to be used in the MEKF     |  SENSORS
      - `ecl`:      Eclipse factor (1.0 = full sun, 0.0 = full eclipse)                      |  Scalar
      - `noise`:    Struct containing the noise used in the sensor. Used for debugging       |  NOISE

"""
function generate_measurements(sat::SATELLITE, alb::ALBEDO, x::STATE, t::Epoch, dt::T; 
             E_am₀ = 1366.9, σB = deg2rad(0.25), σ_gyro = 0.5e-4, 
             σr = 5e3, σ_current = 0.05, use_albedo = true) where {T}

    ᴮQᴵ = quat2rot(x.q)'  # Compute once and pass in 

    # Call each sensor
    Bᴵ, Bᴮ, B̃ᴮ  = mag_measurement( sat, x, ᴮQᴵ, t, dt; σ = σB) 
    w, w̃, ηw    = gyro_measurement(x; σ_scale = σ_gyro) 
    r, r̃, ηr    = pos_measurement( x; σ = σr) 
    sᴵ, sᴮ, ecl = sun_measurement( x, ᴮQᴵ, t)
    I, Ĩ, ηI    = diode_measurement(sat, alb, x, ecl, sᴵ, sᴮ; σ_scale = σ_current, E_am₀ = E_am₀, use_albedo = use_albedo)

    # Store the outputs 
    sensors = SENSORS(B̃ᴮ, Ĩ, w̃, r̃)
    truth   = GROUND_TRUTH(t, Bᴵ, sᴵ, sᴮ, Bᴮ, I)  # w, r are stored in state
    noise   = NOISE(ηI, ηw, ηr)

    return truth, sensors, ecl, noise
end



# I have noise as a rotation, the paper has it as additive with σ = 1e-9 T, which seems to only do 0.25 deg of rotation
"""
      Generates the true magnetic field vector in inertial and body frames. Additionally, 
    uses the satellite's magnetometer parameters, the states bias, and a noise rotation 
    matrix to generate the measured body vector.

    Arguments:
      - `sat`: (True) SATELLITE struct, which is used to generate sensor measurements     |  SATELLITE
                      (J  //  magnetometer  //  diodes)
      - `x`:   struct containing current environment state                                |  STATE 
                      (r  //  v  //  q  //  ω  //  β)
      - `ᴮQᴵ`: Rotation matrix from the inertial to the body frame                        |  [3, 3]  (Static)
      - `t`:   Current time, as an epoch                                                  |  Epoch
      - `dt`:  Time step                                                                  |  Scalar 
      - `σ`:   (Opt) stdev of noise used when generating the measured mag vector          |  Scalar
                in body frame, in radians ('σB' in generate_measurements)
 
    Returns:
      - `Bᴵ`: Magnetic field vector in inertial frame (noiseless)                               | [3,] (SVector)
      - `Bᴮ`: Magnetic field vector in body frame (noiseless)                                   | [3,] (SVector)
      - `B̃ᴮ`: Measured magnetic field vector in body frame (noisy, biased, and un-calibrated)   | [3,] (SVector)
"""
function mag_measurement(sat::SATELLITE, x::STATE, ᴮQᴵ::SMatrix{3, 3, T, 9}, t::Epoch, dt::T; σ = deg2rad(0.25)) where {T}

    Bᴵ = SVector{3, Float64}(IGRF13(x.r, t))   # Mag vector in inertial frame 
    Bᴮ = ᴮQᴵ * Bᴵ 

    # Add in noise
    η_mag = rotation_noise(σ, dt)
    mag_calibration_matrix = get_mag_calibration_matrix(sat)
    B̃ᴮ = (η_mag * mag_calibration_matrix * Bᴮ) + sat.magnetometer.bias 

    return Bᴵ, Bᴮ, B̃ᴮ
end

# no longer have a μ_scale
""" 
    Generates the gyroscope measurement by adding noise to the true values.

    Arguments: 
      - x:  struct containing current environment state                             |  STATE 
              (r  //  v  //  q  //  ω  //  β)
      - σ_scale: (Opt) Scale that modifies magnitude of the gyroscope to get the 
              to get the stdev of the gyro noise, so that noise magnitude scales with 
              gyroscope magnitude.  ('σ_gyro_scale' in generate_measurements)       |  Scalar

    Returns:
      - w: Angular velocity                                                            | [3,] (Static)
          (NOTE that this is just returned for consistency, as it is already in state)
      - w̃: Measured angular velocity (including noise and bias)                        | [3,] (Static)
      - ηw: Gyroscope noise                                                            | [3,] (Static)
          (Note that this is really just tracked for debugging purposes)
"""
function gyro_measurement(x::STATE{T}; σ_scale = 0.005) where {T}

    ηω = rand(Normal(0.0, σ_scale * norm(x.ω)), 3)  # Noise scales with actual value

    ω̃  = x.ω .+ x.β .+ ηω

    return x.ω, SVector{3, T}(ω̃ ), SVector{3, T}(ηω)
end

# (Note that r is tracked in state, but this allows for some adjustments); no longer have a μ
""" 
    Generates the position measurement by adding noise to the true values in state

    Arguments: 
      - x:  struct containing current environment state                           |  STATE 
                (r  //  v  //  q  //  ω  //  β)
      - σ:  (Opt) stdev of noise used when generating the position measurement    |  Scalar 
                ('σr' in generate_measurements)

    Returns:
      - r: Global position                                                             | [3,] (Static)
          (NOTE that this is just returned for consistency, as it is already in state)
      - r̃: Measured position (including noise)                                         | [3,] (Static)
      - ηr: Position noise                                                             | [3,] (Static)
          (Note that this is really just tracked for debugging purposes)
"""
function pos_measurement(x::STATE; σ = 1e4)

    ηr = SVector{3, Float64}(rand(Normal(0.0, σ), 3))
    r̃ = x.r + ηr

    return x.r, r̃, ηr
end

# NOTE i am not scaling sᴵ, sᴮ by ecl anymore, nor am I adding noise to sᴮ
"""
    Generates the sun vector in the body and inertial frames, as well as whether or not the 
    sun is being eclipsed by the Earth (from the satellite's perspective).
      
    Arguments:
      - x:   struct containing current environment state                     |  STATE 
                   (r  //  v  //  q  //  ω  //  β)
      - ᴮQᴵ: Rotation matrix from the inertial to the body frame             |  [3, 3] (Static)
      - t:   Current time, as an epoch                                       |  Epoch
                                                                              
    Returns:
      - sᴵ: Satellite-sun vector in the inertial frame (noiseless)           | [3,] (Static)
      - sᴮ: Unit satellite-sun vector in the body frame (noiseless)          | [3,] (Static)
      - ecl: Eclipse factor (1.0 for full sun, 0.0 for full eclipse)         | Scalar
"""
function sun_measurement(x::STATE, ᴮQᴵ::SMatrix{3, 3, T, 9}, t::Epoch)::Tuple{SVector{3, T}, SVector{3, T}, T} where {T} 

    sᴵₑ = (sun_position(t))   # Earth-Sun vector in inertial frame
    sᴵ  = sᴵₑ - x.r           # Satellite-Sun vector in inertial frame 

    ecl = eclipse_cylindrical(vcat([x.r;]...), sᴵₑ)
    sᴮ  = ᴮQᴵ * (sᴵ / norm(sᴵ))  # Make it unit
    
    return sᴵ, sᴮ, ecl
end

# VERY slow (≈ 3ms); (NOTE - Clips to zero 𝑎𝑓𝑡𝑒𝑟 adding noise, so should do the same on est)
""" 
      Generates the photodiode currents and noisy current measurements.
    (NOTE that current generated cannot be negative so these are clipped at zero)

  Arguments:
    - sat:  (True) SATELLITE struct, which is used to generate sensor measurements     |  SATELLITE
                    (J  //  magnetometer  //  diodes)
    - alb:  struct containing REFL data and cell centers                               |  ALBEDO
                    (REFL  //  cell_centers_ecef)
    - x:    struct containing current environment state                                |  STATE 
                    (r  //  v  //  q  //  ω  //  β)
    - ecl:  Eclipse factor (1.0 for full sun, 0.0 for full eclipse)                    | Scalar
    - sᴵ:   Satellite-sun vector in the inertial frame  (noiseless)                    | [3,] (Static)
    - sᴮ:   Unit satellite-sun vector in the body frame (noiseless)                    | [3,] (Static)
    - σ_scale: (Opt) Scale that modifies ||diode current|| to get stdev of diode noise          |  Scalar 
                           ('σ_current_scale' in generate_measurements)
    - E_am₀:   (Opt) Irradiance of sunlight (TSI - visible & infrared). Default is 1366.0 W/m²  |  Scalar
      
  Returns:
    - I:   Currents generated by each photodiode                                           | [N,]  Scalar
    - Ĩ:   Measured (noisy) currents for each photodiode                                   | [N,]  Scalar
    - ηI:  Current noise added to the measurement                                          | [N,]  Scalar
              (Note that this is really just tracked for debugging purposes)
"""  
function diode_measurement(sat::SATELLITE{N, T}, alb::ALBEDO, x::STATE{T}, ecl::Real, sᴵ::SVector{3, T}, sᴮ::SVector{3, T}; 
                            σ_scale = 0.05, E_am₀ = 1366.9, use_albedo = true)::Tuple{SVector{N, T}, SVector{N, T}, SVector{N, T}} where {N, T}

    if (ecl == 0)   # No need to waste time calculating in eclipse - I = 0, Ĩ is just noise
        η = SVector{N, T}(rand(Normal(0.0, σ_scale * 0.001), N))
        return SVector{N, T}(zeros(N)), η, η
    else 

        C, α, ϵ = sat.diodes.calib_values, sat.diodes.azi_angles, sat.diodes.elev_angles 

        I  = zeros(T, N)   # Currents being generated by each photodiode
        ηI = zeros(T, N)   # Noise being applied to each photodiode
        Ĩ  = zeros(T, N)   # Noisy current measurements

        albedo_matrix = earth_albedo(x.r, sᴵ, alb.refl.data) 

        for i = 1:N 
            # Get photodiode surface normal
            surface_normal = SVector{3, T}((cos(ϵ[i])*cos(α[i])), (cos(ϵ[i])*sin(α[i])), sin(ϵ[i]))   

            diode_albedo = (use_albedo) ? compute_diode_albedo(albedo_matrix, alb.cell_centers_ecef, surface_normal, x.r) : 0.0

            # Calculate current, including noise and Earth's albedo 
            current = C[i] * (dot(surface_normal, sᴮ) + (diode_albedo / E_am₀) )
            current *= ecl
            η = rand(Normal(0.0, abs(σ_scale * current)))

            # Photodiodes don't generate negative current
            I[i] = (current < 1e-8) ? 0.0 : current 
            Ĩ[i] = (current + η < 1e-8) ? 0.0 : (current + η) 
            ηI[i] = η 
        end

        return SVector{N, T}(I), SVector{N, T}(Ĩ), SVector{N, T}(ηI)
    end
end 



####################
# HELPER FUNCTIONS #  Could probably go in a shared file
####################


#  Do I need dt? -> Not for sensor noise ?
"""
      Generates a [3 × 3] noise rotation matrix by generating a random axis of rotation and the 
    angle by which to rotate about the axis. This axis-angle representation is then converted
    to R ∈ SO(3) using Rodrigues' Rotation formula (similar to matrix exponential exp(ωθ))
    (Note that if the standard deviation provided is 0, the identity matrix is returned)

    Arguments:
      - σ:  Stdev of noise, in radians    |  Scalar 
      - dt: Time step                     |  Scalar 

    Returns: 
      - R:  Noise rotation matrix         |  [3, 3]  (Static)
"""
function rotation_noise(σ::T, dt::T)::SMatrix{3, 3, T, 9} where {T}

    if σ == 0
        return SMatrix{3, 3, T, 9}(I(3))
    else 
        v   = rand(Normal(0, σ), 3)
        mag = norm(v)
        v̂   = hat(v / mag) 

        R = I(3) + (v̂) * sin(mag * dt) + (v̂^2) * (1 - cos(mag * dt))  # Equivalent to exp(v̂)

        return SMatrix{3, 3, T, 9}(R)
    end
end

"""
      Generates the calibration matrix that corrects the un-calibrated magnetic
    field vector (in body frame). Accounts for scale factors and non-orthogonality 
    angles. 

    Arguments:
      - sat:  (True) SATELLITE struct, which is used to generate sensor measurements     |  SATELLITE
                      (J  //  magnetometer  //  diodes)
    Returns:
      - T:  Calibration matrix            |  [3, 3]

"""
function get_mag_calibration_matrix(sat::SATELLITE)

    a, b, c = sat.magnetometer.scale_factors 
    ρ, λ, ϕ = sat.magnetometer.non_ortho_angles 

    T = zeros(Float64, 3, 3)
    T[1, 1] = a 
    T[2, 1] = b * sin(ρ) 
    T[2, 2] = b * cos(ρ)
    T[3, 1] = c * sin(λ) 
    T[3, 2] = c * cos(λ) * sin(ϕ)
    T[3, 3] = c * cos(λ) * cos(ϕ)

    # T = [a         0.0              0.0;
    #      b*sin(ρ)  b*cos(ρ)         0.0;
    #      c*sin(λ)  c*sin(ϕ)*cos(λ)  c*cos(ϕ)*cos(λ)]

    return T
end

""" 
      Estimates the effect of Earth's albedo on a specific photodiode (by using the surface normal of that diode)
    (NOTE that this comes from the Springmann paper referenced in the Estimator section about photodiode calibration).
    There is also a duplciate version of this in the Estimators file, kept separate from the Simulator files.
       
      = cell_albedo * surface_normal^T * r_g,
    with r_g as a unit vector in the direction of the grid point on Earth

    Arguments:
      - data:   Albedo values for each cell on the Earth's surface  (calculated with EarthAlbedo)   | [num_lat x num_lon] 
      - cell_centers_ecef:  Center of each grid cell on Earth's surface (in ECEF)                   | [num_lat x num_lon x 3]
      - surface_normal: Photodiode surface normal                                                   | [3,] (Static)
      - sat_pos: Cartesian position of satellite                                                    | [3,] (Static)

    Returns:
      - diode_albedo: Total effect of albedo on specified photodiode              | Scalar
"""  
function compute_diode_albedo(data::Matrix{T}, cell_centers_ecef::Array{T, 3}, surface_normal::SVector{3, T}, sat_pos::SVector{3, T}) where {T}

    Nlat, Nlon = size(data)

    diode_albedo = 0.0
    r_g = zeros(T, 3)

    for r = 1:Nlat
        for c = 1:Nlon

            if data[r,c] != 0
                r_g .= view(cell_centers_ecef, r, c, :) .- sat_pos
                r_g .= r_g ./ norm(r_g)  # Make unit

                cell_albedo = (data[r,c] * dot(surface_normal, r_g))

                if cell_albedo > 0.0    # Can't be negative
                    diode_albedo += cell_albedo 
                end
            end
        end
    end
    
    return diode_albedo
end