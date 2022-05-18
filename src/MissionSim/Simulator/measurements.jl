# [src/MissionSim/Simulator/measurements.jl]

"""
    Contains the code used in the Simulator to generate 
  the measurements for each of the sensors.

  To Do:
    - diode_meas is reaaaaal slow (2.8ms / 19 / 413 KiB)  -> earth_albedo takes up 2.241 ms of that, but if we downsample albedo before...
    - IGRF in mag_meas is also slow ∼ μs
    - Update all the comments
"""

# update comments. These noise parameters are mostly random
function generate_measurements(sat::SATELLITE, alb::ALBEDO, x::STATE, t::Epoch, dt::T; 
             E_am₀ = 1366.9, σβ = deg2rad(3.0), μ_gyro_scale = 0.05, σ_gyro_scale = 0.005, 
             μr = 5e4, σr = 5e3, σ_sun = deg2rad(3.0), σ_current_scale = 0.05) where {T}
  
  """
   (All the noise)
    E_am₀ =  1366.9 # Irradiance of sunlight (TSI - visible & infrared), W/m^2 
  """
  
  ᴮQᴵ = quat2rot(x.q)  # Compute once and pass in 

  # Call each sensor
  Bᴵ, Bᴮ, B̃ᴮ  = mag_measurement(sat, x, ᴮQᴵ, t, dt; σ = σβ) 
  w, w̃, ηw    = gyro_measurement(x; μ_scale = μ_gyro_scale, σ_scale = σ_gyro_scale) 
  r, r̃, ηr    = pos_measurement(x; μ = μr, σ = σr) 
  sᴵ, sᴮ, ecl = sun_measurement(x, ᴮQᴵ, t, dt; σ = σ_sun)
  I, Ĩ, ηI    = diode_measurement(sat, alb, x, ecl, sᴵ, sᴮ; σ_scale = σ_current_scale, E_am₀ = E_am₀)

  # Store the outputs 
  sensors = SENSORS(B̃ᴮ, Ĩ, w̃, r̃)
  truth   = GROUND_TRUTH(t, Bᴵ, sᴵ, sᴮ, Bᴮ, I)  # w, r are stored in state
  noise = NOISE(ηI, ηw, ηr)

  return truth, sensors, ecl, noise
end

# test, alloc, type;    Comment   (Speed up IGRF13? Takes ~μs)
function mag_measurement(sat::SATELLITE, x::STATE, ᴮQᴵ::SMatrix{3, 3, T, 9}, t::Epoch, dt::T; σ = deg2rad(3.0)) where {T}
  """
    Generates the true magnetic field vector in inertial and body frames, 
        as well as the noisy sensor measurement.

    Arguments:
    - pos: Current (x, y, z) position of the satellite in ECEF                                          | [3,] (SVector)
    - t: Current time, as an epoch                                                                      | Epoch
    - sat: TRUE satellite data, which is used to generate sensor measurements                           | SATELLITE
                J (inertia matrix, [3 x 3])          //  magnetometer (calibration values, [3,]),
                diodes (photodiode parameters [3,])  //  state (attitude and bias, [7,])
    - q: Quaternion from inertial frame to body frame                                                   | [4,] (SVector)
    - dt: Time step                                                                                     | Scalar
    
    Returns:
    - Bᴵ: Magnetic field vector in inertial frame (noiseless)                                           | [3,] (SVector)
    - Bᴮ: Magnetic field vector in body frame (noiseless)                                               | [3,] (SVector)
    - B̃ᴮ: Measured magnetic field vector in body frame (noisy, biased, and scaled)                      | [3,] (SVector)
  """
  Bᴵ = SVector{3, Float64}(IGRF13(x.r, t))   # Mag vector in inertial frame 
  Bᴮ = ᴮQᴵ * Bᴵ 

  # Add in noise
  η_mag = rotation_noise(σ, dt)
  mag_calibration_matrix = get_mag_calibration_matrix(sat)
  B̃ᴮ = (mag_calibration_matrix * η_mag * Bᴮ) + sat.magnetometer.bias 

  return Bᴵ, Bᴮ, B̃ᴮ
end

# test, alloc, type;  Verify scale values, update comments; Noise scales with value (?)
function gyro_measurement(x::STATE; μ_scale = 0.05, σ_scale = 0.005) 
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
  ηω = rand(Normal(μ_scale * norm(x.ω), σ_scale * norm(x.ω)), 3)  # Noise scales with actual value

  ω̃  = x.ω .+ x.β .+ ηω

  return x.ω, ω̃ , ηω
end

# test, alloc, type;   Verify noise values, comments; (Note that r is tracked in state, but this allows for some adjustments)
function pos_measurement(x::STATE; μ = 5e4, σ = 5e3)
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

  ηr = SVector{3, Float64}(rand(Normal(μ, σ), 3))
  r̃ = x.r + ηr

  return x.r, r̃, ηr
end

# test, allocs, type;  update comments NOTE i am not scaling sᴵ by ecl anymore; Verify noise values, update comments; dt may not be needed
function sun_measurement(x::STATE, ᴮQᴵ::SMatrix{3, 3, T, 9}, t::Epoch, dt::T; σ = deg2rad(3.0))::Tuple{SVector{3, T}, SVector{3, T}, T} where {T} 
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
  sᴵₑ = (sun_position(t)) 
  sᴵ = sᴵₑ - x.r 

  ecl = eclipse_cylindrical(vcat([x.r;]...), sᴵₑ)

  # Avoid divide-by-zero problem 
  if ecl < 0.001 
    return sᴵ, SVector(0.0, 0.0, 0.0), 0.0

  else

    η_sun = rotation_noise(σ, dt)   
    𝐬ᴵ = sᴵ / norm(sᴵ)
    𝐬ᴮ = (η_sun * (ᴮQᴵ * 𝐬ᴵ))

    return sᴵ, 𝐬ᴮ, ecl
  end
end

# VERY slow (≈ 3ms)
# Verify noise values, comments; noise is scaled by current  (Clips to zero 𝑎𝑓𝑡𝑒𝑟 adding noise, so should do the same on est)
function diode_measurement( sat::SATELLITE{N, T}, alb::ALBEDO, x::STATE{Sₓ, T}, ecl::Real, sᴵ::SVector{3, T}, sᴮ::SVector{3, T}; 
                            σ_scale = 0.05, E_am₀ = 1366.9)::Tuple{SVector{N, T}, SVector{N, T}, SVector{N, T}} where {N, Sₓ, T}
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
    - sᴵ: Satellite-Sun vector in the inertial frame (noiseless)                                      | [3,]
    - sᴮ: Unit satellite-sun vector in the body frame (noiseless)                                     | [3,]
    - CONSTANTS: Various constants needed (Earth radius, μ , E_am₀)                                 | [3,]       
    
    Returns:
    - I: Currents generated by each photodiode                                           | [N,]
    - Ĩ: Measured (noisy) currents for each photodiode                                   | [N,]
    - ηI: Current noise added to the measurement                                        | [N,]
        (Note that this is really just tracked for debugging purposes)
  """  

  if (ecl == 0)   # I = 0, Ĩ is just noise
    η = SVector{N, T}(rand(Normal(0.0, σ_scale * 0.001), N))
    return SVector{N, T}(zeros(N)), η, η
  else 

    C, α, ϵ = sat.diodes.calib_values, sat.diodes.azi_angles, sat.diodes.elev_angles 

    I  = zeros(T, N)   # Currents being generated by each photodiode
    ηI = zeros(T, N)
    Ĩ  = zeros(T, N)

    albedo_matrix = earth_albedo(x.r, sᴵ, alb.refl.data) 

    for i = 1:N 
      surface_normal = SVector{3, T}((cos(ϵ[i])*cos(α[i])), (cos(ϵ[i])*sin(α[i])), sin(ϵ[i]))  # Photodiode surface normal 

      diode_albedo = compute_diode_albedo(albedo_matrix, alb.cell_centers_ecef, surface_normal, x.r)  # N * (0.1 ms, 2 alloc )

      current = C[i] * dot(surface_normal, sᴮ) + (C[i] / E_am₀) * diode_albedo  # Calculate current, including noise and Earth's albedo 
      current *= ecl
      η = rand(Normal(0.0, abs(σ_scale * current)))

      # Photodiodes don't generate negative current
      I[i] = (current < 0) ? 0.0 : current 
      Ĩ[i] = (current + η < 0.0) ? 0.0 : (current + η) 
      ηI[i] = η 
    end

    return SVector{N, T}(I), SVector{N, T}(Ĩ), SVector{N, T}(ηI)
  end
end 




####################
# HELPER FUNCTIONS #
####################
# Could probably go in a shared file

# Tests (?), allocs, type;    Comments; Do I need dt? -> Not for sensor noise 
function rotation_noise(σ::T, dt::T)::SMatrix{3, 3, T, 9} where {T}
  """
      Generates a [3 × 3] noise rotation matrix by generating a random axis of rotation and the 
    angle by which to rotate about the axis. This axis-angle representation is then converted
    to R ∈ SO(3) using Rodrigues' Rotation formula (similar to matrix exponential exp(ωθ))
    (Note that if the standard deviation provided is 0, the identity matrix is returned)
  
    σ is in radians
  """
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

# Tests, type, allocs;   Comments  (look at satDyn.jl comments)
function get_mag_calibration_matrix(sat::SATELLITE)
  """ Generates the calibration matrix that alters the measured magnetic field vector in body frame """
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

# Tests, type, allocs;   Comments
function compute_diode_albedo(data::Matrix{T}, cell_centers_ecef::Array{T, 3}, surface_normal::SVector{3, T}, sat_pos::SVector{3, T}; lat_step::Int = 1, lon_step::Int = 1) where {T}
  """ 
      Estimates the effect of Earth's albedo on a specific photodiode (by using the surface normal of that diode)
          = cell_albedo * surface_normal^T * r_g, with r_g as a unit vector in the direction of the grid point on Earth
          (NOTE this comes from the Springmann paper referenced in the estimator section about photodiode calibration)

      Arguments:
      - data: Albedo values for each cell on the Earth's surface  (calculated with EarthAlbedo)   | [num_lat x num_lon] 
      - surface_normal: Photodiode surface normal                                 | [3,]
      - sat_pos: Cartesian position of satellite                                  | [3,]

      Returns:
      - diode_albedo: Total effect of albedo on specified photodiode              | Scalar
  """    
  Nlat, Nlon = size(data)

  diode_albedo = 0.0
  r_g = zeros(T, 3)

  for r = 1:lat_step:Nlat
      for c = 1:lon_step:Nlon

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