# [src/MissionSim/Simulator/dynamics.jl]

"""
    Contains the code used in the Simulator to update 
  the state of the satellite.

  To Do:
    - Add in drag that is dependent on orientation?
    - Get accurate gyro bias values
    - Speed up accel_perturbations? (it is ≈30x slower than the rest combined)
    - Migrate to using STATE struct (tests; remove old, or keep...?)
"""


"""
    dynamics(J, x, u, t; Rₑ, σβ)
    
  Propagates the state dynamics for a satellite, where the state is defined as 
  [position, velocity, scalar-first quaternion, angular velocity, gyro bias]
    'x = [r⃗, v⃗, (q₀, q⃗), ω, β]'

  Includes gravity from spherical harmonics, uniform drag, SRP, and third body 
  gravity from the sun and moon. 
  (NOTE that there is also a simplified version that takes in inertia rather than a satellite)

  Arguments:
    - J:   Inertia matrix (as a Static Matrix)                                    |  [3, 3]  (Static) 
    - x:   state vector (as listed above)                                         |  [16,]
    - u:   control input (rotation only)                                          |  [3,]
    - t:   Current time (as an epoch)                                             |  Epoch
    - Rₑ:  (Optional) Radius of the Earth (default is meters)                     |  Scalar
    - σβ:  (Optional) Noise for the gyro bias                                     |  Scalar

  Returns:
    - x:   updated state vector                                                   |  [16,]
"""
function dynamics(J::SMatrix{3, 3, T, 9}, x::STATE{S, T}, u::SVector{3, T}, t::Epoch; Rₑ = 6378136.3, σβ = deg2rad(0.1), kwargs...) where {S, T}

  if norm(x.r) < Rₑ                  
      error("ERROR: Impact at time $t")
  end

  ṙ = x.v   
  v̇ = accel_perturbations(t, x.r, x.v; kwargs...)    # <- THIS FUNCTION is reeeeally slow
  q̇ = qdot(x.q, x.ω)  
  ω̇ = J \ (u - cross(x.ω, J * x.ω))   
  β̇ = SVector{3, T}(rand(Normal(0.0, σβ), 3))  

  return STATE(ṙ, v̇, q̇, ω̇ , β̇ )
end

""" Alternate function call that uses an array for state rather than a struct """
function dynamics(J::SMatrix{3, 3, T, 9}, x::SVector{S, T}, u::SVector{3, T}, t::Epoch; Rₑ = 6378136.3, σβ = deg2rad(0.1), kwargs...)::SVector{S, T} where {S, T}

    if norm(@view x[1:3]) < Rₑ                  
        error("ERROR: Impact at time $t")
    end

    r, v, q, ω, β = view(x, 1:3), view(x, 4:6), view(x, 7:10), view(x, 11:13), view(x, 14:16)

    ṙ = v   
    v̇ = accel_perturbations(t, x[1:6]; kwargs...)    # <- THIS FUNCTION is reeeeally slow
    q̇ = qdot(q, ω)  
    ω̇ = J \ (u - cross(ω, J * ω))   
    β̇ = rand(Normal(0.0, σβ), 3)    

    return SVector{S, T}([ṙ; v̇; q̇; ω̇ ; β̇ ])
end

""" Alternate function call for dynamics that allows for a satellite struct to be used, rather than just J """
function dynamics(sat::SATELLITE, x::SVector{S, T}, u::SVector{3, T}, t::Epoch; Rₑ = 6378136.3, σβ = deg2rad(0.1))::SVector{S, T} where {S, T}

    return dynamics(sat.J, x, u, t; Rₑ = Rₑ, σβ = σβ)
end

"""
    accel_perturbations(epc, x; mass, area_drag, coef_drag, area_srp, coef_srp, n_grav, m_grav)

  Generates the acceleration for a spacecraft in LEO. Accounts for a variety of factors, 
  including spherical harmonics, atmospheric drag, SRP, and thirdbody from sun and moon

  ForwardDiff friendly, written by Kevin
"""
function accel_perturbations(epc::Epoch, x::Array{<:Real} ;
    mass::Real=1.0, area_drag::Real=0.01, coef_drag::Real=2.3,
    area_srp::Real=1.0, coef_srp::Real=1.8,
    n_grav::Integer=10, m_grav::Integer=10, third_body::Bool = true )::Vector{Float64}
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
    nu = eclipse_cylindrical(x, r_sun)  # conical doesn't work correctly
    a += nu * accel_srp(x, r_sun, mass, area_srp, coef_srp)

    if third_body
        # third body sun
        a += accel_thirdbody_sun(x, r_sun)

        # third body moon
        a += accel_thirdbody_moon(x, r_moon)
    end

    return a
end
function accel_perturbations(epc::Epoch, r::SVector{3, T}, v::SVector{3, T} ;
  mass::Real=1.0, area_drag::Real=0.01, coef_drag::Real=2.3,
  area_srp::Real=1.0, coef_srp::Real=1.8,
  n_grav::Integer=10, m_grav::Integer=10, third_body::Bool = true )::SVector{3, T} where {T}
  """Accelerations for spacecraft in LEO, ForwardDiff friendly (From Kevin)"""

  # Extract position and velocity
  # r = x[1:3]
  # v = x[4:6]

  # These functions don't like static, so we gotta adjust
  x = vcat([r; v]...)
  r = vcat([r;]...)
  v = vcat([v;]...)

  # Compute ECI to ECEF Transformation -> IAU2010 Theory
  PN = bias_precession_nutation(epc)
  E  = earth_rotation(epc)
  W  = polar_motion(epc)
  R  = W * E * PN

  # Compute sun and moon position
  r_sun  = sun_position(epc)
  r_moon = moon_position(epc)

  # Compute acceleration (eltype(x) makes this forward diff friendly)
  a = zeros(T, 3)

  # spherical harmonic gravity
  a += accel_gravity(x, R, n_grav, m_grav)

  # atmospheric drag
  ρ = density_harris_priester(epc,r)
  a += accel_drag([r;v],ρ,mass, area_drag, coef_drag, Array{Real, 2}(PN))

  # SRP
  nu = eclipse_cylindrical(x, r_sun)  # conical doesn't work correctly
  a += nu * accel_srp(x, r_sun, mass, area_srp, coef_srp)

  if third_body
      # third body sun
      a += accel_thirdbody_sun(x, r_sun)

      # third body moon
      a += accel_thirdbody_moon(x, r_moon)
  end

  return SVector{3, T}(a)
end

"""
    rk4(model, x, u, t, h)

  Modified RK4 function for integrating state.
  (Note there is also a simplied version that takes in J rather than a model)
"""
function rk4(J::SMatrix{3, 3, T, 9}, x::STATE{S, T}, u::SVector{3, T}, t::Epoch, h::Real; kwargs...) where {S, T} 
  k₁ = h * dynamics(J, x, u, t; kwargs...)
  k₂ = h * dynamics(J, x + k₁/2, u, t + h/2; kwargs...)
  k₃ = h * dynamics(J, x + k₂/2, u, t + h/2; kwargs...)
  k₄ = h * dynamics(J, x + k₃,   u, t + h; kwargs...)

  return (x + (1/6) * (k₁ + 2 * k₂ + 2 * k₃ + k₄) )
end

""" Alternate function call for rk4 that takes in a vector for state rather than a struct """
function rk4(J::SMatrix{3, 3, T, 9}, x::SVector{S, T}, u::SVector{3, T}, t::Epoch, h::Real; kwargs...) where {S, T} 
    k₁ = h * dynamics(J, x, u, t; kwargs...)
    k₂ = h * dynamics(J, x + k₁/2, u, t + h/2; kwargs...)
    k₃ = h * dynamics(J, x + k₂/2, u, t + h/2; kwargs...)
    k₄ = h * dynamics(J, x + k₃,   u, t + h; kwargs...)

    return (x + (1/6) * (k₁ + 2 * k₂ + 2 * k₃ + k₄) )
end

""" Alternate function call for rk4 that takes in a satellite struct rather than inertia"""
function rk4(model::SATELLITE, x::SVector{S, T}, u::SVector{3, T}, t::Epoch, h::Real) where {S, T} 
    k₁ = h * dynamics(model, x, u, t)
    k₂ = h * dynamics(model, x + k₁/2, u, t + h/2)
    k₃ = h * dynamics(model, x + k₂/2, u, t + h/2)
    k₄ = h * dynamics(model, x + k₃,   u, t + h)

    return (x + (1/6) * (k₁ + 2 * k₂ + 2 * k₃ + k₄) )
end;


