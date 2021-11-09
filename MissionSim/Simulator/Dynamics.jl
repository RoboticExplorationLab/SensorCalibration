""" Dynamics.jl: 

    Contains all relevant functions to propagate the state dynamics for 
        the satellite.

    TODO:
    - Eliminate magic numbers (dynamics, )
    - Something from ARDS rather than normalizing quaternion in RK4
"""

# TODO - Something fancier than normalizing the quaternion in RK4
function run_dynamics(model, x, u, t, h)
    """ Modified rk4 function to update satellite state """

    k1 = h * dynamics(model, x, u, t)
    k2 = h * dynamics(model, x + k1/2, u, t + h/2)
    k3 = h * dynamics(model, x + k2/2, u, t + h/2)
    k4 = h * dynamics(model, x + k3,   u, t + h)

    x_next = (x + (1/6)*(k1+2*k2+2*k3 + k4))

    x_next[7:10] /= norm(x_next[7:10]) # Normalize quaternion

    return x_next
end

# TODO σβ -> Magic numbers
function dynamics(sat::SATELLITE, x, u, t)
    """ Propagates the state dynamics ([r⃗, v⃗, (q⃗, q₀), ω, β]) """

    if norm(x[1:3]) < _Re # Satellite has crashed!
        error("Impact at time $t")
    end

    ṙ = @view x[4:6]
    v̇ = 0 * accel_perturbations(t, x[1:6]);

    w = @view x[11:13] # angular velocity 
    q̇ = 0.5 * qmult(x[7:10], [w; 0]) # Scalar LAST quaternion

    ẇ = (sat.J^(-1)) * (u - cross(w, (sat.J*w)))

    δβgyro = rand(Normal(0.0, σ_gyro_bias), 3) # β̇  looks funny so I am using δ
    δβmag  = rand(Normal(0.0, σ_mag_bias), 3)

    return [ṙ[:]; v̇[:]; q̇[:]; ẇ[:]; δβgyro[:]; δβmag[:]]
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

