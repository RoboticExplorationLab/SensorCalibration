# module simpleOrbit 
"""
    Generates a simplified orbit for a satellite, with functionality for (1) position only, 
(2) orientation only, or (3) position and orientation. Acceleration can include J₂ if desired, 
and control inputs are only allowed for the orientation.
"""

using LinearAlgebra, StaticArrays
using Plots
using Test
include("quaternions.jl")


######################
# DYNAMICS FUNCTIONS #
######################

# Dont allow for controls, and dont allow for just passing in state... (need a different dynamics interface...? like in energy)
#        OR use static...?

function dynamics(x, u = zeros(3), J::Matrix = zeros(1, 1); use_J₂::Bool = false, μ::Float64 = 3.9860044188e14, J₂::Float64 = 0.0010826269)
    N = size(x, 1)
    if N == 6       # Pure translation 
        return dynamics(x[1:3], x[4:6], use_J₂; μ = μ, J₂ = J₂)
    elseif N == 7   # Pure rotation
        return dynamics(J, x[1:4], x[5:7], u)
    elseif N == 13  # Both 
        return dynamics(x[1:3], x[4:7], x[8:10], x[11:13], J, u, use_J₂; μ = μ, J₂ = J₂)
    else 
        @warn "Invalid state size!"
    end
end

function dynamics(r, v, use_J₂::Bool = false; μ::Float64 = 3.9860044188e14, J₂::Float64 = 0.0010826269, Rₑ = 6378.1363e3)
    """ Translational dynamics ONLY (J₂ defaults to false but can be set to true)"""
    a = (-μ / (norm(r)^3)) * r 
    if use_J₂
        a += -((3.0 * J₂ * μ * (Rₑ^2)) / (2 * (norm(r)^5))) * [(1 - 5 * ((r[3]/norm(r))^2))*r[1]; 
                                                               (1 - 5 * ((r[3]/norm(r))^2))*r[2]; 
                                                               (3 - 5 * ((r[3]/norm(r))^2))*r[3]]
    end

    return [v; a]
end

function dynamics(J, q, ω, u  = zeros(3))
    """ Rotationaly dynamics ONLY (scalar first quaternion) """
    q̇ = 0.5 * L(q) * H * ω   # Scalar FIRST 
    ω̇ = J \ (u - cross(ω, J * ω))

    return [q̇; ω̇ ]
end

function dynamics(r, q, v, ω, J, u = zeros(3), use_J₂::Bool = false; μ::Float64 = 3.9860044188e14, J₂::Float64 = 0.0010826269)
    # Call the other functions to avoid retyping code...

    posDot = dynamics(r, v, use_J₂; μ, J₂)
    rotDot = dynamics(J, q, ω, u)

    # Extract and reorder the derivatives
    ṙ = @view posDot[1:3]
    v̇ = @view posDot[4:6]
    q̇ = @view rotDot[1:4]
    ω̇ = @view rotDot[5:7]

    return [ṙ; q̇; v̇; ω̇ ]
end



####################
# HELPER FUNCTIONS #
####################

# Just do it like [r, 0, 0]; [0, vcosθ, vsinθ]? or have them pass in a unit vector for v direction?
function init_orbit(r₀; θ = 0.0, Rₑ = 6378.1363e3, μ = 3.9860044188e14)
    """
           Initializes an simplified orbit. Allows a user to specify a variety of values, or just 
        leave them as default values.

        Arguments:
          - r₀: Initial position (Cartesian)
          - θ:  Direction of initial velocity, described as yaw about the Earth/Sat vector
          - Rₑ: Radius of the Earth   (defaults to meters)
          - μ:  Standarde gravitational constant
    """

    # Cart coordinates -> how do i define a plane perp to a normal?

    if norm(r₀) < Rₑ
        @warn "Satellite within Earth's radius. Try again..."
        return nothing
    end

    _a  = _r0 # Semi-major axis of (elliptical) orbit (circular orbit if a = r0)
    _v0 = sqrt(_μ*( (2/_r0) - (1/_a))) 

end

function random_orbit()
end

function get_inertia(; m = 1.0, l = 1.0, w = 1.0, h = 1.0)
    """ Calculates inertia matrix for a cuboid """
    Ih = (m / 12.0) * (l^2 + w^2) 
    Iw = (m / 12.0) * (l^2 + h^2)
    Id = (m / 12.0) * (h^2 + w^2)

    return diagm([Ih, Iw, Id])
end


####################
# ENERGY EQUATIONS #
####################

# Not very multiple-dispatch-y (need static vectors or struct tho...)
function energy(x::Vector, m::Float64 = 0.0, J::Matrix = zeros(3, 3); g = 9.8, G = 6.6743e-11, Mₑ = 5.972e24) # kg 

    if size(x, 1) == 6          # [r, v]
        r = @view x[1:3]
        v = @view x[4:6]
        # Eₚ = m * g * norm(r)
        Eₚ = -G * m * Mₑ / norm(r)
        Eₖ = 0.5 * m * v' * v

        return Eₚ + Eₖ        

    elseif size(x, 1) == 7      # [q, ω]
        q = @view x[1:4]
        ω = @view x[5:7]

        return 0.5 * ω' * J * ω

    elseif size(x, 1) == 13     # [r, q, v, ω]

        r = @view x[1:3]
        v =  @view x[8:10]
        ω = @view x[11:13]

        # Eₚ = m * g * norm(r)
        Eₚ = -G * m * Mₑ / norm(r)
        Eₖ = 0.5 * m * v' * v + 0.5 * ω' * J * ω

        return Eₚ + Eₖ

    else 
        @warn "Invalid state size!"
        return nothing
    end
end

function energy(m::Float64, r::Vector, v::Vector; g = 9.8)
    """ Energy for a translational orbit """

    Eₚ = m * g * norm(r)
    Eₖ = 0.5 * m * dot(v, v)

    return Eₚ + Eₖ
end

function energy(J::Matrix, ω::Vector)
    """ Energy for a rotating satellite (no potential energy) """
    return 0.5 * ω' * J * ω 
end

function energy(m::Float64, r::Vector, v::Vector, I::Matrix, ω::Vector; g = 9.8)
    # Eₚ = m * g * norm(r)
    # Eₖ = 0.5 * m * dot(v, v) + 0.5 * ω' * I * ω 
    return energy(m, r, v) + energy(I, ω)  # Sum together translational and rotational energy 
end





# Include alternative integrators for fun? RK4, Implicit Midpoint, etc...?


# CONSTANTS (DON'T EXPLICITELY EXPORT FROM PACKAGE)
_Rₑ = 6378.1363e3
_μ  = 3.9860044188e14

function rk4(f, x, u, h, J::Matrix = zeros(1, 1); kwargs...)
    """ Modified rk4 function """
    k₁ = h * f(x,        u, J; kwargs...)
    k₂ = h * f(x + k₁/2, u, J; kwargs...)
    k₃ = h * f(x + k₂/2, u, J; kwargs...)
    k₄ = h * f(x + k₃,   u, J; kwargs...)

    x⁺ = x + (1/6.0) * (k₁ + 2 * k₂ + 2 * k₃ + k₄)

    # Normalize the quaternion
    if size(x⁺, 1) == 7 
        q, ω = x⁺[1:4], x⁺[5:7]
        q = q / norm(q)
        return [q; ω]
        
    elseif size(x⁺, 1) == 13 
        r, q, v, ω = x⁺[1:3], x⁺[4:7], x⁺[8:10], x⁺[11:13]

        q = q / norm(q)
        return [r; q; v; ω]
    end

    return x⁺
end

function orbit(x0; N = 20000, h = 5.0)
    Xhist = zeros(6, N)
    Xhist[:, 1] = x0 
    for i = 1:N - 1 
        Xhist[:, i + 1] = rk4(dynamics, Xhist[:, i], zeros(3), h; use_J₂ = false)
    end    

    display(plot(Xhist[1:3, :]'))
end

if false 
    @testset "Example (Translation)" begin

        function orbit(x0; N = 20000, h = 5.0, use_J₂ = false, plt = false)
            Xhist = zeros(6, N)
            Xhist[:, 1] = x0 
        
            e = zeros(N)
            e[1] = energy(Xhist[:, 1], _m)
        
            for i = 1:N - 1 
                Xhist[:, i + 1] = rk4(dynamics, Xhist[:, i], zeros(3), h; use_J₂ = use_J₂)
                e[i + 1] = energy(Xhist[:, i + 1], _m)
            end    
        
            μe = sum(e) / N
            if plt
                display(plot(Xhist[1:3, :]'))
                display(plot(e))
            end
            
            err = [abs(e[i] - μe) for i = 1:N]
            @test maximum(err) < abs(μe * 0.002)
        end

        # Example 1 (r in +X, v in +Y, circular)
        _r0 = (550+6371)*(10^(3))   # Distance from two center of masses (100*km) 
        _a  = _r0                   # Semi-major axis of (elliptical) orbit (circular orbit if a = r0)
        _v0 = sqrt(_μ*( (2/_r0) - (1/_a))) 

        _m  = 1.0                   # sat mass (Kg) 
        _r0 = [_r0; 0.0; 0.0]
        _v0 = [0.0; _v0; 0.0]

        x0 = [_r0; _v0]
        orbit(x0)



        # Example 2 (r in +Y, v in +Z, elliptical)
        _r0 = (1000+6371)*(10^(3))   
        _a  = 2.0 * _r0                  
        _v0 = sqrt(_μ*( (2/_r0) - (1/_a))) 

        _m  = 2.2                   
        _r0 = [0.0; _r0; 0.0]
        _v0 = [0.0; 0.0; _v0]

        x0 = [_r0; _v0]
        orbit(x0)


        # Example 3 (r in +Z, v in +X/+Y, elliptical)
        _r0 = (2000+6371)*(10^(3))  
        _a  = 2.0 * _r0              
        _v0 = sqrt(_μ*( (2/_r0) - (1/_a))) 

        _m  = 0.25                   
        _r0 = [0.0; 0.0; _r0]
        _v0 = [_v0 * sin(pi/2); _v0 * cos(pi/2); 0.0]

        x0 = [_r0; _v0]
        orbit(x0)


        # Example 4 (r in +Z, v in +X/+Y, elliptical)
        _r0 = (4000+6371)*(10^(3))  
        _a  = 2.3 * _r0              
        _v0 = sqrt(_μ*( (2/_r0) - (1/_a))) 

        _m  = 11.023                   
        _r0 = [0.0; _r0; 0.0]
        _v0 = [_v0 * sin(pi/8); 0.0; _v0 * cos(pi/8)]

        x0 = [_r0; _v0]
        orbit(x0; use_J₂ = true, N = 100000)

    end

    @testset "Example (Rotation)" begin

        function orbit(x0, J; N = 10000, h = 5.0, plt = false)
            Xhist = zeros(7, N)
            Xhist[:, 1] = x0 
        
            e = zeros(N)
            e[1] = energy(Xhist[:, 1], 0.0, J)
        
            for i = 1:N - 1 
                Xhist[:, i + 1] = rk4(dynamics, Xhist[:, i], zeros(3), h, J)
                e[i + 1] = energy(Xhist[:, i + 1], 0.0, J)
            end    
        
            μe = sum(e) / N
            if plt
                display(plot(Xhist[1:4, :]', title = "Quats"))
                display(plot(Xhist[5:7, :]', title = "ω"))
                display(plot(e, title = "Energy"))
            end
            
            err = [abs(e[i] - μe) for i = 1:N]
            @test maximum(err) < abs(μe * 0.00001)
        end

        # Example 1 
        _q0 = [1, 0, 0, 0]
        _ω0 = [0.0, 1.0, -1.0]
        _J = get_inertia()

        x0 = [_q0; _ω0]
        orbit(x0, _J; plt = true)

        # Example 2
        _q0 = [0, 0, 1, 0]
        _ω0 = [0.5, 1.0, -1.0]
        _J = get_inertia()

        x0 = [_q0; _ω0]
        orbit(x0, _J; N = 100, plt = true)

        # Example 3
        _q0 = [1, 0, 0, 0]
        _ω0 = [0.0, 1.0, 0.0]
        _J = get_inertia()

        x0 = [_q0; _ω0]
        orbit(x0, _J; N = 1000, plt = true)

        # Example 4
        _q0 = [0, sqrt(2)/2, sqrt(2)/2, 0]
        _ω0 = [0.5, 1.0, -1.0]
        _J = get_inertia()

        x0 = [_q0; _ω0]
        orbit(x0, _J; N = 1000, plt = true)

        # Example 5
        _q0 = [0, sqrt(2)/2, sqrt(2)/2, 0]
        _ω0 = [0.5, 1.0, -1.0]
        _J = get_inertia()

        x0 = [_q0; _ω0]
        orbit(x0, _J; N = 40000, plt = true)
    end

    @testset "Example (full state)" begin 

        function orbit(x0, J, m; N = 10000, h = 5.0, plt = false)
            Xhist = zeros(size(x0, 1), N)
            Xhist[:, 1] = x0 
        
            e = zeros(N)
            e[1] = energy(Xhist[:, 1], m, J)
        
            for i = 1:N - 1 
                Xhist[:, i + 1] = rk4(dynamics, Xhist[:, i], zeros(3), h, J)
                e[i + 1] = energy(Xhist[:, i + 1], m, J)
            end    
        
            μe = sum(e) / N
            if plt
                display(plot(Xhist[1:3, :]', title = "Pos"))
                display(plot(Xhist[4:7, :]', title = "Quat"))
                display(plot(e, title = "Energy"))
            end
            
            err = [abs(e[i] - μe) for i = 1:N]
            @test maximum(err) < abs(μe * 0.001)
        end

        _r0 = (550+6371)*(10^(3))   # Distance from two center of masses (100*km) 
        _a  = _r0                   # Semi-major axis of (elliptical) orbit (circular orbit if a = r0)
        _v0 = sqrt(_μ*( (2/_r0) - (1/_a))) 
        _q0 = [1; 0; 0; 0]
        _ω0 = [0.0; 0.0; 0.0]

        _m  = 1.0                   # sat mass (Kg) 

        _r0 = [_r0; 0.0; 0.0]
        _v0 = [0.0; _v0; 0.0]

        x0 = [_r0;  _q0;  _v0;  _ω0]
        J = get_inertia(; m = _m)
        orbit(x0, J, _m; N = 10000, plt = false)

        #####
        _r0 = (550+6371)*(10^(3))   
        _a  = _r0                   
        _v0 = sqrt(_μ*( (2/_r0) - (1/_a))) 
        _q0 = [1; 0; 0; 0]
        _ω0 = [0.1; -0.1; 0.0]

        _m  = 1.5                   # sat mass (Kg) 

        _r0 = [_r0; 0.0; 0.0]
        _v0 = [0.0; 0.0; _v0]

        x0 = [_r0;  _q0;  _v0;  _ω0]
        J = get_inertia(; m = _m)
        orbit(x0, J, _m; N = 10000, plt = true)

        #####
        _r0 = (300+6371)*(10^(3))   
        _a  = _r0                   
        _v0 = sqrt(_μ*( (2/_r0) - (1/_a))) 
        _q0 = [1; 0; 0; 0]
        _ω0 = [0.1; -0.1; 0.0]

        _m  = 1.5                   # sat mass (Kg) 

        _r0 = [_r0; 0.0; 0.0]
        _v0 = [0.0; 0.0; _v0]

        x0 = [_r0;  _q0;  _v0;  _ω0]
        J = get_inertia(; m = _m)
        orbit(x0, J, _m; N = 15000, plt = true)

        #####
        _r0 = (900+6371)*(10^(3))   
        _a  = 2.0 * _r0                   
        _v0 = sqrt(_μ*( (2/_r0) - (1/_a))) 
        _q0 = [0; 0; sqrt(2)/2; sqrt(2)/2]
        _ω0 = [0.1; -0.1; 1.0]

        _m  = 10.0                   # sat mass (Kg) 

        _r0 = [_r0; 0.0; 0.0]
        _v0 = [0.0; 0.0; _v0]

        x0 = [_r0;  _q0;  _v0;  _ω0]
        J = get_inertia(; h = 2.0,  m = _m)
        orbit(x0, J, _m; N = 8000, plt = true)
    end
end


# end


