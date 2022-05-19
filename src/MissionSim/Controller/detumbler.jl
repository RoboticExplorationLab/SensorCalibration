# [src/MissionSim/Controller/detumbler.jl]

# TODO: Calculate the optimal value of κ (Used to be 7e-6, but I estimate it should be large...)

using StaticArrays, LinearAlgebra

"""
      DETUMBLER(ω::SVector, Bᴮ::SVector, dt::Float)

    Simple controller that slows the rotation of the CubeSat.
"""
struct DETUMBLER{S, T}

    ω::SVector{S, T}   # Angular velocity
    Bᴮ::SVector{S, T}  # Measured magnetic field vector (not necessarily unit)
    dt::T              # Time step

    function DETUMBLER(ω::SVector{S, T}, Bᴮ::SVector{S, T}, dt::T) where {S, T} 
        """ Constructor """
        new{S, T}(ω, Bᴮ, dt)
    end

    function DETUMBLER(_ω, _Bᴮ, _dt) 
        """ Constructor for incorrect types """
        # @warn "\tDETUMBLER object made without static vectors!"
        S, T = length(_ω), typeof(_ω[1])
        _ω  = SVector{S, T}(_ω)
        _Bᴮ = SVector{S, T}(_Bᴮ)

        new{S, T}(_ω, _Bᴮ, _dt)
    end

    function DETUMBLER() 
        """ Random Constructor """
        ω  = SVector{3, Float64}(randn(3))
        Bᴮ = SVector{3, Float64}(randn(3))
        dt = round(5.0 * rand(), digits = 2)
    end
end;

"""
    generate_command(ctrl::DETUMBLER)

  Generates an appropriate control command corresponding to the provided 
  controller style. Defaults to B-Cross, but B-Dot is being worked on. 

  (NOTE that this function is here to provide a common interface across  
  all controllers.)
"""
function generate_command(ctrl::DETUMBLER; func = b_cross, kwargs...)
    return func(ctrl.ω, ctrl.Bᴮ; kwargs...)
end;

"""
    b_cross(ω, Bᴮ; κ)

  Provides asymptotic convergence from arbitrary initial conditions to zero angular velocity.
  Requires the spin vector to not be parallel to the magnetic field vector, which (apparently) 
  isn't a problem in practice. Does NOT require full (3dof) actuation, and may require a time-
  varying magnetic field vector. Control gain κ defaults to 7e-6, which was determined through 
  experiment in simulation.

  (NOTE this is based off of "Magnetic Detumbling of a Rigid Spacecraft" (Avanzini).)
"""
function b_cross(ω::SVector{S, T}, Bᴮ::SVector{S, T}; κ::T = 5e-4) where {S, T}   # κ = 7e-6
    """ Simple B-cross controller that commands a torque opposite the direction of rotation """
    B̂ᴮ = Bᴮ / norm(Bᴮ)  # Make unit 
 
    M = -κ * (I - (B̂ᴮ)*(B̂ᴮ)') * ω 

    return M
end;

"""
    b_dot(ω, Bᴮ, dt; κ)

  A B-Dot controller that commands a torque opposite the direction of rotation. 
  Control gain κ defaults to 7e-6, which was what was determined for the B-Cross controller.

  (NOTE this is not currently working)
"""
function b_dot(ω::SVector{S, T}, Bᴮ::SVector{S, T}, dt::T, κ::T = 7e-6) where {S, T}
    @warn "\tWarning! B-Dot controller is still under development"

    # if first_time
    #     global first_time = false 
    #     B_last = (Bᴮ[:] / norm(Bᴮ))
    #     m = [0; 0; 0]
    # else
    #     B̂ᴮ = Bᴮ[:] / norm(Bᴮ)
    #     Ḃᴮ = (B̂ᴮ - B_last[:]) / dt
    #     # Ḃᴮ_unit = Ḃᴮ / norm(Bᴮ)
    #     m = -k * Ḃᴮ #/ norm(Bᴮ)
    #     # @infiltrate
    #     global B_last = B̂ᴮ
    # end
    # # Time derivative of unit vector, divided by norm of vector
    
    # # Ḃᴮ = cross(Bᴮ, ω)
    # # Ḃᴮ_unit = Ḃᴮ / norm(Bᴮ) # norm(Ḃᴮ)

    # # m = -k * Ḃᴮ_unit

    # # b̂ = Bᴮ / norm(Bᴮ)
    # # m =(-k/norm(Bᴮ)) * cross(b̂, ω)

    return 0.0
end;

