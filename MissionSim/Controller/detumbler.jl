####################################################################
#                        DETUMBLER                                 #
####################################################################

struct DETUMBLER 
    ω   # Angular velocity
    Bᴮ  # Measured magnetic field vector
    dt  # Time step
end

function __init__()
    """ Sets up python functions to later be used with PyCall for testing against the Julia implementations """

    py"""
    import numpy as np 
    # import ulab.numpy as np

    def b_cross(omega, Bb):
        k = 7e-6 
        Bb_unit = Bb / np.linalg.norm(Bb)
        m = (-k * (np.eye(3) - (Bb_unit[:,None] @ Bb_unit[None, :])))  @ omega
        return m
    """
end

function generate_command(data::DETUMBLER)
    """ Uses a B-cross controller on a magnetorquer to slow down the tumbling of the satellite """

    ω = data.ω
    
    m = b_cross(ω, data.Bᴮ)
    m_test = py"b_cross"(ω, data.Bᴮ) # Testing the Julia code against the Python
    if !(m ≈ m_test)
        @infiltrate
    end

    return m
end

function b_cross(ω, Bᴮ)
    """ Simple B-cross controller that commands a torque opposite the direction of rotatoin """
    k = 7e-6   

    B̂ᴮ = Bᴮ / norm(Bᴮ)  # Make unit 
    m = -k * (I(3) - (B̂ᴮ)*(B̂ᴮ)') * ω 

    return m
end 

function b_dot(ω, Bᴮ, dt)
    """ B-Dot controller that DOES NOT CURRENTLY WORK """
    k = 7e-6

    if first_time
        global first_time = false 
        B_last = (Bᴮ[:] / norm(Bᴮ))
        m = [0; 0; 0]
    else
        B̂ᴮ = Bᴮ[:] / norm(Bᴮ)
        Ḃᴮ = (B̂ᴮ - B_last[:]) / dt
        # Ḃᴮ_unit = Ḃᴮ / norm(Bᴮ)
        m = -k * Ḃᴮ #/ norm(Bᴮ)
        # @infiltrate
        global B_last = B̂ᴮ
    end
    # Time derivative of unit vector, divided by norm of vector
    
    # Ḃᴮ = cross(Bᴮ, ω)
    # Ḃᴮ_unit = Ḃᴮ / norm(Bᴮ) # norm(Ḃᴮ)

    # m = -k * Ḃᴮ_unit

    # b̂ = Bᴮ / norm(Bᴮ)
    # m =(-k/norm(Bᴮ)) * cross(b̂, ω)

    return m 
end

