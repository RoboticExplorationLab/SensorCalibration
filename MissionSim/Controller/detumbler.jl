####################################################################
#                        DETUMBLER                                 #
####################################################################

println("Make DTUMBLER struct immutable again")
mutable struct DETUMBLER 
    ω 
    Bᴮ
    dt
end

B_last = [0 0 0]
first_time = true
println("Ḃ Detumbler still doesn't work, k is fixed, and could probably account for bias...")

function __init__()
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

function generate_command(data::DETUMBLER, state)

    ω = data.ω
    
    m = b_cross(ω, data.Bᴮ)
    m_test = py"b_cross"(ω, data.Bᴮ)
    if !(m ≈ m_test)
        @infiltrate
    end
    # m = b_dot(ω, data.Bᴮ, data.dt)
    # println("\r$m  \t | $m_test")

    return m
end

# DETERMINE OPTIMAL k !!!
function b_cross(ω, Bᴮ)
    k = 7e-6   #### NEED TO DO BETTER

    B̂ᴮ = Bᴮ / norm(Bᴮ)  # Make unit 
    m = -k * (I(3) - (B̂ᴮ)*(B̂ᴮ)') * ω #[:]

    return m
end 

function b_dot(ω, Bᴮ, dt)
    k = 7e-6

    # if first_time
    #     global first_time = false 
    #     B_last = Bᴮ 
    #     m = [0; 0; 0]
    # else 
    #     Ḃᴮ = (B_last - Bᴮ) / dt 
    #     m = -k * Ḃᴮ
    #     @infiltrate
    #     global B_last = Bᴮ 
    # end

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

