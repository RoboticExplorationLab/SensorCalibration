####################################################################
#                        DETUMBLER                                 #
####################################################################

struct DETUMBLER 
    ω 
    Bᴮ
    dt
end

B_last = [0 0 0]
first_time = true
# IGNORES NOISE AND BIAS RN!!

function generate_command(data::DETUMBLER, state)

    ω = data.ω
    # ω = state[11:13]
    
    m = b_cross(ω, data.Bᴮ)
    # m = b_dot(ω, data.Bᴮ, data.dt)


    if (norm(ω) < 0.01)
        finished = true 
    else 
        finished = false
    end

    return m, finished
end

# DETERMINE OPTIMAL k !!!
function b_cross(ω, Bᴮ)
    k = 7e-6   #### NEED TO DO BETTER

    B̂ᴮ = Bᴮ / norm(Bᴮ)  # Make unit 
    m = -k * (I(3) - (B̂ᴮ)*(B̂ᴮ)') * ω #[:]

    return m
end 

function b_dot(ω, Bᴮ, dt)
    k = 0.01 #7e-6

    if first_time
        global first_time = false 
        B_last = Bᴮ 
        m = [0; 0; 0]
    else 
        Ḃᴮ = B_last - Bᴮ 
        global B_last = Bᴮ 
        m = -Ḃᴮ
    end

    # if first_time
    #     global first_time = false 
    #     B_last = (Bᴮ[:] / norm(Bᴮ))
    #     m = [0; 0; 0]
    # else
    #     B̂ᴮ = Bᴮ[:] / norm(Bᴮ)
    #     Ḃᴮ = (B̂ᴮ - B_last[:]) / dt
    #     global B_last = B̂ᴮ
    #     # Ḃᴮ_unit = Ḃᴮ / norm(Bᴮ)
    #     m = -k * Ḃᴮ / norm(Bᴮ)
    # end




    # Time derivative of unit vector, divided by norm of vector
    
    # Ḃᴮ = cross(Bᴮ, ω)
    # Ḃᴮ_unit = Ḃᴮ / norm(Bᴮ) # norm(Ḃᴮ)

    # m = -k * Ḃᴮ_unit




    # b̂ = Bᴮ / norm(Bᴮ)
    # m =(-k/norm(Bᴮ)) * cross(b̂, ω)

    return m 
end

