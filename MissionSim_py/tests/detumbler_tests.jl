"""
    TEST FILE for comparing detumbler.jl with detumbler.py 
"""

using Infiltrator, PyCall, Test
using LinearAlgebra 


include("detumbler.jl")

function __init_python__()
    py"""
    import sys 
    sys.path.insert(0, "./") # Allows us to include nearby python scripts

    import numpy as np 

    from detumbler import Detumbler
    """
end

### SPECIFIC FUNCTIONS ###

__init_python__()

@testset "Specific Functions" begin

    # DETUMBLER 
    dt = 1.0
    # J = I(3)
    for i = 1:20
        ω = 2 * randn(3)   # [0; -0.1; 0]
        Bᴮ = 5 * randn(3)  # [1; 0; 0]

        data = DETUMBLER(ω, Bᴮ, dt)
        c_jl = generate_command(data)

        obj_py = py"Detumbler"()
        c_py = obj_py.generate_command(ω, Bᴮ)

        @test c_jl ≈ c_py 
    end


    # # DETERMINE OPTIMAL GAIN
    # J = I(3) * 0.001666666666667 
    # inc = 51.6426                           # Degrees
    # ecc = 0.0001717
    # sma = (6378136.3 + 421e3) / (1 + ecc)   # km

    # T 

end


function dynamics(x, u, t)
    w = x 
    ẇ = (J^(-1)) * (u - cross(w, J*w))
    return ẇ
end

function rk4(x_n,u,t_n,h)
    """ 
        Modified rk4 function (from Kevin)
    """

    k1 = h * dynamics(x_n, u, t_n)
    k2 = h * dynamics(x_n + k1/2, u, t_n + h/2)
    k3 = h * dynamics(x_n + k2/2, u, t_n + h/2)
    k4 = h * dynamics(x_n + k3, u, t_n + h)

    return (x_n + (1/6)*(k1+2*k2+2*k3 + k4))
end

const J = I(3)

@testset "Sequence" begin
    N = 5000
    obj_py = py"Detumbler"()

    w₀ = [0.2; -0.5; 0.1]
    B₀ = [0.33; 0.33; 0.33]
    state = w₀
    dt = 1.0

    for i = 1:N 
        t = (i-1) * dt

        B = [B₀[1]; B₀[2]*cos(deg2rad(t)); B₀[3]*sin(deg2rad(t)) ]
        
        data = DETUMBLER(state, B, dt)
        u = generate_command(data)

        u_py = obj_py.generate_command(state, B)

        @test u ≈ u_py

        state = rk4(state, u, t, dt)  
    end

end


