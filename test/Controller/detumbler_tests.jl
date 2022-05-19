# [test/Controller/detumbler_tests.jl]

using Test, BenchmarkTools
include("../src/MissionSim/Controller/detumbler.jl")

@testset "\tDetumbler Tests" begin

    @testset "(P) Struct and Interface" begin 
        """ Make sure the provided function calls work """

        ω  = [1.0, 2.0, 3.0]
        Bᴮ = [1.0, 1.5, 1.0]
        dt = 0.2

        # @info "\tShould throw a warning..."
        ctrl = DETUMBLER(ω, Bᴮ, dt)

        ω  = SVector{3, Float64}([1.0, 2.0, 3.0])
        Bᴮ = SVector{3, Float64}([1.0, 1.5, 1.0])
        dt = 0.1 

        ctrl = DETUMBLER(ω, Bᴮ, dt)

        M = b_cross(ω, Bᴮ)

        # M = b_dot(ω, Bᴮ, dt)

        M = generate_command(ctrl)

        numAllocs = @ballocated generate_command($ctrl; κ = 0.01)

        @test numAllocs == 0
    end

    @testset "(1) Command is perpendiclar to mag field " begin 
        """ M ⟂ b̂ (M = u for this case), so dot(M, B) ≈ 0.0 """
        
        ω  = SVector{3, Float64}([1.0, 2.0, 3.0])
        Bᴮ = SVector{3, Float64}([1.0, 1.5, 1.0])
        dt = 0.1

        ctrl = DETUMBLER(ω, Bᴮ, dt)
        M = generate_command(ctrl; func = b_cross, κ = 1e-2)  # Should work for any κ

        M̂ = M / norm(M)
        B̂ᴮ = Bᴮ / norm(Bᴮ) 

        @test norm(dot(M̂, B̂ᴮ)) ≈ 0.0 atol = 1e-10

        for i = 1:9
            ω  = SVector{3, Float64}(rand(3))
            Bᴮ = SVector{3, Float64}(rand(3))
            ctrl = DETUMBLER(ω, Bᴮ, dt)
            M = generate_command(ctrl) 

            M̂ = M / norm(M)
            
            B̂ᴮ = Bᴮ / norm(Bᴮ) 
            @test norm(dot(M̂, B̂ᴮ)) ≈ 0.0 atol = 1e-10
        end
    end

    # EMPTY - B dot not working
    @testset "(2) Comparison to b dot" begin 
        """ ω × b̂ ∝ b dot """
    end

    include("../../../test/SimpleOrbit.jl")
    @testset "(3) System Energy" begin 
        """ Verify that the net energy in the system decreases """

        # Set up a satellite, verify energy decreases
        _q0 = [1, 0, 0, 0]
        _ω0 = [2.2, 1.0, -1.0]
        J = get_inertia() 

        x0 = [_q0; _ω0]

        N = 15000
        Eₓ = zeros(N)
        Eₓᵤ = zeros(N)
        x  = zeros(7, N); x[:,  1] = x0
        xᵤ = zeros(7, N); xᵤ[:, 1] = x0
        Eₓ[1]  = energy(x0, 0.0, J)
        Eₓᵤ[1] = energy(x0, 0.0, J)
        dt = 5.0   # seconds
        for i = 1:N-1
            x[:, i+1]  = rk4(dynamics, x[:, i],  zeros(3), dt, J) # Shouldn't change

            ω  = xᵤ[5:7, i]
            Bᴮ = [1.0 * sin(i / 1000), 2.0 * cos(i / 1000), -1.0 * sin(i / 5000)]

            ctrl = DETUMBLER(ω, Bᴮ, dt)
            M = generate_command(ctrl; func = b_cross) 
            xᵤ[:, i + 1] = rk4(dynamics, xᵤ[:, i], M, dt, J)

            Eₓ[i + 1] = energy(x[:,  i + 1], 0.0, J)
            Eₓᵤ[i+ 1] = energy(xᵤ[:, i + 1], 0.0, J)
        end

        if false
            display(plot(x[5:7, :]', title = "Unforced"))
            display(plot(xᵤ[5:7, :]', title = "Forced"))

            plot(Eₓ, label = "Unforced")
            display(plot!(Eₓᵤ, label = "Forced"))
        end

        # Verify that the energy decreases 
        @test Eₓᵤ[end] < 0.1
        @test minimum(Eₓ - Eₓᵤ) ≥ 0.0
    end

    @testset "(4) Full State test" begin 
        using SatelliteDynamics
        include("../mag_field.jl")

        ## Generate an orbit
        _r0 = (550+6371)*(10^(3))   # Distance from two center of masses 
        _a  = 1.2 * _r0             # Semi-major axis of (elliptical) orbit (circular orbit if a = r0)
        _v0 = sqrt(_μ*( (2/_r0) - (1/_a))) 
        _q0 = [1; 0; 0; 0]
        _ω0 = [0.2; -1.0; 0.5]

        _m  = 1.2                   # sat mass (Kg) 

        _r0 = [_r0; 0.0; 0.0]
        _v0 = [0.0; _v0; 0.0]

        x0 = [_r0;  _q0;  _v0;  _ω0]
        J = get_inertia(; m = _m)

        ## Set up simulation
        N = 1000
        dt = 5.0   # seconds
        _epc = Epoch(2021, 9, 1, 11, 0, 0, 0.0); # Initial time for sim

        x  = zeros(13, N); x[:,  1] = x0
        xᵤ = zeros(13, N); xᵤ[:, 1] = x0

        # Set up energy vectors
        Eₓ  = zeros(N)
        Eₓᵤ = zeros(N)
        Eₓ[1]  = energy(x0, _m, J)
        Eₓᵤ[1] = energy(x0, _m, J)

        Bᴮs = []

        for i = 1:N-1
            # Reference, uncontrolled orbit
            x[:, i+1]  = rk4(dynamics, x[:, i],  zeros(3), dt, J) # No control input, spin shouldn't change
            Eₓ[i + 1] = energy(x[:,  i + 1], _m, J)

            # Controlled orbit 
            ω  = xᵤ[11:13, i]
            Bᴵ = IGRF13(xᵤ[1:3, i], _epc + ((i-1) * dt))  # Use current position and time (as epoch) to get mag field
            Bᴮ = H' * L(xᵤ[4:7, i]) * R(xᵤ[4:7, i])' * H * Bᴵ

            push!(Bᴮs, Bᴮ)

            ctrl = DETUMBLER(ω, Bᴮ, dt)
            M = generate_command(ctrl; func = b_cross, κ = 0.01)
            xᵤ[:, i + 1] = rk4(dynamics, xᵤ[:, i], M, dt, J)
            Eₓᵤ[i+ 1] = energy(xᵤ[:, i + 1], _m, J)
        end

        if false
            display(plot( hcat(Bᴮs...)', title = "Mag Field"))
            display(plot(x[11:13, :]', title = "Unforced"))
            display(plot(xᵤ[1:3, :]', title = "r (Forced)"))
            display(plot(xᵤ[4:7, :]', title = "q (Forced)"))
            display(plot(xᵤ[8:10, :]', title = "V (Forced)"))
            display(plot(xᵤ[11:13, :]', title = "ω (Forced)"))

            plot(Eₓ, label = "Unforced")
            display(plot!(Eₓᵤ, label = "Forced"))
        end

        # Verify that the energy decreases 
        @test 0.1 * norm(xᵤ[11:13, 1]) > norm(xᵤ[11:13, end])
        @test minimum(Eₓ - Eₓᵤ) ≥ 0.0
    end

    # EMPTY 
    @testset "(5) Bonus tests " begin 
        """ (verify it 𝑑𝑜𝑒𝑠𝑛𝑡 work when parallel?) """

        # """ Test/compare gain κ? """
    end
end;



