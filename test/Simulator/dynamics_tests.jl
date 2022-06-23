# [test/Simulator/dynamics_tests.jl]

""" To Do:
  - Compare rk4 methods (esp quaternion portions)
  - to avoid scope conflicts, I am just calling CustomStructs through Simulator rn...
"""

@testset "Dynamics Tests" begin 

    @testset "  Orbit Time" begin 
        # Checks that the sat position is roughly the same after one orbit
        _μ  = 3.9860044188e14
        G = 6.6743e-11
        Mₑ = 5.972e24      # kg
        Rₑ = 6378.1363e3   # m
        Rₛ = Rₑ + 550e3     # m
 
        # Set up initial state  
        r₀ = SVector{3, Float64}(Rₛ, 0, 0)           # Initial position
        a  = 1.5 * norm(r₀)      # Semi-major axis of (elliptical) orbit (circular orbit if a = r0)
        v₀ = sqrt(_μ*( (2/norm(r₀)) - (1/a)))  * SVector{3, Float64}(0.0, 1.0, 0.0)

        Tₒ = 2 * π * sqrt((a^3) / _μ) # Orbit period, seconds 

        J  = SMatrix{3, 3, Float64, 9}(I(3)) # Doesnt matter, just need for function call
        q₀ = SVector{4, Float64}(1.0, 0, 0, 0)
        ω₀ = SVector{3, Float64}(zeros(3))
        β₀ = SVector{3, Float64}(zeros(3))

        t = Epoch(2021, 9, 1, 11, 0, 0, 0.0)
        dt = 1.0

        # Orbit for T₀ seconds, should end up in about the same sun_position
        x₀ = Simulator.STATE(r₀, v₀, q₀, ω₀, β₀)
        u  = SVector{3, Float64}(zeros(3))

        N = Int(round(Tₒ)) + 1
        X = zeros(3, N)
        X[:, 1] = x₀.r

        x = x₀
        for i = 1:N - 1
            x = rk4(J, x, u, t, dt; coef_drag = 0.0, coef_srp = 0.0, n_grav = 0, m_grav = 0, third_body = false)
            X[:, i + 1] .= x.r
        end

        @test norm(x₀.r - X[1:3, end]) < 0.001 * norm(x₀.r)
    end;

    @testset "  Check energy isnt lost" begin 

        # Set up initial state  
        Rₛ = 6378.1363e3 + 550e3  # m
        r₀ = [Rₛ; 0; 0]           # Initial position
        a  = 1.0 * Rₛ             # Semi-major axis of (elliptical) orbit (circular orbit if a = r0)
        v₀ = sqrt(SimpleOrbit._μ*( (2/Rₛ) - (1/a)))  * [0.0; 1.0; 0.0]

        m = 1.2
        J = SimpleOrbit.get_inertia(; m = m)
        J = SMatrix{3, 3, Float64, 9}(J) # Doesnt matter, just need for function call
        q₀ = [1.0; 0; 0; 0]
        ω₀ = zeros(3)
        β₀ = zeros(3)

        t = Epoch(2021, 9, 1, 11, 0, 0, 0.0)
        dt = 1.0

        x₀ = Simulator.STATE(r₀, v₀, q₀, ω₀, β₀)
        u  = SVector{3, Float64}(zeros(3))

        N = 50000 
        
        E = zeros(N)
        E[1] = SimpleOrbit.energy( vcat([x₀.r; x₀.v; x₀.q; x₀.ω]...), m, reshape([J...], 3, 3))

        x = x₀
        for i = 1:N - 1
            x = rk4(J, x, u, t, dt; coef_drag = 0.0, coef_srp = 0.0, n_grav = 0, m_grav = 0, third_body = false)
            E[i + 1] = SimpleOrbit.energy( vcat([x₀.r; x₀.v; x₀.q; x₀.ω]...), m, reshape([J...], 3, 3))
        end

        # plot(E)
        @test maximum(E) - minimum(E) < 1.0

    end;

    @testset "  Check angular velocity decreases" begin 
        if !(@isdefined generate_command)
            include("../../src/MissionSim/Controller/Controller.jl"); using .Controller
        end

        # Set up initial state  
        Rₛ = 6378.1363e3 + 900e3  # m
        r₀ = SVector{3, Float64}(Rₛ, 0, 0)           # Initial position
        a  = 1.0 * Rₛ             # Semi-major axis of (elliptical) orbit (circular orbit if a = r0)
        v₀ = sqrt(SimpleOrbit._μ*( (2/Rₛ) - (1/a)))  * SVector{3, Float64}(0.0, 0.0, 1.0)

        m  = 2.0
        J  = SimpleOrbit.get_inertia(; m = m)
        J  = SMatrix{3, 3, Float64, 9}(J) # Doesnt matter, just need for function call
        q₀ = SVector{4, Float64}(1.0, 0, 0, 0)
        ω₀ = SVector{3, Float64}(1.0, -2.0, 0.3)
        β₀ = SVector{3, Float64}(zeros(3))

        t = Epoch(2021, 3, 1, 12, 0, 0, 0.0)
        dt = 5.0

        x₀ = Simulator.STATE(r₀, v₀, q₀, ω₀, β₀)
        u  = SVector{3, Float64}(zeros(3))

        N = 2000 

        ωs = []

        x = x₀
        for i = 1:N-1 
            push!(ωs, x.ω)
            Bᴵ = IGRF13(x.r, t + ((i-1) * dt))  # Use current position and time (as epoch) to get mag field
            Bᴮ = H' * L(x.q) * R(x.q)' * H * Bᴵ
            
            
            ctrl = Controller.DETUMBLER(x.ω, SVector{3, Float64}(Bᴮ), dt)
            M = generate_command(ctrl; func = Controller.b_cross) #, κ = 5e-4)
            x = rk4(J, x, M, t, dt; coef_drag = 0.0, coef_srp = 0.0, n_grav = 0, m_grav = 0, third_body = false)
        end

        # plot( hcat(ωs...)' )
        @test norm(ωs[end]) < 0.01
    end;

    @testset "  Check eclipse" begin 
        # Did it get fixed yet? 

        # r_sun = [1.47e11, 0.0, 0.0]
        # r_sat₁ = [1e7, 0.0, 0.0]
        # r_sat₂ = -r_sat₁

        # ν₁ = eclipse_conical(r_sat₁, r_sun) # SHOULD be 1.0 if the function is fixed 
        # ν₂ = eclipse_conical(r_sat₂, r_sun) # SHOULD be 0.0

        # # ν₁ = eclipse_cylindrical(r_sat₁, r_sun)  # <- this works correct, as a Comparison

        # @warn "This testset will probably fail..."
        # @test ν₁ == 1.0
        # @test ν₂ == 0.0

    end;

    @testset "  High/low velocity" begin 
        # ||r|| should constantly grow 

        # Set up initial state  
        Rₛ = 6378.1363e3 + 550e3  # m
        r₀ = SVector{3, Float64}(Rₛ, 0, 0)           # Initial position
        a  = 1.0 * Rₛ             # Semi-major axis of (elliptical) orbit (circular orbit if a = r0)
        
        # Make initial velocity double of stable orbit
        v₀ = 2.0 * sqrt(SimpleOrbit._μ*( (2/Rₛ) - (1/a)))  * SVector{3, Float64}(0.0, 0.0, 1.0)

        m  = 0.1
        J  = SimpleOrbit.get_inertia(; m = m)
        J  = SMatrix{3, 3, Float64, 9}(J) # Doesnt matter, just need for function call
        q₀ = SVector{4, Float64}(1.0, 0, 0, 0)
        ω₀ = SVector{3, Float64}(1.1, -1.1, 0.2)
        β₀ = SVector{3, Float64}(zeros(3))

        t = Epoch(2021, 9, 1, 11, 0, 0, 0.0)
        dt = 3.0

        x₀ = Simulator.STATE(r₀, v₀, q₀, ω₀, β₀)
        u  = SVector{3, Float64}(zeros(3))

        N = 10000 
        
        r = zeros(N)
        r[1] = norm(r₀)

        x = x₀
        for i = 1:N - 1
            x = rk4(J, x, u, t, dt; coef_drag = 0.0, coef_srp = 0.0, n_grav = 0, m_grav = 0, third_body = false)
            r[i + 1] = norm(x.r)
        end

        # display(plot(r))

        @test all(r[2:end] .> r[1:end-1])


        # Make the initial velocity too low now
        N = 500        # Can't go long without crashing...
        Rₛ = (6378.1363e3 + 2000e3)  # m
        r₀ = SVector{3, Float64}(Rₛ, 0, 0)           # Initial position
        v₀ = 0.97 * sqrt(SimpleOrbit._μ*( (2/Rₛ) - (1/a)))  * SVector{3, Float64}(0.0, 0.0, 1.0)
        x₀ = Simulator.STATE(r₀, v₀, q₀, ω₀, β₀)

        r = zeros(N)
        r[1] = norm(r₀)

        x = x₀
        for i = 1:N - 1
            x = rk4(J, x, u, t, dt; coef_drag = 0.0, coef_srp = 0.0, n_grav = 0, m_grav = 0, third_body = false)
            r[i + 1] = norm(x.r)
        end

        # display(plot(r))

        @test all(r[2:end] .< r[1:end-1])
    end;

    @testset "  No Initial Velocity" begin 
        # ||r|| should constantly shrink, ||v|| grow, both along one axis

        # Set up initial state  
        Rₛ = 6378.1363e3 + 800e3  # m
        r₀ = SVector{3, Float64}(Rₛ, 0, 0)           # Initial position
        a  = 1.0 * Rₛ             # Semi-major axis of (elliptical) orbit (circular orbit if a = r0)
        
        # No initial velocity
        v₀ = SVector{3, Float64}( zeros(3) )

        m = 10.0
        J = SimpleOrbit.get_inertia(; m = m)
        J = SMatrix{3, 3, Float64, 9}(J) # Doesnt matter, just need for function call
        q₀ = SVector{4, Float64}(1.0, 0, 0, 0)
        ω₀ = SVector{3, Float64}(1.1, -1.1, 2.2)
        β₀ = SVector{3, Float64}( zeros(3))

        t = Epoch(2021, 1, 5, 2, 0, 0, 0.0)
        dt = 1.0

        x₀ = Simulator.STATE(r₀, v₀, q₀, ω₀, β₀)
        u  = SVector{3, Float64}( zeros(3))

        N = 400 
        
        r = zeros(3, N)
        v = zeros(3, N)
        r[:, 1] = r₀
        v[:, 1] = v₀

        x = x₀
        for i = 1:N - 1
            x = rk4(J, x, u, t, dt; coef_drag = 0.0, coef_srp = 0.0, n_grav = 0, m_grav = 0, third_body = false)
            r[:, i + 1] .= x.r
            v[:, i + 1] .= x.v
        end

        # Y and Z components don't change...
        @test maximum(r[2, :]) ≈ 0.0 atol = 1e-10
        @test maximum(r[3, :]) ≈ 0.0 atol = 1e-10
        @test maximum(v[2, :]) ≈ 0.0 atol = 1e-10
        @test maximum(v[3, :]) ≈ 0.0 atol = 1e-10

        # X component is decreasing for all (velocity increasing in magnitude but neg direction)
        @test all(r[1, 1:end - 1] .≥ r[1, 2:end])
        @test all(v[1, 1:end - 1] .≥ v[1, 2:end])
        
    end;
           
    @testset "  No Angular Velocity" begin
        # 0 ω means no change in q

        # Set up initial state  
        Rₛ = 6378.1363e3 + 2000e3  # m
        r₀ = SVector{3, Float64}(0.0, Rₛ, 0)           # Initial position
        a  = 1.0 * Rₛ             # Semi-major axis of (elliptical) orbit (circular orbit if a = r0)
        
        # No initial velocity
        v₀ = sqrt(SimpleOrbit._μ*( (2/Rₛ) - (1/a)))  * SVector{3, Float64}(1.0, 0.0, 0.0)

        m = 5.0
        J = SimpleOrbit.get_inertia(; m = m, l = 0.1, h = 1.0)
        J = SMatrix{3, 3, Float64, 9}(J) 
        q₀ = SVector{4, Float64}(1.0, 0, 0, 0)
        ω₀ = SVector{3, Float64}( zeros(3))
        β₀ = SVector{3, Float64}( zeros(3))

        t = Epoch(2021, 1, 5, 2, 0, 0, 0.0)
        dt = 1.0

        x₀ = Simulator.STATE(r₀, v₀, q₀, ω₀, β₀)
        u  = SVector{3, Float64}( zeros(3))

        N = 1000 
                
        qs = zeros(4, N)
        ωs = zeros(3, N)
        qs[:, 1] = q₀
        ωs[:, 1] = ω₀

        x = x₀
        for i = 1:N - 1
            x = rk4(J, x, u, t, dt; coef_drag = 0.0, coef_srp = 0.0, n_grav = 0, m_grav = 0, third_body = false)
            qs[:, i + 1] .= x.q
            ωs[:, i + 1] .= x.ω    
        end

        @test all(qs[:, 1:end - 1] .≈ qs[:, 2:end])
        @test all(ωs .≈ 0.0)
    end;

    @testset "  Verify quaternion is unit" begin 
        # @warn "Not implemented!"

        t = Epoch(2021, 9, 1, 11, 0, 0, 0.0)
        dt = 3.0
        x = Simulator.STATE()
        J = Simulator.SATELLITE().J
        u  = SVector{3, Float64}( zeros(3) )

        @test norm(x.q) ≈ 1.0 
        for i = 1:20
            x = rk4(J, x, u, t + (i - 1) * dt, dt)
            @test norm(x.q) ≈ 1.0
        end

    end;

    @testset "  Verify one full rotation" begin 
        # Does my RK4 rotate in a circle correctly?
        q₀ = SVector{4, Float64}(1, 0, 0, 0)
        β  = SVector{3, Float64}(0, 0, 0)
        ω  = SVector{3, Float64}(0.0, deg2rad(36), 0.0)
        J  = SMatrix{3, 3, Float64, 9}([0.2 0 0; 0 0.2 0; 0 0 0.2])
        u  = SVector{3, Float64}(0.0, 0, 0)

        x = STATE(; q = q₀, β = β, ω = ω)
        t = Epoch(2020, 1, 1)
        dt = 0.1
        N = Int(10 / dt)
        for i = 1:N
            x = Simulator.rk4(J, x, u, t, dt; σβ = 0.0)  # No bias
        end

        @test norm(cayley_map(x.q, q₀)) ≈ 0.0 atol = 1e-6
    end;

    # @testset "  Compare rk4 methods" begin
    #     @warn "Not implemented!"
    # end;
end; 
