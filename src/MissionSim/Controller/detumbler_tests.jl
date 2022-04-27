# [test/Controller/detumbler_tests.jl]

include("detumbler.jl")

@testset "\tDetumbler Tests" begin

    @testset "(P) Struct and Interface" begin 
        """ Make sure the provided function calls work """

        ω  = [1.0, 2.0, 3.0]
        Bᴮ = [1.0, 1.5, 1.0]
        dt = 0.2

        @info "\tShould throw a warning..."
        ctrl = DETUMBLER(ω, Bᴮ, dt)

        ω  = SVector{3, Float64}([1.0, 2.0, 3.0])
        Bᴮ = SVector{3, Float64}([1.0, 1.5, 1.0])
        dt = 0.1 

        ctrl = DETUMBLER(ω, Bᴮ, dt)

        M = b_cross(ω, Bᴮ)

        # M = b_dot(ω, Bᴮ, dt)

        M = generate_command(ctrl)

        numAllocs = @ballocated generate_command($ctrl; κ = 1e-8)

        @test numAllocs == 0
    end

    # Adjusted value for κ... should be adjusted to account for Bᴮ
    @testset "(1) M ⟂ B̂" begin 
        """ M ⟂ b̂ (M = u for this case) """
        
        ω  = SVector{3, Float64}([1.0, 2.0, 3.0])
        Bᴮ = SVector{3, Float64}([1.0, 1.5, 1.0])
        dt = 0.1

        ctrl = DETUMBLER(ω, Bᴮ, dt)
        M = generate_command(ctrl; func = b_cross, κ = 1e-8)
        
        B̂ᴮ = Bᴮ / norm(Bᴮ) 
        @test norm(cross(M, B̂ᴮ)) ≈ 0.0 atol = 1e-7

        for i = 1:9
            ω  = SVector{3, Float64}(rand(3))
            Bᴮ = SVector{3, Float64}(rand(3))
            M = generate_command(ctrl; κ = 1e-8)
            
            B̂ᴮ = Bᴮ / norm(Bᴮ) 
            @test norm(cross(M, B̂ᴮ)) ≈ 0.0 atol = 1e-7
        end
    end

    @testset "(2) Comparison to b dot" begin 
        """ ω × b̂ ∝ b dot """

    end

    @testset "(3) System Energy" begin 
        """ Verify that the net energy in the system decreases """
    end

    @testset "(4) Satellite Detumbling " begin 
        """ Verify that we can detmble a satellite """

    end

    @testset "(5) Bonus tests " begin 
        """ (verify it 𝑑𝑜𝑒𝑠𝑛𝑡 work when parallel?) """

        # """ Test/compare gain κ? """

        # """ Test allocations are low? (or 0 or 3?)

    end

end;

println("Done")