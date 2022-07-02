# [test/Estimator/magnetometer_calibration_tests.jl]

""" To Do:
 - 
"""
######## Stuff needed to run this independently #############################
# using Test, BenchmarkTools
# using StaticArrays, Plots, LinearAlgebra, Distributions, ForwardDiff

# include("../../src/CustomStructs.jl");  using .CustomStructs 
# include("../../src/Estimator/magnetometer_calibration.jl")
#############################################################################



@testset "Mag Calib Tests" begin
    import .Estimator.extract_elements, .Estimator.vec_to_matrix_bias, .Estimator.update!
    import .Estimator.residual, .Estimator.gauss_newton, .Estimator.estimate

    function make_T(; a = rand(Normal(1.0, 0.1)), b = rand(Normal(1.0, 0.1)), c = rand(Normal(1.0, 0.1)),
        ρ = rand(Normal(0.0, deg2rad(3.0))), λ = rand(Normal(0.0, deg2rad(3.0))), ϕ = rand(Normal(0.0, deg2rad(3.0))) )

        T = [a          0.0              0.0;
        b*sin(ρ)   b*cos(ρ)         0.0;
        c*sin(λ)   c*cos(λ)*sin(ϕ)  c*cos(λ)*cos(ϕ) ]

        return T
    end

    @testset "Extract Elements" begin 

        # Should be the exact same...
        function extract_parameters_old(T)
            """ Extracts the calibration parameters from a matrix 
                    [a, b, c] (scale factors) and [ρ, λ, ϕ] (non-orthogonality angles) """
        
            a = T[1,1] # Easy 
        
            b = sqrt((T[2,1]^2) + (T[2,2]^2)) # (bsin)^2 + (bcos)^2 = b^2
            ρ = atan(T[2,1], T[2,2]) # sin/cos to maintain signs
        
            c = sqrt((T[3,1]^2) + (T[3,2]^2) + (T[3,3]^2))
            ϕ = atan(T[3,2] / T[3,3])
            λ = atan(  sign(T[3,1]) * sqrt( (T[3,1]^2) ),  
                        sign((T[3,2]^2) + (T[3,3]^2)) * sqrt( (T[3,2]^2) + (T[3,3]^2) ) ) # Not positve this portion of the sign is actually useful
        
            return a, b, c, ρ, λ, ϕ
        end

        function random_params(; a = rand(Normal(1.0, 0.1)), b = rand(Normal(1.0, 0.1)), c = rand(Normal(1.0, 0.1)),
            ρ = rand(Normal(0.0, deg2rad(3.0))), λ = rand(Normal(0.0, deg2rad(3.0))), ϕ = rand(Normal(0.0, deg2rad(3.0))),
            βx = randn(), βy = randn(), βz = randn())
            return a, b, c, ρ, λ, ϕ, βx, βy, βz
        end

        T = make_T()
        β = rand(3)

        # @code_warntype extract_elements(T, rand(3))
        # @btime extract_elements($T, $β);

        @testset "Comparison to old" begin
            N = 1000
            t = zeros(N)
            for i = 1:N
                T = make_T()
                els  = extract_elements(T, randn(3))
                els₂ = extract_parameters_old(T)

                t[i] = (els[1:6] == els₂)
            end

            @test sum(t) == N
        end

        # Ensure a, b, c are positive; ρ, λ, ϕ are small mag 
        for i = 1:10
            T, β = make_T(), randn(3)
            els  = extract_elements(T, β)
            @test all(els[1:3] .≥ 0.0)
            @test all(abs.(els[4:6]) .< pi/3)
        end

        # No angle error 
        T = make_T(ρ = 0.0, λ = 0.0, ϕ = 0.0)
        β = randn(3)
        els = extract_elements(T, β)
        @test all(els[4:6] .≈ 0.0) 

        # No scale offset
        T = make_T(a = 1.0, b = 1.0, c = 1.0)
        β = randn(3)
        els = extract_elements(T, β)
        @test all(els[1:3] .≈ 1.0) 

        # More 
        T = make_T(ρ = deg2rad(45), λ = deg2rad(45), ϕ = deg2rad(45))
        els = extract_elements(T, randn(3))
        @test all(els[4:6] .≈ deg2rad(45)) 

        T = make_T(ρ = deg2rad(45), λ = deg2rad(45), ϕ = deg2rad(-45))
        els = extract_elements(T, randn(3))
        @test all(els[4:5] .≈ deg2rad(45)) 
        @test (els[6] ≈ deg2rad(-45)) 
 
        T = make_T(a = 1.0, b = 0.9, c = 1.1, ρ = 0.0, λ = deg2rad(30), ϕ = deg2rad(-30))
        els = extract_elements(T, randn(3))
        @test [els[1:6]...] ≈ [1.0, 0.9, 1.1, 0.0, deg2rad(30), deg2rad(-30)]
    
    
        N = 1000
        t = zeros(N)
        for i = 1:N
            p = random_params();
            T = make_T(a = p[1], b = p[2], c = p[3], ρ = p[4], λ = p[5], ϕ = p[6]);
            β = p[7:9];
            els  = extract_elements(T, β);
            t[i] = ([p...] ≈ [els...])
        end
        @test sum(t) == N


        # Handwritten one 
        T = [1.05      0.0        0.0;
             0.238114  0.888652   0.0;
            -0.19101   0.113234   1.07735]

        els = extract_elements(T, randn(3))
        @test [els[1:6]...] ≈ [1.05, 0.92, 1.1, deg2rad(15), deg2rad(-10), deg2rad(6)] atol = 1e-3
    end;

    @testset "Vec to Matrix" begin 

        function parameters_to_matrix_bias_old(p)
            """ Reshapes Gauss-Newton output into a calibration matrix and biases """
        
            T = [p[1]   0       0;
                    p[2]   p[4]    0;
                    p[3]   p[5]    p[6]];      # Calibration matrix
        
            β = p[7:9];     # Bias vector 
        
        
            # Account for sign ambiguities in the calibration matrix
            if T[3,3] < 0
                T[3,3] = -T[3,3]
            end
            if T[2,2] < 0 
                T[2,2] = -T[2,2]
                T[3,2] = -T[3,2]
            end
            if T[1,1] < 0
                T[:, 1] = -T[:, 1]
            end
        
            return T, β
        end

        # @code_warntype Estimator.vec_to_matrix_bias(rand(9))
        p = rand(9)
        # @btime Estimator.vec_to_matrix_bias($p);

        @testset "Comparison to old" begin
            N = 1000
            t = zeros(N)
            for i = 1:N
                p = rand(9)
                p1 = vec_to_matrix_bias(p)
                p2 = parameters_to_matrix_bias_old(p)

                t[i] = (p1 == p2)
            end

            @test sum(t) == N
        end

        # Verify T is upper triangular 
        N = 100
        t = zeros(N)
        for i = 1:N
            T, β = vec_to_matrix_bias(rand(9))
            t[i] = (T[1, 2] == 0.0 && T[1, 3] == 0.0 && T[2, 3] == 0.0) 
        end
        @test sum(t) == N

        # Verify diagonal terms are non-zero
        N = 100
        t = zeros(N)
        for i = 1:N
            T, β = vec_to_matrix_bias(rand(9))
            t[i] = (T[1, 1] > 0.0 && T[2, 2] > 0.0 && T[3, 3] > 0.0) 
        end
        @test sum(t) == N

    end;

    @testset "Constructor" begin 

        N = 5
        Bmeas = [1.1, 2, 3]
        Bpred = [1, 2.1, 2.9]
        mag_cal = MAG_CALIBRATOR(N, Bmeas, Bpred)

        # @btime update!($mag_cal, $Bmeas, $Bpred) evals = 10 samples = 10;

        # Test simple constructor 
        @test size(mag_cal.B_meas) == (15,)
        @test size(mag_cal.B_pred) == (15,)
        @test size(mag_cal.A)      == (15 , 9)
        @test mag_cal.idx[1]       == 1
        @test mag_cal.N            == N 

        @test mag_cal.B_meas[1:3] ==  [1.1, 2, 3]
        @test mag_cal.B_meas[4:15] == zeros(12)
        @test mag_cal.B_pred[1:3] == [1, 2.1, 2.9]
        @test mag_cal.B_pred[4:15] == zeros(12)
        @test mag_cal.A[1:3, 7:9] == I(3)
        @test (mag_cal.A[1, 4] == 0.0) && (mag_cal.A[2, 6] == 0.0)

        # and the update stuff too 
        Bmeas = [3, 2, 1.1]
        Bpred = [-1.0, 2.0, 4.4]
        update!(mag_cal, Bmeas, Bpred)

        @test size(mag_cal.B_meas) == (15,)
        @test size(mag_cal.B_pred) == (15,)
        @test size(mag_cal.A)      == (15 , 9)
        @test mag_cal.idx[1]       == 2
        @test mag_cal.N            == N 

        @test mag_cal.B_meas[1:6] ==  [1.1, 2, 3, 3, 2, 1.1]
        @test mag_cal.B_meas[7:15] == zeros(9)
        @test mag_cal.B_pred[1:6] == [1, 2.1, 2.9, -1.0, 2.0, 4.4]
        @test mag_cal.B_pred[7:15] == zeros(9)
        @test mag_cal.A[4:6, 7:9] == I(3)
        @test (mag_cal.A[4, 4] == 0.0) && (mag_cal.A[5, 6] == 0.0)
        @test (mag_cal.A[4, 1] == -1.0) && (mag_cal.A[5, 4] == 2.0)

        update!(mag_cal, Bmeas, Bpred)
        update!(mag_cal, Bmeas, Bpred)
        update!(mag_cal, Bmeas, Bpred)
        @test mag_cal.idx[1]       == N
        @test sum(mag_cal.A[:, 7:9]) == 3 * N

        # Test that it provides a warning and does not keep filling once full 
        @test_logs (:warn, "Mag Calibrator is already full!") update!(mag_cal, Bmeas, Bpred)
    end;

    @testset "Residual" begin 

        function generate_data(; N = 50, T = make_T(), β = randn(3))

            bp = randn(3)
            bm = T * bp + β
            data = MAG_CALIBRATOR(N, bm, bp)

            for i = 2:N
                bp = randn(3)
                bm = T * bp + β
                update!(data, bm, bp)
            end

            return T, β, data, N
        end

        # Verify it runs
        T, β, data, N = generate_data()
        init_guess = ([T[T .!= 0.0]; β]) + 0.01 * randn(9) # Perturb to get initial guess
        r = residual(init_guess, data)
        @test size(r) == (N,)


        # Verify a perfect guess has all zero residual 
        T, β, data, N = generate_data(T = I(3), β = zeros(3))
        T = Matrix(T)
        init_guess = ([1.0; 0.0; 0.0; 1.0; 0.0; 1.0; β])   # Perfect initial guess
        r = residual(init_guess, data)
        @test all(r .== 0.0)


        T, β, data, N = generate_data()
        init_guess = ([T[T .!= 0.0]; β])   # Perfect initial guess
        r = residual(init_guess, data)
        @test isapprox(r, zeros(N), atol = 1e-12)

        Nₜ = 100
        t = zeros(Nₜ)
        for i = 1:Nₜ
            T, β, data, N = generate_data()
            init_guess = ([T[T .!= 0.0]; β])   # Perfect initial guess
            r = residual(init_guess, data)
            t[i] = isapprox(r, zeros(N), atol = 1e-12)
        end
        @test sum(t) == Nₜ


        # Verify that it is the same if T is negative (bc we are looking at magnitude)
        T, β, data, N = generate_data(T = -I(3), β = zeros(3))
        T = Matrix(T)
        init_guess = ([1.0; 0.0; 0.0; 1.0; 0.0; 1.0; β])   # Perfect initial guess
        r = residual(init_guess, data)
        @test all(r .== 0.0)

        T, β, data, N = generate_data(; β = zeros(3))
        init_guess = ([T[T .!= 0.0]; β])   # Perfect initial guess
        data.B_meas .= -data.B_meas
        r = residual(init_guess, data)
        @test isapprox(r, zeros(N), atol = 1e-12)

        # Verify it is non-zero for a noisy guess
        T, β, data, N = generate_data(; N = 1000)
        init_guess = data.A \ data.B_meas

        init_guess =  [T[1, 1]; T[2, 1]; T[2, 2]; T[3, 1]; T[3, 2]; T[3, 3]; β[:]] + 0.1 * randn(9)
        r = residual(init_guess, data)
        @test norm(r) ≉ 0.0


        # @code_warntype residual(init_guess, data)
        # @btime residual($init_guess, $data);
    end;
   
    @testset "Gauss Newton" begin 

        function generate_data(; N = 200, T = make_T(), β = randn(3), σ = 0.001)

            bp = randn(3)
            bm = T * bp + β
            data = MAG_CALIBRATOR(N, bm, bp)

            for i = 2:N
                bp = randn(3)
                bm = T * bp + β + σ * randn(3)
                update!(data, bm, bp)
            end

            return T, β, data, N
        end

        function matrix_bias_to_parameters(T, β)
            return [T[1, 1]; T[2, 1]; T[3, 1]; T[2, 2]; T[3, 2]; T[3, 3]; β[:]]
        end

        # Verify function call works and returns correct thing 
        T, β, data, N = generate_data(N = 50)
        init_guess = data.A \ data.B_meas
        final_est  = gauss_newton(init_guess, data)
        @test size(final_est) == (9,)

        # Verify that the solution ≈ answer 
        for i = 1:10
            T, β, data, N = generate_data(N = 1000, β = -2 * randn(3))
            init_guess = data.A \ data.B_meas
            final_est  = gauss_newton(init_guess, data)
            @test final_est ≈ matrix_bias_to_parameters(T, β) atol = 1e-3
        end


        # Verify residual decreases...
        for i = 1:1000
            T, β, data, N = generate_data(N = 100, β = -2 * randn(3))
            init_guess = data.A \ data.B_meas
            final_est  = gauss_newton(init_guess, data)

            ri = residual(init_guess, data)
            rf = residual(final_est, data)
            @test norm(rf) ≤ norm(ri)
        end
    end;

    @testset "Estimate" begin

        function generate_data(; N = 200, T = make_T(), β = randn(3), σ = 0.005)

            bp = randn(3)
            bm = T * bp + β
            data = MAG_CALIBRATOR(N, bm, bp)

            for i = 2:N
                bp = randn(3)
                bm = T * bp + β + σ * randn(3)
                update!(data, bm, bp)
            end

            return T, β, data, N
        end

        # Verify function call, sat updates, updates correctly (no noise)
        sat = SATELLITE();
        T, β, data, N = generate_data(N = 1000, σ = 0.0);
        sat_new = estimate(sat, data);
        @test sat_new != sat 
        @test sat_new.magnetometer.bias ≈ β atol = 1e-12
        a, b, c, ρ, λ, ϕ, _, _, _ = extract_elements(T, β)
        @test sat_new.magnetometer.scale_factors ≈ [a, b, c] atol = 1e-12
        @test sat_new.magnetometer.non_ortho_angles ≈ [ρ, λ, ϕ] atol = 1e-12
    
        # With noise
        for i = 1:10
            sat = SATELLITE();
            T, β, data, N = generate_data(N = 5000);
            sat_new = estimate(sat, data);
            @test sat_new != sat 
            @test sat_new.magnetometer.bias ≈ β atol = 1e-3
            a, b, c, ρ, λ, ϕ, _, _, _ = extract_elements(T, β)
            @test sat_new.magnetometer.scale_factors ≈ [a, b, c] atol = 1e-3
            @test sat_new.magnetometer.non_ortho_angles ≈ [ρ, λ, ϕ] atol = 1e-3
        end
    end;

end;

