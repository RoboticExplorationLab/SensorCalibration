"""
    TEST FILE for comparing magnetometer_calibration.jl with 
        magnetometerCalibrator.py to verify the Python implementation

    
"""

using Infiltrator, PyCall, Test
using LinearAlgebra 
using StaticArrays

include("../CustomStructs.jl");         using .CustomStructs 
include("magnetometer_calibration.jl")


function __init__()
    py"""
    import sys 
    sys.path.insert(0, "./")

    import numpy as np 
    from magnetometerCalibrator import * #MagnetometerCalibrator, sign
    """
end

__init__()

function generate_random_parameters()
    invalid = true
    while invalid
        a, b, c = 0.5 .+ rand(3)
        rho, lambda, phi = deg2rad.(rand(-15:15, 3))
        
        bx, by, bz = randn(3)

        if (a <= 0 || b <= 0 || c <= 0) || (rho == 0 || lambda == 0 || phi == 0)
            invalid = true 
        else 
            return [a, b, c, rho, lambda, phi, bx, by, bz]
        end
    end
end


########################
# INDIVIDUAL FUNCTIONS #
########################


if false
    @testset "Sign" begin
        for i = 1:1000
            p = 100 * randn()
            @test sign(p) == py"sign"(p)
        end
    end

    @testset "Parameters to Matrix" begin
        obj_py = py"MagnetometerCalibrator"()
        p = [1, 1, 1, 0, 0, 0, 0, 0, 0]
        @test parameters_to_matrix_bias(p) == obj_py.parameters_to_matrix_bias(p)

        for i = 1:500
            p = generate_random_parameters()
            T,  B  = parameters_to_matrix_bias(p)
            T2, B2 = obj_py.parameters_to_matrix_bias(p)

            @test T ≈ T2
            @test B ≈ B2 
        end
    end

    @testset "Matrix to Parameters" begin
        obj_py = py"MagnetometerCalibrator"()
        T = [1 0 0; 0 1 0; 0 0 1]
        @test extract_parameters(T) == obj_py.matrix_to_parameters(T)
        @test collect(extract_parameters(T)) == [1, 1, 1, 0.0, 0, 0]  # Collect converts tuple -> array

        for i = 1:500
            p = generate_random_parameters()
            T, B = parameters_to_matrix_bias(p)

            @test collect(extract_parameters(T)) ≈ collect(obj_py.matrix_to_parameters(T))
        end
    end

    @testset "F" begin
        obj_py = py"MagnetometerCalibrator"()

        for i = 1:500
            B_meas = 10 .* randn(3)
            p = generate_random_parameters()
            obj_py.f(B_meas, p)

            @test f(B_meas, p) ≈ obj_py.f(B_meas, p)
        end
    end

end



########################
#     FULL PROCESS     #
########################

@testset "Full Procedure" begin

    for i = 1:10
        dummy_sat = SATELLITE(randn(3,3), MAGNETOMETER(zeros(3), zeros(3), zeros(3)),
                            DIODES(zeros(6), zeros(6), zeros(6)), zeros(7), zeros(6,6)) 
        obj_py = py"MagnetometerCalibrator"(93)
        data = initialize(MAG_CALIB(0, 0))

        count = 1
        sat, dict = 0, 0
        data = 0
        p = generate_random_parameters()
        T,  B  = parameters_to_matrix_bias(p)
        while !obj_py.has_run
            print("\r$count")
            pred = 10 .* randn(3)
            meas = T * pred + B
            data = MAG_CALIB(meas, pred)

            sat, data = estimate_vals(dummy_sat, data, count == (93))
            dict = obj_py.update_estimate(meas, pred)

            count += 1
        end

        @test sat.magnetometer.scale_factors    ≈ [obj_py.a, obj_py.b, obj_py.c]
        @test sat.magnetometer.non_ortho_angles ≈ [obj_py.rho, obj_py.lamb, obj_py.phi]
        @test sat.magnetometer.bias             ≈ [obj_py.bias_x, obj_py.bias_y, obj_py.bias_z]

    end
end





