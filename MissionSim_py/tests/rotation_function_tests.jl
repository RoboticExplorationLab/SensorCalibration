include("../rotationFunctions.jl")
using PyCall, Test, Infiltrator, LinearAlgebra
using Plots


__init_rot_functions__()

function get_random_quat()
    r = randn(3); r = r / norm(r) # random axis 
    θ = deg2rad(rand(0:360))

    q = [r * sin(θ/2); cos(θ/2)]
    return q
end

@testset "Hat Function" begin
    # Check Hat
    v = [0 0 0]
    @test hat(v) == [0 0 0; 0 0 0; 0 0 0]
    v = [1 2 3]
    @test hat(v) == [0 -3 2; 3 0 -1; -2 1 0]

    # Compare with Python
    for i = 1:10
        v = 10.0 * randn(3)
        @test hat(v) == py"hat"(v)
    end
end

@testset "q from DCM" begin
    function get_random_dcm()
        Q, R = qr(randn(3,3))
        # to avoid only covering a half sphere...
        Q = Matrix(Q)
        Q[:,1] = Q[:,1] * (2 * (rand() > 0.5) - 1);
        Q[:,2] = det(Q) * Q[:, 2]
        return Q
    end

    # Verify Julia implementation
    dcm = I(3)
    @test q_from_DCM(dcm) == [0.0, 0.0, 0.0, 1.0]  # Identity
    dcm = [0 1 0; -1 0 0; 0 0 1]
    @test q_from_DCM(dcm) ≈ [0.0, 0.0, -sqrt(2)/2, sqrt(2)/2]
    dcm = [0 -1 0; 1 0 0; 0 0 1]
    @test q_from_DCM(dcm) ≈ [0.0, 0.0, sqrt(2)/2, sqrt(2)/2]

    # Compare with Python
    for i = 1:10
        dcm = get_random_dcm()
        @test q_from_DCM(dcm) ≈ py"q_from_DCM"(dcm)
    end
end

@testset "DCM from q" begin
    # Verify Julia implementation
    q = [0.0, 0.0, 0.0, 1.0]
    @test dcm_from_q(q) == I(3)
    q = [0.0, 0.0, -sqrt(2)/2, sqrt(2)/2]
    @test dcm_from_q(q) ≈ [0 1 0; -1 0 0; 0 0 1]
    q = [0.0, 0.0, sqrt(2)/2, sqrt(2)/2]
    @test dcm_from_q(q) ≈ [0 -1 0; 1 0 0; 0 0 1]

    # Compare with Python
    for i = 1:10
        q = get_random_quat()
        @test dcm_from_q(q) ≈ py"dcm_from_q"(q)
    end
end 

@testset "Q Mult" begin
    # Verify Julia implementation
    q = [0, 0, 0, 1]; p = [1, 2, 1, 0.5]; p = p / norm(p)
    @test qmult(q, p) == p

    q = [1 0 0 3]; p = [5 1 -2 0];
    @test qmult(q, p) == [15, 5, -5, -5]
    @test qmult(p, q) == [15, 1, -7, -5]

    # Compare with Python
    for i = 1:10
        q, p = get_random_quat(), get_random_quat()
        @test qmult(q, p) ≈ py"qmult"(q, p)
    end
end

@testset "Q Conj" begin
    q = [1, 2, 3, 4];
    @test qconj(q) == [-1, -2, -3, 4]
    q = [0.0, 1.0, 0.0, 0.0]
    @test qconj(q) == [0.0, -1.0, 0.0, 0.0]

    for i = 1:10
        q = get_random_quat()
        @test qconj(q) ≈ py"qconj"(q)
    end
end

println("Finished testing Rotation Functions")