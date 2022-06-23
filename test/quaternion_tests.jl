# [tests/quaternion_tests.jl]

""" To Do:
 - Finish the tests, only a few are included
 - make Quat into a module?
"""

using Test

@testset "Quaternions" begin 
    using LinearAlgebra, BenchmarkTools, StaticArrays
    (@isdefined qdot) ? nothing : include("../src/MissionSim/quaternions.jl");

    function rand_quat() 
        q = randn(4)
        return q / norm(q)
    end


    q = [1, 0, 0, 0]
    @test (H' * L(q) * R(q)' * H) == I(3)
    q = [0, 1, 0, 0]
    @test (H' * L(q) * R(q)' * H) == [1.0 0 0; 0 -1 0; 0 0 -1]
    q = [sqrt(2)/2, 0, 0, sqrt(2)/2.0]
    @test (H' * L(q) * R(q)' * H) ≈ [0.0 -1 0; 1 0 0; 0 0 1]
    
    # Hᵀ L(q) R(q)' H = [3 x 3], SO(3)
    for i = 1:10
        q = rand_quat()
        Q = (H' * L(q) * R(q)' * H) 
        @test (Q' * Q) ≈ I(3)
        @test det(Q) ≈ 1.0
        @test Q == quat2rot(q)
    end

    # hat(v) * v = 0⃗, trace(hat(v)) = 0
    for i = 1:10 
        v = randn(3)
        v̂ = hat(v)
        @test norm(v̂ * v) == 0
        @test tr(v̂) == 0
    end

    # L(Tq) = L(q)' = L(q)⁻¹, same with R
    for i = 1:10
        q = rand_quat()
        @test (L(T * q) == L(q)') && (L(q)' ≈ inv(L(q)))
        @test (R(T * q) == R(q)') && (R(q)' ≈ inv(R(q)))
    end

    # L(q2)q1 == R(q1)q2
    for i = 1:10 
        q₁ = rand_quat()
        q₂ = rand_quat()
        @test L(q₂) * q₁ ≈ R(q₁) * q₂
    end

    # Rotate a vector (using old school method, too)
    for i = 1:10 
        v = randn(3) 
        q = rand_quat()
        v₁ = H' * L(q) * R(q)' * H * v 
        v₂ = ((q ⊙ [0; v]) ⊙ (T * q))[2:4]
        v₃ = hamilton(hamilton(q, [0; v]), T * q)[2:4]
        @test (v₁ ≈ v₂) && (v₂ ≈ v₃) && (v₃ ≈ v₁)
    end
    
    # Type stability and all that 

end;