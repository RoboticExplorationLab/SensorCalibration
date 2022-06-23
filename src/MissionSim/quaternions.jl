# [test/quaternions.jl] 

""" Common functions used when working with quaternions (assumes scalar first!) """
# This really needs to be tidied up and tested


using LinearAlgebra   # For identity matrix 

# TO DO:
#  - Add in functionality for static arrays ?
#  - Test for type stability and allocations? 
#  - Add in derivative, conj, etc...

function hat(v::Vector{T})::Matrix{T} where {T}
    M = zeros(T, 3, 3)
    M[1, 2] = -v[3]
    M[1, 3] =  v[2]
    M[2, 1] =  v[3] 
    M[2, 3] = -v[1]
    M[3, 1] = -v[2]
    M[3, 2] =  v[1]

    # M =  [0.0   -v[3]  v[2];  
    #       v[3]   0.0  -v[1]; 
    #      -v[2]   v[1]  0.0]

    return M 
end

function hat(v::SVector{3, T})::SMatrix{3, 3, T, 9} where {T}
    M = zeros(T, 3, 3)
    M[1, 2] = -v[3]
    M[1, 3] =  v[2]
    M[2, 1] =  v[3] 
    M[2, 3] = -v[1]
    M[3, 1] = -v[2]
    M[3, 2] =  v[1]

    return SMatrix{3, 3, T, 9}(M)
end

function hat(v::SubArray{T})::Matrix{T} where {T}
    M = zeros(T, 3, 3)
    M[1, 2] = -v[3]
    M[1, 3] =  v[2]
    M[2, 1] =  v[3] 
    M[2, 3] = -v[1]
    M[3, 1] = -v[2]
    M[3, 2] =  v[1]

    return M 
end

function hat(v)
    return V = [0    -v[3]  v[2];
                v[3]  0.0  -v[1];
               -v[2]  v[1]  0.0]
end

function lHat(s::T, v::SubArray{T})::Matrix{T} where {T}
    """ qₛ * I + hat(qᵥ), but  faster and fewer allocs"""
    M = zeros(T, 3, 3)
    M[1, 2] = -v[3]
    M[1, 3] =  v[2]
    M[2, 1] =  v[3] 
    M[2, 3] = -v[1]
    M[3, 1] = -v[2]
    M[3, 2] =  v[1]

    M[1, 1] = s 
    M[2, 2] = s 
    M[3, 3] = s

    return M 
end

function rHat(s::T, v::SubArray{T})::Matrix{T} where {T}
    """ qₛ * I - hat(qᵥ), but  faster and fewer allocs"""
    M = zeros(T, 3, 3)
    M[1, 2] =  v[3]
    M[1, 3] = -v[2]
    M[2, 1] = -v[3] 
    M[2, 3] =  v[1]
    M[3, 1] =  v[2]
    M[3, 2] = -v[1]

    M[1, 1] = s 
    M[2, 2] = s 
    M[3, 3] = s

    return M 
end

function L(q::SubArray{T})::Matrix{T} where {T}
    qₛ = q[1]
    qᵥ = view(q, 2:4)      # ≈ 2ns for the pair

    M = zeros(T, 4, 4)     # 25 ns, 1 alloc 200 byte
    M[1, 1] = q[1]         # 1 ns, 0 alloc

    M[1, 2:4] .= -qᵥ       # 30 ns, 1 alloc, 112 bytes  (b/c negative)
    M[2:4, 1] .=  qᵥ       # 20 ns, 0 alloc
    M[2:4, 2:4] .= lHat(qₛ, qᵥ) #  About 150ns faster, 3 fewer alloc, and 400 fewer bytes

    return M
end

function L(q::Vector{T})::Matrix{T} where {T} 
    qₛ = q[1]
    qᵥ = view(q, 2:4)      # ≈ 2ns for the pair

    M = zeros(T, 4, 4)     # 25 ns, 1 alloc 200 byte
    M[1, 1] = q[1]         # 1 ns, 0 alloc

    M[1, 2:4] .= -qᵥ       # 30 ns, 1 alloc, 112 bytes  (b/c negative)
    M[2:4, 1] .=  qᵥ       # 20 ns, 0 alloc
    M[2:4, 2:4] .= lHat(qₛ, qᵥ) #  About 150ns faster, 3 fewer alloc, and 400 fewer bytes
    # M[2:4, 2:4] .= qₛ * I(3) + hat(qᵥ)  # 225, 5 alloc, 750 bytes  (50 ns, 1 alloc, 160 from hat)

    return M
end

function L(q::SVector{4, T})::SMatrix{4, 4, T, 16} where {T} 
    qₛ = q[1]
    qᵥ = view(q, 2:4)      # ≈ 2ns for the pair

    M = zeros(T, 4, 4)     # 25 ns, 1 alloc 200 byte
    M[1, 1] = q[1]         # 1 ns, 0 alloc

    M[1, 2:4] .= -qᵥ       # 30 ns, 1 alloc, 112 bytes  (b/c negative)
    M[2:4, 1] .=  qᵥ       # 20 ns, 0 alloc
    M[2:4, 2:4] .= lHat(qₛ, qᵥ) #  About 150ns faster, 3 fewer alloc, and 400 fewer bytes
    # M[2:4, 2:4] .= qₛ * I(3) + hat(qᵥ)  # 225, 5 alloc, 750 bytes  (50 ns, 1 alloc, 160 from hat)

    return SMatrix{4, 4, T, 16}(M)
end

function R(q::Vector{T})::Matrix{T} where {T}
    qₛ = q[1]
    qᵥ = view(q, 2:4)      # ≈ 2ns for the pair

    M = zeros(T, 4, 4)     # 25 ns, 1 alloc 200 byte
    M[1, 1] = q[1]         # 1 ns, 0 alloc

    M[1, 2:4] .= -qᵥ       # 30 ns, 1 alloc, 112 bytes  (b/c negative)
    M[2:4, 1] .=  qᵥ       # 20 ns, 0 alloc
    M[2:4, 2:4] .= rHat(qₛ, qᵥ) #  About 150ns slower, 3 fewer alloc, and 400 fewer bytes
    # M[2:4, 2:4] .= qₛ * I(3) - hat(qᵥ)

    # # NOTE the below method is the same, but about 10x slower and 10x #allocations
    # qₛ, qᵥ = q[1], q[2:4]
    # M = [qₛ     -qᵥ';                
    #      qᵥ  (qₛ * I - hat(qᵥ)) ]

    return M
end

function R(q::SVector{4, T})::Matrix{T} where {T}
    qₛ = q[1]
    qᵥ = view(q, 2:4)      # ≈ 2ns for the pair

    M = zeros(T, 4, 4)     # 25 ns, 1 alloc 200 byte
    M[1, 1] = q[1]         # 1 ns, 0 alloc

    M[1, 2:4] .= -qᵥ       # 30 ns, 1 alloc, 112 bytes  (b/c negative)
    M[2:4, 1] .=  qᵥ       # 20 ns, 0 alloc
    M[2:4, 2:4] .= rHat(qₛ, qᵥ) #  About 150ns slower, 3 fewer alloc, and 400 fewer bytes

    # # NOTE the below method is the same, but about 10x slower and 10x #allocations
    # qₛ, qᵥ = q[1], q[2:4]
    # M = [qₛ     -qᵥ';                
    #      qᵥ  (qₛ * I - hat(qᵥ)) ]

    return SMatrix{4, 4, T, 16}(M)
end


const Hs = SMatrix{4, 3, Float64, 12}([zeros(1, 3); I(3)])
const H = [zeros(1, 3); I(3)]     # Converts from a 3-element vector to a 4-element vector with 0 scalar part 

const T = [1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 -1]   # Forms the conjugate of q, i.e. q† = Tq   

qconj(q) = [1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 -1] * q

# Add test - L(q₂) q₁ == R(q₁) * q₂
function qmult(q₁, q₂)
    return L(q₁) * q₂
end

⊙(q₁::Vector, q₂::Vector) = qmult(q₁, q₂)
⊙(q₁::SVector{4, T}, q₂::SVector{4, T}) where {T} = qmult(q₁, q₂)
⊙(q₁, q₂) = qmult(q₁, q₂)

# Takes q and ω → q̇
function qdot(q::Vector{T}, ω::Vector{T}) where {T}
    q̇ = 0.5 * L(q) * H * ω

    return q̇
end

function qdot(q::SVector{4, T}, ω::SVector{3, T}) where {T}
    q̇ = 0.5 * L(q) * Hs * ω

    return q̇
end

function qdot(q::SubArray{T}, ω::SubArray{T}) where {T}
    q̇ = 0.5 * L(q) * H * ω

    return q̇
end

function quat2rot(q::Vector)
    return H' * L(q) * R(q)' * H
end

function quat2rot(q::SVector{4, T}) where {T} 
    return SMatrix{3, 3, T, 9}(H' * L(q) * R(q)' * H)
end

function att_jac(q)
    return G(SVector{4, typeof(q)}(q))
end

function G(q::SVector{4, T}) where {T}
    qₛ = q[1] 
    qᵥ = q[2:4]

    G = zeros(T, 4, 3)
    G[1,   :] .= -qᵥ
    G[2:4, :] .= qₛ * I(3) + hat(qᵥ)

    return G
end



function rot2quat(R)
    # Converts a Rotation matrix into a scalar-first quaternion (From Kevin)

    T = R[1,1] + R[2,2] + R[3,3];
    if (T > R[1,1]) && (T > R[2,2]) && (T>R[3,3])
        q4 = .5*sqrt(1+T);
        r  = .25/q4;
        q1 = (R[3,2] - R[2,3])*r;
        q2 = (R[1,3] - R[3,1])*r;
        q3 = (R[2,1] - R[1,2])*r;
    elseif R[1,1]>R[2,2] && R[1,1]>R[3,3]
        q1 = .5*sqrt(1-T + 2*R[1,1]);
        r  = .25/q1;
        q4 = (R[3,2] - R[2,3])*r;
        q2 = (R[1,2] + R[2,1])*r;
        q3 = (R[1,3] + R[3,1])*r;
    elseif R[2,2]>R[3,3]
        q2 = .5*sqrt(1-T + 2*R[2,2]);
        r  = .25/q2;
        q4 = (R[1,3] - R[3,1])*r;
        q1 = (R[1,2] + R[2,1])*r;
        q3 = (R[2,3] + R[3,2])*r;
    else
        q3 = .5*sqrt(1-T + 2*R[3,3]);
        r  = .25/q3;
        q4 = (R[2,1] - R[1,2])*r;
        q1 = (R[1,3] + R[3,1])*r;
        q2 = (R[2,3] + R[3,2])*r;
    end

    q = [q4; q1; q2; q3] 

    if q4 < 0
       q = -q;
    end

    return q
end

function cayley_map(q₁, q₂)
    qₑ = L(q₁)' * q₂
    e  =  (qₑ[2:4] / qₑ[1])
    return e
end

function qErr(q₁, q₂)
    return norm((L(q₁)' * q₂) - [1, 0, 0, 0])
end

""" Hamilton product w/o matrices """
function hamilton(q₁, q₂)
    s₁, v₁ = q₁[1], q₁[2:4]
    s₂, v₂ = q₂[1], q₂[2:4]   

    s = (s₁ * s₂) - dot(v₁, v₂)
    v = (s₁ * v₂ + s₂ * v₁ + cross(v₁, v₂))
    return [s; v]
end;

# @testset "Rot 2 Quat" begin 
#     N = 1000
#     ts = zeros(N)
#     for i = 1:N

#         q = randn(4); q /= norm(q)
#         q̂ = rot2quat(quat2rot(q))
#         ts[i] = ((q ≈ q̂) || (q ≈ -q̂))
#     end
#     @test all(ts .== 1)
# end

