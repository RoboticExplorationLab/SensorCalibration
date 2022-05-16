# [test/quaternions.jl] 

""" Common functions used when working with quaternions (assumes scalar first!) """

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


# Add test - L(q₂) q₁ == R(q₁) * q₂
function qmult(q₁, q₂)
    return L(q₂) * q₁
end

⊙(q₁::Vector, q₂::Vector) = L(q₂) * q₁

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

function quat2rot(q)
    return H' * L(q) * R(q)' * H
end

