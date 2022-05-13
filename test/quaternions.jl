# [test/quaternions.jl] 

""" Common functions used when working with quaternions (assumes scalar first!) """

using LinearAlgebra   # For identity matrix 

# TO DO:
#  - Add in functionality for static arrays ?
#  - Test for type stability and allocations? 
#  - Add in derivative, conj, etc...


function L(q::Vector{T})::Matrix{T} where {T} 
    qₛ, qᵥ = q[1], q[2:4]

    M = [qₛ     -qᵥ'; 
         qᵥ  (qₛ * I + hat(qᵥ)) ]

    return M
end

function R(q::Vector{T})::Matrix{T} where {T}
    qₛ, qᵥ = q[1], q[2:4]

    M = [qₛ     -qᵥ'; 
         qᵥ  (qₛ * I - hat(qᵥ)) ]

    return M
end

function hat(v::Vector{T})::Matrix{T} where {T}
    M =  [0.0   -v[3]  v[2];  
          v[3]   0.0  -v[1]; 
         -v[2]   v[1]  0.0]

    return M 
end

const H = [zeros(1, 3); I(3)]     # Converts from a 3-element vector to a 4-element vector with 0 scalar part 

const T = [1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 -1]   # Forms the conjugate of q, i.e. q† = Tq   


# Add test - L(q₂) q₁ == R(q₁) * q₂
function qmult(q₁, q₂)
    return L(q₂) * q₁
end

⊙(q₁::Vector, q₂::Vector) = L(q₂) * q₁

# Takes q and ω → q̇
function diff(q::Vector, ω::Vector) 
    return 0.5 * L(q) * H * ω
end