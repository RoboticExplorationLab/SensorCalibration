# [src/quaternions.jl] 

""" Utility functions for working with quaternions. In an effort to be fast, there are a large number of multiple dispatch functions """

""" To Do:
 - 
"""

""" Common functions used when working with quaternions (assumes scalar first!) """


using LinearAlgebra   # For identity matrix 

                                        #################
                                        # `hat` methods #
                                        #################

""" Takes in a three-element vector and returns a skew-symmetric cross-product matrix """
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


                                        #######################
                                        # `L` and `R` Methods #
                                        #######################

""" Takes in a scalar-first unit quaternion and converts it to a left-side matrix for multiplication """
function L(q::SubArray{T})::Matrix{T} where {T}
    qₛ = q[1]
    qᵥ = view(q, 2:4)      

    M = zeros(T, 4, 4)    
    M[1, 1] = q[1]         

    M[1, 2:4] .= -qᵥ       
    M[2:4, 1] .=  qᵥ      
    M[2:4, 2:4] .= lHat(qₛ, qᵥ) 

    return M
end

function L(q::Vector{T})::Matrix{T} where {T} 
    qₛ = q[1]
    qᵥ = view(q, 2:4)      

    M = zeros(T, 4, 4)    
    M[1, 1] = q[1]        

    M[1, 2:4] .= -qᵥ       
    M[2:4, 1] .=  qᵥ       
    M[2:4, 2:4] .= lHat(qₛ, qᵥ) 

    return M
end

function L(q::SVector{4, T})::SMatrix{4, 4, T, 16} where {T} 
    qₛ = q[1]
    qᵥ = view(q, 2:4)      

    M = zeros(T, 4, 4)    
    M[1, 1] = q[1]         

    M[1, 2:4] .= -qᵥ      
    M[2:4, 1] .=  qᵥ       
    M[2:4, 2:4] .= lHat(qₛ, qᵥ) 

    return SMatrix{4, 4, T, 16}(M)
end

""" Takes in a scalar-first unit quaternion and converts it to a right-side matrix for multiplication """
function R(q::Vector{T})::Matrix{T} where {T}
    qₛ = q[1]
    qᵥ = view(q, 2:4)      

    M = zeros(T, 4, 4)    
    M[1, 1] = q[1]         

    M[1, 2:4] .= -qᵥ      
    M[2:4, 1] .=  qᵥ       
    M[2:4, 2:4] .= rHat(qₛ, qᵥ) 

    return M
end

function R(q::SVector{4, T})::Matrix{T} where {T}
    qₛ = q[1]
    qᵥ = view(q, 2:4)     

    M = zeros(T, 4, 4)     
    M[1, 1] = q[1]         

    M[1, 2:4] .= -qᵥ      
    M[2:4, 1] .=  qᵥ       
    M[2:4, 2:4] .= rHat(qₛ, qᵥ) 

    return SMatrix{4, 4, T, 16}(M)
end

                                        #######################
                                        # `H` and `T` Methods #
                                        #######################

const Hs = SMatrix{4, 3, Float64, 12}([zeros(1, 3); I(3)])  # Static equivalent
const H = [zeros(1, 3); I(3)]     # Converts from a 3-element vector to a 4-element vector with 0 scalar part 
const T = [1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 -1]   # Forms the conjugate of q, i.e. q† = Tq   
qconj(q) = [1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 -1] * q

                                        ##############
                                        # Arithmatic #
                                        ##############
function qmult(q₁, q₂)
    return L(q₁) * q₂
end

⊙(q₁::Vector, q₂::Vector) = qmult(q₁, q₂)
⊙(q₁::SVector{4, T}, q₂::SVector{4, T}) where {T} = qmult(q₁, q₂)
⊙(q₁, q₂) = qmult(q₁, q₂)

# Takes q and ω → q̇
qdot(q::Vector{T}, ω::Vector{T})         where {T} = 0.5 * L(q) * H * ω
qdot(q::SVector{4, T}, ω::SVector{3, T}) where {T} = 0.5 * L(q) * Hs * ω
qdot(q::SubArray{T}, ω::SubArray{T})     where {T} = 0.5 * L(q) * H * ω


quat2rot(q::Vector) = H' * L(q) * R(q)' * H
quat2rot(q::SVector{4, T}) where {T} =  SMatrix{3, 3, T, 9}(H' * L(q) * R(q)' * H)



""" Attitude Jacobian to convert between the Jacobian of 3 and 4 element attitude representations """
function G(q::SVector{4, T}) where {T}
    qₛ = q[1] 
    qᵥ = q[2:4]

    G = zeros(T, 4, 3)
    G[1,   :] .= -qᵥ
    G[2:4, :] .= qₛ * I(3) + hat(qᵥ)

    return G
end
att_jac(q) = G(SVector{4, typeof(q)}(q))




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
  
""" 
    Cayley map is a convenient approximation of the error between two quaternions. 

    The magnitude of the Cayley map is about one half of the radians, but for angles under 10-15 degrees
    (more robust would be to use a quaternion logrithm)
"""
function cayley_map(q₁, q₂)
    qₑ = L(q₁)' * q₂
    e  =  (qₑ[2:4] / qₑ[1])
    return e
end

function qErr(q₁, q₂)
    return norm((L(q₁)' * q₂)[2:4]) # Need to subtract [±1, 0, 0, 0], but because its unit we just ignore the scalar
end

""" Hamilton product w/o matrices """
function hamilton(q₁, q₂)
    s₁, v₁ = q₁[1], q₁[2:4]
    s₂, v₂ = q₂[1], q₂[2:4]   

    s = (s₁ * s₂) - dot(v₁, v₂)
    v = (s₁ * v₂ + s₂ * v₁ + cross(v₁, v₂))
    return [s; v]
end;

