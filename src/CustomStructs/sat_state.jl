# [src/CustomStructs/sat_state.jl]

"""
    SAT_STATE{N, T} -> q, β

      State of the SATELLITE (not to be confused with the environment state 'state'). Comprised of 
    the attitude 'q' (as a scalar-first unit quaternion) and the gyroscope bias 'β'. 
"""
struct SAT_STATE{T}
    q::SVector{4, T}     # Scalar-first unit quaternion representing orientation of the satellite 
    β::SVector{3, T}     # Gyroscope bias  

    function SAT_STATE(q::SVector{4, T}, β::SVector{3, T}) where {T}
        """ Primary Constructor """

        q = (norm(q) == 1) ? q  : (q / norm(q))
        q = (q[1] < 0)  ?  -q : q                   # Because quaternions double cover, we force the scalar to be positive
        new{T}(q, β)
    end

    function SAT_STATE(; q = nothing, β = nothing, ideal = false)
        """ Randomly generate either a noisy or ideal satellite state """

        function rand_quat()
            t = randn(4)
            return t / norm(t)
        end

        if ideal  # No rotation, no bias
            _q = SVector{4, Float64}(1.0, 0, 0, 0)
            _β = SVector{3, Float64}(0.0, 0.0, 0.0)

            return SAT_STATE(_q, _β)
        else

            _q = isnothing(q)  ?  rand_quat()  :  q 
            _β = isnothing(β)  ?  randn(3)     :  β 
    
            return SAT_STATE(SVector{4, Float64}(_q),  SVector{3, Float64}(_β))
        end
    end
end

"""
    update_state(x; q, β, C, α, ϵ)

      'Updates' portions of a SAT_STATE struct by creating a new struct 
    and copying relevant stuff over. 
"""
function update_state(x::SAT_STATE{T}; q = nothing, β = nothing) where {N, T}

    _q = isnothing(q) ? Array(x.q)  : Array(q) 
    _β = isnothing(β) ? Array(x.β)  : Array(β)
    
    return SAT_STATE(SVector{4, T}(_q), SVector{3, T}(_β))
end
