# [src/MissionSim/CustomStructs/sat_state.jl]

"""
    SAT_STATE{N, T} -> q, β, C, α, ϵ

      State of the SATELLITE (not to be confused with the environment state 'state'). Comprised of 
    the attitude 'q' (as a scalar-first unit quaternion) and the gyroscope bias 'β', as well as the 
    diode calibration parameters: Calibration value 'C', and the azimuth 'α' and elevation 'ϵ' angles 
    for the surface normals of each diode.
"""
struct SAT_STATE{T}
    q::SVector{4, T}     # Scalar-first unit quaternion representing orientation of the satellite 
    β::SVector{3, T}     # Gyroscope bias 
    # C::SVector{N, T}     # Current scale factor for each diode 
    # α::SVector{N, T}     # Azimuth angle corresponding to each diode 
    # ϵ::SVector{N, T}     # Elevation angle corresponding to each diode 

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

        if ideal  # No rotation, no bias, unit scale factors, and perfect 

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
# Base.deepcopy(ss::SAT_STATE) = SAT_STATE(; q = deepcopy(ss.q), β = deepcopy(ss.β))

"""
    update_state(x; q, β, C, α, ϵ)

      'Updates' portions of a SAT_STATE struct by creating a new struct 
    and copying relevant stuff over. 
"""
function update_state(x::SAT_STATE{T}; q = nothing, β = nothing) where {N, T}

    _q = isnothing(q) ? Array(x.q)  : Array(q) 
    _β = isnothing(β) ? Array(x.β)  : Array(β)
    # _q = isnothing(q) ? deepcopy( SVector{4, T}(vcat([x.q;]...))) : q 
    # _β = isnothing(β) ? SVector{3, T}(deepcopy(vcat([x.β;]...))) : SVector{3, T}(deepcopy(vcat([β;]...)))

    return SAT_STATE(SVector{4, T}(_q), SVector{3, T}(_β))
end
