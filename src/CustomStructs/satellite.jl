# [src/CustomStructs/satellite.jl]

"""
    SATELLITE{S, T} -> J, magnetometer, diodes, state, covariance

      Satellite struct that contains the inertia matrix for the CubeSat, as well
    as magnetometers and diodes. Also tracks the SATELLITE state and covariance 
    (not the environment state state).
"""
struct SATELLITE{S, T}
    J::SMatrix{3, 3, T, 9}               # Inertia Matrix of satellite  
    magnetometer::MAGNETOMETER{T}        # Mag calibration values      
    diodes::DIODES{S, T}                 # Diode calibration values 
    state::SAT_STATE{T}                  # Satellite state 
    covariance::Matrix{T}                # Satellite state covariance 

    function SATELLITE(J::SMatrix{3, 3, T, 9}, mag::MAGNETOMETER{T}, dio::DIODES{S, T}, 
                        sta::SAT_STATE{T}, cov::Matrix{T})  where {S, T}
        """ Primary Constructor  """

        new{S, T}(SMatrix{3, 3, T, 9}(J), mag, dio, sta, cov)
    end

    function SATELLITE(; J = nothing, mag::MAGNETOMETER{T} = MAGNETOMETER(), dio::DIODES = DIODES(), 
                            sta::SAT_STATE = SAT_STATE(), cov::Matrix{T} = SAT_COVARIANCE().Σ, ideal::Bool = false) where {T}
        """ Generate a SATELLITE with random parameters """

        if ideal 
            _J = isnothing(J) ? SMatrix{3, 3, Float64, 9}(0.2 * I(3)) : J
            _mag = MAGNETOMETER(; ideal = true)
            _dio = DIODES(; ideal = true)
            _sta = SAT_STATE(; ideal = true)
            _cov = cholesky(Hermitian(Matrix(cov))).U

            return SATELLITE(_J, _mag, _dio, _sta, Matrix(_cov))
        end


        if !isnothing(J)
            return SATELLITE(J, mag, dio, sta, cov)
        else
            m = 1.5 + 0.5 * rand()
            l, w, h = 1 .+ 0.5 * rand(3)
    
            J = zeros(T, 3, 3)
            J[1, 1] = (m / 12.0) * (l^2 + w^2) 
            J[2, 2] = (m / 12.0) * (l^2 + h^2)
            J[3, 3] = (m / 12.0) * (h^2 + w^2)
            
            return SATELLITE(SMatrix{3, 3, T, 9}(J), mag, dio, sta, cov)
        end
    end

end
