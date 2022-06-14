# [src/MissionSim/Estimator/mekf.jl]

""" To Do:

 - A (in prediction) does not appear to be correct..., and maybe Hb and Hi, too
 - current_measurement is slow (~ ms)

 - Read up on IMU stuff to get a fill for the noise W and V and general covariance stuff
 - Read up on setting initial covariances too

 - unit for Bᴵ, Bᴮ, Ĩ, etc...?

 - comment

 - current + diode cal does *worse* with negative middle term???
"""

# Note that I am not doing σ_gyro == σ_orient b/c units; σ_sun, mag are waaaay bigger than the others (0.05 vs 4e-6)
"""
    MEKF_DATA(Wchol, Vchol)

      Struct containing the process and measurement noise matrices for the
    MEKF. Because this is a square-root MEKF, we are using Cholesky-factorized 
    versions of both W and V. 
"""
struct MEKF_DATA{T}

    Wchol::UpperTriangular{T}      # Process noise for state [ϕ, β, C, α, ϵ]
    Vchol::UpperTriangular{T}      # Measurement noise for measurements [B̃ᴮ, Ĩ]

    function MEKF_DATA(Wchol::UpperTriangular{T}, Vchol::UpperTriangular{T}) where {T}
        """ Primary Constructor """
        new{T}(Wchol, Vchol)
    end

    function MEKF_DATA(; N = 6, σ_sun = deg2rad(3.0), σ_mag = deg2rad(3.0),
                            σC = 0.1, σα = deg2rad(3.), σϵ = deg2rad(3.),  # Was 0.1, 3, 3
                            gyro_bias_instability = 0.8,   # in deg/hour  
                            angle_random_walk     = 0.06)  # in deg/sqrt(hour)

        """ Generates a new MEKF_DATA by generating the W and V matrices with the provided values """

        ## Process Noise:
        σ_gyro = deg2rad(gyro_bias_instability) / 3600.0  # Convert (deg/hour) to (rad/sec)
        σ_bias = deg2rad(angle_random_walk) / 60.0        # Convert (deg/sqrt(hour)) to ( rad/sqrt(s) )
        
        W = Diagonal( [σ_gyro * ones(3); σ_bias * ones(3); σC * ones(N); σα * ones(N); σϵ * ones(N)].^2 )
        Wchol = chol(Matrix(W)) 
        
        ## Measurement Noise:
        V = Diagonal( [σ_mag * ones(3); σ_sun * ones(N)].^2 )
        Vchol = chol(Matrix(V))
    
        MEKF_DATA(Wchol, Vchol)
    end

end


# How do i go from noise w, v to covariance P? Right now I just set it to be either 1 or 3σ...
# Is it bad to be re-diagonalizing? Would need to make new constructor for sat_cov
"""
    reset_cov!(sat; reset_calibration = true, σϕ = 10*, σβ = 10*, σC = 0.1, σα = 3*, σϵ = 3*)

      Resets the covariances for a satellite (most likely when coming out of eclipse or 
    transitioning between different states). Completely resets the attitude portion, but for the bias 
    the covariance is just increased byt 50%. Can also reset the covariance corresponding to the 
    calibration values if the flag is set to true.

    Note that this covariance uses 3 parameters for the attitude to prevent dropping rank.
"""
function reset_cov!(sat::SATELLITE{N, T}; reset_calibration = false, σϕ = deg2rad(10), σβ = deg2rad(10), 
                        σC = 0.1, σα = deg2rad(3), σϵ = deg2rad(3)) where {N, T}
    
    # # Update Covariances 
    # Reset q, increase β, pass C, α, and ϵ
    Σϕ = diagm( σϕ * ones(3) )  # Reset (remember this is Cholesky so no σ²)
    Σβ = 1.5 * sat.covariance[4:6, 4:6]    # Increase a bit due to unknown propagation in eclipse
    
    sat.covariance[1:3, 1:3] .= Σϕ
    sat.covariance[4:6, 4:6] .= Σβ

    if reset_calibration 
        i₀ = 7; i₁ = i₀ + N - 1
        sat.covariance[i₀:i₁, i₀:i₁] .= diagm(σC * ones(N))

        i₀ = i₁ + 1; i₁ = i₀ + N - 1
        sat.covariance[i₀:i₁, i₀:i₁] .= diagm(σα * ones(N))

        i₀ = i₁ + 1; i₁ = i₀ + N - 1
        sat.covariance[i₀:i₁, i₀:i₁] .= diagm(σϵ * ones(N))
    end

    return Nothing
end

"""
    estimate(sat, sens, noise, alb, t, dt; use_albedo = true, calibrate_diodes = true)

      Calls a single iteration of the multiplicative extended Kalman filter, updating the 
    satellite attitude and bias estimates. Additionally, if the 'calibrate_diodes' flag is 
    set to 'true', this updates estimates of the diode calibration factors (C, α, and ϵ).
    Note that this function should NOT be called during eclipses.

    Arguments:
      - `sat`:    Struct containing current estimate of satellite values      |  SATELLITE
      - `sens`:   Struct containing the most recent sensor measurements       |  SENSORS 
      - `noise`:  Struct containing the noise matrices V, W (as Cholesky)     |  MEKF_DATA
      - `alb`:    Struct containing the reflectivity data for Earth's Albedo  |  ALBEDO 
      - `t`:      Current time, as an Epoch                                   |  Epoch 
      - `dt`:     Size of each time step                                      |  Scalar
      - `use_albedo`:  (Opt) Whether or not to factor in Earth's albedo  (defaults to 'true')     |  Bool 
      - `calibrate_diodes`:  (Opt) Whether or not to update diode estimates (defaults to 'true')  |  Bool

    Returns:
        - sat:  Updated satellite struct containing new estimates           |  SATELLITE
"""  
function estimate(sat::SATELLITE{N, T}, sens::SENSORS{N, T}, noise::MEKF_DATA{T}, alb::ALBEDO, t::Epoch, dt::Real; use_albedo = true, calibrate_diodes = true) where {N, T}

    # Prepare data 
    sᴵₑ_est = sun_position(t)               # Estimated Sun-Earth vector in inertial frame
    if eclipse_conical(sᴵₑ_est, vcat([sens.pos;]...)) < 0.001 
        @warn "Shouldnt have an eclipse in MEKF!"
    end

    sᴵ_est = sᴵₑ_est - sens.pos             # Estimated sun vector
    Bᴵ_est = SVector{3, T}(IGRF13(sens.pos, t))           # Estimated magnetic field vector in inertial

    # Pass into square-root MEKF 
    Pchol = (calibrate_diodes) ? sat.covariance : sat.covariance[1:6, 1:6]
    Wchol = (calibrate_diodes) ? noise.Wchol    : UpperTriangular(noise.Wchol[1:6, 1:6])
    x⁺, diodes⁺, Pchol⁺ = sqrt_mekf(sat.state, sat.diodes, Pchol, alb, sens.gyro, Bᴵ_est, sᴵ_est,
                             sens.magnetometer, sens.diodes, sens.pos, dt, Wchol, noise.Vchol; 
                             use_albedo = use_albedo, calibrate_diodes = calibrate_diodes)


    # Process result and update sat (by making a new copy)

    # 
    sat⁺ = SATELLITE(deepcopy(sat.J), MAGNETOMETER(sat.magnetometer.scale_factors, sat.magnetometer.non_ortho_angles, sat.magnetometer.bias),
                            diodes⁺, x⁺, copy(sat.covariance) )
    # sat⁺ = SATELLITE(; J = deepcopy(sat.J), dio = deepcopy(sat.diodes), mag = deepcopy(sat.magnetometer), cov = copy(sat.covariance), 
    #                     sta = SAT_STATE(q = x⁺.q, β = x⁺.β, C = x⁺.C, α = x⁺.α, ϵ = x⁺.ϵ) )   
    
    if calibrate_diodes 
        sat⁺.covariance .= Pchol⁺
    else
        sat⁺.covariance[1:6, 1:6] .= Pchol⁺
    end

    return sat⁺
end


# TESTS: Covariance should decrease if meas is good? 
# Does covariance really need to be a struct?
# Test on something real simple to see if it converges?
# W and V are alreaady chol'd
# ASSUMES SMALL ANGLE APPROXIMATION HOLDS, uses 3-param

# Does having Bᴮ not be unit give it heavy weighting bias or is that taken care of with L?
"""
    sqrt_mekf(x, Pchol, alb, ω, Bᴵ, sᴵ, B̃ᴮ, Ĩ, pos, dt, W, V; E_am₀ = 1366.9, use_albedo = true, calibrate_diodes = true)

    
"""
function sqrt_mekf(x::SAT_STATE{T}, diodes::DIODES{N, T}, Pchol::Matrix{T}, alb::ALBEDO, 
                   ω::SVector{3, T}, Bᴵ::SVector{3, T}, sᴵ::SVector{3, T}, B̃ᴮ::SVector{3, T}, 
                   Ĩ::SVector{N, T}, pos::SVector{3, T}, dt::T, W::UpperTriangular, V::UpperTriangular; 
                   E_am₀ = 1366.9, use_albedo = true, calibrate_diodes = true) where {N, T}
    

    # Prediction 
    # x_pred, A = prediction(x, ᴵQᴮ * (ω - x.β) + x.β, dt, N; calibrate_diodes = false)
    x_pred, A = prediction(x, ω, dt, N; calibrate_diodes = calibrate_diodes)
    Pchol_pred = qrᵣ([Pchol * A'; W])  # Update covariance (as Cholesky)

    # Measurement 
    Bᴮ_exp, H_mag = mag_measurement(x_pred, Bᴵ / norm(Bᴵ), N; calibrate_diodes = calibrate_diodes)
    I_exp,  H_cur = current_measurement(x_pred, diodes, sᴵ, pos, alb; E_am₀ = E_am₀,
                                             use_albedo = use_albedo, calibrate_diodes = calibrate_diodes)

    # Innovation 
    z = [(B̃ᴮ / norm(B̃ᴮ)) - Bᴮ_exp; (Ĩ / norm(Ĩ)) - I_exp]

    H_mag *= 2 
    H_cur[:, 1:3] *= 2
    H = [H_mag; H_cur]

    # Kalman Gain (how much we trust the measurement over the dynamics)
    Pchol_yy = qrᵣ([Pchol_pred * H'; V])
    L = (((Pchol_pred' * Pchol_pred * H') / Pchol_yy) / Pchol_yy')

    # Update 
    x⁺, diodes⁺ = update(x_pred, diodes, L, z; calibrate_diodes = calibrate_diodes) 
    Pchol⁺ = qrᵣ([ Pchol_pred * (I - L * H)'; V * L']); 

    @debug (norm(Pchol) < norm(Pchol_pred)) && (norm(Pchol) > norm(Pchol⁺))  # Goes up with prediction, down with measurement

    # #### REMOVE ! 
    # θ_err = acos( (Bᴮ_exp / norm(Bᴮ_exp))' * (B̃ᴮ / norm(B̃ᴮ)) )
    # nP   = norm(Pchol)
    # nPp  = norm(Pchol_pred)
    # nP⁺  = norm(Pchol⁺)
    # nPd  = norm(Pchol[7:end, 7:end])
    # nP⁺d = norm(Pchol⁺[7:end, 7:end])
    # @debug (x_pred !== x) && (x⁺ !== x_pred)
    # @debug (Pchol !== Pchol_pred) && (Pchol_pred !== Pchol⁺)

    # @infiltrate
    # #########

    return x⁺, diodes⁺, Pchol⁺
end

# Add in new tests, comment and clean it up 
# ForwardDIff test does NOT pass (especially the bias term. May be 'good enough'?)
# PANIC - I dont know if we should be using ±hat(γ) b/c Zacs notes say minus. TEST!

# WHAT IF I do the standard q̇ and then use the attitude Jacobian (somehow...?)
"""
function prediction(x::SAT_STATE{T}, ω::SVector{3, T}, dt::T, N::Int; calibrate_diodes::Bool = true) where {T}

    # Generate axis-angle representation
    γ  = ω - x.β     # Corrected angular velocity 
    nγ = norm(γ)     # Magnitude of corrected angular velocity 

    # Predict next orientation (as a quaternion)
    r  = γ / nγ       # Unit axis of rotation
    θ  = nγ * dt      # Angular rotaation about unit axis
    q⁺ = (nγ == 0.0) ? x.q :
                       qmult(x.q, [cos(θ / 2); r*sin(θ / 2)])  # <- right side is axis-angle to unit quaternion

    # q⁺ = (nγ == 0.0) ? x.q :
    #                    qmult([cos(θ / 2); r*sin(θ / 2)], x.q, )  # <- left side is axis-angle to unit quaternion


    @assert (norm(q⁺) ≈ 1.0) "Unit quaternion isn't unit!"    
    x⁺ = update_state(x; q = q⁺)  # Doesn't matter if we are calibrating diodes, only q changes in prediction

    # Calculate Jacobian ∂f/∂x

    # @debug γ̂ is currently positive (Zac's is negative)
    γ̂ = hat(γ / nγ)  # Skew-symmetric matrix 
    R = (nγ == 0.0) ? I(3)  :  # Avoid the divide-by-zero error
                      I(3) + (γ̂ ) * sin(nγ) + ((γ̂ )^2) * (1 - cos(nγ)); # Rodrigues formula
                      
    A = zeros(T, 6, 6)     # Jacobian of f(x) wrt x,   A  =   [R   -dt I₃;  0₃  I₃]
    A[1:3, 1:3] .= R  # or -hat(ω?)
    A[1:3, 4:6] .= -dt * I(3) 
    A[4:6, 4:6] .= I(3)
    # A[4:6, 1:3] .= zeros(3, 3)

    if calibrate_diodes
        # Add in additional states if calibrating (i.e., x = [q, β, C, α, ϵ])
        Ac = zeros(T, 6 + 3 * N, 6 + 3 * N)  # Ac = [A  0; 0 I]
        Ac[1:6, 1:6] .= A 
        Ac[7:end, 7:end] .= I(3 * N)
        
        return x⁺, Ac
    else
        return x⁺, A
    end
end;
"""

# """
    function attitude_jacobian(q)
        G = [-q[2:4]' ; 
            q[1] * I + hat(q[2:4])]
        return SMatrix{4, 3, Float64, 12}(G)
    end
    function prediction(x::SAT_STATE{T}, ω::SVector{3, T}, dt::T, N::Int; calibrate_diodes::Bool = true) where {T}

        function nextq(q, ω, β)
            # Generate axis-angle representation
            γ  = ω - β     # Corrected angular velocity 
            nγ = norm(γ)     # Magnitude of corrected angular velocity 
        
            # Predict next orientation (as a quaternion)
            r  = γ / nγ       # Unit axis of rotation
            θ  = nγ * dt      # Angular rotaation about unit axis
            q⁺ = (nγ == 0.0) ?  q :
                                qmult(q, [cos(θ / 2); r*sin(θ / 2)])

            return q⁺, γ, nγ
        end

        q⁺, γ, nγ = nextq(x.q, ω, x.β)

        x⁺ = update_state(x; q = q⁺)  # Doesn't matter if we are calibrating diodes, only q changes in prediction

        # Calculate Jacobian ∂f/∂x

        γ̂ = hat(γ / nγ)  # Skew-symmetric matrix 
        R = (nγ == 0.0) ? I(3)  :  # Avoid the divide-by-zero error
                        I(3) + (γ̂ ) * sin(nγ) + ((γ̂ )^2) * (1 - cos(nγ)); # Rodrigues formula

        @debug R == exp(nγ * hat(γ / nγ) )
                        
        A = zeros(T, 6, 6)     # Jacobian of f(x) wrt x,   A  =   [R   -dt I₃;  0₃  I₃]


        _nextq_q(_q) = nextq(_q, ω, x.β)[1]
        _nextq_β(_β) = nextq(x.q, ω, _β)[1]

        q⁺ₐ = Array(q⁺) # To prevent FD from somehow overwriting q⁺...
        A[1:3, 1:3] .= attitude_jacobian(q⁺ₐ)' * ForwardDiff.jacobian(_nextq_q, x.q) * attitude_jacobian(x.q)   # df/dq with Attitude Jacobian 
        A[1:3, 4:6] .= attitude_jacobian(q⁺ₐ)' * ForwardDiff.jacobian(_nextq_β, x.β) 
        # A[1:3, 1:3] .= R 
        # A[1:3, 4:6] .= -dt * I(3) 
        A[4:6, 4:6] .= I(3)
        # A[4:6, 1:3] .= zeros(3, 3)

        
        if calibrate_diodes
            # Add in additional states if calibrating (i.e., x = [q, β, C, α, ϵ])
            Ac = zeros(T, 6 + 3 * N, 6 + 3 * N)  # Ac = [A  0; 0 I]
            Ac[1:6, 1:6] .= A 
            Ac[7:end, 7:end] .= I(3 * N)
            
            return x⁺, Ac
        else
            return x⁺, A
        end
    end;
# """


# Clean up and add new tests
# Test - Our H seems to be half of forward diff?
# Good comments in the OG
# Should the Bs be unit?
function mag_measurement(x::SAT_STATE{T}, Bᴵ::SVector{3, T}, N::Int; calibrate_diodes::Bool = true) where {T}

    ## Generate expected measurements
    Bᴮ = quat2rot(x.q)' * Bᴵ;  # Transpose to convert from inertial → body

    ## Generate Jacobian H (mag_measurement wrt state) 

    B̂ᴮ = hat(Bᴮ)  # Skew-symmetric matrix 

    H = (calibrate_diodes) ? zeros(T, 3, 6 + 3 * N) :  # Jacobian Matrix  H = [∂θ ∂β ∂C ∂α ∂ϵ]
                             zeros(T, 3, 6)            # Jacobian Matrix  H = [∂θ ∂β] 

    H[1:3, 1:3] .= B̂ᴮ    
    # None of the other states directly affect measured mag field and are left as zeros

    return Bᴮ, H
end;

# Clean and add new testss
# Test - Our H seems to be half of forward diff? (as does OG)
# Best way to deal with 'corner' in H when in unlit? Make it zero? Real small?
# SUN INERTIAL is UNSCALED NOW!!
# Should the Is be unit?
# SUPER slow; with albedo(1) its 5.5ms /40alloc/416KiB, without it is 1.4μs / 30 / 4KiB; ~1ms comes from earth_albedo()
function current_measurement(x::SAT_STATE{T}, diodes::DIODES{N, T}, sᴵ::SVector{3, T}, pos::SVector{3, T}, alb::ALBEDO; 
                                E_am₀ = 1366.9, use_albedo = true, calibrate_diodes::Bool = true) where {N, T}
    # SUN INERTIAL is UNSCALED NOW!!
    C, α, ϵ = diodes.calib_values, diodes.azi_angles, diodes.elev_angles

    # Generate expected measurements
    sᴵ_unit = sᴵ / norm(sᴵ)  # Unit sun vector in inertial frame
    sᴮ = quat2rot(x.q)' * (sᴵ_unit);  # Transpose to convert from inertial → body
    n = SMatrix{N, 3, Float64, N * 3}([cos.(ϵ).*cos.(α) cos.(ϵ).*sin.(α) sin.(ϵ)]);  # Surface normals for each photodiode
    currents = C .* (n * sᴮ)     # Expected current measurements
    
    albedos  = zeros(T, N)
    if use_albedo
        albedo_matrix = earth_albedo(pos, sᴵ, alb.refl.data)
        for i = 1:N 
            diode_albedo = compute_diode_albedo(albedo_matrix, alb.cell_centers_ecef, n[i, :], pos)  
            albedos[i] = (C[i] / E_am₀) * diode_albedo 
        end
    end

    # Account for eclipses in I, H (a bit problematic for H)
    Is = zeros(T, N)  # This is to deal with the pesky Static problems
    for i = 1:N 
        Is[i] = currents[i] + albedos[i]
    end

    # # Generate Jacobian H (current measurement wrt state)
    ŝᴮ = hat(sᴮ);          # Skew-symmetric 
    ∂θ = (C .* n) * ŝᴮ     # How slight changes in orientation affect expected current field (∂I/∂θ)
    ∂β = zeros(T, N, 3)    # How slight changes in gyro bias affect expected current field (∂I/∂β)
    

    # H = [∂θ ∂β ∂C ∂α ∂ϵ]  :  [∂θ ∂β]
    H = (calibrate_diodes) ? zeros(T, 6, 6 + 3 * N)  : zeros(T, 6, 6)  
    H[:, 1:3] .= ∂θ
    H[:, 4:6] .= ∂β    # Should be just zeros...


    if calibrate_diodes
        # Jacobian Matrix  H = [∂θ ∂β ∂C ∂α ∂ϵ]
        ndα = [(-cos.(ϵ).*sin.(α)) ( cos.(ϵ).*cos.(α)) zeros(size(α))];
        ndϵ = [(-sin.(ϵ).*cos.(α)) (-sin.(ϵ).*sin.(α)) cos.(ϵ)]; # (With negative middle term, differing from the paper)
        
        ∂C = n * sᴮ;            # [i,]
        ∂α = C .* (ndα * sᴮ);   # [i,]
        ∂ϵ = C .* (ndϵ * sᴮ);   # [i,]  

        H[:, 7:end] .= [Diagonal(∂C) Diagonal(∂α) Diagonal(∂ϵ)]
    end

    Is[Is .< 1e-8]   .= 0.0   # Photodiodes don't generate negative current
    H[Is .≤  0.0, :] .= 0.0   #   ^ to match the above 

    return (Is / norm(Is)), H
end;



function update(x::SAT_STATE{T}, diodes::DIODES{N, T}, L::Matrix{T}, z::SVector{N₂, T}; calibrate_diodes::Bool = true) where {N, N₂, T}

    """ Uses kalman gain and innovation to update the best estimate for the satellite state """
    
    δx = L * z 

    # Muliplicative portion (orientation) - use axis-angle again
    δϕ = δx[1:3]  
    θ  = norm(δϕ)  # Magnitude of rotation 

    r  = (θ > 0.0) ? (δϕ / θ)  : [0, 0, 0]  # Unit axis of rotation (dealing with 0 rotation)
    δq = SVector{4, T}([cos(θ / 2); r * sin(θ / 2)]) # axis-angle to scalar-first unit quaternion 
    q⁺ = qmult(x.q, δq);  # Rotate a little 

    # Additive portion
    β⁺ = x.β + δx[4:6]
    
    # Update the calibration states if desired
    i₀ = 7; i₁ = i₀ + N - 1
    C⁺ = (calibrate_diodes) ? diodes.calib_values + δx[i₀:i₁] : diodes.calib_values

    i₀ = i₁ + 1; i₁ = i₀ + N - 1
    α⁺ = (calibrate_diodes) ? diodes.azi_angles + δx[i₀:i₁] : diodes.azi_angles

    i₀ = i₁ + 1; i₁ = i₀ + N - 1
    ϵ⁺ = (calibrate_diodes) ? diodes.elev_angles + δx[i₀:i₁] : diodes.elev_angles 


    x⁺ = SAT_STATE(q⁺, β⁺)
    diodes⁺ = DIODES(C⁺, α⁺, ϵ⁺)

    return x⁺, diodes⁺
end







function compute_diode_albedo(data::Matrix{T}, cell_centers_ecef::Array{T, 3}, surface_normal::SVector{3, T}, sat_pos::SVector{3, T}) where {T}

    Nlat, Nlon = size(data)

    diode_albedo = 0.0
    r_g = zeros(T, 3)

    for r = 1:Nlat
        for c = 1:Nlon

            if data[r,c] != 0
                r_g .= view(cell_centers_ecef, r, c, :) .- sat_pos
                r_g .= r_g ./ norm(r_g)  # Make unit

                cell_albedo = (data[r,c] * dot(surface_normal, r_g))

                if cell_albedo > 0.0    # Can't be negative
                    diode_albedo += cell_albedo 
                end
            end
        end
    end
    
    return diode_albedo
end

function chol(M)
    # return cholesky(Symmetric(M)).U 
    return cholesky(Hermitian(M)).U
end

function qrᵣ(M)
    return qr(M).R 
end
