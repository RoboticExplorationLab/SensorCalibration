# [src/MissionSim/Estimator/magnetometer_calibration.jl]

""" To Do:
      Add:
       - Way to deal with failures

    - Adjust main so this is only called when data is full
    - Remove the check that ensures the angles are right
    - Do I actually 𝑤𝑎𝑛𝑡 to force them all to call/return the same way? Returning data is useless
    - Gauss-Newton is pretty slow but is run only once, so...
    (None of this section needs to be fast, except maybe update!)
    - Do we actually 𝑛𝑒𝑒𝑑 nonlinear least squares...? The relation between meas and pred is affine ya?
"""



# Add alternate constructor with no initial vals? 
"""
    MAG_CALIBRATOR(N, bm, bp)

      Struct containing the necessary matrices to run least-squares and 
    estimate the magnetometer calibration values. Tracks the measured (B̃)
    and predicted (B) magnetic field vectors, as well as a matrix (A) that
    uses the calibration values (ρ) to map measured onto predicted (predicted are 
    treated as 'truth').

      B̃ = T * B + β  (for one instance), where B̃ is measured mag field, B 
      is predicted, T is the calibration matrix, and β is the bias. 

      B̃ = Aρ, where B̃ is a column vector of all measured mag field vectors, 
      ρ are the parameters that form T and β from above, and A is a series of
        
      Aᵢ =  [ Bₓ  0   0   0   0   0   1  0  0;
              0   Bₓ  0   By  0   0   0  1  0;
              0   0   Bₓ  0   By  Bz  0  0  1]

      ρ = [T₁₁ T₂₁ T₃₁ T₂₂ T₃₂ T₃₃ βx βy βz]
  
      (See Eq 1, 2, 3 in the referenced *Attitude-Independent Magnetometer Calibration with Time-Varying Bias* in README)

      Note that, when initialized, this struct allocates all the necessary 
    space in advance and then fills it up, keeping track of the current location 
    with 'idx', which is a vector of size 1 so I can update it in the immutable 
    struct. 
        
"""
struct MAG_CALIBRATOR{T} 
    A::Matrix{T}
    B_meas::Vector{T}   # Measured magnetic field vectors, stacked in a column
    B_pred::Vector{T}   # Predicted magnetic field vectors, stacked in a column
    N::Int              # Max number of samples

    # Index for most recent sample point (b/c we fill this up over time)
    idx::Vector{Int}    # This is a single-element vector so it can be updated

    function MAG_CALIBRATOR(N::Int, bm::Vector{T}, bp::Vector{T}) where {T}
        """ Primary Constructor """
        B_meas = zeros(T, 3 * N)
        B_meas[1:3] .= bm 

        B_pred = zeros(T, 3 * N)
        B_pred[1:3] .= bp 

        A = zeros(T, 3 * N, 9)  # Because we are estimating 9 parameters
        # Arrange the equations into a 3-row section of the matrix
        A[1:3, :] .= [(bp[1] * I(3))  (bp[2] * I(3))[:, 2:3]  (bp[3] * I(3))[:, 3]  I(3)]

        idx = [1]

        new{T}(A, B_meas, B_pred, N, idx)
    end
end

"""
    update!(mag_cal, bm, bp)

      Adds a sample point to the matrices in mag_cal by 
    updating the appropriate rows. Also increments the 
    index. Note that idx is a vector so it can be updated.

    Throws a warning if the matrices have already been filled
    
    Arguments:
      - mag_cal:  Struct containing the matrices to be updated         |  MAG_CALIBRATOR
      - bm:       Measured B-vector to be added to mag_cal             |  [3,]
      - bp:       Predicted B-vector to be added to mag_cal            |  [3,]

    Returns:
      - idx:      Index of newly added sample in the mag_cal struct    |  Int
"""
function update!(mag_cal::MAG_CALIBRATOR{T}, bm::Vector{T}, bp::Vector{T}) where {T} 

    i = mag_cal.idx[1] + 1  # Index for new sample (remember that it is a "vector")

    if i > mag_cal.N
        @warn "Mag Calibrator is already full!"
        return mag_cal.N
    end

    mag_cal.B_meas[(3 * i - 2) : (3 * i)] .= bm 
    mag_cal.B_pred[(3 * i - 2) : (3 * i)] .= bp

    # row = [ bp[1]  0.0    0.0    0.0    0.0    0.0    1.0   0.0   0.0;
    #         0.0    bp[1]  0.0    bp[2]  0.0    0.0    0.0   1.0   0.0;
    #         0.0    0.0    bp[1]  0.0    bp[2]  bp[3]  0.0   0.0   1.0]

    row = zeros(T, 3, 9)
    row[1, 1] = bp[1]
    row[2, 2] = bp[1] 
    row[3, 3] = bp[1]
    row[2, 4] = bp[2]
    row[3, 5] = bp[2]
    row[3, 6] = bp[3]
    row[1, 7] = 1
    row[2, 8] = 1 
    row[3, 9] = 1

    mag_cal.A[(3 * i - 2):(3 * i), :] .= row

    mag_cal.idx[1] += 1 # This line is why idx needs to be a vector in the immutable struct 
end


"""
    estimate(sat, data)

      Estimates the 9 magnetometer calibration values that form the
    lower-triangular calibration matrix T and bias vector β, so that 

                T = [p1  0   0;      β = [p7;
                     p2  p4  0;           p8;
                     p3  p5  p6]          p9]

      This is done by using the provided measured and predicted magnetic 
    field vectors stored in data to first perform linear least squares, 
    and then to improve upon that with nonlinear least squares using 
    Gauss-Newton. 

    NOTE Currently does NOT include time-varying current-induced bias (assume we turn off everything while calibrating) 
        -> Check out initial commit for version with currents included (or the separate folder has it in the past, too)

    Arguments:
      - sat:  Satellite struct that is updated (by copying) with the estimated         |  SATELLITE
                   magnetometer calibration values 
      - data: Struct containing the measured and predicted magnetic field              |  MAG_CALIBRATOR
                   vectors in body frame
"""
function estimate(sat::SATELLITE, data::MAG_CALIBRATOR{T}) where {T}

    # Ensure we have enough data first
    if data.idx[1] < data.N
        @error "Ah! I am not ready for this in mag cal yet!"

    else 
        init_guess = data.A \ data.B_meas       # Linear least squares to get initial guess 
        params = gauss_newton(init_guess, data) # Non-linear least squares to refine 

        # Convert the params into a calibration matrix and bias 
        T̂, β̂  = vec_to_matrix_bias(params)  

        â, b̂, ĉ, ρ̂ , λ̂ , ϕ̂ , β̂x, β̂y, β̂z = extract_elements(T̂, β̂ )

        # TEMPORARY (?) - Check that they are valid
        if (abs(ρ̂ ) > pi/3) || (abs(λ̂ ) > pi/3) || (abs(ϕ̂ ) > pi / 3)  # BAD 
            @error "Error with major non-orthogonal angles!"
        elseif (â < 0.0) || (b̂ < 0.0) || (ĉ < 0.0)
            @error "Error with negative scale values!"
        end

        # Update Satellite Estimate 
        mag_est = MAGNETOMETER(SVector{3, T}(â, b̂, ĉ),
                               SVector{3, T}(ρ̂ , λ̂ , ϕ̂ ),
                               SVector{3, T}(β̂x, β̂y, β̂z) )
                            
        # Because structs are immutable, we just make a new one and copy the relevant stuff over 
        sat_est = SATELLITE(sat.J, mag_est, sat.diodes) 

        return sat_est, data
    end
end

"""
    gauss_newton(x0, data; max_iters = 50, ls_iters = 20)

      Gauss-Newton nonlinear least squares method for batch estimation,
    only slightly modified from code provided by Kevin. Iteratively 
    updates an initial guess for the calibration parameters to minimize
    a residual function, so that 

        x⁺ = x⁻ - 𝐉⁻¹ r(x⁻)
    
    where the jacobian 𝐉 is calculated using ForwardDiff. Also utilizes 
    a line search to help speed up the process. 

    Arguments:
      - x0:   Initial guess for the parameters to be estimated                       |  Vector ([9,] here)
      - data: Measured and predicted mag vectors needed for the residual function    |  MAG_CALIBRATOR 

    Returns:
      - x:  Final guess for the parameters to be estimated                           |  [9,]                
"""
function gauss_newton(x0, data::MAG_CALIBRATOR{T}; max_iters = 50, ls_iters = 20) where {T}

    x = copy(x0) # copy initial guess

    Ds = 0.0
    v = zeros(length(x))
    _residual(x) = residual(x, data)  # Wrapper function for ForwardDiff

    for i = 1:max_iters

        J = ForwardDiff.jacobian(_residual, x)  # ∂r/∂x

        r = _residual(x)

        v = -J\r # solve for Gauss-Newton step (direct, or indirect)

        S_k = dot(r,r) # calculate current cost
 
        α = 1.0 # step size (learning rate)

        # run a simple line search
        for ii = 1:ls_iters
            x_new = x + α*v
            S_new = norm(_residual(x_new))^2

            # this could be updated for strong frank-wolfe conditions
            if S_new < S_k
                x = copy(x_new)
                Ds = S_k - S_new
                break
            else
                α /= 2
            end

            if ii == 20
                # @warn "line search failed"
                Ds = 0
            end
        end

        # depending on problems caling, termination criteria should be updated
        if Ds < 1e-5
            break
        end
    end
    return x
end

# Kinda slow but not used often
"""
    residual(x, data)

      Residual function to be minimized by Gauss-Newton, where the 
    cost function to be minimized is the difference between the magnitudes
    of the predicted mag field and the measured mag field once it has been
    corrected using the provided parameters, so that 

        minₓ  J = 0.5 * (B² - f(B, I, x))ᵀ(B² - f(B, I, x)) 

    Note that we are trying to find the KKT Stationarity condition 
    for this cost (∇J = 0.0) in order to find the local minimum.

      Because this relies only on magnitude, it is independent of satellite
    attitude. 

    Arguments:
      -x:     Parameters used to form the calibration matrix to correct the        |  [9,]
                  measured magnetic field vector (in body frame)
      -data:  Struct containing the measured and predicted magnetic field          |  MAG_CALIBRATOR
                  vectors.

    Returns:
      -r:   Vector of residuals for each sample                                    |  [N,]
"""
function residual(x, data::MAG_CALIBRATOR{T}) where {T} 

    function f(B_meas_i::SubArray, p::Vector) 
        """ 
              Helper function that takes a measured mag field vector and corrects 
            it using the provided calibration parameters 
        """
        T̂, β̂  = vec_to_matrix_bias(p)
        B = T̂ \ (B_meas_i - β̂ )
        return (B[1]^2) + (B[2]^2) + (B[3]^2)
    end

    function get_index(vec::Vector{T}, i::Int) where {T} 
        """ Extracts the three elements from a 1D vector corresponding to the given index"""
        return @view vec[ (3 * i - 2) : (3 * i) ]
    end

    N = data.N
    r = zeros(eltype(x), N)  # ForwardDiff friendly for .jacobian
    
    for i = 1:N
        B_meas_i = get_index(data.B_meas, i)
        B_pred_i = get_index(data.B_pred, i)

        mag_B_expected = (B_pred_i[1]^2) + (B_pred_i[2]^2) + (B_pred_i[3]^2)
        
        r[i] = mag_B_expected - f(B_meas_i, x)
    end

    return r
end

"""
    vec_to_matrix_bias(p) 

      Reshapes the 9-element parameter vector into a lower triangular
    [3 × 3] calibration matrix T and bias vector β. Accounts for sign 
    ambiguities in the calibration matrix by assuming that the non-
    orthogonality angles are less than π/2.

    Arguments:
      - p:  Calibration values determined by least squares to correct the       |  [9,]
                measured magnetic field vectors 
    
    Returns:
      - T:  Lower-triangular calibration matrix corresponding to p              |  [3, 3]
      - β:  Vector of magnetometer bias values corresponding to p               |  [3,]
"""
function vec_to_matrix_bias(p::Vector{ET}) where {ET}

    T = zeros(ET, 3, 3)

    # T = [p[1]   0       0;
    #      p[2]   p[4]    0;
    #      p[3]   p[5]    p[6]];      # Calibration matrix

    β = p[7:9]

    # |λ| or |ϕ| would have to at least π/2 for T[3, 3] < 0.0, which is a major problem anyway
    T[3, 3] = (p[6] ≥ 0.0) ? p[6]  :  -p[6]

    # |ρ| would have to be at least π/2
    if p[4] ≥ 0.0
        T[2, 2] = p[4] 
        T[3, 2] = p[5]
    else 
        T[2, 2] = -p[4] 
        T[3, 2] = -p[5]
    end

    # Scale factors (like a) cant be negative
    if p[1] ≥ 0.0
        T[1, 1] = p[1] 
        T[2, 1] = p[2]
        T[3, 1] = p[3]
    else
        T[1, 1] = -p[1] 
        T[2, 1] = -p[2]
        T[3, 1] = -p[3]
    end

    return T, β
end


"""
    extract_elements(T, β)

      Extracts the calibration elements from a calibration matrix T, 
    which consists of the scale factor along each axis [a, b, c] and 
    the non-orthogonality angles [ρ, λ, ϕ], as well as the bias values
    along each axis. T is lower-triangular and is formed as follows:

            T = [a           0.0              0.0;
                 bsin(ρ)     bcos(ρ)          0.0;
                 csin(λ)     ccos(λ)sin(ϕ)    ccos(λ)cos(ϕ)]

    Arguments:
      - T:  Lower-triangular calibration matrix corresponding to p              |  [3, 3]
      - β:  Vector of magnetometer bias values corresponding to p               |  [3,]
      
    Returns:
      - a:  Magnetometer scale factor for the X axis                            |  Scalar
      - b:  Magnetometer scale factor for the Y axis                            |  Scalar
      - c:  Magnetometer scale factor for the Z axis                            |  Scalar
      - ρ:  Non-orthogonality angle between X and Y axes                        |  Scalar
                (i.e., the angular distance from the 'true' Y to measured Y)
      - λ:  Non-orthogonality angle between X and Z axes                        |  Scalar
                (i.e., the angular distance from the 'true' Z to measured Z)
      - ϕ:  Non-orthogonality angle between Y and Z axes                        |  Scalar
                (i.e., the angular distance from the 'true' Z to measured Z)
      - βx: Magnetometer bias along the X axis                                  |  Scalar
      - βy: Magnetometer bias along the Y axis                                  |  Scalar
      - βz: Magnetometer bias along the Z axis                                  |  Scalar

"""
function extract_elements(T, β)

    a = T[1, 1]

    b = sqrt( (T[2,1]^2) + (T[2,2]^2) )  # (bsin)^2 + (bcos)^2 = b^2
    ρ = atan(T[2,1], T[2,2])             # atan₂ to maintain sign

    # Similarly, ...
    c = sqrt( (T[3,1]^2) + (T[3,2]^2) + (T[3,3]^2) )
    ϕ = atan(T[3,2] / T[3,3])
    λ = atan(  sign( T[3,1]) * sqrt( (T[3,1]^2) ),  
               sign((T[3,2]^2) + (T[3,3]^2)) * sqrt( (T[3,2]^2) + (T[3,3]^2) ) ) 

    βx, βy, βz = β 

    return a, b, c, ρ, λ, ϕ, βx, βy, βz
end 

