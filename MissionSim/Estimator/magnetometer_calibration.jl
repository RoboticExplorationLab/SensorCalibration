####################################################################
#               MAGNETOMETER CALIBRATION                           #
####################################################################

using JLD2, Infiltrator
using Symbolics


function save_global_variables(sat_truth, sat_init_est)
    """
        Temporary function used to save orbit data to test various downsampling rates
    """
    @save "orbit_data_for_mag_calib2.jl" mag_field_meas_hist mag_field_pred_hist A sat_truth sat_init_est
    println("Saved Orbit data!")
end

function generate_mag_calib_matrix(sat::SATELLITE)
    """ Generates the calibration matrix that alters the measured magnetic field vector in body frame """
    a, b, c = sat.magnetometer.scale_factors
    ρ, λ, ϕ = sat.magnetometer.non_ortho_angles

    T = [a        0.0              0.0;
        b*sin(ρ)  b*cos(ρ)         0.0;
        c*sin(λ)  c*sin(ϕ)*cos(λ)  c*cos(ϕ)*cos(λ)]

    return T
end



# Int(round(2 * orbit_period(oe0[1]))/_dt)
# @info "Relies on there being no change to dt and stuff (mag cal)"  # Use an init() function? Better if no globals...
@info "Speed up magnetometer calibration"
# max_idx = Int(round(55780/600))
# cur_idx = 0
mag_field_meas_hist = 0# @MArray zeros(max_idx*3)
mag_field_pred_hist = 0#@MArray zeros(max_idx*3)
A = 0 #@MArray zeros(3*max_idx, 9)
has_run = false

struct MAG_CALIB 
    mag_field_meas
    mag_field_pred
end

function check_if_run()
    """ Simple "getter" function """
    return has_run
end

# Currently does NOT include time-varying current-induced bias (assume we turn off everything while calibrating) -> Check out initial commit for version with currents included 
function estimate_vals(sat::SATELLITE, data::MAG_CALIB, estimate_flag::Bool)
    """
        Accumulates predicted and measured magnetic field vectors. Once sufficiant data are gathered,
            runs a Gauss-Newton solver to estimate the magnetometer calibration matrix

        Arguments:
            - sat:  Struct containing current estimate of satellite values      |  SATELLITE 
            - data: Struct containing predicted/measured mag field vectors      |  MAG_CALIB
            - estimate_flag: Boolean flag signalling when to run Gauss-Newton   |  Bool

        Returns:
            - sat:  Updated satellite struct containing new estimates (if run)  |  SATELLITE
            - data: Unaltered data struct                                       |  MAG_CALIB
    """  

    
    if isempty(size(mag_field_meas_hist))  # Initialize   
        global mag_field_meas_hist = data.mag_field_meas[:]    
        global mag_field_pred_hist = data.mag_field_pred[:]   
        global A = [(data.mag_field_pred[1]*I(3))       (data.mag_field_pred[2]*I(3))[:, 2:3]       (data.mag_field_pred[3]*I(3))[:, 3]    I(3)]
    else 
        global mag_field_meas_hist = [mag_field_meas_hist[:]; data.mag_field_meas[:]] 
        global mag_field_pred_hist = [mag_field_pred_hist[:]; data.mag_field_pred[:]] 
        new_row =  [(data.mag_field_pred[1]*I(3))       (data.mag_field_pred[2]*I(3))[:, 2:3]       (data.mag_field_pred[3]*I(3))[:, 3]    I(3)]

        global A = [A; new_row]
    end
    
    # idx = cur_idx % max_idx # what about when it == 0? it should be max_idx instead yeah?  (cur_idx % max_idx) + 1, cur_idx\_0  = 0 != 1
    # i₀ = (3*idx + 1)
    # i₁ = (3*idx + 3)
    # global mag_field_meas_hist[i₀:i₁] = data.mag_field_meas
    # global mag_field_pred_hist[i₀:i₁] = data.mag_field_pred
    # global A[i₀:i₁, :] .= [(data.mag_field_pred[1]*I(3))    (data.mag_field_pred[2]*I(3))[:, 2:3]    (data.mag_field_pred[3]*I(3))[:, 3]    I(3)]

    # global cur_idx += 1

    if estimate_flag # Only estimates once enough data has been gathered
        return run_gn(sat, data)
    else    
        return sat, data
    end
end

function initialize(data::MAG_CALIB)
    """ Resets system to prepare to calibrate again. Note that the data is returned unadjusted """

    global mag_field_meas_hist = @MArray zeros(max_idx*3)
    global mag_field_pred_hist = @MArray zeros(max_idx*3)
    global A = @MArray zeros(3*max_idx, 9)
    global has_run = false
    global cur_idx = 0

    return data
end


# @info "Gauss Newton portion commented out!"
function run_gn(sat, data)
    """
        Runs a Gauss-Newton solver on accumulated magnetic field vector data to estimate    
            the calibration matrix. Accounts for sign ambiguity 

        Arguments:
            - sat:  Struct containing current estimate of satellite values      |  SATELLITE 
            - data: Struct containing predicted/measured mag field vectors      |  MAG_CALIB

        Returns:
            - sat:  Updated satellite struct containing new estimates (if run)  |  SATELLITE
            - data: Unaltered data struct                                       |  MAG_CALIB
    """ 

    global has_run = true
    params = A \ mag_field_meas_hist 

    params = gauss_newton(params)

    mag_calib_matrix_est, β = parameters_to_matrix_bias(params)
    bx_est, by_est, bz_est = β[:]

    a_est, b_est, c_est, ρ_est, λ_est, ϕ_est = extract_parameters(mag_calib_matrix_est)

    if (abs(ρ_est) > pi/3) || (abs(λ_est) > pi/3) || (abs(ϕ_est) > pi/3)
        println("Error with major non-ortho angles!")
        finished = false
    else
        finished = true
    end


    # UPDATE SATELLITE ESTIMATES
    updated_magnetometer_est = MAGNETOMETER([a_est, b_est, c_est],     # NONE should be negative
                                            [ρ_est, λ_est, ϕ_est],     # NONE should be | | > pi/2 (or even close)
                                            [bx_est, by_est, bz_est] )
    sat.magnetometer = updated_magnetometer_est

    return sat, data 
end

#########################
function new_mag_calib()
    return MAG_CALIB(0.0, 0.0)
end

function extract_parameters(T)
    """ Extracts the calibration parameters from a matrix 
            [a, b, c] (scale factors) and [ρ, λ, ϕ] (non-orthogonality angles) """

    a = T[1,1] # Easy 

    b = sqrt((T[2,1]^2) + (T[2,2]^2)) # (bsin)^2 + (bcos)^2 = b^2
    ρ = atan(T[2,1], T[2,2]) # sin/cos to maintain signs

    c = sqrt((T[3,1]^2) + (T[3,2]^2) + (T[3,3]^2))
    ϕ = atan(T[3,2] / T[3,3])
    λ = atan(  sign(T[3,1]) * sqrt( (T[3,1]^2) ),  
                sign((T[3,2]^2) + (T[3,3]^2)) * sqrt( (T[3,2]^2) + (T[3,3]^2) ) ) # Not positve this portion of the sign is actually useful

    return a, b, c, ρ, λ, ϕ
end

function parameters_to_matrix_bias(p)
    """ Reshapes Gauss-Newton output into a calibration matrix and biases """

    T = [p[1]   0       0;
            p[2]   p[4]    0;
            p[3]   p[5]    p[6]];      # Calibration matrix

    β = p[7:9];     # Bias vector 


    # Account for sign ambiguities in the calibration matrix
    if T[3,3] < 0
        T[3,3] = -T[3,3]
    end
    if T[2,2] < 0 
        T[2,2] = -T[2,2]
        T[3,2] = -T[3,2]
    end
    if T[1,1] < 0
        T[:, 1] = -T[:, 1]
    end

    return T, β
end

function f(bm, p)
    """ Undoes the effect of the calibration matrix/bias vectors """
    T_hat, bias_hat =  parameters_to_matrix_bias(p)

    B = (T_hat^(-1))*(bm - bias_hat)
    B_squared = (B[1]^2) + (B[2]^2) + (B[3]^2)
    return B_squared 
end

@info "Unclear if using correct function for residual (J or dJ/dx?)"
function residual(x)
    """ residual vector for Gauss-Newton. rᵀr = MAP cost function
            Note that x = parameters 
            Meas = [(B_meas, B_pred) x T]
            Loss Function: 
                J = 0.5*(B^2 - f(B,I,x))^T (B^2 - f(B,I,x))
    """
    N = Int(size(mag_field_meas_hist, 1) / 3)  # stored as [x₁ y₁ z₁ x₂ y₂ .....] so we need to divide by 3 
    r = zeros(eltype(x), (N))
    
    for i = 1:N
        B_meas = mag_field_meas_hist[(3*i - 2):(3*i)]
        B_exp_squared = (mag_field_pred_hist[(3*i - 2)]^2) + (mag_field_pred_hist[(3*i - 1)]^2) + (mag_field_pred_hist[(3*i)]^2) # Should be [M x 1] -> Should M > 1?

        J = f(B_meas, x) - B_exp_squared
        r[i] =  J 
    end

    return reshape(r, length(r))
end

function gauss_newton(x0)
    """ Gauss-Newton for batch estimation (From Kevin) """

    # copy initial guess
    x = copy(x0)

    # create sparse jacobian
    # J = spzeros(nx*(T-1) + m*(T),nx*T)

    Ds = 0.0
    v = zeros(length(x))

    # run Gauss-Newton for 100 iterations max
    for i = 1:50

        # ∂r/∂x
        _res(x) = residual(x) # residual(x, data)
        J = ForwardDiff.jacobian(_res, x)


        # calculate residual at x
        r = _res(x)

        # solve for Gauss-Newton step (direct, or indirect)
        v = -J\r
        # lsqr!(v,-J,r)

        # calculate current cost
        S_k = dot(r,r)
        # @show S_k

        # step size (learning rate)
        α = 1.0

        # run a simple line search
        for ii = 1:25
            x_new = x + α*v
            S_new= norm(_res(x_new))^2

            # this could be updated for strong frank-wolfe conditions
            if S_new < S_k
                x = copy(x_new)
                Ds = S_k - S_new
                # @show ii
                break
            else
                α /= 2
            end
            if ii == 25
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
