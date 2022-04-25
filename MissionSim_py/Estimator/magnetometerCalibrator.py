import numpy as np

# TODO add to ulab/numpy helper functions 
def sign(v):
    if (v == 0.0):
        return 0.0
    elif v < 0.0:
        return -1 
    else: 
        return 1

#### ForwardDiff - JACOBIAN 
# NOT FIXED YET!
"""
def gauss_newton(x):
    # "" Gauss-Newton for batch estimation (From Kevin) ""

    Ds = 0.0
    v = np.zeros(x.shape[0])

    # run Gauss-Newton for 100 iterations max
    for i in range(50):

        # ∂r/∂x
        _res(x) = residual(x) # residual(x, data)
        # J = ForwardDiff.jacobian(lambda x = residual(x), x) ###################

        # calculate residual at x
        r = _res(x)

        # solve for Gauss-Newton step (direct, or indirect)
        v = -J / r # -J\r
        # lsqr!(v,-J,r)

        # calculate current cost
        S_k = np.dot(r,r)
        # @show S_k

        # step size (learning rate)
        α = 1.0

        # run a simple line search
        for ii in range(25): 
            x_new = x + α*v
            S_new= norm(_res(x_new))^2

            # this could be updated for strong frank-wolfe conditions
            if S_new < S_k:
                x = copy(x_new)
                Ds = S_k - S_new
                # @show ii
                break
            else:
                α /= 2

            if ii == 25:
                # @warn "line search failed"
                Ds = 0

        # depending on problems caling, termination criteria should be updated
        if Ds < 1e-5:
            break

    return x
"""


print("WARNING - GAUSS_NEWTON NOT IMPLEMENTED IN MAGNETOMETER CALIBRATOR!")
class MagnetometerCalibrator():

    ### MEMBER CONSTANTS ### 
    orbit_length = None 
    max_idx      = None

    ### MEMBER VARIABLES ###
    measurement_history = None
    prediction_history  = None
    A = None 
    
    ready_to_run = False
    has_run      = False 

    current_idx = 0

    def __init__(self, max_idx = int((55780/600)) ):
        self.max_idx = max_idx 

        self.measurement_history = np.zeros(3 * self.max_idx)
        self.prediction_history  = np.zeros(3 * self.max_idx)
        self.A  = np.zeros([3 * self.max_idx, 9])

    # TODO pass in sat so it can be updated? -> run_batch_least...
    def update_estimate(self, meas, pred, estimate_flag = None):
        idx0 = self.current_idx * 3 
        idxf = idx0 + 3
        self.current_idx = (self.current_idx + 1) % self.max_idx

        self.measurement_history[idx0:idxf] = meas 
        self.prediction_history[idx0:idxf]  = pred

        I = np.eye(3)
        row = np.hstack([ (pred[0] * I),
                          (pred[1] * I)[:, 1:3],
                          (pred[2] * I)[:, 2][:, None],
                          I ])

        self.A[idx0:idxf, :] = row

        if self.current_idx == 0: # Means it has been completely filled up and is looping through again
            self.ready_to_run = True 

        if estimate_flag is None: # Allow flag to override everything
            estimate_flag = self.ready_to_run

        if estimate_flag:
            return self.run_batch_least_squares()  # Updated satellite
        else:
            return None 

    ####################
    # HELPER FUNCTIONS #   (Those not called directly)
    ####################

    # TODO check if we can smooth this one out
    def matrix_to_parameters(self, T):
        a = T[0,0]

        b = np.sqrt( (T[1, 0]**2) + (T[1, 1]**2) )
        rho = np.arctan2(T[1, 0], T[1, 1])

        c = np.sqrt( (T[2,0]**2) + (T[2,1]**2) + (T[2,2]**2) )
        phi = np.arctan( T[2,1] / T[2,2])

        # Is there a reason we sqrt and square here, rather than abs()? And no dot?
        lamb = np.arctan2( sign(T[2,0]) * np.sqrt( (T[2,0]**2) ), \
                        sign( (T[2,1]**2) + (T[2,2]**2) ) * np.sqrt( (T[2,1]**2) + (T[2,2]**2) ) )
        
        return a, b, c, rho, lamb, phi

    def parameters_to_matrix_bias(self, p):
        T = np.array([
            [p[0], 0,    0],
            [p[1], p[3], 0],
            [p[2], p[4], p[5]]
        ])

        bias = p[6:9]

        # Account for sign ambiguities in the calibration matrix
        if T[2,2] < 0:  # Scale factor (c) can't be negative
            T[2,2] = -T[2,2]
        
        if T[1,1] < 0:  # Scale factor (b) can't be negative
            T[:, 1] = -T[:, 1]

        if T[0,0] < 0:  # Scale factor (a) can't be negative
            T[:, 0] = -T[:, 0]

        return T, bias

    # TODO rename...?
    def f(self, bm, p):
        """ Undoes the effect of the calibration matrix/bias vectors """
        T, bias = self.parameters_to_matrix_bias(p)

        B = np.linalg.inv(T) @ (bm - bias)
        B_squared = (B[0]**2) + (B[1]**2) + (B[2]**2)

        return B_squared 

    def residual(self, x):
        """ 
            residual vector for Gauss-Newton. rᵀr = MAP cost function
                Note that x = parameters 
                Meas = [(B_meas, B_pred) x T]
                Loss Function: 
                J = 0.5*(B^2 - f(B,I,x))^T(B^2 - f(B,I,x))
        """
        N = int(mag_field_meas_hist.shape[0] / 3)  # stored as [x₁ y₁ z₁ x₂ y₂ .....] so we need to divide by 3 
        r = np.zeros(N)
        
        for i in range(N):
            B_meas = mag_field_meas_hist[(3*i):(3*(i+1))]
            B_exp_squared = (mag_field_pred_hist[(3*i)]^2) + (mag_field_pred_hist[(3*i + 1)]^2) + (mag_field_pred_hist[(3*i + 2)]^2) # Should be [M x 1] -> Should M > 1?

            J = f(B_meas, x) - B_exp_squared
            r[i] =  J 

        return reshape(r, length(r))

    def run_batch_least_squares(self):
        self.has_run = True 

        # Solve using QR Decomposition (TODO UNTESTED!)
        q, r = np.linalg.qr(self.A)
        params = np.linalg.inv(r) @ q.T @ self.measurement_history

        # TODO Run through gauss_newton once able to 
        # params = self.gauss_newton(params)

        T_est, bias_est = self.parameters_to_matrix_bias(params)

        a, b, c, rho, lamb, phi = self.matrix_to_parameters(T_est)

        if ((np.abs(rho) > np.pi/3) or (np.abs(lamb) > np.pi/3) or (np.abs(phi) > np.pi/3)):
            print("Error with non-ortho angles!")

        self.a, self.b, self.c                = a, b, c 
        self.rho, self.lamb, self.phi         = rho, lamb, phi 
        self.bias_x, self.bias_y, self.bias_z = bias_est

        updated_params = {'a': a, 'b': b, 'c': c, 'rho': rho, 'lamb': lamb, 
                            'phi': phi, 'bx': bias_est[0], 'by': bias_est[1], 'bz': bias_est[2]}
        return updated_params







