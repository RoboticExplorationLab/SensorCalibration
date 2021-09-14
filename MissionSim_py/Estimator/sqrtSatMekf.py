import numpy as np
from satelliteDynamics import * 
from rotationFunctions import * 
from earthAlbedo import EarthAlbedo 
from sqrtMekf import SqrtMekf

# TODO List:
# - Move 'diag' func from here and child to support functions (~numpy)
# - Validate the noise constants being used
# - Verify that the Epoch.add_seconds works
# - Comment and clean up 


# TODO Could be moved to support functions
def diag(v): # Copy of numpy, because ulab doesnt have it, but also abbreviated (?)
    n = v.shape[0]
    result = np.zeros((n,n), v.dtype)
    result[:n].flat[0::n+1] = v 
    return result

class SqrtSatMekf(SqrtMekf):

    ### Member Constants ###
    _E_am0 = 1366.9
    
        # Noise constants (TODO Validate) 
    gyro_orient_sigma2 = (np.radians(0.8)**2)/(3600**3)
    gyro_bias_sigma2   = (np.radians(0.06)**2)/3600 
    diode_sigma2       = (1e-5)**2   # Arbitrary-ish
    mag_vec_sigma2     = np.radians(3.0)**2
    curr_sigma2        = 0.05**2     # Also arbitrary-ish 
    q_sigma2           = (10*np.pi/180)**2 # ^
    bias_sigma2        = (10*np.pi/180)**2 # ^

    #### Member Variables ###
        # State
    q      = None 
    bias   = None
    Pchol  = None 
    time   = None

    W      = None 
    V      = None

        # Sensors
    sun_inert = None 
    mag_inert = None 
    sun_body  = None 
    mag_body  = None
    diodes    = None
    gyro      = None
    gps       = None

        # Diode Calibration stuff (Constant once initialized)
    calibration_values = None
    azimuth_angles     = None
    elevation_angles   = None
 
        # System
    dt         = None
    num_diodes = None

    #####################
    # PRIMARY FUNCTIONS #   (Those called directly)
    #####################

    def __init__(self, c, a, e, system, alb, vals):
        """ 
            Initializes a new SqrtSatMekf object by calling the parent template (SqrtMekf) init function
                and then setting up the values needed for this specific satellite problem.
                Note that c, a, and e are CONSTANTS that have already been determined

            Arguments:
                - c: Array of calibration values for each photodiode                |  [i,]
                - a: Array of azimuth angles for each photodiode surface normal     |  [i,]
                - e: Array of elevation angles for each photodiode surface normal   |  [i,]
                - system: Dictionary of various system parameters                   |  [4,]
                        (_dt, _T, _max_sim_length, _num_diodes)   
                - alb: EarthAlbedo class used to estimate Earth's albedo            |  EarthAlbedo object
                - vals: Dictionary containing whatever member values need to be     |  Dict  [13 + ]
                        stored (including state & sensors)
                        (NOTE that spelling matters a lot for this dict!)
        """
        super().__init__()
        self.calibration_values = c 
        self.azimuth_angles     = a 
        self.elevation_angles   = e 

        self.dt = system._dt
        self.num_diodes = system._num_diodes

        self.alb = alb
        self.cell_centers_ecef = alb.get_albedo_cell_centers()

        # Initialize covariance matrix
        p = np.hstack([ (self.q_sigma2) * np.ones(3),
                        (self.bias_sigma2) * np.ones(3) ])
        P = diag(p)
        self.Pchol = self.chol(P)

        # Initialize noise matrices
        w = np.hstack([(self.gyro_orient_sigma2) * np.ones(3),
                       (self.gyro_bias_sigma2) * np.ones(3),
                       (self.diode_sigma2) * np.ones(3 * self.num_diodes) ])

        v = np.hstack([(self.mag_vec_sigma2)*np.ones(3),
                       (self.curr_sigma2)*np.ones(system._num_diodes) ])

        self.W = diag(w)
        self.V = diag(v) 

        self.update_system(vals) # State & Sensors

    def reset(self):
        pass

    def update_estimate(self, d, sat):
        """
            Runs a step of the MEKF to update estimates for state values.
                First updates member values, and then calls the mekf_step function.
                State estimates are updated here and on the satellite, which is returned.

            Arguments:
                - d: Dictionary containing member values to be updated      |  Dict 
                - sat: Current best estimate of the satellite parameters    |  Satellite Object

            Returns:
                - sat: Updated best estimate of the satellite parameters    |  Satellite Object
        """
        self.update_system(d)

        x_n, Pchol_n = self.mekf_step()


        self.Pchol = Pchol_n 
        self.q     = x_n[:4]
        self.bias  = x_n[4:7]

        sat.covariance = Pchol_n 
        sat.state      = x_n 

        return sat

    ####################
    # HELPER FUNCTIONS #   (Those not called directly)
    ####################

    def mekf_step(self):
        """
            Runs a step of the MEKF to predict next state and covariance matrix.
                First predicts next state and covariance, then uses the difference between what
                measurements would be expected with those predictions and what the measurements 
                actually were to update the prediction for next state and covariance.

            Returns:
                - x_next: Best estimate of next state                   |  [7,]
                            ( [q, bias] )
                - Pchol_next: Next covariance matrix (Cholesky factor)  |  [6 x 6]
        """
 
        x_p, Pchol_p = self.predict()

        L, C, z = self.innovate(x_p, Pchol_p)

        x_next, Pchol_next = self.update(x_p, Pchol_p, L, C, z)

        return x_next, Pchol_next

    def predict(self):
        """ Calls the standard "predict" function defined by parent class """
        x_p, A, P_p = super().predict(self.q, self.bias, self.Pchol, self.gyro, self.W, self.dt)
        return x_p, P_p

    def innovate(self, x_p, Pchol_p):
        """
            Innovation step of the MEKF that generates what the expected measurements
                would be based off of the predicted next state, and compares them to what
                the actual measurements were in order to estimate optimal Kalman gain

            Arguments:
                - x_p: Predicted next state                 |  [7,]
                        ( [q, bias] )
                - Pchol_p: Predicted covariance matrix      |  [6 x 6]
                        ( As a Cholesky factor )

            Returns: 
                - L: Kalman gain used to control how much the measurement   
                        affects the original prediction (large L means we 
                        trust the measurements more than the dynamics)

                - C: Sensitivity Matrix (Jacobian of measurements y wrt states x)       |  [3 + i x 6]
                - z: Innovation vector of differences between measured and predicted    |  [3 + i] 
                        measurement values 
        """
        ## Generate measurements
        sun_inert_unit = self.sun_inert / np.linalg.norm(self.sun_inert)
        mag_inert_unit = self.mag_inert / np.linalg.norm(self.mag_inert)
        mag_body_unit  = self.mag_body  / np.linalg.norm(self.mag_body )

            # Magnetometer 
        yp_mag, C_mag = self.mag_measurement(x_p, mag_inert_unit)
        z_mag = mag_body_unit - yp_mag 

            # Photodiodes 
        yp_cur, C_cur = self.current_measurement(x_p, sun_inert_unit)
        z_cur = self.diodes - yp_cur 

        C = np.vstack([ C_mag, C_cur ])
        z = np.hstack([ z_mag, z_cur ])

        ## Compute Kalman Gain
        L = super().innovate(C, z, self.V, Pchol_p)
        return L, C, z 

    def update(self, x_p, Pchol_p, L, C, z):
        """ Updates the predicted state and covariance using the parent's update method """ 
        return super().update(x_p, Pchol_p, L, C, z, self.V) # x_next, Pchol_next 

    # TODO Move "Add_seconds" to Epoch, not Filter
    def update_system(self, data):
        """
            Takes in data from a dictionary and uses it to update member values.
                Note that this method is highly sensitive to spelling errors so 
                verify that you use same parameter names as listed above
                (NOTE that that the time step is automatically updated, except
                    during initialization)

            Arguments:
                - Data: A dictionary containing variable names and updated values
        """
        if not (self.time == None): # Don't want to time step when initializing
            self.__dict__.update(data)
            self.time = self.add_seconds(self.dt)
        else:
            self.__dict__.update(data)
        
    def add_seconds(self, sec):
        """ Method for adding seconds to an Epoch """
        s = self.time.seconds
        d = self.time.days 
        ns = self.time.nanoseconds

        if sec < 1:
            ns += (sec * 1e9)
        else:
            s += sec

        if ns >= 1e9:
            ns -= 1e9 
            s += 1
        if s >= 86400:
            s -= 86400 
            d += 1

        return Epoch(d, s, ns, self.time.tsys)

    def mag_measurement(self, x_p, mag_inert_unit, 
                            c = None, a = None, e = None):

        c = self.calibration_values if np.any(c == None) else c 
        a = self.azimuth_angles     if np.any(a == None) else a 
        e = self.elevation_angles   if np.any(e == None) else e

        q, b = x_p[:4], x_p[4:7]

        bQi = dcm_from_q(q).T 
        B_body = bQi @ mag_inert_unit

        dtheta = hat(B_body)       # dy/dtheta
        dbias  = np.zeros([3, 3])  # dy/dbias = 0 as gyro bias doesn't affect mag field measurement

        H = np.hstack([dtheta, dbias])  # Sensitivity Matrix (Jacobian dy/dx)
        y = B_body 
        return y, H 

    def current_measurement(self, x_p, sun_inert_unit, 
                            c = None, a = None, e = None):

        c = self.calibration_values if np.any(c == None) else c 
        a = self.azimuth_angles     if np.any(a == None) else a 
        e = self.elevation_angles   if np.any(e == None) else e

        q, b = x_p[:4], x_p[4:7]
        bQi = dcm_from_q(q).T 
        
        s_body = bQi @ sun_inert_unit
        s_body_skew = hat(s_body)

        surf_norm = np.array([np.cos(e)*np.cos(a),
                                np.cos(e)*np.sin(a), 
                                np.sin(e)]).T     # [i x 3] 

        dtheta = (c[:, None] * surf_norm) @ s_body_skew # [i x 3] 
        dbias  = np.zeros([self.num_diodes, 3])                       # [i x 3]

        H = np.hstack([dtheta, dbias]) # Sensativity Matrix dy/dx       | [i x 6]
        I = c * (surf_norm @ s_body)   # Expected current measurements  | [i,]

        # Account for Earth's albedo
        s_inert_unscaled = sun_position(self.time) - self.gps
        albedo_matrix, _ = self.alb.albedo(self.gps, s_inert_unscaled, self.alb.data)
        diode_albedos = self.alb.get_diode_albedo(albedo_matrix, self.cell_centers_ecef, surf_norm, self.gps)
        diode_albedos = c * diode_albedos / self._E_am0
        I = I + diode_albedos 

        I[I < 0] = 0 
        H[I <= 0] = 0
        y = I 
        return y, H