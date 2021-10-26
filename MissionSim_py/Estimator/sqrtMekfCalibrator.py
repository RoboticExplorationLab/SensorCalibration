import numpy as np
from satelliteDynamics import * 
from rotationFunctions import * 
from earthAlbedo import EarthAlbedo 
from sqrtSatMekf import SqrtSatMekf

# TODO List:
# - Move 'diag' func from here and parent to support functions (~numpy)
# - Validate noise constants
# - Verify Ulab has h/vstack (or implement...)



# TODO Could be moved to support functions
def diag(v): # Copy of numpy, because ulab doesnt have it, but also abbreviated (?)
    n = v.shape[0]
    result = np.zeros((n,n), v.dtype)
    result[:n].flat[0::n+1] = v 
    return result



class SqrtMekfCalibrator(SqrtSatMekf):

    ### Member Constants ###
    calib_sigma2 = 0.2**2 
    azi_sigma2   = 2.0**2 # degrees    # TODO SHOULD THESE BE RADIANS? -> Probably
    elev_sigma2  = 1.0**2 # degrees

    ### Member Variables ###
    # Note that c, a, e are CHANGING

    #####################
    # PRIMARY FUNCTIONS #   (Those called directly)
    #####################

    def __init__( self, c, a, e, system, alb, vals):
        """ 
            Initializes a new SqrtMekfCalibrator object by calling the parent (SqrtSatMekf) init function
                and then augmenting the covariance matrix to include the calibration values

            Arguments:
                - c: Array of calibration values for each photodiode                |  [i,]
                - a: Array of azimuth angles for each photodiode surface normal     |  [i,]
                - e: Array of elevation angles for each photodiode surface normal   |  [i,]
                - system: Dictionary of various system parameters                   |  [4,]
                        (_dt, _T, _max_sim_length, _num_diodes)   
                - alb: EarthAlbedo class used to estimate Earth's albedo            |  EarthAlbedo object
                - vals: Dictionary containing whatever member values need to be     |  Dict
                        stored (e.g., sun/mag vectors, sensor data, etc.)
                        (NOTE that spelling matters a lot for this dict!)
        """
        super().__init__(c, a, e, system, alb, vals)

        # Append covariance of calibration terms 
        p = np.hstack([ (self.calib_sigma2) * np.ones(self.num_diodes),
                        (self.azi_sigma2)   * np.ones(self.num_diodes),
                        (self.elev_sigma2)  * np.ones(self.num_diodes) ])

        t1 = np.hstack([self.Pchol, np.zeros([6, 3 * self.num_diodes])])
        t2 = np.hstack([np.zeros([3 * self.num_diodes, 6]), self.chol(diag(p))])
        self.Pchol = np.vstack([t1, t2])
        
        # Augment process noise matrix 
        w_diodes = (self.diode_sigma2) * np.ones(3 * self.num_diodes)
        w1 = np.hstack([ self.W, np.zeros([6, 3*self.num_diodes]) ])
        w2 = np.hstack([ np.zeros([3*self.num_diodes, 6]), diag(w_diodes)])
        self.W = np.vstack( [w1, w2])
        


    def reset(self):
        pass

    def update_estimate(self, d, sat):
        """
            Runs a step of the MEKF to update estimates for state and calibration values.
                First updates member values, and then calls the parent mekf_step function.
                Parameter estimates are updated here and on the satellite, which is returned.

            Arguments:
                - d: Dictionary containing member values to be updated      |  Dict 
                - sat: Current best estimate of the satellite parameters    |  Satellite Object

            Returns:
                - sat: Updated best estimate of the satellite parameters    |  Satellite Object
        """
        super().update_system(d)

        x_n, Pchol_n = super().mekf_step()

        self.Pchol = Pchol_n 
        self.q     = x_n[:4]
        self.bias  = x_n[4:7]

        i = self.num_diodes
        self.calibration_values  = x_n[7:(7+i)]
        self.azimuth_angles      = x_n[(7+i):(7+2*i)]
        self.elevation_angles    = x_n[(7+2*i):]


        # This allows it to be run with PyCall, which doesn't seem to let me update Julia structs
        try:
            sat.covariance = Pchol_n
            sat.state      = x_n[:7]

            sat.diodes.calibration_values  = x_n[7:(7+i)]
            sat.diodes.azimuth_angles      = x_n[(7+i):(7+2*i)]
            sat.diodes.elevation_angles    = x_n[(7+2*i):]
        except:
            print("Could not update sat!")

        return sat

    ####################
    # HELPER FUNCTIONS #   (Those not called directly)
    ####################

    def predict(self):
        """         ---DOES NOT CALL PARENT's PREDICT METHOD!---
            Runs the "predict" step of the MEKF by using current state and calibration estimates
                predict the next state and and next covariance matrix.
                - Estimates the angular velocity 
                - Converts axis/angle representation to a quaternion to update orientation
                - Estimates the Jacobian of the state wrt the state (dx/dx)
                (Note that this always predicts no change in calibration values or bias, as expected)

            Returns:
                - x_p:  Predicted state values                  |  [7 + 3i]
                            ( [q, bias, c, a, e] )
                - P_p:  Predicted state covariance matrix       |  [6 + 3i x 6 + 3i]
        """
        gamma = self.gyro - self.bias  # Adjusted angular velocity 
        gamma_mag  = np.linalg.norm(gamma)
        gamma_skew = -hat(gamma) / gamma_mag 

        ## STATE (x_p)
        theta = gamma_mag * self.dt   # Angle of rotation of agent
        r = gamma / gamma_mag    # Axis of rotation (unit)
        p = np.hstack([r*np.sin(theta/2), np.cos(theta/2)])  # Quaternion representation

        q_p = qmult(self.q, p)       # Predicted quaternion

        x_p = np.hstack([q_p, self.bias])    # Predicted state

        ## JACOBIAN (A)
        # Rodrigues
        R = np.eye(3) + (gamma_skew * np.sin(theta)) + (gamma_skew @ gamma_skew) * (1 - np.cos(theta))

        # Jacobian (dx/dx)
        A = np.vstack([ np.hstack([R, -self.dt*np.eye(3)]),
                        np.hstack([np.zeros([3,3]), np.eye(3)]) ])

        i = self.num_diodes 
        A = np.vstack([ np.hstack([A, np.zeros([6, 3*i])]),
                        np.hstack([np.zeros([3*i, 6]), np.eye(3*i)])   ])

        # (No predicted change in c, a, e values)
        x_p = np.hstack([x_p, self.calibration_values, self.azimuth_angles, self.elevation_angles])

        ## COVARIANCE (P_p)
        temp = np.vstack([ (self.Pchol @ A.T), self.chol(self.W)])
        Pchol_p = np.linalg.qr(temp, "r")  # Only returns R matrix

        return x_p, Pchol_p

    def innovate(self, x_p, Pchol_p):
        """ Runs the "innovate" step of the MEKF by calling the parent method """
        return super().innovate(x_p, Pchol_p) # L, the Kalman Gain

    def update(self, x_p, Pchol_p, L, C, z):
        """ Runs the "update" step of the MEKF by calling the parent method """
        return super().update(x_p, Pchol_p, L, C, z)

    # Wrapper around parent method that adds states for calibration variables
    def mag_measurement(self, x_p, mag_inert_unit, 
                            c = None, a = None, e = None):
        """
            Updates the magnetometer portion of the measurements by calling the parent method 
                and then augmenting the sensitivity matrix H (dy/dx) to include the state values 
                corresponding to the calibration states 

            Arguments:
                - x_p: Predicted state 
                        ( [q, bias, c, a, e])
                - mag_inert_unit:  Unit vector corresponding to magnetic             |  [3,]
                        field vector in inertial frame
                - c (Opt): Calibration values for each photodiode. If not provided,  |  [i,]
                        most recent estimate recorded is used 
                - a (Opt): Azimuth angles for each photodiode. If not provided,      |  [i,]
                        most recent estimate recorded is used 
                - e (Opt): Elevation angles for each photodiode. If not provided,    |  [i,]
                        most recent estimate recorded is used 

            Returns: 
                - y: Expected magnetic field vector measurements given predicted state    |  [3,]
                - H: Sensitivity matrix (Jacobian of y wrt x)                             |  [3 x 6 + 3i]
        """

        c = self.calibration_values if np.any(c == None) else c 
        a = self.azimuth_angles     if np.any(a == None) else a 
        e = self.elevation_angles   if np.any(e == None) else e

        y, H = super().mag_measurement(x_p, mag_inert_unit, c, a, e) # [3,], [3 x 6]

        dc = np.zeros([3, self.num_diodes])
        da = np.zeros([3, self.num_diodes])
        de = np.zeros([3, self.num_diodes])

        H = np.hstack([H, dc, da, de])  # [3 x 6 + 3 i]
        return y, H

    # Wrapper around parent method
    def current_measurement(self, x_p, sun_inert_unit, 
                            c = None, a = None, e = None):
        """
            Updates the current portion of the measurements by calling the parent method 
                and then augmenting the sensitivity matrix H (dy/dx) to include the state values 
                corresponding to the calibration states 

            Arguments:
                - x_p: Predicted state 
                        ( [q, bias, c, a, e])
                - sun_inert_unit:  Unit vector corresponding to sun                  |  [3,]
                        vector in inertial frame
                - c (Opt): Calibration values for each photodiode. If not provided,  |  [i,]
                        most recent estimate recorded is used 
                - a (Opt): Azimuth angles for each photodiode. If not provided,      |  [i,]
                        most recent estimate recorded is used 
                - e (Opt): Elevation angles for each photodiode. If not provided,    |  [i,]
                        most recent estimate recorded is used 

            Returns: 
                - y: Expected current measurements given predicted state     |  [i,]
                - H: Sensitivity matrix (Jacobian of y wrt x)                |  [i x 6 + 3i]
        """

        c = self.calibration_values if np.any(c == None) else c 
        a = self.azimuth_angles     if np.any(a == None) else a 
        e = self.elevation_angles   if np.any(e == None) else e
        
        y, H = super().current_measurement(x_p, sun_inert_unit, c, a, e)

        bQi = dcm_from_q(x_p[:4]).T 
        s_body = bQi @ sun_inert_unit 

        surf_norm = np.array([np.cos(e)*np.cos(a),
                        np.cos(e)*np.sin(a), 
                        np.sin(e)]).T           # [i x 3] 

        dnda = np.array([-np.cos(e)*np.sin(a),
                          np.cos(e)*np.cos(a),
                          np.zeros(a.shape)]).T # [6 x 3] 

        dnde = np.array([-np.sin(e)*np.cos(a),
                         -np.sin(e)*np.sin(a), # negative of the equation in the paper
                          np.cos(e)]).T         # [6 x 3] 

        dc     = surf_norm @ s_body   # [i,]
        da     = c * (dnda @ s_body)  # [i,]
        de     = c * (dnde @ s_body)  # [i,]

        H = np.hstack([H, diag(dc), diag(da), diag(de)]) # [i x 6 + 3i]
        H[y <= 0] = 0   # Diodes cant go negative
        return y, H








