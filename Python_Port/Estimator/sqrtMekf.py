""" 
    TEMPLATE Square-root Multiplicative Kalman Filter (MEKF) Class
"""

import numpy as np
from satelliteDynamics import * 
from rotationFunctions import * 
from earthAlbedo import EarthAlbedo 


# TODO List:
# - Comment, clean

# NOTE Currently a template class / interface
class SqrtMekf():

    # Member Constants 
    # None here

    # Member Variables
    # None here

    def __init__(self):
        pass 

    def predict(self, q, bias, Pchol, gyro, W, dt):
        """
            Predicts the next orientation of the agent by applying a small rotation.
                Also returns the Jacobian of the state wrt the state (dx/dx)
                as well as the predicted state covariance matrix P_p
                (NOTE that because we are using quaternions, x will have more parameter than A and P)

            Arguments:
                - q:     Current orientation of the agent as a scalar-last quaternion       |  [4,]
                - bias:  Estimated gyroscope bias                                           |  [3,]
                - gyro:  Measured angular velocity                                          |  [3,]
                - dt:    System time step                                                   |  Scalar

            Returns:
                - x_p:   Predicted state - [quaternion, bias]                               |  [7,]
                - A:     Jacobian of the state wrt the state                                |  [6 x 6]
                - P_p:   Predicted state covariance matrix                                  |  [6 x 6]
        """
        gamma = gyro - bias  # Adjusted angular velocity 
        gamma_mag  = np.linalg.norm(gamma)
        gamma_skew = -hat(gamma) / gamma_mag 

        ## STATE (x_p)
        theta = gamma_mag * dt   # Angle of rotation of agent
        r = gamma / gamma_mag    # Axis of rotation (unit)
        p = np.hstack([r*np.sin(theta/2), np.cos(theta/2)])  # Quaternion representation

        q_p = qmult(q, p)       # Predicted quaternion

        x_p = np.hstack([q_p, bias])    # Predicted state

        ## JACOBIAN (A)
        # Rodrigues
        R = np.eye(3) + (gamma_skew * np.sin(theta)) + (gamma_skew @ gamma_skew) * (1 - np.cos(theta))

        # Jacobian (dx/dx)
        A = np.vstack([ np.hstack([R, -dt*np.eye(3)]),
                        np.hstack([np.zeros([3,3]), np.eye(3)]) ])

        ## COVARIANCE (P_p)
        temp = np.vstack([ (Pchol @ A.T), self.chol(W)])
        Pchol_p = np.linalg.qr(temp, "r") 

        return x_p, A, Pchol_p

    def innovate(self, C, z, V, Pchol_p):
        """
            Determines the optimal Kalman gain for the system by comparing the 
                expected measurements with the actual measurements, factoring in 
                noise and the sensitivity matrix and using QR decomposition
                (NOTE that this does not update measurement values, which must be done elsewhere)

            Arguments:
                - C:       Sensitivity Matrix (Jacobian of outputs y wrt inputs x)
                - z:       Innovation vector (difference between expected and received measurements)
                - V:       Measurement noise matrix 
                - Pchol_p: Predicted state covariance as an upper Cholesky factor
            
            Returns:
                - L:       Kalman Gain
        """

        temp = np.vstack([ (Pchol_p @ C.T), self.chol(V) ])
        Pyy_chol = np.linalg.qr( temp, "r") # Only need R 

        Pyy_chol_inv = np.linalg.inv(Pyy_chol)

        L = ((Pchol_p.T @ Pchol_p @ C.T) @ Pyy_chol_inv) @ (Pyy_chol_inv.T)

        return L

    def update(self, x_p, Pchol_p, L, C, z, V):
        """
            Updates the predicted state and covariance values using the Kalman gain calculated previously.
                Converts 3-parameter angle/axis representation into a quaternion and uses a 
                Hamiltonian product to update. Bias is updated linearly
                (NOTE that dimensions are for the vanilla MEKF, and may not be accurate for a child class implementation)

            Arguments:
                - x_p:      Predicted state                        |  [7,]
                                ( [quaternion, bias] )
                - Pchol_p:  Predicted covariance matrix            |  [6 x 6]
                - L:        Kalman gain
                - C:        Sensitivity matrix (dy/dx)
                - z:        Innovation vector (y_expected - y_measured)
                - V:        Measurement noise

            Returns:
                - x_next:      Predicted next state                   |  [7,]
                                    ( [quaternion, bias] )
                - Pchol_next:  Next state covariance                  |  [6 x 6]
                                    (Upper Cholesky factor)
        """
        dx = L @ z  
        
        # Convert to angle/axis and then to quaternion update
        dphi = dx[:3]
        theta = np.linalg.norm(dphi)    # Angle of rotation
        r = dphi/theta                  # Axis of rotation (unit)

        dq = np.hstack([ r * np.sin(theta/2), np.cos(theta/2) ])

        x_next = np.hstack([ qmult(x_p[:4], dq), dx[3:] + x_p[4:]  ])

        I = np.eye(Pchol_p.shape[0])
        temp = np.vstack([ Pchol_p @ (I - (L @ C)).T, self.chol(V) @ L.T])
        Pchol_next = np.linalg.qr(temp, "r")
        
        return x_next, Pchol_next

    ####################
    # HELPER FUNCTIONS #
    ####################
    def chol(self,M):
        """ Returns an Upper Cholesky factor of the given matrix """
        return np.linalg.cholesky(M).T
