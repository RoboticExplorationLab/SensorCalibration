import numpy as np 
# import ulab.numpy as np 

class Detumbler():
    # MEMBER CONSTANTS 
    k = None  # Optimal gain constant for given inertia matrix
    
    # MEMBER VARIABLES 
    #   (None)

    def __init__(self, k = 7e-6):
        self.k = k

    def determine_optimal_gain(self, T_orbit, inclination, J):
        Jmin    = np.min( [J[0,0], J[1,1], J[2,2]] ) # Min of diagonal terms 
        Omega   = 2 * np.pi / T_orbit

        self.k = 2 * Omega * (1 + np.sin(inclination)) * Jmin

    def generate_command(self, omega, B_body):
        # Using B-Cross 
        B_body_unit = B_body / np.linalg.norm(B_body)

        return (-self.k * (np.eye(3) - (B_body_unit[:,None] @ B_body_unit[None, :])))  @ omega
