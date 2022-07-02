# In theory, this will not ever be run on the satellite ya...?
import numpy as np
from satelliteDynamics import * 
from rotationFunctions import * 


# TODO 
print("Verify how Im generating random numbers - should be Normal with mu and sigma")
# (IDK if np.random.random is included in ulab...)

class Simulator():

    ### MEMBER CONSTANTS
    earth_radius    = 6378136.3  # m
    _E_am0 = 1366.9

        # SYSTEM 
    dt  = None
    sat = None
    num_diodes = None
    alb = None 
    cell_centers_ecef = None
    mag_calibration_matrix = None 
    magnetometer_bias = None

        # NOISE  (TODO are these techincally sigma^2?)
    mag_noise_sigma  = None  # in radians, for rotation
    gps_noise_factor = None  # linear multiplier 
    gyro_noise_sigma  = None
    sun_noise_sigma  = None  # in radians, for rotation
    current_noise_sigma = None

    ### MEMBER VARIABLES

    # Since we use a true satellite that doesnt change, can i just member-variablize it?
    def __init__(self, sat, dt, alb):
        self.dt = dt
        self.mag_calibration_matrix = self.get_mag_calib_matrix(sat)
        self.magnetometer_bias      = sat.magnetometer.bias
        self.alb = alb 
        self.sat = sat
        self.num_diodes = sat.diodes.calib_values.shape[0]
        self.cell_centers_ecef = alb.get_albedo_cell_centers()


    ################
    #   DYNAMICS   #
    ################

    def dynamics(self, J, x, u, t):
        """ 
            Propagates the state dynamics 
                ([r⃗, v⃗, (q⃗, q₀), ω, β]) 
        """
        if np.linalg.norm(x[:3] < self.earth_radius):
            raise Exception("Impact at time {}".format(t))

        rdot = x[3:6]
        vdot = accel_perturbations(t, x[:6])

        w    = x[11:13]   # Angular Velocity 
        qdot = 0.5 * qmult(x[6:10], np.array([w, 0]))  # Scalar-last quaternion

        wdot = np.linalg.inv(J) * (u - np.cross(w, (J * w)))

        # TODO Not actually sure if this is the right value for bias!
        biasDot = self.bias_sigma * np.random.random(3)

        return np.hstack([ rdot, vdot, qdot, wdot, biasDot])

    def accel_perturbations():
        pass 

    def rk4(self, J, x, u, t, h):
        """ Modified rk4 function """
        k1 = h * self.dynamics(J, x, u, t)
        k2 = h * self.dynamics(J, x + k1/2, u, t + h/2)
        k3 = h * self.dynamics(J, x + k2/2, u, t + h/2)
        k4 = h * self.dynamics(J, x + k3,   u, t + h)

        return (x + (1/6) * (k1 + 2*k2 + 2*k3 + k4))


    ################
    # MEASUREMENTS #
    ################
    def generate_measurements(self, x, t, CONSTANTS, dt): # and alb
        #TODO Update comments
        """ UPDATE!
            Generates sensor measurements, including noise.

            Arguments:
            - sim: Used to determine which simulator to use                                 | SIM 
            - sat: TRUE satellite data, which is used to generate sensor measurements       | SATELLITE
                        J (inertia matrix, [3 x 3])          //  magnetometer (calibration values, [3,]),
                        diodes (photodiode parameters [3,])  //  state (attitude and bias, [7,])
            - alb: Albedo struct containing REFL data and cell centers                      | ALBEDO
                        REFL    //   cell_centers_ecef
            - x: Environmental state (pos, vel, att, ang vel, bias)                         | [16,]
            - t: Current time, as an epoch                                                  | Scalar{Epoch}
            - CONSTANTS: Various constants needed (Earth radius, μ , E_am₀)                 | [3,]       
            - dt: Time step                                                                 | Scalar
            
            Returns:
            - truth: TRUTH struct containing true values for current time and vectors in body/inertial frame    | TRUTH
            - sensors: SENSORS struct containing noisy values for each satellite sensor being simulated         | SENSOR
            - ecl: Scale factor used to track how much the satellite is eclipsed by the Earth (∈ [0, 1])        | Scalar
            - noise: NOISE struct containing the amount of noise added to the sensors                           | NOISE
        """  

        bRi = dcm_from_q(x[6:10]).T 
        Bi, Bb, Bb_meas          = self.magnetometer_measurement(x[:3], t, bRi)
        w, w_meas, gyro_noise    = self.gyro_measurement(x)
        gps, gps_meas, gps_noise = self.gps_measurement(x)
        
        si, sb_unit, ecl         = self.update_sun_vectors(x[:3], t, bRi)
        I, I_meas, I_noise    = self.diode_measurement( x[:3], si, sb_unit, ecl)

        sensors = SENSORS(Bb_meas, I_meas, w_meas, gps_meas)
        truth   = GROUND_TRUTH(t, Bi, si, sb_unit, Bb)
        noise   = NOISE(I_noise, gyro_noise, gps_noise)

        return truth, sensors, noise, ecl # NOTE ADJUSTED ORDER!

    # SHOULD IGRF13 be a member method? 
    def magnetometer_measurement(self, pos, time, bRi):
        # Bi = IGRF13(pos, time)

        igrf = igrfclass()  
        mag = np.sqrt( (pos[0]**2) + (pos[1]**2) + (pos[2]**2) )
        lat  = -48.70272645037185 # np.degrees(np.arcsin(pos[2] / mag )) 
        lon  = 5.206967731682121  # np.degrees(np.arctan2(pos[1], pos[0]))
        r_norm_km = mag / 1000 
        a = caldate(time)
        date = 2021.384675934482 # a[0] + (30 * a[1] + a[2]) / 365.2425 # Roughly the year 
        Bi = igrf.ned_igrf(date, lat, lon, r_norm_km)

        print("Lon: {:3f}\t Lat: {:3f}".format(lon, lat))
        print("Rkm: {:3}\t  Date: {}".format(r_norm_km, date))
        print("Bi: ", Bi)
        Bb = bRi @ Bi

        mag_noise = self.noise_matrix(self.mag_noise_sigma)

        Bb_meas = self.mag_calibration_matrix @ mag_noise @ Bb + self.magnetometer_bias
        
        return Bi, Bb, Bb_meas


    # TODO Should I be adding in noise linearly or with a rotation...?
    #   (Not a rotation, so should be ok...? But then i could remove norm)
    def gyro_measurement(self, x):
        w    = x[10:13]
        bias = x[13:16]
        gps_noise = self.gps_noise_factor * np.linalg.norm(w) * np.random.random(3)

        w_meas = w + bias + gps_noise 

        return w, w_meas, gps_noise 

    def gps_measurement(self, x):
        pos = x[:3]
        pos_noise = self.gps_noise_factor * np.random.random(3)
        
        return (pos + pos_noise), pos_noise

    def diode_measurement(self, pos, si, sb_unit, ecl):

        c = self.sat.diodes.calib_values #calibration_values 
        a = self.sat.diodes.azi_angles   #azimuth_angles 
        e = self.sat.diodes.elev_angles  #elevation_angles 
        
        i = self.num_diodes

        surf_norm = np.array([np.cos(e)*np.cos(a),
                                np.cos(e)*np.sin(a), 
                                np.sin(e)]).T     # [i x 3] 

        albedo_matrix, _ = self.alb.albedo(pos, si, self.alb.data) 

        diode_albedos = self.alb.get_diode_albedo(albedo_matrix, self.cell_centers_ecef, surf_norm, pos)
        diode_albedos = c * diode_albedos / self._E_am0

        I = (c * (surf_norm @ sb_unit) + diode_albedos) * ecl
        I_noise = self.current_noise_sigma * np.random.random(i)

        I_meas = (I + I_noise) * ecl
        I[I < 0.0] = 0.0 
        I_meas[I_meas < 0.0] = 0.0

        return I, I_meas, I_noise

    print("Using wrong sign in python update sun!")
    def update_sun_vectors(self, pos, time, bRi):
        si_earth = sun_position(time)
        si_sat   = si_earth - pos 

        ecl = eclipse_conical(-pos, si_earth)
        si_sat *= ecl

        sun_noise = self.noise_matrix(self.sun_noise_sigma)

        # Avoid divide-by-zero error during eclipses
        sb_unit = (sun_noise @ bRi @ (si_sat / np.linalg.norm(si_sat))) if ecl > 0.01 else np.zeros(3)

        return si_sat, sb_unit, ecl 

    def noise_matrix(self, sigma):
        if (sigma != 0.0):
            noise_vec = sigma * np.random.random(3)
            
            noise_norm = np.linalg.norm(noise_vec)
            skew = hat(noise_vec) / noise_norm

            # Rodrigues for matrix exponential (?)
            R = np.eye(3) + skew * np.sin(noise_norm * self.dt) + (skew @ skew) * (1 - np.cos(noise_norm * self.dt))
        else:
            R = np.eye(3)

        return R

    def get_mag_calib_matrix(self, sat):
        """ Generates the calibration matrix that alters the measured magnetic field vector in body frame """
        
        a = sat.magnetometer.scale_factors[0]
        b = sat.magnetometer.scale_factors[1]
        c = sat.magnetometer.scale_factors[2]

        rho = sat.magnetometer.non_ortho_angles[0] 
        lamb = sat.magnetometer.non_ortho_angles[1] 
        phi = sat.magnetometer.non_ortho_angles[2] 

        # a, b, c = sat.magnetometer.scale_factors 
        # rho, lamb, phi = sat.magnetometer.non_ortho_angles 

        T = np.array([  
            [a,              0.0,                        0.0],
            [b*np.sin(rho),  b*np.cos(rho),              0.0],
            [c*np.sin(lamb), c*np.sin(phi)*np.cos(lamb), c*np.cos(phi)*np.cos(lamb)]
        ])

        return T

