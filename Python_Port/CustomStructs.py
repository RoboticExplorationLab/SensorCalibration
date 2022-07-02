
# Class? OR make it a dict?


class Magnetometer():
    def __init__(self, scale_factors, non_ortho_angles, bias):
        self.scale_factors = scale_factors 
        self.non_ortho_angles = non_ortho_angles 
        self.bias = bias

class Diodes():
    def __init__(self, calib_values, azi_angles, elev_angles):
        self.calibration_values = calib_values 
        self.azimuth_angles   = azi_angles 
        self.elevation_angles  = elev_angles

class Satellite():
    def __init__(self, J, magnetometer, diodes, state, covariance):
        self.J = J
        self.magnetometer = magnetometer 
        self.diodes       = diodes 
        self.state        = state 
        self.covariance   = covariance


    def other():
        pass

class Sensors():
    def __init__(self, magnetometer, diodes, gyro, gps):
        self.magnetometer = magnetometer 
        self.diodes       = diodes 
        self.gyro         = gyro 
        self.gps          = gps

class Noise():
    def __init__(self, diodes, gyro, gps, sb_rot, Bb_rot):
        self.diodes = diodes 
        self.gyro   = gyro
        self.gps    = gps 
        self.sb_rot = sb_rot 
        self.Bb_rot = Bb_rot

class GroundTruth():
    def __init__(self, t, Bi, si, Bb, sb):
        self.t_hist  = t 
        self.Bi_hist = Bi 
        self.si_hist = si 
        self.Bb_hist = Bb
        self.sb_hist = sb 

class Flags():
    def __init__(self, in_sun, magnetometer_calibrated,
                    diodes_calibrated, detumbling, calibrating):
        self.magnetometer_calibrated = magnetometer_calibrated
        self.diodes_calibrated       = diodes_calibrated
        self.in_sun                  = in_sun 
        self.detumbling              = detumbling
        self.calibrating             = calibrating





