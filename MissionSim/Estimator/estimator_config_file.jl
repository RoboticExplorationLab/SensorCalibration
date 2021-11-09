
# PHOTODIODES ------------------------------------------------------------------

#NOISE VALUES! 
const σ_q = (10*pi/180) 
const σ_βgyro = (10*pi/180)
const σ_βmag  = (10 * pi / 180)
const σ_c = 0.15 # 0.2 
const σ_α = deg2rad(2.0) # 1.0 #2.0   # 2.0 is the σ used when generating these
const σ_ϵ = deg2rad(2.0) # 0.3 #1.0 


const _E_am0 = 1366.9 # Irradiance of sunlight (TSI - visible & infrared), W/m^2 

##### VERIFY THIS SECTION ########
const estimator_params = (angle_random_walk      = 0.06,   # in deg/sqrt(hour)   
                          gyro_bias_instability  = 0.8,    # Bias instability in deg/hour
                          velocity_random_walk   = 0.014,  # in m/sec/sqrt(hour)
                          accel_bias_instability = 6)      # in microG

const Q_gyro = ((estimator_params[:gyro_bias_instability] * (pi/180)    )^2)/(3600^3)  
const σ_orient = sqrt(Q_gyro);

const Q_bias = ((estimator_params[:angle_random_walk]*(pi/180))^2)/(3600)   # This is super small
const σ_bias_gyro = sqrt(Q_bias)
const σ_bias_mag  = sqrt(Q_bias)  # TODO arbitrary

const Q_diode = 1e-5   # Diode Noise 
const σ_cal = Q_diode
const σ_azi = Q_diode
const σ_ele = Q_diode;

const σ_sunVec = deg2rad(3.0)
const σ_magVec = deg2rad(3.0)
const σ_curr = 0.05 #0.008; #3, 3, 0.005
# NOTE (Not sure why σ_cal != σ_c... )
##################################



