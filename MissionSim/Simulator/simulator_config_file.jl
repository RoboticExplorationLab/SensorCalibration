################
# NOISE VALUES #
################

const σ_mag = deg2rad(2.0)  # radians  

const μ_gyro_scale = 0.05   # Scale factor 
const σ_gyro_scale = 0.005  #   (Both scale with actual value)

const μ_gps = 5e4           # m
const σ_gps = 0.1 * μ_gps   # m 

const σ_sun = deg2rad(2.0)  # radians 

const σ_current_scale = 0.05 # Scale factor 

const σ_gyro_bias = 0.25 * deg2rad(0.22) # TODO this is random
const σ_mag_bias  = 0.5 #1.0  # TODO this is also random

const _Re = 6378136.3                 # Radius of the earth (m))
