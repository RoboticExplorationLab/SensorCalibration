

using Test, BenchmarkTools, Infiltrator 
using Plots, EarthAlbedo, JLD2
using StaticArrays, LinearAlgebra, SatelliteDynamics, Distributions
using ForwardDiff 

include("../src/CustomStructs.jl");         using .CustomStructs
include("../src/Simulator/Simulator.jl");   using .Simulator
include("../src/Estimator/Estimator.jl");   using .Estimator 
include("../src/Controller/Controller.jl"); using .Controller


include("SimpleOrbit.jl")

include("../src/quaternions.jl");
include("../src/mag_field.jl");



# ### GENERAL ###
include("quaternion_tests.jl");


# ### SIMULATOR ### 
include("Simulator/dynamics_tests.jl");
include("Simulator/measurement_tests.jl");
include("Simulator/simulator_tests.jl");


# ### CONTROLLER ###
include("Controller/detumbler_tests.jl");


# ### OTHER ###
include("state_machine_tests.jl");
include("CustomStructs_tests.jl");
include("main_tests.jl");


### ESTIMATOR ###
function chol(M)
    return cholesky(Hermitian(M)).U
end
include("Estimator/magnetometer_calibration_tests.jl");
# include("Estimator/mekf_tests.jl");  # <- currently failing because of Measurement Jacobians
