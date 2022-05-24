# [tests/Estimator_runtests.jl]

""" Runs all tests for files in ../src/MissionSim/Estimator/ """


using Test, BenchmarkTools 
using StaticArrays, LinearAlgebra
using ForwardDiff, Plots, EarthAlbedo, SatelliteDynamics, JLD2
using Distributions

include("../src/MissionSim/Estimator/Estimator.jl");  using .Estimator 

using .Estimator.CustomStructs  # so can just call them directly with Estimator.

include("../src/MissionSim/quaternions.jl")
include("../src/MissionSim/mag_field.jl")

# include("Estimator/magnetometer_calibration_tests.jl");

function chol(M)
    # return cholesky(Symmetric(M)).U 
    return cholesky(Hermitian(M)).U
end

include("Estimator/mekf_tests.jl");


