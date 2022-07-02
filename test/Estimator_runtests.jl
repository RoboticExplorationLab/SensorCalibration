# [tests/Estimator_runtests.jl]

""" Runs all tests for files in ../src/Estimator/ """


using Test, BenchmarkTools 
using StaticArrays, LinearAlgebra
using ForwardDiff, Plots, EarthAlbedo, SatelliteDynamics, JLD2
using Distributions

include("../src/Estimator/Estimator.jl");  using .Estimator 

# Ensure the appropriate scripts are included, but only once (important for struct declarations)
(@isdefined STATE)       ? nothing : using .Estimator.CustomStructs
(@isdefined SimpleOrbit) ? nothing : include("SimpleOrbit.jl");
(@isdefined IGRF13)      ? nothing : include("../src/mag_field.jl");
(@isdefined qdot)        ? nothing : include("../src/quaternions.jl");

function chol(M)
    # return cholesky(Symmetric(M)).U 
    return cholesky(Hermitian(M)).U
end
include("Estimator/magnetometer_calibration_tests.jl");
include("Estimator/mekf_tests.jl");


