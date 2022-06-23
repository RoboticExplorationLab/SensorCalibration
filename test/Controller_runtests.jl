# [tests/Controller_runtests.jl]

""" Runs all tests for files in ../src/MissionSim/Controller/ """

using Test, BenchmarkTools
using LinearAlgebra, StaticArrays, SatelliteDynamics

include("../src/MissionSim/Controller/Controller.jl"); using .Controller

# Ensure the appropriate scripts are included, but only once (important for struct declarations)
(@isdefined STATE)       ? nothing : using .Controller.CustomStructs
(@isdefined SimpleOrbit) ? nothing : include("SimpleOrbit.jl");
(@isdefined IGRF13)      ? nothing : include("../src/MissionSim/mag_field.jl");
(@isdefined qdot)        ? nothing : include("../src/MissionSim/quaternions.jl");

# Run Tests
include("Controller/detumbler_tests.jl");



