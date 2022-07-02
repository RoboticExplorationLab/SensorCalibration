# [tests/Simulator_runtests.jl]

# Make a module for the sake of testing to avoid scope issues again?

""" Runs all tests for files in ../src/Simulator """

using Test, BenchmarkTools
using Plots, EarthAlbedo, JLD2
using StaticArrays, LinearAlgebra, SatelliteDynamics, Distributions 

include("../src/Simulator/Simulator.jl"); using .Simulator

# Ensure the appropriate scripts are included, but only once (important for struct declarations)
(@isdefined STATE)       ? nothing : using .Simulator.CustomStructs
(@isdefined SimpleOrbit) ? nothing : include("SimpleOrbit.jl");
(@isdefined IGRF13)      ? nothing : include("../src/mag_field.jl");
(@isdefined qdot)        ? nothing : include("../src/quaternions.jl");


# Run Tests
include("Simulator/dynamics_tests.jl");
include("Simulator/measurement_tests.jl");
include("Simulator/simulator_tests.jl");

