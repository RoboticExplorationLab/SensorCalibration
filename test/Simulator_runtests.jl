# [tests/Simulator_runtests.jl]

# Make a module for the sake of testing to avoid scope issues again?

""" Runs all tests for files in ../src/MissionSim/Simulator """

using Test, BenchmarkTools
using Plots, EarthAlbedo, JLD2
using StaticArrays, LinearAlgebra, SatelliteDynamics, Distributions 

include("../src/MissionSim/Simulator/Simulator.jl"); using .Simulator

include("SimpleOrbit.jl")

include("../src/MissionSim/quaternions.jl");
include("../src/MissionSim/mag_field.jl");


include("Simulator/dynamics_tests.jl");
include("Simulator/measurement_tests.jl");
include("Simulator/simulator_tests.jl");

