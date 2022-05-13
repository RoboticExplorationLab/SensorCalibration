# [tests/Controller_runtests.jl]

""" Runs all tests for files in ../src/MissionSim/Controller/ """

using Test, BenchmarkTools
using LinearAlgebra, StaticArrays
include("../src/MissionSim/Controller/Controller.jl"); using .Controller

include("Controller/detumbler_tests.jl");



