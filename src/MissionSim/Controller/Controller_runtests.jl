# [tests/Controller_runtests.jl]

""" Runs all tests for files in ../src/MissionSim/Controller/ """

using Test, BenchmarkTools
using LinearAlgebra, StaticArrays
include("Controller.jl"); using .Controller

include("detumbler_tests.jl");



