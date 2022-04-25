# Tests for Triad
using Random, LinearAlgebra, Plots

include("../../CustomStructs.jl");  using .CustomStructs
include("../../rotationFunctions.jl")
include("../diode_calibration.jl")


errs = zeros(100)
for i = 1:100
    # Randomly generate two Inertial vectors 
    r₁ᴵ = rand(1:1000) * randn(3)
    r₂ᴵ = rand(1:1000) * randn(3)

    # Randomly generate a rotation matrix/quaternion 
    q = randn(4); q = q / norm(q)
    R = dcm_from_q(q)

    # Apply, pass in
    r₁ᴮ = R * r₁ᴵ
    r₂ᴮ = R * r₂ᴵ

    # Compare results

    q̂, R̂ = triad(r₁ᴵ, r₂ᴵ, r₁ᴮ, r₂ᴮ)
    R̂ = R̂'

    errs[i] = sum(abs.(R - R̂))
    # println(sum(abs.(R - R̂)))
end

plot(errs)