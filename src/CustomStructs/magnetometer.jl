# [src/CustomStructs/magnetometer.jl]

"""
    MAGNETOMETER{T} -> scale_factors, non_ortho_angles, bias

      Magnetometer sensor with that includes scale factors along each axis, non-orthogonality angles
    (angular distance from perfectly orthogonal axes), and a constant bias. 

      Includes a primary constructor that allows for static and non-static inputs, as well as a  
    generator that can make an ideal magnetometer or a random noisy one. 
"""
struct MAGNETOMETER{T}
    scale_factors::SVector{3, T}         # Linear scale factors for soft iron materials (a, b, c)    |   [3,]
    non_ortho_angles::SVector{3, T}      # Offset from purely orthogonal vectors (ρ, λ, ϕ)           |   [3,]
    bias::SVector{3, T}                  # Constant bias on each axis (DC offset) (x₀, y₀, z₀)       |   [3,] 

    function MAGNETOMETER(s::SVector{3, T}, noa::SVector{3, T}, b::SVector{3, T}) where {T}
        """ Primary Constructor """
        new{T}(s, noa, b)
    end

    function MAGNETOMETER(s, noa, b) 
        """ Converts to static """
        T = typeof(s[1]) 
        MAGNETOMETER(SVector{3, T}(s), SVector{3, T}(noa), SVector{3, T}(b))
    end

    function MAGNETOMETER(; ideal = false) 
        """ Randomly generate """

        if ideal # Generate a perfect magnetometer
            sf  = SVector{3, Float64}(ones(3))   # No scale 
            noa = SVector{3, Float64}(zeros(3))  # No offset 
            b   = SVector{3, Float64}(zeros(3))  # No Bias
            return MAGNETOMETER(sf, noa, b);
        else # Make a noisy one
            sf  = rand(Normal(1.0, 0.2), 3) 
            noa = rand(Normal(0.0, deg2rad(3.0)), 3) 
            b   = rand(Normal(0.0, 1.0), 3)
            return MAGNETOMETER(SVector{3, Float64}(sf), SVector{3, Float64}(noa), SVector{3, Float64}(b))
        end
    end
end