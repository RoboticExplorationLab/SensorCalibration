# [src/MissionSim/CustomStructs/diodes.jl]

"""
    DIODES{S, T} -> calib_values, azi_angles, elev_angles

      Calibration values (~ scale factor) and surface normal (in azimuth and elevation angle) 
    for a set of sun-sensing photodiodes on the CubeSat. Allows for any number of diodes, but 
    the default is six, with one located on each side of the CubeSat. 
"""
struct DIODES{S, T}
    calib_values::SVector{S, T}       #   Calibration value (How much 1 unit of current is scaled by)
    azi_angles::SVector{S, T}         #   Azimuth angle for each photodiode surface normal
    elev_angles::SVector{S, T}        #   Elevation angle for each photodiode surface normal

    function DIODES(cv::SVector{S, T}, aa::SVector{S, T}, ea::SVector{S, T}) where {S, T}

        # Wrap α ∈ (-pi, pi] and ϵ ∈ (-pi/2, pi/2]
        if any(ea .≤ -pi/2) || any(ea .> pi/2)
            # Convert to cartesian, and convert back to stay in range

            x, y, z = sph2cart(aa, ea)
            aa, ea = cart2sph(x, y, z)
            aa = SVector{S, T}(aa)
            ea = SVector{S, T}(ea)
        elseif any(aa .≤ -pi) || any(aa .> pi)
            aa = wrap(aa)
        end
        
        """ Primary Constructor """
        new{S, T}(cv, aa, ea)
    end

    function DIODES(cv, aa, ea)
        """ Constructor for non-static types """
        @warn "Improper declaration of DIODES; converting to static..."

        S = length(cv)
        T = typeof(cv[1])
        DIODES(SVector{S, T}(cv), SVector{S, T}(aa), SVector{S, T}(ea))
    end

    function DIODES(; ideal = false, N = 6)
        """ 
              Randomly generates a set of diodes. If using six diodes, the default orientation is
            set up and perturbed slightly. If some other number of diodes is used, each is assigned
            a completely random surface normal. Also allows for an 'ideal' set of diodes with perfect 
            scale and surface normals exactly as we would want.
        """

        if ideal 
            cv = ones(6)
            ea = [(-pi/4); (pi/4);  0.0;    0.0;    (pi/4); (-pi/4)] 
            aa = [0.0;      pi;     (pi/2); (-pi/2); 0.0;    pi]  

            # 
            # aa = [0.0,  pi, pi/2, -pi/2,  0.0,   0.0]
            # ea = [0.0, 0.0,  0.0,   0.0, pi/2, -pi/2]

            return DIODES(SVector{N, Float64}(cv), SVector{N, Float64}(aa), SVector{N, Float64}(ea))
        end

        # cv = rand(Normal(2.0, 0.3), N)   # On-Orbit... paper has 0.5 or 0.2 
        cv = rand(Normal(1.0, 0.2), N)

        if N == 6
            aa = [0.0;        pi;   (pi/2); (-pi/2);      0.0;      pi]  
            ea = [(-pi/4); (pi/4);     0.0;     0.0;   (pi/4); (-pi/4)] 


            # aa = [0.0,  pi, pi/2, -pi/2,  0.0,   0.0]
            # ea = [0.0, 0.0,  0.0,   0.0, pi/2, -pi/2]

            ea = ea + rand(Normal(0.0, deg2rad(5.0)), N)   # On-Orbit paper has 10 or 2 deg
            aa = aa + rand(Normal(0.0, deg2rad(5.0)), N)
        else
            ea = deg2rad.(rand(-90:90, N))
            aa = deg2rad.(rand(-179:180, N))
        end

        return DIODES(SVector{N, Float64}(cv), SVector{N, Float64}(aa), SVector{N, Float64}(ea))
    end
end

function wrap(a::SVector{N, T}) where {N, T}
    v = atan.(sin.(a), cos.(a))
    return SVector{N, T}(v)
end

function sph2cart(α, ϵ)
    θ = (pi/2) .- ϵ  # Convert to from elevation to inclination or whatever
    x = cos.(α) .* sin.(θ)
    y = sin.(α) .* sin.(θ)
    z = cos.(θ)
    return x, y, z
end

function cart2sph(x, y, z)
    α = atan.(y, x)
    θ = atan.(sqrt.( (x.^2) .+ (y.^2)), z)  # Assumes unit 
    ϵ = pi/2 .- θ
    return α, ϵ
end
