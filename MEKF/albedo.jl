# Calculation of albedo for a given satellite and sun location
#   and specified reflectivity data
# 
# (Ported over from Matlab's Earth Albedo Toolbox)

using MAT
using CoordinateTransformations # For Cartesian/Spherical conversions
using StaticArrays              # For speed
using LinearAlgebra

# CONSTANTS 
_Re = 6371.01e3;   # Average radius of the Earth, m
_Am0 = 1366.9;     # Irradiance of the sun, W/m^2 
_d2r = pi / 180.0; # Standard degrees to radians conversion

struct refl_struct
    """ REFL Structure for holding reflectivity data"""
    data
    type
    start_time
    stop_time
end

function albedo(sat, sun, refl)
    """
        Determines the Earth's albedo at a satellite for a given set of reflectivity data.
            - Divides Earth into a set of cells
            - Determines which cells are visible by both the satellite and the Sun
            - Determines how much sunlight is reflected by the Earth towards the satellite 
            (Primary function)

        Inputs:
            - sat: vector from center of Earth to satellite                                   |  [3,]
            - sun: vector from center of Earth to Sun                                         |  [3,]
            - refl: Refl struct containing the averaged reflectivity data for each cell       | Custom struct

        Outputs:
            - cell_albedos: Matrix with each element corresponding to the albedo coming from 
                                the corresponding cell on the surface of the earth            | [num_lat x num_lon]
    """
    
    # Adjust input dimenions to support both row and column vectors
    sat = sat[:];
    sun = sun[:];

    # 
    num_lat, num_lon = size(refl.data); 

    # Verify no satellite collision has occured 
    if norm(sat) < _Re
        error("albedo.m: The satellite has crashed into Earth!");
    end

    # Convert from Cartesian -> Spherical (r θ ϕ)
    sat_sph_t = SphericalFromCartesian()(sat);
    sun_sph_t = SphericalFromCartesian()(sun);

    # Adjust order to match template code to (θ ϵ r), with ϵ = (π/2) - ϕ
    sat_sph = [sat_sph_t.θ; (pi/2) - sat_sph_t.ϕ; sat_sph_t.r ];
    sun_sph = [sun_sph_t.θ; (pi/2) - sun_sph_t.ϕ; sun_sph_t.r ];

    # REFL Indices 
    sat_i, sat_j = rad2idx(sat_sph[1], sat_sph[2], num_lat, num_lon);
    sun_i, sun_j = rad2idx(sun_sph[1], sun_sph[2], num_lat, num_lon);

    # SKIP GENERATING PLOTS FOR NOW

    satellite_fov_cells = earthfov(sat_sph, num_lat, num_lon) # Cells in the satellite field-of-view
    sunlit_cells     =    earthfov(sun_sph, num_lat, num_lon) # Cells lit by the sun
    union = (satellite_fov_cells .!= 0) .& (sunlit_cells .!= 0) # Cells lit by the sun AND visible by the satellite

    cell_albedos = zeros(num_lat, num_lon)
    for i = 1:num_lat 
        for j = 1:num_lon 
            if union[i,j] # If a cell is illuminated by the sun AND visible by the satellite...
                ϕ_in = gridangle(i, j, sun_i, sun_j, num_lat, num_lon); # angle of incident solar Irradiance

                if ϕ_in > pi/2 # Account for numerical inaccuracies
                    ϕ_in = pi/2
                end

                E_inc = _Am0 * cellarea(i, j, num_lat, num_lon) * cos(ϕ_in) # Incident power  (Pc?)

                grid_theta, grid_phi = idx2rad(i, j, num_lat, num_lon); # Distance to satellite from grid center

                grid_spherical = Spherical(_Re, grid_theta, (pi/2) - grid_phi)
                grid = CartesianFromSpherical()(grid_spherical);    

                sat_dist = norm(sat - grid);  # Unit vector pointing from grid center to satellite

                ϕ_out = acos( ((sat - grid)/sat_dist)' * grid / norm(grid) ); # Angle to sat from grid

                P_out = E_inc * refl.data[i, j] * cos(ϕ_out) / (pi * sat_dist^2);

                cell_albedos[i, j] = P_out;
            
            end
        end
    end

    return cell_albedos;

end

function rad2idx(θ, ϕ, sy, sx)
    """ Transforms location (in radians) to TOMS REFL matrix indices """
    dx = 2 * pi / sx; # 360* / number of longitude cells 
    dy = pi / sy;     # 180* / number of latitude cells 

    i = round( (pi - dy/2 - ϕ)/dy ) + 1;
    j = round( (θ + pi - dx/2 )/dx ) + 1;

    # Adjust so that 180/-180 is included in the interval 
    if i == 0
        i = 1;
    end
    if j == 0;
        j = 1;
    end

    return i, j
end

function idx2rad(i, j, sy, sx)
    """ Transforms TOMS REFL matrix indices to radians """

    dx = 2 * pi / sx;
    dy = pi / sy;

    ϕ = pi - dy/2 - (i-1)*dy;
    θ = (j-1) * dx - pi + dx/2;

    return θ, ϕ
end

function earthfov(sat_sph, sy, sx)
    """ 
        Determines the field of view on earth using spherical coordinates.

        Input:
            - sat_sph: vector to the satellite in ECEF frame using spherical coordinates (from Earth or Sun) | [3,]
            - sy, sx: Number of latitude and longitude cells, respectively                                   | Scalars
    """

    default_outvalue = 1; # Circle value (?)

    if sat_sph[3] < _Re  # LEO shortcut (?)
        sat_sph[3] = sat_sph[3] + _Re
    end

    dx = 2 * pi / sx; 
    dy = pi / sy;

    result = zeros(sy, sx)

    # Small circle center (?)
    θ_0 = sat_sph[1];
    ϕ_0 = sat_sph[2];

    ρ = acos(_Re / sat_sph[3]) # FOV on Earth 
    
    for i = 1:sy
        for j = 1:sx
            θ, ϕ = idx2rad(i, j, sy, sx)
            rd = acos( sin(ϕ_0)*sin(ϕ)*cos(θ_0-θ) + cos(ϕ_0)*cos(ϕ) ); # Radial Distance
            
            if rd <= ρ 
                result[i,j] = default_outvalue;
            end

        end 
    end

    return result
end

function gridangle(i1, j1, i2, j2, sy, sx)
    """ Calculate the angle between two grid index pairs """

    θ1, ϕ1 = idx2rad(i1, j1, sy, sx)
    θ2, ϕ2 = idx2rad(i2, j2, sy, sx);

    ρ = acos( sin(ϕ1)*sin(ϕ2)*cos(θ1 - θ2) + cos(ϕ1)*cos(ϕ2) );

    return ρ
end

function cellarea(i, j, sy, sx)
    """ Calculate area of TOMS cell. Input is in refl matrix indices """

    θ, ϕ = idx2rad(i, j, sy, sx) # Convert to angles (radians)
    dϕ = (180.0 / sy) * _d2r;
    dθ = (360.0 / sx) * _d2r;

    # Get diagonal points
    ϕ_max = ϕ + dϕ/2;
    ϕ_min = ϕ - dϕ/2;

    A = (_Re^2) * dθ * (cos(ϕ_min) - cos(ϕ_max));

    return A
end



""" TESTING """
# sun = [  2.7965379648052444e10; -1.3250066780843158e11; -5.7446070038960976e10];
# sat = [ 68264.69562539333; 6.92779997710015e6;  0.0];

# refl_dict = matread("refl.mat")
# refl = refl_struct(refl_dict["data"], refl_dict["type"], refl_dict["start_time"], refl_dict["stop_time"])
# a = albedo(-sat, sun, refl);

# # COMPARED WITH 
# using MATLAB
# mat"""
#     addpath TEMP;
#     addpath TEMP/private_like
#     $a = albedo(-$sat, $sun, $refl)
# """