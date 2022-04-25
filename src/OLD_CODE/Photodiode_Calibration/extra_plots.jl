###############################################################################################
# Extra data visualization plots - to be run after fullDataGeneration.jl
###############################################################################################


#---------------------------------------------------------------------------------------------#
# Plot cell centers
#---------------------------------------------------------------------------------------------#
"""
for i = 1:1
    pyplot()

    δa = 5; δb = 5;

    tt_x = cell_centers_ecef[1:δa:end,1:δb:end,1] * dscale;
    tt_y = cell_centers_ecef[1:δa:end,1:δb:end,2] * dscale;
    tt_z = cell_centers_ecef[1:δa:end,1:δb:end,3] * dscale;

    plot_sun = 0.0002 * sun;

    t = scatter3d(reshape(tt_x, length(tt_x)), reshape(tt_y, length(tt_y)), reshape(tt_z, length(tt_z)), color = :blue, label = "Grid centers")

    t = scatter3d!([pos[1,1]], [pos[2,1]], [pos[3,1]], color = :red, label = "Satellite")

    t = scatter3d!([plot_sun[1,1]], [plot_sun[2,1]], [plot_sun[3,1]], color = :yellow, markersize = 40, label = false)

    sat_pos = 1.0 * pos
    i = 1
    C_sat = [sat_pos[1,i], sat_pos[2,i], sat_pos[3,i]]
    vx, vy, vz = 30 * rN1[:,i] * dscale;
    F = [C_sat[1] + vx, C_sat[2] + vy, C_sat[3] + vz];  # WHY not NEGATIVE NOW???
    plot3d!([C_sat[1], F[1]], [C_sat[2], F[2]], [C_sat[3], F[3]], label = false)

    lims = 2e7
    t = scatter3d!(xlim = [-lims, lims], ylim = [-lims, lims], zlim = [-lims, lims])

    display(t)
end

#---------------------------------------------------------------------------------------------# """



#---------------------------------------------------------------------------------------------#
# (Simpler) Plot Earth with sun/sat-lit cells highlighted
#---------------------------------------------------------------------------------------------#
"""
display(for i = 1:1 
    pyplot()
    δa = 7;
    δb = 7;

    lims = 1.2e7

    t_x = cell_centers_ecef[1:δa:end,1:δb:end,1] * dscale;
    t_y = cell_centers_ecef[1:δa:end,1:δb:end,2] * dscale;
    t_z = cell_centers_ecef[1:δa:end,1:δb:end,3] * dscale;

    # scatter3d(reshape(t_x, length(t_x)), reshape(t_y, length(t_y)), reshape(t_z, length(t_z)), color = :blue, xlim = [-lims, lims], ylim = [-lims, lims], zlim = [-lims, lims],  label = false)
    scatter3d([0.0], [0.0], [0.0], color = :black, xlim = [-lims, lims], ylim = [-lims, lims], zlim = [-lims, lims],  label = "Grid centers")

    sat_pos = 1.5 .* pos
    T = size(pos,2)
    @gif for i = 1:T
        scatter3d(reshape(t_x, length(t_x)), reshape(t_y, length(t_y)), reshape(t_z, length(t_z)), markersize = 1, color = :gray, xlim = [-lims, lims], ylim = [-lims, lims], zlim = [-lims, lims],  label = false)
        tt_x_lit = cell_centers_ecef[:,:,1] .* reshape(unions[:,i], (180, 288)) .* dscale;
        tt_y_lit = cell_centers_ecef[:,:,2] .* reshape(unions[:,i], (180, 288)) .* dscale;
        tt_z_lit = cell_centers_ecef[:,:,3] .* reshape(unions[:,i], (180, 288)) .* dscale;

        tt_x = tt_x_lit[1:δa:end, 1:δb:end]
        tt_y = tt_y_lit[1:δa:end, 1:δb:end]
        tt_z = tt_z_lit[1:δa:end, 1:δb:end]

        scatter3d!(reshape(tt_x, length(tt_x)), reshape(tt_y, length(tt_y)), reshape(tt_z, length(tt_z)), color = :yellow, label = false)
        scatter3d!([sat_pos[1,i]], [sat_pos[2,i]], [sat_pos[3,i]], markersize = 6, color = :red, xlim = [-lims, lims], ylim = [-lims, lims], zlim = [-lims, lims], label = false)
        C_sat = [sat_pos[1,i], sat_pos[2,i], sat_pos[3,i]]
        
        vx, vy, vz = 10 * rN1[:,i] * dscale;
        F = [C_sat[1] + vx, C_sat[2] + vy, C_sat[3] + vz];  # WHY not NEGATIVE NOW???
        plot3d!([C_sat[1], F[1]], [C_sat[2], F[2]], [C_sat[3], F[3]], label = "Sun (N)")

    end every 30
end)
#---------------------------------------------------------------------------------------------# """



#---------------------------------------------------------------------------------------------#
# Plot Earth with sun/sat-lit cells highlighted, sun/mag fields
#---------------------------------------------------------------------------------------------#
""" 
for i = 1:1
    pyplot()
    δa = 5; δb = 5;
    lims = 1.25e7 #1.2e7

    t_x = cell_centers_ecef[1:δa:end,1:δb:end,1] * dscale;
    t_y = cell_centers_ecef[1:δa:end,1:δb:end,2] * dscale;
    t_z = cell_centers_ecef[1:δa:end,1:δb:end,3] * dscale;

    scatter3d([0.0], [0.0], [0.0], color = :black, xlim = [-lims, lims], ylim = [-lims, lims], zlim = [-lims, lims], markersize = 5, label = "Grid centers")

    sat_pos = 1.0 .* pos
    plot_sun = 0.00015 * sun;


    T = size(pos,2)
    @gif for i = 1:T
        scatter3d(reshape(t_x, length(t_x)), reshape(t_y, length(t_y)), reshape(t_z, length(t_z)), markersize = 1, markerstrokewidth = 0, color = :gray, xlim = [-lims, lims], ylim = [-lims, lims], zlim = [-lims, lims],  label = false)
        tt_x_lit = cell_centers_ecef[:,:,1] .* reshape(unions[:,i], (180, 288)) .* dscale;
        tt_y_lit = cell_centers_ecef[:,:,2] .* reshape(unions[:,i], (180, 288)) .* dscale;
        tt_z_lit = cell_centers_ecef[:,:,3] .* reshape(unions[:,i], (180, 288)) .* dscale;

        tt_x = tt_x_lit[1:δa:end, 1:δb:end]
        tt_y = tt_y_lit[1:δa:end, 1:δb:end]
        tt_z = tt_z_lit[1:δa:end, 1:δb:end]

        scatter3d!(reshape(tt_x, length(tt_x)), reshape(tt_y, length(tt_y)), reshape(tt_z, length(tt_z)), color = :yellow, label = false)
        scatter3d!([sat_pos[1,i]], [sat_pos[2,i]], [sat_pos[3,i]], markersize = 6, color = :red, xlim = [-lims, lims], ylim = [-lims, lims], zlim = [-lims, lims], label = false)
        
        
        scatter3d!([plot_sun[1,i]], [plot_sun[2,i]], [plot_sun[3,i]], color = :yellow, markersize = 40, label = false)
        
        C_sat = [sat_pos[1,i], sat_pos[2,i], sat_pos[3,i]]
        
        vx, vy, vz = 10 * rN1[:,i] * dscale; # Sun Vector 
        F = [C_sat[1] + vx, C_sat[2] + vy, C_sat[3] + vz];  
        plot3d!([C_sat[1], F[1]], [C_sat[2], F[2]], [C_sat[3], F[3]], linewidth = 2, label = false)

        vx, vy, vz = 8 * rN2[:,i] * dscale; # Mag Vector
        F = [C_sat[1] + vx, C_sat[2] + vy, C_sat[3] + vz];  
        plot3d!([C_sat[1], F[1]], [C_sat[2], F[2]], [C_sat[3], F[3]], color = :blue, label = false)

        scatter3d!(xlabel = "x", ylabel = "y", zlabel = "z", title = "Satellite Orbit")
    end every 30
end
#---------------------------------------------------------------------------------------------# """





###########################################################################################
# Extra MEKF Plots - To be run after main.jl
###########################################################################################

#---------------------------------------------------------------------------------------------#
# Covariance plots
#---------------------------------------------------------------------------------------------#
# """ 
# for i = 1:1
#     using Dates
#     pyplot()

#     function sphere(center, ngrid=25)
#         # Modified from https://discourse.julialang.org/t/drawing-an-ellipsoid-to-visualize-a-tensor/31286/2

#         # Set of all spherical angles:
#         u = range(0, 2pi, length=ngrid);
#         v = range(0, pi, length=ngrid);

#         # Cartesian coordinates that correspond to the spherical angles:
#         # (this is the equation of an ellipsoid):
#         x = [x * y for (x, y) in  Iterators.product(cos.(u), sin.(v))];
#         y = [x * y for (x, y) in Iterators.product(sin.(u), sin.(v))];
#         z = [x * y for (x, y) in Iterators.product(ones(length(u)), cos.(v))];
        
#         return x .+ center[1], y .+ center[2], z .+ center[3];
#     end

#     function ellipsoid(center, coefs, ngrid=25)
#         # From https://discourse.julialang.org/t/drawing-an-ellipsoid-to-visualize-a-tensor/31286/2
        
#         # Radii corresponding to the coefficients:
#         rx, ry, rz = sqrt.(coefs);

#         # Set of all spherical angles:
#         u = range(0, 2pi, length=ngrid);
#         v = range(0, pi, length=ngrid);

#         # Cartesian coordinates that correspond to the spherical angles:
#         # (this is the equation of an ellipsoid):
#         x = [rx * x * y for (x, y) in Iterators.product(cos.(u), sin.(v))];
#         y = [ry * x * y for (x, y) in Iterators.product(sin.(u), sin.(v))];
#         z = [rz * x * y for (x, y) in Iterators.product(ones(length(u)), cos.(v))];
#         return x .+ center[1], y .+ center[2], z .+ center[3];
#     end

#     function simple_covariance_ellipsoid_gif(eig_vals)
#         C = [0 0 0]
#         surface(sphere(C)) # Reset axes
#         lim = 0.2

#         @gif for i ∈ 1:ell 
#             x, y, z = ellipsoid(C, eig_vals[:,i])
#             surface(x, y, z, aspect_ratio=:equal, colorbar = false, title = "Covariance vs Time\n$i", xlim = [-lim, lim], ylim = [-lim, lim], zlim = [-lim, lim])

#         end every 50
#     end

#     function newton_vectors_and_ellipsoid_gif(eig_vals, rN1, rN2)
#         C = [0 0 0]
#         surface(sphere(C)) # Reset axes

#         @gif for i ∈ 1:ell
#             x, y, z = ellipsoid(C, 0.1 * eig_vals[:,i])
#             surface(x, y, z, aspect_ratio=:equal, colorbar = false, title = "Covariance vs Time\n$i", xlim = [-.05, .05], ylim = [-.05, .05], zlim = [-.05, .05])

#             Bx, By, Bz = rN2[:,i] * 0.05 # Mag vectors, scaled purely for visual reasons
#             plot!([0.0, Bx], [0.0, 0.0], [0.0, 0.0], label = "Bx", lw = 1)
#             plot!([0.0, 0.0], [0.0, By], [0.0, 0.0], label = "By", lw = 1)
#             plot!([0.0, 0.0], [0.0, 0.0], [0.0, Bz], label = "Bz", lw = 1)
        
#             Sx, Sy, Sz = rN1[:,i] * 0.05 # Sun vectors, scaled for aesthetics
#             plot!([0.0, Sx], [0.0, 0.0], [0.0, 0.0], label = "Sx", lw = 2, linestyle = :dash, xlim = [-0.05, 0.05], ylim = [-0.05, 0.05], zlim = [-0.05, 0.05])
#             plot!([0.0, 0.0], [0.0, Sy], [0.0, 0.0], label = "Sy", lw = 2, linestyle = :dash)
#             plot!([0.0, 0.0], [0.0, 0.0], [0.0, Sz], label = "Sz", lw = 2, linestyle = :dash)

#         end every 50
#     end

#     function orbiting_gif(eig_vals, rN1, rN2)
#         C = [0 0 0]
#         lim = 1.5 # axis limits
#         _Re = 6378136.3

#         pos_scaled = 1.25 * (pos ./ _Re) # Position normalized by Earth radius (and scaled up a bit for aesthetics)
#         surface(sphere(C)) # Reset axes

#         @gif for i ∈ 1:5700 #ell
#             # Draw the earth
#             x, y, z = sphere(C);
#             surface(x, y, z, aspect_ratio=:equal, c = :grays, colorbar = false, xlim = [-lim, lim], ylim = [-lim, lim], zlim = [-lim, lim])
            
#             # Draw the satellite
#             local C_sat = [pos_scaled[1,i], pos_scaled[2,i], pos_scaled[3,i]] # Center of satellite
#             x, y, z = ellipsoid((C_sat), eig_vals[:,i]); 
#             surface!(x, y, z, aspect_ratio=:equal, colorbar = false)

#             # Plot the Sun vectors
#             vx, vy, vz = 1.5 * rN1[:,i];
#             F = [C_sat[1] - vx, C_sat[2] - vy, C_sat[3] - vz];  # WHY NEGATIVE???   
#             plot3d!([C_sat[1], F[1]], [C_sat[2], F[2]], [C_sat[3], F[3]], label = "Sun (N)")
        
#             # Plot the mag vectors
#             vx, vy, vz = 1.5 * rN2[:,i]; 
#             F = [C_sat[1] + vx, C_sat[2] + vy, C_sat[3] + vz];
#             plot3d!([C_sat[1], F[1]], [C_sat[2], F[2]], [C_sat[3], F[3]], label = "Mag (N)")
        
#             # plot3d!(title = "Combined")
#         end every 50

#     end


#     # Generate eigenvalues and eigenvectors from covariance matrix
#     ell = size(Phist, 3)

#     quaternion_cov = Phist[1:3, 1:3, :]
#     eig_vecs = zeros(3,3,ell)
#     eig_vals = zeros(3, ell)
#     norm_vals = zeros(ell)
#     for i = 1:ell
#         p_quat = quaternion_cov[:,:,i]
#         F = eigen(p_quat)
#         e_vals = F.values;
#         e_vecs = F.vectors

#         e_vecs[:,1] = e_vecs[:,1] * e_vals[1];
#         e_vecs[:,2] = e_vecs[:,2] * e_vals[2];
#         e_vecs[:,3] = e_vecs[:,3] * e_vals[3];

#         eig_vecs[:,:,i] = e_vecs 
#         eig_vals[:,i]  = e_vals 
#         norm_vals[i] = norm(e_vals)
#     end
#     println("Finished generating eigen stuff")

#     # println("Starting #1 at ", now())
#     # display(simple_covariance_ellipsoid_gif(eig_vals))

#     display(newton_vectors_and_ellipsoid_gif(eig_vals, rN1, rN2))

#     # display(orbiting_gif(eig_vals, rN1, rN2))
# end
# #---------------------------------------------------------------------------------------------# """





###############################################################################################
# Other - no other file must be run first
###############################################################################################

#---------------------------------------------------------------------------------------------#
# Plot magnetic field vectors (in testing)
#---------------------------------------------------------------------------------------------#
"""
for i = 1:1
    using Plots
    using SatelliteDynamics
    include("mag_field.jl")

    function get_albedo_cell_centers(lat_step = 1, lon_step = 1.25)
        """
            # Returns the cell centers for the grid covering the surface of the Earth in Cartesian ECEF, to be used in later estimations of Earth's albedo,
            #     by looping through each cell's LLA coordinate and converting to ECEF 

            # Arguments:
            # - lat_step: (Optional) The step size (in degrees) to take across the latitude domain. Defaults to 1*        | Scalar 
            # - lon_step: (Optional) The step size (in degrees) to take across the longitude domain. Defaults to 1.25*    | Scalar

            # Returns:
            # - cells_ecef: Matrix containing [x,y,z] coordinate for each latitude, longitude point.
            #                 Of form [lat, lon, [x,y,z]]                                                                 | [num_lat x num_lon x 3]
        """
        alt = 0#1000e3 # Assume all cells are on surface of earth
        num_lat = Int(round((180 - lat_step) / lat_step) + 1)
        num_lon = Int(round((360 - lon_step) / lon_step) + 1)

        lon_offset = lon_step + (360 - lon_step) / 2   # Centers at 0 (longitude: [1.25, 360] => [-179.375, 179.375])
        lat_offset = lat_step + (180 - lat_step) / 2   # Centers at 0 (latitude:  [1.00, 180] => [-89.5, 89.5])

        cells_ecef = zeros(num_lat, num_lon, 3) # Lat, Lon, [x,y,z]
        for lat = 1:num_lat 
            for lon = 1:num_lon
                geod = [(lon * lon_step - lon_offset), (lat * lat_step - lat_offset), alt]
                ecef = sGEODtoECEF(geod, use_degrees = true)

                cells_ecef[Int(lat), Int(lon), :] = ecef
            end
        end

        return cells_ecef 
    end

    function sphere(center, ngrid=25)
        # Modified from https://discourse.julialang.org/t/drawing-an-ellipsoid-to-visualize-a-tensor/31286/2

        # Set of all spherical angles:
        u = range(0, 2pi, length=ngrid);
        v = range(0, pi, length=ngrid);

        # Cartesian coordinates that correspond to the spherical angles:
        # (this is the equation of an ellipsoid):
        x = [x * y for (x, y) in Iterators.product(cos.(u), sin.(v))];
        y = [x * y for (x, y) in Iterators.product(sin.(u), sin.(v))];
        z = [x * y for (x, y) in Iterators.product(ones(length(u)), cos.(v))];

        
        return x .+ center[1], y .+ center[2], z .+ center[3];
    end

    center = [0 0 0]
    lims = 0.7e7
    plt = surface(6371e3 .* sphere(center), color = :grays)# aspect_ratio = :equal)
    plt = scatter3d!(xlim = [-lims, lims], ylim = [-lims, lims], zlim = [-lims, lims])

    cells = get_albedo_cell_centers()

    δa = 10;
    δb = 10;

    lims = 1.25e7 #1.2e7

    t_x = cells[1:δa:end,1:δb:end,1]; t_x = reshape(t_x, length(t_x));
    t_y = cells[1:δa:end,1:δb:end,2]; t_y = reshape(t_y, length(t_y));
    t_z = cells[1:δa:end,1:δb:end,3]; t_z = reshape(t_z, length(t_z));

    num_points = length(t_x)
    t = Epoch(2019, 1, 1, 12, 0, 0, 0.0);

    @gif for i = 10:num_points
        plt = surface(6371e3 .* sphere(center), color = :grays)# aspect_ratio = :equal)
        plt = scatter3d!(xlim = [-lims, lims], ylim = [-lims, lims], zlim = [-lims, lims])

        for j = (i-9):i
            C = [t_x[j], t_y[j], t_z[j]]
            mag = 80000 .* IGRF13(C, t) 
            plt = plot3d!([C[1], C[1] + mag[1]], [C[2], C[2] + mag[2]], [C[3], C[3] + mag[3]], label = false, color = :blue, linewidth = 1)
        end
        plt = surface!(6371e3 .* sphere(center), color = :grays)

        # if (true) # == 0
        #     plt = plot3d!([C[1], C[1] + mag[1]], [C[2], C[2] + mag[2]], [C[3], C[3] + mag[3]], label = false, color = :blue, linewidth = 1)
        # else
        #     plt = plot3d!([C[1], C[1] + mag[1]], [C[2], C[2] + mag[2]], [C[3], C[3] + mag[3]], label = false, color = :red, linewidth = 1)
        # end
    end every 1
end
#---------------------------------------------------------------------------------------------# """





#---------------------------------------------------------------------------------------------#
# Eclipse Test (Simple test to verify that the eclipse_conical is working)
#---------------------------------------------------------------------------------------------#
# for i = 1:1
    using Plots
    using SatelliteDynamics
    pyplot()

    T = 100
    θs = range(0, 2*pi, length = T);
    sun_0 = [0; 1.451e11; 0];
    sat_0 = [0; 6.928e6;  0];

    Earth = [0; 0; 0];

    sun = zeros(3, T)
    sat = zeros(3, T)
    ecl_hist = zeros(T)

    sun[:, 1] = sun_0;
    sat[:, 1] = sat_0;
    ecl_hist[1] = eclipse_conical(-sat[:,1], sun[:,1])

    for i = 2:T 
        θ = θs[i];
        R = [cos(θ) sin(θ) 0; -sin(θ) cos(θ) 0; 0 0 1];
        sun[:,i] = sun_0;
        sat[:,i] = R * sat_0;
        ecl_hist[i] = eclipse_conical(-sat[:,i], sun[:,i])
    end

    sun = sun * 0.001; # For plotting convenience

    llims = 1.5e7
    rlims = 1.5e8

    @gif for i = 1:T    
        scatter3d([0], [0], [0], color = :green, label = "Earth", markersize = 7)
        scatter3d!([sun[1,i]], [sun[2,i]], [sun[3,i]], color = :yellow, label = "Sun", markersize = 10)
        if ecl_hist[i] == 1
            scatter3d!([sat[1,i]], [sat[2,i]], [sat[3,i]], color = :red, label = "Sat", markersize = 4)
        else
            scatter3d!([sat[1,i]], [sat[2,i]], [sat[3,i]], color = :black, label = "Sat", markersize = 4)
        end
        scatter3d!(xlim = [-llims, rlims], ylim = [-llims, rlims], zlim = [-llims, rlims], camera = (0, 90))

    end every 1
# end
#---------------------------------------------------------------------------------------------# """