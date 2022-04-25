import numpy as np
from satelliteDynamics import *

class EarthAlbedo():

    def __init__(self, data, data_type, start_time, stop_time,
                    radius = 6371.01e3, irradiance = 1366.9):
        self.data = data 
        self.type = data_type
        self.start_time = start_time 
        self.stop_time = stop_time 
        self.earth_radius = radius 
        self.sun_irradiance = irradiance 

    def rad2idx(self, theta, eps, sy, sx):
        """ 
            Transforms location (in radians) to TOMS REFL matrix indices.
            
            Arguments:
            - θ, ϵ: Location of cell center in spherical coordinates (radians)                          | Scalars 
            - sy, sx: Number of latitude, longitude cells that the surface is being divided into        | Scalars 

            Returns: 
            - i, j: TOMS REFL matrix indices of cell corresponding to given spherical coordinates       | Scalars 
        """

        dx = 2 * np.pi / sx # 360* / number of longitude cells 
        dy = np.pi / sy     # 180* / number of latitude cells 

        i = np.around( (np.pi - dy/2 - eps)/dy ) 
        j = np.around( (theta + np.pi - dx/2 )/dx )

        return i, j

    def idx2rad(self, i, j, sy, sx):
        """ 
            Transforms TOMS REFL matrix indices to radians.

            Arguments:
            - i, j: TOMS REFL matrix indices of desired cell                                         | Scalars
            - sy, sx: Number of latitude, longitude cells that the surface is being divided into     | Scalars

            Returns:
            - θ, ϵ: Location of cell center in spherical coordinates (radians)                       | Scalars
        """
        dx = 2 * np.pi / sx
        dy = np.pi / sy

        eps = np.pi - dy/2 - i*dy   
        theta = j * dx - np.pi + dx/2

        return theta, eps

    def gridangle(self, i1, j1, i2, j2, sy, sx):
        """ 
            Calculate the angle between two grid index pairs 

            Arguments:
            - i1, j1: First grid index pair                                       | Scalars
            - i2, j2: Second grid index pair                                      | Scalars
            - sy, sx: Number of latitude and longitude cells, respectively        | Scalars

            Returns:
            - ρ: Angle between the two grid index pairs                           | Scalars
        """

        theta1, phi1 = self.idx2rad(i1, j1, sy, sx)
        theta2, phi2 = self.idx2rad(i2, j2, sy, sx)

        angle = np.sin(phi1) * np.sin(phi2) * np.cos(theta1 - theta2) + np.cos(phi1) * np.cos(phi2)

        # Protect against numerical instability 
        if hasattr(angle, "__len__"): # If it is an array...
            angle[(angle > 1.0) & (angle < 1.0001)] = 1.0
            angle[(angle < -1.0) & (angle > -1.0001)] = -1.0
        else:  # scalar 
            if (angle > 1.0) and (angle < 1.0001):
                angle = 1.0   # Protect against numerical instability
            elif (angle < -1.0) and (angle > -1.0001):
                angle = -1.0

        rho = np.arccos(angle)

        return rho

    def earthfov(self, pos_sph, sy, sx):
        """ 
            Determines the field of view on earth using spherical coordinates.

            Arguments:
            - pos_sph: vector from Earth to the object in question (satellite, sun) 
                            in ECEF frame using spherical coordinates                                      | [3,]
            - sy, sx: Number of latitude and longitude cells, respectively                                 | Scalars
        
            Returns:
            - fov: Cells on Earth's surface that are in the field of view of the given object (sat or sun) | [sy x sx]
        """
        if (pos_sph[2] < self.earth_radius): # LEO shortcut (?)
            pos_sph[2] = pos_sph[2] + self.earth_radius
    
        dx = 2 * np.pi / sx
        dy = np.pi / sy

        # Small circle center 
        theta0 = pos_sph[0]
        phi0 = pos_sph[1]

        rho = np.arccos(self.earth_radius / pos_sph[2]) # FOV on Earth 
        fov = np.zeros([sy, sx])

        js = np.array([j for j in range(sx)])
        for i in range(sy):
            thetas, phis = self.idx2rad(i, js, sy, sx)
            rds = np.arccos(np.sin(phi0) * np.sin(phis) * np.cos(theta0 - thetas) + np.cos(phi0) * np.cos(phis)) # Radial Distance
            fov[i, (rds <= rho)] = 1.0

        return fov

    def cellarea(self, i, j, sy, sx):
        """ 
            Calculate area of TOMS cell for use in albedo calculation.
                
            Arguments:
            - i, j: TOMS REFL matrix indices of desired cell                                      | Scalars
            - sy, sx: Number of latitude and longitude cells, respectively                        | Scalars
            
            Returns:
            - A: Area of cell                                                                     | Scalars
        """
        _d2r = np.pi / 180.0; # Standard degrees to radians conversion

        theta, phi = self.idx2rad(i, j, sy, sx) # Convert to angles (radians)
        dphi = (180.0 / sy) * _d2r
        dtheta = (360.0 / sx) * _d2r

        # Get diagonal points
        phi_max = phi + dphi/2
        phi_min = phi - dphi/2

        A = ((self.earth_radius)**2) * dtheta * (np.cos(phi_min) - np.cos(phi_max))

        return A

    # NOT VECTORIZED
    def get_albedo_cell_centers(self, lat_step = 1, lon_step = 1.25):
        """
            Returns the cell centers for the grid covering the surface of the Earth in Cartesian ECEF, to be used in later estimations of Earth's albedo,
                by looping through each cell's LLA coordinate and converting to ECEF 

            Arguments:
            - lat_step: (Optional) The step size (in degrees) to take across the latitude domain. Defaults to 1*        | Scalar 
            - lon_step: (Optional) The step size (in degrees) to take across the longitude domain. Defaults to 1.25*    | Scalar

            Returns:
            - cells_ecef: Matrix containing [x,y,z] coordinate for each latitude, longitude point.
                            Of form [lat, lon, [x,y,z]]                                                                 | [num_lat x num_lon x 3]
        """
        alt = 0.0 # Assume all cells are on surface of earth
        num_lat = int(np.around((180 - lat_step) / lat_step) + 1)
        num_lon = int(np.around((360 - lon_step) / lon_step) + 1)

        lon_offset = (360 - lon_step) / 2   # Centers at 0 (longitude: [1.25, 360] => [-179.375, 179.375])
        lat_offset = (180 - lat_step) / 2   # Centers at 0 (latitude:  [1.00, 180] => [-89.5, 89.5])
        
        cells_ecef = np.zeros([num_lat, num_lon, 3]) # Lat, Lon, [x,y,z]
        for lat in range(num_lat): 
            for lon in range(num_lon): 
                geod = np.array([(lon * lon_step - lon_offset), (lat * lat_step - lat_offset), alt])
                ecef = sGEODtoECEF(geod, use_degrees = True)

                cells_ecef[int(lat), int(lon), :] = ecef

        return cells_ecef 

    # NOT sure that this one works with scalars now... just use an if statement...?
    def get_diode_albedo(self, albedo_matrix, cell_centers_ecef, surface_normal, sat_pos):
        """ 
            Estimates the effect of Earth's albedo on a specific photodiode (by using the surface normal of that diode)
                = cell_albedo * surface_normal^T * r_g, with r_g as a unit vector in the direction of the grid point on Earth

            Arguments:
            - albedo_matrix: Albedo values for each cell on the Earth's surface         | [num_lat x num_lon] 
            - surface_normal: Photodiode surface normal                                 | [3,]    or    [3, n]
            - sat_pos: Cartesian position of satellite                                  | [3,]

            Returns:
            - diode_albedo: Total effect of albedo on specified photodiode              | Scalar
        """    
        cell_albedos = np.zeros(albedo_matrix.shape)

        rows, cols = albedo_matrix.shape
        diode_albedo = 0.0

        r_gs = cell_centers_ecef - sat_pos 
        r_gs = r_gs / np.linalg.norm(r_gs, axis = 2)[:, :, None]

        cell_albedo = albedo_matrix[:, :, None] * (r_gs @ surface_normal.T)
        cell_albedo[cell_albedo < 0.0] = 0.0

        result = np.sum(cell_albedo, axis = (0, 1))
        return result

    def albedo(self, sat, sun, refl_data): # self.data):
        """
            Determines the Earth's albedo at a satellite for a given set of reflectivity data.
            - Divides Earth into a set of cells
            - Determines which cells are visible by both the satellite and the Sun
            - Determines how much sunlight is reflected by the Earth towards the satellite 

            Arguements:
            - sat: vector from center of Earth to satellite, [m]                              |  [3,]
            - sun: vector from center of Earth to Sun, [m]                                    |  [3,]
            - refl: Refl struct containing the averaged reflectivity data for each cell       | Custom struct

            Returns:
            - cell_albedos: Matrix with each element corresponding to the albedo coming from 
                                the corresponding cell on the surface of the earth            | [num_lat x num_lon]
            - union: Matrix containing indicators to which cells are illuminated by the sun 
                        AND visible by the satellite                                          | [num_lat x num_lon]
        """

        num_lat, num_lon = refl_data.shape 

        # Verify no satellite collision has occured 
        if np.linalg.norm(sat) < self.earth_radius: #CONSTANTS.EARTH_RADIUS
            raise Exception("albedo.m: The satellite has crashed into Earth!")

        # Convert from Cartesian -> Spherical (r θ ϕ)
        sat_sph_t = cartesian2spherical(sat)
        sun_sph_t = cartesian2spherical(sun)

        # Adjust order to match template code to (θ ϵ r), with ϵ = (π/2) - ϕ
        sat_sph = np.array([sat_sph_t[1], (np.pi/2) - sat_sph_t[2], sat_sph_t[0]])
        sun_sph = np.array([sun_sph_t[1], (np.pi/2) - sun_sph_t[2], sun_sph_t[0]])

        # REFL Indices 
        sat_i, sat_j = self.rad2idx(sat_sph[0], sat_sph[1], num_lat, num_lon)
        sun_i, sun_j = self.rad2idx(sun_sph[0], sun_sph[1], num_lat, num_lon)

        # SKIP GENERATING PLOTS FOR NOW
        satellite_fov_cells = self.earthfov(sat_sph, num_lat, num_lon)   # Cells in the satellite field-of-view
        sunlit_cells     =    self.earthfov(sun_sph, num_lat, num_lon)   # Cells lit by the sun

        union = satellite_fov_cells * sunlit_cells # Cells lit by the sun AND visible by the satellite

        cell_albedos = np.zeros([num_lat, num_lon])

        js = np.array([j for j in range(num_lon)])
        for i in range(num_lat):  
            phi_ins = self.gridangle(i, js, sun_i, sun_j, num_lat, num_lon)
            phi_ins[phi_ins > np.pi/2] = np.pi/2 

            E_incs = self.sun_irradiance * self.cellarea(i, js, num_lat, num_lon) * np.cos(phi_ins) 

            grid_thetas, grid_phis = self.idx2rad(i, js, num_lat, num_lon) # Distance to satellite from grid center

            grid_sphericals = np.array([np.ones(num_lon)*self.earth_radius, grid_thetas, np.ones(num_lon)*(np.pi/2) - grid_phis])
            grids = spherical2cartesian(grid_sphericals)

            sat_dists = np.linalg.norm(sat - grids.T, axis = 1) 
            r_norm = ((-grids + sat[:,None])/sat_dists)
            grids_norm = (grids / np.linalg.norm(grids, axis = 0))
            angles = np.sum( r_norm * grids_norm, axis = 0) 
            phi_outs = np.arccos(angles) # Angle to sat from grid

            P_outs = E_incs * refl_data[i, :] * np.cos(phi_outs) / (np.pi * sat_dists**2)

            cell_albedos[i,:] = P_outs * union[i,:] # Zero out if not in union

        return cell_albedos, union

