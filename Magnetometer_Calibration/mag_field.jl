using LinearAlgebra, SatelliteDynamics

const SD = SatelliteDynamics


function IGRF13(r_eci,epc::Epoch)
    """IGRF 13 model for NED magnetic field vector in nT

    Arguments:
        r_eci: position vector from eci to sc (m)
        epc: time (Epoch)

    Returns:
        B_eci_T: magnetic field in ECI, Teslas (T)

    Comments:
        I used to use SatelliteToolbox for their IGRF, but I wrote my own
    """

    # get decimal date
    date = SD.caldate(epc)
    year = date[1]
    day = SD.day_of_year(epc)
    decimal_date = year + day/365.2425

    # eci and ecef location
    ecef_Q_eci = SD.rECItoECEF(epc)
    eci_Q_ecef = transpose(ecef_Q_eci)

    # get ecef location
    r_ecef = ecef_Q_eci*r_eci

    # long lat geod
    longitude,latitude,altitude = SD.sECEFtoGEOC(r_ecef,use_degrees=true)

    # IGRF
    # SatelliteToolbox v0.7.1
    # B_ned_nT = igrf(decimal_date, norm(r_ecef), latitude, longitude, Val(:geocentric))
    # SatelliteToolbox v0.6.3
    # B_ned_nT = igrf12(decimal_date, norm(r_ecef), latitude, longitude, Val{:geocentric})
    # my own IGRF function
    B_ned_nT = my_igrf_13(decimal_date,norm(r_ecef)/1000,latitude,longitude,13)

    # NED and ECEF DCM
    ecef_Q_ned = ecef_Q_ned_mat(deg2rad(longitude),deg2rad(latitude))

    # conver to eci
    B_eci_nT = eci_Q_ecef*ecef_Q_ned*B_ned_nT

    # convert from nT to T
    return B_eci_nT*1e-9
end




function my_igrf_13(date,alt,lat,elong,order)
"""Truncated IGRF model.

Arguments:
   gh: truncated coefficients
   date: decimal date
   alt: radius from center of earth (km)
   lat: latitude (degrees)
   elong: east longitude (degrees)
   order: order of IGRF model
"""
gh= [
     -29404.8, -1450.9,  4652.5, -2499.6,  2982.0, -2991.6,    # 2020
       1677.0,  -734.6,  1363.2, -2381.2,   -82.1,  1236.2,    # 2020
        241.9,   525.7,  -543.4,   903.0,   809.5,   281.9,    # 2020
         86.3,  -158.4,  -309.4,   199.7,    48.0,  -349.7,    # 2020
       -234.3,   363.2,    47.7,   187.8,   208.3,  -140.7,    # 2020
       -121.2,  -151.2,    32.3,    13.5,    98.9,    66.0,    # 2020
         65.5,   -19.1,    72.9,    25.1,  -121.5,    52.8,    # 2020
        -36.2,   -64.5,    13.5,     8.9,   -64.7,    68.1,    # 2020
         80.6,   -76.7,   -51.5,    -8.2,   -16.9,    56.5,    # 2020
          2.2,    15.8,    23.5,     6.4,    -2.2,    -7.2,    # 2020
        -27.2,     9.8,    -1.8,    23.7,     9.7,     8.4,    # 2020
        -17.6,   -15.3,    -0.5,    12.8,   -21.1,   -11.7,    # 2020
         15.3,    14.9,    13.7,     3.6,   -16.5,    -6.9,    # 2020
         -0.3,     2.8,     5.0,     8.4,   -23.4,     2.9,    # 2020
         11.0,    -1.5,     9.8,    -1.1,    -5.1,   -13.2,    # 2020
         -6.3,     1.1,     7.8,     8.8,     0.4,    -9.3,    # 2020
         -1.4,   -11.9,     9.6,    -1.9,    -6.2,     3.4,    # 2020
         -0.1,    -0.2,     1.7,     3.6,    -0.9,     4.8,    # 2020
          0.7,    -8.6,    -0.9,    -0.1,     1.9,    -4.3,    # 2020
          1.4,    -3.4,    -2.4,    -0.1,    -3.8,    -8.8,    # 2020
          3.0,    -1.4,     0.0,    -2.5,     2.5,     2.3,    # 2020
         -0.6,    -0.9,    -0.4,     0.3,     0.6,    -0.7,    # 2020
         -0.2,    -0.1,    -1.7,     1.4,    -1.6,    -0.6,    # 2020
         -3.0,     0.2,    -2.0,     3.1,    -2.6,    -2.0,    # 2020
         -0.1,    -1.2,     0.5,     0.5,     1.3,     1.4,    # 2020
         -1.2,    -1.8,     0.7,     0.1,     0.3,     0.8,    # 2020
          0.5,    -0.2,    -0.3,     0.6,    -0.5,     0.2,    # 2020
          0.1,    -0.9,    -1.1,     0.0,    -0.3,     0.5,    # 2020
          0.1,    -0.9,    -0.9,     0.5,     0.6,     0.7,    # 2020
          1.4,    -0.3,    -0.4,     0.8,    -1.3,     0.0,    # 2020
         -0.1,     0.8,     0.3,     0.0,    -0.1,     0.4,    # 2020
          0.5,     0.1,     0.5,     0.5,    -0.4,    -0.5,    # 2020
         -0.4,    -0.4,    -0.6,                               # 2020
          5.7,     7.4,   -25.9,   -11.0,    -7.0,   -30.2,    # 2022
         -2.1,   -22.4,     2.2,    -5.9,     6.0,     3.1,    # 2022
         -1.1,   -12.0,     0.5,    -1.2,    -1.6,    -0.1,    # 2022
         -5.9,     6.5,     5.2,     3.6,    -5.1,    -5.0,    # 2022
         -0.3,     0.5,     0.0,    -0.6,     2.5,     0.2,    # 2022
         -0.6,     1.3,     3.0,     0.9,     0.3,    -0.5,    # 2022
         -0.3,     0.0,     0.4,    -1.6,     1.3,    -1.3,    # 2022
         -1.4,     0.8,     0.0,     0.0,     0.9,     1.0,    # 2022
         -0.1,    -0.2,     0.6,     0.0,     0.6,     0.7,    # 2022
         -0.8,     0.1,    -0.2,    -0.5,    -1.1,    -0.8,    # 2022
          0.1,     0.8,     0.3,     0.0,     0.1,    -0.2,    # 2022
         -0.1,     0.6,     0.4,    -0.2,    -0.1,     0.5,    # 2022
          0.4,    -0.3,     0.3,    -0.4,    -0.1,     0.5,    # 2022
          0.4,     0.0,zeros(115)...                         # 2022
        ]

# colat from lat
colat = 90-lat


# Declaration of variables
fn    = 0
gn    = 0
kmx   = 0
ll    = 0
nc    = 0
nmx   = 0
x     = 0.0
y     = 0.0
z     = 0.0
t     = 0.0
tc    = 0.0

# since date is after 2020
t  = date - 2020
tc = 1.0

# gh = gh[3256:end]
# ll  = 3255
ll = 0
nmx = order

nc  = Int(nmx*(nmx+2))
# nc = 195

kmx = Int((nmx+1)*(nmx+2)/2)
# kmx = 105

# allocate
cl          = zeros(eltype(alt),nmx)
sl          = zeros(eltype(alt),nmx)
p           = zeros(eltype(alt),kmx)
q           = zeros(eltype(alt),kmx)


r  = alt
ct          = cos(colat*pi/180)
st          = sin(colat*pi/180)
cl[1]       = cos(elong*pi/180)
sl[1]       = sin(elong*pi/180)
Cd = 1.0
sd = 0.0
l      = 1
m      = 1
n     = 0

ratio = 6371.2/r
  rr    = ratio^2

  # Computation of Schmidt quasi-normal coefficients p and x(=q).
  # =============================================================

  p[1] = 1.0
  p[3] = st
  q[1] = 0.0
  q[3] = ct

  for k = 2:kmx
      # There is no need to check bounds here. The code guarantees that
      # everything is inside the right bounds. This increased the source code
      # performance by 13%.
      @inbounds begin
          if n < m
              m  = 0
              n  = n+1
              rr = rr*ratio
              fn = n
              gn = n-1
          end

          fm = m

          if (m == n)
              if k != 3
                  one   = sqrt(1 - 0.5/fm)
                  j     = k - n - 1
                  p[k]  = one*st*p[j]
                  q[k]  = one*(st*q[j] + ct*p[j])
                  cl[m] = cl[m-1]*cl[1] - sl[m-1]*sl[1]
                  sl[m] = sl[m-1]*cl[1] + cl[m-1]*sl[1]
              end
          else
              gmm   = m^2
              one   = sqrt(fn^2 - gmm)
              two   = sqrt(gn^2 - gmm)/one
              three = (fn + gn)/one
              i     = k - n
              j     = i - n + 1
              p[k]  = three*ct*p[i] - two*p[j]
              q[k]  = three*(ct*q[i] - st*p[i]) - two*q[j]
          end

          # Synthesis of x, y, and z in geocentric coordinates.
          lm  = ll + l
          one = (tc*gh[lm] + t*gh[lm+nc])*rr

          if m != 0
              two   = (tc*gh[lm+1] + t*gh[lm+nc+1])*rr
              three = one*cl[m] + two*sl[m]
              x     = x + three*q[k]
              z     = z - (fn + 1)*three*p[k]

              if st != 0
                  y = y + (one*sl[m] - two*cl[m])*fm*p[k]/st
              else
                  y = y + (one*sl[m] - two*cl[m])*q[k]*ct
              end

              l = l + 2
          else
              x = x + one*q[k]
              z = z - (fn + 1)*one*p[k]
              l = l + 1
          end

          m = m + 1
      end
  end

  # Conversion to coordinate system specified by itype.
  # ===================================================
  one   = copy(x)
  x     = x*Cd +   z*sd
  z     = z*Cd - one*sd

  return 1e6*[x;y;z]
end


function ecef_Q_ned_mat(longitude,latitude)

    # ϕ = geoc[2]
    # λ = geoc[1]
    ϕ = latitude
    λ = longitude

    ecef_Q_ned = [-sin(ϕ)*cos(λ)  -sin(λ)  -cos(ϕ)*cos(λ);
                  -sin(ϕ)*sin(λ)   cos(λ)  -cos(ϕ)*sin(λ);
                   cos(ϕ)          0.0     -sin(ϕ)]

    return ecef_Q_ned
end

function geoc_from_ecef(r_ecef)
    x,y,z = ecef

    lat = atan(z, sqrt(x*x + y*y))
    lon = atan(y, x)

    return lat, lon
end


# check diffability
function diffme()
    # initial time for sim
    epc = Epoch(2020, 1, 1, 12, 0, 0, 0.0)
    r_eci = 1.05*R_EARTH*normalize(randn(3))

    # @show IGRF13(r_eci,epc::Epoch)
    B_fx(x) = IGRF13(x,epc::Epoch)

    # @show FD.jacobian(B_fx,r_eci)

end

diffme();
