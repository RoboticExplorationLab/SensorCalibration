import numpy as np
from numpy import sin, cos
import pysofa2  # For testing

# CONSTANTS
MJD2000  = 51544.0   # Modified Julian Date of Jan 2, 2000
MJD_ZERO = 2400000.5   # Offset of modified Julian Days wrt Julian Days
DAYSEC = 86400
DJMIN = -68569.5 
DJMAX = 1e9

# CoordinateTransformations.jl, not SatelliteDynamics
def cartesian2spherical(p):
    # p: [x, y, z] => [r, theta, phi] | radians

    r = np.sqrt(np.dot(p, p))
    theta = np.arctan2(p[1], p[0])
    phi = np.arctan2(p[2], np.sqrt((p[0]**2) + (p[1]**2)))

    return np.array([r, theta, phi])

def spherical2cartesian(p):
    # p: [r theta phi] => [x, y, z] | radians

    sTheta = np.sin(p[1])
    cTheta = np.cos(p[1])

    sPhi = np.sin(p[2])
    cPhi = np.cos(p[2])

    r = p[0]

    return np.array([ r*cTheta*cPhi,
                        r*sTheta*cPhi,
                        r*sPhi   ])


# In SatelliteDynamics
class Epoch():
    # Python implementation of parts of SatelliteDynamics.jl    

    def align_epoch_data(self, days, seconds, nanoseconds):
        # Align nanoseconds to [0.0, 1.0e9)
        while nanoseconds < 0.0:
            nanoseconds += 1.0e9 
            seconds -= 1 

        while nanoseconds >= 1.0e9:
            nanoseconds -= 1.0e9 
            seconds += 1 

        # Align seconds to [0, 86400) 
        while seconds < 0:
            seconds += 86400 
            days -= 1 
        
        while seconds >= 86400:
            seconds -= 86400 
            days += 1 

        return days, seconds, nanoseconds

    def __init__(self, day, sec, nano, tsys):
        self.days = day 
        self.seconds = sec 
        self.nanoseconds = nano
        self.tsys = tsys 

    @classmethod  # Allows this to be run before an object is created (AKA overloading __init__)
    def ymd(cls, year, month, day, hour = 0, minute = 0, second = 0, nanosecond = 0.0, tsys = "TAI"):
        """ Initializes an Epoch from year, month, and day """

        if not valid_time_system(tsys):
            raise Exception("Invalid time system $(String(tsys))")

        # Convert date into days and fractional days in desired time system
                # _sofa.Dtf2d
        status, jd, fd = iauDtf2d("TAI", year, month, day, hour, minute, 0)
        if status < 0:
            print("Year: {}, Month: {}, Day:{}, Hour: {}, Minute: {}".format(year, month, day, hour, minute))
            raise Exception("jd < 0!")
        jd2, fd2 = pysofa2.Dtf2d("TAI", year, month, day, hour, minute, 0)
        if (jd != jd2) or (fd != fd2):
            print("ERROR with Dtf2d!")
            print("JD: {}  vs  {}".format(jd, jd2))
            print("FD: {}  vs  {}".format(fd, fd2))

        # Get time system offset based on days and fractional days using SOFA
        tsys_offset     = time_system_offset(jd, fd, tsys, "TAI")
        if tsys_offset != 0:
            print("Error with offset!")
        foffset, offset = modf(tsys_offset)

        # Ensure days is an integer number, add entire fractional component to the fractional days variable
        fdays, days = modf(jd)
        fd          = fd + fdays

        # Convert fractional days into total seconds still retaining fractional part
        seconds = fd * 86400.0
        fsecs, secs = modf(seconds)

        # Now trip second of the fractional part
        fsecond, second = modf(second)

        # Add the integer parts together
        seconds = secs + second + offset

        nanoseconds = nanosecond + fsecs*1e9 + fsecond*1e9 + foffset*1e9
        nanoseconds = float(nanoseconds)

        # Ensure type consistency of input:
        days        = int(days)
        seconds     = int(seconds)
        nanoseconds = float(nanoseconds)

        days, seconds, nanoseconds = cls.align_epoch_data(cls, days, seconds, nanoseconds)

        # Parse inputs into Epoch internal storage
        return cls(days, seconds, nanoseconds, tsys)

    # DONT THINK I DID THIS RIGHT...?
    # TODO do i need a similar method for not epc + epc but epc + sec?
    def __add__(self, epc):
        if not (epc.tsys == self.tsys):
            raise Exception("Different T systems used!")

        ns = self.nanoseconds + epc.nanoseconds 
        if ns >= 1e9:
            ns -= 1e9
            s = self.seconds + epc.seconds + 1
        else:
            s = self.seconds + epc.seconds

        if s >= 86400:
            s -= 86400 
            d = self.days + epc.days + 1 
        else:
            d = self.days + epc.days 

        return Epoch(d, s, ns, tsys)

    def add_seconds(self, sec):
        s = self.seconds
        d = self.days 
        ns = self.nanoseconds

        if sec < 1:
            ns += (sec * 1e9)
        else:
            s += sec

        if ns >= 1e9:
            ns -= 1e9 
            s += 1
        if s >= 86400:
            s -= 86400 
            d += 1

        return Epoch(d, s, ns, self.tsys)

    # Havent implemented UTC conversion yet
    # def __str__(self):
    #     year, month, day, hour, minute, second, nanodseconds = caldate(self, tsys="UTC")
    #     s = "Epoch({:2}-{:2}-{:2}T{:2}:{:2}:{:2}.{:3f}Z)".format(year, month, day, hour, minute, second, np.floor(nanodseconds / 1.0e6))
    #     return s


######################
#  PRIMARY FUNCTIONS #
######################

def sun_position(epc):
    # Constants
    mjd_tt  = mjd(epc, tsys="TT")      # MJD of epoch in TT
    epsilon = 23.43929111*np.pi/180.0     # Obliquity of J2000 ecliptic
    T       = (mjd_tt-MJD2000)/36525.0 # Julian cent. since J2000

    # Variables
    # Mean anomaly, ecliptic longitude and radius
    M = 2.0*np.pi * modf(0.9931267 + 99.9973583*T)[0]                 # [rad]
    L = 2.0*np.pi * modf(0.7859444 + M/(2.0*np.pi) + (6892.0*np.sin(M)+72.0*np.sin(2.0*M)) / 1296.0e3)[0] # [rad]
    r = 149.619e9 - 2.499e9*np.cos(M) - 0.021e9*np.cos(2*M)           # [m]

    # Equatorial position vector
    p_sun = Rx(-epsilon) @ np.array([r*np.cos(L), r*np.sin(L), 0.0])

    return p_sun
 
def sGEODtoECEF(geod, use_degrees = False):

    WGS84_a = 6378137.0  # WGS-84 semi-major axis (_Re?)
    WGS84_f = 1.0/298.257223563   # WGS-84 flattening 
    ECC2 = WGS84_f * (2.0 - WGS84_f) # Square of eccentricity
    # Extract lat and lon
    lon = geod[0]
    lat = geod[1]
    
    # Handle non-explict use-degrees
    if geod.shape[0] == 3:
        alt = geod[2]
    else:
        alt = 0.0

    # Convert input to radians
    if use_degrees:
        lat = lat*np.pi/180.0
        lon = lon*np.pi/180.0
    

    # Check validity of input
    if (lat < -np.pi/2) or (lat > np.pi/2):
        raise Exception("Lattiude, $lat, out of range. Must be between -90 and 90 degrees.")

    # Compute Earth-fixed position vector
    N = WGS84_a / np.sqrt(1.0 - ECC2*(np.sin(lat)**2))
    x = (N + alt) * np.cos(lat) * np.cos(lon)
    y = (N + alt) * np.cos(lat) * np.sin(lon)
    z = ((1.0 - ECC2) * N + alt) * np.sin(lat)
    
    return np.array([x, y, z])

# UNVERIFIED
def sECEFtoGEOC(ecef, use_degrees = False):
    print("sECEFtoGEOC not implemented!")
    pass


########################
# reference_systems.jl #
########################

#######################################################
def iauNut00b():
    pass 

def iauPr00():
    pass 

def iauObl80():
    pass 

def iauBp00():
    pass 

def iauCr():
    pass 

def iauNumat():
    pass 

def iauRxr():
    pass 

def iauPn00(date1, date2, dpsi, deps):
    
    # IAU 2000 precession-rate adjustments. 
    dpsipr, depspr = iauPr00(date1, date2)

    # Mean obliquity, consistent with IAU 2000 precession-nutation. */
    epsa = iauObl80(date1, date2) + depspr

    # Frame bias and precession matrices and their product.
    iauBp00(date1, date2, rb, rp, rbpw)
    iauCr(rbpw, rbp)

    # Nutation matrix.
    iauNumat(*epsa, dpsi, deps, rnw)
    iauCr(rnw, rn)

    # Bias-precession-nutation matrix (classical)
    iauRxr(rnw, rbpw, rbpn)

def iaupn00b(date1, date2):
    iauNut00b(date1, date2, dpsi, deps)

def iauPnm00b(date1, date2):
    iauPn00b(date1, date2)

def iauBpn2xy(rpbn):
    return rpbn[2, 0], rpbn[2, 1]

def iauS00(date1, date2, x, y):
    pass 

def iauXys00b(date1, date2):
    rbpn = iauPnm00b(date1, date2)
    x, y = iauBpn2xy(rbpn)
    s, x, y = iauS00(date1, date2)

    return s, x, y

def iauC2ixys(x, y, s):
    r2 = x*x + y*y 
    e  = np.arctan2(y, x) if (r2 > 0.0) else 0.0 
    d  = np.arctan( np.sqrt(r2 / (1.0 - r2)) )

    r = iauIr()
    r = iauRz(e, r)
    r = iauRy(d, r)
    r = iauRz( -(e+s), r)

    return r

def bias_precession_nutation(epc):
    """
        Computes the Bias-Precession-Nutation matrix transforming the GCRS to the 
        CIRS intermediate reference frame. This transformation corrects for the 
        bias, precession, and nutation of Celestial Intermediate Origin (CIO) with
        respect to inertial space.

        Arguments:
        - `epc::Epoch`: Epoch of transformation

        Returns:
        - `rc2i::Matrix{<:Real}`: 3x3 Rotation matrix transforming GCRS -> CIRS
    """
    # Constants of IAU 2006A transofrmation
    DMAS2R =  4.848136811095359935899141e-6 / 1.0e3
    dx06   =  0.0001750 * DMAS2R
    dy06   = -0.0002259 * DMAS2R

    # Compute X, Y, s terms using low-precision series terms
    x, y, s = iauXys00b(MJD_ZERO, mjd(epc, tsys="TT"))

    # Apply IAU2006 Offsets
    x += dx06
    y += dy06

    # Compute transformation and return
    rc2i = iauC2ixys(x, y, s)

    return rc2i

def iauAnp(a):
    w = a % D2PI 
    if (w < 0):
        w += D2PI

    return w

# TODO DJ00, D2PI
def iauEra00(dj1, dj2):

    # Days since fundamental epoch
    if (dj1 < dj2):
        d1 = dj1
        d2 = dj2 
    else:
        d1 = dj2 
        d2 = dj1 
    
    t = d1 + d2 - DJ00

    # Fractional part of T (days)
    f = (d1 % 1.0) + (d2 % 1.0)

    # Earth rotation angle at this UT1 
    theta = iauAnp(D2PI * (f + 0.7790572732640 + 0.00273781191135448 * t))

    return theta

def iauRz(psi, r):
    s, c = np.sin(theta), np.cos(theta)

    r[0, 0] =  c*r[0][0] + s*r[1][0]
    r[0, 1] =  c*r[0][1] + s*r[1][1]
    r[0, 2] =  c*r[0][2] + s*r[1][2]
    r[1, 0] = -s*r[0][0] + c*r[1][0]
    r[1, 1] = -s*r[0][1] + c*r[1][1]
    r[1, 2] = -s*r[0][2] + c*r[1][2]

    return r

# TODO need to implement UT1 in mjd :( 
def earth_rotation(epc):
    """
        Computes the Earth rotation matrix transforming the CIRS to the TIRS
        intermediate reference frame. This transformation corrects for the Earth
        rotation.

        Arguments:
        - `epc::Epoch`: Epoch of transformation

        Returns:
        - `r::Matrix{<:Real}`: 3x3 Rotation matrix transforming CIRS -> TIRS
    """
    # Compute Earth rotation angle

    # era = iauEra00(MJD_ZERO, mjd(epc, tsys="UT1"))

    # # Rotate Matrix and return
    # r = iauRz(era, Matrix{Float64}(I, 3, 3))

    return r

def iauIr(r, dim = 3):
    """ returns an identity matrix """
    return np.eye(dim)

# These three just seem to be standard rotation matrices...
def iauRx(phi, r):
    s, c = np.sin(phi), np.cos(phi)

    r = np.zeros([3,3])

    r[1, 0] = c * r[1, 0] + s * r[2, 0]
    r[1, 1] = c * r[1, 1] + s * r[2, 1]
    r[1, 2] = c * r[1, 2] + s * r[2, 2]
    r[2, 0] = s * r[1, 0] + c * r[2, 0] 
    r[2, 1] = s * r[1, 1] + c * r[2, 1]
    r[2, 2] = s * r[1, 2] + c * r[2, 2]

    return r 

def iauRy(theta, r):
    s, c = np.sin(theta), np.cos(theta)

    r[0, 0] = c * r[0, 0] - s * r[2, 0]
    r[0, 1] = c * r[0, 1] - s * r[2, 1]
    r[0, 2] = c * r[0, 2] - s * r[2, 2]
    r[2, 0] = s * r[0][0] + c * r[2][0]
    r[2, 1] = s * r[0][1] + c * r[2][1]
    r[2, 2] = s * r[0][2] + c * r[2][2]

    return r

def iauPom00(xp, yp, sp):
    rpom = iauIr() 
    rpom = iauRz(sp,  rpom)
    rpom = iauRy(-xp, rpom)
    rpom = iauRx(-yp, rpom)
    
    return rpom

def POLE_LOCATOR():
    pass 

def polar_motion(epc):
    """
        Computes the Earth rotation matrix transforming the TIRS to the ITRF reference 
        frame.

        # Arguments
        - `epc::Epoch`: Epoch of transformation

        # Returns
        - `rpm::Matrix{<:Real}`: 3x3 Rotation matrix transforming TIRS -> ITRF
    """
    xp, yp = POLE_LOCATOR(mjd(epc, tsys="UTC"))

    # Compute transformation and return
    rpm = iauPom00(xp, yp, iauSp00(MJD_ZERO, mjd(epc, tsys="TT")))

    return rpm
#########################################


# NOTE - UNTESTED & simplified by only using z rotation
def rECItoECEF(epc):
    """
        Computes the combined rotation matrix from the ECI (Earth Centered Inertial) to the ECEF (Earth-Centered Earth-Fixed)
        reference frame. Applies correction for Earth-rotation.

        # Arguments
        - `epc::Epoch`: Epoch of transformation

        # Returns
        - `r::Matrix{<:Real}`: 3x3 Rotation matrix transforming ECI -> ECEF
    """
    # Based on Kevin's SpacecraftSim code 
    mjd_current = mjd(epc, tsys = "UTC")
    theta = 1.0047517553493037 + 6.300387486754831 * mjd_current
    theta = theta % (2 * np.pi)
    # rotation 
    return np.array([[cos(theta),   sin(theta), 0],
                     [-sin(theta),  cos(theta), 0],
                     [0,            0,          1]])


######################
#  HELPTER FUNCTIONS #
######################

def valid_time_system(tsys):
        return True 

def mjd(epc, tsys= "TT"): #epc.tsys):
        jd = epc.days 
        fd = (epc.seconds + epc.nanoseconds/1.0e9)/ 86400.0
        offset = time_system_offset(jd, fd, "TAI", tsys)

        return (epc.days + (epc.seconds + epc.nanoseconds/1.0e9 + offset)/86400.0) - MJD_ZERO

# ONLY INCLUDES NECESSARY CONVERSIONS!
def time_system_offset(jd, fd, tsys_src, tsys_dest):
    GPS_TAI  = -19.0    # Offset of GPS wrt TAI
    TAI_GPS  = -GPS_TAI # Offset of TAI wrt GPS
    TT_TAI   = 32.184   # Offset of TT wrt TAI
    TAI_TT   = -TT_TAI  # Offset of TAI wrt TT

    if tsys_src == tsys_dest:
        return 0.0 # No offset 

    offset = 0.0

    # # Convert to TAI 
        # if tsys_src == "GPS":
        #     offset += self.TAI_GPS
        # elif tsys_src == "TT":
        #     offset += self.TAI_TT 
    if tsys_src == "UTC":
        status, iy, im, id, ihmsf = iauD2dtf("UTC", 6, jd, fd) # Returns TAI-UTC
        status, dutc = iauDat(iy, im, id, (ihmsf[0]*3600 + ihmsf[1]*60 + ihmsf[2] + ihmsf[3]/1e6)/86400.0)
        offset += dutc
        # elif tsys_src == "UT1":
        #     # Convert UT1 -> UTC
        #     offset -= UT1_UTC((jd - MJD_ZERO) + fd)

        #     # Convert UTC -> TAI
        #     status, iy, im, id, ihmsf = iauD2dtf("UTC", 6, jd, fd + offset) # Returns TAI-UTC
        #     status, dutc = iauDat(iy, im, id, (ihmsf[0]*3600 + ihmsf[1]*60 + ihmsf[2] + ihmsf[3]/1e6)/86400.0)
        #     offset += dutc
    if tsys_src == "TAI":
        # Do nothing in this case
        offset = 0.0 
    else:
        raise Exception("Only current valid source is TAI")

    # Covert from TAI to source
        # if tsys_dest == "GPS":
        #     offset += self.GPS_TAI
        # elif tsys_dest == "TT":
        #     offset += self.TT_TAI
        # elif tsys_dest == "UTC":
        #     # Initial UTC guess
        #     u1, u2 = jd, fd + offset/86400.0

        #     # Iterate to get the UTC time
        #     for i in range(3):
        #         status, d1, d2 = iauUtctai(u1, u2)

        #         # Adjust UTC guess
        #         u1 += jd - d1
        #         u2 += fd + offset/86400.0 - d2

        #     # Compute Caldate from two-part date
        #     status, iy, im, id, ihmsf = iauD2dtf("UTC", 6, u1, u2)

        #     status, dutc = iauDat(iy, im, id, (ihmsf[1]*3600 + ihmsf[2]*60 + ihmsf[3] + ihmsf[4]/1e6)/86400.0)
        #     offset -= dutc

        # elif tsys_dest == "UT1":
        #     # Initial UTC guess
        #     u1, u2 = jd, fd + offset/86400.0

        #     # Iterate to get the UTC time
        #     for i in range(3):
        #         status, d1, d2 = iauUtctai(u1, u2)

        #         # Adjust UTC guess
        #         u1 += jd - d1
        #         u2 += fd + offset/86400.0 - d2
        #     # Compute Caldate from two-part date
        #     status, iy, im, id, ihmsf = iauD2dtf("UTC", 6, u1, u2)

        #     status, dutc = iauDat(iy, im, id, (ihmsf[1]*3600 + ihmsf[2]*60 + ihmsf[3] + ihmsf[4]/1e6)/86400.0)
        #     offset -= dutc

        #     # Convert UTC to UT1
        #     offset += UT1_UTC(u1 + u2 + offset/86400.0 - MJD_ZERO)
        # elif tsys_dest == "TAI":
        #     asdf = 0
        # else:
        #     raise Exception("Invalid option!")
    if tsys_dest == "TT":
        offset += TT_TAI 
    else:
        raise Exception("Only (currently) valid tsys destination is TT")

    return offset


def Rx(angle, use_degrees = False):
    # Rotation matrix about x-axis 
    if use_degrees:
        angle *= np.pi / 180.0 
    
    c = np.cos(angle)
    s = np.sin(angle)

    return np.array([ [1.0, 0.0, 0.0],
                        [0.0,   c,   s],
                        [0.0,   -s,  c]])

def modf(x):
    intgr = np.floor(np.abs(x))
    dcml = np.abs(x) - intgr

    return (-dcml, -intgr) if x < 0.0 else (dcml, intgr)

def iauCal2jd(year: int, month: int, day: int):
    # from IAU's SOFA

    YEAR_MIN = -4799 # no earlier than 4800 BC
    DJM0 = 2400000.5
    mLengths = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])

    j = 0
    djm0 = DJM0
    djm = 0
    # Validate year, month
    if (year < YEAR_MIN):
        return -1, djm0, djm 
    if (month < 1) or (month > 12):
        return -2, djm0, djm 
    
    leap_year = True if ((month == 2) and not (year % 4) and ((year % 100) or not (year % 400))) else False 

    # Validate day, accounting for leap years 
    if ((day < 1) or (day > (mLengths[month-1] + leap_year))):
        j = -3 

    # Return result -> NOTE that there are lots of forced types so as to match the SOFA values
    my = int((month - 14) / 12.0)
    iypmy = int(year + my)  
    djm0 = DJM0

    djm = int(1461 * (iypmy + 4800) / 4) + int((367 * (month - 2 - 12 * my)) / 12) \
          - int((3 * (iypmy + 4900) / 100) / 4) + int(day - 2432076)

    return j, djm0, djm

def iauDtf2d(tsys, year, month, day, hour, minute, second):
    # from IAU's SOFA
    
    # Today's Julian Day Number
    js, dj, w = iauCal2jd(year, month, day)
    dj2, w2   = pysofa2.Cal2jd(year, month, day)

    if (dj != dj2) or (w != w2):
        print("Error with Cal2jd!")
        print("W : {}  vs  {}".format(w,  w2))


    if js != 0:
        raise Exception("Error with iauCal2jd")
    
    dj += w
    day = DAYSEC 
    seclim = 60.0 

    """if not (tsys == "UTC"):

        # TAI - UTC at 0h today
        js, dat0 = iauDat(year, month, day, 0.0)
        if js:
            return js 

        # TAI - UTC at 12h today (to detect drift)
        js, dat12 = iauDat(year, month, day, 0.5)
        if js < 0:
            return js 
        
        # TAI - UTC at 0h tomorrow (to detect jumps)
        js, year2, month2, day2, w = iauJd2cal(dj, 1.5, w)
        if js:
            return js 
        
        js, dat24 = iauDat(year2, month2, day2, 0.0)
        if js < 0:
            return js 

        # Check for a sudden change in TAI-UTC between today and tomorrow
        dleap = dat24 - (2.0 * dat12 - dat0)

        # If leap second day, correct
        day += dleap
        if (hour == 23) and (minute == 59):
            seclim += dleap
    else:
        raise Exception(("Invalid tsys in iauDtf2d"))"""

    # Validate the time 
    if (hour >= 0) and (hour <= 23):
        if (minute >= 0) and (minute <= 59):
            if second >= 0:
                if second >= seclim:
                    js += 2 
            else:
                js = -6 
        else:
            js = -5
    else:
        js = -4

    if js < 0:
        return js, 0.0, 0.0

    time = (60.0 * (60.0 * hour + minute) + second) / day

    return js, dj, time

def iauD2tf(ndp: int, days: float):
    sign = "+" if days >= 0.0 else "-"
    a = DAYSEC * np.abs(days) # Interval in seconds 

    # Pre-round if resolution is coarser than 1s 
    if ndp < 0:
        nrs = 1 
        for n in range(-ndp):
            nrs = (nrs * 6) if (n == 2 or n == 4) else (nrs * 10)
        rs = float(nrs)
        w = a / rs 
        a = rs * np.round(w)

    nrs = 1
    for n in range(ndp):
        nrs *= 10
    rs = float(nrs)
    rm = rs * 60.0 
    rh = rm * 60.0

    a = np.round(rs * a)

    ah = int(a / rh) 
    a -= (ah * rh)

    am = int(a / rm)
    a -= (am * rm)

    ash = int(a / rs) 
    af = a - ash * rs 

    hmsf = np.array([int(ah), int(am), int(ash), int(af)])

    return sign, hmsf

def iauJd2cal(dj1, dj2):
    dj = dj1 + dj2 
    if (dj < DJMIN) or (dj > DJMAX): 
        return -1, 0.0, 0.0, 0.0, 0.0

    if (dj1 >= dj2):
        d1 = dj1 
        d2 = dj2  
    else:
        d1 = dj2 
        d2 = dj1 
    d2 -= 0.5

    # Separate day and fraction 
    f1 = d1 % 1.0 
    f2 = d2 % 1.0
    f = (f1 + f2) % 1.0 
    if f < 0.0:
        f += 1.0 
    
    d = np.round(d1 - f1) + np.round(d2 - f2) + np.round(f1 + f2 - f)
    jd = int(np.round(d) + 1)

    # Express in Gregorian calendar (trying to match SOFA by following types)
    l = int(jd + 68569)
    n = int((4 * l) / 146097)
    l -= int((146097 * n + 3) / 4)
    i = int((4000 * (l + 1)) / 1461001)
    l -= int((1461 * i) / 4 - 31)
    k = int((80 * l) / 2447)
    ids = int(l - int((2447 * k) / 80))
    l = int(k / 11)
    im = int(k + 2 - 12 * l)
    iy = int(100 * (n - 49) + i + l)
    fd = f

    return 0.0, iy, im, ids, fd

# Skips all the UTC stuff
def iauD2dtf(scale, ndp, d1, d2):
    # From IAU's SOFA

    error = False # Flag to track errors 

    a1, b1 = d1, d2

    # Provisional calendar date
    js, y1, m1, d1, fd = iauJd2cal(a1, b1)
    if js:
        error = True

    # Is this a leap second day?
    leap = 0
    if (scale != "UTC"):

        # TAI-UTC at 0h today 
        js, dat0  = iauDat(y1, m1, d1, 0.0)
        if js < 0:
            error = True

        # TAI-UTC at 12h today (to detect drift)
        js, dat12 = iauDat(y1, m1, d1, 0.5)
        if js < 0:
            error = True 

        # TAI-UTC at 0h tomorrow (to detect jumps)
        js, y2, m2, d2, w = iauJd2cal(a1+1.5, b1-fd)
        if js < 0:
            error = True 
        
        js, dat24 = iauDat(y2, m2, d2, 0.0)
        if js < 0:
            error = True 

        # Any sudden change in TAI-UTC (seconds) 
        dleap = dat24 - (2.0 * dat12 - dat0)

        # If leap second day, scale the fraction of a day into SI 
        leap = (dleap != 0.0)
        if leap:
            fd += fd * dleap / DAYSEC 


    # Provisional Time of Day
    _, hmsf= iauD2tf(ndp, fd)

    # Has (rounded) time gone past 24h?
    if hmsf[0] > 23:
        js, y2, m2, d2, w = iauJd2cal(a1+1.5, b1-fd)
        if js:
            return -1, 0, 0, 0, 0

        # is today a leap second day?
        if not leap:
            # No - use 0h tomorrow 
            y1 = y2
            m1 = m2 
            d1 = d2 
            hmsf[0] = 0 
            hmsf[1] = 0 
            hmsf[2] = 0
        
        else:
            # Yes - are we past the leap second itself?

            if hmsf[2] > 0:
                # Yes. Use tomorrow but allow for leap second 
                y1 = y2 
                m1 = m2 
                d1 = d2 
                hmsf[0] = 0 
                hmsf[1] = 0 
                hmsf[2] = 0

            else:
                # No - use 23 59 60 ... today
                hmsf[0] = 23
                hmsf[1] = 59 
                hmsf[2] = 60

            # If rounding to 10s or coarser always go up to new day 
            if ((ndp < 0) and (hmsf[2] == 60)):
                y1 = y2 
                m1 = m2 
                d1 = d2 
                hmsf[0] = 0
                hmsf[1] = 0 
                hmsf[2] = 0 
                
        # Assume no leap because no UTC...
        year = y1 # y2?
        month = m1 
        day = d1 
        # hmsf[0:3] = np.zeros(3)

        # More assumptions because no UTC...
    else:
        year, month, day = y1, m1, d1 

    return js, year, month, day, hmsf



def caldate(epc, tsys = "TAI"): #epc.tsys):
    # Returns Gregorian Calendar date
    jd = epc.days 
    fd = (epc.seconds + epc.nanoseconds / 1.0e9) / 86400.0
    offset = time_system_offset(jd, fd, "TAI", tsys)

    status, year, month, day, hmsf = iauD2dtf("TAI", 9, epc.days, (epc.seconds + offset + epc.nanoseconds/1.0e9)/86400.0)

    return year, month, day, hmsf[0], hmsf[1], hmsf[2], hmsf[3]

def day_of_year(epc, tsys = "TAI"): #epc.tsys):
    year, month, day, hour, minute, second, _ = caldate(epc, tsys = tsys)
    mjd0 = caldate_to_mjd(year, 1, 1) 

    jd = epc.days 
    fd = (epc.seconds + epc.nanoseconds / 1.0e9) / 86400.0
    offset = time_system_offset(jd, fd, "TAI", tsys)

    mjd = (epc.days + (epc.seconds + epc.nanoseconds/1.0e9 + offset)/86400.0) - MJD_ZERO

    doy = mjd - mjd0 + 1.0

    return doy

def caldate_to_mjd(year, month, day, hour = 0, minute = 0, second = 0.0, nanoseconds = 0.0):
    status, jd, fd = iauDtf2d("TAI", year, month, day, hour, minute, second + nanoseconds/1.0e9)

    mjd = (jd - MJD_ZERO) + fd
    return mjd


######################
#  IGRF13  FUNCTIONS #
######################
# From Kevin
# date MUST be after 2020

def IGRF13(r_eci, epc):
    """IGRF 13 model for NED magnetic field vector in nT

    Arguments:
        r_eci: position vector from eci to sc (m)
        epc: time (Epoch)

    Returns:
        B_eci_T: magnetic field in ECI, Teslas (T)

    Comments:
        I used to use SatelliteToolbox for their IGRF, but I wrote my own - Kevin
    """
    date = caldate(epc)
    year = date[0]
    day = day_of_year(epc)
    decimal_date = year + day / 365.2425

    # ECI & ECEF Location 
    ecef_Q_eci = rECItoECEF(epc)
    eci_Q_ecef = ecef_Q_eci.T 

    # get ECEF location 
    r_ecef = ecef_Q_eci * r_eci 
    
    # Long lat geod 
    longitude, latitude, altitude = sECEFtoGEOC(r_ecef, ues_degrees = True)

    # IGRF
    # SatelliteToolbox v0.7.1
    # B_ned_nT = igrf(decimal_date, norm(r_ecef), latitude, longitude, Val(:geocentric))
    # SatelliteToolbox v0.6.3
    # B_ned_nT = igrf12(decimal_date, norm(r_ecef), latitude, longitude, Val{:geocentric})
    # my own IGRF function - Kevin
    B_ned_nT = my_igrf_13(decimal_date, np.linalg.norm(r_ecef)/1000,latitude,longitude,13)

    # NED & ECEF DCM 
    ecef_Q_ned = ecef_Q_ned_mat(np.radians(longitude), np.radians(latitude))

    # Convert to ECI 
    B_eci_nT = eci_Q_ecef * ecef_Q_ned * B_ned_nT 

    # convert from nT to T 
    return B_eci_nT * 1e-9

def my_igrf_13(date, alt, lat, elong, order):
    """Truncated IGRF model.

    Arguments:
    gh: truncated coefficients
    date: decimal date
    alt: radius from center of earth (km)
    lat: latitude (degrees)
    elong: east longitude (degrees)
    order: order of IGRF model
    """
    gh = np.array([-29404.8, -1450.9,  4652.5, -2499.6,  2982.0, -2991.6,   # 2020
                1677.0,  -734.6,  1363.2, -2381.2,   -82.1,  1236.2,    # 2020
                241.9,   525.7,  -543.4,   903.0,   809.5,   281.9,     # 2020
                86.3,  -158.4,  -309.4,   199.7,    48.0,  -349.7,      # 2020
                -234.3,   363.2,    47.7,   187.8,   208.3,  -140.7,    # 2020
                -121.2,  -151.2,    32.3,    13.5,    98.9,    66.0,    # 2020
                65.5,   -19.1,    72.9,    25.1,  -121.5,    52.8,      # 2020
                -36.2,   -64.5,    13.5,     8.9,   -64.7,    68.1,     # 2020
                80.6,   -76.7,   -51.5,    -8.2,   -16.9,    56.5,      # 2020
                2.2,    15.8,    23.5,     6.4,    -2.2,    -7.2,       # 2020
                -27.2,     9.8,    -1.8,    23.7,     9.7,     8.4,     # 2020
                -17.6,   -15.3,    -0.5,    12.8,   -21.1,   -11.7,     # 2020
                15.3,    14.9,    13.7,     3.6,   -16.5,    -6.9,      # 2020
                -0.3,     2.8,     5.0,     8.4,   -23.4,     2.9,      # 2020
                11.0,    -1.5,     9.8,    -1.1,    -5.1,   -13.2,      # 2020
                -6.3,     1.1,     7.8,     8.8,     0.4,    -9.3,      # 2020
                -1.4,   -11.9,     9.6,    -1.9,    -6.2,     3.4,      # 2020
                -0.1,    -0.2,     1.7,     3.6,    -0.9,     4.8,      # 2020
                0.7,    -8.6,    -0.9,    -0.1,     1.9,    -4.3,       # 2020
                1.4,    -3.4,    -2.4,    -0.1,    -3.8,    -8.8,       # 2020
                3.0,    -1.4,     0.0,    -2.5,     2.5,     2.3,       # 2020
                -0.6,    -0.9,    -0.4,     0.3,     0.6,    -0.7,      # 2020
                -0.2,    -0.1,    -1.7,     1.4,    -1.6,    -0.6,      # 2020
                -3.0,     0.2,    -2.0,     3.1,    -2.6,    -2.0,      # 2020
                -0.1,    -1.2,     0.5,     0.5,     1.3,     1.4,      # 2020
                -1.2,    -1.8,     0.7,     0.1,     0.3,     0.8,      # 2020
                0.5,    -0.2,    -0.3,     0.6,    -0.5,     0.2,       # 2020
                0.1,    -0.9,    -1.1,     0.0,    -0.3,     0.5,       # 2020
                0.1,    -0.9,    -0.9,     0.5,     0.6,     0.7,       # 2020
                1.4,    -0.3,    -0.4,     0.8,    -1.3,     0.0,       # 2020
                -0.1,     0.8,     0.3,     0.0,    -0.1,     0.4,      # 2020
                0.5,     0.1,     0.5,     0.5,    -0.4,    -0.5,       # 2020
                -0.4,    -0.4,    -0.6,                                 # 2020
                5.7,     7.4,   -25.9,   -11.0,    -7.0,   -30.2,       # 2022
                -2.1,   -22.4,     2.2,    -5.9,     6.0,     3.1,      # 2022
                -1.1,   -12.0,     0.5,    -1.2,    -1.6,    -0.1,      # 2022
                -5.9,     6.5,     5.2,     3.6,    -5.1,    -5.0,      # 2022
                -0.3,     0.5,     0.0,    -0.6,     2.5,     0.2,      # 2022
                -0.6,     1.3,     3.0,     0.9,     0.3,    -0.5,      # 2022
                -0.3,     0.0,     0.4,    -1.6,     1.3,    -1.3,      # 2022
                -1.4,     0.8,     0.0,     0.0,     0.9,     1.0,      # 2022
                -0.1,    -0.2,     0.6,     0.0,     0.6,     0.7,      # 2022
                -0.8,     0.1,    -0.2,    -0.5,    -1.1,    -0.8,      # 2022
                0.1,     0.8,     0.3,     0.0,     0.1,    -0.2,       # 2022
                -0.1,     0.6,     0.4,    -0.2,    -0.1,     0.5,      # 2022
                0.4,    -0.3,     0.3,    -0.4,    -0.1,     0.5,       # 2022
                0.4,     0.0 ]) #                        # 2022
    gh = np.hstack([gh, np.zeros(115)])

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

    nc  = int(nmx * (nmx + 2))
    # nc = 195

    kmx = int((nmx+1) * (nmx + 2) / 2)
    # kmx = 105

    # allocate
    cl          = np.zeros(nmx)
    sl          = np.zeros(nmx)
    p           = np.zeros(kmx)
    q           = np.zeros(kmx)

    r  = alt
    ct          = np.cos(colat*np.pi/180)
    st          = np.sin(colat*np.pi/180)
    cl[0]       = np.cos(elong*np.pi/180)
    sl[0]       = np.sin(elong*np.pi/180)
    Cd = 1.0
    sd = 0.0
    l  = 1
    m  = 1
    n  = 0

    ratio = 6371.2/r
    rr    = ratio**2

    # Computation of Schmidt quasi-normal coefficients p and x(=q).
    # =============================================================

    p[0] = 1.0
    p[2] = st
    q[0] = 0.0
    q[2] = ct

    for k in range(2, kmx+1): 
        # There is no need to check bounds here. The code guarantees that
        # everything is inside the right bounds. This increased the source code
        # performance by 13%.
        if n < m:
            m  = 0
            n  = n+1
            rr = rr*ratio
            fn = n
            gn = n-1

        fm = m

        if (m == n):
            if k != 3:
                one   = np.sqrt(1 - 0.5/fm)
                j     = k - n - 2
                p[k-1]  = one*st*p[j]
                q[k-1]  = one*(st*q[j] + ct*p[j])
            # if ((k != 3)):
                cl[m-1] = cl[m-2]*cl[0] - sl[m-2]*sl[0]
                sl[m-1] = sl[m-2]*cl[0] + cl[m-2]*sl[0]
        else:
            gmm   = m**2
            one   = np.sqrt(fn**2 - gmm)
            two   = np.sqrt(gn**2 - gmm)/one
            three = (fn + gn)/one
            i     = k - n - 1
            j     = i - n + 1
            p[k-1]  = three*ct*p[i] - two*p[j]
            q[k-1]  = three*(ct*q[i] - st*p[i]) - two*q[j]

        # Synthesis of x, y, and z in geocentric coordinates.
        lm  = ll + l -1
        one = (tc*gh[lm] + t*gh[lm+nc])*rr

        if m != 0:
            two   = (tc*gh[lm+1] + t*gh[lm+nc+1])*rr
            three = one*cl[m-1] + two*sl[m-1]
            x     = x + three*q[k-1]
            z     = z - (fn + 1)*three*p[k-1]

            if st != 0:
                y = y + (one*sl[m-1] - two*cl[m-1])*fm*p[k-1]/st
            else:
                y = y + (one*sl[m-1] - two*cl[m-1])*q[k-1]*ct

            l = l + 2
        else:
            x = x + one*q[k-1]
            z = z - (fn + 1)*one*p[k-1]
            l = l + 1

        m = m + 1

    # Conversion to coordinate system specified by itype.
    # ===================================================
   
    one   = x 
    x     = x*Cd +   z*sd
    z     = z*Cd - one*sd

    return 1e6*np.array([x, y, z])

# UNVERIFIED 
def ecef_Q_ned_mat(lon, lat):
    sp, cp = np.sin(lat), np.cos(lat)
    sl, cl = np.sin(lon), np.cos(lon)

    return np.array([  [-sp*cl, -sl, -cp*cl],
                         [-sp*sl,  cl, -cp*sl],
                         [cp,     0.0, -sp]      ])

def geoc_from_ecef(r_ecef):
    # I dont think this is used because it seems to have an error in mag_field.jl...
    print("Using GEOC_from_ECEF")
    x, y, z = r_ecef 
    lat = np.atan2(z, np.sqrt(x*x + y*y))
    lon = np.atan2(y, x)
    return lat, lon


#######################
#  IGRF13_5 FUNCTIONS #
#######################
# From Kevin (For circuit python)

# clean the lists of arrays before using them
def reset_list(input_list):
    n = len(input_list[0])
    for i in range(len(input_list)):
        input_list[i] = np.zeros(n)

def reset_array(input_array):
    for i in range(len(input_array)):
        input_array[i] = 0.0

def igrf13_5(gh, date, latitude_degrees, elongitude_degrees, r_norm_km, cl, sl, p ,q):

    # reset the lists that are passed by reference
    reset_array(cl)
    reset_array(sl)
    reset_array(p)
    reset_array(q)

    # colatitude
    colat = 90 - latitude_degrees

    # Declaration of variables
    fn = 0
    gn = 0
    kmx = 0
    ll = 0
    nc = 0
    nmx = 0
    x = 0.0
    y = 0.0
    z = 0.0
    t = 0.0
    tc = 0.0

    t = date - 2020
    tc = 1.0

    ll = 0
    nmx = 5

    # nc = int(nmx * (nmx + 2))
    nc = 35

    # kmx = int((nmx + 1) * (nmx + 2) / 2)
    kmx = 21

    # allocate
    #cl = [0.0 for _ in range(nmx)]
    #sl = [0.0 for _ in range(nmx)]
    #p = [0.0 for _ in range(kmx)]
    #q = [0.0 for _ in range(kmx)]

    r = r_norm_km
    ct = np.cos(colat * np.pi / 180)
    st = np.sin(colat * np.pi / 180)
    cl[0] = np.cos(elongitude_degrees * np.pi / 180)
    sl[0] = np.sin(elongitude_degrees * np.pi / 180)
    Cd = 1.0
    sd = 0.0
    l = 1
    m = 1
    n = 0

    ratio = 6371.2 / r
    rr = ratio ** 2

    p[0] = 1.0
    p[2] = st
    q[0] = 0.0
    q[2] = ct

    for k in range(2, kmx + 1):
        if n < m:
            m = 0
            n = n + 1
            rr = rr * ratio
            fn = n
            gn = n - 1

        fm = m

        if m == n:
            if k != 3:
                one = np.sqrt(1 - 0.5 / fm)
                j = k - n - 1
                p[k - 1] = one * st * p[j - 1]
                q[k - 1] = one * (st * q[j - 1] + ct * p[j - 1])
                cl[m - 1] = cl[m - 2] * cl[0] - sl[m - 2] * sl[0]
                sl[m - 1] = sl[m - 2] * cl[0] + cl[m - 2] * sl[0]
        else:
            gmm = m ** 2
            one = np.sqrt(fn ** 2 - gmm)
            two = np.sqrt(gn ** 2 - gmm) / one
            three = (fn + gn) / one
            i = k - n
            j = i - n + 1
            p[k - 1] = three * ct * p[i - 1] - two * p[j - 1]
            q[k - 1] = three * (ct * q[i - 1] - st * p[i - 1]) - two * q[j - 1]

        lm = ll + l
        one = (tc * gh[lm - 1] + t * gh[lm + nc - 1]) * rr

        if m != 0:
            two = (tc * gh[lm] + t * gh[lm + nc]) * rr
            three = one * cl[m - 1] + two * sl[m - 1]
            x = x + three * q[k - 1]
            z = z - (fn + 1) * three * p[k - 1]

            if st != 0:
                y = y + (one * sl[m - 1] - two * cl[m - 1]) * fm * p[k - 1] / st
            else:
                y = y + (one * sl[m - 1] - two * cl[m - 1]) * q[k - 1] * ct

            l = l + 2

        else:
            x = x + one * q[k - 1]
            z = z - (fn + 1) * one * p[k - 1]
            l = l + 1

        m = m + 1

    one = x
    x = x * Cd + z * sd
    z = z * Cd - one * sd

    return np.array([x, y, z])

class igrfclass:
    gh = [-29404.8,     -1450.9,    4652.5,     -2499.6,    2982.0,      -2991.6,  # 2020
            1677.0,      -734.6,    1363.2,     -2381.2,     -82.1,       1236.2,  # 2020
             241.9,       525.7,    -543.4,       903.0,     809.5,        281.9,  # 2020
              86.3,      -158.4,    -309.4,       199.7,      48.0,       -349.7,  # 2020
            -234.3,       363.2,      47.7,       187.8,     208.3,       -140.7,  # 2020
            -121.2,      -151.2,      32.3,        13.5,      98.9,         66.0,  # 2020
              65.5,       -19.1,      72.9,        25.1,    -121.5,         52.8,  # 2020
             -36.2,       -64.5,      13.5,         8.9,     -64.7,         68.1,  # 2020
              80.6,       -76.7,     -51.5,        -8.2,     -16.9,         56.5,  # 2020
               2.2,        15.8,      23.5,         6.4,      -2.2,         -7.2,  # 2020
             -27.2,         9.8,      -1.8,        23.7,       9.7,          8.4,  # 2020
             -17.6,       -15.3,      -0.5,        12.8,       ]

    def __init__(self):
        self.cl = [0.0 for _ in range(5)]
        self.sl = [0.0 for _ in range(5)]
        self.p = [0.0 for _ in range(21)]
        self.q = [0.0 for _ in range(21)]

    def ned_igrf(self,date, latitude_degrees, elongitude_degrees, r_norm_km):
        return (igrf13_5(self.gh, date, latitude_degrees, elongitude_degrees, r_norm_km,self.cl,self.sl,self.p,self.q))


# igrf = igrfclass()

# R_EARTH = 6.378136300e6  # m
# date = 2020.3
# latitude_degrees = 34.567
# elongitude_degrees = 45.678
# r_norm_km = 1.05 * R_EARTH / 1000

# a = igrf.ned_igrf(date, latitude_degrees, elongitude_degrees, r_norm_km)








######################
#  Orbit_dynamics.jl #
######################

# NOTE that this does not have the sign fixed so it will match Julias eclipse_conical(x, r_sun)
print("Need to fix eclipse_conical sign in python AND julia")
def eclipse_conical(x, r_sun):
    R_SUN = 6.957e8
    R_EARTH = 6.3781363e6

    r = x[:3]  # Inertial position of satellite 

    # Occultation Geometry 
    a = np.arcsin(R_SUN / np.linalg.norm(r_sun - r))
    b = np.arcsin(R_EARTH / np.linalg.norm(r))
    c = np.arccos(np.dot( r, r_sun - r) / (np.linalg.norm(r) * np.linalg.norm(r_sun - r)) )

    e_sun = r_sun / np.linalg.norm(r_sun)

    # Test Occulatation Conditions 
    nu = 0.0 
    if (np.abs(a - b) < c) and (c < (a + b)):
        # Partial Occultation 
        xx = ((c**2) + (a**2) - (b**2)) / (2*c)
        yy = np.sqrt( (a**2) - (xx**2) )
        A  = (a**2) * np.arccos(xx/a) + (b**2) * np.arccos( (c-xx)/b) - (c*yy)

        nu = 1 - A/(np.pi * (a**2))

    elif (a + b) <= c:
        # No Occultation 
        nu = 1.0 

    else:
        # Full Occultation 
        nu = 0.0

    return nu








