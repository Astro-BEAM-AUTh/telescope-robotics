import numpy as np
from datetime import datetime, timezone
from astropy.time import Time
from astropy.coordinates import get_sun, EarthLocation, AltAz
from astropy.coordinates import get_sun
import astropy.units as u
from astropy.coordinates import solar_system_ephemeris, EarthLocation
from astropy.coordinates import get_body_barycentric, get_body

import turtle

# !!! constant coordinate format (alt,az) with np arrays !!!

# initializing tables
sun_coords_ephimeris_previous = np.array([0,0])
sun_coords_ephimeris_now = np.array([0,0])

# defining basic parameters

# degrees (for now, we will keep the same step for alt and az)
d_az = 10
d_alt = 10

alt_dim = 100
az_dim = 100

alt_count = int(alt_dim/d_alt) # grid rows count
az_count = int(az_dim/d_az) # grid columns count

# telescope earth coordinates (the following are obviously not correct and only for experimental purposes)
telescopeLongitude = 22.9638
telescopeLatitude = 40.6401
elevation = 0 * u.m

# turtle color
turtle.color("green")
screen = turtle.Screen()
screen.screensize(400, 400)
turtle.speed(1)

# matrix that will contain all the data [(alt_sun_ephimeris,az_sun_ephimeris),(alt_sun_history,az_sun_history),signal_strength]
matrix = np.empty([alt_count,az_count], dtype=object) # will be altCount x azCount (10x10)
# set initial values to [-999.9999,-999.9999,-999.9999,-999.9999]
four_minus_999 = np.array([-999.9999,-999.9999,-999.9999,-999.9999])

for i in range(alt_count):
    for j in range(az_count):
        matrix[i,j] = four_minus_999
# this way we will know if a value has not been changed

def getImuCoords():
    rng = np.random.default_rng()

    rand_az = rng.uniform(0, 360)
    rand_alt = rng.uniform(-180, 180)

    return np.array([rand_alt, rand_az])

def getSunCoordsEphimeris(time):
    # uses astropy to get the coordinates, based on the sun_tracking programm

    # time = Time(datetime.now(timezone.utc))

    location = EarthLocation(lat=telescopeLatitude * u.deg,
                             lon=telescopeLongitude * u.deg,
                             height=elevation)

    with solar_system_ephemeris.set('jpl'):
        sun_coord = get_body('sun', time, location)
    
    #------------from JPL ephimerides--------------------------

    print("--------from astropy jpl ephimerides--------")

    altaz_frame = AltAz(obstime=time, location=location)
    sun_altaz = sun_coord.transform_to(altaz_frame)

    alt_deg = sun_altaz.alt.degree
    az_deg = sun_altaz.az.degree

    print(f"sun's altitude in degrees: {alt_deg}")
    print(f"sun's azimuth in degrees: {az_deg}")
    print(f"time: {time}")

    return np.array([alt_deg,az_deg])

def moveRtAltAz(d_alt_az):
    # will move the steppers
    # moveAlt(dAltAz(0))
    # moveAz(dAltAz(1))
    currentX,currentY = turtle.position()

    turtle.goto(currentX+d_alt_az[1],currentY+d_alt_az[0])

    return

def getSdrSignalStrength():
    # ...
    rng = np.random.default_rng()

    rand_signal_strength = rng.uniform(0, 100)
    return np.array([rand_signal_strength])

def main():
    # get initial approximate coordinates from IMU and point radiotelescope based on them
    # SPEED = MINIMUM

    initial_coords = getImuCoords()
    print(f"initial_coords={initial_coords}")

    theoretical_coords = initial_coords

    d_coords_initial = getSunCoordsEphimeris(Time(datetime.now(timezone.utc))) - initial_coords

    moveRtAltAz(d_coords_initial) # moves the RT to an initial approximate position, the scan will have this position as center

    multiplier = -1

    # starting scan
    # SPEED = MAXIMUM

    print("-----------------START OF LOG-----------------\n")

    turtle.color("blue")

    moveRtAltAz(np.array([-az_dim/2,-alt_dim/2]))

    turtle.color("red")

    for i in range(alt_count):
        multiplier = -multiplier

        for j in range(az_count):

            # Each time the coordinates will be incremented by d_coords_movement and by d_coords_scan

            moveRtAltAz(np.array([0,multiplier*d_az]))

            theoretical_coords[1] = theoretical_coords[1] + multiplier*d_az # based on the coords of the initial point, that was found with the IMU

            signal_strength = getSdrSignalStrength()

            t = Time(datetime.now(timezone.utc))

            t_matrix = np.array([t], dtype=object)

            matrix_to_append = np.array([np.array([theoretical_coords[0]]),np.array([theoretical_coords[1]]),signal_strength,t_matrix])

            print(f"matrix_to_append={matrix_to_append}")

            print("\n")

            matrix[i,j] = matrix_to_append
        
        # change alt
        if i==alt_count-1:
            print("end")
        else:
            moveRtAltAz(np.array([d_alt,0]))

            theoretical_coords[0] = theoretical_coords[0] + d_alt # based on the coords of the initial point, that was found with the IMU

    # finding maximum signal strength

    i_max = -1
    j_max = -1
    max_signal = -999
    element_ij2 = 0

    for i in range(alt_count):
        for j in range(az_count):
            element_ij2 = matrix[i,j][2]
            if element_ij2 > max_signal:
                max_signal = element_ij2
                i_max = i
                j_max = j

    print(f"max_signal={max_signal}, i_max={i_max}, jmax={j_max}")

    print(f"max_signal={max_signal}, i_max={i_max}, jmax={j_max}  \n")

    # finding the corrections

    alt_max = matrix[i_max,j_max][0][0]
    az_max = matrix[i_max,j_max][1][0]
    t_max = matrix[i_max,j_max][3][0]

    # the theoretical coordinates of the sun, when maximum signal strength was found
    alt_az_max = np.array([alt_max,az_max])

    # the real coordinates of the sun at that same time, from the ephimerides
    alt_az_max_ephimeris = getSunCoordsEphimeris(t_max)

    corrections = alt_az_max_ephimeris - alt_az_max

    print(f"corrections={corrections}")

    print("-----------------END OF LOG-----------------\n")
    # will be 1 x 2 matrix

    # so whenever we will need to find a point we will use the initial IMU coord, but corrected with the corrections
    # keeping track of the movements of the radiotelescope, we will know where it points

    input("type sth to exit")

if __name__ == "__main__":

    main()
