import numpy as np
from datetime import datetime, timezone
from astropy.time import Time
from astropy.coordinates import get_sun, EarthLocation, AltAz
#from astropy.coordinates import get_sun
import astropy.units as u
from astropy.coordinates import solar_system_ephemeris, EarthLocation
from astropy.coordinates import get_body #, get_body_barycentric
import matplotlib.pyplot as plt
#import matplotlib.animation as animation
from mpl_interactions import ioff, panhandler, zoom_factory
import matplotlib.pyplot as plt
import time

fpath = "FILE_PATH\\scan_log.txt"

fpath2 = "FILE_PATH\\scan_matrix.txt"

# !!! constant coordinate format (alt,az) with np arrays !!!

# initializing tables
sun_coords_ephimeris_previous = np.array([0,0])
sun_coords_ephimeris_now = np.array([0,0])

# defining basic parameters

# degrees (for now, we will keep the same step for alt and az)
d_az = 0.1
d_alt = 0.1

# degrees
alt_dim = 5
az_dim = 5

alt_count = int(alt_dim/d_alt) # grid rows count
az_count = int(az_dim/d_az) # grid columns count

# telescope earth coordinates (the following are obviously not correct and only for experimental purposes)
telescopeLongitude = 22.9638
telescopeLatitude = 40.6401
elevation = 0 * u.m

# matrix that will contain all the data [(alt_sun_ephimeris,az_sun_ephimeris),(alt_sun_history,az_sun_history),signal_strength]
matrix = np.empty([alt_count,az_count], dtype=object) # will be altCount x azCount (10x10)
# set initial values to [-999.9999,-999.9999,-999.9999,-999.9999]
four_minus_999 = np.array([-999.9999,-999.9999,-999.9999,-999.9999])

for i in range(alt_count):
    for j in range(az_count):
        matrix[i,j] = four_minus_999
# this way we will know if a value has not been changed

def getImuCoords():
    #rng = np.random.default_rng()

    #rand_az = rng.uniform(0, 360)
    #rand_alt = rng.uniform(-180, 180)

    return np.array([60, 53])#np.array([rand_alt, rand_az])

def getSunCoordsEphimeris(time):
    # uses astropy to get the coordinates, based on the sun_tracking programm

    # time = Time(datetime.now(timezone.utc))

    location = EarthLocation(lat=telescopeLatitude * u.deg,
                             lon=telescopeLongitude * u.deg,
                             height=elevation)

    with solar_system_ephemeris.set('jpl'):
        sun_coord = get_body('sun', time, location)
    
    #------------from JPL ephimerides--------------------------

    ##print("--------from astropy jpl ephimerides--------")

    altaz_frame = AltAz(obstime=time, location=location)
    sun_altaz = sun_coord.transform_to(altaz_frame)

    alt_deg = sun_altaz.alt.degree
    az_deg = sun_altaz.az.degree

    ##print(f"sun's altitude in degrees: {alt_deg}")
    ##print(f"sun's azimuth in degrees: {az_deg}")
    ##print(f"time: {time}")

    return np.array([alt_deg,az_deg])

def moveRtAltAz(d_alt_az):
    # will move the steppers
    # moveAlt(dAltAz(0))
    # moveAz(dAltAz(1))

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
    #print(f"initial_coords={initial_coords}")

    theoretical_coords = initial_coords # theoretical coords of the RT

    d_coords_initial = getSunCoordsEphimeris(Time(datetime.now(timezone.utc))) - initial_coords - np.array([az_dim/2,alt_dim/2])

    theoretical_coords = theoretical_coords + d_coords_initial

    moveRtAltAz(d_coords_initial) # moves the RT to an initial approximate position, the scan will have this position as center

    multiplier = -1

    # starting scan
    # SPEED = MAXIMUM

    f=open(fpath,"w")
    f.write("-----------------START OF LOG-----------------\n")

    prev_sun_alt_az = np.array([0.0,0.0])

    X_DATA = np.empty((0))
    Y_DATA = np.empty((0))

    #print(X_DATA)
    #print(Y_DATA)

    #print(theoretical_coords)
    print("Scanning, this will take a while...")

    for i in range(alt_count):
        multiplier = -multiplier

        current_sun_coords = getSunCoordsEphimeris(Time(datetime.now(timezone.utc)))

        for j in range(az_count):

            if i==0 and j==0:
                d_coords_total = np.array([0.0,0.0])
            elif j==0 and i != 0:
                current_sun_coords = getSunCoordsEphimeris(Time(datetime.now(timezone.utc)))
                d_coords_total = current_sun_coords - prev_sun_alt_az
            else:
                current_sun_coords = getSunCoordsEphimeris(Time(datetime.now(timezone.utc)))

                #d_alt_az_movement = current_sun_coords - prev_sun_alt_az

                #d_alt_az_scan = np.array([0,multiplier*d_az])

                d_coords_total = current_sun_coords - prev_sun_alt_az + np.array([0.0,multiplier*d_az])

            # Each time the coordinates will be incremented by d_coords_movement and by d_coords_scan

            moveRtAltAz(d_coords_total)

            theoretical_coords = theoretical_coords + d_coords_total # based on the coords of the initial point, that was found with the IMU

            signal_strength = getSdrSignalStrength()

            t = Time(datetime.now(timezone.utc))

            t_matrix = np.array([t], dtype=object)

            matrix_to_append = np.array([np.array([theoretical_coords[0]]),np.array([theoretical_coords[1]]),signal_strength,t_matrix])

            f.write(f"matrix_to_append={matrix_to_append}")

            f.write("\n")

            matrix[i,j] = matrix_to_append

            X_DATA = np.append(X_DATA,np.array([theoretical_coords[1]]),axis=0)
            Y_DATA = np.append(Y_DATA,np.array([theoretical_coords[0]]),axis=0)

            prev_sun_alt_az = current_sun_coords

            time.sleep(0.1)
        
        # change alt
        if i==alt_count-1:
            print("end")
        else:
            moveRtAltAz(np.array([d_alt,0]))

            theoretical_coords[0] = theoretical_coords[0] + d_alt # based on the coords of the initial point, that was found with the IMU

            time.sleep(0.1)
    
    print("Scan is over")

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

    f.write(f"max_signal={max_signal}, i_max={i_max}, jmax={j_max}  \n")

    # finding the corrections

    alt_max = matrix[i_max,j_max][0][0]
    az_max = matrix[i_max,j_max][1][0]
    t_max = matrix[i_max,j_max][3][0]

    alt_az_max = np.array([alt_max,az_max])
    alt_az_max_ephimeris = getSunCoordsEphimeris(t_max)

    corrections = alt_az_max_ephimeris - alt_az_max

    print(f"corrections={corrections}")
    f.write(f"corrections={corrections} \n")

    f.write("-----------------END OF LOG-----------------\n")
    f.close()
    # will be 1 x 2 matrix

    # so whenever we will need to find a point we will use the initial IMU coord, but corrected with the corrections
    # keeping track of the movements of the radiotelescope, we will know where it points

    f2 = open(fpath2,"w")
    for row in range(len(matrix)):
        for column in range(len(matrix[0])):
            f2.write(str(matrix[row,column][0]) + str(matrix[row,column][1]) + str(matrix[row,column][2]) + str(matrix[row,column][3]) + "\n")
    
    f2.close()

    num_points = alt_count*az_count

    
    '''
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            X_DATA = np.append(X_DATA,matrix[i,j][0],axis=0)
            Y_DATA = np.append(Y_DATA,matrix[i,j][1],axis=0)
    '''
    #print("xdata")
    #print(X_DATA)
    #print(Y_DATA)

    fig, ax = plt.subplots()

    #ax.scatter(X_DATA, Y_DATA)
    line, = ax.plot(X_DATA, Y_DATA, marker="o", picker=5)

    ax.scatter(initial_coords[1],initial_coords[0], color="red", s=100, marker ="x")
    zz = initial_coords + d_coords_initial
    ax.scatter(zz[1],zz[0], color="green", s=100, marker ="x")

    ax.grid()
    ax.set_aspect("equal", adjustable='box')
    ax.set_xlim(min(X_DATA), max(X_DATA))
    ax.set_ylim(min(Y_DATA), max(Y_DATA))
    disconnect_zoom = zoom_factory(ax)
    # Enable scrolling and panning with the help of MPL
    # Interactions library function like panhandler.
    pan_handler = panhandler(fig)
    #for xy in zip(X_DATA, Y_DATA):                                       # <--
    #    ax.annotate('(%.6f, %.6f)' % xy, xy=xy, textcoords='data',fontsize=8) # <--
    # Create an initial annotation object, hidden by default
    

    annot = ax.annotate("", xy=(0., 0.), xytext=(20, 20),
                    textcoords="offset points",
                    bbox=dict(boxstyle="round,pad=0.5", fc="yellow", alpha=0.5),
                    arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=0.2"))
    annot.set_visible(False)

    # 3. Define the callback function to handle the click event
    def on_pick(event):
        if event.artist != line:
            return

        # Get the index of the closest point
        ind = event.ind[0]
        thisx = event.artist.get_xdata()[ind]
        thisy = event.artist.get_ydata()[ind]

        # Update the annotation position (xy) and text
        annot.xy = (thisx, thisy)
        text = f"x={thisx:.6f}, y={thisy:.6f}"
        annot.set_text(text)
        
        # Make the annotation visible and redraw the canvas
        annot.set_visible(True)
        fig.canvas.draw_idle() # Redraw the figure without blocking the GUI loop

    # 4. Connect the event handler
    fig.canvas.mpl_connect('pick_event', on_pick)

    # 5. Add a second event to hide the annotation on a general mouse click if desired
    def on_click_anywhere(event):
        # Hide annotation if the user clicks anywhere else in the plot area
        if event.artist != line and annot.get_visible():
            annot.set_visible(False)
            fig.canvas.draw_idle()

    fig.canvas.mpl_connect('button_press_event', on_click_anywhere)

    plt.show()

    # Assume data is shape (M, N) and each element is a 1D array of length 4
    M, N = matrix.shape

    x = np.zeros(M * N)
    y = np.zeros(M * N)
    z = np.zeros(M * N)
    t = np.empty(M * N, dtype=object)   # Time objects

    k = 0
    for i in range(M):
        for j in range(N):
            element = matrix[i, j]   # shape (4,)
            x[k] = element[0]
            y[k] = element[1]
            z[k] = element[2]
            t[k] = element[3]      # Time object
            k += 1
    plt.scatter(x, y, c=z)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Scatter plot with z as color dimension")
    plt.colorbar(label="z")
    plt.show()

    #print(X_DATA[0:10])
    #print(Y_DATA[0:10])

    #print(np.where(X_DATA < 1)[0])

if __name__ == "__main__":
    main()
