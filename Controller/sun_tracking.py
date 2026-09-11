# motors and sensor
import RPi.GPIO as GPIO
import time
import board
import adafruit_bno055

# skyfield and plots
from skyfield.api import load, wgs84
import numpy as np

STEPS_PER_DEGREE_ALT = 17.7
STEPS_PER_DEGREE_AZ = 17.7

# --- GPIO SETUP ---
DIR1, PUL1 = 21, 20
DIR2, PUL2 = 23, 22

GPIO.setmode(GPIO.BCM)
for pin in [DIR1, PUL1, DIR2, PUL2]:
    GPIO.setup(pin, GPIO.OUT)
    GPIO.output(pin, GPIO.LOW)
    
# setting motor speed
dt = 0.005 # intermediate speed (for high speed set dt = 0.001)

# BNO055 sensor setup
class Mode:
    CONFIG_MODE = 0x00
    ACCONLY_MODE = 0x01
    MAGONLY_MODE = 0x02
    GYRONLY_MODE = 0x03
    ACCMAG_MODE = 0x04
    ACCGYRO_MODE = 0x05
    MAGGYRO_MODE = 0x06
    AMG_MODE = 0x07
    IMUPLUS_MODE = 0x08
    COMPASS_MODE = 0x09
    M4G_MODE = 0x0A
    NDOF_FMC_OFF_MODE = 0x0B
    NDOF_MODE = 0x0C

i2c = board.I2C()  # uses board.SCL and board.SDA
# i2c = board.STEMMA_I2C()  # For using the built-in STEMMA QT connector on a microcontroller
sensor = adafruit_bno055.BNO055_I2C(i2c)

# setting sensor mode
# sensor.mode = Mode.ACCMAG_MODE #accelerometer and magnetometer
# sensor.mode = Mode.NDOF_MODE  # Set the sensor to NDOF_MODE # not needed, it is default
# sensor.mode = Mode.AMG_MODE #accelerometer, magnetometer and gyroscope
# If you are going to use UART uncomment these lines
# uart = board.UART()
# sensor = adafruit_bno055.BNO055_UART(uart)

last_val = 0xFFFF

magnetic_declination = 5.31 # degrees

# Tracking setup
az_list = []
alt_list = []

UPDATE_INTERVAL = 1.0  # slower updates are fine

# coordinates
# (Wellington - New Zealand)
LAT = 41.2924
LON = 174.7787
ELEV = 50

ts = load.timescale()

# Load planetary data (includes Sun)
eph = load('de421.bsp')
sun = eph['sun']
earth = eph['earth']

observer = earth + wgs84.latlon(LAT, LON, elevation_m=ELEV)

running = True

def get_dish_dir_sensor():
    euler_angle = sensor.euler # format: (heading, roll, pitch) (where heading is azimuth)
    print(f"Euler angle: {euler_angle}")
    #current_alt = 0.0
    #current_az = 0.0
    current_alt = euler_angle[2]
    current_az = euler_angle[0] + magnetic_declination
    return current_alt, current_az

def get_sensor_temperature():
    global last_val  # noqa: PLW0603
    result = sensor.temperature
    if abs(result - last_val) == 128:
        result = sensor.temperature
        if abs(result - last_val) == 128:
            return 0b00111111 & result
    last_val = result
    return result

def motor_control(alt_steps, az_steps, dt):
    print(f"Moving: Alt={alt_steps}, Az={az_steps}, dt={dt}")
    
    # Set Motor Directions
    GPIO.output(DIR1, GPIO.HIGH if alt_steps >= 0 else GPIO.LOW)
    GPIO.output(DIR2, GPIO.HIGH if az_steps >= 0 else GPIO.LOW)
    
    # Take the absolute value for the loop count
    alt_remaining = abs(alt_steps)
    az_remaining = abs(az_steps)
    max_steps = max(alt_remaining, az_remaining)

    for _ in range(max_steps):        
        if alt_remaining > 0:
            GPIO.output(PUL1, GPIO.HIGH)
        if az_remaining > 0:
            GPIO.output(PUL2, GPIO.HIGH)

        time.sleep(dt)
        
        GPIO.output(PUL1, GPIO.LOW)
        GPIO.output(PUL2, GPIO.LOW)
        
        time.sleep(dt)
        
        alt_remaining -= 1
        az_remaining -= 1

def on_key(event):
    global running
    if event.key == 'q':
        print("Q pressed → stopping", flush=True)
        running = False

def get_sun_altaz_at(t):
    astrometric = observer.at(t).observe(sun)
    apparent = astrometric.apparent()
    
    alt, az, distance = apparent.altaz()

    print("alt = ")
    print(alt.degrees)
    print("az = ")
    print(az.degrees)
    
    return alt.degrees, az.degrees

# prev_alt = None
# prev_az = None
prev_alt, prev_az = get_dish_dir_sensor()

def get_sun_state():
    global prev_alt, prev_az

    alt, az = get_sun_altaz_at(ts.now())

    if prev_alt is None:
        prev_alt = alt
        prev_az = az
        return alt, az, 0.0, 0.0

    d_alt = alt - prev_alt
    d_az  = az - prev_az

    if d_az > 180:
        d_az -= 360
    elif d_az < -180:
        d_az += 360

    prev_alt = alt
    prev_az = az

    return alt, az, d_alt, d_az

def track_sun():
    global running

    while running:

        dir_from_sensor = get_dish_dir_sensor()
        alt_sensor = dir_from_sensor[0]
        if alt_sensor<0 or alt_sensor>180:
            break
    
        alt, az, d_alt, d_az = get_sun_state()

        d_alt_steps = int(d_alt * STEPS_PER_DEGREE_ALT)
        d_az_steps = int(d_az * STEPS_PER_DEGREE_AZ)

        motor_control(d_alt_steps, d_az_steps)

    print("Stopping tracking...", flush=True)

if __name__ == "__main__":
    track_sun()