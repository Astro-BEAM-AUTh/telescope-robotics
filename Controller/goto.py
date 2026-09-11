import RPi.GPIO as GPIO
import time
import board
import adafruit_bno055

# --- CONFIGURATION / CONSTANTS ---
STEPS_PER_DEGREE_ALT = 17.7
STEPS_PER_DEGREE_AZ = 17.7

# GPIO Pin Mapping (BCM numbering)
DIR1, PUL1 = 21, 20  # Altitude Motor
DIR2, PUL2 = 23, 22  # Azimuth Motor

# --- INITIALIZATION ---
GPIO.setmode(GPIO.BCM)
for pin in [DIR1, PUL1, DIR2, PUL2]:
    GPIO.setup(pin, GPIO.OUT)
    GPIO.output(pin, GPIO.LOW)

print("Telescope GOTO Control System Initialized.")

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

    print(f"Moving {alt_steps} Alt steps, {az_steps} Az steps...")
    print("!! Press 'q' then 'Enter' to Emergency Stop !!")

    for _ in range(max_steps):
        # Check if the radiotelescope's altitude is still inside limits.
        dir_from_sensor = get_dish_dir_sensor()
        alt_sensor = dir_from_sensor[0]
        if alt_sensor<0 or alt_sensor>180:
            break
        
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

def goto():
    # Validate Target
    # 1. Get initial position from sensor
    curr_alt, curr_az = get_dish_dir_sensor()
    print("-" * 40)
    print(f"SENSOR READOUT: Alt={curr_alt:.2f}°, Az={curr_az:.2f}°")
    
    try:
        t_alt = float(input("Enter Target Altitude (0 to 180): "))
        t_az = float(input("Enter Target Azimuth: "))
        t_dt = input("Enter dt (Speed) [Default 0.001]: ")
        dt = float(t_dt) if t_dt.strip() else 0.001
    except ValueError:
        print("Invalid input. Please enter numeric values.")
        return

    # --- CONSTRAINT 1: Altitude [0, 180] ---
    if not (0 <= t_alt <= 180):
        print(f"LIMIT REJECTED: Target Alt {t_alt}° is outside allowed range [0, 180].")
        return

    # --- CONSTRAINT 2: Azimuth Delta within +-360 ---
    az_diff_deg = t_az - curr_az
    if abs(az_diff_deg) > 360:
        print(f"LIMIT REJECTED: Target move of {az_diff_deg:.2f}° exceeds 360° safety limit.")
        print("Reset your azimuth base or check your coordinates.")
        return

    # 2. Convert angle differences to motor pulses
    alt_diff_deg = t_alt - curr_alt
    d_alt = int(alt_diff_deg * STEPS_PER_DEGREE_ALT)
    d_az = int(az_diff_deg * STEPS_PER_DEGREE_AZ)

    # 3. Execute
    print("---------------------------------")
    print("d_alt = ")
    print(d_alt)
    print("d_az = ")
    print(d_az)
    print("---------------------------------")
    success = motor_control(d_alt, d_az, dt)
    
    if success:
        print("Target Reached Successfully.")
    else:
        print("Movement Interrupted.")

def loop_goto():
    try:
        while True:
            goto()
            choice = input("\nPerform another move? (y/n): ")
            if choice.lower() != 'y':
                break
    finally:
        GPIO.cleanup()
        print("\nGPIO Cleaned up. Program Terminated.")

if __name__ == "__main__":
    loop_goto()