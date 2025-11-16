# SPDX-FileCopyrightText: 2021 ladyada for Adafruit Industries
# SPDX-License-Identifier: MIT

import time

import board

import adafruit_bno055

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


def temperature():
    global last_val  # noqa: PLW0603
    result = sensor.temperature
    if abs(result - last_val) == 128:
        result = sensor.temperature
        if abs(result - last_val) == 128:
            return 0b00111111 & result
    last_val = result
    return result


while True:
    #print(f"Temperature: {sensor.temperature} degrees C")
    """
    print(
        "Temperature: {} degrees C".format(temperature())
    )  # Uncomment if using a Raspberry Pi
    """
    #print(f"Accelerometer (m/s^2): {sensor.acceleration}")
    #print(f"Magnetometer (microteslas): {sensor.magnetic}")
    #print(f"Gyroscope (rad/sec): {sensor.gyro}")
    print(f"Euler angle: {sensor.euler}")
    #print(f"Quaternion: {sensor.quaternion}")
    #print(f"Linear acceleration (m/s^2): {sensor.linear_acceleration}")
    #print(f"Gravity (m/s^2): {sensor.gravity}")
    print()

    time.sleep(1)
