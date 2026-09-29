# For the motors
from time import sleep
import RPi.GPIO as gpio
# For the GUI
import tkinter as tk
from tkinter import ttk, messagebox
# threading
import threading

stop_flag = False

# Setting up drivers and motors

# motor
direction_pin1 = 21
pulse_pin1 = 20

# directions:
# clockwise (cw) -> LOW
# counterclockwise (ccw) -> HIGH
cw_direction = gpio.LOW
ccw_direction = gpio.HIGH

# setting GPIO pin numbering system -> BCM (Broadcom SOC channel)
gpio.setmode(gpio.BCM)

# setting pins for the motor
gpio.setup(direction_pin1, gpio.OUT)
gpio.setup(pulse_pin1, gpio.OUT)


def degreesToSteps(angle, steps_per_rev):
    steps = float(20.0 * (angle / 360.0) * steps_per_rev) # multiplying by 20, because we have a drive systme with 20:1 ratio
    return int(steps)

def step_motor(steps_per_rev, dt, direction, angle):
    global stop_flag

    try:
        if direction == "CW":
            gpio.output(direction_pin1, cw_direction)
        elif direction == "CCW":
            gpio.output(direction_pin1, ccw_direction)
        else:
            messagebox.showerror("Error", "Invalid direction")
            return

        maxSteps = degreesToSteps(angle, steps_per_rev)

        for currentStep in range(maxSteps):
            if stop_flag:
                print("Motor stopped")
                break

            gpio.output(pulse_pin1, gpio.HIGH)
            sleep(dt)
            gpio.output(pulse_pin1, gpio.LOW)

    except KeyboardInterrupt:
        gpio.cleanup()


def run_motor():
    global stop_flag
    stop_flag = False

    try:
        steps_per_rev = int(steps_entry.get())
        dt = float(dt_entry.get())
        angle = float(angle_entry.get())
        direction = direction_var.get()

        thread = threading.Thread(
            target=step_motor,
            args=(steps_per_rev, dt, direction, angle)
        )
        thread.start()

    except ValueError:
        messagebox.showerror("Input Error", "Please enter valid numeric values.")

def stop_motor():
    global stop_flag
    stop_flag = True

# GUI

# Main window
root = tk.Tk()
root.title("Stepper Motor Control")
root.geometry("450x400")

# Style font
FONT = ("Arial", 15)

# Steps per revolution
ttk.Label(root, text="Steps per revolution:", font=FONT).pack(pady=(15, 0))
steps_entry = ttk.Entry(root, font=FONT)
steps_entry.insert(0, "6400")
steps_entry.pack()

# dt (speed)
ttk.Label(root, text="dt (speed):", font=FONT).pack(pady=(15, 0))
dt_entry = ttk.Entry(root, font=FONT)
dt_entry.insert(0, "0.001")
dt_entry.pack()

# Angle
ttk.Label(root, text="Angle to step:", font=FONT).pack(pady=(15, 0))
angle_entry = ttk.Entry(root, font=FONT)
angle_entry.pack()

# Direction radio buttons
direction_var = tk.StringVar(value="CW")

ttk.Label(root, text="Direction:", font=FONT).pack(pady=(15, 0))

frame = ttk.Frame(root)
frame.pack()

ttk.Radiobutton(frame, text="Clockwise", variable=direction_var, value="CW").pack(side="left", padx=10)
ttk.Radiobutton(frame, text="Counter Clockwise", variable=direction_var, value="CCW").pack(side="left", padx=10)

button_frame = ttk.Frame(root)
button_frame.pack(pady=25)

ttk.Button(button_frame, text="Run Motor", command=run_motor).pack(side="left", padx=10)
ttk.Button(button_frame, text="STOP", command=stop_motor).pack(side="left", padx=10)

# Run app
root.mainloop()