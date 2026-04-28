import curses
import RPi.GPIO as GPIO
import time

def main_control(stdscr, up_down_step, left_right_step, dt):
    # Set getch() to non-blocking mode (waits 100ms then moves on)
    curses.halfdelay(1) 
    stdscr.keypad(True)
    curses.noecho()
    curses.curs_set(0) # Hide the blinking cursor

    while True:
        try:
            key = stdscr.getch()
            
            # Clear lines to refresh status
            stdscr.move(0, 0)
            stdscr.clrtoeol() 

            if key == ord('q'):
                # Hard low on pulses for safety
                GPIO.output(PUL1, GPIO.LOW)
                GPIO.output(PUL2, GPIO.LOW)
                break
            elif key == curses.KEY_UP:
                stdscr.addstr(0, 0, f"Status: MOVING UP (Step: {up_down_step})")
                motor_worker(up_down_step,0,dt)
            elif key == curses.KEY_DOWN:
                stdscr.addstr(0, 0, f"Status: MOVING DOWN (Step: {up_down_step})")
                motor_worker(-up_down_step,0,dt)
            elif key == curses.KEY_LEFT:
                stdscr.addstr(0, 0, f"Status: MOVING LEFT (Step: {left_right_step})")
                motor_worker(0,left_right_step,dt)
            elif key == curses.KEY_RIGHT:
                stdscr.addstr(0, 0, f"Status: MOVING RIGHT (Step: {left_right_step})")
                motor_worker(0,-left_right_step,dt)
            else:
                stdscr.addstr(0, 0, "Status: IDLE - Ready for Motor Input")
                GPIO.output(PUL1, GPIO.LOW)
                GPIO.output(PUL2, GPIO.LOW)

            stdscr.addstr(2, 0, f"Settings: dt={dt} | Step Up/Down={up_down_step} | Step Left/Right={left_right_step}")
            stdscr.addstr(3, 0, "Press 'q' to quit.")
            stdscr.refresh()
            
        except KeyboardInterrupt:
            break

def motor_worker(ud_steps, lr_steps, dt):
    """Background thread to pulse motors."""
    print(f"Moving: UD={ud_steps}, LR={lr_steps}, dt={dt}")
    
    # Set Directions (Assuming positive is one way, negative the other)
    GPIO.output(DIR1, GPIO.HIGH if ud_steps >= 0 else GPIO.LOW)
    GPIO.output(DIR2, GPIO.HIGH if lr_steps >= 0 else GPIO.LOW)
    
    # Take the absolute value for the loop count
    ud_remaining = abs(ud_steps)
    lr_remaining = abs(lr_steps)
    max_steps = max(ud_remaining, lr_remaining)

    for _ in range(max_steps):        
        # Pulse Motor 1 if steps remain
        if ud_remaining > 0:
            GPIO.output(PUL1, GPIO.HIGH)
        # Pulse Motor 2 if steps remain
        if lr_remaining > 0:
            GPIO.output(PUL2, GPIO.HIGH)

        time.sleep(dt)
        
        GPIO.output(PUL1, GPIO.LOW)
        GPIO.output(PUL2, GPIO.LOW)
        
        time.sleep(dt)
        
        ud_remaining -= 1
        lr_remaining -= 1

# --- GPIO SETUP ---
DIR1, PUL1 = 21, 20
DIR2, PUL2 = 23, 22

GPIO.setmode(GPIO.BCM)
for pin in [DIR1, PUL1, DIR2, PUL2]:
    GPIO.setup(pin, GPIO.OUT)

# --- DATA COLLECTION (Outside Curses) ---
print("""Enter movement step for each axis, in motor microsteps. 
The value must be an integer. 
Default value: 10 for both axes.
      """)

u_in = input("up/down step [10]: ")
l_in = input("left/right step [10]: ")

print("""Enter dt in seconds -> speed
(dt is the time between two motor microsteps -> the lower the value of dt, the hogher the speed)
Default value: 0.001
      """)

d_in = input("dt seconds [0.001]: ")

# Simple logic to handle empty inputs (defaults)
up_down_step = int(u_in) if u_in.strip() else 10
left_right_step = int(l_in) if l_in.strip() else 10
dt = float(d_in) if d_in.strip() else 0.001

# --- START CURSES ---
# We pass the variables into the function
curses.wrapper(main_control, up_down_step, left_right_step, dt)