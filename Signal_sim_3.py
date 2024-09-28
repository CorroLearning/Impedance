import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import Button

# Constants for the impedance values
Rel = 1  # Ohms
Cdl = 1e0  # Farads
Rp = 10  # Ohms

# Function to calculate the impedance
def calculate_impedance(omega):
    Z_parallel = (1 / Rp + 1j * omega * Cdl) ** -1
    Z_total = Rel + Z_parallel
    return Z_total

# Function to calculate current phasor given the voltage phasor
def calculate_current_phasor(V, Z):
    return V / Z

# User input for angular frequency (omega) and voltage amplitude
omega = float(input("Enter the angular frequency (omega) in rad/s: "))
amplitude = float(input("Enter the voltage amplitude: "))

# Calculate the total impedance for the given omega
Z_total = calculate_impedance(omega)

# Calculate the maximum current amplitude for axis scaling
max_current_amplitude = amplitude / abs(Z_total)

# Create a figure with four subplots: 2 rows x 2 columns
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Assign axes
# Left column: Time plots
ax_voltage_time = axes[0, 0]    # Top left
ax_current_time = axes[1, 0]    # Bottom left

# Right column: Phasor plots
ax_voltage_phasor = axes[0, 1]  # Top right
ax_current_phasor = axes[1, 1]  # Bottom right

# Set up the voltage vs time plot (ax_voltage_time)
ax_voltage_time.set_xlim(0, 2 * np.pi / omega)  # One period of the wave
ax_voltage_time.set_ylim(-1.5 * amplitude, 1.5 * amplitude)
ax_voltage_time.grid(True)
ax_voltage_time.set_xlabel('Time (s)')
ax_voltage_time.set_ylabel('Amplitude')
ax_voltage_time.set_title('Voltage vs. Time')

# Initialize time-domain line for voltage
time_data = []
voltage_real_data = []
voltage_time_line, = ax_voltage_time.plot([], [], 'r-', lw=2)

# Set up the current vs time plot (ax_current_time)
ax_current_time.set_xlim(0, 2 * np.pi / omega)  # One period of the wave
ax_current_time.set_ylim(-1.5 * max_current_amplitude, 1.5 * max_current_amplitude)
ax_current_time.grid(True)
ax_current_time.set_xlabel('Time (s)')
ax_current_time.set_ylabel('Amplitude')
ax_current_time.set_title('Current vs. Time')

# Initialize time-domain line for current
current_real_data = []
current_time_line, = ax_current_time.plot([], [], 'b-', lw=2)

# Set up the voltage phasor plot (ax_voltage_phasor)
phasor_limit_voltage = 1.5 * amplitude
ax_voltage_phasor.set_xlim(-phasor_limit_voltage, phasor_limit_voltage)
ax_voltage_phasor.set_ylim(-phasor_limit_voltage, phasor_limit_voltage)
ax_voltage_phasor.set_aspect('equal')
ax_voltage_phasor.grid(True)
ax_voltage_phasor.set_title('Voltage Phasor')
ax_voltage_phasor.set_xlabel('Real Part')
ax_voltage_phasor.set_ylabel('Imaginary Part')

# Initialize the voltage phasor line
voltage_phasor_line, = ax_voltage_phasor.plot([], [], 'r-', lw=2)

# Set up the current phasor plot (ax_current_phasor) with autoscaling
ax_current_phasor.set_aspect('equal')
ax_current_phasor.grid(True)
ax_current_phasor.set_title('Current Phasor')
ax_current_phasor.set_xlabel('Real Part')
ax_current_phasor.set_ylabel('Imaginary Part')

# Initialize the current phasor line
current_phasor_line, = ax_current_phasor.plot([], [], 'b-', lw=2)

# Initialize the animation function
def init():
    voltage_time_line.set_data([], [])
    current_time_line.set_data([], [])
    voltage_phasor_line.set_data([], [])
    current_phasor_line.set_data([], [])
    return voltage_time_line, current_time_line, voltage_phasor_line, current_phasor_line

# Animation update function
def update(frame):
    # Calculate the current phase angle based on frame number
    time = frame / fps  # in seconds
    phase_angle = omega * time

    # Voltage phasor as a rotating vector
    V = amplitude * np.exp(1j * phase_angle)
    voltage_phasor_line.set_data([0, V.real], [0, V.imag])

    # Current phasor based on voltage phasor and impedance
    I = calculate_current_phasor(V, Z_total)
    current_phasor_line.set_data([0, I.real], [0, I.imag])

    # Autoscale current phasor plot
    phasor_limit_current = 1.1 * np.abs(I)
    ax_current_phasor.set_xlim(-phasor_limit_current, phasor_limit_current)
    ax_current_phasor.set_ylim(-phasor_limit_current, phasor_limit_current)

    # Update time-domain data
    time_data.append(time)
    voltage_real_data.append(V.real)
    current_real_data.append(I.real)

    # Update voltage vs time plot
    voltage_time_line.set_data(time_data, voltage_real_data)

    # Update current vs time plot
    current_time_line.set_data(time_data, current_real_data)

    # Keep the time-domain plots within one period for better visualization
    period = 2 * np.pi / omega
    if time > 6 * period:
        # Stop the animation after 6 periods
        ani.event_source.stop()
    elif time > period:
        ax_voltage_time.set_xlim(time - period, time)
        ax_current_time.set_xlim(time - period, time)

    return voltage_time_line, current_time_line, voltage_phasor_line, current_phasor_line

# Calculate total frames needed for at least 6 periods
fps = 30  # frames per second
period = 2 * np.pi / omega
total_time = 6 * period
total_frames = int(np.ceil(total_time * fps))

# Define a function to start the animation
def start_animation(event):
    global ani
    ani = animation.FuncAnimation(
        fig, update, frames=total_frames,
        init_func=init, blit=False, interval=1000 / fps
    )
    plt.draw()

# Add a start button
ax_button = plt.axes([0.4, 0.01, 0.2, 0.05])  # Position for the button [left, bottom, width, height]
button = Button(ax_button, 'Start Simulation')
button.on_clicked(start_animation)

# Adjust layout to make sure plots occupy the same space
plt.subplots_adjust(wspace=0.3, hspace=0.3)

# Display the figure with the start button
plt.suptitle('Voltage and Current Phasor Animation and Time-Domain Plots')
plt.tight_layout(rect=[0, 0.05, 1, 0.95])
plt.show()
