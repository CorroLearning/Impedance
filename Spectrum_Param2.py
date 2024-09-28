import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

st.title("Circuit Impedance Visualization")

# Initialize session state for axis limits
if 'axis_limits_fig1' not in st.session_state:
    st.session_state.axis_limits_fig1 = None
if 'axis_limits_fig2' not in st.session_state:
    st.session_state.axis_limits_fig2 = None

st.header("Input Parameters")
st.write("Use the sliders below to adjust the component values. The sliders are mapped logarithmically.")

# Arrange sliders in a row
slider_col1, slider_col2, slider_col3 = st.columns(3)
with slider_col1:
    # Rel slider
    log_Rel = st.slider('Log10(Rel) [Ohms]', min_value=-1.0, max_value=3.0, value=1.0, step=0.1)
    Rel = 10 ** log_Rel
with slider_col2:
    # Rpol slider
    log_Rpol = st.slider('Log10(Rpol) [Ohms]', min_value=1.0, max_value=6.0, value=4.0, step=0.1)
    Rpol = 10 ** log_Rpol
with slider_col3:
    # Cdl slider
    log_Cdl = st.slider('Log10(Cdl) [Farads]', min_value=-9.0, max_value=-3.0, value=-6.0, step=0.1)
    Cdl = 10 ** log_Cdl

# Autoscale checkbox for Bode plots
autoscale = st.checkbox('Autoscale Bode Plots', value=True)

# Frequency range
f = np.logspace(-3, 5, num=1000)  # From 0.001 Hz to 100,000 Hz
omega = 2 * np.pi * f

# Calculate impedance components
Z_Cdl = 1 / (1j * omega * Cdl)
Z_parallel = 1 / (1 / Rpol + 1 / Z_Cdl)
Z_total = Rel + Z_parallel

# Calculate modulus and phase
Z_modulus = np.abs(Z_total)
Z_phase_degrees = np.degrees(np.angle(Z_total))  # Convert phase to degrees

# Prepare Bode plot (Modulus)
fig1, ax1 = plt.subplots(figsize=(5, 4))
ax1.plot(np.log10(f), np.log10(Z_modulus))
ax1.set_xlabel('Log Frequency (Hz)')
ax1.set_ylabel('Log Modulus of Impedance (Ω)')
ax1.set_title('Log Modulus of Impedance vs. Log Frequency')

# Handle axis limits for Bode Modulus plot
if not autoscale and st.session_state.axis_limits_fig1 is not None:
    ax1.set_xlim(st.session_state.axis_limits_fig1[0])
    ax1.set_ylim(st.session_state.axis_limits_fig1[1])
elif autoscale:
    st.session_state.axis_limits_fig1 = (ax1.get_xlim(), ax1.get_ylim())

# Prepare Bode plot (Phase)
fig2, ax2 = plt.subplots(figsize=(5, 4))
ax2.plot(np.log10(f), -Z_phase_degrees)  # Plot phase in degrees
ax2.set_xlabel('Log Frequency (Hz)')
ax2.set_ylabel('Negative Phase of Impedance (Degrees)')
ax2.set_title('Negative Phase vs. Log Frequency')

# Handle axis limits for Bode Phase plot
if not autoscale and st.session_state.axis_limits_fig2 is not None:
    ax2.set_xlim(st.session_state.axis_limits_fig2[0])
    ax2.set_ylim(st.session_state.axis_limits_fig2[1])
elif autoscale:
    st.session_state.axis_limits_fig2 = (ax2.get_xlim(), ax2.get_ylim())

# Prepare Nyquist plot
fig3, ax3 = plt.subplots(figsize=(5, 4))  # Matching size with other plots
ax3.plot(Z_total.real, -Z_total.imag)
ax3.set_xlabel('Real Part of Impedance (Ω)')
ax3.set_ylabel('Negative Imaginary Part of Impedance (Ω)')
ax3.set_title('Nyquist Plot')

# Ensure Nyquist plot always autoscales
ax3.autoscale()  # Always autoscale the Nyquist plot
ax3.set_aspect('equal', adjustable='datalim')  # Keeps aspect ratio equal

# Display plots in a row
st.header("Plots")
col1, col2, col3 = st.columns(3)
with col1:
    st.pyplot(fig1)
with col2:
    st.pyplot(fig2)
with col3:
    st.pyplot(fig3)

# Display component values in a row below the graphs
st.header("Component Values")
value_col1, value_col2, value_col3 = st.columns(3)
with value_col1:
    st.write(f'Rel = {Rel:.4g} Ω')
with value_col2:
    st.write(f'Rpol = {Rpol:.4g} Ω')
with value_col3:
    st.write(f'Cdl = {Cdl:.4g} F')
