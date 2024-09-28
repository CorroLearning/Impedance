import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import time

st.markdown("<h1 style='font-size:24px;'>This is a Smaller Header</h1>", unsafe_allow_html=True)


def initialize_session_state():
    # Initialize session state for axis limits and accumulated data
    if 'axis_limits_fig1' not in st.session_state:
        st.session_state.axis_limits_fig1 = None
    if 'axis_limits_fig2' not in st.session_state:
        st.session_state.axis_limits_fig2 = None
    if 'axis_limits_fig3' not in st.session_state:
        st.session_state.axis_limits_fig3 = None
    if 'freq_list' not in st.session_state:
        st.session_state.freq_list = []
    if 'modulus_list' not in st.session_state:
        st.session_state.modulus_list = []
    if 'phase_list' not in st.session_state:
        st.session_state.phase_list = []
    if 'real_list' not in st.session_state:
        st.session_state.real_list = []
    if 'imag_list' not in st.session_state:
        st.session_state.imag_list = []

def input_parameters():
    st.header("Input Parameters")
    st.write("Use the sliders below to adjust the component values. The sliders are mapped logarithmically.")

    # Rel slider
    log_Rel = st.slider('Log10(Rel) [Ohms]', min_value=-1.0, max_value=3.0, value=1.0, step=0.1)
    Rel = 10 ** log_Rel
    st.write(f'Rel = {Rel:.4g} Ω')

    # Rpol slider
    log_Rpol = st.slider('Log10(Rpol) [Ohms]', min_value=1.0, max_value=6.0, value=4.0, step=0.1)
    Rpol = 10 ** log_Rpol
    st.write(f'Rpol = {Rpol:.4g} Ω')

    # Cdl slider
    log_Cdl = st.slider('Log10(Cdl) [Farads]', min_value=-9.0, max_value=-3.0, value=-6.0, step=0.1)
    Cdl = 10 ** log_Cdl
    st.write(f'Cdl = {Cdl:.4g} F')

    # Autoscale checkbox
    autoscale = st.checkbox('Autoscale', value=True)

    return Rel, Rpol, Cdl, autoscale

def voltage_signal_parameters():
    st.subheader("Voltage Signal Parameters")
    st.write("Adjust the amplitude and frequency of the applied voltage signal.")

    # Amplitude slider (linear mapping)
    V = st.slider('Amplitude V [Volts]', min_value=0.0, max_value=0.05, value=0.005, step=0.001)

    # Frequency slider (logarithmic mapping)
    log_F = st.slider('Log10(Frequency) [Hz]', min_value=-3.0, max_value=5.0, value=1.0, step=0.1)
    F = 10 ** log_F
    st.write(f'Frequency F = {F:.4g} Hz')

    return V, F

def simulation_parameters():
    st.subheader("Simulation Parameters")
    start_freq = st.number_input('Start Frequency (Hz)', value=100000.0, min_value=0.001, max_value=100000.0)
    end_freq = st.number_input('End Frequency (Hz)', value=0.001, min_value=0.001, max_value=100000.0)
    simulate = st.button('Simulate')
    return start_freq, end_freq, simulate

def create_plot_placeholders():
    st.markdown("<h1 style='font-size:16px;'>Frequency Domain</h1>", unsafe_allow_html=True)
    col1, col2, col3 = st.columns(3)
    with col1:
        plot_area1 = st.empty()
    with col2:
        plot_area2 = st.empty()
    with col3:
        plot_area3 = st.empty()

    st.markdown("<h1 style='font-size:16px;'>Time Domain</h1>", unsafe_allow_html=True)
    col4, col5, col6 = st.columns(3)
    with col4:
        plot_area4 = st.empty()
    with col5:
        plot_area5 = st.empty()
    with col6:
        plot_area6 = st.empty()

    st.markdown("<h1 style='font-size:16px;'>Phasors</h1>", unsafe_allow_html=True)
    col7, col8, col9 = st.columns(3)
    with col7:
        plot_area7 = st.empty()
    with col8:
        plot_area8 = st.empty()
    with col9:
        plot_area9 = st.empty()

    return (plot_area1, plot_area2, plot_area3, plot_area4, plot_area5, plot_area6, plot_area7, plot_area8, plot_area9)

def calculate_and_plot(F_sim, Rel, Rpol, Cdl, V, autoscale, plot_areas, simulate=False):
    # Unpack plot areas
    plot_area1, plot_area2, plot_area3, plot_area4, plot_area5, plot_area6, plot_area7, plot_area8, plot_area9 = plot_areas

    # Frequency range for frequency domain plots
    f = np.logspace(-3, 5, num=1000)  # From 0.001 Hz to 100,000 Hz
    omega = 2 * np.pi * f

    # Calculate impedance components over frequency range
    Z_Cdl = 1 / (1j * omega * Cdl)
    Z_parallel = 1 / (1 / Rpol + 1 / Z_Cdl)
    Z_total = Rel + Z_parallel

    # Calculate modulus and phase over frequency range
    Z_modulus = np.abs(Z_total)
    Z_phase = np.angle(Z_total, deg=False)  # In radians

    # Compute Z_total at single frequency F_sim
    omega_single = 2 * np.pi * F_sim
    Z_Cdl_single = 1 / (1j * omega_single * Cdl)
    Z_parallel_single = 1 / (1 / Rpol + 1 / Z_Cdl_single)
    Z_total_single = Rel + Z_parallel_single

    Z_modulus_single = np.abs(Z_total_single)
    Z_phase_single = np.angle(Z_total_single)  # in radians

    # Append data to accumulated lists if simulating
    if simulate:
        st.session_state.freq_list.append(F_sim)
        st.session_state.modulus_list.append(np.log10(Z_modulus_single))
        st.session_state.phase_list.append(-Z_phase_single)
        st.session_state.real_list.append(Z_total_single.real)
        st.session_state.imag_list.append(-Z_total_single.imag)

        # For the frequency domain plots, plot the accumulated data
        fig1, ax1 = plt.subplots()
        ax1.plot(np.log10(st.session_state.freq_list), st.session_state.modulus_list, 'o-')
        ax1.set_xlabel('Log Frequency (Hz)')
        ax1.set_ylabel('Log Modulus of Impedance (Ω)')
        ax1.set_title('Log Modulus of Impedance vs. Log Frequency')
        # Add vertical line at current frequency
        ax1.axvline(x=np.log10(F_sim), color='red', linestyle='--', label=f'F = {F_sim:.4g} Hz')
        ax1.legend()

        # Handle axis limits for fig1
        if not autoscale and st.session_state.axis_limits_fig1 is not None:
            ax1.set_xlim(st.session_state.axis_limits_fig1[0])
            ax1.set_ylim(st.session_state.axis_limits_fig1[1])
        elif autoscale:
            st.session_state.axis_limits_fig1 = (ax1.get_xlim(), ax1.get_ylim())

        # Similarly for fig2
        fig2, ax2 = plt.subplots()
        ax2.plot(np.log10(st.session_state.freq_list), st.session_state.phase_list, 'o-')
        ax2.set_xlabel('Log Frequency (Hz)')
        ax2.set_ylabel('Negative Phase of Impedance (Radians)')
        ax2.set_title('Negative Phase vs. Log Frequency')
        # Add vertical line at current frequency
        ax2.axvline(x=np.log10(F_sim), color='red', linestyle='--', label=f'F = {F_sim:.4g} Hz')
        ax2.legend()

        # Handle axis limits for fig2
        if not autoscale and st.session_state.axis_limits_fig2 is not None:
            ax2.set_xlim(st.session_state.axis_limits_fig2[0])
            ax2.set_ylim(st.session_state.axis_limits_fig2[1])
        elif autoscale:
            st.session_state.axis_limits_fig2 = (ax2.get_xlim(), ax2.get_ylim())

        # For fig3 (Nyquist plot), plot the accumulated data
        fig3, ax3 = plt.subplots()
        ax3.plot(st.session_state.real_list, st.session_state.imag_list, 'o-')
        ax3.set_xlabel('Real Part of Impedance (Ω)')
        ax3.set_ylabel('Negative Imaginary Part of Impedance (Ω)')
        ax3.set_title('Nyquist Plot')

        # Set aspect ratio first
        ax3.set_aspect('equal', adjustable='box')

        # Ensure that x and y axes have the same range
        xlims = ax3.get_xlim()
        ylims = ax3.get_ylim()
        x_range = xlims[1] - xlims[0]
        y_range = ylims[1] - ylims[0]
        max_range = max(x_range, y_range)
        x_middle = (xlims[0] + xlims[1]) / 2
        y_middle = (ylims[0] + ylims[1]) / 2
        ax3.set_xlim(x_middle - max_range / 2, x_middle + max_range / 2)
        ax3.set_ylim(y_middle - max_range / 2, y_middle + max_range / 2)

    else:
        # Not simulating, plot the full frequency response
        # Prepare frequency domain plots
        fig1, ax1 = plt.subplots()
        ax1.plot(np.log10(f), np.log10(Z_modulus))
        ax1.set_xlabel('Log Frequency (Hz)')
        ax1.set_ylabel('Log Modulus of Impedance (Ω)')
        ax1.set_title('Log Modulus of Impedance vs. Log Frequency')
        # Add vertical line at current frequency
        ax1.axvline(x=np.log10(F_sim), color='red', linestyle='--', label=f'F = {F_sim:.4g} Hz')
        ax1.legend()

        # Handle axis limits for fig1
        if not autoscale and st.session_state.axis_limits_fig1 is not None:
            ax1.set_xlim(st.session_state.axis_limits_fig1[0])
            ax1.set_ylim(st.session_state.axis_limits_fig1[1])
        elif autoscale:
            st.session_state.axis_limits_fig1 = (ax1.get_xlim(), ax1.get_ylim())

        fig2, ax2 = plt.subplots()
        ax2.plot(np.log10(f), -Z_phase)
        ax2.set_xlabel('Log Frequency (Hz)')
        ax2.set_ylabel('Negative Phase of Impedance (Radians)')
        ax2.set_title('Negative Phase vs. Log Frequency')
        # Add vertical line at current frequency
        ax2.axvline(x=np.log10(F_sim), color='red', linestyle='--', label=f'F = {F_sim:.4g} Hz')
        ax2.legend()

        # Handle axis limits for fig2
        if not autoscale and st.session_state.axis_limits_fig2 is not None:
            ax2.set_xlim(st.session_state.axis_limits_fig2[0])
            ax2.set_ylim(st.session_state.axis_limits_fig2[1])
        elif autoscale:
            st.session_state.axis_limits_fig2 = (ax2.get_xlim(), ax2.get_ylim())

        fig3, ax3 = plt.subplots()
        ax3.plot(Z_total.real, -Z_total.imag)
        ax3.set_xlabel('Real Part of Impedance (Ω)')
        ax3.set_ylabel('Negative Imaginary Part of Impedance (Ω)')
        ax3.set_title('Nyquist Plot')

        # Set aspect ratio first
        ax3.set_aspect('equal', adjustable='box')

        # Ensure that x and y axes have the same range
        xlims = ax3.get_xlim()
        ylims = ax3.get_ylim()
        x_range = xlims[1] - xlims[0]
        y_range = ylims[1] - ylims[0]
        max_range = max(x_range, y_range)
        x_middle = (xlims[0] + xlims[1]) / 2
        y_middle = (ylims[0] + ylims[1]) / 2
        ax3.set_xlim(x_middle - max_range / 2, x_middle + max_range / 2)
        ax3.set_ylim(y_middle - max_range / 2, y_middle + max_range / 2)

    # Time domain calculations
    # Time vector for 3 periods
    t_max = 3 / F_sim  # Total time for 3 periods
    t = np.linspace(0, t_max, num=1000)  # Time vector

    # Voltage and current as functions of time
    V_t = V * np.sin(2 * np.pi * F_sim * t)
    I_magnitude = V / Z_modulus_single
    I_t = I_magnitude * np.sin(2 * np.pi * F_sim * t + Z_phase_single)

    # Convert current to microamps for time-domain plot
    I_t_uA = I_t * 1e6  # Convert to microamps

    # Prepare time domain plots
    fig4, ax4 = plt.subplots()
    ax4.plot(t, V_t)
    ax4.set_xlabel('Time (s)')
    ax4.set_ylabel('Voltage (V)')
    ax4.set_title('Applied Voltage vs. Time')
    ax4.grid(True)

    fig5, ax5 = plt.subplots()
    ax5.plot(t, I_t_uA)
    ax5.set_xlabel('Time (s)')
    ax5.set_ylabel('Current (µA)')
    ax5.set_title('Resulting Current vs. Time')
    ax5.grid(True)

    fig6, ax6 = plt.subplots()
    ax6.plot(V_t, I_t_uA)
    ax6.set_xlabel('Voltage (V)')
    ax6.set_ylabel('Current (µA)')
    ax6.set_title('Lissajous Figure (Current vs. Voltage)')
    ax6.grid(True)

    # Prepare phasor plots
    # Voltage phasor in millivolts
    V_mV = V * 1000  # Convert to mV
    fig7, ax7 = plt.subplots()
    ax7.arrow(0, 0, V_mV, 0, head_width=V_mV * 0.05, head_length=V_mV * 0.1, fc='blue', ec='blue')
    ax7.set_xlabel('Real (mV)')
    ax7.set_ylabel('Imaginary (mV)')
    ax7.set_title('Voltage Phasor (mV)')
    ax7.grid(True)

    # Set axis limits for voltage phasor
    V_max = abs(V_mV) * 1.2
    if V_max == 0:
        V_max = 1  # Avoid zero range
    ax7.set_xlim(-V_max, V_max)
    ax7.set_ylim(-V_max, V_max)
    ax7.set_aspect('equal', adjustable='box')

    # Current phasor in microamps
    I_magnitude_uA = I_magnitude * 1e6  # Convert to microamps
    I_real_uA = I_magnitude_uA * np.cos(Z_phase_single)
    I_imag_uA = -I_magnitude_uA * np.sin(Z_phase_single)
    fig8, ax8 = plt.subplots()
    I_arrow_length_uA = np.hypot(I_real_uA, I_imag_uA)

    # Set head_width and head_length proportional to arrow length with limits
    I_max_uA = max(abs(I_real_uA), abs(I_imag_uA)) * 1.2
    if I_max_uA == 0:
        I_max_uA = 1  # Avoid zero range
    min_head_size = I_max_uA * 0.02
    max_head_size = I_max_uA * 0.1
    head_width = I_arrow_length_uA * 0.05
    head_length = I_arrow_length_uA * 0.1
    head_width = np.clip(head_width, min_head_size, max_head_size)
    head_length = np.clip(head_length, min_head_size, max_head_size)

    ax8.arrow(0, 0, I_real_uA, I_imag_uA, head_width=head_width, head_length=head_length, fc='red', ec='red')
    ax8.set_xlabel('Real (µA)')
    ax8.set_ylabel('Imaginary (µA)')
    ax8.set_title('Current Phasor (µA)')
    ax8.grid(True)

    # Set axis limits for current phasor
    ax8.set_xlim(-I_max_uA, I_max_uA)
    ax8.set_ylim(-I_max_uA, I_max_uA)
    ax8.set_aspect('equal', adjustable='box')

    # Impedance phasor
    Z_real = Z_total_single.real
    Z_imag = Z_total_single.imag
    fig9, ax9 = plt.subplots()
    Z_arrow_length = np.hypot(Z_real, Z_imag)

    # Set head_width and head_length proportional to arrow length with limits
    Z_max = max(abs(Z_real), abs(Z_imag)) * 1.2
    if Z_max == 0:
        Z_max = 1  # Avoid zero range
    min_head_size_Z = Z_max * 0.02
    max_head_size_Z = Z_max * 0.1
    head_width_Z = Z_arrow_length * 0.05
    head_length_Z = Z_arrow_length * 0.1
    head_width_Z = np.clip(head_width_Z, min_head_size_Z, max_head_size_Z)
    head_length_Z = np.clip(head_length_Z, min_head_size_Z, max_head_size_Z)

    ax9.arrow(0, 0, Z_real, Z_imag, head_width=head_width_Z, head_length=head_length_Z, fc='green', ec='green')
    ax9.set_xlabel('Real')
    ax9.set_ylabel('Imaginary')
    ax9.set_title('Impedance Phasor')
    ax9.grid(True)

    # Set axis limits for impedance phasor
    ax9.set_xlim(-Z_max, Z_max)
    ax9.set_ylim(-Z_max, Z_max)
    ax9.set_aspect('equal', adjustable='box')

    # Update plots
    with plot_area1:
        st.pyplot(fig1)
    with plot_area2:
        st.pyplot(fig2)
    with plot_area3:
        st.pyplot(fig3)
    with plot_area4:
        st.pyplot(fig4)
    with plot_area5:
        st.pyplot(fig5)
    with plot_area6:
        st.pyplot(fig6)
    with plot_area7:
        st.pyplot(fig7)
    with plot_area8:
        st.pyplot(fig8)
    with plot_area9:
        st.pyplot(fig9)

# Main application code
st.title("Circuit Impedance Visualization")

initialize_session_state()

Rel, Rpol, Cdl, autoscale = input_parameters()

V, F = voltage_signal_parameters()

start_freq, end_freq, simulate = simulation_parameters()

plot_areas = create_plot_placeholders()

if simulate:
    # Reset accumulated data
    st.session_state.freq_list = []
    st.session_state.modulus_list = []
    st.session_state.phase_list = []
    st.session_state.real_list = []
    st.session_state.imag_list = []

    # Ensure start_freq is greater than end_freq
    if start_freq < end_freq:
        st.error('Start Frequency must be greater than End Frequency.')
    else:
        # Calculate number of decades
        decades = np.log10(start_freq) - np.log10(end_freq)
        num_freqs = int(decades * 5) + 1  # 5 values per decade
        frequencies = np.logspace(np.log10(start_freq), np.log10(end_freq), num=num_freqs)

        progress_bar = st.progress(0)
        progress_step = 100 / len(frequencies)

        for idx, F_sim in enumerate(frequencies):
            calculate_and_plot(F_sim, Rel, Rpol, Cdl, V, autoscale, plot_areas, simulate=True)

            # Update progress bar
            progress_bar.progress(int((idx + 1) * progress_step))

            # Introduce a short delay (optional)
            time.sleep(0.1)
else:
    calculate_and_plot(F, Rel, Rpol, Cdl, V, autoscale, plot_areas, simulate=False)
