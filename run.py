import sm as sm_model
from matplotlib import pyplot as plt
import numpy as np
import dialectric

'''
Execute the model/simulation and display the associated figures and animations.

'''

# Sentinel-1 Altitude in meters arbitrary as long as it's a magnitude larger than the wavelength
satellite_altitude = 50

t_final = 250
envelope = 10
dx = .075
sample_time = 100

delta_mvs = np.zeros(len(dialectric.mv))
delta_eps = np.zeros(len(dialectric.mv))
delta_phis = np.zeros(len(dialectric.mv))

# np.linspace(.1, .6, 10): # use a single value or remove comment to try different depths
for depth in [0.5]:
    for i in range(0, len(dialectric.mv)):
        depth = np.round(depth, 2)

        # # get a sample pair of soil moistures
        mv1 = dialectric.mv[0]
        mv2 = dialectric.mv[i]

        # # Related Permittivities
        epsilon_1 = dialectric.epsilon_r[0]
        epsilon_2 = dialectric.epsilon_r[i]

        # Define soil mediums

        soil0 = sm_model.medium(
            satellite_altitude, depth, epsilon_2)
        soil1 = sm_model.medium(
            satellite_altitude, depth, epsilon_1)

        # Setup Simulation

        wave0 = sm_model.microwave(
            soil0, envelope=envelope, tfinal=t_final, dx=dx)
        wave1 = sm_model.microwave(
            soil1, envelope=envelope, tfinal=t_final, dx=dx)

        # Solve wave equation
        wave0.solve()
        wave1.solve()

        sm_model.animate_multiple([wave0, wave1], 20)

        # Find a time to sample at
        index = int(sample_time / wave1.dt)

        # Sample waveforms at time t
        waveform0 = np.array(wave0.uarr[index])
        waveform1 = np.array(wave1.uarr[index])

        # Calculate a wrapped offset in multiples of pi
        offset = wave0.x[np.argmax(waveform0)] - wave0.x[np.argmax(waveform1)]
        print(f'Offset in cm: {offset}')
        offset *= (2 / sm_model.sentinel1_wavelength)
        print(f'Offset in pi radins: {offset}')
        offset = (offset % 2)
        offset = np.round(offset, 2)

        # Plot waveforms
        plt.plot(wave0.x, np.array(wave0.uarr[index]), label=f'{mv1}mv')
        plt.plot(wave1.x, np.array(wave1.uarr[index]), label=f'{mv2}mv')
        plt.legend(loc='upper right')
        plt.xlabel('x (cm)')
        plt.ylabel('E Field Magnitude')
        plt.title(
            f'??1 = {epsilon_1}, ??2 = {epsilon_2}, Phase Offset: {offset} ??')
        plt.show()

        # get plottable values
        delta_mv = np.abs(mv1 - mv2)
        delta_ep = np.abs(epsilon_1 - epsilon_2)
        delta_phi = offset

        print(f'Delta MV: {delta_mv}, Delta Phi: {delta_phi}')

        delta_eps[i] = delta_ep
        delta_mvs[i] = delta_mv
        delta_phis[i] = delta_phi

    delta_phis[delta_phis < 0] += 2

    plt.plot(delta_mvs, delta_phis, label=f'Delta phi (depth = {depth})')

plt.plot(delta_mvs, delta_eps / delta_eps.max(),
         label='Delta Epsilon / Epsilon')

plt.legend(bbox_to_anchor=(0, 1), loc='upper left', ncol=1)
plt.xlabel('Change in % Water Volume')
plt.ylabel('Interferometric Phase')

plt.show()
