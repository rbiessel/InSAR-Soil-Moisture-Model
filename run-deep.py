import sm as sm_model
from matplotlib import pyplot as plt
import numpy as np
import dialectric

'''
This file is for demonstration purposes and highlights the effect a change in dielectric has on a
single wave. It should be observed that there is reflection at the boundary as well as the transmission.

'''

# Sentinel-1 Altitude in meters arbitrary as long as it's a magnitude larger than the wavelength
satellite_altitude = 10

t_final = 250
envelope = 2

delta_mvs = np.zeros(len(dialectric.mv))
delta_eps = np.zeros(len(dialectric.mv))
delta_phis = np.zeros(len(dialectric.mv))

for depth in [5]:  # np.linspace(.1, .6, 10):
    for i in range(0, len(dialectric.mv)):
        depth = np.round(depth, 2)
        # # get a sample pair of soil moistures
        mv1 = dialectric.mv[0]
        mv2 = dialectric.mv[i]

        # # Related Permittivities
        epsilon_1 = dialectric.epsilon_r[0]
        epsilon_2 = dialectric.epsilon_r[i]

        epsilon_1 = 3
        epsilon_2 = 1

        # Define soil mediums

        soil0 = sm_model.medium(
            satellite_altitude, depth, epsilon_2)
        soil1 = sm_model.medium(
            satellite_altitude, depth, epsilon_1)

        # Setup Simulation

        dx = .05

        wave0 = sm_model.microwave(
            soil0, envelope=envelope, tfinal=t_final, dx=dx)
        wave1 = sm_model.microwave(
            soil1, envelope=envelope, tfinal=t_final, dx=dx)

        wave1.omega_0 = 2 * np.pi / 1

        # Solve wave equation
        wave0.solve()
        wave1.solve()

        sm_model.animate_multiple([wave1], 20)
        break
        # Find a time to sample at
        index = int(100 / wave1.dt)
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
        if True:
            plt.plot(wave0.x, np.array(wave0.uarr[index]), label='mv1')
            plt.plot(wave1.x, np.array(wave1.uarr[index]), label='mv2')
            plt.legend(loc='upper right')
            plt.xlabel('x (cm)')
            plt.ylabel('E Field Magnitude')
            plt.title(
                f'ε1 = {epsilon_1}, ε2 = {epsilon_2}, Phase Offset: {offset} π')
            plt.show()

        delta_mv = np.abs(mv1 - mv2)
        delta_ep = np.abs(epsilon_1 - epsilon_2)
        delta_phi = offset

        print(f'Delta MV: {delta_mv}, Delta Phi: {delta_phi}')

        delta_eps[i] = delta_ep
        delta_mvs[i] = delta_mv
        delta_phis[i] = delta_phi

    delta_phis[delta_phis < 0] += 2

    plt.plot(delta_mvs, delta_phis, label=f'Delta phi (depth = {depth})')
