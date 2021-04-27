import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from matplotlib import rc
rc('animation', html='jshtml')

'''
PHYS220 Computational Physics Final Project
by Rowan Biessel

Modeling the effect of a dielectric change on microwave propogation in
order to illustrate the effect of soil-moisture changes on repeat-pass
radar interferometric phase.

This file contains all of the core code for the soil-moisture model.
'''


sentinel1_wavelength = 5.66  # cm


class medium:
    def __init__(self, air_length, soil_length, dielectric):
        self.air_length = air_length  # m
        self.soil_length = soil_length  # m
        self.dielectric = dielectric


class microwave:
    def __init__(self, medium, tfinal=1000, dx=.1, envelope=1):

        # Wave Parameters

        self.omega_0 = 2 * np.pi / sentinel1_wavelength
        self.envelope_width = envelope  # envelope width of pulse
        self.medium = medium

        # Numerical Scheme Setup and Initial Conditions

        self.tfinal = tfinal

        self.nx = 400.0  # number of points along the line
        self.dx = dx  # spatial grid size
        self.x = np.arange(0, medium.air_length + medium.soil_length,
                           self.dx)  # generate spatial grid

        # Wave Initial Conditions
        self.u = np.zeros(len(self.x))
        self.psi = np.zeros(len(self.x))

        # setup the s-array to be dependend on position, find an appropriate delta t
        s1 = 1
        s2 = s1 / self.medium.dielectric

        s = np.ones(self.x.shape) * s1

        self.medium_boundary = int(self.medium.air_length / self.dx)

        s[self.medium_boundary::] = s2
        self.dt = 1 * dx**2 / s1**2

        print(f'dt: {self.dt}')

        self.s = s

        # setup arrays to hold solution
        self.uarr = []
        self.tarr = []

    def boundary_conditions(self, u, psip, t):
        h = self.envelope_width
        envelope = np.exp(-(t - (3*h))**2 / h**2)

        # Oscillate boundary to produce singular wave packet
        u[0] = np.sin(self.omega_0 * t) * envelope
        psip[0] = np.cos(self.omega_0 * t) * envelope

        u[-1] = 0
        psip[-1] = 0

        return u, psip

    def solve(self):
        t = 0
        psip = self.psi
        u = self.u
        while t < self.tfinal:
            u, psip = self.boundary_conditions(u, psip, t)

            psip[1:-1] = self.psi[1:-1] + self.s[1:-1] * \
                (u[2:] - 2 * u[1:-1] + u[:-2])

            # Use the previous value of psi to progress the wave displacements
            u = u + self.dt * self.psi

            self.psi = psip  # update psi
            self.uarr.append(u)  # store displacements and time
            self.tarr.append(t)
            t += self.dt

    def animate_func(self, nstep):
        # plt.style.use('dark_background')
        fig, ax = plt.subplots(figsize=(14, 6))

        line, = ax.plot([], [])
        ax.axvspan(0, self.medium.air_length, color='gray', alpha=0.2)

        ax.set_ylim(-2, 2)
        ax.set_xlim(self.x.min(), self.x.max())
        ax.set_xlabel('x (cm)')
        ax.set_ylabel('E_Field(x)')
        annotation = ax.annotate('Frame', xy=(2, 1.8))
        annotation.set_animated(True)

        def init():
            return line, annotation

        def update_line(i):  # Get frame to animate
            line.set_ydata(data[i])
            line.set_xdata(self.x)
            annotation.set_text(f'Time: {np.round(i * nstep * self.dt)}')
            return line, annotation

        data = []

        for i in range(0, len(self.tarr), nstep):
            data.append(self.uarr[i])

        nfrm = int(len(data))
        ani = animation.FuncAnimation(
            fig, update_line, init_func=init, frames=nfrm, interval=1, blit=True, repeat=True, cache_frame_data=True)

        plt.show()
        return ani


def animate_multiple(waves, nstep):

    # plt.style.use('dark_background')
    fig, ax = plt.subplots(figsize=(14, 6))

    lines = [ax.plot([], [])[0] for _ in range(len(waves))]  # lines to animate

    ax.axvspan(0, waves[0].medium.air_length, color='gray', alpha=0.2)

    ax.set_ylim(-2, 2)
    ax.set_xlim(waves[0].x.min(), waves[0].x.max())
    ax.set_xlabel('x (cm)')
    ax.set_ylabel('E_Field(x)')
    annotation = ax.annotate('Frame', xy=(2, 1.8))
    annotation.set_animated(True)

    def init():
        for line in lines:
            line.set_data([], [])
        return lines + [annotation]

    def update_line(i):  # Get frame to animate
        for j, line in enumerate(lines):
            line.set_ydata(data[i][j])
            line.set_xdata(np.array(waves[0].x))

        annotation.set_text(f'Time: {np.round(i * nstep * waves[0].dt)}')
        return lines + [annotation]

    data = []

    for i in range(0, len(waves[0].tarr), nstep):
        snapshot = []
        for wave in waves:
            snapshot.append(wave.uarr[i])

        data.append(snapshot)

    nfrm = int(len(data))
    ani = animation.FuncAnimation(
        fig, update_line, init_func=init, frames=nfrm, interval=1, blit=True, repeat=True, cache_frame_data=True)

    plt.show()
    return ani
