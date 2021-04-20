import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from matplotlib import rc
rc('animation', html='jshtml')


def initial_wave(x, x0, A, h):
    u = np.zeros(len(x))
    psi = np.zeros(len(x))
    return u, psi


def boundary_conditions(u, psip, t, omega=1):

    h = 30
    envelope = np.exp(-(t - (3*h))**2 / h**2)

    u[0] = np.cos(omega * t) * envelope
    u[-1] = 0
    psip[0] = np.sin(omega * t) * envelope
    psip[-1] = 0

    return u, psip


def solve_wave(u, psi, tfinal, dt):
    t = 0  # Intialize time and empty arrays to gold time and displacements
    uarr = []
    tarr = []
    psip = psi
    while t < tfinal:
        u, psip = boundary_conditions(u, psip, t)

        psip[1:-1] = psi[1:-1] + s[1:-1] * (u[2:] - 2 * u[1:-1] + u[:-2])

        u = u + dt * psi  # Use the previous value of psi to progress the wave displacements

        # if t < 20 and t > 1:
        #     u[1] = np.sin(4 * np.pi * t)

        psi = psip  # update psi
        uarr.append(u)  # store displacements and time
        tarr.append(t)
        t += dt
    return uarr, tarr


def solve_wave2(u0, psi0, tfinal, dt):
    tarr = np.arange(0, tfinal, dt)
    uarr = np.zeros((tarr.size, psi0.size))
    psis = np.zeros((tarr.size, psi0.size))

    # set initial conditions
    psis[0, :] = psi0
    uarr[0, :] = u0

    for t in range(0, len(tarr) - 1):  # loop through each time step
        # update the displacements
        uarr[t + 1] = uarr[t] + dt * psis[t]

        for x in range(0, uarr.shape[1] - 1):  # loop through each spatial step
            # find next psi
            # forward time, centered space

            psis[t + 1, x] = psis[t, x] + \
                (s[x - 1] * uarr[t, x - 1] + s[x] * 2 *
                 uarr[t, x] + s[x + 1] * uarr[t, x + 1])

    return uarr, tarr


def animate_func(x, uarr, tarr, A, nstep):
    plt.style.use('dark_background')
    fig, ax = plt.subplots(figsize=(14, 6))

    line, = ax.plot([], [])
    ax.axvspan(0, 300, color='gray', alpha=0.2)

    ax.set_ylim(-A, A)
    ax.set_xlim(x.min(), x.max())
    ax.set_xlabel('x')
    ax.set_ylabel('f(x)')

    def update_line(i):  # Get frame to animate
        line.set_ydata(data[i])
        line.set_xdata(x)
        return line,

    data = []

    for i in range(0, len(tarr), nstep):
        data.append(uarr[i])

    nfrm = int(len(data))
    ani = animation.FuncAnimation(
        fig, update_line, frames=nfrm, interval=1, blit=True, repeat=True, cache_frame_data=True)

    plt.show()
    return ani
# -------------------------------------------------------------------------------------------------------------------------


nx = 400.0  # number of points along the line
dx = 0.1  # spatial grid size
x = np.arange(0, nx, dx)  # generate spatial grid
x0 = 100  # Initial position of wave
h = 30.0  # width of wave "standard deviation" of the Gaussian
A = 1.0  # wave amplitude
v = 1.0  # wave velocity m/s
tfinal = 1000

# setup the s-array to be dependend on position
s1 = 3.0
s2 = 0.1
s = np.ones(x.shape)

s[600::] = s2

# s[-1] = 0
# s[0] = 0

dt = 1 * dx**2 / v**2

u, psi = initial_wave(x, x0, A, h)
uarr, tarr = solve_wave(u, psi, tfinal, dt=dt)

exit
cnt = 20

fig, ax = plt.subplots(figsize=(10, 5))

ax.plot(x, np.sum(np.array(uarr), axis=0))
plt.show()

nstep = 10
ani = animate_func(x, uarr, tarr, A, nstep)
ani
