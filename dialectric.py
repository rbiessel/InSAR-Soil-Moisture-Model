import numpy as np
from matplotlib import pyplot as plt
import scipy


mv = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55])
epsilon_r = np.array([4, 6, 7.5, 10, 13, 16.5, 19, 23, 26.5, 30, 34, 37.5])

# def interpolate(x0, y0, padding=0.5):
#     # Define a set of dense x values to apply the interpolated function on
#     x = np.linspace(x0[0] - padding, x0[-1] + padding, x0.size * 100)
#     # Create a multi-dimensional array to hold the interpolaion for each y value)
#     pjs = np.zeros((x0.size, x.size))
#     for j in range(y0.size):
#         pj = np.ones(x.shape) * y0[j]  # Create an individual row
#         for k in range(y0.size):
#             if j != k:
#                 # Iteratively perform the multiplication for each x value
#                 pj *= (x - x0[k]) / (x0[j] - x0[k])
#         pjs[j, :] = pj

#     # Sum down each column to get the total interpolated function
#     return x, np.sum(pjs, axis=0)


# plt.plot(mv, epsilon_r)
# plt.show()
