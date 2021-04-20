from matplotlib import pyplot as plt
import numpy as np


x = np.arange(0, 4 * np.pi, np.pi/64)

f1 = np.sin(np.pi * x)
f2 = np.sin(2 * np.pi * x)
f3 = (1/3) * np.sin(3 * np.pi * x)


plt.plot(x, f1, label='sin(pi * x)')
plt.plot(x, f2, label='sin(2 * pi * x)')
plt.plot(x, f3, label='1/3 * sin(3 * pi * x)')
plt.plot(x, f1 + f2 + f3, label='SUM')
plt.legend(loc='lower left')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.show()
