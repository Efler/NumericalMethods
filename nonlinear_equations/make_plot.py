import math

import numpy as np
import matplotlib.pyplot as plt


def function_1(x: float) -> float:
    return math.exp(x)


def function_2(x: float) -> float:
    return 2 * x + 2


if __name__ == '__main__':
    x_vals = np.arange(-1.5, 2, 0.05)
    plt.plot(x_vals, [function_1(i) for i in x_vals], label=r'$e^x$')
    plt.plot(x_vals, [function_2(i) for i in x_vals], label=r'$2x + 2$')
    plt.xlabel(r'$x$', fontsize=14)
    plt.ylabel(r'$y$', fontsize=14)
    plt.grid(True)
    plt.legend(loc='best', fontsize=12)
    plt.show()
