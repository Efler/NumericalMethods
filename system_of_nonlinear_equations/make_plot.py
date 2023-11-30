import math

import numpy as np
import matplotlib.pyplot as plt


# region test functions
# def function_1(x: float) -> float:
#     return ((0.3 - 0.1 * math.pow(x, 2) - x) / 0.2) ** 0.5
#
#
# def function_2(x: float) -> float:
#     return (0.7 - 0.2 * math.pow(x, 2)) / (1 - 0.1 * x)
# endregion


def function_1(x: float) -> float:
    return math.acos(x - 1)


def function_2(x: float) -> float:
    return math.log(x + 1, 10) + 3


if __name__ == '__main__':
    x_1_vals = np.linspace(0.0, 2.0, 1000)
    x_2_vals = np.linspace(-0.9, math.pi / 2, 1000)
    plt.plot(x_1_vals, [function_1(i) for i in x_1_vals], label=r'$x_1 - cos(x_2) = 1$')
    plt.plot(x_2_vals, [function_2(i) for i in x_2_vals], label=r'$x_2 - lg(x_1 + 1) = 3$')
    plt.xlabel(r'$x$', fontsize=14)
    plt.ylabel(r'$y$', fontsize=14)
    plt.grid(True)
    plt.legend(loc='best', fontsize=12)
    plt.show()


# 0.00943991
# 3.00408047
