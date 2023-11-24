import math
import numpy as np
import matplotlib.pyplot as plt
import itertools as it


def function(x: float) -> float:
    return math.exp(x)


def lagrange_interpolation_polynomial(f: function, x_vals: list[float]) -> list[float]:
    n = len(x_vals)
    y_factors = []
    for i in range(n):
        y = f(x_vals[i])
        for j in range(n):
            if i != j:
                y /= (x_vals[i] - x_vals[j])
        y_factors.append(y)
    return y_factors


def calculate_lagrange_polynomial(x: float, y_factors: list[float], x_vals: list[float]) -> float:
    result = 0
    n = len(x_vals)
    for i in range(n):
        tmp = y_factors[i]
        for j in range(n):
            if i != j:
                tmp *= (x - x_vals[j])
        result += tmp
    return result


def newton_interpolation_polynomial(f: function, x_vals: list[float]) -> list[float]:
    n = len(x_vals)
    h = x_vals[1] - x_vals[0]
    delta_y = []
    for x in x_vals:
        delta_y.append(f(x))
    y_factors = [delta_y[0]]
    for k in range(1, n):
        delta_y = [delta_y[i + 1] - delta_y[i] for i in range(len(delta_y) - 1)]
        y_factors.append(delta_y[0] / (math.factorial(k) * math.pow(h, k)))
    return y_factors


def calculate_newton_polynomial(x: float, y_factors: list[float], x_vals: list[float]) -> float:
    n = len(x_vals)
    result = y_factors[0]
    for i in range(1, n):
        tmp = y_factors[i]
        for j in range(i):
            tmp *= (x - x_vals[j])
        result += tmp
    return result


def computation_error(f: function, x: float, res: float) -> tuple[float, float]:
    if f(x) != 0:
        return abs(f(x) - res), abs(f(x) - res) / f(x) * 100
    return abs(res), 100 if res != 0 else 0


def make_plot() -> None:
    x = np.arange(x_values_lagrange[0], x_values_lagrange[3], 0.01)
    plt.plot(x, [function(i) for i in x], label=r'$e^x$')
    plt.plot(x, [calculate_lagrange_polynomial(i, y_coefficients_lagrange, x_values_lagrange) for i in x],
             label=r'$Lagrange$')
    plt.plot(x, [calculate_newton_polynomial(i, y_coefficients_newton, x_values_newton) for i in x],
             label=r'$Newton$',
             )
    all_values = set(it.chain(x_values_lagrange, x_values_newton))
    plt.scatter(list(all_values), [function(i) for i in all_values])
    plt.xlabel(r'$x$', fontsize=14)
    plt.ylabel(r'$y$', fontsize=14)
    plt.grid(True)
    plt.legend(loc='best', fontsize=12)
    plt.show()


if __name__ == '__main__':

    x_values_lagrange = [-2., -1., 0., 1.]
    x_values_newton = [-2., -1., 0., 1.]
    x_aster = -0.5
    rounding = 4

    # lagrange
    y_coefficients_lagrange = lagrange_interpolation_polynomial(function, x_values_lagrange)
    value_lagrange = calculate_lagrange_polynomial(x_aster, y_coefficients_lagrange, x_values_lagrange)
    absolute_error_lagrange, relative_error_lagrange = computation_error(function, x_aster, value_lagrange)
    print(
        '\n--- Lagrange polynomial coefficients ---',
        '  '.join(str(round(y, rounding)) for y in y_coefficients_lagrange),
        f'\n--- Value of the Lagrange polynomial at the point x = {x_aster} ---',
        round(value_lagrange, rounding),
        f'\n--- Function at the point x = {x_aster} ---',
        round(function(x_aster), rounding),
        '\n--- Absolute error ---',
        round(absolute_error_lagrange, rounding),
        '\n--- Relative error ---',
        f'{round(relative_error_lagrange, rounding)} %',
        sep='\n'
    )

    print('\n\n######################################\n')

    # newton
    y_coefficients_newton = newton_interpolation_polynomial(function, x_values_newton)
    value_newton = calculate_newton_polynomial(x_aster, y_coefficients_newton, x_values_newton)
    absolute_error_newton, relative_error_newton = computation_error(function, x_aster, value_newton)
    print(
        '\n--- Newton polynomial coefficients ---',
        '  '.join(str(round(y, rounding)) for y in y_coefficients_newton),
        f'\n--- Value of the Newton polynomial at the point x = {x_aster} ---',
        round(value_newton, rounding),
        f'\n--- Function at the point x = {x_aster} ---',
        round(function(x_aster), rounding),
        '\n--- Absolute error ---',
        round(absolute_error_newton, rounding),
        '\n--- Relative error ---',
        f'{round(relative_error_newton, rounding)} %',
        sep='\n'
    )

    make_plot()
