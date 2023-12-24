import math

import numpy as np
import pylab
from tabulate import tabulate

from mnk import mnk
from thermodynamic_properties_of_substances import ThermodynamicPropertiesOfSubstances as TPoS


def make_plot(a: float, b: float, props: list, width: float, height: float) -> None:
    props_titles = {
        'Cp': 'Cp',
        'Fi': 'Ф',
        'S': 'S',
        'delta_H': 'ΔH'
    }
    temps = my_props.temperatures_for_element(current_element)
    temperature_segment = temps[temps.index(a):temps.index(b) + 1]
    x_vals = np.arange(a, b, 1)
    for index, (prop, func) in enumerate(props, 1):
        pylab.subplot(2, 2, index)
        pylab.plot(x_vals, [func(i) for i in x_vals], label=fr'${props_titles[prop]}(T)$')
        pylab.scatter(temperature_segment,
                      [my_props.get_property(current_element, x, prop) for x in temperature_segment],
                      c='black')
        pylab.xlabel(r'$T$', fontsize=10)
        pylab.ylabel(fr'${props_titles[prop]}$', fontsize=10)
        pylab.grid(True)
        pylab.legend(loc='best', fontsize=12)
    pylab.gcf().set_size_inches(width, height)
    pylab.gcf().suptitle(current_element, fontsize=20)
    pylab.show()


def log_coefficients(coefficients: list) -> str:
    table = tabulate((coefficients,),
                     tablefmt='rounded_outline',
                     stralign='center',
                     numalign='center')
    header = '\n'.join([line.center(table.find('\n') + 1)
                        for line in tabulate((('Polynomial coefficients MNK',),),
                                             tablefmt='heavy_grid',
                                             stralign='center',
                                             numalign='center').splitlines()])
    return header + '\n' + table


def log_properties(title: str, temps, theoretical, practical) -> str:
    columns = ['T', 'theoretical value', 'practical value']
    data = list(zip(temps, theoretical, practical))
    table = tabulate(data,
                     headers=columns,
                     tablefmt='rounded_outline',
                     stralign='center',
                     numalign='center')
    header = '\n'.join([line.center(table.find('\n') + 1)
                        for line in tabulate(((title,),),
                                             tablefmt='heavy_grid',
                                             stralign='center',
                                             numalign='center').splitlines()])
    return header + '\n' + table


if __name__ == '__main__':
    function_array = [lambda x: 1, lambda x: math.log(x),
                      lambda x: 1 / (x ** 2), lambda x: 1 / x,
                      lambda x: x, lambda x: x ** 2, lambda x: x ** 3]
    my_props = TPoS()
    # ---------------------
    current_element = 'H'
    # ---------------------
    temperatures = my_props.temperatures_for_element(current_element)
    # temperatures = temperatures[temperatures.index(500):temperatures.index(6000) + 1]  # todo!
    polynomial_coefficients = mnk(list(map(lambda x: x * (1 / (10 ** 4)), temperatures)),
                                  [my_props.get_property(current_element, t, 'Fi') for t in temperatures],
                                  function_array)

    def fi(t: float) -> float:
        t_adapt = t * (1 / (10 ** 4))
        return (polynomial_coefficients[0] + polynomial_coefficients[1] * math.log(t_adapt)
                + polynomial_coefficients[2] * (1 / (t_adapt ** 2)) + polynomial_coefficients[3] * 1 / t_adapt
                + polynomial_coefficients[4] * t_adapt + polynomial_coefficients[5] * (t_adapt ** 2)
                + polynomial_coefficients[6] * (t_adapt ** 3))

    def g(t: float) -> float:
        return (my_props.get_constant(current_element, 'f_H_298.15')
                - my_props.get_property(current_element, 298.15, 'delta_H')
                - (fi(t) * t))

    def g_derivative(t: float) -> float:
        return (((-(2 * (math.pow(10, 8)) * polynomial_coefficients[2] / math.pow(t, 3))
                  - (math.pow(10, 4) * polynomial_coefficients[3] / math.pow(t, 2))
                  + (3 * polynomial_coefficients[6] * math.pow(t, 2) / math.pow(10, 12))
                  + (polynomial_coefficients[5] * t / (5 * math.pow(10, 7)))
                  + (polynomial_coefficients[1] / t)
                  + (polynomial_coefficients[4] / math.pow(10, 4)))
                 * (-t))
                - (polynomial_coefficients[6] * math.pow(t, 3) / math.pow(10, 12))
                - (polynomial_coefficients[5] * math.pow(t, 2) / math.pow(10, 8))
                - (polynomial_coefficients[4] * t / math.pow(10, 4))
                - (polynomial_coefficients[3] * math.pow(10, 4) / t)
                - (polynomial_coefficients[2] * math.pow(10, 8) / math.pow(t, 2))
                - (polynomial_coefficients[1] * math.log(t / math.pow(10, 4)))
                - polynomial_coefficients[0])

    def g_second_derivative(t: float) -> float:
        return ((-1 * math.pow(10, -12) * math.pow(t, -3) *
                 (4 * polynomial_coefficients[2] * math.pow(10, 20)
                  + polynomial_coefficients[3] * math.pow(10, 16) * t
                  + polynomial_coefficients[4] * math.pow(10, 8) * math.pow(t, 3)
                  + 4 * polynomial_coefficients[5] * math.pow(10, 4) * math.pow(t, 4)
                  + 9 * polynomial_coefficients[6] * math.pow(t, 5)))
                - (math.pow(10, -12) * math.pow(t, -3) *
                   (polynomial_coefficients[1] * math.pow(10, 12) * math.pow(t, 2)
                    - 2 * polynomial_coefficients[2] * math.pow(10, 20)
                    - polynomial_coefficients[3] * math.pow(10, 16) * t
                    + polynomial_coefficients[4] * math.pow(10, 8) * math.pow(t, 3)
                    + 2 * polynomial_coefficients[5] * math.pow(10, 4) * math.pow(t, 4)
                    + 3 * polynomial_coefficients[6] * math.pow(t, 5))))

    def s(t: float) -> float:
        return -g_derivative(t)

    def c_p(t: float) -> float:
        return -t * g_second_derivative(t)

    def delta_h(t: float) -> float:
        return (g(t) + t * s(t)) * math.pow(10, -3)


    print(log_coefficients(polynomial_coefficients))
    print(log_properties('Ф', temperatures,
                         [my_props.get_property(current_element, t, 'Fi') for t in temperatures],
                         [fi(t) for t in temperatures]))
    print(log_properties('S', temperatures,
                         [my_props.get_property(current_element, t, 'S') for t in temperatures],
                         [s(t) for t in temperatures]))
    print(log_properties('Cp', temperatures,
                         [my_props.get_property(current_element, t, 'Cp') for t in temperatures],
                         [c_p(t) for t in temperatures]))
    print(log_properties('ΔH', temperatures,
                         [my_props.get_property(current_element, t, 'delta_H') for t in temperatures],
                         [delta_h(t) for t in temperatures]))

    make_plot(100, 10000, [
        ('Fi', fi), ('S', s), ('Cp', c_p), ('delta_H', delta_h)
    ], 16, 8)
