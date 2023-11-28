import math

import numpy as np


# region test functions
# def function(x: float) -> float:
#     return math.exp(2 * x) + 3 * x - 4
#
#
# def first_derivative(x: float) -> float:
#     return 2 * math.exp(2 * x) + 3
#
#
# def second_derivative(x: float) -> float:
#     return 4 * math.exp(2 * x)
#
#
# def phi(x: float) -> float:
#     return math.log(4 - 3 * x) / 2
#
#
# def phi_derivative(x: float) -> float:
#     return -3. / (2 * (4 - 3 * x))
# endregion


def function(x: float) -> float:
    return math.exp(x) - 2 * x - 2


def first_derivative(x: float) -> float:
    return math.exp(x) - 2


def second_derivative(x: float) -> float:
    return math.exp(x)


def phi_1(x: float) -> float:
    return math.log(2 * x + 2)


def phi_derivative_1(x: float) -> float:
    return 2 / (2 * x + 2)


def phi_2(x: float) -> float:
    return (math.exp(x) - 2) / 2


def phi_derivative_2(x: float) -> float:
    return math.exp(x) / 2


def check_conditions_newton(a: float, b: float, x_0: float) -> bool:
    x_vals = np.arange(a, b + 0.01, 0.01)
    first_derivative_sign = False
    for i in x_vals:
        if i == a:
            if first_derivative(i) < 0:
                first_derivative_sign = False
            elif first_derivative(i) > 0:
                first_derivative_sign = True
            else:
                return False
        else:
            if first_derivative(i) == 0:
                return False
            if ((first_derivative(i) < 0 and first_derivative_sign) or
                    (first_derivative(i) > 0 and not first_derivative_sign)):
                return False
    second_derivative_sign = False
    for i in x_vals:
        if i == a:
            if second_derivative(i) < 0:
                second_derivative_sign = False
            elif second_derivative(i) > 0:
                second_derivative_sign = True
            else:
                return False
        else:
            if second_derivative(i) == 0:
                return False
            if ((second_derivative(i) < 0 and second_derivative_sign) or
                    (second_derivative(i) > 0 and not second_derivative_sign)):
                return False
    if function(x_0) * second_derivative(x_0) <= 0:
        return False
    return True


def newton_method(a: float, b: float, eps: float) \
        -> tuple[list[float], list[float], list[float], list[float]]:
    if check_conditions_newton(a, b, a):
        x_0 = a
    elif check_conditions_newton(a, b, b):
        x_0 = b
    else:
        raise ValueError('Error: Conditions are not met for Newton method!')
    x = [x_0]
    f = []
    f_der = []
    diff = []
    while len(x) == 1 or abs(x[-1] - x[-2]) >= eps:
        f_k = function(x[-1])
        f_der_k = first_derivative(x[-1])
        diff_k = -f_k / f_der_k
        x.append(x[-1] + diff_k)
        f.append(f_k)
        f_der.append(f_der_k)
        diff.append(diff_k)
    return x, f, f_der, diff


def check_conditions_simple_iteration(a: float, b: float, phi_func: function, phi_der_func: function) \
        -> float:
    x_vals = np.arange(a, b + 0.01, 0.01)
    try:
        for x in x_vals:
            phi_func(x)
    except ValueError:
        return -1.
    q = float('-inf')
    for x in x_vals:
        q = max(q, abs(phi_der_func(x)))
    if q >= 1:
        return -1.
    return q


def simple_iteration_method(a: float, b: float, phi_func: function, phi_der_func: function, eps: float) \
        -> tuple[list[float], list[float]]:
    if (q := check_conditions_simple_iteration(a, b, phi_func, phi_der_func)) < 0:
        raise ValueError('Error conditions are not met for Simple iteration method!')
    x = [(a + b) / 2]
    phi_vals = []
    while len(x) == 1 or ((q / (1 - q)) * abs(x[-1] - x[-2])) > eps:
        phi_k = phi_func(x[-1])
        x.append(phi_k)
        phi_vals.append(phi_k)
    return x, phi_vals


if __name__ == '__main__':
    # TODO !!!!!!!!!!!!!!!!!! --------------------------------
    epsilon = 0.001
    try:
        for values in newton_method(1.5, 2.0, epsilon):
            print([round(el, -round(math.log(epsilon, 10)) + 1) for el in values])
        print()
        for values in newton_method(-1.0, -0.5, epsilon):
            print([round(el, -round(math.log(epsilon, 10)) + 1) for el in values])
    except ValueError as ex:
        print(ex.args[0])

    print('\n-----\n')

    try:
        for values in simple_iteration_method(1.5, 2.0, phi_1, phi_derivative_1, epsilon):
            print([round(el, -round(math.log(epsilon, 10)) + 1) for el in values])
        print()
        for values in simple_iteration_method(-1.0, -0.5, phi_2, phi_derivative_2, epsilon):
            print([round(el, -round(math.log(epsilon, 10)) + 1) for el in values])
    except ValueError as ex:
        print(ex.args[0])
    # TODO !!!!!!!!!!!!!!!!!! --------------------------------
