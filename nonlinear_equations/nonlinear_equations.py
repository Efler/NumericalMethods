import math
import sys

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


def check_conditions_simple_iteration(a: float, b: float, phi_func, phi_der_func) \
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


def simple_iteration_method(a: float, b: float, phi_func, phi_der_func, eps: float) \
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


def create_table_newton(a: float, b: float, eps: float) -> str:
    table_newton = []
    width_newton = 6 + 5 + 15 * 3 + 17
    line_newton = '-' * width_newton
    title_newton = ('|' + 'Newton method'.center(width_newton - 2, ' ')
                    + '|\n|'
                    + f'[{a}, {b}], epsilon = {eps}'.center(width_newton - 2, ' ')
                    + '|')
    header_newton = '|' + '|'.join([
        'k'.center(5, ' '),
        'x'.center(15, ' '),
        'f(x)'.center(15, ' '),
        "f'(x)".center(15, ' '),
        "-f(x) / f'(x)".center(17, ' '),
    ]) + '|'
    try:
        x_list, f_list, f_der_list, diff_list = newton_method(a, b, eps)
        table_newton.append(line_newton)
        table_newton.append(title_newton)
        table_newton.append(line_newton)
        table_newton.append(header_newton)
        table_newton.append(line_newton)
        for j in range(len(x_list) - 1):
            tmp = '|' + '|'.join([
                f'{j}'.center(5, ' '),
                f'{round(x_list[j], -round(math.log(eps, 10)) + 1)}'
                .center(15, ' '),
                f'{round(f_list[j], -round(math.log(eps, 10)) + 1)}'
                .center(15, ' '),
                f'{round(f_der_list[j], -round(math.log(eps, 10)) + 1)}'
                .center(15, ' '),
                f'{round(diff_list[j], -round(math.log(eps, 10)) + 1)}'
                .center(17, ' '),
            ]) + '|'
            table_newton.append(tmp)
        table_newton.append('|' + '|'.join([
            f'{len(x_list) - 1}'.center(5, ' '),
            f'[{round(x_list[len(x_list) - 1], -round(math.log(eps, 10)) + 1)}]'
            .center(15, ' '),
            f''.center(15, ' '),
            f''.center(15, ' '),
            f''.center(17, ' '),
        ]) + '|')
        table_newton.append(line_newton)
        table_newton.append('|' + f'x = {round(x_list[len(x_list) - 1], -round(math.log(eps, 10)))}'
                            .center(width_newton - 2) + '|')
        table_newton.append(line_newton)
        table_newton_str = '\n'.join(table_newton)
    except ValueError as ex:
        table_newton_str = ex.args[0]
    return table_newton_str


def create_table_simple_iteration(a: float, b: float, eps: float, phi_func, phi_der_func) -> str:
    table_simple_iteration = []
    width_simple_iteration = 4 + 15 * 3
    line_simple_iteration = '-' * width_simple_iteration
    title_simple_iteration = ('|' + 'Simple iteration method'.center(width_simple_iteration - 2, ' ')
                              + '|\n|'
                              + f'[{a}, {b}], epsilon = {eps}'.center(width_simple_iteration - 2, ' ')
                              + '|')
    header_simple_iteration = '|' + '|'.join([
        'k'.center(15, ' '),
        'x'.center(15, ' '),
        'phi(x)'.center(15, ' ')
    ]) + '|'
    try:
        x_list_2, phi_list_2 = simple_iteration_method(a, b, phi_func, phi_der_func, eps)
        table_simple_iteration.append(line_simple_iteration)
        table_simple_iteration.append(title_simple_iteration)
        table_simple_iteration.append(line_simple_iteration)
        table_simple_iteration.append(header_simple_iteration)
        table_simple_iteration.append(line_simple_iteration)
        for j in range(len(x_list_2) - 1):
            tmp = '|' + '|'.join([
                f'{j}'.center(15, ' '),
                f'{round(x_list_2[j], -round(math.log(eps, 10)) + 1)}'
                .center(15, ' '),
                f'{round(phi_list_2[j], -round(math.log(eps, 10)) + 1)}'
                .center(15, ' ')
            ]) + '|'
            table_simple_iteration.append(tmp)
        table_simple_iteration.append('|' + '|'.join([
            f'{len(x_list_2) - 1}'.center(15, ' '),
            f'[{round(x_list_2[len(x_list_2) - 1], -round(math.log(eps, 10)) + 1)}]'
            .center(15, ' '),
            f''.center(15, ' ')
        ]) + '|')
        table_simple_iteration.append(line_simple_iteration)
        table_simple_iteration.append('|' + f'x = {round(x_list_2[len(x_list_2) - 1], -round(math.log(eps, 10)))}'
                                      .center(width_simple_iteration - 2) + '|')
        table_simple_iteration.append(line_simple_iteration)
        table_simple_iteration_str = '\n'.join(table_simple_iteration)
    except ValueError as ex:
        table_simple_iteration_str = ex.args[0]
    return table_simple_iteration_str


if __name__ == '__main__':
    with open(sys.argv[1], 'r', encoding='UTF-8') as fin:
        epsilon = float(fin.readline().rstrip('\n'))
        a1, b1 = [float(el) for el in fin.readline().rstrip('\n').split()]
        a2, b2 = [float(el) for el in fin.readline().rstrip('\n').split()]

    # [a1, b1] Newton
    table_newton_string_1 = create_table_newton(a1, b1, epsilon)
    # [a2, b2] Newton
    table_newton_string_2 = create_table_newton(a2, b2, epsilon)
    # [a1, b1] Simple iteration
    table_simple_iteration_string_1 = create_table_simple_iteration(a1, b1, epsilon, phi_1, phi_derivative_1)
    # [a2, b2] Simple iteration
    table_simple_iteration_string_2 = create_table_simple_iteration(a2, b2, epsilon, phi_2, phi_derivative_2)

    with open(sys.argv[2], 'w', encoding='UTF-8') as f_out:
        all_tables = '\n\n'.join([
            table_newton_string_1,
            table_newton_string_2,
            table_simple_iteration_string_1,
            table_simple_iteration_string_2
        ])
        f_out.write(all_tables)
    print(all_tables)
