import math
import sys

import my_linear_algebra as lin_al
import tabulate as tbl
import itertools as it


# region -------- test functions --------
# def function_1(x_1: float, x_2: float) -> float:
#     return 0.1 * math.pow(x_1, 2) + x_1 + 0.2 * math.pow(x_2, 2) - 0.3
#
#
# def function_2(x_1: float, x_2: float) -> float:
#     return 0.2 * math.pow(x_1, 2) + x_2 - 0.1 * x_1 * x_2 - 0.7
#
#
# def function_1_dx1(x_1: float, x_2: float) -> float:
#     return 0.2 * x_1 + 1 + (x_2 * 0)
#
#
# def function_1_dx2(x_1: float, x_2: float) -> float:
#     return 0.4 * x_2 + (x_1 * 0)
#
#
# def function_2_dx1(x_1: float, x_2: float) -> float:
#     return 0.4 * x_1 - 0.1 * x_2
#
#
# def function_2_dx2(x_1: float, x_2: float) -> float:
#     return 1 - 0.1 * x_1 + (x_2 * 0)
#
#
# def phi_1(x_1: float, x_2: float) -> float:
#     return 0.3 - 0.1 * math.pow(x_1, 2) - 0.2 * math.pow(x_2, 2)
#
#
# def phi_2(x_1: float, x_2: float) -> float:
#     return 0.7 - 0.2 * math.pow(x_1, 2) + 0.1 * x_1 * x_2
#
#
# def phi_1_dx1(x_1: float, x_2: float) -> float:
#     return -0.2 * x_1 + (x_2 * 0)
#
#
# def phi_1_dx2(x_1: float, x_2: float) -> float:
#     return -0.4 * x_2 + (x_1 * 0)
#
#
# def phi_2_dx1(x_1: float, x_2: float) -> float:
#     return -0.4 * x_1 + 0.1 * x_2
#
#
# def phi_2_dx2(x_1: float, x_2: float) -> float:
#     return 0.1 * x_1 + (x_2 * 0)
# endregion -------- test functions --------


# region -------- variant math functions --------
def function_1(x_1: float, x_2: float) -> float:
    return x_1 - math.cos(x_2) - 1


def function_2(x_1: float, x_2: float) -> float:
    return x_2 - math.log(x_1 + 1, 10) - 3


def function_1_dx1(x_1: float, x_2: float) -> float:
    return 1 + (x_1 * 0) + (x_2 * 0)


def function_1_dx2(x_1: float, x_2: float) -> float:
    return math.sin(x_2) + (x_1 * 0)


def function_2_dx1(x_1: float, x_2: float) -> float:
    return -1 / (math.log(10) * x_1 + math.log(10)) + (x_2 * 0)


def function_2_dx2(x_1: float, x_2: float) -> float:
    return 1 + (x_1 * 0) + (x_2 * 0)


def phi_1(x_1: float, x_2: float) -> float:
    return math.cos(x_2) + 1 + (x_1 * 0)


def phi_2(x_1: float, x_2: float) -> float:
    return math.log(x_1 + 1, 10) + 3 + (x_2 * 0)


def phi_1_dx1(x_1: float, x_2: float) -> float:
    return 0 * x_1 * x_2


def phi_1_dx2(x_1: float, x_2: float) -> float:
    return -math.sin(x_2) + (x_1 * 0)


def phi_2_dx1(x_1: float, x_2: float) -> float:
    return 1 / (math.log(10) * x_1 + math.log(10)) + (x_2 * 0)


def phi_2_dx2(x_1: float, x_2: float) -> float:
    return 0 * x_1 * x_2
# endregion -------- variant math functions --------


def make_j(xk: tuple[float, float]) -> list[list[float]]:
    return [
        [function_1_dx1(*xk), function_1_dx2(*xk)],
        [function_2_dx1(*xk), function_2_dx2(*xk)]
    ]


def make_a1(xk: tuple[float, float]) -> list[list[float]]:
    return [
        [function_1(*xk), function_1_dx2(*xk)],
        [function_2(*xk), function_2_dx2(*xk)]
    ]


def make_a2(xk: tuple[float, float]) -> list[list[float]]:
    return [
        [function_1_dx1(*xk), function_1(*xk)],
        [function_2_dx1(*xk), function_2(*xk)]
    ]


def newton_method(a_1: float, b_1: float, a_2: float, b_2: float, eps: float) \
        -> list[list[tuple[float, float]]]:
    x0_1 = (b_1 + a_1) / 2
    x0_2 = (b_2 + a_2) / 2
    x_list = [(x0_1, x0_2)]
    func_vals_list = []
    func_dx1_vals_list = []
    func_dx2_vals_list = []
    det_a1_list = []
    det_a2_list = []
    det_j_list = []
    while len(x_list) < 2 or lin_al.vector_norm(lin_al.subtract_vectors(list(x_list[-1]), list(x_list[-2]))) > eps:
        func_vals_list.append((function_1(*x_list[-1]), function_2(*x_list[-1])))
        func_dx1_vals_list.append((function_1_dx1(*x_list[-1]), function_2_dx1(*x_list[-1])))
        func_dx2_vals_list.append((function_1_dx2(*x_list[-1]), function_2_dx2(*x_list[-1])))
        det_a1 = lin_al.determinant(make_a1(x_list[-1]))
        det_a2 = lin_al.determinant(make_a2(x_list[-1]))
        det_j = lin_al.determinant(make_j(x_list[-1]))
        if round(det_j, -int(math.log(eps, 10)) + 1) == 0:
            raise ValueError('Newton method ERROR --> det(J) is 0!')
        det_a1_list.append(det_a1)
        det_a2_list.append(det_a2)
        det_j_list.append(det_j)
        x_list.append((
            x_list[-1][0] - (det_a1 / det_j),
            x_list[-1][1] - (det_a2 / det_j)
        ))
    return [x_list, func_vals_list, func_dx1_vals_list, func_dx2_vals_list,
            det_a1_list, det_a2_list, det_j_list]


def get_q(a_1: float, b_1: float, a_2: float, b_2: float):
    phi_derivative = [
        [max(abs(phi_1_dx1(a_1, a_2)),
             abs(phi_1_dx1(b_1, a_2)),
             abs(phi_1_dx1(a_1, b_2)),
             abs(phi_1_dx1(b_1, b_2))),
         max(abs(phi_1_dx2(a_1, a_2)),
             abs(phi_1_dx2(b_1, a_2)),
             abs(phi_1_dx2(a_1, b_2)),
             abs(phi_1_dx2(b_1, b_2)))],
        [max(abs(phi_2_dx1(a_1, a_2)),
             abs(phi_2_dx1(b_1, a_2)),
             abs(phi_2_dx1(a_1, b_2)),
             abs(phi_2_dx1(b_1, b_2))),
         max(abs(phi_2_dx2(a_1, a_2)),
             abs(phi_2_dx2(b_1, a_2)),
             abs(phi_2_dx2(a_1, b_2)),
             abs(phi_2_dx2(b_1, b_2)))],
    ]
    return max(sum(t) for t in phi_derivative)


def check_phi(borders, g):
    for c in it.combinations(borders, r=2):
        if abs(phi_1(*c)) >= g or abs(phi_2(*c)) >= g:
            return False
    return True


def simple_iteration_method(a_1: float, b_1: float, a_2: float, b_2: float, eps: float) \
        -> list[list[tuple[float, float]]]:
    g = (((b_1 - a_1) / 4) + ((b_2 - a_2) / 4)) / 2
    x0_1 = (b_1 + a_1) / 2
    x0_2 = (b_2 + a_2) / 2
    q = get_q(a_1, b_1, a_2, b_2)
    if q >= 1:
        raise ValueError(f'Simple iteration method ERROR --> q is greater than 1 ({q})!')
    if check_phi([a_1, b_1, a_2, b_2], g):
        raise ValueError(f'ERROR!')
    else:
        print('its ok')
    x_list = [(x0_1, x0_2)]
    phi_list = []
    while (len(x_list) < 2 or
           lin_al.vector_norm(lin_al.subtract_vectors(list(x_list[-1]), list(x_list[-2]))) * (q / (1 - q)) > eps):
        phi_k_1 = phi_1(*x_list[-1])
        phi_k_2 = phi_2(*x_list[-1])
        if abs(phi_k_1 - x0_1) > g or abs(phi_k_2 - x0_2) > g:
            raise ValueError('Simple iteration method ERROR --> value not in G area!')
        x_list.append((phi_k_1, phi_k_2))
        phi_list.append((phi_k_1, phi_k_2))
    return [x_list, phi_list, q]


def round_up(arr: list[tuple[float]], eps: float) -> list[tuple[float]]:
    res = []
    for el in arr:
        if isinstance(el, tuple):
            res.append((round(el[0], -int(math.log(eps, 10)) + 2), round(el[1], -int(math.log(eps, 10)) + 2)))
        elif isinstance(el, float):
            res.append(round(el, -int(math.log(eps, 10)) + 2))
        else:
            res.append(el)
    return res


if __name__ == '__main__':
    with open(sys.argv[1], 'r', encoding='UTF-8') as fin:
        epsilon = float(fin.readline().rstrip('\n'))
        a1, b1 = [float(el) for el in fin.readline().rstrip('\n').split()]
        a2, b2 = [float(el) for el in fin.readline().rstrip('\n').split()]

    try:
        (x_vals, func_vals, func_dx1_vals, func_dx2_vals,
         det_a1_vals, det_a2_vals, det_j_vals) = newton_method(a1, b1, a2, b2, epsilon)

        header_newton = ['k', 'x1, x2', 'f1, f2', 'df1/dx1, df2/dx1',
                         'df1/dx2, df2/dx2', 'det A1', 'det A2', 'det J']
        k_newton = [*range(len(x_vals))]
        table_newton = tbl.tabulate(
            map(lambda x: round_up(x, epsilon),
                it.zip_longest(
                    k_newton, x_vals, func_vals, func_dx1_vals, func_dx2_vals,
                    det_a1_vals, det_a2_vals, det_j_vals
                )), headers=header_newton, tablefmt='rounded_outline', stralign='center')
        title_newton = f'   Newton method (epsilon = {epsilon})'
        result_newton = (f'   x1 = {round(x_vals[-1][0], -int(math.log(epsilon, 10)) + 1)}, '
                         f'x2 = {round(x_vals[-1][1], -int(math.log(epsilon, 10)) + 1)}')
        all_strings_newton = '\n'.join([title_newton, table_newton, result_newton])
    except ValueError as ex:
        all_strings_newton = ex.args[0]

    try:
        x_vals_simple_it, phi_vals, q_val = simple_iteration_method(a1, b1, a2, b2, epsilon)

        header_simple_it = ['k', 'x1, x2', 'phi1, phi2']
        k_simple_it = [*range(len(x_vals_simple_it))]
        table_simple_it = tbl.tabulate(
            map(lambda x: round_up(x, epsilon),
                it.zip_longest(
                    k_simple_it, x_vals_simple_it, phi_vals
                )), headers=header_simple_it, tablefmt='rounded_outline', stralign='center')
        title_simple_it = f'   Simple iteration method (epsilon = {epsilon})'
        result_simple_it = (f'   x1 = {round(x_vals_simple_it[-1][0], -int(math.log(epsilon, 10)) + 1)}, '
                            f'x2 = {round(x_vals_simple_it[-1][1], -int(math.log(epsilon, 10)) + 1)}\n'
                            f'   q = {round(q_val, -int(math.log(epsilon, 10)) + 1)}')
        all_strings_simple_it = '\n'.join([title_simple_it, table_simple_it, result_simple_it])
    except ValueError as ex:
        all_strings_simple_it = ex.args[0]

    all_strings = all_strings_newton + '\n\n\n' + all_strings_simple_it
    with open(sys.argv[2], 'w', encoding='UTF-8') as f_out:
        f_out.write(all_strings)
        print()
        print(all_strings)
