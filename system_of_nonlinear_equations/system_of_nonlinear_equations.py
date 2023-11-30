import math
import sys

import my_linear_algebra as lin_al
import tabulate as tbl
import itertools as it


# region test functions
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

# endregion


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


def newton_method(x0_1: float, x0_2: float, eps: float) \
        -> list[list[tuple[float, float]]]:
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
        x1_0, x2_0 = [float(el) for el in fin.readline().rstrip('\n').split()]

    (x_vals, func_vals, func_dx1_vals, func_dx2_vals,
     det_a1_vals, det_a2_vals, det_j_vals) = newton_method(x1_0, x2_0, epsilon)

    header_newton = ['k', 'x1, x2', 'f1, f2', 'df1/dx1, df2/dx1',
                     'df1/dx2, df2/dx2', 'det A1', 'det A2', 'det J']
    k_newton = [*range(len(x_vals))]
    table_newton = tbl.tabulate(
        map(lambda x: round_up(x, epsilon),
            it.zip_longest(
                k_newton, x_vals, func_vals, func_dx1_vals, func_dx2_vals,
                det_a1_vals, det_a2_vals, det_j_vals
            )), headers=header_newton, tablefmt='rounded_outline', stralign='center')
    width_newton = table_newton.find('\n') + 1
    title_newton = f'Newton method (epsilon = {epsilon})'.center(width_newton, ' ')
    result_newton = (f'x1 = {round(x_vals[-1][0], -int(math.log(epsilon, 10)) + 1)}, '
                     f'x2 = {round(x_vals[-1][1], -int(math.log(epsilon, 10)) + 1)}'
                     .center(width_newton, ' '))

    with open(sys.argv[2], 'w', encoding='UTF-8') as f_out:
        f_out.write('\n'.join([title_newton, table_newton, result_newton]))
        print(title_newton)
        print(table_newton)
        print(result_newton)
