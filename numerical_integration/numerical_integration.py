import math
import sys
from scipy import integrate


def get_by_step(h: float, x_0: float, x_k: float) -> list[float]:
    res = []
    while x_0 <= x_k:
        res.append(x_0)
        x_0 += h
    return res


def function(x: float) -> float:
    return x / ((2 * x + 7) * (3 * x + 4))


def rectangle_method(y_vals: list[float], h: float) -> float:
    return sum(y_vals) * h


def trapezoid_method(y_vals: list[float], h: float) -> float:
    n = len(y_vals)
    if n < 2:
        return 0.
    return sum(y_vals[i] / 2 if (i == 0 or i == n - 1) else y_vals[i] for i in range(n)) * h


def simpsons_method(y_vals: list[float], h: float) -> float:
    n = len(y_vals)
    res = 0.
    if n < 2:
        return res
    for i in range(n):
        if i == 0 or i == n - 1:
            res += y_vals[i]
        else:
            if i % 2 == 0:
                res += y_vals[i] * 2
            else:
                res += y_vals[i] * 4
    return res * h / 3


def handle_integration(h: float, x_0: float, x_k: float, f: function, method='rectangle', rounding=5) -> list[float]:
    x_vals = get_by_step(h, x_0, x_k)
    n = len(x_vals)
    match method:
        case 'rectangle':
            y_vals = [f((x_vals[i] + x_vals[i + 1]) / 2) for i in range(n - 1)]
        case 'trapezoid' | 'simpson':
            if method == 'simpsons' and n % 2 == 0:
                raise ValueError('Error: N is not even!')
            y_vals = [f((x_vals[i])) for i in range(n)]
        case _:
            raise ValueError('Error: invalid method name!')
    res = []
    match method:
        case 'rectangle':
            for i in range(n):
                res.append(round(rectangle_method(y_vals[:i], h), rounding))
        case 'trapezoid':
            for i in range(1, n + 1):
                res.append(round(trapezoid_method(y_vals[:i], h), rounding))
        case 'simpson':
            for i in range(1, n + 1):
                if i % 2 != 0:
                    res.append(round(simpsons_method(y_vals[:i], h), rounding))
                else:
                    res.append(None)
    return res


def runge_romberg_richardson(f_h: float, f_kh: float, k: float, p=2, rounding=5) -> float:
    return round(f_h + ((f_h - f_kh) / (math.pow(k, p) - 1)), rounding)


if __name__ == '__main__':
    with open(sys.argv[1], 'r', encoding='UTF-8') as fin:
        x0, xk, h1, h2 = [float(val) for val in fin.readline().rstrip('\n').split()]
    func = function

    rectangle_h1 = handle_integration(h1, x0, xk, func, method='rectangle')
    trapezoid_h1 = handle_integration(h1, x0, xk, func, method='trapezoid')
    simpson_h1 = handle_integration(h1, x0, xk, func, method='simpson')

    rectangle_h2 = handle_integration(h2, x0, xk, func, method='rectangle')
    trapezoid_h2 = handle_integration(h2, x0, xk, func, method='trapezoid')
    simpson_h2 = handle_integration(h2, x0, xk, func, method='simpson')

    rectangle_r_r_r = runge_romberg_richardson(rectangle_h2[-1], rectangle_h1[-1], h1 / h2)
    trapezoid_r_r_r = runge_romberg_richardson(trapezoid_h2[-1], trapezoid_h1[-1], h1 / h2)
    simpson_r_r_r = runge_romberg_richardson(simpson_h2[-1], simpson_h1[-1], h1 / h2)
    exact_value = round(integrate.quad(func, x0, xk)[0], 5)

    x_val_h1 = get_by_step(h1, x0, xk)
    x_val_h2 = get_by_step(h2, x0, xk)
    y_val_h1 = [round(func(x), 5) for x in x_val_h1]
    y_val_h2 = [round(func(x), 5) for x in x_val_h2]
    count_h1 = len(x_val_h1)
    count_h2 = len(x_val_h2)

    width = 5 + 15 * 5 + 7
    width_2 = 0 + 15 * 4 + 5
    header = ('|' + '|'.join(["i".center(5, " "),
                              "x".center(15, " "),
                              "y".center(15, " "),
                              "rectangle".center(15, " "),
                              "trapezoid".center(15, " "),
                              "simpson".center(15, " ")])
              + '|')
    header_2 = ('|' + '|'.join(["exact value".center(15, " "),
                                "rectangle".center(15, " "),
                                "trapezoid".center(15, " "),
                                "simpson".center(15, " ")])
                + '|')

    output_h1 = ['-' * width, header, '-' * width]
    for index, x, y, rect, trap, simp in zip(range(count_h1), x_val_h1, y_val_h1,
                                             rectangle_h1, trapezoid_h1, simpson_h1):
        line = ('|' + '|'.join([str(index).center(5, " "),
                                str(x).center(15, " "),
                                str(y).center(15, " "),
                                str(rect).center(15, " "),
                                str(trap).center(15, " "),
                                str(simp).center(15, " ")
                                if simp is not None
                                else ''.center(15, " ")])
                + '|')
        output_h1.append(line)
    output_h1.append('-' * width)
    output_str_h1 = '\n'.join(output_h1)

    output_h2 = ['-' * width, header, '-' * width]
    for index, x, y, rect, trap, simp in zip(range(count_h2), x_val_h2, y_val_h2,
                                             rectangle_h2, trapezoid_h2, simpson_h2):
        line = ('|' + '|'.join([str(index).center(5, " "),
                                str(x).center(15, " "),
                                str(y).center(15, " "),
                                str(rect).center(15, " "),
                                str(trap).center(15, " "),
                                str(simp).center(15, " ")
                                if simp is not None
                                else ''.center(15, " ")])
                + '|')
        output_h2.append(line)
    output_h2.append('-' * width)
    output_str_h2 = '\n'.join(output_h2)

    output_r_r_r = ['-' * width_2, header_2, '-' * width_2, '|'
                    + '|'.join([str(exact_value).center(15, " "),
                                str(rectangle_r_r_r).center(15, " "),
                                str(trapezoid_r_r_r).center(15, " "),
                                str(simpson_r_r_r).center(15, " ")])
                    + '|',
                    '|'
                    + '|'.join([''.center(15, " "),
                                str(round(abs(rectangle_r_r_r - exact_value), 5)).center(15, " "),
                                str(round(abs(trapezoid_r_r_r - exact_value), 5)).center(15, " "),
                                str(round(abs(simpson_r_r_r - exact_value), 5)).center(15, " ")])
                    + '|',
                    '-' * width_2]
    output_str_r_r_r = '\n'.join(output_r_r_r)

    all_table_str = output_str_h1 + '\n\n' + output_str_h2 + '\n\n' + output_str_r_r_r
    with open(sys.argv[2], 'w', encoding='UTF-8') as f_out:
        f_out.write(all_table_str)
        print(all_table_str)
