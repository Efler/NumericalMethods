import sys
from copy import deepcopy
from math import sqrt, pow, log


def beautify_matrix(matrix: [[float]]) -> [[str]]:
    str_matrix = []
    res_matrix = []
    for row in matrix:
        for el in row:
            str_matrix.append(str(el))
    max_width = max(map(lambda x: len(x), str_matrix))
    for row in matrix:
        res_matrix.append([str(x).rjust(max_width, ' ') for x in row])
    return res_matrix


def glue_str_matrix(matrix: [[str]], row_sep='\n', col_sep='  ') -> str:
    return row_sep.join(col_sep.join(row) for row in matrix)


def round_up(matrix: [[float]], n: int = 2) -> [[float]]:
    res = deepcopy(matrix)
    for i in range(len(res)):
        for j in range(len(res[i])):
            res[i][j] = round(res[i][j], n)
    return res


def multiply_matrix(a: [[float]], b: [[float]], eps: int = None) -> [[float]]:
    n, m = len(a), len(a[0])
    n2, m2 = len(b), len(b[0])
    if m != n2:
        raise TypeError('Invalid sizes of matrices to multiply')
    result = [[0] * n for _ in range(m2)]
    for r in range(n):
        for j in range(m2):
            for k in range(n2):
                result[r][j] += a[r][k] * b[k][j]
            result[r][j] = result[r][j]
            if eps is not None:
                result[r][j] = round(result[r][j], 2)
    return result


def multiply_vectors_to_matrix(a: [[float]], b: [float]) -> [[float]]:
    if len(a) != len(b):
        raise TypeError('Invalid sizes of vectors to multiply')
    res = []
    for k in range(len(a)):
        row = []
        for j in range(len(b)):
            row.append(a[k][0] * b[j])
        res.append(row)
    return res


def multiply_vectors_to_value(a: [float], b: [[float]]) -> float:
    if len(a) != len(b):
        raise TypeError('Invalid sizes of vectors to multiply')
    return sum(x * y for x, y in zip(a, [row[0] for row in b]))


def transpose(vector: list[float]) -> list[list[float]]:
    return [[v] for v in vector]


def eye(n: int) -> [[float]]:
    res = []
    for k in range(n):
        res.append([0.0] * n)
    for k in range(n):
        res[k][k] = 1.0
    return res


def sign(x: float) -> int:
    if x < 0:
        return -1
    elif x == 0:
        return 0
    else:
        return 1


def vector_euclidean_norm(vector: list[float], k=0) -> float:
    return sqrt(sum(pow(v, 2) for v in vector[k:]))


def get_v_vector(a_vector: [float], k=0) -> list[float]:
    v_vector = [0] * k
    for n in range(k, len(a_vector)):
        if n == k:
            v_vector.append(a_vector[n] + sign(a_vector[n]) * vector_euclidean_norm(a_vector, k))
        else:
            v_vector.append(a_vector[n])
    return v_vector


def householder_matrix(v_vector: list[float], eps: int = None) -> list[list[float]]:
    cf = 2 / multiply_vectors_to_value(v_vector, transpose(v_vector))
    vv_matrix = multiply_vectors_to_matrix(transpose(v_vector), v_vector)
    eye_matrix = eye(len(vv_matrix))
    for j in range(len(vv_matrix)):
        for k in range(len(vv_matrix)):
            eye_matrix[j][k] -= vv_matrix[j][k] * cf
            if eps is not None:
                eye_matrix[j][k] = round(eye_matrix[j][k], eps)
    return eye_matrix


def decompose(matrix: [[float]], eps: int = None) -> tuple[[[float]], [[float]]]:
    if len(matrix) != len(matrix[0]):
        raise TypeError('Invalid matrix')
    a_matrix = deepcopy(matrix)
    h_matrices = []
    for k in range(len(matrix) - 1):
        h_k = householder_matrix(get_v_vector([row[k] for row in a_matrix], k))
        h_matrices.append(h_k)
        a_matrix = multiply_matrix(h_k, a_matrix, eps if k == len(matrix) - 2 else None)
    q_matrix = h_matrices[0]
    for k in range(1, len(h_matrices)):
        q_matrix = multiply_matrix(q_matrix, h_matrices[k],
                                   eps if k == len(h_matrices) - 1 else None)
    return q_matrix, a_matrix


def complex_abs(x: float, y: float) -> float:
    return sqrt(pow(x, 2) + pow(y, 2))


def convert_epsilon(eps: float | None) -> int | None:
    if eps is None:
        return None
    res = -round(log(eps, 10))
    return res if res >= 0 else 0


def find_complex_conjugate_pairs(a11: int, a12: int, a21: int, a22: int)\
        -> tuple[tuple[float, float] | None, tuple[float, float] | None]:
    b = -a11 - a22
    c = (a11 * a22) - (a12 * a21)
    d = pow(b, 2) - (4 * c)
    if d >= 0:
        return None, None
    else:
        return (-b / 2, sqrt(abs(d)) / 2), (-b / 2, -sqrt(abs(d)) / 2)


def check_matrix_blocks(matrix: [[float]], pairs: [[tuple[float, float]]], eps: float, first=False)\
        -> [tuple[int, [[tuple[float, float]]]]]:
    res = []
    for k in range(len(matrix) - 1):
        h1, h2 = find_complex_conjugate_pairs(matrix[k][k], matrix[k][k + 1],
                                              matrix[k + 1][k], matrix[k + 1][k + 1])
        if first:
            if h1 is not None and h2 is not None:
                pairs.append([h1, h2])
            else:
                pairs.append(None)
        else:
            if h1 is not None and h2 is not None:
                if pairs[k] is not None:
                    if (complex_abs(h1[0] - pairs[k][0][0], h1[1] - pairs[k][0][1]) <= eps
                       and complex_abs(h2[0] - pairs[k][1][0], h2[1] - pairs[k][1][1]) <= eps):
                        res.append((k, [h1, h2]))
                    pairs[k] = [h1, h2]
            else:
                pairs[k] = None
    res.sort(key=lambda x: x[0])
    return res


def check_eigenvalues(matrix: [[float]], errors: [float], eps: float, first=False) -> [tuple[int, float, float]]:
    res = []
    for k in range(len(matrix)):
        if first:
            if k != len(matrix) - 1:
                checker = vector_euclidean_norm([row[k] for row in matrix], k + 1)
                errors.append(checker)
                if checker <= eps and errors[k] is not None:
                    res.append((k, matrix[k][k], checker))
            else:
                res.append((k, matrix[k][k], None))
        else:
            if k != len(matrix) - 1:
                checker = vector_euclidean_norm([row[k] for row in matrix], k + 1)
                if errors[k] is not None and checker <= errors[k]:
                    errors[k] = checker
                else:
                    errors[k] = None
                if checker <= eps and errors[k] is not None:
                    res.append((k, matrix[k][k], checker))
            else:
                res.append((k, matrix[k][k], None))
    res.sort(key=lambda x: x[0])
    return res


if __name__ == '__main__':
    input_matrix = []
    for line in sys.stdin:
        line = line.rstrip('\n')
        if line == '':
            break
        else:
            input_matrix.append([float(x) for x in line.split()])
    epsilon = float(sys.stdin.readline().rstrip('\n'))
    try:
        header_logs = f'------------- Logs of iterations (epsilon = {epsilon}) -------------'
        line_logs = ''.ljust(len(header_logs), '#')
        header_results = f'------------------ RESULTS ------------------'
        line_results = ''.ljust(len(header_results), '#')
        logs = [line_logs, header_logs, line_logs, '']

        a_matrix = deepcopy(input_matrix)
        final_result_values = []
        complex_conjugate_pairs = []
        errors_under_diagonal = []
        iterations_count = 0
        iterations_flag = True
        while iterations_flag:
            header_iteration = f'Iteration {iterations_count}'
            header_a = f'A({iterations_count}) matrix'
            line_iteration = ''.ljust(len(header_iteration), '=')
            line_a = ''.ljust(len(header_a), '-')
            logs.append(line_iteration)
            logs.append(header_iteration)
            logs.append(line_iteration)
            logs.append('')
            logs.append(header_a)
            logs.append(line_a)
            logs.append(glue_str_matrix(
                            beautify_matrix(
                                round_up(a_matrix, convert_epsilon(epsilon) + 1))))

            q, r = decompose(a_matrix)
            a_matrix = multiply_matrix(r, q)
            eigenvalues_checker = check_eigenvalues(a_matrix, errors_under_diagonal, epsilon,
                                                    True if iterations_count == 0 else False)
            blocks_checker = check_matrix_blocks(a_matrix, complex_conjugate_pairs, epsilon,
                                                 True if iterations_count == 0 else False)
            results = [None] * len(a_matrix)
            for i in range(len(results)):
                if i != len(results) - 1:
                    tmp = None
                    for j in eigenvalues_checker:
                        if j[0] == i:
                            tmp = j[1]
                    if tmp is not None:
                        results[i] = tmp
                    else:
                        for j in blocks_checker:
                            if j[0] == i:
                                tmp = j[1][0]
                        if tmp is not None:
                            results[i] = tmp
                        else:
                            for j in blocks_checker:
                                if j[0] == i - 1:
                                    tmp = j[1][1]
                            if tmp is not None:
                                results[i] = tmp
                else:
                    tmp = None
                    for j in blocks_checker:
                        if j[0] == i - 1:
                            tmp = j[1][1]
                    if tmp is not None:
                        results[i] = tmp
                    else:
                        for j in eigenvalues_checker:
                            if j[0] == i:
                                tmp = j[1]
                        if tmp is not None:
                            results[i] = tmp
            none_flag = False
            for i in results:
                if i is None:
                    none_flag = True
            if not none_flag:
                iterations_flag = False
                final_result_values = deepcopy(results)
            iterations_count += 1

            header_q = f'Q({iterations_count - 1}) matrix'
            header_r = f'R({iterations_count - 1}) matrix'
            header_a_next = f'A({iterations_count}) = R({iterations_count - 1}) * Q({iterations_count - 1})'
            line_q = ''.ljust(len(header_q), '-')
            line_r = ''.ljust(len(header_r), '-')
            line_a_next = ''.ljust(len(header_a_next), '-')
            logs.append('')
            logs.append(header_q)
            logs.append(line_q)
            logs.append(glue_str_matrix(
                            beautify_matrix(
                                round_up(q, convert_epsilon(epsilon) + 1))))
            logs.append('')
            logs.append(header_r)
            logs.append(line_r)
            logs.append(glue_str_matrix(
                            beautify_matrix(
                                round_up(r, convert_epsilon(epsilon) + 1))))
            logs.append('')
            logs.append(header_a_next)
            logs.append(line_a_next)
            logs.append(glue_str_matrix(
                            beautify_matrix(
                                round_up(a_matrix, convert_epsilon(epsilon) + 1))))
            logs.append('')
            logs.append(line_results)

        logs.append('')
        logs.append(line_results)
        logs.append(header_results)
        logs.append(line_results)
        logs.append('')
        for index, result in enumerate(final_result_values, 1):
            if isinstance(result, float):
                str_res = f'λ{index} = {round(result, convert_epsilon(epsilon))}'
            elif isinstance(result, tuple):
                x, y = result
                if y < 0:
                    str_sign = '-'
                    y = abs(y)
                else:
                    str_sign = '+'
                str_res = (f'λ{index} = {round(x, convert_epsilon(epsilon))}'
                           f' {str_sign}'
                           f' {round(y, convert_epsilon(epsilon))}i')
            else:
                raise TypeError('Invalid result output', result)
            logs.append(str_res)
        logs.append('')
        logs.append(f'epsilon: {epsilon}')
        logs.append(f'total iterations: {iterations_count}')
        logs.append('')
        logs.append(line_results)
        logs_string = '\n'.join(logs)
        print(logs_string)

    except OSError as ex:
        print(f'OSError! Message: {"; ".join(map(str, ex.args))}')
    except TypeError as ex:
        print("; ".join(map(str, ex.args)))
