import math
import sys
from copy import deepcopy


def multiply_matrix(a: [[float]], b: [[float]]) -> [[float]]:
    n, m = len(a), len(a[0])
    n2, m2 = len(b), len(b[0])
    if m != n2:
        raise TypeError('Invalid sizes of matrices to multiply')
    result = [[0] * n for _ in range(m2)]
    for i in range(n):
        for j in range(m2):
            for k in range(n2):
                result[i][j] += a[i][k] * b[k][j]
            result[i][j] = result[i][j]
    return result


def multiply_matrix_vector(matrix: [[float]], vector: [float]) -> [float]:
    res = []
    for i in range(len(matrix)):
        res.append(sum(list(map(lambda x, y: x * y, matrix[i], vector))))
    return res


def transport_matrix(matrix: [[float]]) -> [[float]]:
    res_matrix = deepcopy(matrix)
    for i in range(len(matrix)):
        for j in range(i + 1, len(matrix)):
            res_matrix[i][j], res_matrix[j][i] = res_matrix[j][i], res_matrix[i][j]
    return res_matrix


def multiply_vector(a: [float], b: [float]) -> float:
    return sum(list(map(lambda x, y: x * y, a, b)))


def divide_vector(vector: [float], divider: float) -> [float]:
    return [i / divider for i in vector]


def vector_norm(vector: [float]) -> float:
    return max(vector)


def check_matrix(in_matrix: [[float]]) -> bool:
    n = len(in_matrix)
    for row in in_matrix:
        if len(row) != n:
            return False
    return True


def t(matrix: [[float]]) -> float:
    res = 0
    for i in range(len(matrix) - 1):
        for j in range(i + 1, len(matrix)):
            res += matrix[i][j] * matrix[i][j]
    return res ** 0.5


def max_non_diagonal(matrix: [[float]]) -> [int, int]:
    max_num = matrix[0][1]
    max_ij = [0, 1]
    for i in range(len(matrix)):
        for j in range(i + 1, len(matrix)):
            if abs(max_num) < abs(matrix[i][j]):
                max_ij = [i, j]
    return max_ij


def phi(matrix: [[float]], i: int, j: int) -> float:
    if matrix[i][i] == matrix[j][j]:
        return math.pi / 4
    return 0.5 * (math.atan((2 * matrix[i][j]) / (matrix[i][i] - matrix[j][j])))


def rotation_matrix(phi_num: float, size: int, m_i: int, m_j: int) -> [[float]]:
    res = [[1.0 if i == j else 0.0 for j in range(size)] for i in range(size)]
    res[m_i][m_j] = -math.sin(phi_num)
    res[m_j][m_i] = math.sin(phi_num)
    res[m_i][m_i] = math.cos(phi_num)
    res[m_j][m_j] = math.cos(phi_num)
    return res


def check_orthogonal(res_x: [[float]], eps: float) -> bool:
    if (multiply_vector(res_x[0], res_x[1]) > eps or
            multiply_vector(res_x[0], res_x[2]) > eps or
            multiply_vector(res_x[1], res_x[2]) > eps):
        return False
    return True


def rotation_iterations(matrix: [[float]], eps: float) -> [[float], [[float]], int]:
    iter_count = 0
    res_a_matrix = deepcopy(matrix)
    res_u_matrix = [[1 if i == j else 0 for i in range(len(matrix))] for j in range(len(matrix))]
    while t(res_a_matrix) > eps:
        iter_count += 1
        max_i, max_j = max_non_diagonal(res_a_matrix)
        u_matrix = rotation_matrix(phi(res_a_matrix, max_i, max_j), len(res_a_matrix), max_i, max_j)
        res_u_matrix = multiply_matrix(res_u_matrix, u_matrix)
        res_a_matrix = multiply_matrix(multiply_matrix(transport_matrix(u_matrix), res_a_matrix), u_matrix)
    res_x = [[res_u_matrix[j][i] for j in range(len(res_u_matrix))] for i in range(len(res_u_matrix))]
    if not check_orthogonal(res_x, eps):
        raise ValueError("Eigenvectors of matrix are not orthogonal")
    return [res_a_matrix[i][i] for i in range(len(res_a_matrix))], res_x, iter_count


def powers_method(matrix: [[float]], eps: float) -> (float, [float], int):
    j = 0
    last_pair = []
    y = [1] * len(matrix)
    while len(last_pair) != 2:
        tmp = multiply_matrix_vector(matrix, y)
        last_pair.append(tmp[j] / y[j])
        y = divide_vector(tmp, vector_norm(tmp))
    iter_count = 2
    while abs(last_pair[1] - last_pair[0]) >= eps:
        iter_count += 1
        tmp = multiply_matrix_vector(matrix, y)
        last_pair[0], last_pair[1] = last_pair[1], tmp[j] / y[j]
        y = divide_vector(tmp, vector_norm(tmp))
    return last_pair[1], y, iter_count


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Error: wrong number of arguments!')
    else:

        epsilon = 0.01 if len(sys.argv) < 4 else float(sys.argv[3])
        input_matrix = []
        try:
            with open(sys.argv[1], encoding="UTF-8") as file_in:
                for line in file_in.readlines():
                    if len(line) == 0:
                        print('Error: invalid input!')
                        exit(1)
                    input_matrix.append([int(x) for x in line.rstrip('\n').split()])

                if not check_matrix(input_matrix):
                    print('Error: invalid matrix size!')
                else:

                    eigenvalues, eigenvectors, rotation_iters = rotation_iterations(input_matrix, epsilon)
                    my_lambda, y_vector, powers_iters = powers_method(input_matrix, epsilon)

                    output = []
                    print('\n########## LOGS ##########')
                    output.append(f'{" Jacobi rotation method ".center(50, "/")}\n\n')
                    output.append('----- eigenvalues -----\n\n')
                    for num, lmb in enumerate(eigenvalues, 1):
                        if lmb == abs(lmb):
                            lmb = abs(lmb)
                        output.append(f'λ{num} = {round(lmb, 4)}\n')
                    else:
                        output.append('\n')

                    output.append('----- eigenvectors -----\n\n')
                    for num, h in enumerate(eigenvectors, 1):
                        nums = f'[ {"  ".join(str(round(i, 4)) for i in h)} ]T'
                        output.append(f'h{num} = {nums}\n')
                    else:
                        output.append('\n')

                    output.append('----- iterations -----\n\n')
                    output.append(f'{str(rotation_iters)}\n\n')

                    output.append('----- epsilon -----\n\n')
                    output.append(f'{str(epsilon)}\n\n')

                    output.append(f'\n{" Powers method ".center(50, "/")}\n\n')

                    output.append('----- eigenvalue -----\n\n')
                    output.append(f'λ = {round(my_lambda, 4)}\n\n')

                    output.append('----- eigenvector -----\n\n')
                    output.append(f'y = [ {"  ".join(str(round(i, 4)) for i in y_vector)} ]T\n\n')

                    output.append('----- iterations -----\n\n')
                    output.append(f'{str(powers_iters)}\n\n')

                    output.append('----- epsilon -----\n\n')
                    output.append(f'{str(epsilon)}')

                    output_str = ''.join(output)
                    print(output_str)
                    with open(sys.argv[2], 'w', encoding="UTF-8") as file_out:
                        file_out.write(output_str)

        except OSError:
            print('Error: could not open file!')
        except ValueError as ex:
            print('\n'.join(ex.args))
