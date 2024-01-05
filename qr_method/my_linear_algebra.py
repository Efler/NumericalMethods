from copy import deepcopy
from math import sqrt, pow


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
