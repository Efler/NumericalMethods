import copy


def to_triangular(in_matrix, *args):
    height = len(in_matrix)
    res_matrix = copy.deepcopy(in_matrix)
    res_args = copy.deepcopy(args)
    perm_count = 0

    for column in range(height):
        res_matrix, perm_count, *res_args = leading_element(res_matrix, column, perm_count, *res_args)
        lead_el = res_matrix[column][column]
        for row in range(column + 1, height):
            if not lead_el == 0:
                multiplier = -(res_matrix[row][column] / lead_el)
                for k in range(column, height):
                    res_matrix[row][k] += res_matrix[column][k] * multiplier
                for arg in res_args:
                    if isinstance(arg[column], int | float):
                        arg[row] += arg[column] * multiplier
                    else:
                        for k in range(len(arg[column])):
                            arg[row][k] += arg[column][k] * multiplier

    if len(args) > 0:
        return res_matrix, perm_count, *res_args
    else:
        return res_matrix, perm_count


def leading_element(matrix, index, permutations_count, *args):
    res_matrix = copy.deepcopy(matrix)
    res_args = copy.deepcopy(args)
    max_index = index
    max_el = res_matrix[index][index]

    for i in range(index + 1, len(res_matrix)):
        if res_matrix[i][index] > max_el:
            max_index = i
            max_el = res_matrix[i][index]
    if max_index != index:
        res_matrix[index], res_matrix[max_index] = res_matrix[max_index], res_matrix[index]
        for arg in res_args:
            arg[index], arg[max_index] = arg[max_index], arg[index]
        permutations_count += 1

    if len(res_args) > 0:
        return res_matrix, permutations_count, *res_args
    else:
        return res_matrix, permutations_count


def find_solutions(in_matrix, b):
    height = len(in_matrix)
    res = [None for _ in range(height)]
    triangular_matrix, tmp, b = to_triangular(in_matrix, b)

    for row in range(height - 1, -1, -1):
        result = b[row]
        for column in range(row + 1, height):
            result -= triangular_matrix[row][column] * res[column]
        res[row] = result / triangular_matrix[row][row]

    return res
