import copy
import sys


def beautify_matrix(in_matrix):
    str_matrix = []
    res_matrix = []

    for row in in_matrix:
        if isinstance(row, int | float):
            str_matrix.append(str(row))
        else:
            for el in row:
                str_matrix.append(str(el))
    max_width = max(map(lambda x: len(x), str_matrix))
    for row in in_matrix:
        if isinstance(row, int | float):
            res_matrix.append(str(row).rjust(max_width, ' '))
        else:
            res_matrix.append([str(x).rjust(max_width, ' ') for x in row])

    return res_matrix


def check_matrix(in_matrix):
    n = len(in_matrix)
    for row in in_matrix:
        if len(row) != n:
            return False
    return True


def round_up(*args, n=2):
    res_args = []

    for arg in args:
        res_arg = copy.deepcopy(arg)
        for i in range(len(res_arg)):
            if isinstance(res_arg[i], int | float):
                res_arg[i] = round(res_arg[i], n)
            else:
                for j in range(len(res_arg[i])):
                    res_arg[i][j] = round(res_arg[i][j], n)
        res_args.append(res_arg)

    return res_args


def make_unit_matrix(n):
    matrix = []

    for i in range(n):
        row = []
        for j in range(n):
            if i == j:
                row.append(1)
            else:
                row.append(0)
        matrix.append(row)

    return matrix


def inverse_matrix(in_matrix):
    height = len(in_matrix)
    unit_matrix = make_unit_matrix(len(in_matrix))
    res_matrix, tmp, unit_matrix = to_triangular(in_matrix, unit_matrix)

    for column in range(height - 1, -1, -1):
        lead_el = res_matrix[column][column]
        for row in range(column - 1, -1, -1):
            multiplier = -(res_matrix[row][column] / lead_el)
            for k in range(column, -1, -1):
                res_matrix[row][k] += res_matrix[column][k] * multiplier
            if isinstance(unit_matrix[column], int | float):
                unit_matrix[row] += unit_matrix[column] * multiplier
            else:
                for k in range(len(unit_matrix[column]) - 1, -1, -1):
                    unit_matrix[row][k] += unit_matrix[column][k] * multiplier

    for column in range(height):
        multiplier = 1 / res_matrix[column][column]
        res_matrix[column][column] *= multiplier
        if isinstance(unit_matrix[column], int | float):
            unit_matrix[column] *= multiplier
        else:
            for k in range(len(unit_matrix[column])):
                unit_matrix[column][k] *= multiplier

    return unit_matrix


def determinant_from_triangular(triangular_matrix, permutations_count):
    product = 1
    height = len(triangular_matrix)

    for i in range(height):
        product *= triangular_matrix[i][i]

    return product * ((-1) ** permutations_count)


def to_triangular(in_matrix, *args):
    height = len(in_matrix)
    res_matrix = copy.deepcopy(in_matrix)
    res_args = copy.deepcopy(args)
    perm_count = 0

    for column in range(height):
        res_matrix, perm_count, *res_args = leading_element(res_matrix, column, perm_count, *res_args)
        lead_el = res_matrix[column][column]
        for row in range(column + 1, height):
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


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Error: wrong number of arguments!')
    else:
        input_matrix = []
        b_vector = None
        try:
            with open(sys.argv[1], encoding="UTF-8") as file_in:
                while (line := file_in.readline()) != '\n':
                    input_matrix.append([int(x) for x in line.rstrip('\n').split()])
                b_vector = [int(x.rstrip('\n')) for x in file_in.readlines()]

                if not check_matrix(input_matrix) or len(input_matrix) != len(b_vector):
                    print('Error: invalid matrix or b vector size!')
                else:
                    triangular, permutations = to_triangular(input_matrix)
                    determinant = determinant_from_triangular(triangular, permutations)
                    if determinant == 0:
                        print('Error: singular matrix (determinant equals zero)!')
                    else:
                        solutions = find_solutions(input_matrix, b_vector)
                        inverse = inverse_matrix(input_matrix)
                        solutions, triangular, inverse = round_up(solutions, triangular, inverse, n=3)

                        output = []
                        print('\n########## LOGS ##########')
                        output.append('----- input matrix -----\n\n')
                        str_input_matrix = beautify_matrix(input_matrix)
                        for r in str_input_matrix:
                            output.append(f"{' '.join(r)}\n")
                        else:
                            output.append('\n')

                        output.append('----- b vector -----\n\n')
                        str_b_vector = beautify_matrix(b_vector)
                        for r in str_b_vector:
                            output.append(f"{r}\n")
                        else:
                            output.append('\n')

                        output.append('----- triangular matrix -----\n\n')
                        str_triangular_matrix = beautify_matrix(triangular)
                        for r in str_triangular_matrix:
                            output.append(f"{' '.join(r)}\n")
                        else:
                            output.append('\n')

                        output.append('----- permutations -----\n\n')
                        output.append(f'{str(permutations)}\n\n')

                        output.append('----- solutions -----\n\n')
                        for num, x in enumerate(solutions, 1):
                            output.append(f'x{num} = {x}\n')
                        else:
                            output.append('\n')

                        output.append('----- determinant -----\n\n')
                        output.append(f'{str(determinant)}\n\n')

                        output.append('----- inverse matrix -----\n\n')
                        str_inverse_matrix = beautify_matrix(inverse)
                        for r in str_inverse_matrix:
                            output.append(f"{' '.join(r)}\n")

                        output_str = ''.join(output)
                        print(output_str)
                        with open(sys.argv[2], 'w', encoding="UTF-8") as file_out:
                            file_out.write(output_str)
        except OSError:
            print('Error: could not open file!')
