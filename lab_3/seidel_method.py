import copy
import sys


def check_matrix(in_matrix):
    n = len(in_matrix)
    for row in in_matrix:
        if len(row) != n:
            return False
    return True


def convergence_condition(matrix):
    for i in range(len(matrix)):
        if abs(matrix[i][i]) <= sum(abs(x) for x in matrix[i] if x != matrix[i][i]):
            return False
    return True


def make_iteration(matrix, b, x):
    result_x = copy.deepcopy(x)
    for i in range(len(result_x)):
        tmp = b[i]
        for j in range(len(matrix[i])):
            if i != j:
                tmp -= matrix[i][j] * result_x[j]
        result_x[i] = (tmp / matrix[i][i])
    return result_x


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Error: wrong number of arguments!')
    else:

        epsilon = 0.01 if len(sys.argv) < 4 else float(sys.argv[3])
        input_matrix = []
        b_vector = None
        try:
            with open(sys.argv[1], encoding="UTF-8") as file_in:
                while (line := file_in.readline()) != '\n':
                    if len(line) == 0:
                        print('Error: invalid input!')
                        exit(1)
                    input_matrix.append([int(x) for x in line.rstrip('\n').split()])
                b_vector = [int(x.rstrip('\n')) for x in file_in.readlines()]

                if not check_matrix(input_matrix) or len(input_matrix) != len(b_vector):
                    print('Error: invalid matrix or b vector size!')
                else:

                    if not convergence_condition(input_matrix):
                        print('Error: the convergence condition is not satisfied!')
                        exit(1)

                    counter = 0
                    accuracy_condition = False
                    prev_x = [0 for _ in range(len(input_matrix[0]))]
                    while not accuracy_condition:
                        counter += 1
                        new_x = make_iteration(input_matrix, b_vector, prev_x)
                        max_accuracy = max(abs(x - y) for x, y in zip(new_x, prev_x))
                        if max_accuracy < epsilon:
                            accuracy_condition = True
                        prev_x = new_x

                    output = []
                    print('\n########## LOGS ##########')
                    output.append('----- solutions -----\n\n')
                    for num, x in enumerate(new_x, 1):
                        if x == abs(x):
                            x = abs(x)
                        output.append(f'x{num} = {x}\n')
                    else:
                        output.append('\n')

                    output.append('----- iterations -----\n\n')
                    output.append(f'{counter}\n\n')

                    output.append('----- epsilon -----\n\n')
                    output.append(str(epsilon))

                    output_str = ''.join(output)
                    print(output_str)
                    with open(sys.argv[2], 'w', encoding="UTF-8") as file_out:
                        file_out.write(output_str)

        except OSError:
            print('Error: could not open file!')
