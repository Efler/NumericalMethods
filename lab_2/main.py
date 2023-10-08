import sys


def check_matrix(in_matrix):
    n = len(in_matrix)
    for row in in_matrix:
        if len(row) != n:
            return False
    return True


def get_p_and_q(in_matrix, d_vector, acc=3):
    p = []
    q = []
    errors = []
    a_index = -1
    b_index = 0
    c_index = 1

    for (i, row), d in zip(enumerate(in_matrix, 1), d_vector):
        a = 0 if a_index < 0 else row[a_index]
        b = row[b_index]
        c = 0 if c_index >= len(row) else row[c_index]

        # ???
        if a == 0 and i != 1:
            errors.append(f'row #{i}: condition [a != 0] ---> False\n')
        if c == 0 and i != len(row):
            errors.append(f'row #{i}: condition [c != 0] ---> False\n')
        if abs(b) < (abs(a) + abs(c)):
            errors.append(f'row #{i}: condition [|b| >= |a| + |c|] ---> False\n')

        if i == 1:
            q.append(round(d / b, acc))
            p.append(round(-c / b, acc))
        else:
            denominator = b + a * p[-1]
            q.append(round((d - a * q[-1]) / denominator, acc))
            p.append(0 if i == len(row) else round(-c / denominator, acc))
        a_index += 1
        b_index += 1
        c_index += 1

    return p, q, errors


def find_solutions(p, q, acc=3):
    solutions = [None for _ in range(len(p))]

    for i in range(len(p) - 1, -1, -1):
        if i == len(p) - 1:
            solutions[i] = round(q[i], acc)
        else:
            solutions[i] = round(p[i] * solutions[i + 1] + q[i], acc)
    return solutions


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Error: wrong number of arguments!')
    else:
        accuracy = 3 if len(sys.argv) < 4 else int(sys.argv[3])
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
                    p_values, q_values, cautions = get_p_and_q(input_matrix, b_vector)
                    x_values = find_solutions(p_values, q_values, acc=accuracy)

                    output = []
                    print('\n########## LOGS ##########')
                    output.append('----- solutions -----\n\n')
                    for num, x in enumerate(x_values, 1):
                        if x == abs(x):
                            x = abs(x)
                        output.append(f'x{num} = {x}\n')
                    else:
                        output.append('\n')

                    output.append('----- P and Q -----\n\n')
                    for k, (p_val, q_val) in enumerate(zip(p_values, q_values), 1):
                        output.append(f'P{k} = {p_val}  Q{k} = {q_val}\n')

                    if len(cautions) != 0:
                        output.append('\n----- cautions -----\n\n')
                        for k in cautions:
                            output.append(k)

                    output_str = ''.join(output)
                    print(output_str)
                    with open(sys.argv[2], 'w', encoding="UTF-8") as file_out:
                        file_out.write(output_str)

        except OSError:
            print('Error: could not open file!')
