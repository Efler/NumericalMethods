import sys


def numerical_differentiation(points, x):
    index = 0
    for i in range(1, len(points) - 2):
        if points[i][0] <= x <= points[i + 1][0]:
            index = i
            break
    left_handed_derivative = ((points[index + 1][1] - points[index][1])
                              / (points[index + 1][0] - points[index][0]))
    right_handed_derivative = ((points[index + 2][1] - points[index + 1][1])
                               / (points[index + 2][0] - points[index + 1][0]))
    first_der = (left_handed_derivative
                 + ((right_handed_derivative - left_handed_derivative)
                    / (points[index + 2][0] - points[index][0]))
                 * (2. * x - points[index][0] - points[index + 1][0]))
    second_der = (2. * (right_handed_derivative - left_handed_derivative)
                  / (points[index + 2][0] - points[index][0]))
    return left_handed_derivative, right_handed_derivative, first_der, second_der


if __name__ == '__main__':
    rounding = 4
    input_points = []
    with open(sys.argv[1], encoding="UTF-8") as file_in:
        while (line := file_in.readline()) != '\n':
            if len(line) == 0:
                print('Error: invalid input!')
                exit(1)
            input_points.append([float(x) for x in line.rstrip('\n').split()])
        x_aster = float(file_in.readline())

    left_handed, right_handed, first_derivative, second_derivative = numerical_differentiation(input_points, x_aster)

    output = [f'left handed derivative -> {round(left_handed, rounding)}',
              f'right handed derivative -> {round(right_handed, rounding)}',
              f'first order derivative -> {round(first_derivative, rounding)}',
              f'second order derivative -> {round(second_derivative, rounding)}']

    output_str = '\n'.join(output)
    print(output_str)
    with open(sys.argv[2], 'w', encoding="UTF-8") as file_out:
        file_out.write(output_str)
