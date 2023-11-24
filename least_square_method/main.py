import sys
import least_square_method


if __name__ == '__main__':
    rounding = 4
    points = []
    with open(sys.argv[1], encoding="UTF-8") as file_in:
        for line in file_in.readlines():
            if len(line) == 0:
                print('Error: invalid input!')
                exit(1)
            points.append([float(x) for x in line.rstrip('\n').split()])

    lsm = least_square_method.LSM(points)

    output = ['First order polynomial LSM:',
              f'{lsm.get_first_order_polynomial(4)}',
              'Sum of squared errors:',
              f'{round(lsm.quadratic_error_first_order(), rounding)}\n',
              'Second order polynomial LSM:',
              f'{lsm.get_second_order_polynomial(4)}',
              'Sum of squared errors:',
              f'{round(lsm.quadratic_error_second_order(), rounding)}']

    output_str = '\n'.join(output)
    print(output_str)
    with open(sys.argv[2], 'w', encoding="UTF-8") as file_out:
        file_out.write(output_str)

    lsm.make_plot()
