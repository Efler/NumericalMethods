import sys
import spline as spl


rounding = 4
points = []
with open(sys.argv[1], encoding="UTF-8") as file_in:
    while (line := file_in.readline()) != '\n':
        if len(line) == 0:
            print('Error: invalid input!')
            exit(1)
        points.append([float(x) for x in line.rstrip('\n').split()])
    x_aster = float(file_in.readline())

my_spline = spl.Spline(points)

output = []
a, b, c, d = my_spline.get_table()
head = ''.join(['i'.center(12, ' '), 'segment'.center(12, ' '),
                'a'.center(12, ' '), 'b'.center(12, ' '),
                'c'.center(12, ' '), 'd'.center(12, ' ')])
output.append(' Table '.center(len(head), '-'))
output.append(head)
output.append('-' * len(head))
for i in range(len(points) - 1):
    output.append(''.join([
        f'{i + 1}'.center(12, ' '),
        f'[{points[i][0]}, {points[i + 1][0]}]'.center(12, ' '),
        f'{round(a[i], rounding)}'.center(12, ' '),
        f'{round(b[i], rounding)}'.center(12, ' '),
        f'{round(c[i], rounding)}'.center(12, ' '),
        f'{round(d[i], rounding)}'.center(12, ' ')
    ]))

output.append(f'\n>>  Value of the Newton polynomial at the point x = {x_aster} '
              f'equals {round(my_spline.get_value(x_aster), rounding)}')

output_str = '\n'.join(output)
print(output_str)
with open(sys.argv[2], 'w', encoding="UTF-8") as file_out:
    file_out.write(output_str)

my_spline.get_plot(x_aster)
