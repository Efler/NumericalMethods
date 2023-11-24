import gauss
import numpy as np
import matplotlib.pyplot as plt


class LSM:

    __points = []
    __first_order_polynomial = []
    __second_order_polynomial = []

    def __normal_system_lsm(self):
        first_order_matrix = [[0.] * 2 for _ in range(2)]
        first_order_vector = [0.] * 2
        second_order_matrix = [[0.] * 3 for _ in range(3)]
        second_order_vector = [0.] * 3
        first_order_matrix[0][0] = len(self.__points)
        second_order_matrix[0][0] = len(self.__points)
        sum_x = 0.
        sum_x_pow_2 = 0.
        sum_x_pow_3 = 0.
        sum_x_pow_4 = 0.
        sum_y = 0.
        sum_y_x = 0.
        sum_y_x_pow_2 = 0.
        for i in range(len(self.__points)):
            sum_x += self.__points[i][0]
            sum_x_pow_2 += self.__points[i][0] ** 2
            sum_x_pow_3 += self.__points[i][0] ** 3
            sum_x_pow_4 += self.__points[i][0] ** 4
            sum_y += self.__points[i][1]
            sum_y_x += self.__points[i][0] * self.__points[i][1]
            sum_y_x_pow_2 += (self.__points[i][0] ** 2) * self.__points[i][1]

        (first_order_matrix[0][1], first_order_matrix[1][0], first_order_matrix[1][1],
         first_order_vector[0], first_order_vector[1]) = sum_x, sum_x, sum_x_pow_2, sum_y, sum_y_x

        (second_order_matrix[0][1], second_order_matrix[0][2], second_order_matrix[1][0],
         second_order_matrix[1][1], second_order_matrix[1][2], second_order_matrix[2][0],
         second_order_matrix[2][1], second_order_matrix[2][2]) = (
            sum_x, sum_x_pow_2, sum_x,
            sum_x_pow_2, sum_x_pow_3, sum_x_pow_2,
            sum_x_pow_3, sum_x_pow_4)
        second_order_vector[0], second_order_vector[1], second_order_vector[2] = sum_y, sum_y_x, sum_y_x_pow_2
        return first_order_matrix, first_order_vector, second_order_matrix, second_order_vector

    def __init__(self, points: list[list[float]]):
        self.__points = gauss.copy.deepcopy(points)
        first_order_m, first_order_v, second_order_m, second_order_v = self.__normal_system_lsm()
        self.__first_order_polynomial = gauss.find_solutions(first_order_m, first_order_v)
        self.__second_order_polynomial = gauss.find_solutions(second_order_m, second_order_v)

    def calculate_first_order_polynomial(self, x: float):
        return (self.__first_order_polynomial[0]
                + self.__first_order_polynomial[1] * x)

    def calculate_second_order_polynomial(self, x: float):
        return (self.__second_order_polynomial[0]
                + self.__second_order_polynomial[1] * x
                + self.__second_order_polynomial[2] * x ** 2)

    def quadratic_error_first_order(self):
        error = 0
        for i in range(len(self.__points)):
            error += (self.calculate_first_order_polynomial(self.__points[i][0]) - self.__points[i][1]) ** 2
        return error

    def quadratic_error_second_order(self):
        error = 0
        for i in range(len(self.__points)):
            error += (self.calculate_second_order_polynomial(self.__points[i][0]) - self.__points[i][1]) ** 2
        return error

    def get_first_order_polynomial(self, rounding: int):
        return (f'f (x)  =  {round(self.__first_order_polynomial[0], rounding)}'
                f'  +  {round(self.__first_order_polynomial[1], rounding)} * x')

    def get_second_order_polynomial(self, rounding: int):
        return (f'f (x)  =  {round(self.__second_order_polynomial[0], rounding)}'
                f'  +  {round(self.__second_order_polynomial[1], rounding)} * x'
                f'  +  {round(self.__second_order_polynomial[2], rounding)} * x^2')

    def make_plot(self):
        x_vals = np.arange(self.__points[0][0] - 2, self.__points[-1][0] + 2, 0.01)
        plt.plot(x_vals, [self.calculate_first_order_polynomial(i) for i in x_vals], label='first order', c='lime')
        plt.plot(x_vals, [self.calculate_second_order_polynomial(i) for i in x_vals], label='second order', c='blue')
        plt.scatter(*list(zip(*self.__points)), c='black')
        plt.xlabel(r'$x$', fontsize=14)
        plt.ylabel(r'$y$', fontsize=14)
        plt.grid(True)
        plt.legend(loc='best', fontsize=12)
        plt.show()
