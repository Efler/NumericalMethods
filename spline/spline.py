import copy

import gauss
import numpy as np
import matplotlib.pyplot as plt


class Spline:

    __a = []
    __b = []
    __c = []
    __d = []
    __h = []
    __points = []

    def __get_index(self, x: float) -> int:
        index = 0
        for i in range(1, len(self.__points)):
            if self.__points[i - 1][0] <= x <= self.__points[i][0]:
                index = i
                break
        if index == 0:
            raise ValueError
        return index

    def __make_linear_c(self):
        matrix = [[0.] * (len(self.__points) - 2) for _ in range(len(self.__points) - 2)]
        vector = [0.] * (len(self.__points) - 2)
        matrix[0][0] = (2 * (self.__h[0] + self.__h[1]))
        matrix[0][1] = self.__h[1]
        vector[0] = 3 * (
                (self.__points[2][1] - self.__points[1][1]) / self.__h[1]
                - (self.__points[1][1] - self.__points[0][1]) / self.__h[0]
        )
        for j in range(1, len(matrix) - 1):
            matrix[j][j - 1] = self.__h[j]
            matrix[j][j] = (self.__h[j] + self.__h[j + 1]) * 2
            matrix[j][j + 1] = self.__h[j + 1]
            vector[j] = 3 * (
                    (self.__points[j + 2][1] - self.__points[j + 1][1]) / self.__h[j + 1]
                    - ((self.__points[j + 1][1] - self.__points[j][1]) / self.__h[j])
            )
        matrix[len(matrix) - 1][len(matrix) - 2] = self.__h[len(matrix) - 1]
        matrix[len(matrix) - 1][len(matrix) - 1] = 2 * (
                self.__h[len(matrix) - 1] + self.__h[len(matrix)]
        )
        vector[len(matrix) - 1] = 3 * (
                (self.__points[len(matrix) + 1][1] - self.__points[len(matrix)][1])
                / self.__h[len(matrix)]
                - (self.__points[len(matrix)][1] - self.__points[(len(matrix) - 1)][1])
                / self.__h[(len(matrix) - 1)]
        )
        return matrix, vector

    def __init__(self, points: list[list[float]]):
        self.__points = copy.deepcopy(points)
        self.__h = [self.__points[i + 1][0] - self.__points[i][0] for i in range(len(self.__points) - 1)]

        c_matrix, c_vector = self.__make_linear_c()
        self.__c = gauss.find_solutions(c_matrix, c_vector)
        self.__c = [0] + self.__c

        for i in range(1, len(points)):
            self.__a.append(points[i - 1][1])

        for i in range(1, len(points) - 1):
            self.__b.append(
                (points[i][1] - points[i - 1][1]) / self.__h[i - 1]
                - (self.__h[i - 1] * (self.__c[i] + 2 * self.__c[i - 1]) / 3)
            )
            self.__d.append((self.__c[i] - self.__c[i - 1]) / (3 * self.__h[i - 1]))
        self.__b.append(
            (points[len(points) - 1][1] - points[len(points) - 2][1]) / self.__h[len(points) - 2]
            - (2 * self.__h[len(points) - 2] * self.__c[len(points) - 2]) / 3
        )
        self.__d.append(-1 * self.__c[len(points) - 2] / (3 * self.__h[len(points) - 2]))

    def get_value(self, x: float):
        index = self.__get_index(x)
        return (self.__a[index - 1]
                + self.__b[index - 1] * (x - self.__points[index - 1][0])
                + self.__c[index - 1] * ((x - self.__points[index - 1][0]) ** 2)
                + self.__d[index - 1] * ((x - self.__points[index - 1][0]) ** 3))

    def get_table(self):
        return self.__a, self.__b, self.__c, self.__d

    def get_plot(self, x: float):
        x_vals = np.arange(self.__points[0][0], self.__points[-1][0], 0.01)
        plt.plot(x_vals, [self.get_value(i) for i in x_vals], label=r'$spline$')
        plt.scatter(*list(zip(*self.__points)))
        plt.scatter(x, self.get_value(x))
        plt.xlabel(r'$x$', fontsize=14)
        plt.ylabel(r'$y$', fontsize=14)
        plt.grid(True)
        plt.legend(loc='best', fontsize=12)
        plt.show()
