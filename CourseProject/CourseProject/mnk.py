import my_linear_algebra as lin_al


def mnk(x_arr: list[float], y_arr: list[float], func_array) -> list[float]:

    def composition_y(x_array, y_array, function):
        res = 0
        for k in range(len(x_array)):
            res += function(x_array[k]) * y_array[k]
        return res

    def composition_x(x_array, function1, function2):
        res = 0
        for k in range(len(x_array)):
            res += function1(x_array[k]) * function2(x_array[k])
        return res

    matrix = []
    vector = []
    for i in range(len(func_array)):
        prom = []
        for j in range(len(func_array)):
            prom.append(composition_x(x_arr, func_array[i], func_array[j]))
        matrix.append(prom)
        vector.append(composition_y(x_arr, y_arr, func_array[i]))
    coefficients = lin_al.find_solutions(matrix, vector)

    return coefficients
