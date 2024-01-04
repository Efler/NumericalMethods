from copy import deepcopy
import my_linear_algebra as la


def get_v_vector(a_vector: [float], k=0) -> list[float]:
    v_vector = [0] * k
    for n in range(k, len(a_vector)):
        if n == k:
            v_vector.append(a_vector[n] + la.sign(a_vector[n]) * la.vector_euclidean_norm(a_vector, k))
        else:
            v_vector.append(a_vector[n])
    return v_vector


def householder_matrix(v_vector: list[float], eps: int = None) -> list[list[float]]:
    cf = 2 / la.multiply_vectors_to_value(v_vector, la.transpose(v_vector))
    vv_matrix = la.multiply_vectors_to_matrix(la.transpose(v_vector), v_vector)
    eye_matrix = la.eye(len(vv_matrix))
    for j in range(len(vv_matrix)):
        for k in range(len(vv_matrix)):
            eye_matrix[j][k] -= vv_matrix[j][k] * cf
            if eps is not None:
                eye_matrix[j][k] = round(eye_matrix[j][k], eps)
    return eye_matrix


def decompose(matrix: [[float]], eps: int = None) -> tuple[[[float]], [[float]]]:
    if len(matrix) != len(matrix[0]):
        raise TypeError('Invalid matrix')
    a_matrix = deepcopy(matrix)
    h_matrices = []
    for k in range(len(matrix) - 1):
        h_k = householder_matrix(get_v_vector([row[k] for row in a_matrix], k))
        h_matrices.append(h_k)
        a_matrix = la.multiply_matrix(h_k, a_matrix, eps if k == len(matrix) - 2 else None)
    q_matrix = h_matrices[0]
    for k in range(1, len(h_matrices)):
        q_matrix = la.multiply_matrix(q_matrix, h_matrices[k],
                                      eps if k == len(h_matrices) - 1 else None)
    return q_matrix, a_matrix
