import numpy as np


def jacobi(A, eps, max_iter):
    A = np.array(A, dtype=float)
    n = len(A)
    
    
    H = np.eye(n)

    for k in range(max_iter):
        max_value = 0
        i_max = 0
        j_max = 1

        for i in range(n):
            for j in range(i + 1, n):

                if abs(A[i, j]) > max_value:
                    max_value = abs(A[i, j])
                    i_max = i
                    j_max = j
        i = i_max
        j = j_max

        if abs(A[i, i] - A[j, j]) == 0:
            phi = np.pi / 4
        else:
            phi = 0.5 * np.arctan2(2 * A[i, j],A[i, i] - A[j, j])

        c = np.cos(phi)
        s = np.sin(phi)

        U = np.eye(n)

        U[i, i] = c
        U[j, j] = c

        U[i, j] = -s
        U[j, i] = s

        A = U.T @ A @ U

        H = H @ U

        t = 0
        
        for i1 in range(n):
            for j1 in range(n):
        
                if i1 != j1:
                    t += A[i1, j1] ** 2
        
        t = np.sqrt(t)
        
        if t <= eps:
            break

    else:
        print("Метод не сошёлся за заданное число итераций.")

    eigenvalues = np.diag(A)

    eigenvectors = H

    return eigenvalues, eigenvectors

A = [
    [5, -3, -4],
    [-3, -3, 4],
    [-4, 4, 0]
]

eigenvalues, eigenvectors = jacobi(A, 0.01, 100)


print("\nСобственные значения:")
for i in range(len(eigenvalues)):
    print("λ", i + 1, "=", eigenvalues[i])

print("\nСобственные векторы:")
for i in range(len(eigenvalues)):
    print("x", i + 1, "=", eigenvectors[:, i])