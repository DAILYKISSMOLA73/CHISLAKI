import numpy as np

def simple_iteration(A, b, eps, max_iter):
    n = len(b)
    x = np.zeros(n)

    alpha = np.zeros((n, n))
    beta = np.zeros(n)

    for i in range(n):
        beta[i] = b[i]/A[i][i]

        for j in range(n):
            if i != j:
                alpha[i][j] = - A[i][j]/A[i][i]

    q = np.max(np.sum(np.abs(alpha), axis = 1))
    if q >= 1:
        print("Условие сходимости не выполняется")
        return None

    x = beta.copy()

    for k in range(max_iter):
        x_new = alpha @ x + beta
        error = q / (1 - q) * np.max(np.abs(x_new - x))
        if error < eps:
            return x_new
        x = x_new
    print("Метод простых итераций не сошелся за заданное количество итераций")
    return None

def seidel(A, b, eps, max_iter):
    n = len(b)
    x = np.zeros(n)
    
    alpha = np.zeros((n, n))
    beta = np.zeros(n)
    
    for i in range(n):
        beta[i] = b[i]/A[i][i]
    
        for j in range(n):
            if i != j:
                alpha[i][j] = - A[i][j]/A[i][i]
    
    q = np.max(np.sum(np.abs(alpha), axis = 1))
    if q >= 1:
        print("Условие сходимости не выполняется")
        return None

    x = beta.copy()
    for k in range(max_iter):
        x_old = x.copy()
        for i in range(n):
            x[i] = beta[i]
            for j in range(n):
                if j < i:
                    x[i] += alpha[i, j] * x[j]
                elif j > i:
                    x[i] += alpha[i, j] * x_old[j]

        
        error = np.max(np.abs(x - x_old))
        if error < eps:
            return x
    print("Метод Зейделя не сошелся за заданное количество итераций")
    return None


A = np.array([
    [-19,  2, -1, -8],
    [  2, 14,  0, -4],
    [  6, -5, -20, -6],
    [ -6,  4, -2, 15]
], dtype=float)

b = np.array([38, 20, 52, 43], dtype=float)

print("Решиние методом итераций: ", simple_iteration(A, b, 0.01, 100))
print("Решение методом Зейделя: ", seidel(A, b, 0.01, 100))
