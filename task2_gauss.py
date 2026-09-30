
def gauss_solve(A, f):
    n = len(A)
    M = [row[:] + [f[i]] for i, row in enumerate(A)]

    # Прямой ход
    for j in range(n - 1):
        for i in range(j + 1, n):
            mu = M[i][j] / M[j][j]
            for k in range(j, n + 1):
                M[i][k] -= M[j][k] * mu

    # Обратный ход
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = M[i][n]
        for j in range(i + 1, n):
            s -= M[i][j] * x[j]
        x[i] = s / M[i][i]

    return x


def gauss_det(A):
    n = len(A)
    M = [row[:] for row in A]

    det = 1.0
    for j in range(n - 1):  
        if M[j][j] == 0:
            return 0.0
        det *= M[j][j]
        for i in range(j + 1, n):
            mu = M[i][j] / M[j][j]
            for k in range(j, n):
                M[i][k] -= M[j][k] * mu

    det *= M[n - 1][n - 1]
    return det


def gauss_inverse(A):
    n = len(A)
    M = [row[:] + [1.0 if i == j else 0.0 for j in range(n)]
         for i, row in enumerate(A)]

    for j in range(n):
        piv = M[j][j]
        for k in range(j, 2 * n):
            M[j][k] /= piv

        for i in range(n):
            if i == j:
                continue
            mu = M[i][j]
            for k in range(j, 2 * n):
                M[i][k] -= M[j][k] * mu

    return [row[n:] for row in M]


A = [
    [8,  8, -5, -8],
    [8, -5,  9, -8],
    [5, -4, -6, -2],
    [8,  3,  6,  6],
]
f = [13, 38, 14, -95]

print("Решение:", gauss_solve([r[:] for r in A], f))
print("Определитель:", gauss_det(A))
print("Обратная матрица:")
for row in gauss_inverse(A):
    print([round(v, 6) for v in row])