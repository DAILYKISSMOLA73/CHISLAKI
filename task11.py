
import numpy as np

def finite_difference(N):

    a = 0
    b = 1

    h = (b - a) / N
    x = np.linspace(a, b, N + 1)

    A = np.zeros((N + 1, N + 1))
    B = np.zeros(N + 1)

    A[0, 0] = -3
    A[0, 1] = 4
    A[0, 2] = -1

    for k in range(1, N):

        p = (x[k] - 3) / (x[k]**2 - 1)
        q = -1 / (x[k]**2 - 1)
        f = 0

        A[k, k-1] = 1 - p*h/2
        A[k, k] = -2 + q*h**2
        A[k, k+1] = 1 + p*h/2

        B[k] = h**2 * f


    A[N, N-2] = 1
    A[N, N-1] = -4
    A[N, N] = 3 + 2*h

    B[N] = -1.5*h

    y = np.linalg.solve(A, B)

    return x, y

N = 10

x, y = finite_difference(N)

y_exact = x - 3 + 1/(x + 1)

error = np.abs(y - y_exact)

print("      x       y численное       y точное       Погрешность")

for k in range(N + 1):
    print(f"{x[k]:.2f}    {y[k]:.8f}    "
          f"{y_exact[k]:.8f}    {error[k]:.8f}")


