import numpy as np
import sympy as sp

def peregonka(a, b, c, f):
    n = len(f)
    alpha = [0.0]*n
    beta = [0.0]*n
    alpha[0] = -c[0]/b[0]
    beta[0] = f[0]/b[0]

    for i in range(1, n):
        den = b[i] + a[i]*alpha[i-1]
        alpha[i] = -c[i]/den
        beta[i] = (f[i]-a[i]*beta[i-1])/den

    x = [0.0]*n
    x[-1] = beta[-1]

    for i in range(n-2, -1, -1):
        x[i] = alpha[i]*x[i+1] + beta[i]
    return x


def seidel(A, b, eps, max_iter):

    n = len(b)
    x = np.zeros(n)

    alpha = np.zeros((n, n))
    beta = np.zeros(n)

    for i in range(n):
        beta[i] = b[i]/A[i][i]

        for j in range(n):
            if i != j:
                alpha[i][j] = -A[i][j]/A[i][i]

    x = beta.copy()

    for k in range(max_iter):

        x_old = x.copy()

        for i in range(n):

            x[i] = beta[i]

            for j in range(n):

                if j < i:
                    x[i] += alpha[i][j]*x[j]

                elif j > i:
                    x[i] += alpha[i][j]*x_old[j]

        error = np.max(np.abs(x - x_old))

        if error < eps:
            return x

    print("Метод Зейделя не сошелся за заданное количество итераций")
    return None



def spline(x, y):
    n = len(x)
    h = [x[i+1]-x[i] for i in range(n-1)]
    m = n-1
    C = [0.0]*n

    if m > 1:
        a = [0.0]*(m-1)
        b = [0.0]*(m-1)
        c = [0.0]*(m-1)
        f = [0.0]*(m-1)

        for i in range(1, m):
            k = i-1
            b[k] = 2*(h[i-1]+h[i])
            if k > 0:
                a[k] = h[i-1]
            if k < m-2:
                c[k] = h[i]
            f[k] = 3*((y[i+1]-y[i])/h[i] - (y[i]-y[i-1])/h[i-1])

        C_inner = peregonka(a, b, c, f)

        for i in range(1, m):
            C[i] = C_inner[i-1]

    A = []
    B = []
    D = []
    for i in range(m):
        A.append(y[i])
        B.append((y[i+1]-y[i])/h[i] - h[i]*(2*C[i]+C[i+1])/3)
        D.append((C[i+1]-C[i])/(3*h[i]))
        
    return A, B, C[:-1], D


def spline_value(x0, x, A, B, C, D):
    n = len(x)
    i = np.searchsorted(x, x0) - 1
    if i < 0: i = 0
    if i >= n-1: i = n-2
    dx = x0-x[i]
    return (A[i] + B[i]*dx + C[i]*dx**2 + D[i]*dx**3)

def mnk(x, y, m):

    n = len(x)

    S = [0.0]*(2*m+1)
    t = [0.0]*(m+1)

    for k in range(2*m+1):

        for i in range(n):
            S[k] += x[i]**k

    for k in range(m+1):

        for i in range(n):
            t[k] += (x[i]**k)*y[i]

    A = [[0.0]*(m+1) for i in range(m+1)]

    for i in range(m+1):

        for j in range(m+1):
            A[i][j] = S[i+j]

    b = t

    a = seidel(A, b, 0.000001, 10000)

    phi = 0

    for i in range(n):

        y_calc = 0

        for j in range(m+1):
            y_calc += a[j]*x[i]**j

        phi += (y_calc-y[i])**2

    return a, phi

# сплайн
x = [-2.0, -1.0, 0.0, 1.0, 2.0]
y = [-1.8647, -0.63212, 1.0, 3.7183, 9.3891]

A, B, C, D = spline(x, y)

X = sp.symbols('x')

print("\nСплайн\n")
for i in range(len(x)-1):

    dx = X - x[i]

    S = (A[i] + B[i]*dx
         + C[i]*dx**2 + D[i]*dx**3)

    print(f"{x[i]} <= x <= {x[i+1]}:")
    print("S(x) =", sp.expand(S))
    print()

x0 = -0.5 
result = spline_value(x0, x, A, B, C, D)
print(f"S({x0}) = {result:.6f}")

# МНК

x2 = [-3, -2, -1, 0, 1, 2]
y2 = [-2.9502, -1.8647, -0.63212, 1.0, 3.7183, 9.3891]

n = len(x2)

a1, phi1 = mnk(x2, y2, 1)

a2, phi2 = mnk(x2, y2, 2)

print("\nМНК")
print("\nМногочлен первой степени:")
print(f"Q1(x) = {a1[0]:.6f} + ({a1[1]:.6f})*x")
print(f"Сумма квадратов ошибок Phi = {phi1:.6f}")

print("\nМногочлен второй степени:")
print(f"Q2(x) = {a2[0]:.6f} + ({a2[1]:.6f})*x + ({a2[2]:.6f})*x**2")
print(f"Сумма квадратов ошибок Phi = {phi2:.6f}")