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

x = [-2.0, -1.0, 0.0, 1.0, 2.0]
y = [-1.8647, -0.63212, 1.0, 3.7183, 9.3891]

A, B, C, D = spline(x, y)

X = sp.symbols('x')

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
