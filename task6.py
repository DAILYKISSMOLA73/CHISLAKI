import numpy as np
import sympy as sp


def simple_iteration(phi, var, x0, eps, max_iter):
    n = len(var)

    phi_func = sp.lambdify(var, phi, 'numpy')

    x = np.array(x0, dtype=float)

    for k in range(max_iter):
        x_new = np.array(phi_func(*x), dtype=float)

        if np.max(np.abs(x_new - x)) < eps:
            return x_new

        x = x_new

    print("Метод итераций не сошёлся за заданное число итераций")
    return x


x, y = sp.symbols('x y')

phi = [
    sp.cos(y)/3,
    (sp.E**x)/3
]

x0 = [0.3, 0.3]

result = simple_iteration(phi, [x, y], x0, 0.001, 1000)

print("Метод простых итераций:")
for i in range(len(result)):
    print("x" + str(i + 1), "=", result[i])