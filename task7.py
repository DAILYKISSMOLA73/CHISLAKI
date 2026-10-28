import sympy as sp

def lagrange(nodes, f_expr, x):
    n = len(nodes)
    P = 0

    for i in range(n):
        L = 1

        for j in range(n):
            if i != j:
                L = L * (x - nodes[j]) / (nodes[i] - nodes[j])

        P = P + f_expr.subs(x, nodes[i]) * L

    return sp.expand(P)


def newton(nodes, f_expr, x):
    n = len(nodes)

    y = [f_expr.subs(x, xi) for xi in nodes]

    differences = [y.copy()]

    for k in range(1, n):
        prev = differences[k - 1]
        current = []

        for i in range(n - k):
            value = (prev[i + 1] - prev[i]) / (
                nodes[i + k] - nodes[i]
            )
            current.append(value)

        differences.append(current)

    P = differences[0][0]
    product = 1

    for k in range(1, n):
        product = product * (x - nodes[k - 1])
        P = P + differences[k][0] * product

    return sp.expand(P)

x = sp.symbols('x')
f_expr = sp.exp(x)+x
f = sp.lambdify(x, f_expr, 'numpy')

nodes = [-2, -1, 0, 1]
x_check = -0.5

L_lagrange = lagrange(nodes, f_expr, x)

P_newton = newton(nodes, f_expr, x)

print("Исходная функция:")
print("f(x) =", f_expr)

print("\nМногочлен Лагранжа:")
sp.pprint(L_lagrange)

print("\nМногочлен Ньютона:")
sp.pprint(P_newton)

L_func = sp.lambdify(x, L_lagrange, 'numpy')
P_func = sp.lambdify(x, P_newton, 'numpy')

true_value = float(f_expr.subs(x, x_check))
lagrange_value = float(L_func(x_check))
newton_value = float(P_func(x_check))

error_lagrange = abs(true_value - lagrange_value)
error_newton = abs(true_value - newton_value)

print("Абсолютная погрешность Лагранжа:", error_lagrange)

print("Абсолютная погрешность Ньютона:", error_newton)
