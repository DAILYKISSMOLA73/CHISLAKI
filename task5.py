#начальное приближения брал из десмоса

import sympy as sp

def newton(f, var, x0, eps, max_iter):


    x = var
    df = sp.diff(f, x)

    f_func = sp.lambdify(x, f, 'math')
    df_func = sp.lambdify(x, df, 'math')

    x = x0

    for k in range(max_iter):
        x_new = x - f_func(x) / df_func(x)

        if abs(x_new - x) < eps:
            return x_new

        x = x_new

def simple_iter(phi, var, x0, eps, max_iter):
    
    x = var
    phi_func = sp.lambdify(x, phi, 'math')
    
    x = x0

    for k in range(max_iter):
        x_new = phi_func(x)

        if abs(x_new-x) < eps:
            return x_new
        x = x_new
    print("Метод не сошёлся за заданное число итераций.")








x = sp.symbols('x')
f = 4**x - 5*x - 2
phi = sp.log(5*x + 2) / sp.log(4)


print("Метод Ньютона: x=", (newton(f, x, 2, 0.001, 100)))
print("Метод простых итераций: x= ", simple_iter(phi, x, 2, 0.001, 100))

     