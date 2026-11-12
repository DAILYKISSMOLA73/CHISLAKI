import numpy as np

a, b = -2, 2
h1, h2 = 1.0, 0.5

def f(x):
    return 1 / (256 - x**4)


def integrate(method, h):
    x = np.arange(a, b + h / 2, h)
    y = f(x)
    n = len(x) - 1

    if method == "rectangles":            
        xm = (x[:-1] + x[1:]) / 2
        return h * np.sum(f(xm))

    if method == "trapezoids":              
        return h * (y[0] / 2 + np.sum(y[1:-1]) + y[-1] / 2)

    if method == "simpson":                
        return h / 3 * (y[0] + y[-1]
                       + 4 * np.sum(y[1:-1:2])
                       + 2 * np.sum(y[2:-1:2]))


methods = [
    ("Средние прямоугольники", "rectangles", 2),
    ("Трапеции", "trapezoids", 2),
    ("Симпсона", "simpson", 4)
]

for name, method, p in methods:
    I1 = integrate(method, h1)
    I2 = integrate(method, h2)

    I = I2 + (I2 - I1) / (2**p - 1)
    R = abs(I2 - I1) / (2**p - 1)

    print(name)
    print(f"  При h1 = {h1}: {I1:.10f}")
    print(f"  При h2 = {h2}: {I2:.10f}")
    print(f"  Уточнённое значение: {I:.10f}")

xi = np.array([-0.2, 0, 0.2, 0.4, 0.6])
yi = np.array([-0.40136, 0.0, 0.40136, 0.81152, 1.2435])

x_star = 0.2

i = np.where(np.isclose(xi, x_star))[0][0]

h = xi[1] - xi[0]

first_derivative = (yi[i+1] - yi[i-1]) / (2*h)

second_derivative = (
    yi[i+1] - 2*yi[i] + yi[i-1]
) / h**2

print("Дифференцирование")
print(f"x* = {x_star}, h = {h}")

print(f"Первая производная y'({x_star}) = "
      f"{first_derivative:.5f}")

print(f"Вторая производная y''({x_star}) = "
      f"{second_derivative:.5f}")
