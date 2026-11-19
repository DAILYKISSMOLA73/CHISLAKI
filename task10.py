
import math

def f(x, y):
    return (-y + x * math.sin(x)) / x

def y_exact(x):
    return (math.sin(x) - x * math.cos(x)) / x


a = 1.0
b = 2.0
h = 0.1
y0 = 0.3011687

N = int((b - a) / h)

x = [a + k * h for k in range(N + 1)]

y_euler = [0.0] * (N + 1)
y_improved = [0.0] * (N + 1)
y_runge = [0.0] * (N + 1)

y_euler[0] = y0
y_improved[0] = y0
y_runge[0] = y0


for k in range(N):
    y_euler[k + 1] = y_euler[k] + h * f(x[k], y_euler[k])


for k in range(N):
    x_half = x[k] + h / 2
    y_half = y_improved[k] + h / 2 * f(x[k], y_improved[k])

    y_improved[k + 1] = y_improved[k] + h * f(x_half, y_half)

for k in range(N):
    K1 = h * f(x[k], y_runge[k])

    K2 = h * f(x[k] + h / 2, y_runge[k] + K1 / 2)

    K3 = h * f(x[k] + h / 2, y_runge[k] + K2 / 2)

    K4 = h * f(x[k] + h, y_runge[k] + K3)

    delta_y = (K1 + 2 * K2 + 2 * K3 + K4) / 6

    y_runge[k + 1] = y_runge[k] + delta_y

print(f"{'x':>6} {'Эйлер':>14} {'Улучш. Эйлер':>16} "
      f"{'Рунге-Кутта':>14} {'Точное':>14}")

for k in range(N + 1):
    print(f"{x[k]:6.1f} "
          f"{y_euler[k]:14.8f} "
          f"{y_improved[k]:16.8f} "
          f"{y_runge[k]:14.8f} "
          f"{y_exact(x[k]):14.8f}")



