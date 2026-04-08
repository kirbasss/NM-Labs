import numpy as np
import matplotlib.pyplot as plt

# 1) x2 = 8 / (x1^2 + 4)
def curve1(x):
    return 8.0 / (x**2 + 4.0)

# 2) (x1 - 1)^2 + (x2 - 1)^2 = 4
#    x2 = 1 ± sqrt(4 - (x1 - 1)^2)
def circle_upper(x):
    inside = 4.0 - (x - 1.0)**2
    return 1.0 + np.sqrt(np.maximum(inside, 0.0))

def circle_lower(x):
    inside = 4.0 - (x - 1.0)**2
    return 1.0 - np.sqrt(np.maximum(inside, 0.0))

# Диапазон для первой кривой
x_all = np.linspace(-3.0, 5.0, 1200)

# Диапазон для окружности: (x - 1)^2 <= 4  =>  x in [-1, 3]
x_circle = np.linspace(-1.0, 3.0, 800)

y1 = curve1(x_all)
y2u = circle_upper(x_circle)
y2l = circle_lower(x_circle)

plt.figure(figsize=(8, 6))

plt.plot(x_all, y1, label=r"$x_2=\frac{8}{x_1^2+4}$")
plt.plot(x_circle, y2u, label=r"$(x_1-1)^2+(x_2-1)^2=4$ (верхняя ветвь)")
plt.plot(x_circle, y2l, label=r"$(x_1-1)^2+(x_2-1)^2=4$ (нижняя ветвь)")

plt.axhline(0, linewidth=1)
plt.axvline(0, linewidth=1)
plt.grid(True)

plt.xlim(-3, 5)
plt.ylim(-2, 4)

plt.xlabel(r"$x_1$")
plt.ylabel(r"$x_2$")
plt.title("Графическое решение системы нелинейных уравнений")
plt.legend()

plt.tight_layout()
plt.savefig("plot_2_2.png", dpi=150)
plt.show()