import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return 2**x - x**2 - 0.5

x = np.linspace(-3, 5, 1000)
y = f(x)

plt.figure(figsize=(8, 5))
plt.axhline(0, linewidth=1, color="black")
plt.axvline(0, linewidth=1, color="black")
plt.plot(x, y, label=r"$2^x - x^2 - 0.5$")

plt.grid(True)
plt.legend()
plt.xlabel("x")
plt.ylabel("f(x)")
plt.title("График уравнения 2^x - x^2 - 0.5 = 0")

plt.savefig("plot_2_1.png", dpi=150, bbox_inches="tight")
plt.show()