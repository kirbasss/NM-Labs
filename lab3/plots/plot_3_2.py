import numpy as np
import matplotlib.pyplot as plt

# точки
points = np.loadtxt("output/3.2.csv", delimiter=",", skiprows=1)
xp = points[:, 0]
yp = points[:, 1]

# сплайн
spline = np.loadtxt("output/3.2_spline.csv", delimiter=",", skiprows=1)
xs = spline[:, 0]
ys = spline[:, 1]

plt.figure()

# точки
plt.scatter(xp, yp, label="Исходные точки")

# сплайн
plt.plot(xs, ys, label="Кубический сплайн")

x_true = np.linspace(0, 4, 200)
y_true = np.sin(x_true)

plt.legend()
plt.grid(True)
plt.title("Кубический сплайн (задача 3.2)")

plt.savefig("plot_3_2.png", dpi=150, bbox_inches="tight")

plt.show()