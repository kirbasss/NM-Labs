import sys
import numpy as np
import matplotlib.pyplot as plt


def parse_output_file(filename):
    x = None
    y = None
    c1 = None
    c2 = None

    with open(filename, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue

            key = parts[0]
            if key == "X":
                x = np.array(list(map(float, parts[1:])))
            elif key == "Y":
                y = np.array(list(map(float, parts[1:])))
            elif key == "C1":
                c1 = np.array(list(map(float, parts[1:])))
            elif key == "C2":
                c2 = np.array(list(map(float, parts[1:])))

    if x is None or y is None or c1 is None or c2 is None:
        raise ValueError("Не удалось прочитать PLOT_DATA из файла результата")

    return x, y, c1, c2


def eval_poly(coeffs, x):
    result = np.zeros_like(x, dtype=float)
    power = np.ones_like(x, dtype=float)
    for a in coeffs:
        result += a * power
        power *= x
    return result


def main():
    if len(sys.argv) < 2:
        print("Использование: python plot_3_3.py output/3.3.txt")
        return

    filename = sys.argv[1]
    x, y, c1, c2 = parse_output_file(filename)

    x_min = np.min(x) - 0.5
    x_max = np.max(x) + 0.5
    xs = np.linspace(x_min, x_max, 500)

    y1 = eval_poly(c1, xs)
    y2 = eval_poly(c2, xs)

    plt.figure(figsize=(9, 6))
    plt.plot(xs, y1, label="МНК, степень 1")
    plt.plot(xs, y2, label="МНК, степень 2")
    plt.scatter(x, y, label="Табличные данные")

    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Приближение функции методом наименьших квадратов")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig("plot_3_3.png", dpi=150, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()