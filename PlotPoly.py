import numpy as np
import matplotlib.pyplot as plt

def PlotPoly():
    # Read data from file
    data = np.loadtxt('poly.data')
    x = data[:, 0]
    y = data[:, 1]

    # Plotting
    plt.figure(figsize=(10, 6))
    plt.plot(x, y, linestyle='dashed', linewidth=2, marker='o', color='red', markersize=6)
    plt.xlim(-1.0, 1.0)
    plt.xticks(np.linspace(-1.0, 1.0, num=5))
    plt.grid(True)
    plt.xlabel("x-axis", fontsize=14)
    plt.ylabel("y-axis", fontsize=14)
    plt.title("Polynomial Plot", fontsize=16)
    plt.savefig('polynomial_plot.png', dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    PlotPoly()
