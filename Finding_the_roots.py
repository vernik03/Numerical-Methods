import random
from traceback import format_tb
import numpy as np
import matplotlib.pyplot as plt


class system:

    f1 = lambda arr: arr[0]**2 / arr[1]**2 - np.cos(arr[1]) - 2
    dxf1 = lambda arr: (2*arr[0])/arr[1]**2
    dyf1 = lambda arr: (-2*arr[0]**2)/arr[1]**3 + np.sin(arr[1])

    f2 = lambda arr: arr[0]**2 + arr[1]**2 - 6
    dxf2 = lambda arr: 2*arr[0]
    dyf2 = lambda arr: 2*arr[1]

class Calculator:


    def __init__(self, eps, system):
        self.system = system
        self.e = eps

    def jacobian(self, f, x, y):
        return np.array([[f.dxf1([x, y]), f.dyf1([x, y])], [f.dxf2([x, y]), f.dyf2([x, y])]])

    def newton(self, x0, y0):
        x = x0
        y = y0
        while True:
            A = self.jacobian(self.system, x, y)
            F = np.array([self.system.f1([x, y]), self.system.f2([x, y])])
            delta = np.linalg.solve(A, F)
            x -= delta[0]
            y -= delta[1]
            if np.linalg.norm(delta) <= self.e:
                return x, y

    def draw(self, a, b):
        x = np.linspace(a, b, 100)
        y = np.linspace(a, b, 100)
        X, Y = np.meshgrid(x, y)
        
        x_y = []
        for i in range(len(x)):
            for j in range(len(y)):
                x_y += [[x[i], y[j]]]
        print(x_y[0])
        Z1 = np.array([self.system.f1(x_y[i]) for i in range(len(x_y))])
        Z1 = Z1.reshape(X.shape)
        Z2 = np.array([self.system.f2(x_y[i]) for i in range(len(x_y))])
        Z2 = Z2.reshape(X.shape)
        plt.contour(X, Y, Z1, [0], colors='red')
        plt.contour(X, Y, Z2, [0], colors='blue')
        plt.grid()
        plt.show()

if __name__ == '__main__':
    calc = Calculator(0.0001, system)
    print(calc.newton(1, 1))
    # # print(calc.newton(5, 5))
    # # print(calc.newton(-0.6, -0.3))
    # # print(calc.newton(1.0, 1.0))
    calc.draw(-5,5)

