import random
from traceback import format_tb
import numpy as np
import matplotlib.pyplot as plt
import time


class system:
    graph = "https://www.geogebra.org/calculator/qp7ynmja"
    n = 0

    f1 = lambda arr: arr[0]**2 / arr[1]**2 - np.cos(arr[1]) - 2
    dxf1 = lambda arr: (2*arr[0])/arr[1]**2
    dyf1 = lambda arr: (-2*arr[0]**2)/arr[1]**3 + np.sin(arr[1])

    f2 = lambda arr: arr[0]**2 + arr[1]**2 - 6
    dxf2 = lambda arr: 2*arr[0]
    dyf2 = lambda arr: 2*arr[1]

class system2:
    def __init__(self, n):
        self.n = n
        
    def f(self, arr, i):
        result = 0
        for j in range(len(arr)):
            if i == j:
                result += arr[j]**3.0 - (j+1)**3.0
            else:
                result += arr[j]**2.0 - (j+1)**2.0
        return result
    
    def df(self, arr, i): 
        result = []
        for j in range(len(arr)):
            if i == j:
                result.append(3.0*arr[j]**2.0)
            else:
                result.append(2.0*arr[j])  
        return result
        

class Calculator:
    def __init__(self, eps, system):
        self.system = system
        self.e = eps

    def jacobian(self, f, arr):
        # print("n = ", f.n)
        if f.n > 0:
            for i in range(f.n):
                if i == 0:
                    J = np.array([f.df(arr, i)])
                else:
                    J = np.append(J, [f.df(arr, i)], axis=0)
            return J
        else:
            return np.array([[f.dxf1(arr), f.dyf1(arr)], [f.dxf2(arr), f.dyf2(arr)]])

    def newton(self, arr): 
        while True:
            A = self.jacobian(self.system, arr)
            print("A = ", A)
            if self.system.n > 0:
                F = np.array([self.system.f(arr, i) for i in range(self.system.n)])
            else:
                F = np.array([self.system.f1(arr), self.system.f2(arr)]) 
            # print("F = ", F)
            delta = np.linalg.solve(A, F)
            for i in range(len(arr)):
                arr[i] -= delta[i]
            if np.linalg.norm(delta) <= self.e:
                return arr

    def newton_mod(self, arr):
        A = self.jacobian(self.system, arr)
        print("A = ", A)
        while True:          
            if self.system.n > 0:
                F = np.array([self.system.f(arr, i) for i in range(self.system.n)])
            else:
                F = np.array([self.system.f1(arr), self.system.f2(arr)]) 
            delta = np.linalg.solve(A, F)
            for i in range(len(arr)):
                arr[i] -= delta[i]
            print(arr)
            time.sleep(1)
            if np.linalg.norm(delta) <= self.e:
                return arr

    def simple_iteration(self, arr):
        while True:
            if self.system.n > 0:
                for i in range(self.system.n):
                    arr[i] = self.system.f(arr, i)
            else:
                print(arr)
                time.sleep(1)
                arr[0] = self.system.f1(arr)
                arr[1] = self.system.f2(arr)  
            if np.linalg.norm(arr) <= self.e:
                return arr


    def draw(self, a, b):
        if self.system.n > 2:
            return print("Can't draw in dimensions higher than 2")
        x = np.linspace(a, b, 100)
        y = np.linspace(a, b, 100)
        Y, X = np.meshgrid(x, y)        
        x_y = []
        for i in range(len(x)):
            for j in range(len(y)):
                x_y += [[x[i], y[j]]]
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
    print(calc.newton([1, 1]))
    print(calc.newton_mod([1, 1]))
    # print(calc.simple_iteration([1, 1]))
    # print(calc.newton([-0.6, -0.3]))
    # print(calc.newton([1, 2.5]))
    # print(calc.newton([-1.0, -2.5]))
    # calc.draw(-20,20)
    
    calc2 = Calculator(0.0001, system2(5))
    print(calc2.newton([1, 1, 1, 1, 1]))
    print(calc2.newton_mod([1, 1, 1, 1, 1]))

