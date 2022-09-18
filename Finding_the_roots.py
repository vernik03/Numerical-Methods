import math
import struct
from time import sleep
import numpy as np
import matplotlib.pyplot as plt
import threading

class function1:
    def __init__(self):
        self.a = 5.0
        self.b = -2.0
        self.c = -2.0/5.0

    f = lambda self, x: self.a*x**3 + self.b*np.sin(x)*x**2 + self.c
    df = lambda self, x: self.a*3*x**2 + self.b*np.sin(x)*2*x + self.b*5*np.cos(x)*x

class function2:
    def __init__(self):
        self.a = -4.0

    f = lambda self, x: x**2 + self.a
    df = lambda self, x: 2*x

class Calculator:
    def __init__(self, e=0.00001, segment = (-2.0,1.1), start = 5, func = function1()):
        self.e = e

        self.start = start
        self.segment = segment

        self.func = func
        
        self.real_x = np.linspace(self.segment[0], self.segment[1], 100)
        self.real_y = self.func.f(self.real_x)
        self.der_y = self.func.df(self.real_x)

    def dichotomy_method(self):
        a, b = self.segment
        while abs(b-a) > self.e:
            x = (a+b)/2
            if abs(self.func.f(x)) <= self.e:
                return x
            if self.func.f(a)*self.func.f(x) < 0:
                b = x
            else:
                a = x
        return x

    def newton_method(self):
        x = self.func.f(self.start)/self.func.df(self.start)
        while abs(self.func.f(x)) > self.e:
            x = x - self.func.f(x)/self.func.df(x)
        return x

    def relaxation_method(self):
        x = (self.segment[0] + self.segment[1]) / 2
        while abs(self.func.f(x)) > self.e:
            t = -1/self.func.df(x)
            x += t * self.func.f(x)
        return x

    def draw(self):
        plt.plot(self.real_x, self.real_y, label='f(x)')
        plt.plot(self.real_x, self.der_y, label='df(x)')
        plt.grid(True)
        plt.show()


if __name__ == '__main__':
    calculator = Calculator()
    print(calculator.dichotomy_method())
    print(calculator.newton_method())
    print(calculator.relaxation_method())
    calculator.draw()
    