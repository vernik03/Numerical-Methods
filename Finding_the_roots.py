from cmath import sin
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

class function3:
    def __init__(self):
        self.a = -4.0

    f = lambda self, x: x**2 + np.sin(x) - 12.0*x - 0.25

class function4:
    def __init__(self):
        self.a = 3.0
        self.b = 1.0
        
    f = lambda self, x: self.a*x + np.cos(x) + self.b
    df = lambda self, x: self.a - np.sin(x)
    phi = lambda self, x: (-np.cos(x)-self.b)/self.a

class result:
    def __init__(self, x, i):
        self.x = x
        self.i = i
   
class Calculator:
    def __init__(self, e=0.0001, segment = (-3,0), start = -5, func = function4()):
        self.e = e

        self.start = start
        self.segment = segment

        self.func = func
        
        self.real_x = np.linspace(self.segment[0], self.segment[1], 100)
        self.real_y = self.func.f(self.real_x)
        self.der_y = self.func.df(self.real_x)
        self.phi_y = []

    def dichotomy_method(self):
        a, b = self.segment
        i = 0
        while abs(b-a) > self.e:
            x = (a+b)/2
            i+=1
            if abs(self.func.f(x)) <= self.e:
                return result(x, i)
            if self.func.f(a)*self.func.f(x) < 0:
                b = x
            else:
                a = x
        plt.scatter(x, 0, color='orange', s=40, marker='o')
        return result(x, i)

    def newton_method(self):
        x = self.func.f(self.start)/self.func.df(self.start)
        x_prev = 0
        i = 0
        while abs(x_prev - x) > self.e:
            x_prev = x
            x -= self.func.f(x)/self.func.df(x)
            i += 1
        plt.scatter(x, 0, color='orange', s=40, marker='o')
        return result(x, i)

    def relaxation_method(self):
        x = (self.segment[0] + self.segment[1]) / 2
        t = -2/(self.der_y.max() + self.der_y.min())
        i = 0
        while abs(self.func.f(x)) > self.e:
            x += t * self.func.f(x)
            i += 1
        plt.scatter(x, 0, color='orange', s=40, marker='o')
        return result(x, i)

    def iteration_method(self):
        self.phi_y = self.func.phi(self.real_x)
        x = (self.segment[0] + self.segment[1]) / 2
        x_prev = 0
        i = 0
        while abs(x_prev - x) > self.e:
            x_prev = x
            x = self.func.phi(x)
            i += 1
        plt.scatter(x, 0, color='orange', s=40, marker='o')
        return result(x, i)

    def draw(self):
        plt.plot(self.real_x, self.real_y, label='f(x)')
        plt.plot(self.real_x, self.der_y, label='df(x)')
        if len(self.phi_y) > 0:
            plt.plot(self.real_x, self.phi_y, label='phi(x)')
        plt.grid(True)
        plt.show()

def print_result(result):
    print("Результат: " + str(result.x) + "\tКількость кроків: " + str(result.i))

if __name__ == '__main__':
    calculator = Calculator()
    print_result(calculator.dichotomy_method())
    print_result(calculator.newton_method())
    print_result(calculator.relaxation_method())
    print_result(calculator.iteration_method())
    calculator.draw()
    

