import math


class Calculator:
    def __init__(self, e=0.00001, a=-10, b=10):
        #ax^2 + bx + c = 0
        #5x^3 - 2x^2*sin(x) - 2/5 = 0
        self.a = 5.0
        self.b = -2.0
        self.c = -2.0/5.0
        self.e = e
        self.segment = (a,b)

    def f(self, x):
        return self.a*x**3 + self.b*math.sin(x)*x**2 + self.c

    def f_prime(self, x):
        return self.a*3*x**2 + self.b*math.sin(x)*2*x + self.b*math.cos(x)*x**2

    def dichotomy_method(self):
        a, b = self.segment
        while abs(b-a) > self.e:
            x = (a+b)/2
            if self.f(a)*self.f(x) < 0:
                b = x
            else:
                a = x
        return x

    def newton_method(self):
        x = self.segment[0]
        while abs(self.f(x)) > self.e:
            x = x - self.f(x)/self.f_prime(x)
        return x

   


if __name__ == '__main__':
    calculator = Calculator()
    print(calculator.dichotomy_method())
    print(calculator.newton_method())
    

