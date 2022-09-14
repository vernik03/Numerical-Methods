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

    f = lambda self, x: self.a*x**3 + self.b*math.sin(x)*x**2 + self.c

    f_derivative = lambda self, x: self.a*3*x**2 + self.b*math.sin(x)*2*x + self.b*5*math.cos(x)*x
    
    f_derivative_second = lambda self, x: self.a*6*x + self.b*math.sin(x)*(x**2 - 2) + self.b*20*math.cos(x)*x

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
            x = x - self.f(x)/self.f_derivative(x)
        return x

    # def relaxation_method(self):
    #     x = self.segment[0]
    #     while abs(self.f(x)) > self.e:
    #         x = x - self.f(x)/self.f_derivative_second(x)
    #     return x



if __name__ == '__main__':
    calculator = Calculator()
    print(calculator.dichotomy_method())
    print(calculator.newton_method())

    

