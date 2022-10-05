import math
import copy
from enum import Enum

class Matrix(Enum):
    RANDOM = 1
    HILBERT = 2

class Method(Enum):
    LU = 1
    JACOBI = 2
    SEIDEL = 3
   
class Calculator:
    def __init__(self, A, b, e):
        self.A = A
        self.b = b
        self.e = e

    def print(self, type, toOutput):
        if type == Method.LU:
            print('L:')
            self.matrix(toOutput[0])
            print('U:')
            self.matrix(toOutput[1])
            print()
        elif type == Method.JACOBI:
            print('Jacobi:')
            self.list(toOutput)
            print()
        elif type == Method.SEIDEL:
            print('Seidel:')
            self.list(toOutput)
            print()
        else:
            print('Error')


    def maxDigit(self, matrix):
        max = 0
        for row in matrix:
            for elem in row:
                if int(elem) > 0:
                    digits = int(math.log10(int(elem)))+1
                elif int(elem) == 0:
                    digits = 1
                    if elem < 0:
                        digits = 2
                elif int(elem) < 0:
                    digits = int(math.log10(-int(elem)))+2
                if digits > max:
                    max = digits
        if max < 2:
            max = 2
        return max

    def matrix(self, matrix):
        digits = int(math.log10(int(1/self.e)))+1
        for row in matrix:
            for elem in row:
                separator = '  ' if elem >= 0 else ' '
                if int(elem) > 0:
                    sep_digits = int(math.log10(int(elem)))+1
                elif int(elem) == 0:
                    sep_digits = 1
                    if elem < 0:
                        sep_digits = 1
                elif int(elem) < 0:
                    sep_digits = int(math.log10(-int(elem)))+1
                for i in range(self.maxDigit(matrix) - sep_digits) :
                    separator += ' '
                print(separator + '{:.{}f}'.format(round(elem, digits), digits) , end='')  
            print()
        

    def list(self, list):
        digits = int(math.log10(int(1/self.e)))+1
        for i in range(len(list)):
            print('x' + str(i) + ' = ' + '{:.{}f}'.format(round(list[i], digits), digits))    
            
            

    def randMatrix(self, n, isCorect = 0):
        import random
        self.A = [[random.randint(1, 100) for i in range(n)] for j in range(n)]
        self.b = [random.randint(1, 100) for i in range(n)]
        if isCorect:
            for i in range(n):
                self.A[i][i] = sum(self.A[i]) + 1
    

    def HilbertMatrix(self, n):
        self.A = [[0] * n for i in range(n)]
        self.b = [1] * n      
        for i in range(n):
            for j in range(n):
                self.A[i][j] = 1/(i + j + 1)


    def isCorrectArray(self):
        for row in range(0, len(self.A)):
            if( len(self.A[row]) != len(self.b) ):
                print('Не соответствует размерность')
                return False

        for row in range(0, len(self.A)):
            if( self.A[row][row] == 0 ):
                print('Нулевые элементы на главной диагонали')
                return False
        return True

    # Решение систем линейныъ уравнений методом LU-разложения с вібором главного элемента
    def LU(self, n = 0, type = Matrix.RANDOM):
        if n > 0:
            if type == Matrix.RANDOM:
                self.randMatrix(n)
            elif type == Matrix.HILBERT:
                self.HilbertMatrix(n)
        n = len(self.A)
        L = [[0] * n for i in range(n)]
        U = [[0] * n for i in range(n)]
        for i in range(n):
            L[i][i] = 1
        for i in range(n):
            for j in range(n):
                if i <= j:
                    s1 = sum(U[k][j] * L[i][k] for k in range(i))
                    U[i][j] = self.A[i][j] - s1
                if i > j:
                    s2 = sum(U[k][j] * L[i][k] for k in range(j))
                    L[i][j] = (self.A[i][j] - s2) / U[j][j]
        for i in range(n):
            for i in range(j + 1):
                s1 = sum(U[k][j] * L[i][k] for k in range(i))
                U[i][j] = self.A[i][j] - s1
            for i in range(j, n):
                s2 = sum(U[k][j] * L[i][k] for k in range(j))
                L[i][j] = (self.A[i][j] - s2) / U[j][j]
        return L, U
    
    # Решение систем линейныъ уравнений методом Якоби
    def Jacobi(self, n = 0, type = Matrix.RANDOM):
        if n > 0:
            if type == Matrix.RANDOM:
                self.randMatrix(n, 1)
            elif type == Matrix.HILBERT:
                self.HilbertMatrix(n)
        if self.isCorrectArray():
            n = len(self.A)
            x = [1] * n
            x1 = [1] * n
            while True:
                for i in range(n):
                    s = 0
                    for j in range(n):
                        if i != j:
                            s += self.A[i][j] * x[j]
                    x1[i] = (self.b[i] - s) / self.A[i][i]
                if math.sqrt(sum([(x1[i] - x[i]) ** 2 for i in range(n)])/sum([(x1[i]) ** 2 for i in range(n)])) < self.e:
                    break
                x = copy.deepcopy(x1)
            return x
        else:
            return None

        # Решение систем линейныъ уравнений методом Зейделя
    def Seidel(self, n = 0, type = Matrix.RANDOM):
        if n > 0:
            if type == Matrix.RANDOM:
                self.randMatrix(n, 1)
            elif type == Matrix.HILBERT:
                self.HilbertMatrix(n)
        if self.isCorrectArray():
            n = len(self.A)
            x = [1] * n
            x1 = [1] * n
            while True:
                for i in range(n):
                    s = 0
                    for j in range(i):
                        s += self.A[i][j] * x1[j]
                    for j in range(i + 1, n):
                        s += self.A[i][j] * x[j]
                    x1[i] = (self.b[i] - s) / self.A[i][i]
                if math.sqrt(sum([(x1[i] - x[i]) ** 2 for i in range(n)])/sum([(x1[i]) ** 2 for i in range(n)])) < self.e:
                    break
                x = copy.deepcopy(x1)
            return x
        else:
            return None


if __name__ == '__main__':
    A = [[10, 1, -1],
     [1, 10, -1],
     [-1, 1, 10]]
     
    b = [11, 10, 10]
    calc = Calculator(A, b, 0.0001)
    calc.print(Method.LU, calc.LU(3, Matrix.RANDOM))
    calc.print(Method.JACOBI, calc.Jacobi(3, Matrix.RANDOM))
    calc.print(Method.SEIDEL, calc.Seidel(3, Matrix.RANDOM))
    
    
    

