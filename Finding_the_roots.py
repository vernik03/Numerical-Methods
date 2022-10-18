import numpy as np
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
    def __init__(self, e, A = 0, b = 0, n = 0, type = Matrix.RANDOM):
        self.e = e
        if A != 0 and b != 0:
            self.A = A
            self.b = b
        elif n > 0:
            self.A = [[0] * n for i in range(n)]
            self.b = [1] * n     
            if type == Matrix.RANDOM:
                self.randMatrix(n, 1)
            elif type == Matrix.HILBERT:
                self.HilbertMatrix(n)

    def print(self, type, toOutput):
        if toOutput != None:
            if type == Method.LU:
                print('L:')
                self.matrix(toOutput[0])
                print('U:')
                self.matrix(toOutput[1])
                print('result:')
                self.list(toOutput[2])
            elif type == Method.JACOBI:
                print('Jacobi:')
                self.list(toOutput)
            elif type == Method.SEIDEL:
                print('Seidel:')
                self.list(toOutput)
            else:
                print('Error')
            print()
        else:
            if type == Method.LU:
                print('in LU')
            elif type == Method.JACOBI:
                print('in Jacobi')
            elif type == Method.SEIDEL:
                print('in Seidel')
            else:
                print('Error')
            print()

    def maxDigit(self, matrix):
        max = 0
        for row in matrix:
            for elem in row:
                if int(elem) > 0:
                    digits = int(np.log10(int(elem)))+1
                elif int(elem) == 0:
                    digits = 1
                    if elem < 0:
                        digits = 2
                elif int(elem) < 0:
                    digits = int(np.log10(-int(elem)))+2
                if digits > max:
                    max = digits
        if max < 2:
            max = 2
        return max

    def matrix(self, matrix):
        digits = int(np.log10(int(1/self.e)))
        for row in matrix:
            for elem in row:
                separator = '  ' if elem >= 0 else ' '
                if int(elem) > 0:
                    sep_digits = int(np.log10(int(elem)))+1
                elif int(elem) == 0:
                    sep_digits = 1
                    if elem < 0:
                        sep_digits = 1
                elif int(elem) < 0:
                    sep_digits = int(np.log10(-int(elem)))+1
                for i in range(self.maxDigit(matrix) - sep_digits) :
                    separator += ' '
                print(separator + '{:.{}f}'.format(round(elem, digits), digits) , end='')  
            print()
        
    def list(self, list):
        digits = int(np.log10(int(1/self.e)))
        for i in range(len(list)):
            print('x' + str(i+1) + ' = ' + '{:.{}f}'.format(round(list[i], digits), digits))    
                        
    def randMatrix(self, n, isCorect = 0):
        import random
        self.A = [[random.randint(1, 10) for i in range(n)] for j in range(n)]
        self.b = [random.randint(1, 10) for i in range(n)]
        if isCorect:
            for i in range(n):
                self.A[i][i] = sum(self.A[i]) + 1
    
    def HilbertMatrix(self, n):
        self.A = [[0] * n for i in range(n)]
        self.b = [1] * n      
        for i in range(n):
            for j in range(n):
                self.A[i][j] = 1/(i + j + 1)

    def isCorrectArray(self, Method = Method.LU):
        for row in range(0, len(self.A)):
            if( len(self.A[row]) != len(self.b) ):
                print('Size does not match')
                return False

        if Method == Method.JACOBI or Method == Method.SEIDEL:
            for row in range(0, len(self.A)):
                if( self.A[row][row] == 0 ):
                    print('Zero elements on the main diagonal')
                    return False
            if not self.conditionJacobi():
                print('The method does not converge')
                return False
            
        return True

    def conditionJacobi(self):
        for i in range(len(self.A)):
            sum = 0
            for j in range(len(self.A)):
                if i != j:
                    sum += abs(self.A[i][j])
            if abs(self.A[i][i]) < sum:
                return False
        return True
        
    def isSymmetric(self):
        for i in range(len(self.A)):
            for j in range(len(self.A)):
                if self.A[i][j] != self.A[j][i]:
                    return False
        return True

    def transpose(self, matrix):
        return list(map(list, zip(*matrix)))

    def det(self, A):
        if len(A) == 1:
            return A[0][0]
        elif len(A) == 2:
            return A[0][0] * A[1][1] - A[0][1] * A[1][0]
        else:
            det = 0
            for i in range(len(A)):
                det += (-1)**i * A[0][i] * self.det(self.minor(A, 0, i))
            return det
    
    def minor(self, A, i, j):
        return [row[:j] + row[j+1:] for row in (A[:i]+A[i+1:])]

    def LU(self, n = 0, type = Matrix.RANDOM):
        if n > 0:
            if type == Matrix.RANDOM:
                self.randMatrix(n)
            elif type == Matrix.HILBERT:
                self.HilbertMatrix(n)
        n_y = len(self.A)
        n_x = len(self.A[0])
        L = [[0] * n_x for i in range(n_y)]
        U = [[0] * n_x for i in range(n_y)]
        for i in range(n_y):
            L[i][i] = 1
        # decomposition of matrix
        for i in range(n_y):
            for j in range(n_x):
                if i <= j:
                    U[i][j] = self.A[i][j] - sum([L[i][k] * U[k][j] for k in range(i)])
                else:
                    L[i][j] = (self.A[i][j] - sum([L[i][k] * U[k][j] for k in range(j)])) / U[j][j]
        # lu = L+U-I
        # find solution of Ly = b
        y = [0] * n_y
        for i in range(n_y):
            y[i] = self.b[i]
            for j in range(i):
                y[i] -= L[i][j] * y[j]
            y[i] /= L[i][i]        
        # find solution of Ux = y
        x = [0] * n_y
        for i in range(n_y-1, -1, -1):
            x[i] = y[i]
            for j in range(i+1, n_y):
                x[i] -= U[i][j] * x[j]
            x[i] /= U[i][i]
        return L, U, x
        
    def findNorm(self, array):
        return max([abs(i) for i in array])       
    
    def Jacobi(self, n = 0, type = Matrix.RANDOM):
        if n > 0:
            if type == Matrix.RANDOM:
                self.randMatrix(n, 1)
            elif type == Matrix.HILBERT:
                self.HilbertMatrix(n)
        if self.isCorrectArray(Method.JACOBI):
            n = len(self.A)
            x = [1] * n
            x1 = [1] * n
            steps = 0
            while steps < 500:
                for i in range(n):
                    s = 0
                    for j in range(n):
                        if i != j:
                            s += self.A[i][j] * x[j]
                    x1[i] = (self.b[i] - s) / self.A[i][i]
                # if np.sqrt(sum([(x1[i] - x[i]) ** 2 for i in range(n)])/sum([(x1[i]) ** 2 for i in range(n)])) < self.e:
                if self.findNorm([(x1[i] - x[i]) for i in range(n)]) < self.e:
                    break 
                x = copy.deepcopy(x1)
                steps += 1
            return x
        else:
            return None

    def Seidel(self, n = 0, type = Matrix.RANDOM):
        if n > 0:
            if type == Matrix.RANDOM:
                self.randMatrix(n, 1)
            elif type == Matrix.HILBERT:
                self.HilbertMatrix(n)
        if self.isCorrectArray(Method.SEIDEL):
            n = len(self.A)
            x = [1] * n
            x1 = [1] * n
            steps = 0
            while steps < 500:
                for i in range(n):
                    s = 0
                    for j in range(i):
                        s += self.A[i][j] * x1[j]
                    for j in range(i + 1, n):
                        s += self.A[i][j] * x[j]
                    x1[i] = (self.b[i] - s) / self.A[i][i]
                if self.findNorm([(x1[i] - x[i]) for i in range(n)]) < self.e:
                    break
                x = copy.deepcopy(x1)
                steps += 1
            return x
        else:
            return None
        

if __name__ == '__main__':
    calc = Calculator( 0.01, 0, 0, 5, Matrix.RANDOM)
    print("A:")
    calc.matrix(calc.A)
    print()
    print("b:")
    calc.list(calc.b)
    print()
    print("Determinant: " + str(calc.det(calc.A)))
    print()
    calc.print(Method.LU, calc.LU())
    calc.print(Method.JACOBI, calc.Jacobi())
    calc.print(Method.SEIDEL, calc.Seidel())
    
    
    

