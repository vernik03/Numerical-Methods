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
    CHOLESKY = 4
   
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
            elif type == Method.JACOBI:
                print('Jacobi:')
                self.list(toOutput)
            elif type == Method.SEIDEL:
                print('Seidel:')
                self.list(toOutput)
            elif type == Method.CHOLESKY:
                print('Cholesky:')
                print('L:')
                self.matrix(toOutput)
                print('L-tranpose:')
                self.matrix(self.transpose(toOutput))
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
            elif type == Method.CHOLESKY:
                print('in Cholesky')
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
        digits = int(np.log10(int(1/self.e)))+1
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
        digits = int(np.log10(int(1/self.e)))+1
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

        if Method == Method.CHOLESKY:
            if not self.isSymmetric():
                print('The matrix is not symmetrical')
                return False
            if self.A != self.transpose(self.A):
                print('Matrix is not equal to transposed')
                return False
            if self.determinant(self.A) <= 0:
                print('Determinant <= 0')
                return False

        if Method == Method.JACOBI or Method == Method.SEIDEL:
            # if np.linalg.norm(self.A) >= 1:
            #     print('Не сходится')
            #     return False
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

    def determinant(self, matrix):
        if len(matrix) == 2:
            return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
        else:
            det = 0
            for i in range(len(matrix)):
                det += ((-1) ** i) * matrix[0][i] * self.determinant([row[:i] + row[i+1:] for row in (matrix[1:])])
            return det


    # найти определитель матрицы A
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
                if np.sqrt(sum([(x1[i] - x[i]) ** 2 for i in range(n)])/sum([(x1[i]) ** 2 for i in range(n)])) < self.e:
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
                if np.sqrt(sum([(x1[i] - x[i]) ** 2 for i in range(n)])/sum([(x1[i]) ** 2 for i in range(n)])) < self.e:
                    break
                x = copy.deepcopy(x1)
                steps += 1
            return x
        else:
            return None
     
    def Cholesky(self, n = 0, type = Matrix.RANDOM):
        if n > 0:
            if type == Matrix.RANDOM:
                self.randMatrix(n, 1)
            elif type == Matrix.HILBERT:
                self.HilbertMatrix(n)
        if self.isCorrectArray(Method.CHOLESKY):
            n = len(self.A)
            L = [[0] * n for i in range(n)]
            for i in range(n):
                for j in range(i + 1):
                    s = sum(L[i][k] * L[j][k] for k in range(j))
                    if i == j:
                        L[i][j] = np.sqrt(self.A[i][i] - s)
                    else:
                        L[i][j] = (self.A[i][j] - s) / L[j][j]
            return L
        else:
            return None
        
        


if __name__ == '__main__':
    # A = [[1, 2, 3], [2, 5, 5], [3, 5, 6]]
    # b = [1, 2, 3]
    A = [[10, 1, 0, 0, 0],
         [1, 30, 1, 0, 0],
         [0, 1, 20, 1, 0],
         [0, 0, 1, 20, 2],
         [0, 0, 0, 2, 10]]
    b = [1, 3, 2, 2, 3]
    calc = Calculator( 0.01, A, b, 5, Matrix.RANDOM)
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
    calc.print(Method.CHOLESKY, calc.Cholesky())
    
    
    

