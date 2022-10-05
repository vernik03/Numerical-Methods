import math
   
class Calculator:
    def __init__(self):
        self.x = 0

    # Решение систем линейныъ уравнений методом LU-разложения
    def LU(self, A):
        n = len(A)
        L = [[0] * n for i in range(n)]
        U = [[0] * n for i in range(n)]
        for i in range(n):
            L[i][i] = 1
        for j in range(n):
            for i in range(j + 1):
                s1 = sum(U[k][j] * L[i][k] for k in range(i))
                U[i][j] = A[i][j] - s1
            for i in range(j, n):
                s2 = sum(U[k][j] * L[i][k] for k in range(j))
                L[i][j] = (A[i][j] - s2) / U[j][j]
        return L, U
        

        # Решение систем линейныъ уравнений методом Якоби
    def Jacobi(self, A, b, e):
        n = len(A)
        x = [0] * n
        x1 = [0] * n
        while True:
            for i in range(n):
                s = 0
                for j in range(n):
                    if i != j:
                        s += A[i][j] * x[j]
                x1[i] = (b[i] - s) / A[i][i]
            if math.sqrt(sum([(x1[i] - x[i]) ** 2 for i in range(n)])) < e:
                break
            x = x1
        return x

        # Решение систем линейныъ уравнений методом Зейделя
    def Seidel(self, A, b, e):
        n = len(A)
        x = [0] * n
        x1 = [0] * n
        while True:
            for i in range(n):
                s = 0
                for j in range(i):
                    s += A[i][j] * x1[j]
                for j in range(i, n):
                    s += A[i][j] * x[j]
                x1[i] = (b[i] - s) / A[i][i]
            if math.sqrt(sum([(x1[i] - x[i]) ** 2 for i in range(n)])) < e:
                break
            x = x1
        return x


if __name__ == '__main__':
    calculator = Calculator()
    A = [[10, -7, 0], [-3, 6, 2], [5, -1, 5]]
    print(calculator.LU(A))
    
    A = [[-0.1, 0.1, 1.1], [-0.1, 0.1, 1.0], [0.1, -0.1, 1.0]]
    b = [1.1, 1.0, 1.0]
    print(calculator.Jacobi(A, b, 0.0001))
    print(calculator.Seidel(A, b, 0.0001))
    #print(calculator.Jacobi(A, b, 0.0001))
    #print(calculator.Seidel(A, b, 0.0001))
    
    
    

