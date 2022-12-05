import numpy as np
import numpy.polynomial.polynomial as polynom
import matplotlib
import matplotlib.pyplot as plt

f1 = lambda x: np.cos(x) + x * np.sin(x)

class Interpolation:
    def __init__(self, func = f1, x_start= -10, x_end = 10, n=10, points = 1000):
        self.x = np.linspace(x_start, x_end, points)
        self.x_interval = np.linspace(x_start, x_end, n)
        self.fx_interval = f1(self.x_interval)
        self.values = np.array([self.x_interval, self.fx_interval]).T

    def calculate(self):
        self.lagrange_value = self.lagrange_interpolation()
        self.newton_value = self.newton_interpolation()
        self.spline_value = self.spline_interpolation()
        
        print(f"Lagrange polynom: {self.lagrange_value}\n")
        print(f"Newton polynom: {self.newton_value}\n")
        print(f"Spline interpolation: {self.spline_value}\n")

    def draw(self):
        plt.plot(self.x_interval, self.fx_interval,'ko')
        plt.plot(self.x, f1(self.x),'black', label='cos(x)+x*sin(x)')
        plt.plot(self.x, self.lagrange_value(self.x), 'pink', label='Lagrange')
        plt.plot(self.x, self.newton_value(self.x), 'red', label='Newton')
        plt.plot(self.x, self.get_spline_values(self.spline_value), 'purple', label='Spline')
        plt.legend()
        plt.show()

    def lagrange_interpolation(self):
        n = self.values.shape[0]
        l_n = polynom.Polynomial([0])

        for i in range(n):
            l_i = polynom.Polynomial([1])
            for j in range(n):
                if i != j:
                    l_i *= polynom.Polynomial([-self.values[j][0], 1]) / (self.values[i][0] - self.values[j][0])

            l_n += self.values[i][1] * l_i


        return l_n

    def newton_interpolation(self):
        l_n = polynom.Polynomial([self.values[0][1]])
        n = self.values.shape[0]
        diff_table = self.get_diff_table()
        poly_factor = polynom.Polynomial([1])

        for i in range(1, n):
            poly_factor *= polynom.Polynomial([-self.values[i - 1][0], 1])
            l_n += diff_table[0][i] * poly_factor

        return l_n


    def get_diff_table(self):
        n = self.values.shape[0]
        table = np.zeros([n, n])
        table[:, 0] = np.copy(self.values[:, 1])

        for i in range(1, n):
            for j in range(n - i):
                x_j = self.values[j][0]
                x_k = self.values[j + i][0]
                diff_x = x_k - x_j
                diff_f = table[j + 1][i - 1] - table[j][i - 1]
                table[j][i] = diff_f / diff_x
                
        return table


    def spline_interpolation(self):
        n = self.values.shape[0]
        intervals = self.get_h_vector()
        c_values = self.get_c_vector(intervals, n)
        s_n = self.get_s_n_for_all_intervals(intervals, c_values, n)

        return s_n


    def get_h_vector(self):
        n = self.values.shape[0] - 1
        interval = np.zeros(n)

        for i in range(n):
            interval[i] = self.values[i + 1][0] - self.values[i][0]

        return interval


    def get_c_vector(self, intervals, n):
        intervals_count = n - 1
        c_matrix = np.zeros([intervals_count - 1, intervals_count - 1])
        right_val = np.zeros(intervals_count - 1)
        c = np.zeros(n)

        for i in range(intervals_count - 1):
            if i > 0:
                c_matrix[i][i - 1] = intervals[i]

            c_matrix[i][i] = 2 * (intervals[i] + intervals[i + 1])

            if i < intervals_count - 2:
                c_matrix[i][i + 1] = intervals[i + 1]

            right_val[i] = 6 * (
                    (self.values[i + 2][1] - self.values[i + 1][1]) / intervals[i + 1] -
                    (self.values[i + 1][1] - self.values[i][1]) / intervals[i])

        c[1:n - 1] = np.linalg.solve(c_matrix, right_val)

        return c


    def get_s_n_for_all_intervals(self, intervals, c_values, n):
        intervals_count = n - 1
        s_n = np.empty([intervals_count], polynom.Polynomial)

        for i in range(1, intervals_count + 1):
            a_i = self.values[i][1]
            d_i = (c_values[i] - c_values[i - 1]) / intervals[i - 1]
            b_i = intervals[i - 1] * c_values[i] / 2 - intervals[i - 1] ** 2 * d_i / 6 + (self.values[i][1] - self.values[i - 1][1]) / intervals[i - 1]

            poly_factor = polynom.Polynomial([a_i])
            poly_factor += b_i * polynom.Polynomial([-self.values[i][0], 1])
            poly_factor += (c_values[i] / 2) * polynom.Polynomial([-self.values[i][0], 1]) ** 2
            poly_factor += (d_i / 6) * polynom.Polynomial([-self.values[i][0], 1]) ** 3

            s_n[i - 1] = poly_factor

        return s_n


    def get_spline_values(self, s_n):
        total_points = self.x.shape[0]
        total_spline = s_n.shape[0]
        s_x = np.zeros(total_points)
        number_of_spline = 0

        for i in range(total_points):
            if number_of_spline < total_spline - 1 and \
                    self.x[i] > self.x_interval[number_of_spline + 1]:
                number_of_spline += 1

            s_x[i] = polynom.polyval(self.x[i], s_n[number_of_spline].coef)

        return s_x


if __name__ == '__main__':
    calc = Interpolation()
    calc.calculate()
    calc.draw()