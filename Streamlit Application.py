import streamlit as st
import numpy as np
from sympy import symbols, Eq, solve

def bisection_fn(func_str, a, b, max_iterations):
    func = lambda x: eval(func_str)
    count = 0

    while count <= max_iterations:
        c = (a + b) / 2
        a1 = func(a)
        c1 = func(c)
        b1 = func(b)

        if c1 * b1 < 0:
            a = c
        if a1 * c1 < 0:
            b = c

        count = count + 1

    return c

def least_squares_approximation(x, y):
    n = len(x)
    sumX = np.sum(x)
    sumY = np.sum(y)
    sumXY = np.sum(x * y)
    sumX_square = np.sum(np.square(x))

    a, b = symbols('a b')
    eq1 = Eq(n * a + b * sumX, sumY)
    eq2 = Eq(a * sumX + b * sumX_square, sumXY)

    result = solve((eq1, eq2), (a, b))
    a_value, b_value = result[a], result[b]
    return a_value, b_value

def lagrange_interpolation(x, y, xp):
    n = len(x)
    yp = 0

    for i in range(n):
        p = 1

        for j in range(n):
            if i != j:
                p = p * (xp - x[j]) / (x[i] - x[j])

        yp = yp + p * y[i]

    return yp

def simpson13(func, x0, xn, n):
    h = (xn - x0) / n
    integration = func(x0) + func(xn)

    for i in range(1, n):
        k = x0 + i * h
        if i % 2 == 0:
            integration += 2 * func(k)
        else:
            integration += 4 * func(k)

    integration = integration * h / 3
    return integration

def simpson38(func, x0, xn, n):
    h = (xn - x0) / n
    integration = func(x0) + func(xn)

    for i in range(1, n):
        k = x0 + i * h
        if i % 3 == 0:
            integration = integration + 2 * func(k)
        else:
            integration = integration + 3 * func(k)

    integration = integration * 3 * h / 8
    return integration

def birge_vieta_solver(coefficients, initial_approx_root, num_iterations):
    p0 = initial_approx_root

    count = 1
    while count <= num_iterations:
        n = len(coefficients) - 1
        b = [0] * (n + 1)
        c = [0] * (n + 1)

        b[0] = coefficients[0]
        c[0] = coefficients[0]

        for i in range(1, n + 1):
            b[i] = coefficients[i] + p0 * b[i - 1]
            c[i] = b[i] + p0 * c[i - 1]

        p1 = p0 - (b[n] / c[n - 1])

        st.write("Iteration:", count)
        st.write("New root is:", p1)

        p0 = p1
        count += 1

    return p0

def rk4_solver(x0, y0, xn, n, equation):
    # Function to be solved
    f = lambda x, y: eval(equation)

    # Calculating step size
    h = (xn - x0) / n

    st.write('\n--------SOLUTION--------')
    st.write('-------------------------')
    st.write('x0\ty0\tyn')
    st.write('-------------------------')
    
    for i in range(n):
        k1 = h * f(x0, y0)
        k2 = h * f(x0 + h / 2, y0 + k1 / 2)
        k3 = h * f(x0 + h / 2, y0 + k2 / 2)
        k4 = h * f(x0 + h, y0 + k3)
        k = (k1 + 2 * k2 + 2 * k3 + k4) / 6
        yn = y0 + k
        st.write('%.4f\t%.4f\t%.4f' % (x0, y0, yn))
        st.write('-------------------------')
        y0 = yn
        x0 = x0 + h

    st.write('\nAt x=%.4f, y=%.4f' % (xn, yn))


def main():
    st.title("Numerical Methods Application")

    method = st.selectbox("Select Method", ["Bisection Method", "Least Squares Approximation", "Lagrange Interpolation","Simpson's 1/3 Method", "Simpson's 3/8 Method", "Birge-Vieta Method","4th Order Runge-Kutta"])
    
    if method == "Bisection Method":
        func_str = st.text_input("Enter the function (in terms of x):")
        a = st.number_input("Enter value of a:")
        b = st.number_input("Enter value of b:")
        max_iterations = st.number_input("Enter maximum number of iterations:")

        if st.button("Calculate"):
            root = bisection_fn(func_str, a, b, max_iterations)
            st.write(f"The final value of the root after all iterations is: {root:.6f}")

    elif method == "Least Squares Approximation":
        n = int(st.number_input("Enter the Total number of Data Points:"))
        x = np.zeros(n)
        y = np.zeros(n)

        st.write("Enter the values of x")
        for i in range(n):
            v = st.number_input(f"x{i}:", key=f"x_{i}")
            x[i] = v

        st.write("Enter the values of y")
        for i in range(n):
            v = st.number_input(f"y{i}:", key=f"y_{i}")
            y[i] = v

        if st.button("Calculate Least Squares Approximation"):
            result_a, result_b = least_squares_approximation(x, y)
            st.write(f"The result of least squares approximation: a = {result_a.evalf():.6f}, b = {result_b.evalf():.6f}")

    elif method == "Lagrange Interpolation":
        n = int(st.number_input("Enter the number of data points:"))
        x = np.zeros(n)
        y = np.zeros(n)

        st.write("Enter data values:")
        for i in range(n):
            x[i] = st.number_input(f'x[{i}]:', key=f'x_{i}')
            y[i] = st.number_input(f'y[{i}]:', key=f'y_{i}')

        xp = st.number_input("Enter interpolation point:")
        
        if st.button("Calculate Lagrange Interpolation"):
            interpolated_value = lagrange_interpolation(x, y, xp)
            st.write(f"Interpolated value at {xp:.3f} is {interpolated_value:.3f}.")
            
    elif method == "Simpson's 1/3 Method":
        func_str = st.text_input("Enter the function (in terms of x):")
        lower_limit = st.number_input("Enter lower limit of integration:")
        upper_limit = st.number_input("Enter upper limit of integration:")
        sub_interval = int(st.number_input("Enter number of sub intervals:"))
        
        if st.button("Calculate Simpson's 1/3 Method"):
            func = lambda x: eval(func_str)
            result = simpson13(func, lower_limit, upper_limit, sub_interval)
            st.write(f"Integration result: {result:.6f}")
            
    elif method == "Simpson's 3/8 Method":
        func_str = st.text_input("Enter the function (in terms of x):")
        lower_limit = st.number_input("Enter lower limit of integration:")
        upper_limit = st.number_input("Enter upper limit of integration:")
        sub_interval = int(st.number_input("Enter number of sub intervals:"))

        if st.button("Calculate Simpson's 3/8 Method"):
            func = lambda x: eval(func_str)
            result = simpson38(func, lower_limit, upper_limit, sub_interval)
            st.write(f"Integration result: {result:.6f}")
            
    elif method == "Birge-Vieta Method":
        equation_str = st.text_input("Enter the coefficients of the polynomial equation (e.g., '1 1 5 4 4'):")
        coefficients = list(map(float, equation_str.split()))
        initial_approx_root = st.number_input("Enter the initial approximate root p0:")
        num_iterations = st.number_input("Enter the number of iterations:")

        if st.button("Calculate Birge-Vieta Method"):
            root = birge_vieta_solver(coefficients, initial_approx_root, num_iterations)
            st.write(f"The final value of the root after all iterations is: {root:.6f}")
            
    elif method == "4th Order Runge-Kutta":
        x0 = st.number_input('Enter initial x (x0): ')
        y0 = st.number_input('Enter initial y (y0): ')
        xn = st.number_input('Enter final x (xn): ')
        n = int(st.number_input('Enter number of steps (n): '))
        equation_str = st.text_input('Enter the function (in terms of x and y): ')

        if st.button("Calculate 4th Order Runge-Kutta"):
            rk4_solver(x0, y0, xn, n, equation_str)
            

if __name__ == "__main__":
    main()
