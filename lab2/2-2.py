import numpy as np
from mpmath import diff
from scipy.optimize import root

eps = 1e-8

def damp_Newton_Raphson(f, x):
    step = 0
    while True:
        tmp = f(x) / diff(f, x)
        lambd = 1
        while lambd > 0 and np.abs(f(x - lambd * tmp)) >= np.abs(f(x)):
            lambd /= 2
        x -= lambd * tmp
        step += 1
        print(f'Step {step}: x = {x}, f(x) = {f(x)}')
        if lambd <= 0:
            break
    return x

def basc_Newton_Raphson(f, x):
    step = 0
    last_x = x
    while np.abs(f(x)) > eps or np.abs(x - last_x) > eps:
        last_x = x
        x -= f(x) / diff(f, x)
        step += 1
        print(f'Step {step}: x = {x}, f(x) = {f(x)}')
    return x

def solve(f, x0):
    print('damping Newton Raphson:')
    damp_r = damp_Newton_Raphson(f, x0)
    print('\nbasic Newton Raphson:')
    basic_r = basc_Newton_Raphson(f, x0)
    scipy_r = root(f, x0).x[0]
    print(f'\nBasic Newton Raphson: {basic_r}')
    print(f'Newton Raphson with damp: {damp_r}')
    print(f'SciPy: {scipy_r}')
    print(f'Newton error: {np.float64((basic_r - scipy_r) / scipy_r):.8%}\nNewton with damp error: {np.float64((damp_r - scipy_r) / scipy_r):.8%}\n')

solve(lambda x: x ** 3 - x - 1, 0.6)
solve(lambda x: - x ** 3 + 5 * x, 1.35)