import numpy as np
import matplotlib.pyplot as plt

n = 1
sum = np.float32(0)

while True:
    tmp = np.float32(sum + np.float32(1) / np.float32(n))
    if tmp == sum:
        break
    sum = tmp
    n += 1

print(f'Stop at n = {n}, sum = {sum}')
sum32 = sum
n32 = n


n = 1
sum = np.float64(0)
eps = 6e-8

while True:
    tmp = np.float64(np.float64(1) / np.float64(n))
    if tmp <= eps * sum / 2:
        break
    sum += tmp
    n += 1
print(f'Theoretical value: n = {n}')


sum64 = np.float64(0)
for _n in range(1, n32+1):
    sum64 = np.float64(sum64 + np.float64(1) / np.float64(_n))
error = np.abs(sum64 - sum32)

print(f'Absolute error = {error}')
print(f'Relative error = {error / sum64:.3%}')


from scipy.optimize import fsolve

def f(n):
    gamma = 0.57721566490153286060651209008240243104215933593992
    return 1e-16 * (np.log(n) + gamma + 1 / (2 * n)) / 2 - 1 / n

n = int(fsolve(f, [1.])[0])
print(f'With float64, n will be {n}')
print(f'Will stop in {n * 2 / (80 * 1e9) / 3600:.3f} hours')