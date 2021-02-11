import numpy as np

# kochenderfer/wheeler, algorithm B.1
def ackley(x, a=20, b=0.2, c=2*np.pi):
    d = len(x)
    return -a*np.exp(-b*np.sqrt(np.sum(x**2)/d)) \
           - np.exp(np.sum(np.cos(c*x))/d) + a + np.exp(1)

# x is 2d
def booth(x):
    return (x[0] + 2*x[1] - 7)**2 + (2*x[0] + x[1] - 5)**2

# x is 2d
def branin(x, a=1, b=5.1/(4*np.pi**2), c=5/np.pi, r=6, s=10, t=1/(8*np.pi)):
    return a*(x[1] - b*x[0]**2 + c*x[0] - r)**2 + s*(1-t)*np.cos(x[0]) + s

# x is 2d
def flower(x, a=1, b=1, c=4):
    return a*np.linalg.norm(x) + b*np.sin(c*np.arctan2(x[1], x[0]))


