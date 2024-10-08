import numpy as np
import pylab as plt

def p1(x,):
    return np.exp(x)/(np.exp(1) - 1)

def p2(x,):
    ## This is the line
    y_0 = p1(0)
    return x + y_0

def p3(x):
    ## This is the x^2 thing
    y_0 = p1(0) + 0.1
    return x**2  + y_0


