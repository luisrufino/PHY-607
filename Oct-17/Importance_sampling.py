import numpy as np
import pylab as plt
from scipy import integrate

def func(x):
    return np.exp(x)
## Now I want to sample from a uniform distribution
def uniform(xmin, xmax):
    return np.random.uniform(xmin, xmax)

## Now I want to sample from some probability distribution

def probability(x): ## THis is my q bar.
    ## This is from 0 to 1
    return (np.exp(x) + 1)/(np.exp(1) - 1)

def probability1(x):
    return ((x/0.8)**2 +5)/5.52083

### so now I want to integrate func(x) with the probability distribution


def test(func, probability, xmin, xmax, nsteps=10000):
    list = []
    for i in range(nsteps):
        x_uniform = uniform(xmin, xmax)
        ## Compute the 'ratio'
        a = func(x_uniform)/probability(x_uniform)
        list.append(a)
    b = np.mean(list)
    res_prob_int = integrate.quad(probability, 0, 1)


    print(b*res_prob_int[0])

test(func, probability1, 0, 1, nsteps=100_000)

