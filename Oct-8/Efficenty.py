import numpy as np
import pylab as plt


def rejection(target, proposal):
    x_acpt = []
    reject = 0
    accept = 0
    while len(x_acpt) < 1000:
        ## First choose a random x
        x = np.random.uniform(0, 1)
        u = np.random.uniform(0, proposal(x))
        ## Calculate the probability of accepting x
        accept_prob = target(x) / proposal(x)
        if accept_prob >= u:
            ## If the accept probability is greater than or equal to a random number between 0 and 1, accept x
            accept += 1
            x_acpt.append(x)
        else:
            reject += 1
    a = f"Number of accepted: {accept}, Number of rejected: {reject}"
    r = 1 - reject/accept
    print(a)
    print(f"ratio of rejection/accept: {r}")

    return a

if __name__ == "__main__":
    print(f'hello world')
    def p1(x):
        return np.exp(x) / (np.exp(1) - 1)


    def p2(x):
        ## This is the line
        y_0 = p1(0) + 0.1
        return x + y_0


    def p3(x):
        ## This is the x^2 thing
        y_0 = p1(0) + 0.1
        return x ** 2 + y_0
    print(f'numbers for the line')
    a = rejection(p1, p2)
    print(f"\nnumbers for the x^2 function")
    a = rejection(p1, p3)
