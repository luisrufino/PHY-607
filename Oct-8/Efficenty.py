import numpy as np
import pylab as plt


def rejection(target, proposal):
    x_acpt = []
    reject = 0
    accept = 0
    while len(x_acpt) < 1000:
        ## First choose a random x
        x = np.random.uniform(0, 1)
        u = np.random.uniform(0, 1)
        ## Calculate the probability of accepting x
        accept_prob = target(x) / proposal(x)
        if accept_prob >= u:
            ## If the accept probability is greater than or equal to a random number between 0 and 1, accept x
            accpet += 1
            x_acpt.append(x)
        else:
            reject += 1
    a = f"accept reject ratio: {accept/reject}"
    print(a)
    return a

if __name__ == "__main":
    print(f'hello world')
    def p1(x, ):
        return np.exp(x) / (np.exp(1) - 1)


    def p2(x, ):
        ## This is the line
        y_0 = p1(0) + 0.1
        return x + y_0


    def p3(x):
        ## This is the x^2 thing
        y_0 = p1(0) + 0.1
        return x ** 2 + y_0

    a = rejection(p1, p2)
    print("hello world")
