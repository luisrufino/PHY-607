import numpy as np
import pylab as plt


def rejection(target, proposal):
    x_acpt = []
    u_acpt = []
    reject = 0
    accept = 0
    while len(x_acpt) < 10_000:
        ## First choose a random x
        x = np.random.uniform(0, 1)
        u = np.random.uniform(0, proposal(x))
        ## Calculate the probability of accepting x
        ## Now check if u is above or below the value at x for the target
        if u < target(x):
            ## If the accept probability is greater than or equal to a random number between 0 and 1, accept x
            accept += 1
            x_acpt.append(x)
            u_acpt.append(u)
        else:
            reject += 1
    a = f"Number of accepted: {accept}, Number of rejected: {reject}"
    r = 1 - reject/accept
    print(a)
    print(f"ratio of rejection/accept: {r}")

    return x_acpt, u_acpt

if __name__ == "__main__":
    print(f'hello world')
    def p1(x):
        return np.exp(x) / (np.exp(1) - 1)


    def p2(x):
        ## This is the line
        y_0 = p1(0) + 0.5
        return x + y_0


    def p3(x):
        ## This is the x^2 thing
        y_0 = p1(0) + 0.5
        return x ** 2 + y_0
    print(f'numbers for the line')
    a, b = rejection(p1, p2)
    plt.scatter(a, b)
    plt.title("proposal g(x) = Line")
    plt.show()
    print(f"\nnumbers for the x^2 function")
    a, b = rejection(p1, p3)
    plt.scatter(a, b)
    plt.title("proposal g(x) = quadratic")
    plt.show()