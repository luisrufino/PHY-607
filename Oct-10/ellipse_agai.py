import numpy as np
import pylab as plt
def arc_length_ellipse(a, b, n_samples,Sweeps):
    arc_length_tot = [] ## Array that contains the
    for n in range(1, n_samples):
        delta = 0.01
        arc_length = [] ## arc length for n samples
        for i in range(Sweeps):
            x = np.random.uniform(-a, a, n)
            y = np.random.uniform(-b, b, n)
            plus_ellipse = (x/(a+delta))**2 + (y/(b+delta))**2 < 1
            plus_area = a * b * 4 * np.mean(plus_ellipse)
            minus_ellipse = (x/(a-delta))**2 + (y/(b-delta))**2 < 1
            minus_area = a * b * 4 * np.mean(minus_ellipse)
            ## Now compare the areas, and this should give me some length unit which is related to the arc length
            ## Now compare
            arc = (plus_area - minus_area)/(2*delta)
            arc_length.append(arc)
        std_arc = np.std(arc_length)
        mean_arc = np.mean(arc_length)
        arc_length_tot.append([mean_arc,std_arc])

    return arc_length_tot

def random_int(target, proposal):
    start = 0
    end = 1
    n = 1_000
    if proposal == None:
        x_rand = np.random.uniform(0, 1, n)
        f_rand = target(x_rand)
        res = np.sum(f_rand)/n
    else:
        x_rand = np.random.uniform(proposal(start), proposal(end), n)
        f_rand = target(x_rand)
        #cont = 1/np.sum(proposal(x_rand))
        prob_x = 1/proposal(x_rand)
        res = np.sum(f_rand*prob_x)
    print(f"Sum: {res}")
    return res

if __name__ == "__main__":
    def f1(x):
        return -x**2 + 10
    def f2(x):
        return -x + 11
    random_int(f1, f2)
    # a = 5
    # b = 2
    # n_samples = 1000
    # sweeps = 1000
    # arc_lengths = arc_length_ellipse(a, b, n_samples, sweeps)
    # arc_lengths = np.array(arc_lengths)
    # ns = np.arange(1,n_samples)
    # print(len(ns))
    # print(arc_lengths.shape)
    # plt.plot(ns, arc_lengths[:,0])
    # plt.title("arc length vs number of samples")
    # plt.xlabel('Number of Samples')
    # plt.ylabel('Mean Arc Length')
    # plt.show()