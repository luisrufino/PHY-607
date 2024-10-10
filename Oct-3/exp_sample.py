import numpy as np
import matplotlib.pyplot as plt

# Inverse CDF Sampling for e^x on [0, 1]
def inverse_cdf_sampling(n_samples):
    # Generate uniform samples
    uniform_samples = np.random.uniform(0, 1, n_samples)

    # Apply the inverse CDF transformation for f(x) = e^x on [0, 1]
    exponential_samples = np.log((uniform_samples * (np.exp(1) - 1)) + 1)

    return exponential_samples

# Rejection Sampling for e^x on [0, 1]
def rejection_sampling_exp(n_samples, M=np.exp(1)):
    samples = []
    while len(samples) < n_samples:
        # Sample from uniform distribution
        x_proposal = np.random.uniform(0, 1)  # Proposal distribution
        u = np.random.uniform(0, 1)  # Uniform sample for rejection

        # Target distribution (e^x) and scaled proposal
        if u <= np.exp(x_proposal) / M:
            samples.append(x_proposal)
    return np.array(samples)

# Monte Carlo Integration: Circumference of an Ellipse
def monte_carlo_ellipse_circumference(a, b, n_samples):
    # Sample theta uniformly between 0 and 2pi
    theta_samples = np.random.uniform(0, 2 * np.pi, n_samples)

    # Compute the integrand (arc length element)
    arc_lengths = np.sqrt((a * np.sin(theta_samples))**2 + (b * np.cos(theta_samples))**2)
    #print(np.mean(arc_lengths))

    # Estimate the circumference using the average arc length times the range of theta
    circumference_estimate = (2 * np.pi) * np.mean(arc_lengths)

    return circumference_estimate

def arc_length_ellipse(a, b, n_samples):
    arc_length = []
    for i in range(n_samples):
        delta = 0.0001
        x = np.random.uniform(-a, a, n_samples)
        y = np.random.uniform(-b, b, n_samples)
        plus_ellipse = (x/(a+delta))**2 + (y/(b+delta))**2 < 1
        plus_area = a * b * 4 * np.mean(plus_ellipse)
        minus_ellipse = (x/(a-delta))**2 + (y/(b-delta))**2 < 1
        minus_area = a * b * 4 * np.mean(minus_ellipse)
        ## Now compare the areas, and this should give me some length unit which is related to the arc length
        ## Now compare
        arc = (plus_area - minus_area)/(2*delta)
        #print(arc)
        arc_length.append(arc)
    std_arc = np.std(arc_length)
    mean_arc = np.mean(arc_length)

    return mean_arc, std_arc

## This is the best aproximation to the arc length of an ellipse

def approx(a, b):
    t = ((a-b)/(a+b))**2
    return np.pi*(a+b)*(1 + 3*t/(10 + np.sqrt(4 - 3*t)))

# Parameters
n_samples = 100_000
a = 2  # Semi-major axis
b = 5  # Semi-minor axis

# 1. Inverse CDF Sampling
# inverse_samples = inverse_cdf_sampling(n_samples)
# plt.hist(inverse_samples, bins=50, density=True, label='Inverse CDF Samples', alpha=0.7)
# plt.title('Sampling from $e^x$ using Inverse CDF (0 to 1)')
# plt.xlabel('x')
# plt.ylabel('Density')
# plt.legend()
# plt.show()
#
# # 2. Rejection Sampling
# rejection_samples = rejection_sampling_exp(n_samples)
# plt.hist(rejection_samples, bins=50, density=True, label='Rejection Sampling', alpha=0.7)
# plt.title('Sampling from $e^x$ using Rejection Sampling (0 to 1)')
# plt.xlabel('x')
# plt.ylabel('Density')
# plt.legend()
# plt.show()

# 3. Monte Carlo Integration for Ellipse Circumference
mean, std = arc_length_ellipse(a, b, n_samples)
print(f"Estimated Circumference of the Ellipse (a = {a}, b = {b}): mean:{mean}, std:{std}\n"
      f"Using Ramanujan aproximation of the circumference of the elipse: {approx(a,b)}")
