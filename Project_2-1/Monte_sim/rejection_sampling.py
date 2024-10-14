import numpy as np

def rejection_sampling(pdf, x_min, x_max, n_samples):
    """
    Perform rejection sampling for a given probability distribution function.

    Parameters
    ----------
    pdf : function
        The probability distribution function to sample from.
    x_min : float
        Minimum value of the domain.
    x_max : float
        Maximum value of the domain.
    n_samples : int
        Number of samples to generate.

    Returns
    -------
    np.ndarray
        Samples generated using rejection sampling.
    """
    samples = []
    while len(samples) < n_samples:
        x = np.random.uniform(x_min, x_max)
        y = np.random.uniform(0, pdf(x_max))
        if y < pdf(x):
            samples.append(x)
    return np.array(samples)
