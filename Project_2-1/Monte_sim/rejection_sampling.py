import numpy as np


def rejection_sampling(pdf, x_min, x_max, n_samples, y_max):
    """
    Perform rejection sampling for a given probability distribution function (PDF).

    Parameters
    ----------
    pdf : function
        The probability distribution function to sample from.
    x_min : float
        Minimum value of the domain for sampling.
    x_max : float
        Maximum value of the domain for sampling.
    n_samples : int
        Number of samples to generate.
    y_max : float
        The maximum value of the PDF in the sampling range.

    Returns
    -------
    samples : np.ndarray
        Array of accepted samples.
    """
    samples = []

    while len(samples) < n_samples:
        # Step 1: Sample x from a uniform distribution over [x_min, x_max]
        x = np.random.uniform(x_min, x_max)

        # Step 2: Sample y from a uniform distribution over [0, y_max]
        y = np.random.uniform(0, y_max)

        # Step 3: Check if y is under the curve of the PDF
        if y <= pdf(x):
            samples.append(x)

    return np.array(samples)


# Example target PDF (normalized Gaussian)
def gaussian_pdf(x, mean=0, std_dev=1):
    """
    Gaussian probability density function (PDF).

    Parameters
    ----------
    x : float
        The point at which to evaluate the PDF.
    mean : float
        Mean of the Gaussian distribution.
    std_dev : float
        Standard deviation of the Gaussian distribution.

    Returns
    -------
    float
        The value of the PDF at x.
    """
    return (1 / (std_dev * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mean) / std_dev) ** 2)
