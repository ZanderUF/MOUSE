# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED

import numpy as np

def create_lognormal_sampler(low_cost, high_cost, class3_cost):
    """Draw one sample from a lognormal fitted to three cost estimates.

    Fits mu = mean(ln(low), ln(class3), ln(high)) and sigma = population std (ddof=0)
    of those logs, then returns a single LogNormal(mu, sigma) draw in $.
    """
    # Calculate the natural logarithms of the given costs
    ln_low_cost = np.log(low_cost)
    ln_high_cost = np.log(high_cost)
    ln_class3_cost = np.log(class3_cost)

    # Calculate the mean (mu) and standard deviation (sigma) of the logarithms
    mu = np.mean([ln_low_cost, ln_high_cost, ln_class3_cost])
    sigma = np.std([ln_low_cost, ln_high_cost, ln_class3_cost], ddof=0) # Population std deviation

    # Define the sampler function
    def sampler():
        return np.random.lognormal(mean=mu, sigma=sigma)

    return sampler()

def truncated_normal_sample(mean, std, lower_bound, upper_bound):
    """Draw one sample from N(mean, std) truncated to [lower_bound, upper_bound].

    Uses rejection sampling: repeatedly draws from the normal until a value falls
    within the bounds. Returns one draw (dimensionless, e.g. a scaling exponent).
    """
    while True:
        sample = np.random.normal(mean, std)
        if lower_bound <= sample <= upper_bound:
            return sample

def uniform_sample(low, high):
    """Draw one sample from a uniform distribution U(low, high).

    Returns a single draw; units match those of low/high (e.g. $).
    """
    return np.random.uniform(low, high)

def sampler(distribution, **kwargs):
    """Dispatch to a distribution sampler by name and return one draw.

    Routes 'Lognormal' to create_lognormal_sampler(low_cost, high_cost, class3_cost),
    'Truncated Normal' to truncated_normal_sample(mean, std, lower_bound, upper_bound),
    and 'Uniform' to uniform_sample(low, high) via kwargs; raises ValueError otherwise.
    Returns one draw (units match the chosen distribution, e.g. $).
    """
    if distribution == "Lognormal":
        return create_lognormal_sampler(kwargs['low_cost'], kwargs['high_cost'], kwargs['class3_cost'])
    elif distribution == "Truncated Normal":
        return truncated_normal_sample(kwargs['mean'], kwargs['std'], kwargs['lower_bound'], kwargs['upper_bound'])
    elif distribution == "Uniform":
        return uniform_sample(kwargs['low'], kwargs['high'])    
    else:
        raise ValueError("Unavailable Distribution")
