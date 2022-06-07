import torch
import math

import numpy as np
import pandas as pd

from hypothesis_tests import get_test, NEWCORR, HSIC

# Parameters.
N = 100
NOISE_LEVELS = 10
NOISE_SCALE = 3
ALPHA = 0.05
# Number of samples to use to estimate the power (number of tests).
NUM_SIMS_TEST = 500
# Number of samples to use to estimate the distribution of the test statistic.
NUM_SIMS_DIST = 500

# Relationships.
LINEAR = "linear"
STEP_FUNC = "step_function"
W_SHAPED = "w_shaped"
SINUSOID = "sinusoid"
CIRCULAR = "circular"
HETERO = "heteroskedastic"


def get_data(relationship, noise_level):
    if relationship == LINEAR:
        return get_linear(noise_level)
    elif relationship == STEP_FUNC:
        return get_step_function(noise_level)
    elif relationship == W_SHAPED:
        return get_w_shaped(noise_level)
    elif relationship == SINUSOID:
        return get_sinusoid(noise_level)
    elif relationship == CIRCULAR:
        return get_circular(noise_level)
    elif relationship == HETERO:
        return get_heteroskedastic(noise_level)
    else:
        raise ValueError(f"Unrecognized function type {relationship}!")


def get_linear(noise_level):
    x = 2 * torch.rand(N) - 1
    noise = torch.normal(torch.zeros(N), torch.ones(N))
    y = 0.5 * x + NOISE_SCALE * noise_level / NOISE_LEVELS * noise
    return x, y


# y = -3*(x< -.5) + 2*(x >= -.5 & x<0) - 4*(x>=0 & x < .5) - 3*(x>=.5) + 10*(l/num.noise)*rnorm(n)


def get_step_function(noise_level):
    x = 2 * torch.rand(N) - 1
    noise = torch.normal(torch.zeros(N), torch.ones(N))
    y = (
        -3.0 * (x < -0.5)
        + 2.0 * (torch.logical_and(x >= -0.5, x < 0))
        - 4.0 * (torch.logical_and(x >= 0, x < 0.5))
        - 3.0 * (x >= 0.5)
    )
    y += 10.0 * noise_level / NOISE_LEVELS * noise
    return x, y


def get_w_shaped(noise_level):
    x = 2 * torch.rand(N) - 1
    noise = torch.normal(torch.zeros(N), torch.ones(N))
    y = torch.abs(x + 0.5) * (x < 0) + torch.abs(x - 0.5) * (x >= 0)
    y += 0.25 * NOISE_SCALE * noise_level / NOISE_LEVELS * noise
    return x, y


def get_sinusoid(noise_level):
    x = 2 * torch.rand(N) - 1
    noise = torch.normal(torch.zeros(N), torch.ones(N))
    y = torch.cos(8 * math.pi * x) + NOISE_SCALE * noise_level / NOISE_LEVELS * noise
    return x, y


def get_circular(noise_level):
    x = 2 * torch.rand(N) - 1
    noise = torch.normal(torch.zeros(N), torch.ones(N))
    rademacher = 2 * torch.bernoulli(0.5 * torch.ones(N)) - 1
    y = rademacher * torch.sqrt(1 - x ** 2)
    y += 0.3 * NOISE_SCALE * noise_level / NOISE_LEVELS * noise
    return x, y


def get_heteroskedastic(noise_level):
    x = 2 * torch.rand(N) - 1
    noise = torch.normal(torch.zeros(N), torch.ones(N))
    y = NOISE_SCALE * (
        3 * (torch.abs(x) < 0.5) * (1 - noise_level / NOISE_LEVELS)
        + noise_level / NOISE_LEVELS
    )
    y *= noise
    return x, y


def save_distributions(relationship, marginal=False, num_sims=NUM_SIMS_DIST):
    for j, noise_level in enumerate(range(NOISE_LEVELS + 1)):
        data = torch.zeros(N, 2 * num_sims)
        for i in range(num_sims):
            x, y = get_data(relationship, noise_level)
            if marginal:
                # Resimulate to get the null scenario.
                x = 2 * torch.rand(N) - 1
            data[:, 2 * i] = x
            data[:, 2 * i + 1] = y
        df = pd.DataFrame(data.numpy())
        if marginal:
            df.to_csv(
                f"data/power/{relationship}_noise_level_{j}_marginal.csv", index=False
            )
        else:
            df.to_csv(
                f"data/power/{relationship}_noise_level_{j}_joint.csv", index=False
            )


def compute_critical_vals(
    relationship=LINEAR, test_name=NEWCORR, num_sims=NUM_SIMS_DIST
):
    test = get_test(test_name)
    c_vals = torch.zeros(NOISE_LEVELS + 1)
    for noise_level in range(NOISE_LEVELS + 1):
        test_stats = torch.zeros(num_sims)
        for i in range(num_sims):
            x, y = get_data(relationship, noise_level)
            # Resimulate to get the null scenario.
            x = 2 * torch.rand(N) - 1
            test_stats[i] = test(x, y)
        c_vals[noise_level] = torch.quantile(test_stats, 1 - ALPHA)

    torch.save(
        c_vals, f"results/power/{relationship}_{test_name}_critical_values.pt",
    )


def compute_powers(relationship=LINEAR, test_name=NEWCORR, num_sims=NUM_SIMS_TEST):
    test = get_test(test_name)
    powers = torch.zeros(NOISE_LEVELS + 1)
    c_vals = torch.load(f"results/power/{relationship}_{test_name}_critical_values.pt")
    for i, noise_level in enumerate(range(NOISE_LEVELS + 1)):
        rejects = 0
        for _ in range(num_sims):
            x, y = get_data(relationship, noise_level)
            if test(x, y, compute_pvalue=False) >= c_vals[i]:
                rejects += 1
        powers[i] = rejects / num_sims
    torch.save(powers, f"results/power/{relationship}_{test_name}_powers.pt")


def compute_true_powers(relationship=LINEAR, test_name=NEWCORR, num_sims=NUM_SIMS_TEST):
    test = get_test(test_name)
    powers = torch.zeros(NOISE_LEVELS + 1)
    for i, noise_level in enumerate(range(NOISE_LEVELS + 1)):
        rejects = 0
        for _ in range(num_sims):
            x, y = get_data(relationship, noise_level)
            stat, pval = test(x, y, compute_pvalue=True)
            if pval <= ALPHA:
                rejects += 1
        powers[i] = rejects / num_sims
    torch.save(powers, f"results/power/{relationship}_{test_name}_true_powers.pt")


if __name__ == "__main__":
    for relationship in [STEP_FUNC, W_SHAPED, SINUSOID, CIRCULAR, HETERO]:
        # for relationship in [LINEAR]:
        compute_true_powers(relationship=relationship, test_name=HSIC)

