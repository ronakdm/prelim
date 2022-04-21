import torch

# import numpy as np

from hypothesis_tests import get_test, NEWCORR

N = 100
NOISE_LEVELS = 10
NOISE_SCALE = 3
LINEAR = "linear"
ALPHA = 0.05

# Number of samples to use to estimate the power (number of tests).
NUM_SIMS_TEST = 500
# Number of samples to use to estimate the distribution of the test statistic.
NUM_SIMS_DIST = 500


def get_data(relationship, noise_level):
    if relationship == LINEAR:
        return get_linear(noise_level)
    else:
        raise ValueError(f"Unrecognized function type {relationship}!")


def get_linear(noise_level):
    x = 2 * torch.rand(N) - 1
    noise = torch.normal(torch.zeros(N), torch.ones(N))
    y = 0.5 * x + NOISE_SCALE * noise_level / NOISE_LEVELS * noise
    return x, y


# TODO: parallelize.
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


# TODO: parallelize.
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


if __name__ == "main":
    compute_critical_vals(num_sims=500)
    compute_powers(num_sims=500)

