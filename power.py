import numpy as np

from hypothesis_tests import get_test, NEWCORR

N = 100
NOISE_LEVELS = 10
NOISE_SCALE = 3
NUM_SIMS_TEST = 500  # Number of samples to use to estimate the power (number of tests).
NUM_SIMS_DIST = (
    500  # Number of samples to use to estimate the distribution of the test statistic.
)
LINEAR = "linear"
ALPHA = 0.05


def get_data(relationship, noise_level):
    if relationship == LINEAR:
        return get_linear(noise_level)
    else:
        raise ValueError(f"Unrecognized function type {relationship}!")


def get_linear(noise_level):
    x = np.random.uniform(-1, 1, N)
    y = 0.5 * x + NOISE_SCALE * noise_level / NOISE_LEVELS * np.random.normal(size=N)
    return x, y


# TODO: parallelize.
def compute_critical_vals(
    relationship=LINEAR, test_name=NEWCORR, num_sims=NUM_SIMS_DIST
):
    test = get_test(test_name)
    c_vals = []
    for noise_level in range(NOISE_LEVELS + 1):
        test_stats = []
        for _ in range(num_sims):
            x, y = get_data(relationship, noise_level)
            test_stats.append(test(x, y, compute_pvalue=False))
        test_stats = np.array(test_stats)
        c_vals.append(np.quantile(test_stats, 1 - ALPHA))

    np.save(
        f"results/power/{relationship}_{test_name}_critical_values.npy",
        np.array(c_vals),
    )


# TODO: parallelize.
def compute_powers(relationship=LINEAR, test_name=NEWCORR, num_sims=NUM_SIMS_TEST):
    test = get_test(test_name)
    powers = []
    c_vals = np.load(f"results/power/{relationship}_{test_name}_critical_values.npy",)
    for i, noise_level in enumerate(range(NOISE_LEVELS + 1)):
        rejects = 0
        for _ in range(num_sims):
            x, y = get_data(relationship, noise_level)
            if test(x, y, compute_pvalue=False) >= c_vals[i]:
                rejects += 1
        powers.append(rejects / num_sims)
    np.save(f"results/power/{relationship}_{test_name}_powers.npy", np.array(powers))


# compute_dists(num_sims=500)
compute_powers(num_sims=500)

