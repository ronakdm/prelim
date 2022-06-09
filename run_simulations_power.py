from hyppo.tools import SIMULATIONS
from hypothesis_tests import get_test, NEWCORR, BOOSTED

from joblib import Parallel, delayed
import pickle

import torch
import numpy as np


def compute_power(sim_name, n, test_name, n_sims=100):
    test_func = get_test(test_name)
    sim = SIMULATIONS[sim_name]

    def worker():
        x, y = sim(n, 1, noise=True)
        x, y = torch.tensor(x).reshape(-1), torch.tensor(y).reshape(-1)
        # perm = torch.randperm(n)
        # y = y[perm]
        stat, pvalue = test_func(x, y, compute_pvalue=True)
        if pvalue <= 0.05:
            return 1
        return 0

    rejects = np.array(Parallel(n_jobs=-3)(delayed(worker)() for _ in range(n_sims)))
    return rejects.mean()


sample_sizes = np.logspace(2, 4, 30).astype(np.int64)
pickle.dump(sample_sizes, open("results/boosting/sample_sizes.p", "wb"))

for n in sample_sizes:
    for sim_name in list(SIMULATIONS.keys())[0:-2]:
        print(f"n = {n} sim = {sim_name}")
        test_name = BOOSTED
        n_sims = 500
        power = compute_power(sim_name, int(n), test_name, n_sims=n_sims)
        pickle.dump(
            power,
            open(
                f"results/boosting/sim_{sim_name}_test_{test_name}_n_{n}_power.p", "wb"
            ),
        )
