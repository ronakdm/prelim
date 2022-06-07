import torch
import numpy as np

from hyppo.tools import SIMULATIONS
from hypothesis_tests import get_test, NEWCORR, HSIC


def compute_power(sim_name, n, test_name, n_sims=100):
    test_func = get_test(test_name)
    sim = SIMULATIONS[sim_name]
    rejects = 0
    for _ in range(n_sims):
        x, y = sim(n, 1, noise=True)
        x, y = torch.tensor(x).reshape(-1), torch.tensor(y).reshape(-1)
        perm = torch.randperm(n)
        y = y[perm]
        stat, pvalue = test_func(x, y, compute_pvalue=True)
        if pvalue <= 0.05:
            rejects += 1
    return rejects / n_sims


for sim_name in list(SIMULATIONS.keys())[0:-2]:
    n = 30
    test_name = NEWCORR
    n_sims = 10
    print(compute_power(sim_name, n, test_name, n_sims=n_sims))

