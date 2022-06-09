import torch

from hyppo.tools import SIMULATIONS
from hypothesis_tests import get_test, NEWCORR, BOOSTED

sim_name = "linear"
sim = SIMULATIONS[sim_name]
n = 100

test_func = get_test(BOOSTED)
x, y = sim(n, 1, noise=True)
x, y = torch.tensor(x).reshape(-1), torch.tensor(y).reshape(-1)
stat, pvalue = test_func(x, y, compute_pvalue=True)

print("lol")
