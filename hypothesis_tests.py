import torch
import math
from scipy.stats import norm

NEWCORR = "new_correlation"
COR = "correlation"
MAXCOR = "maxcor"
DCOR = "dcor"
HSIC = "hsic"
MIC = "mic"
TIC = "tic"
HHG = "hhg"


def get_test(test_name):
    if test_name == NEWCORR:
        return chatterjee
    else:
        raise ValueError(f"Unrecognized test {test_name}!")


def correlation(x, y, compute_pvalue=False):
    x = x - x.mean()
    y = y - y.mean()

    return (x @ y) / math.sqrt((x ** 2).sum() * (y ** 2).sum())


def chatterjee(x, y, compute_pvalue=False):
    n = len(x)

    # Ties are broken at random. For the "argsort" command,
    # they are decided based on which index comes first, so we shuffle.
    perm = torch.randperm(n)
    x, y = x[perm], y[perm]

    # Compute rankings.
    idx = torch.argsort(x)
    y_ord = y[idx]
    # TODO: Optimize ranking.
    r = torch.tensor([(y <= y_ord[i]).sum() for i in range(n)])
    l = torch.tensor([(y >= y_ord[i]).sum() for i in range(n)])

    # Compute statistic.
    num = n * torch.abs(r[1:] - r[:-1]).sum()
    denom = 2 * (l * (n - l)).sum()
    stat = (1 - num / denom).item()

    # Compute p-value.
    if compute_pvalue:
        # TODO: optimize ranking with argsort(argsort).
        R = torch.tensor([(y <= y[i]).sum() for i in range(n)])
        L = torch.tensor([(y >= y[i]).sum() for i in range(n)])
        u = torch.sort(R)[0]
        v = torch.cumsum(u, 0)

        ind = torch.arange(n) + 1
        a = ((2 * n - 2 * ind + 1) * (u ** 2)).sum() / n ** 4
        b = ((v + (n - ind) * u) ** 2).sum() / n ** 5
        c = ((2 * n - 2 * ind + 1) * u).sum() / n ** 3
        d = (L * (n - L)).sum() / n ** 3

        tau2 = ((a - 2 * b + c ** 2) / d ** 2).item()
        pvalue = 1 - norm.cdf(stat * math.sqrt(n / tau2))
        return stat, pvalue

    return stat


def benjamini_hochberg(pvalues, fdr=0.05):
    m = len(pvalues)
    sorted_vals, idx = torch.sort(pvalues)
    thres = fdr * (torch.arange(m) + 1) / m
    candidates = torch.arange(m)[sorted_vals < thres]
    if len(candidates) > 0:
        num_rejects = torch.max(candidates).item()
        return idx[0:num_rejects]
    else:
        return []
