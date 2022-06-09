import torch
import math
from scipy.stats import norm

from hyppo.independence import Hsic

NEWCORR = "new_correlation"
BOOSTED = "boosted_chatterjee"
COR = "correlation"
# MAXCOR = "maxcor"
DCOR = "dcor"
HSIC = "hsic"
MIC = "mic"
TIC = "tic"
HHG = "hhg"


def get_test(test_name):
    if test_name == NEWCORR:
        return chatterjee
    elif test_name == HSIC:
        return hsic
    elif test_name == BOOSTED:
        return boosted_chatterjee
    else:
        raise ValueError(f"Unrecognized test {test_name}!")


def correlation(x, y, compute_pvalue=False):
    x = x - x.mean()
    y = y - y.mean()

    return (x @ y) / math.sqrt((x ** 2).sum() * (y ** 2).sum())


def hsic(x, y, compute_pvalue=False):
    hsic = Hsic()
    x = x.numpy().reshape(-1, 1)
    y = y.numpy().reshape(-1, 1)
    if compute_pvalue:
        stat, pvalue = hsic.test(x, y)
        return stat, pvalue
    else:
        return hsic.statistic(x, y)


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
        a = ((2 * n - 2 * ind + 1) * (u ** 2)).sum()
        b = ((v + (n - ind) * u) ** 2).sum() / n
        c = ((2 * n - 2 * ind + 1) * u).sum() / n
        d = (L * (n - L)).sum() / n

        tau2 = ((a - 2 * b + c ** 2) / d ** 2).item()
        pvalue = 1 - norm.cdf(stat * math.sqrt(n / tau2))
        if math.isnan(pvalue):
            pvalue = 1e-12
        return stat, pvalue

    return stat


def boosted_chatterjee_stat(rank, n, M):
    denom = (n + 1) * (n * M + M * (M + 1) / 4)
    num = 0
    for i in range(n):
        for j in range(1, M + 1):
            if i + j >= n:
                num += rank[i]
            else:
                num += min(rank[i], rank[i + j])

    return -2 + 6 * num / denom


def boosted_chatterjee(x, y, compute_pvalue=False):
    n = len(x)
    M = int(math.sqrt(n))

    # Ties are broken at random. For the "argsort" command,
    # they are decided based on which index comes first, so we shuffle.
    perm = torch.randperm(n)
    x, y = x[perm], y[perm]

    # Compute rankings.
    idx = torch.argsort(x)
    y_ord = y[idx]
    r = torch.tensor([(y <= y_ord[i]).sum() for i in range(n)])
    r2 = get_ranks(y)[idx]  # TODO: Check that these are the same.

    stat1 = boosted_chatterjee_stat(r, n, M)
    stat2 = boosted_chatterjee_stat(n - r + 1, n, M)
    stat = max(stat1, stat2)

    if compute_pvalue:
        n_sim = 1000
        idx = torch.cat([torch.randperm(n).unsqueeze(0) for i in range(n_sim)])
        R = r[idx]
        # null_dist = torch.zeros(n_sim)
        # for b in range(n_sim):
        #     # r = r[torch.randperm(n)]
        #     r = R[b]
        #     # R[b, :] = r.clone()
        #     stat1 = boosted_chatterjee_stat(r, n, M)
        #     stat2 = boosted_chatterjee_stat(n - r + 1, n, M)
        #     null_dist[b] = max(stat1, stat2)

        denom = (n + 1) * (n * M + M * (M + 1) / 4)
        stats1 = torch.zeros(n_sim)
        stats2 = torch.zeros(n_sim)
        for j in range(1, M + 1):
            U = torch.stack([R[:, 0 : n - j], R[:, j:n]])
            V = torch.min(U, dim=0)[0]
            stats1 += torch.sum(V, dim=1)
            stats1 += torch.sum(R[:, n - j : n], axis=1)

            L = n - R + 1
            U = torch.stack([L[:, 0 : n - j], L[:, j:n]])
            V = torch.min(U, dim=0)[0]
            stats2 += torch.sum(V, dim=1)
            stats2 += torch.sum(L[:, n - j : n], axis=1)

        stats1 = -2 + 6 * stats1 / denom
        stats2 = -2 + 6 * stats2 / denom

        null_dist = torch.max(torch.stack([stats1, stats2]), dim=0)[0]
        pvalue = torch.sum((null_dist >= stat).int()) / n_sim

        return stat, pvalue

    return stat


def benjamini_hochberg(pvalues, fdr=0.05):
    m = len(pvalues)
    sorted_vals, idx = torch.sort(pvalues)
    thres = fdr * (torch.arange(m) + 1) / m
    candidates = torch.arange(m)[sorted_vals <= thres]
    if len(candidates) > 0:
        num_rejects = torch.max(candidates).item()
        return torch.sort(idx[0 : (num_rejects + 1)].int())[0]
    else:
        return []


def get_ranks(z):
    uq_vals, inverse_ixs, counts = torch.unique(
        z, sorted=True, return_inverse=True, return_counts=True
    )

    cumsum_counts = counts.cumsum(dim=0)
    possible_ranks = torch.zeros_like(counts)
    possible_ranks[1:] = cumsum_counts[:-1]
    ranks = possible_ranks[inverse_ixs]

    return cumsum_counts[inverse_ixs]
