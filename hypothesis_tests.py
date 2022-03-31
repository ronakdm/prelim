import numpy as np
from scipy.stats import norm


def correlation(x, y):
    x = x - x.mean()
    y = y - y.mean()

    return (x @ y) / np.sqrt((x ** 2).sum() * (y ** 2).sum())


def new_correlation(x, y, compute_pvalue=True):
    n = len(x)

    # Ties are broken at random. For the "argsort" command,
    # they are decided based on which index comes first, so we shuffle.
    perm = np.random.permutation(n)
    x, y = x[perm], y[perm]

    # Compute rankings.
    idx = np.argsort(x)
    y_ord = y[idx]
    r = np.array([(y <= y_ord[i]).sum() for i in range(n)])
    l = np.array([(y >= y_ord[i]).sum() for i in range(n)])

    # Compute statistic.
    num = n * np.array([np.abs(r[i + 1] - r[i]) for i in range(n - 1)]).sum()
    denom = 2 * (l * (n - l)).sum()
    stat = 1 - num / denom

    # Compute p-value.
    pvalue = None
    if compute_pvalue:
        R = np.array([(y <= y[i]).sum() for i in range(n)])
        L = np.array([(y >= y[i]).sum() for i in range(n)])
        u = np.sort(R)
        v = np.array([u[0:i].sum() for i in range(n)])

        ind = np.arange(n)
        a = ((2 * n - 2 * ind + 1) * (u ** 2)).sum() / n ** 4
        b = ((v + (n - ind) * u) ** 2).sum() / n ** 5
        c = ((2 * n - 2 * ind + 1) * u).sum() / n ** 3
        d = (L * (n - L)).sum() / n ** 3

        tau2 = (a - 2 * b + c ** 2) / d ** 2
        pvalue = 1 - norm.cdf(stat * np.sqrt(n / tau2))

    return stat, pvalue
