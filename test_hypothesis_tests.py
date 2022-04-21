import unittest
import pandas as pd
import torch

from hypothesis_tests import chatterjee, benjamini_hochberg


class TestChatterjee(unittest.TestCase):
    def test_galton_pea_statistic(self):
        """
        Matches statistic reported in Galton's peas example (Section 3).
        """
        df = pd.read_csv("data/galton_peas_data.csv", header=0)
        x = torch.tensor(df["parent"].values)
        y = torch.tensor(df["child"].values)

        n_sims = 1000

        torch.manual_seed(0)

        stat_xy = torch.tensor([chatterjee(x, y) for _ in range(n_sims)]).mean().item()
        stat_yx = torch.tensor([chatterjee(y, x) for _ in range(n_sims)]).mean().item()

        msg = "xi_n(X, Y) do not agree on Galton's peas dataset."
        self.assertAlmostEqual(stat_xy, 0.11, 2, msg)

        msg = "xi_n(Y, X) do not agree on Galton's peas dataset."
        self.assertAlmostEqual(stat_yx, 0.92, 2, msg)

    def test_bernoulli_statistic(self):
        """
        Converges to population value for Bernoulli example (Section 4.2).
        """
        p = 0.4
        p_ = 0.5

        pop_value = p_ * (1 - p) / (1 - p * p_)
        msg = "Population value for dependent Bernoulli does not match."
        self.assertAlmostEqual(pop_value, 0.375, 5, msg)

        torch.manual_seed(0)

        n = 100000
        x = torch.bernoulli(p * torch.ones(n))
        z = torch.bernoulli(p_ * torch.ones(n))
        y = x * z

        stat = chatterjee(x, y)

        msg = "Statistic for dependent Bernoulli does not match for high n."
        # This is a very rough test, longer one takes much more time.
        self.assertAlmostEqual(stat, pop_value, 2, msg)

    def test_gene_asympt_pvalue(self):
        """
        Selects and does not select same genes as paper.
        """
        spellman = pd.read_csv("data/spellman_gene_expr_data.csv", header=0)
        x = torch.tensor(spellman["time"].values)
        genes = list(spellman.columns[1:])
        pvals = torch.zeros(len(genes))
        for i, gene in enumerate(genes):
            y = torch.tensor(spellman[gene].values)
            pvals[i] = chatterjee(x, y, compute_pvalue=True)[1]
        rejects_idx = benjamini_hochberg(pvals)
        lookup = {}
        for i in rejects_idx:
            lookup[genes[i]] = True

        picked = ["YOR308C", "YJL115W", "YGR177C", "YPR119W", "YKL127W", "YHR143W"]
        for gene in picked:
            msg = gene + " selected by paper, but not by my implementation."
            self.assertTrue(gene in lookup, msg)

        not_picked = [
            "YOR140W",
            "YCL021W",
            "YJR086W",
            "YLR406C",
            "YLR283W",
            "YOR378W",
        ]
        for gene in not_picked:
            msg = gene + " not selected by paper, but selected by my implementation."
            self.assertFalse(gene in lookup, msg)


if __name__ == "__main__":
    unittest.main()
