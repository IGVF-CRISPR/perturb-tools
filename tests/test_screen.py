import anndata as ad
import numpy as np
import pandas as pd
import perturb_tools as pt
import pytest


@pytest.fixture
def screen():
    return ad.AnnData(
        X=np.random.randint(low=0, high=1000, size=(4, 3)),
        dtype=int,
        var=pd.DataFrame(
            {"target": ["t1", "t1", "t2"]},
            index=["g1", "g2", "g3"],
        ),
        obs=pd.DataFrame(
            {"replicate": [1, 2, 1, 2], "condition": [0, 0, 1, 1]},
            index=["c1", "c2", "c3", "c4"],
        ),
    )


class TestClass:
    def test_subset(self, screen):
        screen_sub = screen[:, ["g1", "g2"]]
        assert screen_sub.shape == (4, 2)
        assert len(screen_sub.var) == 2

    def test_select(self, screen):
        screen_sel = screen[["c1", "c3"], :]
        assert screen_sel.shape == (2, 3)
        assert len(screen_sel.obs) == 2

    def test_lognorm(self, screen):
        pt.pp.log_norm(screen)

    def test_lfcs(self, screen):
        pt.pp.log_norm(screen)
        pt.pp.log_fold_change(screen, "c1", "c3")

    def test_lfc_reps(self, screen):
        pt.pp.log_norm(screen)
        pt.pp.log_fold_change_reps(screen, 1, 0, compare_col="condition")
        pt.pp.log_fold_change_aggregate(screen, 1, 0, compare_col="condition")

    def test_write(self, screen):
        pt.pp.log_norm(screen)
        pt.io.to_Excel(screen, "test.xlsx")
        pt.io.to_mageck_input(screen, "test.txt", target_column="target")
