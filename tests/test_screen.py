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
            np.random.rand(3, 1),
            columns=[
                "a",
            ],
            index=["g1", "g2", "g3"],
        ),
        obs=pd.DataFrame(index=["c1", "c2", "c3", "c4"]),
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

    