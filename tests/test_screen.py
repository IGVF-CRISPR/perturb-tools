import anndata as ad
import numpy as np
import pandas as pd
import perturb_tools as pt
import pytest


@pytest.fixture
def screen():
    return pt.Screen(
        X=np.random.rand(3, 4),
        guides=pd.DataFrame(
            np.random.rand(3, 1),
            columns=[
                "a",
            ],
            index=["g1", "g2", "g3"],
        ),
        condit=pd.DataFrame(index=["c1", "c2", "c3", "c4"]),
    )


@pytest.fixture
def adata():
    return ad.AnnData(
        X=np.random.rand(3, 4),
        obs=pd.DataFrame(
            np.random.rand(3, 1),
            columns=[
                "a",
            ],
            index=["g1", "g2", "g3"],
        ),
        var=pd.DataFrame(index=["c1", "c2", "c3", "c4"]),
    )


class TestClass:
    def test_read(self, screen):
        screen.write("tmp")
        pt.read_h5ad("tmp")

    def test_subset(self, screen):
        screen_sub = screen[["g1", "g2"], :]
        assert screen_sub.shape == (2, 4)
        assert len(screen_sub.obs) == 2

    def test_select(self, screen):
        screen_sel = screen[:, ["c1", "c3"]]
        assert screen_sel.shape == (3, 2)
        assert len(screen_sel.var) == 2

    def test_from_adata(self, adata, screen):
        screen = pt.Screen.from_adata(adata)
        assert (screen.guides.index == adata.obs.index).all()
        assert (screen.condit.index == adata.var.index).all()
