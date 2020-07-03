# AUTOGENERATED! DO NOT EDIT! File to edit: nbs/data.ipynb (unless otherwise specified).

__all__ = ['load_poincare_maps', 'load_index_file', 'load_poincare_index_pair', 'TSDataChaos']

# Cell
from fastcore.all import *
from timeseries.all import *
import pandas as pd
import numpy as np
from .utils import df_slicer

# Cell
def load_poincare_maps(fname):
    "Load the data from a Poincare map.\
    Returns a numpy array with a shape (n_orbits, 2, sequence_length). \
    The time column is removed from the data, since it is always multiples of 2*pi"
    df = pd.read_table(fname, sep='\s+', names=['time', 'x', 'y'])
    seq_length = set(df.time).__len__()
    df = df.drop('time', axis=1)
    nparr = df_slicer(df, w=seq_length, s=seq_length)
    return nparr.transpose([0, 2, 1])

# Cell
def load_index_file(fname):
    return pd.read_table(fname,
                         sep='\s+',
                         header=None,
                         usecols=[7],
                         squeeze=True).values

# Cell
def load_poincare_index_pair(fname_poincare, fname_index):
    "Load the x data from a Poincare file and the y data from the index file.\
    Returns a tuple of 2 numpy arrays: "
    "x : array with a shape (n_samples, n_channels, sequence_length)"
    "y : array with a shape (n_samples)"
    "for the Poincare maps, n_channels is 2 (x and y)"
    return load_poincare_maps(fname_poincare), load_index_file(fname_index)

# Cell
class TSDataChaos(TSData):
    @classmethod
    def from_poincare_and_index_files(cls, fnames):
        "`fnames` is a list of pairs (poincare_file, index_file), or a single pair."
        self = cls(fnames)
        self.x = []
        self.y = []
        self.dsname = []
        self.fnames = []
        xs,ys = [],[]
        if isinstance(fnames, list):
            for fn_poincare, fn_index in fnames:
                x, y = load_poincare_index_pair(fn_poincare, fn_index)
                xs.append(x)
                ys.append(y)
                self.fnames.append((fn_poincare, fn_index))
                data.dsname.append((fn_poincare.stem, fn_index.stem))
            self.x = np.concatenate(xs)
            self.y = np.concatenate(ys)
        else:
            fn_poincare, fn_index = fnames
            self.fnames.append(fnames)
            self.dsname.append(fn_poincare.parent.stem)
            self.x, self.y = load_poincare_index_pair(fn_poincare, fn_index)
        return self