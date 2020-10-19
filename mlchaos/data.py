# AUTOGENERATED! DO NOT EDIT! File to edit: nbs/data.ipynb (unless otherwise specified).

__all__ = ['load_poincare_maps', 'load_index_file', 'load_poincare_index_pair', 'TSDataChaos', 'get_motion_items',
           'show_labelled_ic_map']

# Cell
from fastcore.all import *
from timeseries.all import *
import pandas as pd
import numpy as np
from .utils import df_slicer
import seaborn as sns

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
def load_index_file(fname, index_col=7, uncertainty_index=-1.):
    "Returns the index of an index file. In case the argument `index_col` has more \
    than one value, the value in `uncertainty_index` will be set."
    indices = pd.read_table(fname,
                         sep='\s+',
                         header=None,
                         usecols=L(index_col),
                         squeeze=True).values
    # indices is a numpy array
    if len(indices.shape) > 1:
        indices = array([x[0] if np.all(x == x[0]) else uncertainty_index for x in indices])
    return indices

# Cell
@delegates(to=load_index_file, but=['fname'])
def load_poincare_index_pair(fname_poincare, fname_index, **kwargs):
    "Load the x data from a Poincare file and the y data from the index file.\
    Returns a tuple of 2 numpy arrays: "
    "x : array with a shape (n_samples, n_channels, sequence_length)"
    "y : array with a shape (n_samples)"
    "for the Poincare maps, n_channels is 2 (x and y)"
    return load_poincare_maps(fname_poincare), load_index_file(fname_index, **kwargs)

# Cell
class TSDataChaos(TSData):
    @classmethod
    @delegates(to=load_poincare_index_pair)
    def from_poincare_and_index_files(cls, fnames, **kwargs):
        "`fnames` is a list of pairs (poincare_file, index_file), or a single pair."
        self = cls(fnames)
        self.x = []
        self.y = []
        self.dsname = []
        self.fnames = []
        self.ds = []
        xs,ys, dss = [],[],[]
        if isinstance(fnames, list):
            for i, (fn_poincare, fn_index) in enumerate(fnames):
                dsname = fn_poincare.parent.name
                x, y = load_poincare_index_pair(fn_poincare, fn_index, **kwargs)
                xs.append(x)
                ys.append(y)
                dss.append(np.repeat(i, len(y))) # Map each item with the dataset it belongs to
                self.fnames.append((fn_poincare, fn_index))
                self.dsname.append(dsname)
            self.x = np.concatenate(xs)
            self.y = np.concatenate(ys)
            self.ds = np.concatenate(dss).squeeze()
        else:
            fn_poincare, fn_index = fnames
            dsname = fn_poincare.parent.name
            self.fnames.append(fnames)
            self.dsname.append(dsname)
            self.x, self.y = load_poincare_index_pair(fn_poincare, fn_index, **kwargs)
            self.ds = np.repeat(0, len(y))
        return self

    @classmethod
    def from_poincare_maps(cls, fnames):
        "`fnames` is a list of paths to Poincare maps. No ys are provided"
        self = cls(fnames)
        self.x = []
        self.y = None
        self.dsname = []
        self.fnames = []
        xs = []
        if isinstance(fnames, list):
            for fn_poincare in fnames:
                x = load_poincare_maps(fn_poincare)
                xs.append(x)
                self.fnames.append((fn_poincare, fn_index))
                self.dsname.append(fn_poincare.parent.name)
            self.x = np.concatenate(xs)
        else:
            fn_poincare = fnames
            self.fnames.append(fn_poincare)
            self.dsname.append(fn_poincare.parent.name)
            self.x = load_poincare_maps(fn_poincare)
        return self

# Cell
def get_motion_items(fnames):
    "get_ts_items return list of tuples. Each tuple corresponds \
    to a timeserie (nump.ndarray) and a label (string). fnames is not \
    an array of paths, but an array of tuples of paths. Each tuple has the \
    path of the Poincare map and the path of the index file"
    data = TSDataChaos.from_poincare_and_index_files(fnames)
    return data.get_items()

# Cell
@delegates(sns.scatterplot, but=['data', 'x', 'y', 'hue', 'marker'])
def show_labelled_ic_map(poinc_maps:Union[torch.Tensor, np.ndarray]=None,
                         lbls:list=None, **kwargs):
    "Show a scatter plot with the initial conditions (x0, y0) of each Poincare map in \
    `poinc_maps`, coloured with the labels given in `lbls`.'. The argument `poinc_maps` \
    expects a (n_items, n_channels, n_points) numpy array or tensor"
    initial_conditions = poinc_maps[:,:,0]
    qux = pd.DataFrame(initial_conditions, columns=['x0', 'y0'])
    qux['class'] = lbls
    return sns.scatterplot(data=qux, x='x0', y='y0', hue='class', marker='.', **kwargs)