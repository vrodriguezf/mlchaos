# AUTOGENERATED! DO NOT EDIT! File to edit: nbs/utils.ipynb (unless otherwise specified).

__all__ = ['df_slicer', 'show_attribution_maps_at', 'ClassificationInterpretationAugmented',
           'log_last_plt_as_wandb_img', 'log_plotnine_as_wandb_img']

# Cell
from .imports import *
from fastai.basics import *
import numpy as np
from fastcore.all import *
from fastai.interpret import *
from fastai.data.all import *
import wandb
import matplotlib.pyplot as plt

# Cell
def df_slicer(df, w, s=1, padding=False, padding_value=0, return_as='ndarray'):
    "Transform a numeric dataframe `df` into slices (np arrays) of `w` \
    rows and the same number of columns than the original dataframe. The \
    distance between each slice is given by the stride `s`. If `padding` is \
    equals to True, the last slices which have less than `w` points are filled \
    with the value marked in the argument `padding_value`. Otherwise, those \
    slices are removed from the result."
    aux = [df.iloc[x:x+w] for x in range(0, len(df), s)]
    if padding:
        with_padding = [x.append(pd.DataFrame(
            np.full((w - len(x), len(df.columns)), padding_value),
            columns=df.columns.values)) if len(x) < w else x for x in aux]
    else:
        with_padding = [x for x in aux if len(x) == w]
    return np.rollaxis(np.dstack([x.values for x in with_padding]), -1)

# Cell
@patch
def __getitem__(self:Interpretation, idxs):
    "Get the inputs, preds, targets, decoded outputs, and losses at `idx`"
    if not is_listy(idxs): idxs = [idxs]
    attrs = 'inputs,preds,targs,decoded,losses'
    res = L([getattr(self, attr)[idxs] for attr in attrs.split(',')])
    return res

@patch
def show_at(self:Interpretation, idxs, **kwargs):
    inp, _, targ, dec, _ = self[idxs]
    self.dl.show_results((inp, dec), targ, **kwargs)

# Cell
def show_attribution_maps_at(idxs,  **kwargs):
    return None

# Cell
class ClassificationInterpretationAugmented(ClassificationInterpretation):
    def top_losses(self, k=None, largest=True, predicted=None, actual=None):
        r"""
        `k` largest(/smallest) losses and indexes, defaulting to all losses (sorted by `largest`)."
        `predicted` and `actual` arguments are passed as decoded labels
        """
        if predicted is None and actual is None:
            # Default behaviour
            return self.losses.topk(ifnone(k, len(self.losses)), largest=largest)
        else:
            # Subset losses by the conditions given in predicted and actual arguments
            cond_preds = (self.decoded == self.vocab.o2i[predicted]) if predicted else tensor(True)
            cond_actuals = (self.targs == self.vocab.o2i[actual]) if actual else tensor(True)
            idxs = (cond_preds & cond_actuals).nonzero().squeeze()
            loss_subset = self.losses[idxs].topk(ifnone(k, len(idxs)), largest=largest)
            # The indices in loss_subset are relative to the object `idxs`. We have to
            # return the aboluste idxs with respect to the `self` object.
            # TODO: It's returning a pair instead of a topk object
            return (loss_subset[0], idxs[loss_subset[1]])


    def plot_top_losses(self, k, largest=True, predicted=None,
                        actual=None, **kwargs):
        losses,idx = self.top_losses(k, largest, predicted, actual)
        if not isinstance(self.inputs, tuple): tupled_inputs = (self.inputs,)
        if isinstance(tupled_inputs[0], torch.Tensor): inps = tuple(o[idx] for o in tupled_inputs)
        else: inps = self.dl.create_batch(self.dl.before_batch([tuple(o[i] for o in tupled_inputs) for i in idx]))
        b = inps + tuple(o[idx] for o in (self.targs if is_listy(self.targs) else (self.targs,)))
        x,y,its = self.dl._pre_show_batch(b, max_n=k)
        b_out = inps + tuple(o[idx] for o in (self.decoded if is_listy(self.decoded) else (self.decoded,)))
        x1,y1,outs = self.dl._pre_show_batch(b_out, max_n=k)
        if its is not None:
            plot_top_losses(x, y, its, outs.itemgot(slice(len(inps), None)), self.preds[idx], losses,  **kwargs)

    def get_error_idxs(self):
        r"""
        Get the indices (relative to the validation set) of the items wrongly
        classified by the learner, regardless the type of error.
        """
        mc = self.most_confused()
        error_idxs = [self.top_losses(predicted=x[1], actual=x[0])[1].numpy() for x in mc]
        return np.concatenate(error_idxs).squeeze()

# Cell
def log_last_plt_as_wandb_img(plt_title):
    r"""
    Log the last plot made by matplotlib as a wandb.Image object. This is useful
    when you want to log plots that are not correctly transformed by wandb.log.
    Otherwise, it is preferrable to use just wandb.log and create an interactive plot
    """
    plt.savefig('./tmp.png')
    wandb.log({plt_title: wandb.Image('./tmp.png')})
    os.system("rm tmp.png")

# Cell
def log_plotnine_as_wandb_img(plot, title):
    r"""Log plotnine object as a wandb.Image"""
    plot.save(filename='./tmp.png', verbose=False)
    wandb.log({title: wandb.Image('./tmp.png')})
    os.system("rm tmp.png")