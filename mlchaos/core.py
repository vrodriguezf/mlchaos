# AUTOGENERATED! DO NOT EDIT! File to edit: nbs/core.ipynb (unless otherwise specified).

__all__ = ['StandardizeGlobal', 'StandardizeItem', 'StandardizeNoDecode', 'TensorMotion', 'ToTensorMotion',
           'MotionBlock', 'TensorTSMotion', 'ToTensorTSMotion', 'prettify_plot_TensorTSMotion']

# Cell
from timeseries.all import *
from fastai.basics import *
from fastai.vision.data import get_grid
import numpy as np

# Cell
class StandardizeGlobal(ItemTransform):
    def setups(self, items):
        self.mean, self.std = get_mean_std(torch.stack(list(items)), scale_subtype='all_sam')
        print(self.mean, self.std)
    def encodes(self, x):
        return (x-self.mean)/self.std
    def decodes(self, x):
        f = to_cpu if x.device.type=='cpu' else noop
        return (x*f(self.std) + f(self.mean))

# Cell
class StandardizeItem(ItemTransform):
    "Standardize item per sample per channel"
    def encodes(self, x):
        self.mean = x.mean(axis=1).unsqueeze(1).repeat(1, x.shape[1])
        self.std = x.std(axis=1).unsqueeze(1).repeat(1, x.shape[1])
        return (x - self.mean)/ self.std
    def decodes(self, x):
        f = to_cpu if x.device.type=='cpu' else noop
        return (x*f(self.std) + f(self.mean))

# Cell
class StandardizeNoDecode(Standardize):
    def decodes(self, x:TensorTS):
        return x

# Cell
@typedispatch
def plot_top_losses(x:TensorTS, y:TensorCategory, samples, outs, raws, losses, nrows=None, ncols=None, figsize=None, **kwargs):
    axs = get_grid(len(samples), nrows=nrows,
                   ncols=ncols, add_vert=1,
                   figsize=figsize, title='Prediction/Actual/Loss/Probability')
    for ax,s,o,r,l in zip(axs, samples, outs, raws, losses):
        s[0].show(ctx=ax, **kwargs)
        ax.set_title(f'{o[0]}/{s[1]} / {l.item():.2f} / {r.max().item():.2f}')

# Cell
class TensorMotion(TensorTS):
    def show(self, ctx=None, title=None, label=None, chs=None,
             leg=True, ylim=None, return_fig=False, mode=None, **kwargs):
        "Show method with different modes"
        if mode=='ts':
            return TensorTS.show(self, ctx, title, chs, leg, **kwargs)
            if return_fig: print('Figures cannot be returned in mode ts')
        elif mode=='stacked':
            ret = self.show_stacked(ctx, title, label, chs, leg, **kwargs)
            if return_fig: return ret
        else:
            ret = self.show_poincmap(ctx, title, label, leg, ylim, return_fig, **kwargs)
            if return_fig: return ret

    def show_poincmap(self, ctx=None, title=None, label=None, leg=True,
                      ylim=None, return_fig=False, **kwargs):
        "Display poincare map for a motion"
        if ctx is None: fig, ctx = plt.subplots()
        t = range(self.shape[1])
        ctx.scatter(self[0][1:], self[1][1:])
        # The initial conditions (x0, y0) are plotted in a different colour
        ctx.scatter(self[0][0], self[1][0])
        ctx.set(xlim=[0, 360])
        if ylim: ctx.set(ylim=ylim)
        if title: ctx.set_title(title)
        if return_fig: return fig

    def show_stacked(self, ctx=None, title=None, label=None, chs=None, leg=True, color=None, **kwargs):
        "Display timeseries plots for a motion"
        t = range(self.shape[1])
        if ctx is None:
            fig , ctx = plt.subplots(nrows=2, ncols=1)
            ctx[0].plot(t, self[0], label='x', color='C0')
            ctx[0].scatter(t, self[0], color='C0', marker='.')
            ctx[1].plot(t, self[1], label='y', color='C1')
            ctx[1].scatter(t, self[1], color='C1', marker='.')
            if leg:
                ctx[0].legend(loc='upper right', ncol=2, framealpha=0.5)
                ctx[1].legend(loc='upper right', ncol=2, framealpha=0.5)
                if title: ax0.set_title(title)
            return fig
        else:
            # Call from show_batch, only plot one channel. The title is plotted in the show batch
            ch = chs[0]
            ctx.plot(t, self[ch], label=label, color=color)
            if leg: ctx.legend(loc='upper right', ncol=2, framealpha=0.5)
            if title: ctx.set_title(title)

# Cell
class ToTensorMotion(ItemTransform):
    # "x : 2D numpy array"
    def encodes(self, x): return TensorMotion(x)
    #def decodes(self, x): return x.numpy()

# Cell
def MotionBlock():
    "`TransformBlock` for orbit motions. Transform np array to TensorMotion type"
    return TransformBlock(type_tfms=[ToTensorTS(), ToTensorMotion()])

# Cell
@typedispatch
def plot_top_losses(x:TensorMotion, y:TensorCategory, samples, outs, raws, losses, nrows=None, ncols=None, figsize=None, **kwargs):
    axs = get_grid(len(samples), nrows=nrows,
                   ncols=ncols, add_vert=1,
                   figsize=figsize, title='Prediction/Actual/Loss/Probability')
    for ax,s,o,r,l in zip(axs, samples, outs, raws, losses):
        s[0].show(ctx=ax, **kwargs)
        ax.set_title(f'{o[0]}/{s[1]} / {l.item():.2f} / {r.max().item():.2f}')

# Cell
class TensorTSMotion(TensorTS):
    "This class could be deprecated soon. Please see `TensorMotion`"
    def show(self, ctx=None, title=None, label=None, chs=None, leg=True, color=None, **kwargs):
        "Display timeseries plots for a motion"
        t = range(self.shape[1])
        if ctx is None:
            _, ctx = plt.subplots(nrows=2, ncols=1)
            ctx[0].plot(t, self[0], label='x', color='C0')
            ctx[1].plot(t, self[1], label='y', color='C1')
            if leg:
                ctx[0].legend(loc='upper right', ncol=2, framealpha=0.5)
                ctx[1].legend(loc='upper right', ncol=2, framealpha=0.5)
                if title: ax0.set_title(title)
        else:
            # Call from show_batch, only plot one channel. The title is plotted in the show batch
            ch = chs[0]
            ctx.plot(t, self[ch], label=label, color=color)
            if leg: ctx.legend(loc='upper right', ncol=2, framealpha=0.5)
            if title: ctx.set_title(title)

# Cell
class ToTensorTSMotion(ItemTransform):
    # "x : 2D numpy array"
    def encodes(self, x): return TensorTSMotion(x)
    def decodes(self, x): return x.numpy()

# Cell
def prettify_plot_TensorTSMotion(ctxs, n, nchannels, nrows, ncols):
    "Configure how to paint axes with motion plots. ctxs is a 2xN array of axes"
    for i in range(n):
        ctxs[0 + nchannels*int(i/ncols)][i%ncols].get_xaxis().set_visible(False)
        ctxs[0 + nchannels*int(i/ncols)][i%ncols].spines['bottom'].set_visible(False)
        ctxs[1 + nchannels*int(i/ncols)][i%ncols].spines['top'].set_visible(False)

# Cell
@typedispatch
def show_batch(x:TensorTSMotion, y, samples, ctxs=None, max_n=3, nrows=1, figsize=(15,6), **kwargs):
    "Show batch for TensorTS objects"
    n = min(len(samples), max_n)
    ncols = int(np.ceil(n/nrows))
    nchannels = 2
    if ctxs is None: fig, ctxs = plt.subplots(nchannels*nrows, ncols, figsize=figsize)
    for b, tar, i in zip(samples.itemgot(0), samples.itemgot(1), range(n)):
        b.show(ctx=ctxs[0 + nchannels*int(i/ncols)][i%ncols], title = tar, label='x', chs=[0], color='C0')
        b.show(ctx=ctxs[1 + nchannels*int(i/ncols)][i%ncols], label='y', chs=[1], color='C1')
    prettify_plot_TensorTSMotion(ctxs, n, nchannels, nrows, ncols)
    fig.tight_layout()
    return ctxs

# Cell
@typedispatch
def show_results(x:TensorTSMotion, y, samples, outs, ctxs=None, max_n=3, nrows=1, figsize=(15,6), **kwargs):
    "Show results for TensorTSMotion objects"
    n = min(len(samples), max_n)
    ncols = int(np.ceil(n/nrows))
    nchannels = 2
    if ctxs is None: fig, ctxs = plt.subplots(nchannels*nrows, ncols, figsize=figsize)
    outs = [detuplify(o) for o in outs]
    for b, o, i in zip(samples, outs, range(max_n)):
        b[0].show(ctx=ctxs[0 + nchannels*int(i/ncols)][i%ncols], title=f'{o} / {b[1]}', label='x', chs=[0], color='C0')
        b[0].show(ctx=ctxs[1 + nchannels*int(i/ncols)][i%ncols], label='y', chs=[1], color='C1')
    prettify_plot_TensorTSMotion(ctxs, n, nchannels, nrows, ncols)
    fig.tight_layout()
    return ctxs

# Cell
@typedispatch
def plot_top_losses(x:TensorTSMotion, y:TensorCategory, samples, outs, raws, losses, nrows=1, figsize=None, **kwargs):
    #axs = get_grid(len(samples), nrows=nrows, ncols=ncols, add_vert=1, figsize=figsize, title='Prediction/Actual/Loss/Probability')
    n = len(samples)
    ncols = int(np.ceil(n/nrows))
    nchannels = 2
    fig, ctxs = plt.subplots(nchannels*nrows, ncols, figsize=figsize)
    for s,o,r,l,i in zip(samples, outs, raws, losses, range(len(samples))):
        s[0].show(ctx=ctxs[0 + nchannels*int(i/ncols)][i%ncols], label='x', chs=[0], color='C0', **kwargs)
        s[0].show(ctx=ctxs[1 + nchannels*int(i/ncols)][i%ncols], label='y', chs=[1], color='C1', **kwargs)
        ctxs[0 + nchannels*int(i/ncols)][i%ncols].set_title(f'{o[0]}/{s[1]} / {l.item():.2f} / {r.max().item():.2f}')
    prettify_plot_TensorTSMotion(ctxs, n, nchannels, nrows, ncols)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.suptitle('Prediction/Actual/Loss/Probability')