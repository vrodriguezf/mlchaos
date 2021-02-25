# Notes
docker build -t mlchaos_image --build-arg USER_ID=$(id -u) --build-arg GROUP_ID=$(id -g) .
docker run -d  --name mlchaos --gpus all -v $(pwd):/home/user/work -v {{your_path_to_data}}:/home/user/data -p {{your_local_port}}:8888 mlchaos_image
docker run -d  --name mlchaos --gpus all -v $(pwd):/home/user/work -v /home/victor/rclone_local/gdrive/Universidad/Proyectos/2020/MLchaos/:/home/user/data -p 8801:8888 mlchaos_image

- Normalization
I don't know what ype of normalization I should apply. In terms of the data, does it make sense to normalize by all samples? Are the coordinates x, y part of the same environment? What I am doing right now is to use just batch normalization per sample (i.e., not distinguishing between channels)

The problem is that, if I try any item-based normalization, I ge NaNs in the validation loss.

ACTION: try to normalize directly in the TSData object, without any ItemTransform.

-------------

28/07/2020: The training results in the third dataset are slightly worse than in the second one.

09/08/2020: 

* `plot_top_losses` is not compatible with a Normalization by_sample, because in a `by_sample` normalization, you need a call to `encodes` right before the `decodes` to make it work. Otherwise you don't have the information from the original samples. In the batch plotted in `plot_top_losses`, samples are not encoded as in `show_batch` and `show_results`. Instead, they are just taken from the attribute `interp.inputs`, which in turn, comes from a call to `learn.get_preds(with_input=True)`, which returns the inputs transformed.

* oguiza's does not seem compatible with the inclusion of a vocab in the dataset. One good way to combine the two repos is to adapt the transforms of oguiza (class TSTensor) to the TensorTS from ai-fast-track

* oguiza's `show_batch` works five times faster than `ai-fast-track`. The bottleneck must be in the creation of the `Datasets` object. oguiza's datasets take 5 seconds to be created, while 
ai-fast-track takes 100 milliseconds.

15/09/2020

cam.py: 498:             plt.scatter(t, tseries, cmap=cmap, c = acts)

Replace tseries by tseries[ch], where ch is the channel of interest

25/02/21
Update docker-compose to 1.28