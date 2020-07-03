FROM vrodriguezf/jupyterlab-cuda

# Python packages
RUN pip install nbdev fastai2
RUN pip install git+https://github.com/timeseriesAI/timeseriesAI.git
