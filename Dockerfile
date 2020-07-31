FROM vrodriguezf/jupyterlab-cuda

# Python packages
# Git
RUN pip install git+https://github.com/fastai/fastai2.git
RUN pip install git+https://github.com/ai-fast-track/timeseries.git
# PyPi
RUN pip install nbdev fastinference seaborn papermill

# Add non-root user (call this with the specific UID and GIO of the host, to share permissions)
ARG USER_ID=1000
ARG GROUP_ID=1000
ARG USER=user

RUN addgroup --gid $GROUP_ID $USER
RUN adduser --disabled-password --gecos '' --uid $USER_ID --gid $GROUP_ID $USER
USER $USER
WORKDIR /home/$USER

