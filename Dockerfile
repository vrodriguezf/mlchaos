FROM vrodriguezf/jupyterlab-cuda

# Add non-root user (call this with the specific UID and GIO of the host, to share permissions)
ARG USER_ID=1000
ARG GROUP_ID=1000
ARG USER=user

RUN addgroup --gid $GROUP_ID $USER
RUN adduser --disabled-password --gecos '' --uid $USER_ID --gid $GROUP_ID $USER

###
# Python packages
###

# PyPi
RUN pip install nbdev fastai2 seaborn papermill wandb

# Git
#RUN pip install git+https://github.com/fastai/fastai2.git
RUN pip install git+https://github.com/muellerzr/fastinference.git

# Editable packages
RUN mkdir /home/$USER/lib
RUN cd /home/$USER/lib \
    && git clone https://github.com/vrodriguezf/timeseries.git \
    && cd timeseries \ 
    && pip install -e .

# Environmental variables for wandb
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

### Install Github command line (TODO: Move this to the parent image)
RUN cd ~ \
    && curl -L https://github.com/cli/cli/releases/download/v0.11.1/gh_0.11.1_linux_amd64.deb -o gh.deb \
    && apt install ./gh.deb

# Change the ownership of the editable installs within the lib folder
RUN chown -R $USER:$USER /home/$USER/lib

# change default user to $USER
USER $USER
WORKDIR /home/$USER

