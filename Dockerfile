FROM vrodriguezf/jupyterlab-cuda

# Add non-root user (call this with the specific UID and GIO of the host, to share permissions)
ARG USER_ID=1000
ARG GROUP_ID=1000
ARG USER=user

RUN addgroup --gid $GROUP_ID $USER
RUN adduser --disabled-password --gecos '' --uid $USER_ID --gid $GROUP_ID $USER

# Copy the default jupyterlab settings user settings to the new user folder
RUN cp -r /.jupyter /home/$USER/.jupyter
RUN chown -R $USER_ID:$GROUP_ID /home/$USER/.jupyter

###
# Python packages
###

# PyPi
RUN pip install --upgrade nbdev seaborn papermill plotly wandb papersweep plotnine
RUN pip3 install torch==1.9.0+cu111 torchvision==0.10.0+cu111 -f https://download.pytorch.org/whl/torch_stable.html

# Git
ENV LANG C.UTF-8
RUN pip install git+https://github.com/fastai/fastai \ 
                git+https://github.com/muellerzr/fastinference.git \
                git+https://github.com/vrodriguezf/papersweep.git

# Editable packages
RUN mkdir /home/$USER/lib
RUN cd /home/$USER/lib \
    && git clone https://github.com/timeseriesAI/tsai.git \
    && cd tsai \ 
    && pip install -e .
RUN cd /home/$USER/lib \
    && git clone https://github.com/vrodriguezf/timeseries.git \
    && cd timeseries \ 
    && pip install -e .
    
# Environmental variables for wandb
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

# Change the ownership of the editable installs within the lib folder
RUN chown -R $USER:$USER /home/$USER/lib

# change default user to $USER
USER $USER
WORKDIR /home/$USER

