# MLchaos
> Studying the use of machine learning techniques to characterize chaotic orbits

## Deploy

To run the notebooks, install `docker` and `docker-compose` in your system. 
Then, create a new *.env* file in the root of the project following the structure:
```
# The name of the docker-compose project
COMPOSE_PROJECT_NAME=your_project_name
# The user ID you are using to run docker-compose
USER_ID=your_numeric_id
# The group ID you are using to run docker-compose (you can get it with id -g in a terminal)
GROUP_ID=your_numeric_id
# The user name assigned to the user id
USER_NAME=your_user_name
# The port from which you want to access Jupyter lab
JUPYTER_PORT=XXXX
# The path to your data files to train/test the models
LOCAL_DATA_PATH=/path/to/your/data
# The W&B personal API key (see https://wandb.ai/authorize)
WANDB_API_KEY=your_wandb_api_key
# List of comma separated GPU indices that will be available in the container (by default only 0, the first one)
CUDA_VISIBLE_DEVICES=0
```

You'll also need to have a `.gitconfig` file in your home folder. It can be an empty file that you create manually, or it can contain your git global configuration. For the latter case, run:
- `git config --global user.name "YOUR NAME IN THIS GITLAB INSTANCE"`
- `git config --global user.email "YOUR EMAIL IN THIS GITLAB INSTANCE"`

This will automatically create the `~/.gitconfig` file in your home folder.

Finally, in a terminal located in the root of this repository, run:

```docker-compose up -d --build```

then go to `localhost:{{JUPYTER_PORT}}` to run the notebooks under the `nbs` folder. In case you are working in a remote server, replace `localhost` with the IP of your remote server.

> Note: This is meant to be deployed in a system with at least one GPU. If your system does not have any, you'll have to change the previous ``docker-compose`` command to: ``docker-compose -f docker-compose-cpu.yml up -d --build``

Finally, to get the access token for JupyterLab, execute `docker-compose logs jupyter` and copy-paste the token that will be prompted as part of the container logs.

## Cite

TBA, paper accepted at Nature Scientific Reports
