# Notes
docker build -t mlchaos_image --build-arg USER_ID=$(id -u) --build-arg GROUP_ID=$(id -g) .
docker run -d  --name mlchaos --gpus all -v $(pwd):/home/user/work -v {{your_path_to_data}}:/home/user/data -p {{your_local_port}}:8888 mlchaos_image

- Normalization
I don't know what ype of normalization I should apply. In terms of the data, does it make sense to normalize by all samples? Are the coordinates x, y part of the same environment? What I am doing right now is to use just batch normalization per sample (i.e., not distinguishing between channels)

The problem is that, if I try any item-based normalization, I ge NaNs in the validation loss.

ACTION: try to normalize directly in the TSData object, without any ItemTransform.

-------------

28/07/2020: The trianing results in the third dataset are slightly worse thanin the second one.
