# docker run script

parent_dir=$(cd $(dirname $0)/..;pwd)
#sudo docker run --runtime=nvidia --rm -it -v $parent_dir:/share tensorflow/tensorflow:latest-gpu-py3
sudo docker run --runtime=nvidia --rm -it -v $parent_dir:/share ann_reaction
