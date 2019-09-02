#FROM nvidia/cuda:10.1-base-ubuntu16.04
FROM tensorflow/tensorflow:latest-gpu-py3

RUN apt-get update -y && \
    apt-get install -y --no-install-recommends \
    gfortran \
    libgfortran3 \
    libifport.so.5

#RUN mkdir -p /usr/src/app
#WORKDIR /usr/src/app

#COPY requirements.txt ./

RUN pip install --upgrade pip
RUN pip install \
    scikit-learn \
    pandas \
    matplotlib

#COPY . /usr/src/app

