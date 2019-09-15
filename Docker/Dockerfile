#FROM nvidia/cuda:10.1-base-ubuntu16.04
FROM tensorflow/tensorflow:latest-gpu-py3

RUN apt-get update -y && \
    apt-get install -y --no-install-recommends \
    gfortran \
    libgfortran3 \
    neovim

#RUN mkdir -p /usr/src/app
#WORKDIR /usr/src/app

#COPY requirements.txt ./

RUN pip install --upgrade pip
RUN pip install \
    scikit-learn \
    pandas \
    matplotlib \
    ipython \
    keras

RUN cat << EOF > FILE ~/.keras/keras.jason
{
    "image_data_format": "channels_last",
    "epsilon": 1e-07,
    "floatx": "float64",
    "backend": "tensorflow"
}
EOF

#COPY . /usr/src/app

