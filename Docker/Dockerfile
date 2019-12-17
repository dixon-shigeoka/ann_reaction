#FROM nvidia/cuda:10.1-base-ubuntu16.04
#FROM tensorflow/tensorflow:1.14.0-gpu-py3
FROM tensorflow/tensorflow:1.14.0-py3

RUN apt-get update -y && \
    apt-get install -y --no-install-recommends \
    gfortran \
    libgfortran3 \
    neovim \
    cmake \
    gcc \
    g++ \
    make \
    unzip \
    wget \
    build-essential \
    zlib1g-dev \
    libgl1-mesa-glx\
    libsm6

#RUN mkdir -p /usr/src/app
#WORKDIR /usr/src/app

#COPY requirements.txt ./

RUN pip install --upgrade pip
RUN pip install \
    scikit-learn \
    pandas \
    matplotlib \
    ipython \
    keras \
    scipy \
    tqdm \
    opencv-python \
    vtk \
    pyautogui

RUN mkdir /root/.keras
COPY Docker-keras.json /root/.keras/keras.json

#RUN ln -s /usr/include/libv4l1-videodev.h /usr/include/linux/videodev.h
#RUN mkdir ~/tmp
#RUN cd ~/tmp && wget https://github.com/Itseez/opencv/archive/3.1.0.zip && unzip 3.1.0.zip
#RUN cd ~/tmp/opencv-3.1.0 && cmake CMakeLists.txt -DWITH_TBB=ON \
#                                                  -DINSTALL_CREATE_DISTRIB=ON \
#                                                  -DWITH_FFMPEG=OFF \

#                                                  -DWITH_IPP=OFF \
#                                                  -DCMAKE_INSTALL_PREFIX=/usr/local
#RUN cd ~/tmp/opencv-3.1.0 && make -j2 && make install
