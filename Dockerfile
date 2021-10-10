FROM ubuntu:20.04
ENV DEBIAN_FRONTEND=noninteractive
SHELL ["/bin/bash", "-c"]
ARG username=ubuntu
ARG password=ubuntu
RUN apt-get update &&\
    apt-get -y install \
    build-essential \
    libcairo2-dev \
    libgsl-dev \
    libgtk-3-dev \
    automake \
    libtool \
    curl \
    sudo &&\
    adduser --gecos "" --disabled-password $username &&\
    chpasswd <<<"$username:$password" &&   usermod -aG sudo ubuntu &&   su - ubuntu &&\
    cd ~/ &&\
    curl -O http://ibsimu.sourceforge.net/libibsimu-1.0.6dev_9cbf04.tar.gz &&\
    tar zxvf libibsimu-1.0.6dev_9cbf04.tar.gz &&\
    cd libibsimu-1.0.6dev_9cbf04 &&\
    ./configure --prefix=/home/ubuntu &&\
    make &&\
    make check &&\ 
    make install 
COPY --chown=$username . /home/$username/



