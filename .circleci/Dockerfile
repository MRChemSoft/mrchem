FROM ubuntu:18.04

# Set an encoding to make things work smoothly.
ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8

# Add a timestamp for the build. Also, bust the cache.
ADD https://now.httpbin.org/when/now /opt/docker/etc/timestamp

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get --yes -qq update \
 && apt-get --yes -qq upgrade \
 && apt-get --yes -qq install \
                      bzip2 \
                      cpio \
                      curl \
                      g++ \
                      gcc \
                      git \
                      libblas-dev \
                      libboost-all-dev \
                      libeigen3-dev \
                      libidn11-dev \
                      liblapack-dev \
                      libopenmpi-dev \
                      libpython3-all-dev \
                      openmpi-bin \
                      python3-dev \
                      python3-pip \
 && rm -rf /var/lib/apt/lists/*

# Run common commands
COPY run_commands /opt/docker/bin/run_commands
RUN /opt/docker/bin/run_commands

ENV PATH $PATH:/root/.local/bin:/opt/cmake/bin

ENV LD_LIBRARY_PATH /root/.local/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}

SHELL ["/bin/bash"]
