FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
 && apt-get -y install \
    autoconf \
    automake \
    bison \
    fastqc \
    flex \
    git \
    libfreetype6-dev \
    libxml2-dev \
    libtool \
    m4 \
    pkg-config \
    python3 \
    python3-dev \
    python3-pip \
    python3-venv \
    zlib1g-dev \
    curl \
 && rm -rf /var/cache/apt/*

## Python packages
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
COPY docker/requirements.txt /opt/requirements.txt
RUN pip3 install -r /opt/requirements.txt \
 && rm -rf /root/.cache/pip

RUN pip3 install git+https://github.com/hubmapconsortium/panhumanpy@ef95d738f615856b474a385c734edf9087b745c7

RUN curl -L -o /opt/pan-human-azimuth-crosswalk.csv https://cdn.humanatlas.io/digital-objects/ctann/pan-human-azimuth/v1.0/assets/pan-human-azimuth-crosswalk.csv

COPY bin/ /opt

WORKDIR /opt

CMD ["/bin/bash"]
