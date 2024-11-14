FROM ubuntu:22.04

RUN apt-get update --fix-missing

# Installing Software
RUN apt-get install -y \
  wget \
  perl \
  infernal \
  infernal-doc \
  gcc \
  cmake

# Setting Up

WORKDIR /usr/local

RUN wget https://github.com/UCSC-LoweLab/tRNAscan-SE/archive/refs/tags/v2.0.12.tar.gz \
  && tar -zxvf v2.0.12.tar.gz \
  && rm v2.0.12.tar.gz \
    && cd tRNAscan-SE-2.0.12 \
  && ./configure \
  && make \
  && make install

WORKDIR /work
