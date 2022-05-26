FROM ubuntu:22.04

MAINTAINER rdemko2332@gmail.com

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
RUN wget http://trna.ucsc.edu/software/trnascan-se-2.0.5.tar.gz \
  && tar -zxvf trnascan-se-2.0.5.tar.gz \
  && rm trnascan-se-2.0.5.tar.gz \
  && cp /bin/cm* /usr/local/bin/ \
  && cd tRNAscan-SE-2.0 \
  && ./configure \
  && make \
  && make install \
  && mv tRNAscan-SE /usr/local/bin/ \
  && cd /usr/local/bin/ chmod +x tRNAscan-SE

WORKDIR /work