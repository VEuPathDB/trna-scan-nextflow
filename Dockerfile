FROM ubuntu:22.04
MAINTAINER rdemko2332@gmail.com
WORKDIR /work

RUN apt-get -qq update --fix-missing
RUN apt-get install -y wget cmake gcc perl
RUN wget http://trna.ucsc.edu/software/trnascan-se-2.0.5.tar.gz   && tar -zxvf trnascan-se-2.0.5.tar.gz   && cd tRNAscan-SE-2.0   && ./configure   && make   && make install   && rm /work/trnascan-se-2.0.5.tar.gz

RUN mv /work/tRNAscan-SE-2.0/tRNAscan-SE /usr/local/bin/
RUN cd /usr/local/bin/ chmod +x tRNAscan-SE
RUN cd /work
RUN wget http://eddylab.org/infernal/infernal-1.1.3.tar.gz   && tar -zxvf infernal-1.1.3.tar.gz   && cd infernal-1.1.3   && ./configure   && make   && make install   && rm /work/infernal-1.1.3.tar.gz
