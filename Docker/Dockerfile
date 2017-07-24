FROM ubuntu:16.04
ENV HOME /root

RUN apt-get update
RUN apt-get -y install build-essential python3-dev
RUN apt-get -y upgrade

RUN pip3 install --user whatshap
RUN echo "export PATH=$HOME/.local/bin:$PATH" >> ~/.bashrc

RUN git clone https://github.com/AlgoLab/HapCHAT
RUN mkdir data
