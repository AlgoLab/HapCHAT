FROM ubuntu:16.04
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    g++ \
    git-core \
    python3-biopython-sql \
    python3-dev \
    python3-networkx \
    virtualenv \
    zlib1g-dev

VOLUME ["/data"]

RUN git clone https://github.com/AlgoLab/HapCHAT.git && \
    cd HapCHAT && \
    ./setup.sh


ENTRYPOINT ["./HapCHAT.py"]

CMD [""]
