FROM ubuntu:xenial
RUN apt-get update
RUN apt-get install libboost-all-dev python3-all-dev libhdf5-dev libhdf5-cpp-11 python3-pip git build-essential --yes
RUN git clone --recurse-submodules https://github.com/chanzuckerberg/cellxgene-rest-api.git
RUN make -C cellxgene-rest-api/ExpressionMatrix2/Release-ubuntu16-python3
RUN pip3 install -r cellxgene-rest-api/requirements.txt
RUN pip3 install nose
RUN python3 application.py &

