FROM ubuntu:xenial
ARG GOOGLE_CLIENT_ID
ARG GOOGLE_CLIENT_SECRET
ARG SECRET_KEY
ENV GOOGLE_CLIENT_ID=$GOOGLE_CLIENT_ID
ENV GOOGLE_CLIENT_SECRET=$GOOGLE_CLIENT_SECRET
ENV SECRET_KEY=$SECRET_KEY
RUN apt-get update
RUN apt-get install libboost-all-dev python3-all-dev libhdf5-dev libhdf5-cpp-11 python3-pip git build-essential --yes
RUN git clone --recurse-submodules https://github.com/chanzuckerberg/cellxgene-rest-api.git
RUN make -C cellxgene-rest-api/ExpressionMatrix2/Release-ubuntu16-python3
RUN pip3 install -r cellxgene-rest-api/requirements.txt
RUN pip3 install nose
EXPOSE 5000


