FROM csweaver/cellxgene
ARG SECRET_KEY
ENV SECRET_KEY=$SECRET_KEY
RUN git clone --recurse-submodules https://github.com/chanzuckerberg/cellxgene-rest-api.git
RUN make -C cellxgene-rest-api/ExpressionMatrix2/Release-ubuntu16-python3 -j 3
EXPOSE 5000:5000


