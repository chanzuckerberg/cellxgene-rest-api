FROM csweaver/cellxgene
ARG GOOGLE_CLIENT_ID
ARG GOOGLE_CLIENT_SECRET
ARG SECRET_KEY
ENV GOOGLE_CLIENT_ID=$GOOGLE_CLIENT_ID
ENV GOOGLE_CLIENT_SECRET=$GOOGLE_CLIENT_SECRET
ENV SECRET_KEY=$SECRET_KEY
RUN git clone --recurse-submodules https://github.com/chanzuckerberg/cellxgene-rest-api.git
RUN make -C cellxgene-rest-api/ExpressionMatrix2/Release-ubuntu16-python3 -j 3
EXPOSE 5000:5000


