FROM csweaver/cellxgene
RUN mkdir -p /app
WORKDIR /app
COPY . ./
RUN git submodule update --init ExpressionMatrix2
RUN make -C ExpressionMatrix2/Release-ubuntu16-python3 -j 3
EXPOSE 5000:5000
CMD ["python3", "application.py"]