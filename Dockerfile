FROM csweaver/cellxgene
RUN mkdir -p /app
WORKDIR /app
COPY . ./
RUN git submodule update --init ExpressionMatrix2
RUN npm install --prefix cellxgene/ cellxgene
RUN npm run --prefix cellxgene build
RUN mv cellxgene/build/* templates
RUN mkdir -p ExpressionMatrix2/Release-ubuntu16-python3
RUN cd ExpressionMatrix2/Release-ubuntu16-python3 && cmake ../src && make -j 3
EXPOSE 5000:5000
CMD ["python3", "application.py"]