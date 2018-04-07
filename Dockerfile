FROM csweaver/cellxgene

# Env vars expected when running containers from this image:
# * ENV - staging, etc.
# * CXG_API_BASE
# * CONFIG_FILE

# It also needs to be able to authenticate to AWS. In ECS this will come in the form of a
# task profile. When running elsewhere you should pass in some credentials via environment
# variables. You can do this with `docker run -e AWS_ACCESS_KEY_ID -e AWS_SECRET_ACCESS_KEY ...`
# if you have the variables in your env.

RUN mkdir -p /app
WORKDIR /app
COPY . ./
RUN wget https://github.com/chanzuckerberg/ExpressionMatrix2/releases/download/0.5.0/ExpressionMatrix2-0.5.0.tar
RUN tar -xf  ExpressionMatrix2-0.5.0.tar
RUN mv ExpressionMatrix2-0.5.0/build/Release-ubuntu16-python3/ExpressionMatrix2.so .
RUN rm -rf ExpressionMatrix2-0.5.0.tar ExpressionMatrix2-0.5.0
RUN git submodule update --init cellxgene
RUN npm install --prefix cellxgene/ cellxgene
RUN npm run --prefix cellxgene build
RUN mv cellxgene/build/* templates
EXPOSE 5000:5000

# Install chamber, for pulling secrets into the container.
# Needs a `secret_key` setting.
ADD https://github.com/segmentio/chamber/releases/download/v1.16.0/chamber-v1.16.0-linux-amd64 /bin/chamber
RUN chmod +x /bin/chamber

CMD chamber exec stp-$ENV-cellxgene -- python3 application.py