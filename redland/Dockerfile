FROM python:3.7.2-slim

#RUN apk add gnupg libxml2-dev swig perl

RUN apt-get update && apt-get install -y \
    gcc \
    make \
    pkg-config \
    wget \
    gnupg \
    libxml2-dev \
    swig \
    perl

WORKDIR /build

ADD resources/build-redland.sh /build

ENTRYPOINT ["/build/build-redland.sh"]
