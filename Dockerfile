FROM ubuntu:16.04 

WORKDIR /app
ADD . /app

RUN \
  apt-get update \
  && apt-get install -y build-essential

ENTRYPOINT /bin/bash
