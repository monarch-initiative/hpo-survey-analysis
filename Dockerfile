FROM python:3.6-alpine

EXPOSE 8888

RUN apk update && apk upgrade
RUN apk add --update gcc g++ && rm -rf /var/cache/apk/*

ADD . /hpo-subset
WORKDIR /hpo-subset

RUN pip install -r requirements.txt

ENV PATH="/hpo-subset/scripts/:/hpo-subset/:$PATH"
