FROM maven:3.6.0-jdk-8-alpine

EXPOSE 8080

WORKDIR /data

RUN apk --update add git

RUN git clone https://github.com/monarch-initiative/owlsim-v3.git /data/owlsim3/
RUN cd /data/owlsim3 && git checkout 8cce7d78 && mvn package -DskipTests -DskipITs

ADD resources/configuration.yaml /data/owlsim3/conf/
ADD resources/annotations.tsv /data/owlsim3/conf/
ADD resources/mondo_labels.tsv /data/owlsim3/conf/
ADD resources/hp.owl /data/owlsim3/conf/

CMD java -jar /data/owlsim3/owlsim-services/target/owlsim-services-3.0-SNAPSHOT.jar server /data/owlsim3/conf/configuration.yaml
