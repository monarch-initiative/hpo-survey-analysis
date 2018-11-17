    docker build . --tag owlsimv3:latest

    docker run -e JAVA_OPTS='-Xmx2g' -p 9000:8080 owlsimv3
    
    Go to http://localhost:9000/api/docs/
    
    Test query: http://localhost:9000/api/match/naive-bayes-fixed-weight-two-state?id=HP%3A0009909&id=HP%3A0004590&id=HP%3A0030087&id=HP%3A0100025&id=HP%3A0100154&id=HP%3A0006872&id=HP%3A0006291&id=HP%3A0012086&id=HP%3A0001888&id=HP%3A0002234&id=HP%3A0100037&id=HP%3A0011821&id=HP%3A0006682&id=HP%3A0007185&id=HP%3A0001810