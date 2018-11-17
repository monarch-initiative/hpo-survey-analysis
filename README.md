### HPO survey analysis
A library for analyzing surveys utilizing the HPO


#### What's Included?
- phenom: a python library for analyzing hpo surveys
- Jupyter notebooks demoing various tasks
- Scripts for running analyses

#### Getting Started

##### virtualenv

    virtualenv venv -p /usr/bin/python3.6
    source venv/bin/activate
    pip install -r requirements.txt
    export PYTHONPATH=.:$PYTHONPATH

##### Running notebooks with docker

    docker build . --tag hpo-subset:latest
    docker run --rm -it -p 8888:8888 --hostname localhost -v $PWD/notebooks/:/hpo-subset/notebooks/ hpo-subset jupyter notebook notebooks --ip=0.0.0.0 --allow-root --no-browser
    Copy the URL and replace "(localhost or 127.0.0.1)" with either "localhost" or "127.0.0.1" (without quotes)


#### Features
- Utilities for generating derived profiles from gold standard annotations, either best match or noisy
- Analyzing semantic similarity, disease classification of derived HPO profiles
- Disease and/or phenotype class enrichment
- Utilities for generating subsets of terms
- Simulate survey responses (yes, no, absent)


#### Example Surveys of Interest
- Phenotypr: https://phenotypr.com/
- GenomeConnect: https://www.genomeconnect.org/
