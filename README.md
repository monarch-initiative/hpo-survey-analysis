### HPO survey analysis
A library for analyzing surveys utilizing the HPO


#### What's Included?
- phenom: a python library with modules for querying monarch remote services, working with ontologies,
seand performing semantic similarity and distance computations
- Jupyter notebooks demoing various tasks
- Scripts for generating synthetic patient data, clustering, and evaluation

#### Getting Started

##### virtualenv

    virtualenv venv -p /usr/bin/python3.7
    source venv/bin/activate
    pip install -r requirements.txt
    export PYTHONPATH=.:$PYTHONPATH
    
##### Adding Redland Python Bindings

A dockerfile is provided for building the Redland C dynamic library for operating
systems running glibc and the python bindings.  See http://librdf.org/

To build:

    cd redland
    # Build Docker image
    docker build . --tag redland:latest

    # Build C binary and python lib from source 
    docker run --name build-redland redland

    # Copy binary and python src files
    docker cp build-redland:/build/redland-bindings.tar.gz ./
    
    # Extract and copy into virtual environment
    tar -zxvf redland-bindings.tar.gz
    mv redland-bindings/* ../venv/lib/python3.7/site-packages/


##### Running notebooks with docker

    docker build . --tag hpo-subset:latest
    docker run --rm -it -p 8888:8888 --hostname localhost -v $PWD/notebooks/:/hpo-subset/notebooks/ hpo-subset jupyter notebook notebooks --ip=0.0.0.0 --allow-root --no-browser
    Copy the URL and replace "(localhost or 127.0.0.1)" with either "localhost" or "127.0.0.1" (without quotes)


#### Features
- Implementations of common semantic similarity algorithms
- Disease and/or phenotype class enrichment (fisher exact)
- Clustering phenotype profiles
- Utilities for generating derived profiles from hpo rare disease annotations
- Simulate patient data


#### Example Surveys of Interest
- Phenotypr: https://phenotypr.com/
- GenomeConnect: https://www.genomeconnect.org/
