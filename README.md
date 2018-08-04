### HPO subset analysis
A library for analyzing subsets of the HPO


#### What's Included?
- phenom: a python library for analyzing hpo subsets
- Jupyter notebooks demoing various tasks
- Scripts for running analyses

#### Getting Started

##### virtualenv

    virtualenv venv -p /usr/bin/python3.6
    source venv/bin/activate
    pip install -r requirements.txt
    export PYTHONPATH=.:$PYTHONPATH


#### Features
- Utilities for generating derived profiles from gold standard annotations, either best match or noisy
- Analyzing semantic similarity, disease classification of derived HPO profiles
- Disease and/or phenotype class enrichment
- Utilities for generating subsets of terms
- Simulate survey responses (yes, no, absent)

#### Driving Use Cases
- Prioritize new terms to add to a subset
- Analyze clinical utility of a survey that exposes some subset of the HPO

#### Example Subsets of Interest
- Lay person: Set of terms containing one or more lay person synonyms
- GenomeConnect: Set of terms exposed by the GenomeConnect survey
- OWL Axiomized: Set of terms that contain logical definitions
