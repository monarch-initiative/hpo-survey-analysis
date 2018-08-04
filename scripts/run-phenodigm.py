from phenom.similarity.phenodigm import Phenodigm
from rdflib import Graph

hp_graph = Graph()
hp_graph.parse('http://purl.obolibrary.org/obo/hp.owl', format='xml')
root = "HP:0000118"

ic_fh = open("/path/to/ic-cache.tsv", "r")
ic_map = {}

for line in ic_fh.readlines():
    hpo_id, ic = line.rstrip("\n").split("\t")
    ic_map[hpo_id] = float(ic)

phenodigm = Phenodigm(hp_graph, root, ic_map)

pheno_profile1 = ["HP:0001595", "HP:0002360"]
pheno_profile2 = ["HP:0002219", "HP:0002360"]

print(phenodigm.symmetric_phenodigm(pheno_profile1, pheno_profile2))
print(phenodigm.phenodigm_compare(pheno_profile1, pheno_profile2))
print(phenodigm.phenodigm_compare(pheno_profile2, pheno_profile1))

phenodigm = Phenodigm(hp_graph, root, ic_map, True, 'ic')

print(phenodigm.symmetric_phenodigm(pheno_profile1, pheno_profile2))
print(phenodigm.phenodigm_compare(pheno_profile1, pheno_profile2))
print(phenodigm.phenodigm_compare(pheno_profile2, pheno_profile1))
