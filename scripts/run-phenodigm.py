from phenom.similarity.phenodigm import Phenodigm
from phenom.utils import owl_utils
from phenom.similarity import metric
from rdflib import Graph, RDFS

pheno_profile1 = ["HP:0001595", "HP:0002360"]
pheno_profile2 = ["HP:0002219", "HP:0002360"]


hp_graph = Graph()
hp_graph.parse('http://purl.obolibrary.org/obo/hp.owl', format='xml')
root = "HP:0000118"

ic_fh = open("/home/kshefchek/ic-cache.tsv", "r")
ic_map = {}

for line in ic_fh.readlines():
    hpo_id, ic = line.rstrip("\n").split("\t")
    ic_map[hpo_id] = float(ic)

phenodigm = Phenodigm(hp_graph, root, ic_map)


print(phenodigm.symmetric_phenodigm(pheno_profile1, pheno_profile2))
print(phenodigm.phenodigm_compare(pheno_profile1, pheno_profile2))
print(phenodigm.phenodigm_compare(pheno_profile2, pheno_profile1))

phenodigm = Phenodigm(hp_graph, root, ic_map, True, 'ic')

print(phenodigm.symmetric_phenodigm(pheno_profile1, pheno_profile2))
print(phenodigm.phenodigm_compare(pheno_profile1, pheno_profile2))
print(phenodigm.phenodigm_compare(pheno_profile2, pheno_profile1))


print(metric.profile_jaccard(pheno_profile1, pheno_profile2, hp_graph, root))

print([owl_utils.get_closure(hp_graph, pheno, RDFS['subClassOf'], root) for pheno in pheno_profile1])
