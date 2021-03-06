from phenom.similarity.semantic_sim import SemanticSim
from phenom.similarity.semantic_dist import SemanticDist
from rdflib import Graph

#pheno_profile1 = ["HP:0001595", "HP:0002360", "-HP:0002814"]
#pheno_profile2 = ["HP:0002219", "HP:0002360", "-HP:0007340"]

pheno_profile1 = ["HP:0001595", "HP:0002360", "HP:0002814"]
pheno_profile2 = ["HP:0002219", "HP:0002360", "HP:0007340"]


hp_graph = Graph()
hp_graph.parse('http://purl.obolibrary.org/obo/hp.owl', format='xml')
root = "HP:0000118"

ic_fh = open("../examples/ic-cache.tsv", "r")
ic_map = {}

for line in ic_fh.readlines():
    hpo_id, ic = line.rstrip("\n").split("\t")
    ic_map[hpo_id] = float(ic)

sem_sim = SemanticSim(hp_graph, root, ic_map)
sem_dist = SemanticDist(hp_graph, root, ic_map)


print("Symmetric phenodigm: ", sem_sim.phenodigm_compare(
    pheno_profile1, pheno_profile2, is_symmetric=True))
print("Symmetric owlsim2: ", sem_sim.phenodigm_compare(
    pheno_profile1, pheno_profile2, is_symmetric=True, sim_measure='ic'))
print("jaccard sim: ", sem_sim.jaccard_sim(pheno_profile1, pheno_profile2))
print("sim gic: ", sem_sim.sim_gic(pheno_profile1, pheno_profile2))
print("resnik sim: ", sem_sim.resnik_sim(pheno_profile1, pheno_profile2))
print("resnik sim normalized: ", sem_sim.resnik_sim(
    pheno_profile1, pheno_profile2, is_normalized=True))
print("resnik sim symmetric norm: ", sem_sim.resnik_sim(
    pheno_profile1, pheno_profile2, is_normalized=True, is_symmetric=True))
print("resnik sim symmetric norm avg: ", sem_sim.resnik_sim(
    pheno_profile1, pheno_profile2, matrix_metric='avg', is_normalized=True, is_symmetric=True))
print("resnik sim symmetric norm max: ", sem_sim.resnik_sim(
    pheno_profile1, pheno_profile2, matrix_metric='max', is_normalized=True, is_symmetric=True))
print("cosine sim (neg weight .05): ", sem_sim.cosine_sim(pheno_profile1, pheno_profile2, negative_weight=.05))
print("cosine sim (neg weight .01): ", sem_sim.cosine_sim(pheno_profile1, pheno_profile2, negative_weight=.01))
print("cosine sim weighted: ", sem_sim.cosine_sim(pheno_profile1, pheno_profile2, ic_weighted=True))
print("euclidean distance: ", sem_dist.euclidean_distance(pheno_profile1, pheno_profile2))
print("euclidean matrix distance: ", sem_dist.euclidean_matrix(pheno_profile1, pheno_profile2))
print("euclidean matrix distance, jc dist: ", sem_dist.euclidean_matrix(pheno_profile1, pheno_profile2, distance_measure='jin_conrath'))
