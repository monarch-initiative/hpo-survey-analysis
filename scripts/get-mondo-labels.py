from phenom.utils import owl_utils
from rdflib import Graph
import logging
import gzip

output = "./mondo_labels.tsv"
output_file = open(output, 'w')

mondo_graph = Graph()

# Previous cache made with 2018-08-03 version of mondo
logging.info("Loading MONDO")
mondo_graph.parse(gzip.open("../data/mondo.owl.gz", 'rb'), format='xml')
root = "MONDO:0000001"

logging.info("Getting classes")
all_diseases = owl_utils.get_descendants(mondo_graph, root)
for disease in all_diseases:
    output_file.write("{}\t{}\n".format(disease, owl_utils.label(disease, mondo_graph)))
