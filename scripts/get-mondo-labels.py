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
all_classes = set()
for disease in mondo_graph.subjects():
    if str(disease).startswith("http://purl.obolibrary.org/obo/MONDO_"):
        all_classes.add(disease)

for uri in all_classes:
    curie = str(uri).replace("http://purl.obolibrary.org/obo/MONDO_", "MONDO:")
    output_file.write("{}\t{}\n".format(curie, mondo_graph.label(uri)))
