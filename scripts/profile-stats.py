from statistics import mean, median
import argparse
import csv

parser = argparse.ArgumentParser(
        description='Given a subset of HPO terms and diseases, generates '
                    'derived annotations from the HPO disease to phenotype '
                    'annotations')
parser.add_argument('--annotations', '-a', type=str, required=True,
                    help='Cached gold standard disease phenotype annotations')


args = parser.parse_args()
disease2phen = {}

with open(args.annotations, 'r') as cache_file:
    reader = csv.reader(cache_file, delimiter='\t', quotechar='\"')
    for row in reader:
        if row[0].startswith('#'): continue
        (mondo_id, phenotype_id) = row[0:2]
        if mondo_id in disease2phen:
            disease2phen[mondo_id].append(phenotype_id)
        else:
            disease2phen[mondo_id] = [phenotype_id]

print( "median: ",
    median([
        len(pheno) for pheno in disease2phen.values()
    ])
)

print( "mean: ",
    mean([
        len(pheno) for pheno in disease2phen.values()
    ])
)