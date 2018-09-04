import numpy as np
import argparse
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.cluster import hierarchy
import logging
from scipy.spatial.distance import squareform

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main():
    """
    Cluster distance matrix with scipy.cluster.hierarchy
    """
    parser = argparse.ArgumentParser(description='description')
    parser.add_argument('--input', '-i', type=str, required=True,
                        help='Location of input file'
                             ' that contains the sim matrix as json')
    parser.add_argument('--label', '-l', type=str, required=False,
                        help='Location of id-label mapping file')
    parser.add_argument('--output', '-o', required=False, help='output directory')
    args = parser.parse_args()

    logger.info("loading matrix")
    matrix = np.loadtxt(args.input, delimiter=",")
    labels = [line.rstrip('\n').split('\t')[0] for line in open(args.label, 'r')]


    logger.info("clustering")
    Z = linkage(squareform(matrix), 'ward')

    logger.info("generating flat clusters")
    clusters = fcluster(Z, 2, 'distance')

    # Output clusters
    output = open(args.output, 'w')
    for disease_id, cluster_id in zip(labels, clusters):
        output.write("{}\t{}\n".format(disease_id, cluster_id))

    # Draw dendrogram
    plt.figure()

    dn = hierarchy.dendrogram(Z)

    hierarchy.set_link_color_palette(['m', 'c', 'y', 'k'])
    fig, axes = plt.subplots(1, 2, figsize=(8, 3))
    dn1 = hierarchy.dendrogram(Z, ax=axes[0], above_threshold_color='y',
                               orientation='top')
    dn2 = hierarchy.dendrogram(Z, ax=axes[1], above_threshold_color='#bcbddc',
                               orientation='right')
    hierarchy.set_link_color_palette(None)  # reset to default after use
    plt.show()

if __name__ == "__main__":
    main()