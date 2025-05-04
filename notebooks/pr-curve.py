import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import auc
import gzip
from typing import Tuple, List
import matplotlib as mpl

# Set a professional style
plt.style.use('seaborn-v0_8-whitegrid')
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans']

annotation_file = '../data/synthetic/gold-standard.tsv.gz'
classes = set()
with gzip.open(annotation_file, 'rt') as annot_file:
    for line in annot_file:
        row = line.rstrip().split('\t')
        classes.add(row[0])
n_classes = len(classes)


def get_precision_recall(file) -> Tuple[List[float], List[float]]:
    precision = {disease: [] for disease in classes}
    recall = {disease: [] for disease in classes}

    # iterate over confusion matrix and calculate per disease
    with gzip.open(file, 'rt') as confusion_file:
        for line in confusion_file:
            disease, score, true_pos, false_pos, false_neg, true_neg = line.rstrip().split('\t')
            try:
                precision[disease].append(int(true_pos) / (int(true_pos) + int(false_pos)))
            except ZeroDivisionError:
                precision[disease].append(1.0)  # Convention: precision=1 when TP=FP=0

            try:
                recall[disease].append(int(true_pos) / (int(true_pos) + int(false_neg)))
            except ZeroDivisionError:
                recall[disease].append(0)

    micro_precision = [sum(x) / len(x) for x in zip(*precision.values())]
    micro_recall = [sum(x) / len(x) for x in zip(*recall.values())]

    # Sort by recall (typically x-axis in PR curves)
    sorted_points = sorted(zip(micro_recall, micro_precision))
    micro_recall, micro_precision = zip(*sorted_points)

    return list(micro_recall), list(micro_precision)


gc_confusion = '../data/confusion_matrices/phenodigm/gc-confusion.tsv.gz'
gold_confusion = '../data/confusion_matrices/phenodigm/gold-confusion.tsv.gz'
lay_confusion = '../data/confusion_matrices/phenodigm/lay-confusion.tsv.gz'

sim_recall, sim_precision = get_precision_recall(gold_confusion)
gc_recall, gc_precision = get_precision_recall(gc_confusion)
lay_recall, lay_precision = get_precision_recall(lay_confusion)

sim_auc = auc(sim_recall, sim_precision)
gc_auc = auc(gc_recall, gc_precision)
lay_auc = auc(lay_recall, lay_precision)

# Create figure with high resolution
fig, ax = plt.subplots(figsize=(12, 10), dpi=300)

# Define a professional color palette
colors = {
    'gold': '#E69F00',  # Gold
    'gc': '#1F77B4',  # Blue
    'lay': '#009E73',  # Green
    'random': '#D55E00',  # Orange-Red
    'grid': '#CCCCCC',  # Light Gray
    'background': '#FFFFFF'  # Very Light Gray
}

line_width = 3
alpha = 0.9
random_chance = 1 / 7344  # 0.000136

# Set background color
ax.set_facecolor(colors['background'])
fig.patch.set_facecolor('white')

# Plot the curves with markers at data points
ax.plot(sim_recall, sim_precision, '-', color=colors['gold'], linewidth=line_width, alpha=alpha,
        marker='o', markersize=8, markevery=max(1, len(sim_recall) // 10),
        label=f'Gold Simulated (AUPRC = {sim_auc:.3f})')

ax.plot(lay_recall, lay_precision, '-', color=colors['lay'], linewidth=line_width, alpha=alpha,
        marker='s', markersize=8, markevery=max(1, len(lay_recall) // 10),
        label=f'Lay Subset (AUPRC = {lay_auc:.3f})')

ax.plot(gc_recall, gc_precision, '-', color=colors['gc'], linewidth=line_width, alpha=alpha,
        marker='^', markersize=8, markevery=max(1, len(gc_recall) // 10),
        label=f'Genome Connect (AUPRC = {gc_auc:.3f})')

# Add random chance line with dashed style
ax.axhline(y=random_chance, color=colors['random'], linestyle='--', linewidth=line_width,
           label=f'Random Chance = 1/7,344 ({random_chance:.6f})')

# Improved title and axis labels
ax.set_title('Precision-Recall Performance Comparison', fontsize=22, fontweight='bold', pad=15)
ax.set_xlabel('Recall', fontsize=18, fontweight='bold', labelpad=10)
ax.set_ylabel('Precision', fontsize=18, fontweight='bold', labelpad=10)

# Set axis limits with some padding
ax.set_xlim([-0.02, 1.02])
ax.set_ylim([-0.02, 1.02])

# Better ticks
ax.set_xticks(np.arange(0, 1.01, 0.1))
ax.set_yticks(np.arange(0, 1.01, 0.1))
ax.tick_params(axis='both', which='major', labelsize=12)

# Add grid for readability
# ax.grid(True, linestyle='--', alpha=0.7, color=colors['grid'])

# Better legend with shadow and rounded corners
legend = ax.legend(loc='upper right', prop={'size': 14, 'weight': 'bold'},
                   frameon=True, fancybox=True, framealpha=0.95, shadow=True)
legend.get_frame().set_edgecolor('gray')

# Add a caption/note
plt.figtext(0.5, 0.01, f"Based on analysis of {n_classes} disease classes",
            ha="center", fontsize=12, fontstyle='italic')

# Adjust layout and save with higher quality
plt.tight_layout(rect=[0, 0.03, 1, 0.97])
plt.savefig('./precision-recall-curve.png', bbox_inches='tight', pad_inches=0.2, dpi=300)
plt.show()
