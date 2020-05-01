import sys

from mad_results import mad_results
from bestDTLRoot import findBestDTLRooting

#geneTreePath = str(sys.argv[1])
#speciesTreePath = str(sys.argv[2])
print("Welcome to Briar-Rooting!")

print('Gene Tree:', geneTreePath)
print('Species Tree:', speciesTreePath)

print("Compiling Mad Scores...")

target_trees = mad_results(geneTreePath)

print("Done.")

print("Number of trees to check DTL Scores: ", len(target_trees) )
print("Compiling DTL scores...")

optimalGeneTrees = findBestDTLRooting(speciesTreePath, target_trees)
