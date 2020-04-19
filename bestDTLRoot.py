#Python function to find best DTL rooting among a list of rooting possibilities.

#Input: file to species tree, file to list of gene trees
#Output: the gene tree with the most optimal rooting in newick format
#Requires: Ranger or Optroot?

import os
import subprocess
from getDTLCost import getDTLCost

def findBestDTLRooting(speciesTree, geneTrees):

    #Initialize gene tree data structures. geneCosts contains tuples with gene Tree and rec cost
    geneCosts = []

    for geneTree in geneTrees:
        geneCosts.append([geneTree, -1])

    for i in range(len(geneCosts)):
        geneTree = geneCosts[i][0]
        geneCosts[i][1] = getDTLCost(speciesTree, geneTree)
    
    #initialize Cost and return tree data structure
    minCost = -1
    optimalGeneTrees = []

    for i in range(len(geneCosts)):
        if geneCosts[i][1] < minCost or minCost == -1:
            minCost = geneCosts[i][1]
            optimalGeneTrees = [geneCosts[i][0]]
        elif geneCosts[i][1] == minCost:
            optimalGeneTrees.append(geneCosts[i][0])
    
    #Need to edit this. If multiple trees with an equal optimal cost, return the tree with the best MAD score. Else, return all optimal trees.
    print("Number of Optimal Trees: ", len(optimalGeneTrees))
    print("Reconciliation Cost: ", minCost)
    print("Now Printing Trees.")
    print("------------------------")
    for i in optimalGeneTrees:
        print(i)



findBestDTLRooting("speciesTreeTest.txt", ["(((((a,b),c),d),b),(e,f));", "((((a,b),c),d),(b,(e,f)));", "((((((a,b),c),d),b),e),f);"])

    