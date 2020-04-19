import os
import subprocess

#Python function to get a DTL cost

#input: path to species tree, gene tree in string format
#output: cost of reconciliation

#assumes Ranger-DTL is a.out in Ranger-DTL folder

def getDTLCost(specPath, geneTree):
    spec_tree_newick = ""
    with open(specPath) as specFile:
        spec_tree_newick = specFile.readline()
    specFile.close()
    tempTree = open("tempTree.newick", "w+")
    tempTree.write(spec_tree_newick)
    tempTree.write("\n")
    tempTree.write(geneTree)
    tempTree.close()

    os.system('Ranger-DTL/a.out -i tempTree.newick > out.txt')

    os.system('rm tempTree.newick')

    with open('out.txt') as recFile:
        recLine = recFile.readline()
        while recLine:
            if "The minimum reconciliation cost" in recLine:
                cost = recLine[recLine.index(":") + 2:recLine.index("(") - 1]
                cost = int(cost)
            recLine = recFile.readline()
    recFile.close()

    os.system('rm out.txt')

    return cost

#print(getDTLCost("speciesTreeTest.txt", "((((((a,b),c),d),b),e),f);"))
