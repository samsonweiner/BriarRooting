import sys
import os
import dendropy as dp
from dendropy.calculate import treecompare as tc

def rf_distance_from_orig(filename):
    """
    :param filename: file formatted as the output from bestDTLRoot.py
    :return: rf distance between output of dual brothers rooting and original tree. 0 means correct root.
    """
    with open(filename, "r") as fp:
        file_contents = fp.read()
    fp.close()

    info, data = file_contents.split("------------------------")
    info = info.split("\n")
    gene_tree_loc = info[1]

    # emily's original trees in slighly different location than samson's relative to what's in the sample file
    #gene_tree_loc = gene_tree_loc.replace("../", "../final-project-src/")

    # we're using the formatted for optroot trees as the originals here
    gene_tree_loc = gene_tree_loc.replace("/relaxed_trees/", "/formatted_for_optroot/")

    gene_tree_loc = gene_tree_loc.replace("Gene Tree: ", "")

    # read in the original gene tree.
    oTree = dp.Tree.get_from_path(gene_tree_loc, schema="newick")
    oTree.is_rooted = True

    # make tree from output
    cTree = dp.Tree.get_from_string(data, schema="newick")
    cTree.is_rooted = True

    # calculate rf distance. taxon_namespaces have to be the same between the two trees.
    cTree.migrate_taxon_namespace(oTree.taxon_namespace) # hopefully will not throw error
    return tc.symmetric_difference(oTree, cTree)

# running the script with a directory as input will run rf_distance_from_orig for all files in that tree and
# return stats about how often the trees were rooted correctly.

# USAGE:
# python3 rf_root_compare.py [directory]
# where directory only includes files that have been outputted by dual brothers rooting
# and look like the samples sam sent emily

dir = sys.argv[1]
file_list = os.listdir(dir)
n = len(file_list)
rf_dist_list = [0 for i in range(n)]
for i in range(n):
    rf_dist_list[i] = rf_distance_from_orig(dir+"/"+file_list[i])

output = list(zip(file_list, rf_dist_list))
print(output)

# statistics: percent correct out of number of output files in folder
print(rf_dist_list.count(0)/n)