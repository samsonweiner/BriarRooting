import dendropy as dp
from dendropy.calculate import treecompare as tc

def rf_distance_dualbros_orig(dualbros_file, original_file):
    """
    :param dualbros_file: output file from dualbros rooting, contains output newick string rooted
    :param original_file: original file with correct rooting, must have same taxon names as dualbros output. newick.
    :return: rf distance between the two trees. 0 signifies correct root, >1 means incorrect root. >2 means something
                went pretty wrong. usually is 1 or 2, as i've seen.
    """
    with open(dualbros_file, "r") as fp:
        file_contents = fp.read()
    fp.close()

    info, data = file_contents.split("------------------------")

    # read in the original gene tree.
    o_tree = dp.Tree.get_from_path(original_file, schema="newick")
    o_tree.is_rooted = True

    # make tree from output
    c_tree = dp.Tree.get_from_string(data, schema="newick")
    c_tree.is_rooted = True

    # calculate rf distance. taxon_namespaces have to be the same between the two trees.
    c_tree.migrate_taxon_namespace(o_tree.taxon_namespace) # hopefully will not throw error
    return tc.symmetric_difference(o_tree, c_tree)


# testing out the function.
#print(rf_distance_dualbros_orig("sample_output/tree4.txt", "../final-project-src/AllSimulatedDatasets/R-025-HI-NR/formatted_for_optroot/tree4"))