import dendropy
import subprocess
import sys
import json


# a filter to find the nodes whose label contains "CANDIDATE", as well as removes the "CANDIDATE" string from them.
def candidate_filter(n):

    if n is None:
        return False

    elif n.label is None:
        if n.taxon is not None:
            if "CANDIDATE" in n.taxon.label:
                n.taxon.label = n.taxon.label.replace("CANDIDATE", "")
                return True

    else:
        if "CANDIDATE" in n.label:
            n.label = n.label.replace("CANDIDATE", "")
            return True

    return False


def mad_results(filename):
    args = ["python3", "mad.py", filename, "-f"]
    subprocess.call(args)

    # this will NOT write normal mad output to the original location of the tree.rooted UNLESS uncomment line 655 of mad.py
    # Will write another file which contains all of the different trees (newick strings) with "CANDIDATE" at
    # the top 10% of root locations. I didn't include the branch length of the root node in the newick string results.
    # writes to "top_mad_trees.json"

    with open("top_mad_trees.json", "r") as fp:
        target_trees = json.load(fp)
    fp.close()

    for current_tree in target_trees:

        nwstr = current_tree["newick"]

        # create the tree object from dendropy
        ctree = dendropy.Tree.get_from_string(src=nwstr, schema="newick", rooting="force-rooted")
        ctree.is_rooted=True

        # find the target root node
        target_root = ctree.find_node(candidate_filter)

        # re-root the tree at the candidate node
        if target_root is not None: # the case at the top AD score.
            ctree.reroot_at_edge(
                target_root.edge,
                target_root.edge_length/2,
                target_root.edge_length/2,
                update_bipartitions=True
            )

        # THIS IS PROCESSING TO CHANGE LEAF LABELS TO MATCH SPECIES TREE # Should be removed for general case
        for node in ctree.leaf_node_iter():
            cname = node.taxon.label
            numbers = cname[1:]
            numbers = numbers.split(" ")
            cname = "H" + numbers[1]+" "+numbers[0]
            node.taxon.label = cname

        # we have the re-rooted tree! yay. overwrite the old newick string in the dictionary
        current_tree["newick"] = ctree.as_string("newick")

    return target_trees
