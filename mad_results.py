import dendropy
import subprocess
import sys

# first run mad rooting on the input tree. run with the -f or -g options.
# subprocess with args: python3 mad.py [input tree] -f, input tree is a path to a file containing a newick string.

# we'll have to do some debugging and writing instructional stuff here. but anyways, here's the code without any error protection.
args = ["python3", "mad.py", sys.argv[1], "-f"]
subprocess.call(args)

# this will write two files: the normal mad output to the original location of the tree.rooted
# also will write another file which contains all of the different trees (newick strings) with "CANDIDATE" at
# the top 10% of root locations. I didn't include the branch length of the root node in the newick string results.


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


with open("top_mad_trees.txt", "r") as fp:
    root_candidates = fp.read().splitlines()
fp.close()

for i in range(len(root_candidates)):

    nwstr = root_candidates[i]

    # create the tree object from dendropy
    ctree = dendropy.Tree.get_from_string(src=nwstr, schema="newick", rooting="force-rooted")
    ctree.is_rooted=True

    # find the target root node
    target_root = ctree.find_node(candidate_filter)

    # re-root the tree at the candidate node
    if target_root is not None: # the case at the top AD score.
        ctree.reroot_at_edge(target_root.edge, target_root.edge_length/2, target_root.edge_length/2, update_bipartitions=True)

    # we have the re-rooted tree! yay. overwrite the old newick string
    root_candidates[i] = ctree.as_string("newick")

# it's likely these will need a little more manipulation before they go into DTL-Ranger. Will find out more later.
with open("rerooted_mad_trees.txt", "w") as fp:
    fp.write(";\n".join(root_candidates) + ";")
fp.close()