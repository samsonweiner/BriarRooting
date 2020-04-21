/*
 * *   Copyright (C) 2017 Mukul S. Bansal (mukul.bansal@uconn.edu).
 * *   Based on code originally written by Andre Wehe and Mukul S. Bansal
 * *
 * *   This program is free software: you can redistribute it and/or modify
 * *   it under the terms of the GNU General Public License as published by
 * *   the Free Software Foundation, either version 3 of the License, or
 * *   (at your option) any later version.
 * *
 * *   This program is distributed in the hope that it will be useful,
 * *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 * *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * *   GNU General Public License for more details.
 * *
 * *   You should have received a copy of the GNU General Public License
 * *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * */




#ifndef DTL_TREESET_H
#define DTL_TREESET_H


// ------------------------------------------------------------------------------------------------------------------
// primary mapping (unrooted)
template<class T>
class PrimaryMappingUnrooted {
protected:
	// mapping
	T* mapping[3];

public:
	PrimaryMappingUnrooted() {
		resetMapping();
	}

	// get mapping in direction i
	inline T*& getMapping(const int &i) {
		return mapping[i];
	}

	// set mapping in direction i
	template<class T2>
	inline void setMapping(int &i, T2* &node) {
		mapping[i] = node;
	}

	// reset mappings for all directions
	inline void resetMapping() {
		for (int i=0; i<3; i++) mapping[i] = NULL;
	}

	template<class T2>
	friend ostream & operator << (ostream & os, PrimaryMappingUnrooted<T2> & m);
};

// output the node into a string stream
template<class T>
ostream & operator << (ostream & os, PrimaryMappingUnrooted<T> & m) {
	os << "1st-map:{" << m.mapping[0] << ',' << m.mapping[1] << ','  << m.mapping[2] << '}';
	return os;
}

// ------------------------------------------------------------------------------------------------------------------
// primary mapping (rooted)
template<class T>
class PrimaryMappingRooted {
protected:
	// mapping
	T* mapping;

public:
	PrimaryMappingRooted() {
		resetMapping();
	}

	// get mapping in direction i
	inline T*& getMapping() {
		return mapping;
	}

	// set mapping in direction i
	template<class T2>
	inline void setMapping(T2* &node) {
		mapping = node;
	}

	// reset mappings for all directions
	inline void resetMapping() {
		mapping = NULL;
	}

	template<class T2>
	friend ostream & operator << (ostream & os, PrimaryMappingRooted<T2> & m);
};

// output the node into a string stream
template<class T>
ostream & operator << (ostream & os, PrimaryMappingRooted<T> & m) {
	os << "1st-map:" << m.mapping;
	return os;
}

// ------------------------------------------------------------------------------------------------------------------
// secondary mapping
template<class T>
class SecondaryMapping {
public:
	// mapping
	T* secondarymapping;

	SecondaryMapping() {
		secondarymapping = NULL;
	}

	// check if contained in Gamma-Tree
	inline bool belongs2GammaTree() {
		return secondarymapping != NULL;
	}

	template<class T2>
	friend ostream & operator << (ostream & os, SecondaryMapping<T2> & m);
};

// output the node into a string stream
template<class T>
ostream & operator << (ostream & os, SecondaryMapping<T> & m) {
	os << "2nd-map:" << m.secondarymapping;
	return os;
}

// ------------------------------------------------------------------------------------------------------------------
// order number
class Order {
public:
	unsigned int no, begin, end;

	Order() {
		no = 0;
		begin = 0;
		end = 0;
	}

	friend ostream & operator << (ostream & os, Order & m);
};

// output the node into a string stream
ostream & operator << (ostream & os, Order & m) {
	os << "order:(" << m.begin << '<' << m.no << '<' << m.end << ')';
	return os;
}

// ------------------------------------------------------------------------------------------------------------------
// order number
class GeneDuplication {
public:
	unsigned int gain, lost[2];
	unsigned int genedup, tempgenedup;

	GeneDuplication() {
		gain = 0;
		for (int i=0; i<2; i++) lost[i] = 0;
		genedup = 0;
	}

	friend ostream & operator << (ostream & os, GeneDuplication & m);
};

// output the node into a string stream
ostream & operator << (ostream & os, GeneDuplication & m) {
	os << "genedup:" << m.genedup << '(' << m.lost[0] << '[' << m.gain << ']' << m.lost[1] << ')';
	return os;
}

// ==================================================================================================================
// ==================================================================================================================
// ==================================================================================================================

// ------------------------------------------------------------------------------------------------------------------
class SepciesNode;
class NamedSepciesNode;
class GeneNodeRooted;
class GeneNodeUnrooted;
class NamedGeneNodeUnrooted;
class SpeciesTree;
class GeneTree;

// ------------------------------------------------------------------------------------------------------------------
// species node
class SpeciesNode : public TreeNodeRooted<SpeciesNode>, public Order, public GeneDuplication {
public:
	int idx;
	unsigned int constraint;
	
	int nodeDepth;
	double branchDepth;
//	int POidx; // index for storing post order based index 

//	int numleaves; // number of leaves in the forest at that species node.

	set<GeneNodeRooted*> MGDset; 
		
	SpeciesNode(SpeciesNode  *parent = NULL) : TreeNodeRooted<SpeciesNode>(parent) {
		constraint = 0;
	}

	friend ostream & operator << (ostream & os, SpeciesNode & m);
};

// output a node into a string stream
ostream & operator << (ostream & os, SpeciesNode & m) {
	os << (TreeNodeRooted<SpeciesNode> &)m << ", "
	   << "idx:" << m.idx << ", "
	   << (Order &)m << ", "
	   << (GeneDuplication &)m;
	return os;
}

// ------------------------------------------------------------------------------------------------------------------
// species node with name
class NamedSpeciesNode : public SpeciesNode, public NodeID {
public:
	NamedSpeciesNode(const string &id, SpeciesNode *parent = NULL) : SpeciesNode(parent), NodeID(id) {}

	friend ostream & operator << (ostream & os, NamedSpeciesNode & m);
};

// output a node into a string stream
ostream & operator << (ostream & os, NamedSpeciesNode & m) {
	os << (NodeID &)m << ", "
	   << (SpeciesNode &)m;
	return os;
}

// ------------------------------------------------------------------------------------------------------------------
// gene node (rooted)
class GeneNodeRooted : public TreeNodeRooted<GeneNodeRooted>, public PrimaryMappingRooted<SpeciesNode>, public SecondaryMapping<SpeciesNode> {
public:
	int idx;
	unsigned int treeid;
//	int POidx; // index for storing post order based index 

	GeneNodeRooted *leftMost, *rightMost; 

	using PrimaryMappingRooted<SpeciesNode>::getMapping;

	GeneNodeRooted(GeneNodeRooted *parent = NULL) : TreeNodeRooted<GeneNodeRooted>(parent) {}

	// get LCA mapping in rooted direction
	inline SpeciesNode*& getMapping() {
		return mapping;
	}

	// set LCA mapping in rooted direction
	template<class T2>
	inline void setMapping(T2* &node) {
		mapping = node;
	}

	// true if this node is a gene duplication
	bool isGeneduplication() {
		GeneNodeRooted *gn[2];
		getChildren(gn);
		for (int i = 0; i < 2; i++) {
			if (child(i) != NULL) {
				if (child(i)->getMapping() == getMapping()) return true;
			}
		}
		return false;
	}

	friend ostream & operator << (ostream & os, GeneNodeRooted & m);
};

// output the node into a string stream
ostream & operator << (ostream & os, GeneNodeRooted & m) {
	os << (TreeNodeRooted<GeneNodeRooted> &)m << ", "
	   << (PrimaryMappingRooted<SpeciesNode> &)m << ", "
	   << (SecondaryMapping<SpeciesNode> &)m;
	return os;
}

// ------------------------------------------------------------------------------------------------------------------
// gene node with name
class NamedGeneNodeRooted : public GeneNodeRooted, public NodeID {
public:
	
	string geneName;  // to store the gene name for that node when input gene tree leaf labels are of the form speciesname_genename.
	NamedGeneNodeRooted(const string &id, GeneNodeRooted *parent = NULL) : GeneNodeRooted(parent), NodeID(id) {}

	friend ostream & operator << (ostream & os, NamedGeneNodeRooted & m);
};

// output a node into a string stream
ostream & operator << (ostream & os, NamedGeneNodeRooted & m) {
	os << (NodeID &)m << ", " << (GeneNodeRooted &)m;
	return os;
}

// ------------------------------------------------------------------------------------------------------------------
// gene node (unrooted)
class GeneNodeUnrooted : public TreeNodeUnrooted<GeneNodeUnrooted>, public PrimaryMappingUnrooted<SpeciesNode>, public SecondaryMapping<SpeciesNode> {
public:
	using PrimaryMappingUnrooted<SpeciesNode>::getMapping;

	GeneNodeUnrooted(GeneNodeUnrooted *parent = NULL) : TreeNodeUnrooted<GeneNodeUnrooted>(parent) {}

	// get LCA mapping in rooted direction
	inline SpeciesNode*& getMapping() {
		return mapping[parentno];
	}

	// set LCA mapping in rooted direction
	template<class T2>
	inline void setMapping(T2* &node) {
		mapping[parentno] = node;
	}

	// return the LCA mappings (unrooted)
	inline SpeciesNode* &getMappingDirected(GeneNodeUnrooted* &parent) {
		for (int i=0; i<3; i++) {
			if (relative[i] == parent) {
				return mapping[i];
			}
		}
		EXCEPTION("direction not found in getMappingDirected");
		return mapping[0];
	}

	// set the LCA mapping (unrooted)
	template<class T2>
	inline void setMappingDirected(GeneNodeUnrooted* &parent, T2* &node) {
		for (int i=0; i<3; i++) {
			if (relative[i] == parent) {
				mapping[i] = node;
				return;
			}
		}
		EXCEPTION("direction not found in setMappingDirected");
	}

	friend ostream & operator << (ostream & os, GeneNodeUnrooted & m);
};

// output the node into a string stream
ostream & operator << (ostream & os, GeneNodeUnrooted & m) {
	os << (TreeNodeUnrooted<GeneNodeUnrooted> &)m << ", "
	   << (PrimaryMappingUnrooted<SpeciesNode> &)m << ", "
	   << (SecondaryMapping<SpeciesNode> &)m;
	return os;
}

// ------------------------------------------------------------------------------------------------------------------
// gene node with name
class NamedGeneNodeUnrooted : public GeneNodeUnrooted, public NodeID {
public:
	
	string geneName;  // to store the gene name for that node when input gene tree leaf labels are of the form speciesname_genename.
	NamedGeneNodeUnrooted(const string &id, GeneNodeUnrooted *parent = NULL) : GeneNodeUnrooted(parent), NodeID(id) {}

	friend ostream & operator << (ostream & os, NamedGeneNodeUnrooted & m);
};

// output a node into a string stream
ostream & operator << (ostream & os, NamedGeneNodeUnrooted & m) {
	os << (NodeID &)m << ", " << (GeneNodeUnrooted &)m;
	return os;
}

// ==================================================================================================================
// ==================================================================================================================
// ==================================================================================================================

// ------------------------------------------------------------------------------------------------------------------
class SpeciesTree : public TreeIO<Tree<SpeciesNode, NamedSpeciesNode> > {
protected: unsigned int namecounter;
public:
	SpeciesTree() {
		namecounter = 1;
	}

	virtual ~SpeciesTree() {
	}

	// create a species node from the input - analyses comments for [&CONSTRAINT]
	unsigned int constraintcounter;
	void createNode(TempNode &n) {
		// named node or not
		if (n.name.empty()) {
			ostringstream os;
			os << 'n' << namecounter++;
			string name = os.str();
			n.node = new NAMEDNODE(name, NULL);
		} else {
			n.node = new NAMEDNODE(n.name, NULL);
		}

		// checking for constraints
		unsigned int c = constraintcounter;
		for (int i=0, last=n.comments.size(); i<last; i++) {
			string &str = n.comments[i];
			if (str == "[&CONSTRAINT]") {
				constraintcounter = c + 1;
				n.node->constraint = constraintcounter;
			}
		}
	}

// 	// output a named tree node in newick format
// 	void namednode2newick(ostream &os, NamedSpeciesNode &node) {
// 		os << name2newick(node.getName());
// 		if (node.constraint) os << "[&CONSTRAINT"<<node.constraint<<"]";
// 	}
//
// 	// output a unnamed tree node in newick format
// 	void node2newick(ostream &os, SpeciesNode &node) {
// 		if (node.constraint) os << "[&CONSTRAINT"<<node.constraint<<"]";
// 	}

	// color each node of a constrained subtree
	inline void colorSpeciesTreeByConstraints() {
		colorSpeciesTreeByConstraintsDFS(root, root->constraint);
	}
	inline void colorSpeciesTreeByConstraintsDFS(SpeciesNode *node, unsigned int &color) {
		if (node == NULL) return;
		if (node->constraint == 0) {
			node->constraint = color;
		}
		for (int i=0, last=2; i<last; i++) {
			colorSpeciesTreeByConstraintsDFS(node->child(i), node->constraint);
		}
	}

	// set the gene duplication triple to 0
	void resetGeneDubTriple() {
		for (vector<SpeciesNode*>::iterator itr = nodes.begin(); itr != nodes.end(); itr++) {
			SpeciesNode &t = **itr;
			t.gain = 0;
			t.lost[0] = 0;
			t.lost[1] = 0;
		}
	}

	// assigns each node an index
	void assignIndex() {
		for (int i = 0; i < nodes.size(); i++) nodes[i]->idx = i;
	}

	// node numbering + range of subtree
	inline void establishOrder() {
		int pos = 0;
		establishOrderDFS(root, pos);
	}
	inline void establishOrderDFS(SpeciesNode *node, int &pos) {
		node->begin = pos;
		if (node->child(0) != NULL) {
			establishOrderDFS(node->child(0), pos);
			pos++;
		}
		node->no = pos;
		if (node->child(1) != NULL) {
			pos++;
			establishOrderDFS(node->child(1), pos);
		}
		node->end = pos;
	}

};

// ------------------------------------------------------------------------------------------------------------------
// a rooted binary gene tree
class GeneTreeRooted : public TreeIO<Tree<GeneNodeRooted, NamedGeneNodeRooted> > {
public:
	
	bool rooting; // to specify if the tree is to be interpreted as rooted or unrooted by the reconciliation procedure
	
	inline bool isRooted() {
		return true;
	}

	// assigns each node an index
	void assignIndex() {
		for (int i = 0; i < nodes.size(); i++) nodes[i]->idx = i;
	}
};

// a quasi unrooted binary gene tree
// (rooted gene tree, but it allows rerooting)
class GeneTreeUnrooted : public TreeIO<Tree<GeneNodeUnrooted, NamedGeneNodeUnrooted> > {
public:
	inline bool isRooted() {
		return false;
	}

	// move the root into a different edge
	inline void reroot(GeneNodeUnrooted *u, GeneNodeUnrooted *v) {
		disconnectRoot();
		connectRoot(u, v);
		rerootDFS(root->child(0), root);
		rerootDFS(root->child(1), root);
	}
	inline void rerootDFS(GeneNodeUnrooted *&node, GeneNodeUnrooted *&parent) {
		if (node == NULL) return;
		if (node->parent() == parent) return;
		node->assignParent(parent);
		rerootDFS(node->child(0), node);
		rerootDFS(node->child(1), node);
	}
};

// ==================================================================================================================
// ==================================================================================================================
// ==================================================================================================================

// ------------------------------------------------------------------------------------------------------------------
// the container for a whole set of trees (1 species tree and several gene trees)
class TreeSet {
public:
	// data container for the species tree
	SpeciesTree *speciestree;

	// data container for the gene trees
	vector<GeneTreeRooted*> genetree_rooted;
	vector<GeneTreeUnrooted*> genetree_unrooted;

	TreeSet() {
		#ifdef DEBUG
		cout << "TreeSet created" << endl;
		#endif
		speciestree = NULL;
	}

	virtual ~TreeSet() {
		for(vector<GeneTreeRooted*>::iterator itr=genetree_rooted.begin(); itr!=genetree_rooted.end(); itr++) delete *itr;
		for(vector<GeneTreeUnrooted*>::iterator itr=genetree_unrooted.begin(); itr!=genetree_unrooted.end(); itr++) delete *itr;
		if (speciestree != NULL) delete speciestree;
		#ifdef DEBUG
		cout << "TreeSet destroyed" << endl;
		#endif
	}

	// ------------------------------------------------------------------------------------------------------
	// read all trees from the input
	void readTrees(Input &input) {
		// read the species tree
		speciestree = new SpeciesTree;
		speciestree->constraintcounter = 0;
		if (speciestree->stream2tree(input)) {
			msgout << "input species tree of " << speciestree->leafnodes.size() << " taxa" << endl;
		} else EXCEPTION("missing input for species tree " << input.getLastPos());
		speciestree->colorSpeciesTreeByConstraints();

		// read all gene trees
		unsigned int treecount = 0;
		for (;;) {
			string str;
			bool rooted = true;
			while (!(str=input.getComment()).empty()) {
				if (str == "[&U]") rooted = false;
				if (str == "[&R]") rooted = true;
			}

			if (rooted) {
				GeneTreeRooted *tree = new GeneTreeRooted;
				if (tree->stream2tree(input)) {
					for (int i=0, l=tree->nodes.size(); i<l; i++) tree->nodes[i]->treeid = treecount;
					genetree_rooted.push_back(tree);
					tree->rooting = true;
					msgout << "rooted input gene tree of " << tree->leafnodes.size() << " taxa" << endl;
				} else {
					delete tree;
					break;
				}
			} else {
				
				// print error message and quit if any of the gene trees are unrooted.
				EXCEPTION("Input contains unrooted gene tree(s). Please root all gene trees.");
			
			}
			treecount++;
		}
		if (genetree_rooted.size() + genetree_unrooted.size() < 1)
			EXCEPTION("missing input for gene trees " << input.getLastPos()); // ERROR no gene trees

		msgout << genetree_rooted.size() + genetree_unrooted.size() << " input gene trees total" << endl;
	}

	// ------------------------------------------------------------------------------------------------------
	// establish the initial primary mapping
	// mappings between the leaf nodes of gene trees and the species tree; according to their names
	void createLeafMapping() {
		// get the leaf nodes of the species tree
		// the leaf nodes of the species tree have to be a full set of all species and be unique
		vector<NamedSpeciesNode*> &s = speciestree->leafnodes;
		// sort these node by their name
		map<const string*, NamedSpeciesNode*> sm;
		for(vector<NamedSpeciesNode*>::iterator itr=s.begin(); itr!=s.end(); itr++) {
			NamedSpeciesNode &species = **itr;
			if (sm.find(&species.name->first) != sm.end()) // leafe node is not unique
				EXCEPTION('"' << species.name->first << "\" is not a unique species");
			sm[&species.name->first] = &species;
		}

		bool isPresent = false;
		size_t found;
//		string tempstring;
		// create the leaf node mapping for all rooted gene trees
		for(int i=0; i<genetree_rooted.size(); i++) {
			// get the leaf nodes of the gene tree
			vector<NamedGeneNodeRooted*> &g = genetree_rooted[i]->leafnodes;
			// create links for each leaf node to the species tree
			for(vector<NamedGeneNodeRooted*>::iterator itr=g.begin(); itr!=g.end(); itr++) {
				NamedGeneNodeRooted *genenode = *itr;
				
			// new code for handling the speciesName_geneName leaf name format for gene trees	
				found= genenode->name->first.find("_");
				if (found!=string::npos)
				{
					genenode->geneName = genenode->name->first.substr(0, found);

				}
				else
				{					
					genenode->geneName = genenode->name->first;
				}
					
//				map<const string*, NamedSpeciesNode*>::iterator smt = sm.find(&genenode.geneName);


				for(vector<NamedSpeciesNode*>::iterator itr=s.begin(); itr!=s.end(); itr++) 
				{
					isPresent = false;
					NamedSpeciesNode &species = **itr;
//					cout << species.name->first << endl;
					if (species.name->first == genenode->geneName)
					{
						isPresent= true;
						map<const string*, NamedSpeciesNode*>::iterator smt = sm.find(&species.name->first);
						genenode->setMapping(smt->second);
//						genenode->setMapping(species);
						break;
					}
					
				}
				if(isPresent == false)
					EXCEPTION("species \"" << genenode->geneName << "\" is not covered in the species tree");
				
				
				
				
				
//				map<const string*, NamedSpeciesNode*>::iterator smt = sm.find(&genenode->name->first);
//				if (smt == sm.end()) // leaf node of genetree is not in the species tree!
//				EXCEPTION("species \"" << genenode->geneName << "\" is not covered in the species tree");
				// establish mapping
//				genenode->setMapping(smt->second);
			}
		}
		
		// create the leaf node mapping for all unrooted gene trees
		for(int i=0; i<genetree_unrooted.size(); i++) {
			// get the leaf nodes of the gene tree
			vector<NamedGeneNodeUnrooted*> &g = genetree_unrooted[i]->leafnodes;
			// create links for each leaf node to the species tree
			for(vector<NamedGeneNodeUnrooted*>::iterator itr=g.begin(); itr!=g.end(); itr++) {
				NamedGeneNodeUnrooted *genenode = *itr;
				
				// new code for handling the speciesName_geneName leaf name format for gene trees
				found= genenode->name->first.find("_");
				if (found!=string::npos)
				{
					genenode->geneName = genenode->name->first.substr(0, found);

				}
				else
				{					
					genenode->geneName = genenode->name->first;
				}
					
//				map<const string*, NamedSpeciesNode*>::iterator smt = sm.find(&genenode.geneName);


				for(vector<NamedSpeciesNode*>::iterator itr=s.begin(); itr!=s.end(); itr++) 
				{
					isPresent = false;
					NamedSpeciesNode &species = **itr;
//					cout << species.name->first << endl;
					if (species.name->first == genenode->geneName)
					{
						isPresent= true;
						map<const string*, NamedSpeciesNode*>::iterator smt = sm.find(&species.name->first);
						genenode->setMapping(smt->second);
//						genenode->setMapping(species);
						break;
					}
					
				}
				if(isPresent == false)
					EXCEPTION("species \"" << genenode->geneName << "\" is not covered in the species tree");
				
				
/*				
				map<const string*, NamedSpeciesNode*>::iterator smt = sm.find(&genenode->name->first);
				if (smt == sm.end()) // leaf node of genetree is not in the species tree!
					EXCEPTION("species \"" << genenode->name->first << "\" is not covered in the species tree");
				// establish mapping
				genenode->setMapping(smt->second);
*/				
			}
		}

		
/*		
		// check for complete leaf checking (all species nodes have to map to at least 1 node in any gene tree)
		for(int i=0; i<genetree_rooted.size(); i++) {
			vector<NamedGeneNodeRooted*> &g = genetree_rooted[i]->leafnodes;
			for(vector<NamedGeneNodeRooted*>::iterator itr=g.begin(); itr!=g.end(); itr++) {
				NamedGeneNodeRooted *genenode = *itr;
				sm.erase(&genenode->name->first);
			}
		}
		for(int i=0; i<genetree_unrooted.size(); i++) {
			vector<NamedGeneNodeUnrooted*> &g = genetree_unrooted[i]->leafnodes;
			for(vector<NamedGeneNodeUnrooted*>::iterator itr=g.begin(); itr!=g.end(); itr++) {
				NamedGeneNodeUnrooted *genenode = *itr;
				sm.erase(&genenode->name->first);
			}
		}
		if (!sm.empty()) {
			for(map<const string*, NamedSpeciesNode*>::iterator itr=sm.begin(); itr!=sm.end(); itr++) {
				const string *name = itr->first;
				NamedSpeciesNode *&node = itr->second;
				WARNING("species \"" << *name << "\" is not covered by any gene tree - species ignored");
				speciestree->deleteLeafNode(node);
			}
		}
		
*/		
	}


	// ------------------------------------------------------------------------------------------------------
	// outputs the trees into a string stream
	friend ostream & operator << (ostream &os, TreeSet &t);
};

// outputs the trees into a string stream
ostream & operator << (ostream &os, TreeSet &t) {
	for(int i = 0; i < t.genetree_rooted.size(); i++ ) os << "gene tree" << i << ':' << endl << *t.genetree_rooted[i] << endl;
	for(int i = 0; i < t.genetree_unrooted.size(); i++ ) os << "gene tree" << i << ':' << endl << *t.genetree_unrooted[i] << endl;
	os << "species tree:" << endl << *t.speciestree << endl;
	return os;
}

#endif
