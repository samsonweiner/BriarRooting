/*
 * *   Copyright (C) 2017 Mukul S. Bansal (mukul.bansal@uconn.edu).
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




#ifndef DTL_ALGORITHM_H
#define DTL_ALGORITHM_H

#include <set>
#include <algorithm>
#include <iterator>


//#include <boost/random/mersenne_twister.hpp>
//#include <boost/random/uniform_int_distribution.hpp>
//#include <boost/cstdint.hpp>
#include <boost/multiprecision/cpp_int.hpp> 
#include <boost/multiprecision/random.hpp>

class DTLreconciliation :  public TreeSet {
	
protected:
	
	// define boost's random number generator
	boost::random::mt19937 gen; 
	

	
	vector< vector<int> > DParray; //create the dynamic programming table. (variable c(.,.) in paper).
//	vector< vector<int> > in;
	vector< vector<boost::multiprecision::cpp_int> > sol;  // to store the total number of optimal solutions for given g, s pair. Use for outputing solutions at random.

	
	int c_spec, c_dup, c_trans, minIn, minIn_second, minOut, minFinal, numDup, numTrans;
	int totalCost, totalDup, totalTrans, totalLoss;
	boost::multiprecision::cpp_int totalOptSol; // to store the total number of optimal solutions
	boost::multiprecision::cpp_int index;
	
	struct rootMapping{   // to store the best mappings from the gene tree root
		int bestmapindex;
		boost::multiprecision::cpp_int bestmapnum;
	};
	vector<rootMapping> bestMapIdx; /// // to store the best mappings from the gene tree root
	rootMapping tempBestMap;
	
	GeneNodeRooted *bestRoot; // to store the best rooting for unrooted gene trees.
	int bestScore; // to store the bestScore found so far when trying alternative roots for the gene tree.
	
	struct data{
		SpeciesNode *mapping0; // left child's mapping
		SpeciesNode *mapping1; // right child's mapping
		SpeciesNode *recipient;  // to store the recipient in case of transfer
		int event;		// 0 = speciation, 1 = duplication, 2 = transfer
		
	};
	
	vector< vector<data> > backTrack;    // table for back tracking
//	vector<data> optQueue;  // create the temporary structure that will store all equally optimal scenarios encountered during nested post order


	struct datarand{
		int event;		// 0 = speciation, 1 = duplication, 2 = transfer
		int casenumber; // to specify which specific case in Spec, Dup, or Transfer is being stored.
		boost::multiprecision::cpp_int numSolutions; // to store the number of solutions in that particular case
		
	};	
		vector<datarand> optQueue;  // create the temporary structure that will store all equally optimal scenarios encountered during nested post order

		
	struct dataTrack{  // datastructure for use with the queue used for tracking best in and out mappings.
		boost::multiprecision::cpp_int numSol;
		SpeciesNode *chosenmapping;
		};
		
		vector <dataTrack> trackingQueue; // Queue used for storing best mappings and number of solutins when assigning in and out mappings.
		
		
	// structure for storing the optimal mapping for the in, out, and in_second counters.
	struct datainout{
		SpeciesNode *optmapping; // using a struct for this simply for extensibility
	};
	
	vector< vector<datainout> > inTrack;   // tables for storing optimal mappings for the in, out, and in_second counters.
	vector< vector<datainout> > outTrack;
	 		
	
	vector<SpeciesNode*> bestInList; // vector to store all nodes of the relevant subtree in the getBestIn() procedure
	
	data chosentemp; // temporary instance of structure to store the finally chosen case and corresponding child mappings
	datarand temp; // create temporary instance of structure for storing intermediate optimal cases and event type.
	datainout tempinout; // create temporary instance of inout structure for storing intermediate optimal mapping.	
	dataTrack tempTrack; // temp instance of dataTrack for use in getBestInTrack and getBestOutTrack.
	
public:
	DTLreconciliation() {
		
		

	}


	
/*
	void SetGeneTreePOidx(GeneNodeRooted *currnode, int index) 
	{
		if (currnode == NULL) return;
		SetGeneTreePOidx(currnode->child(0), index);
		SetGeneTreePOidx(currnode->child(1), index);
		
		cout << currnode->idx << endl;
	}
	
	
	void SetSpeciesTreePOidx(SpeciesNode *currnode, int index) 
	{
		if (currnode == NULL) return;
		SetSpeciesTreePOidx(currnode->child(0), index);
		SetSpeciesTreePOidx(currnode->child(1), index);
		
		cout << currnode->idx << endl;

	}
	
*/	
	
	// set the leftmost decendant and rightmost descendants for the nodes of the gene tree.
	void leftRightDescendants(GeneNodeRooted *currnode) 
	{
		
		if (currnode == NULL) return;
		
		leftRightDescendants(currnode->child(0));
		leftRightDescendants(currnode->child(1));
		
		if (currnode->isLeaf() == true)
		{
			currnode->leftMost = currnode;
			currnode->rightMost = currnode;
			return;
		}
		else
		{
			currnode->leftMost = currnode->child(0)->leftMost;
			currnode->rightMost = currnode->child(1)->rightMost;
		}
		
	}
	
	
	
		// do a preOrder traversal of the species tree to label nodes with their depth (used for computing distances between nodes)
	void doNodePreOrder(SpeciesNode *node, int depthCount) {

			node->nodeDepth = depthCount;
			depthCount++;
		
		if (node->child(0) != NULL) {
			doNodePreOrder(node->child(0), depthCount);
		}
		if (node->child(1) != NULL) {
			doNodePreOrder(node->child(1), depthCount);
		}
	}

	
	
	
	
	
	// do a preOrder traversal of the species tree to label nodes with their depth in terms of branch lengths (for computing transfers)
	void doBranchPreOrder(SpeciesNode *node) {
		if (node == speciestree->root)
		{
			node->branchDepth = node->weight;
		}
		else
		{
			node->branchDepth = node->weight + node->parent()->branchDepth;
		}
		
		if (node->child(0) != NULL) {
			doBranchPreOrder(node->child(0));
		}
		if (node->child(1) != NULL) {
			doBranchPreOrder(node->child(1));
		}
	}

	
	
	// Perform nested post order traversal (main part of DTL-recociliation algo)
	void PostOrderGeneTree(GeneNodeRooted *currnode) 
	{
		if (currnode == NULL) return;
		PostOrderGeneTree(currnode->child(0));
		PostOrderGeneTree(currnode->child(1));
		
		// perform post order traversal of the species tree
		PostOrderSpeciesTree(currnode, speciestree->root);
		
		// update the out[][] counters as well as the outTrack table 
//		PreOrderSpeciesTree(currnode, speciestree->root);
		
		
	}
	
	
	void PostOrderSpeciesTree(GeneNodeRooted *currnode, SpeciesNode *CurrSpeciesNode) 
	{


		if (CurrSpeciesNode == NULL) return;
		PostOrderSpeciesTree(currnode, CurrSpeciesNode->child(0));
		PostOrderSpeciesTree(currnode, CurrSpeciesNode->child(1));

		
		// build up the DP table
		if(currnode->isLeaf() == false)
		{
			
			optQueue.clear(); // clear the queue that stores all optimal solutions
//			data temp; // create temporary instance of structure for storing optimal mappings and event type.
//			datainout tempinout; // create temporary instance of inout structure for storing optimal mapping.
			
			if(CurrSpeciesNode->isLeaf() == true)
			{
				
				
				// We always divide INT_MAX by 10 just to prevent overflow of the int value when we are adding two large integers.
				
				minFinal = INT_MAX / 10;  // to keep track of min value among speciation, duplication and transfer.
				
				//speciation 
				c_spec = INT_MAX / 10;
				
				
				// duplication 
				c_dup = DCOST + DParray[currnode->child(0)->idx][CurrSpeciesNode->idx] + DParray[currnode->child(1)->idx][CurrSpeciesNode->idx];
				
				if (c_dup < minFinal) 
				{
					minFinal = c_dup;
					optQueue.clear();

					temp.event = 1;
					temp.casenumber = 0;
					temp.numSolutions = sol[currnode->child(0)->idx][CurrSpeciesNode->idx] * sol[currnode->child(1)->idx][CurrSpeciesNode->idx];
					optQueue.push_back(temp);
				}
				else if ((c_dup == minFinal) || (minFinal >= INT_MAX / 10)) 
				{
					temp.event = 1;
					temp.casenumber = 0;
					temp.numSolutions = sol[currnode->child(0)->idx][CurrSpeciesNode->idx] * sol[currnode->child(1)->idx][CurrSpeciesNode->idx];
					optQueue.push_back(temp);
				}
				
				
				
				// transfer
				c_trans = INT_MAX / 10;
				int T1, T2;
				if(CurrSpeciesNode != speciestree->root)
				{
					// first transfer scenario
					T1 = getBestIn(currnode->child(0), CurrSpeciesNode) + getBestOut(currnode->child(1), CurrSpeciesNode);
					
					c_trans = TCOST + T1;
				
					if (c_trans < minFinal) 
					{
						minFinal = c_trans;
						optQueue.clear();
						temp.event = 2;
						temp.casenumber = 1;
						temp.numSolutions = getBestInCount(currnode->child(0), CurrSpeciesNode) * getBestOutCount(currnode->child(1), CurrSpeciesNode);
						optQueue.push_back(temp);
					}
					else if ((c_trans == minFinal) || (minFinal >= INT_MAX / 10)) 
					{
						temp.event = 2;
						temp.casenumber = 1;
						temp.numSolutions = getBestInCount(currnode->child(0), CurrSpeciesNode) * getBestOutCount(currnode->child(1), CurrSpeciesNode);
						optQueue.push_back(temp);
					}
					
					
					
					// second transfer scenario
					T2 = getBestIn(currnode->child(1), CurrSpeciesNode) + getBestOut(currnode->child(0), CurrSpeciesNode);
					
					c_trans = TCOST + T2;
					
					if (c_trans < minFinal) 
					{
						minFinal = c_trans;
						optQueue.clear();
						temp.event = 2;
						temp.casenumber = 2;
						temp.numSolutions = getBestInCount(currnode->child(1), CurrSpeciesNode) * getBestOutCount(currnode->child(0), CurrSpeciesNode);
						optQueue.push_back(temp);
					}
					else if ((c_trans == minFinal) || (minFinal >= INT_MAX / 10)) 
					{
						temp.event = 2;
						temp.casenumber = 2;
						temp.numSolutions = getBestInCount(currnode->child(1), CurrSpeciesNode) * getBestOutCount(currnode->child(0), CurrSpeciesNode);
						optQueue.push_back(temp);
					}

				}
				
				
				boost::multiprecision::cpp_int numSolSum = 0;
				// pick a random optimal scenario from the queue based on its probability.
				if (!optQueue.empty()) 
				{
					int x;
					
					// compute the total number of different solutions in all cases.
					for (x = 0; x < optQueue.size(); x++)
					{
						numSolSum = numSolSum + optQueue[x].numSolutions;
					}
					
					
//cout << "value 1 = " << numSolSum << " " << minFinal << " " << optQueue.size() << endl;
					// get the random number between 1 and numSolSum.
					if (numSolSum >= 1)
					{						
						index = getRandom(numSolSum);
//cout << "ended 1" << endl;
//					index = 1;
					
						// pick the case to be chosen based on the random number.
						boost::multiprecision::cpp_int tempsum = 0;
						for (x = 0; x < optQueue.size(); x++)
						{
							if ((index > tempsum) && (index <= (tempsum + optQueue[x].numSolutions)))
								break;
							tempsum = tempsum + optQueue[x].numSolutions;
						}
					
						// process the chosen case
						if (optQueue[x].event == 1)  // if duplication
						{
							chosentemp.mapping0 =  CurrSpeciesNode;
							chosentemp.mapping1 = CurrSpeciesNode;
							chosentemp.event = 1;
							chosentemp.recipient = NULL;
						}
					
						else if (optQueue[x].event == 2)  // if transfer
						{
							if (optQueue[x].casenumber == 1)
							{
								getBestInTrack(currnode->child(0), CurrSpeciesNode);
								chosentemp.mapping0 =  inTrack[currnode->child(0)->idx][CurrSpeciesNode->idx].optmapping;
								getBestOutTrack(currnode->child(1), CurrSpeciesNode);
								chosentemp.mapping1 = outTrack[currnode->child(1)->idx][CurrSpeciesNode->idx].optmapping;
								chosentemp.event = 2;
								chosentemp.recipient = chosentemp.mapping1;
							}
						
							else if (optQueue[x].casenumber == 2)
							{
								getBestInTrack(currnode->child(1), CurrSpeciesNode);
								chosentemp.mapping1 =  inTrack[currnode->child(1)->idx][CurrSpeciesNode->idx].optmapping;
								getBestOutTrack(currnode->child(0), CurrSpeciesNode);
								chosentemp.mapping0 = outTrack[currnode->child(0)->idx][CurrSpeciesNode->idx].optmapping;
								chosentemp.event = 2;
								chosentemp.recipient = chosentemp.mapping0;
							}
						
						
						}
					}
						
					// save the total number of solutions in the relevant array.															
					sol[currnode->idx][CurrSpeciesNode->idx] = numSolSum;
					
						
					
				}
				else
				{
					cout << "ERROR!! Empty Queue." << endl;
					exit(1);
				}	
				
				
				// set DParray entry 
				DParray[currnode->idx][CurrSpeciesNode->idx] = minFinal;
				
				if (numSolSum >=1)
				{
				// update the backtracking array
					backTrack[currnode->idx][CurrSpeciesNode->idx].mapping0 =  chosentemp.mapping0;
					backTrack[currnode->idx][CurrSpeciesNode->idx].mapping1 =  chosentemp.mapping1;
					backTrack[currnode->idx][CurrSpeciesNode->idx].recipient =  chosentemp.recipient;
					backTrack[currnode->idx][CurrSpeciesNode->idx].event =  chosentemp.event;
				}
				else 
				{
					backTrack[currnode->idx][CurrSpeciesNode->idx].mapping0 =  NULL;
					backTrack[currnode->idx][CurrSpeciesNode->idx].mapping1 =  NULL;
					backTrack[currnode->idx][CurrSpeciesNode->idx].recipient =  NULL;
					backTrack[currnode->idx][CurrSpeciesNode->idx].event =  -1;
				}
					
				
/*				
				
				// update in[][] entry as well as inTrack
				in[currnode->idx][CurrSpeciesNode->idx] = DParray[currnode->idx][CurrSpeciesNode->idx];
				inTrack[currnode->idx][CurrSpeciesNode->idx].optmapping = CurrSpeciesNode;
*/				

/*				
if((currnode->idx ==1) && (CurrSpeciesNode->idx == 1))
cout << endl << "BEGIN " << DParray[currnode->child(0)->idx][CurrSpeciesNode->idx] << "  " <<  out[currnode->child(1)->idx][CurrSpeciesNode->idx] << " " << DParray[currnode->idx][CurrSpeciesNode->idx];
*/				
			}
			else
			{
				
				minFinal = INT_MAX / 10;  // to keep track of min value among speciation, duplication and transfer.
				
				// speciation
				int S1 = getBestIn(currnode->child(0), CurrSpeciesNode->child(0)) + getBestIn(currnode->child(1), CurrSpeciesNode->child(1));
				
				c_spec = S1;
				
				if (c_spec < minFinal) 
				{
					minFinal = c_spec;
					optQueue.clear();
					temp.event = 0;
					temp.casenumber = 1;
					temp.numSolutions = getBestInCount(currnode->child(0), CurrSpeciesNode->child(0)) * getBestInCount(currnode->child(1), CurrSpeciesNode->child(1));
					optQueue.push_back(temp);
					
					
				}
				else if ((c_spec == minFinal) || (minFinal >= INT_MAX / 10)) 
				{
					temp.event = 0;
					temp.casenumber = 1;
					temp.numSolutions = getBestInCount(currnode->child(0), CurrSpeciesNode->child(0)) * getBestInCount(currnode->child(1), CurrSpeciesNode->child(1));
					optQueue.push_back(temp);
				}
				
				
				
				int S2 = getBestIn(currnode->child(0), CurrSpeciesNode->child(1)) + getBestIn(currnode->child(1), CurrSpeciesNode->child(0));

				c_spec = S2;
				
				if (c_spec < minFinal) 
				{
					minFinal = c_spec;
					optQueue.clear();
					temp.event = 0;
					temp.casenumber = 2;
					temp.numSolutions = getBestInCount(currnode->child(0), CurrSpeciesNode->child(1)) * getBestInCount(currnode->child(1), CurrSpeciesNode->child(0));
					optQueue.push_back(temp);
				}
				else if ((c_spec == minFinal) || (minFinal >= INT_MAX / 10)) 
				{
					temp.event = 0;
					temp.casenumber = 2;
					temp.numSolutions = getBestInCount(currnode->child(0), CurrSpeciesNode->child(1)) * getBestInCount(currnode->child(1), CurrSpeciesNode->child(0));
					optQueue.push_back(temp);
				}
				
				

				
				
				
				// duplication
//				minD = INT_MAX / 10;
				int D1 = LCOST + DParray[currnode->child(0)->idx][CurrSpeciesNode->idx] + getBestIn(currnode->child(1), CurrSpeciesNode->child(1));
					
				c_dup = DCOST + D1;
								
				if (c_dup < minFinal) 
				{
					minFinal = c_dup;
					optQueue.clear();
					temp.event = 1;
					temp.casenumber = 1;
					temp.numSolutions = sol[currnode->child(0)->idx][CurrSpeciesNode->idx] * getBestInCount(currnode->child(1), CurrSpeciesNode->child(1));
					optQueue.push_back(temp);
				}
				else if ((c_dup == minFinal) || (minFinal >= INT_MAX / 10)) 
				{
					temp.event = 1;
					temp.casenumber = 1;
					temp.numSolutions = sol[currnode->child(0)->idx][CurrSpeciesNode->idx] * getBestInCount(currnode->child(1), CurrSpeciesNode->child(1));
					optQueue.push_back(temp);
				}

				
				int D2 = LCOST + DParray[currnode->child(0)->idx][CurrSpeciesNode->idx] + getBestIn(currnode->child(1), CurrSpeciesNode->child(0));
				
				c_dup = DCOST + D2;
								
				if (c_dup < minFinal) 
				{
					minFinal = c_dup;
					optQueue.clear();
					temp.event = 1;
					temp.casenumber = 2;
					temp.numSolutions = sol[currnode->child(0)->idx][CurrSpeciesNode->idx] * getBestInCount(currnode->child(1), CurrSpeciesNode->child(0));
					optQueue.push_back(temp);
				}
				else if ((c_dup == minFinal) || (minFinal >= INT_MAX / 10)) 
				{
					temp.event = 1;
					temp.casenumber = 2;
					temp.numSolutions = sol[currnode->child(0)->idx][CurrSpeciesNode->idx] * getBestInCount(currnode->child(1), CurrSpeciesNode->child(0));
					optQueue.push_back(temp);
				}
				
				int D3 = LCOST + DParray[currnode->child(1)->idx][CurrSpeciesNode->idx] + getBestIn(currnode->child(0), CurrSpeciesNode->child(1));
				
				c_dup = DCOST + D3;
								
				if (c_dup < minFinal) 
				{
					minFinal = c_dup;
					optQueue.clear();
					temp.event = 1;
					temp.casenumber = 3;
					temp.numSolutions = sol[currnode->child(1)->idx][CurrSpeciesNode->idx] * getBestInCount(currnode->child(0), CurrSpeciesNode->child(1));
					optQueue.push_back(temp);
				}
				else if ((c_dup == minFinal) || (minFinal >= INT_MAX / 10)) 
				{
					temp.event = 1;
					temp.casenumber = 3;
					temp.numSolutions = sol[currnode->child(1)->idx][CurrSpeciesNode->idx] * getBestInCount(currnode->child(0), CurrSpeciesNode->child(1));
					optQueue.push_back(temp);
				}
				
				int D4 = LCOST + DParray[currnode->child(1)->idx][CurrSpeciesNode->idx] + getBestIn(currnode->child(0), CurrSpeciesNode->child(0));
				
				c_dup = DCOST + D4;
								
				if (c_dup < minFinal) 
				{
					minFinal = c_dup;
					optQueue.clear();
					temp.event = 1;
					temp.casenumber = 4;
					temp.numSolutions = sol[currnode->child(1)->idx][CurrSpeciesNode->idx] * getBestInCount(currnode->child(0), CurrSpeciesNode->child(0));
					optQueue.push_back(temp);
				}
				else if ((c_dup == minFinal) || (minFinal >= INT_MAX / 10)) 
				{
					temp.event = 1;
					temp.casenumber = 4;
					temp.numSolutions = sol[currnode->child(1)->idx][CurrSpeciesNode->idx] * getBestInCount(currnode->child(0), CurrSpeciesNode->child(0));
					optQueue.push_back(temp);
				}
				
				int D5 = DParray[currnode->child(0)->idx][CurrSpeciesNode->idx] + DParray[currnode->child(1)->idx][CurrSpeciesNode->idx];
				
				c_dup = DCOST + D5;
								
				if (c_dup < minFinal) 
				{
					minFinal = c_dup;
					optQueue.clear();
					temp.event = 1;
					temp.casenumber = 5;
					temp.numSolutions = sol[currnode->child(0)->idx][CurrSpeciesNode->idx] * sol[currnode->child(1)->idx][CurrSpeciesNode->idx];
					optQueue.push_back(temp);
				}
				else if ((c_dup == minFinal) || (minFinal >= INT_MAX / 10)) 
				{
					temp.event = 1;
					temp.casenumber = 5;
					temp.numSolutions = sol[currnode->child(0)->idx][CurrSpeciesNode->idx] * sol[currnode->child(1)->idx][CurrSpeciesNode->idx];
					optQueue.push_back(temp);
				}
				
				
				int D6 = (2 * LCOST) + getBestIn(currnode->child(0), CurrSpeciesNode->child(0)) + getBestIn(currnode->child(1), CurrSpeciesNode->child(1));
				
				c_dup = DCOST + D6;
								
				if (c_dup < minFinal) 
				{
					minFinal = c_dup;
					optQueue.clear();
					temp.event = 1;
					temp.casenumber = 6;
					temp.numSolutions = getBestInCount(currnode->child(0), CurrSpeciesNode->child(0)) * getBestInCount(currnode->child(1), CurrSpeciesNode->child(1));
					optQueue.push_back(temp);
				}
				else if ((c_dup == minFinal) || (minFinal >= INT_MAX / 10)) 
				{
					temp.event = 1;
					temp.casenumber = 6;
					temp.numSolutions = getBestInCount(currnode->child(0), CurrSpeciesNode->child(0)) * getBestInCount(currnode->child(1), CurrSpeciesNode->child(1));
					optQueue.push_back(temp);
				}
				
				int D7 = (2 * LCOST) + getBestIn(currnode->child(0), CurrSpeciesNode->child(1)) + getBestIn(currnode->child(1), CurrSpeciesNode->child(0));
				
				c_dup = DCOST + D7;
								
				if (c_dup < minFinal) 
				{
					minFinal = c_dup;
					optQueue.clear();
					temp.event = 1;
					temp.casenumber = 7;
					temp.numSolutions = getBestInCount(currnode->child(0), CurrSpeciesNode->child(1)) * getBestInCount(currnode->child(1), CurrSpeciesNode->child(0));
					optQueue.push_back(temp);
				}
				else if ((c_dup == minFinal) || (minFinal >= INT_MAX / 10)) 
				{
					temp.event = 1;
					temp.casenumber = 7;
					temp.numSolutions = getBestInCount(currnode->child(0), CurrSpeciesNode->child(1)) * getBestInCount(currnode->child(1), CurrSpeciesNode->child(0));
					optQueue.push_back(temp);
				}
				
				int D8 = (2 * LCOST) + getBestIn(currnode->child(0), CurrSpeciesNode->child(0)) + getBestIn(currnode->child(1), CurrSpeciesNode->child(0));
				
				c_dup = DCOST + D8;
								
				if (c_dup < minFinal) 
				{
					minFinal = c_dup;
					optQueue.clear();
					temp.event = 1;
					temp.casenumber = 8;
					temp.numSolutions = getBestInCount(currnode->child(0), CurrSpeciesNode->child(0)) * getBestInCount(currnode->child(1), CurrSpeciesNode->child(0));
					optQueue.push_back(temp);
				}
				else if ((c_dup == minFinal) || (minFinal >= INT_MAX / 10)) 
				{
					temp.event = 1;
					temp.casenumber = 8;
					temp.numSolutions = getBestInCount(currnode->child(0), CurrSpeciesNode->child(0)) * getBestInCount(currnode->child(1), CurrSpeciesNode->child(0));
					optQueue.push_back(temp);
				}
				
				int D9 = (2 * LCOST) + getBestIn(currnode->child(0), CurrSpeciesNode->child(1)) + getBestIn(currnode->child(1), CurrSpeciesNode->child(1));
				
				c_dup = DCOST + D9;
								
				if (c_dup < minFinal) 
				{
					minFinal = c_dup;
					optQueue.clear();
					temp.event = 1;
					temp.casenumber = 9;
					temp.numSolutions = getBestInCount(currnode->child(0), CurrSpeciesNode->child(1)) * getBestInCount(currnode->child(1), CurrSpeciesNode->child(1));
					optQueue.push_back(temp);
				}
				else if ((c_dup == minFinal) || (minFinal >= INT_MAX / 10)) 
				{
					temp.event = 1;
					temp.casenumber = 9;
					temp.numSolutions = getBestInCount(currnode->child(0), CurrSpeciesNode->child(1)) * getBestInCount(currnode->child(1), CurrSpeciesNode->child(1));
					optQueue.push_back(temp);
				}
				

				
				
				// transfer
				c_trans = INT_MAX / 10; 
				if(CurrSpeciesNode != speciestree->root)
				{
					int T1 = getBestIn(currnode->child(0), CurrSpeciesNode) + getBestOut(currnode->child(1), CurrSpeciesNode);
					
					
					c_trans = TCOST + T1;
					
					if (c_trans < minFinal) 
					{
						minFinal = c_trans;
						optQueue.clear();
						temp.event = 2;
						temp.casenumber = 1;
						temp.numSolutions = getBestInCount(currnode->child(0), CurrSpeciesNode) * getBestOutCount(currnode->child(1), CurrSpeciesNode);
						optQueue.push_back(temp);
					}
					else if ((c_trans == minFinal) || (minFinal >= INT_MAX / 10)) 
					{
						temp.event = 2;
						temp.casenumber = 1;
						temp.numSolutions = getBestInCount(currnode->child(0), CurrSpeciesNode) * getBestOutCount(currnode->child(1), CurrSpeciesNode);
						optQueue.push_back(temp);
					}
					
					
					int T2 = getBestIn(currnode->child(1), CurrSpeciesNode) + getBestOut(currnode->child(0), CurrSpeciesNode);
					
					c_trans = TCOST + T2;
					
					if (c_trans < minFinal) 
					{
						minFinal = c_trans;
						optQueue.clear();
						temp.event = 2;
						temp.casenumber = 2;
						temp.numSolutions = getBestInCount(currnode->child(1), CurrSpeciesNode) * getBestOutCount(currnode->child(0), CurrSpeciesNode);
						optQueue.push_back(temp);
					}
					else if ((c_trans == minFinal) || (minFinal >= INT_MAX / 10)) 
					{
						temp.event = 2;
						temp.casenumber = 2;
						temp.numSolutions = getBestInCount(currnode->child(1), CurrSpeciesNode) * getBestOutCount(currnode->child(0), CurrSpeciesNode);
						optQueue.push_back(temp);
					}
					
				}
					
				
				boost::multiprecision::cpp_int numSolSum = 0;
// pick a random optimal scenario from the queue based on its probability.
				if (!optQueue.empty()) 
				{
					int x;
					
					
					// compute the total number of different solutions in all cases.
					for (x = 0; x < optQueue.size(); x++)
					{
						numSolSum = numSolSum + optQueue[x].numSolutions;
					}
					
				
//cout << "value 2 = " << numSolSum << endl;					
					// get the random number between 1 and numSolSum.
					if (numSolSum >=1)
					{						
						index = getRandom(numSolSum);
//						index = 1;
//cout << "end 2" << endl;
					
						// pick the case to be chosen based on the random number.
						boost::multiprecision::cpp_int tempsum = 0;
						for (x = 0; x < optQueue.size(); x++)
						{
							if ((index > tempsum) && (index <= (tempsum + optQueue[x].numSolutions)))
								break;
							tempsum = tempsum + optQueue[x].numSolutions;
						}
					
						// process the chosen case
					
						if (optQueue[x].event == 0)  // if speciation
						{
							if (optQueue[x].casenumber == 1)
							{
								getBestInTrack(currnode->child(0), CurrSpeciesNode->child(0));
								chosentemp.mapping0 = inTrack[currnode->child(0)->idx][CurrSpeciesNode->child(0)->idx].optmapping;
								getBestInTrack(currnode->child(1), CurrSpeciesNode->child(1));
								chosentemp.mapping1 = inTrack[currnode->child(1)->idx][CurrSpeciesNode->child(1)->idx].optmapping;
								chosentemp.recipient = NULL;
								chosentemp.event = 0;
							}
						
							else if (optQueue[x].casenumber == 2)
							{
								getBestInTrack(currnode->child(0), CurrSpeciesNode->child(1));
								chosentemp.mapping0 = inTrack[currnode->child(0)->idx][CurrSpeciesNode->child(1)->idx].optmapping;
								getBestInTrack(currnode->child(1), CurrSpeciesNode->child(0));
								chosentemp.mapping1 = inTrack[currnode->child(1)->idx][CurrSpeciesNode->child(0)->idx].optmapping;
								chosentemp.recipient = NULL;
								chosentemp.event = 0;
							}
						
						}
					
				
				
					
						else if (optQueue[x].event == 1)  // if duplication
						{
							if (optQueue[x].casenumber == 1)
							{
								chosentemp.mapping0 = CurrSpeciesNode;
								getBestInTrack(currnode->child(1), CurrSpeciesNode->child(1));
								chosentemp.mapping1 = inTrack[currnode->child(1)->idx][CurrSpeciesNode->child(1)->idx].optmapping;
								chosentemp.event = 1;
								chosentemp.recipient = NULL;
							}
						
							else if (optQueue[x].casenumber == 2)
							{
								chosentemp.mapping0 = CurrSpeciesNode;
								getBestInTrack(currnode->child(1), CurrSpeciesNode->child(0));
								chosentemp.mapping1 = inTrack[currnode->child(1)->idx][CurrSpeciesNode->child(0)->idx].optmapping;
								chosentemp.recipient = NULL;
								chosentemp.event = 1;
							}
							
							else if (optQueue[x].casenumber == 3)
							{
								getBestInTrack(currnode->child(0), CurrSpeciesNode->child(1));
								chosentemp.mapping0 = inTrack[currnode->child(0)->idx][CurrSpeciesNode->child(1)->idx].optmapping;
								chosentemp.mapping1 = CurrSpeciesNode;
								chosentemp.recipient = NULL;
								chosentemp.event = 1;
							}	
							
							else if (optQueue[x].casenumber == 4)
							{
								getBestInTrack(currnode->child(0), CurrSpeciesNode->child(0));
								chosentemp.mapping0 = inTrack[currnode->child(0)->idx][CurrSpeciesNode->child(0)->idx].optmapping;
								chosentemp.mapping1 = CurrSpeciesNode;
								chosentemp.recipient = NULL;
								chosentemp.event = 1;
							}	
							
							else if (optQueue[x].casenumber == 5)
							{
								chosentemp.mapping0 = CurrSpeciesNode;
								chosentemp.mapping1 = CurrSpeciesNode;
								chosentemp.recipient = NULL;
								chosentemp.event = 1;
							}
							
							else if (optQueue[x].casenumber == 6)
							{
								getBestInTrack(currnode->child(0), CurrSpeciesNode->child(0));
								chosentemp.mapping0 = inTrack[currnode->child(0)->idx][CurrSpeciesNode->child(0)->idx].optmapping;
								getBestInTrack(currnode->child(1), CurrSpeciesNode->child(1));
								chosentemp.mapping1 = inTrack[currnode->child(1)->idx][CurrSpeciesNode->child(1)->idx].optmapping;
								chosentemp.recipient = NULL;
								chosentemp.event = 1;
							}
							
							else if (optQueue[x].casenumber == 7)
							{
								getBestInTrack(currnode->child(0), CurrSpeciesNode->child(1));
								chosentemp.mapping0 = inTrack[currnode->child(0)->idx][CurrSpeciesNode->child(1)->idx].optmapping;
								getBestInTrack(currnode->child(1), CurrSpeciesNode->child(0));
								chosentemp.mapping1 = inTrack[currnode->child(1)->idx][CurrSpeciesNode->child(0)->idx].optmapping;
								chosentemp.recipient = NULL;
								chosentemp.event = 1;
							}
							
							else if (optQueue[x].casenumber == 8)
							{
								getBestInTrack(currnode->child(0), CurrSpeciesNode->child(0));
								chosentemp.mapping0 = inTrack[currnode->child(0)->idx][CurrSpeciesNode->child(0)->idx].optmapping;
								getBestInTrack(currnode->child(1), CurrSpeciesNode->child(0));
								chosentemp.mapping1 = inTrack[currnode->child(1)->idx][CurrSpeciesNode->child(0)->idx].optmapping;
								chosentemp.recipient = NULL;
								chosentemp.event = 1;
							}	
							
							else if (optQueue[x].casenumber == 9)
							{
								getBestInTrack(currnode->child(0), CurrSpeciesNode->child(1));
								chosentemp.mapping0 = inTrack[currnode->child(0)->idx][CurrSpeciesNode->child(1)->idx].optmapping;
								getBestInTrack(currnode->child(1), CurrSpeciesNode->child(1));
								chosentemp.mapping1 = inTrack[currnode->child(1)->idx][CurrSpeciesNode->child(1)->idx].optmapping;
								chosentemp.recipient = NULL;
								chosentemp.event = 1;
							}						
							
							
						}
						
						else if (optQueue[x].event == 2)  // if transfer
						{
							if (optQueue[x].casenumber == 1)
							{
								getBestInTrack(currnode->child(0), CurrSpeciesNode);
								chosentemp.mapping0 =  inTrack[currnode->child(0)->idx][CurrSpeciesNode->idx].optmapping;
								getBestOutTrack(currnode->child(1), CurrSpeciesNode);
								chosentemp.mapping1 = outTrack[currnode->child(1)->idx][CurrSpeciesNode->idx].optmapping;
								chosentemp.event = 2;
								chosentemp.recipient = chosentemp.mapping1;
							}
							
							else if (optQueue[x].casenumber == 2)
							{
								getBestInTrack(currnode->child(1), CurrSpeciesNode);
								chosentemp.mapping1 =  inTrack[currnode->child(1)->idx][CurrSpeciesNode->idx].optmapping;
								getBestOutTrack(currnode->child(0), CurrSpeciesNode);
								chosentemp.mapping0 = outTrack[currnode->child(0)->idx][CurrSpeciesNode->idx].optmapping;
								chosentemp.event = 2;
								chosentemp.recipient = chosentemp.mapping0;
							}
							
						}
					
					}
					// save the total number of solutions in the relevant array.															
					sol[currnode->idx][CurrSpeciesNode->idx] = numSolSum;																
							
				}
				else
				{
					cout << "ERROR!! Empty Queue." << endl;
					exit(1);
				}	
				
				
				// set DParray entry
				DParray[currnode->idx][CurrSpeciesNode->idx] = minFinal;
				
				// update the backtracking array
				if(numSolSum >=1)
				{
					backTrack[currnode->idx][CurrSpeciesNode->idx].mapping0 =  chosentemp.mapping0;
					backTrack[currnode->idx][CurrSpeciesNode->idx].mapping1 =  chosentemp.mapping1;
					backTrack[currnode->idx][CurrSpeciesNode->idx].recipient =  chosentemp.recipient;
					backTrack[currnode->idx][CurrSpeciesNode->idx].event =  chosentemp.event;
				}
				else
				{
					backTrack[currnode->idx][CurrSpeciesNode->idx].mapping0 =  NULL;
					backTrack[currnode->idx][CurrSpeciesNode->idx].mapping1 =  NULL;
					backTrack[currnode->idx][CurrSpeciesNode->idx].recipient =  NULL;
					backTrack[currnode->idx][CurrSpeciesNode->idx].event =  -1;
				}
					
				
/*				// update in[][] entry and the inTrack table
				
				minIn = DParray[currnode->idx][CurrSpeciesNode->idx];
				tempinout.optmapping = CurrSpeciesNode; 
				if ( (LCOST + in[currnode->idx][CurrSpeciesNode->child(0)->idx]) < minIn)
				{
					minIn = (LCOST + in[currnode->idx][CurrSpeciesNode->child(0)->idx]);
					tempinout.optmapping = inTrack[currnode->idx][CurrSpeciesNode->child(0)->idx].optmapping;
				}
				if ( (LCOST + in[currnode->idx][CurrSpeciesNode->child(1)->idx]) < minIn)
				{
					minIn = (LCOST + in[currnode->idx][CurrSpeciesNode->child(1)->idx]);
					tempinout.optmapping = inTrack[currnode->idx][CurrSpeciesNode->child(1)->idx].optmapping;
				}
				
				in[currnode->idx][CurrSpeciesNode->idx] = minIn;
				inTrack[currnode->idx][CurrSpeciesNode->idx].optmapping = tempinout.optmapping;
 */
			}
				
				
			
		}
	}
	
	
	// get random number using boosts random number generator
	boost::multiprecision::cpp_int getRandom(boost::multiprecision::cpp_int limit)
	{
		boost::random::uniform_int_distribution<boost::multiprecision::cpp_int> distrib(1, limit);
		return distrib(gen);
	}
	
	
	
	// function for computing the best "in" value by brute force.  
	int getBestIn(GeneNodeRooted *currnode, SpeciesNode *node)
	{
		minIn = INT_MAX / 10;
		
		bool valid;
		int num_nodes = speciestree->nodes.size();
		int distance;
		SpeciesNode* temp_node;
		SpeciesNode* lcanode;
		
		// get the list of nodes in the subtree rooted at the current species node.
		bestInList.clear();
		bestInPreorder(node);
		
		int listsize = bestInList.size();
		
		for (int i = 0; i< listsize; i++)
		{
			
			distance = bestInList[i]->nodeDepth - node->nodeDepth;
			
			if (((LCOST * distance) + DParray[currnode->idx][bestInList[i]->idx]) < minIn)
			{
				minIn = (LCOST * distance) + DParray[currnode->idx][bestInList[i]->idx];
//				inTrack[currnode->idx][node->idx].optmapping = bestInList[i];
			}
		}
		
		return minIn;
		
	}
	
	
	boost::multiprecision::cpp_int getBestInCount(GeneNodeRooted *currnode, SpeciesNode *node)
	{
		minIn = INT_MAX / 10;
		
		bool valid;
		int num_nodes = speciestree->nodes.size();
		int distance;
		SpeciesNode* temp_node;
		SpeciesNode* lcanode;
		boost::multiprecision::cpp_int count = 0;
		
		// get the list of nodes in the subtree rooted at the current species node.
		bestInList.clear();
		bestInPreorder(node);
		
		int listsize = bestInList.size();
		
		for (int i = 0; i< listsize; i++)
		{
			
			distance = bestInList[i]->nodeDepth - node->nodeDepth;
			
			if (((LCOST * distance) + DParray[currnode->idx][bestInList[i]->idx]) < minIn)
			{
				minIn = (LCOST * distance) + DParray[currnode->idx][bestInList[i]->idx];
				count = sol[currnode->idx][bestInList[i]->idx];
//				inTrack[currnode->idx][node->idx].optmapping = bestInList[i];
			}
			else if (((LCOST * distance) + DParray[currnode->idx][bestInList[i]->idx]) == minIn)
			{
				count = count + sol[currnode->idx][bestInList[i]->idx];
			}
			
		}
		
		return count;
		
	}	
	
	// to randomly return a best 'in' mapping.
	void getBestInTrack(GeneNodeRooted *currnode, SpeciesNode *node)
	{
		minIn = INT_MAX / 10;
		
		bool valid;
		int num_nodes = speciestree->nodes.size();
		int distance;
		SpeciesNode* temp_node;
		SpeciesNode* lcanode;
		boost::multiprecision::cpp_int count = 0;
		boost::multiprecision::cpp_int Index = 0; 
		// clear the queue;
		trackingQueue.clear();
		
		// get the list of nodes in the subtree rooted at the current species node.
		bestInList.clear();
		bestInPreorder(node);
		
		int listsize = bestInList.size();
		
		for (int i = 0; i< listsize; i++)
		{
			
			distance = bestInList[i]->nodeDepth - node->nodeDepth;
			
			if (((LCOST * distance) + DParray[currnode->idx][bestInList[i]->idx]) < minIn)
			{
				minIn = (LCOST * distance) + DParray[currnode->idx][bestInList[i]->idx];
				count = sol[currnode->idx][bestInList[i]->idx];
				trackingQueue.clear();
				tempTrack.numSol = sol[currnode->idx][bestInList[i]->idx];
				tempTrack.chosenmapping = bestInList[i];
				trackingQueue.push_back(tempTrack);
								
//				inTrack[currnode->idx][node->idx].optmapping = bestInList[i];
			}
			else if (((LCOST * distance) + DParray[currnode->idx][bestInList[i]->idx]) == minIn)
			{
				count = count + sol[currnode->idx][bestInList[i]->idx];
				tempTrack.numSol = sol[currnode->idx][bestInList[i]->idx];
				tempTrack.chosenmapping = bestInList[i];
				trackingQueue.push_back(tempTrack);
			}
			
		}
		
		
//cout << "value 3 = " << count << endl;		
		// get the random number between 1 and count.
		Index = getRandom(count);
//		Index =  1;				
					
//cout << "end 3" << endl;	
		// pick the case to be chosen based on the random number.
		boost::multiprecision::cpp_int tempsum = 0;
		int x = 0;
		for (x = 0; x < trackingQueue.size(); x++)
		{
			if ((Index > tempsum) && (Index <= (tempsum + trackingQueue[x].numSol)))
				break;
			tempsum = tempsum + trackingQueue[x].numSol;
		}
		
		
		inTrack[currnode->idx][node->idx].optmapping = trackingQueue[x].chosenmapping;
		
		
	}
	
	
	
	// function for returning the list of all species tree nodes that are in the subtree rooted at the input node.
	void bestInPreorder(SpeciesNode *node)
	{
		bestInList.push_back(node);
		
		if (node->child(0) != NULL)
			bestInPreorder(node->child(0));
		if (node->child(1) != NULL)
			bestInPreorder(node->child(1));
		
		return;
	}
	
	
	
	
	// This function computes the best receiver for a transfer (by brute force for distance dependant transfer costs). 
	int getBestOut(GeneNodeRooted *currnode, SpeciesNode *node)
	{
		minOut = INT_MAX / 10;
			
		bool valid;
		int transferAdd;   // additional cost of transfer based on distance
		int num_nodes = speciestree->nodes.size();
		int distance;
		SpeciesNode* temp_node;
		SpeciesNode* lcanode;
		
		// iterate over all non-root nodes
		for (int i = 1; i< num_nodes; i++)
		{
			
//if (speciestree->nodes[i]->child(0) == NULL)
//cout << speciestree->nodes[i]->branchDepth << endl;
			
			if (speciestree->nodes[i] == node)
				continue;
			
			
			valid = true;
			
			// disallow transfers that go from descendant to ancestor or vice versa. All other transfers are OK.
			if ((node->no >= speciestree->nodes[i]->begin) && (node->no <= speciestree->nodes[i]->end))
				valid = false;
			if ((speciestree->nodes[i]->no >= node->begin) && (speciestree->nodes[i]->no <= node->end))
				valid = false;
			
			if (valid == true)
			{
				
				temp_node = node;
				while(temp_node->parent() != NULL)
				{
					temp_node = temp_node->parent();
					if((speciestree->nodes[i]->no >= temp_node->begin) && (speciestree->nodes[i]->no <= temp_node->end))
						break;
				}
				
				lcanode = temp_node;
				
				
				distance = (node->nodeDepth - lcanode->nodeDepth) + (speciestree->nodes[i]->nodeDepth - lcanode->nodeDepth);
				
				
				
				
				
//				if((node->no >= speciestree->nodes[i]->begin) && (node->no <= speciestree->nodes[i]->end))
//					cout << "ERROR!!!" << endl;
//				if((speciestree->nodes[i]->no >= node->begin) && (speciestree->nodes[i]->no <= node->end))
//					cout << "ERROR!!!" << endl;
				
				
				if (TYPE == 0)   // if not using variable transfer costs
				{
					if (DParray[currnode->idx][speciestree->nodes[i]->idx] < minOut)
					{
						minOut = DParray[currnode->idx][speciestree->nodes[i]->idx];
//						outTrack[currnode->idx][node->idx].optmapping = speciestree->nodes[i];
					}
				}
				if (TYPE == 1)  // if using the simple transfer cost
				{
					
					if (distance >= THRESHOLD)
						transferAdd = ADD;
					else
						transferAdd = 0;
					
					if ((DParray[currnode->idx][speciestree->nodes[i]->idx] + transferAdd) < minOut)
					{
						minOut = (DParray[currnode->idx][speciestree->nodes[i]->idx] + transferAdd);
//						outTrack[currnode->idx][node->idx].optmapping = speciestree->nodes[i];
					}
					
				}
				
				
				if (TYPE == 2)   // if using the more general transfer cost
				{
					
					if (distance >= THRESHOLD)
						transferAdd = (distance/THRESHOLD) * ADD;
					else
						transferAdd = 0;
					
					if ((DParray[currnode->idx][speciestree->nodes[i]->idx] + transferAdd) < minOut)
					{
						minOut = (DParray[currnode->idx][speciestree->nodes[i]->idx] + transferAdd);
//						outTrack[currnode->idx][node->idx].optmapping = speciestree->nodes[i];
					}
				}
				
				

				
				
			}

		}
 		
		return minOut;
	}
	
	
	
	
	// This function computes the best receiver for a transfer (by brute force). 
	boost::multiprecision::cpp_int getBestOutCount(GeneNodeRooted *currnode, SpeciesNode *node)
	{
		minOut = INT_MAX / 10;
			
		bool valid;
		int transferAdd;   // additional cost of transfer based on distance
		int num_nodes = speciestree->nodes.size();
		int distance;
		SpeciesNode* temp_node;
		SpeciesNode* lcanode;
		
		boost::multiprecision::cpp_int count  = 0;
		
		// iterate over all non-root nodes
		for (int i = 1; i< num_nodes; i++)
		{
			
//if (speciestree->nodes[i]->child(0) == NULL)
//cout << speciestree->nodes[i]->branchDepth << endl;
			
			if (speciestree->nodes[i] == node)
				continue;
			
			
			valid = true;
			
			// disallow transfers that go from descendant to ancestor or vice versa. All other transfers are OK.
			if ((node->no >= speciestree->nodes[i]->begin) && (node->no <= speciestree->nodes[i]->end))
				valid = false;
			if ((speciestree->nodes[i]->no >= node->begin) && (speciestree->nodes[i]->no <= node->end))
				valid = false;
			
			if (valid == true)
			{
				
				temp_node = node;
				while(temp_node->parent() != NULL)
				{
					temp_node = temp_node->parent();
					if((speciestree->nodes[i]->no >= temp_node->begin) && (speciestree->nodes[i]->no <= temp_node->end))
						break;
				}
				
				lcanode = temp_node;
				
				
				distance = (node->nodeDepth - lcanode->nodeDepth) + (speciestree->nodes[i]->nodeDepth - lcanode->nodeDepth);
				
				
				
				
				
//				if((node->no >= speciestree->nodes[i]->begin) && (node->no <= speciestree->nodes[i]->end))
//					cout << "ERROR!!!" << endl;
//				if((speciestree->nodes[i]->no >= node->begin) && (speciestree->nodes[i]->no <= node->end))
//					cout << "ERROR!!!" << endl;
				
				
				if (TYPE == 0)   // if not using variable transfer costs
				{
					if (DParray[currnode->idx][speciestree->nodes[i]->idx] < minOut)
					{
						minOut = DParray[currnode->idx][speciestree->nodes[i]->idx];
//						outTrack[currnode->idx][node->idx].optmapping = speciestree->nodes[i];
						count  = sol[currnode->idx][speciestree->nodes[i]->idx];
					}
					else if (DParray[currnode->idx][speciestree->nodes[i]->idx] == minOut)
					{
						count  =  count + sol[currnode->idx][speciestree->nodes[i]->idx];	
					}
				}
				
				if (TYPE == 1)  // if using the simple transfer cost
				{
					
					if (distance >= THRESHOLD)
						transferAdd = ADD;
					else
						transferAdd = 0;
					
					if ((DParray[currnode->idx][speciestree->nodes[i]->idx] + transferAdd) < minOut)
					{
						minOut = (DParray[currnode->idx][speciestree->nodes[i]->idx] + transferAdd);
//						outTrack[currnode->idx][node->idx].optmapping = speciestree->nodes[i];
						count  = sol[currnode->idx][speciestree->nodes[i]->idx];
					}
					else if ((DParray[currnode->idx][speciestree->nodes[i]->idx] + transferAdd) == minOut)
					{

						count  = count + sol[currnode->idx][speciestree->nodes[i]->idx];
					}
					
				}
				
				
				if (TYPE == 2)   // if using the more general transfer cost
				{
					
					if (distance >= THRESHOLD)
						transferAdd = (distance/THRESHOLD) * ADD;
					else
						transferAdd = 0;
					
					if ((DParray[currnode->idx][speciestree->nodes[i]->idx] + transferAdd) < minOut)
					{
						minOut = (DParray[currnode->idx][speciestree->nodes[i]->idx] + transferAdd);
//						outTrack[currnode->idx][node->idx].optmapping = speciestree->nodes[i];
						count  = sol[currnode->idx][speciestree->nodes[i]->idx];
					}
					if ((DParray[currnode->idx][speciestree->nodes[i]->idx] + transferAdd) == minOut)
					{
						count  =  count + sol[currnode->idx][speciestree->nodes[i]->idx];

					}
				}

				
			}

		}
 		
		return count;
	}	
	
	
	
	// randomly choose a best 'out' mapping
	void getBestOutTrack(GeneNodeRooted *currnode, SpeciesNode *node)
	{
		minOut = INT_MAX / 10;
			
		bool valid;
		int transferAdd;   // additional cost of transfer based on distance
		int num_nodes = speciestree->nodes.size();
		int distance;
		SpeciesNode* temp_node;
		SpeciesNode* lcanode;
		
		boost::multiprecision::cpp_int count  = 0;
		
		boost::multiprecision::cpp_int Index = 0; 
		// clear the queue;
		trackingQueue.clear();
		
		// iterate over all non-root nodes
		for (int i = 1; i< num_nodes; i++)
		{
			
//if (speciestree->nodes[i]->child(0) == NULL)
//cout << speciestree->nodes[i]->branchDepth << endl;
			
			if (speciestree->nodes[i] == node)
				continue;
			
			
			valid = true;
			
			// disallow transfers that go from descendant to ancestor or vice versa. All other transfers are OK.
			if ((node->no >= speciestree->nodes[i]->begin) && (node->no <= speciestree->nodes[i]->end))
				valid = false;
			if ((speciestree->nodes[i]->no >= node->begin) && (speciestree->nodes[i]->no <= node->end))
				valid = false;
			
			if (valid == true)
			{
				
				temp_node = node;
				while(temp_node->parent() != NULL)
				{
					temp_node = temp_node->parent();
					if((speciestree->nodes[i]->no >= temp_node->begin) && (speciestree->nodes[i]->no <= temp_node->end))
						break;
				}
				
				lcanode = temp_node;
				
				
				distance = (node->nodeDepth - lcanode->nodeDepth) + (speciestree->nodes[i]->nodeDepth - lcanode->nodeDepth);
				
				
				
				
				
//				if((node->no >= speciestree->nodes[i]->begin) && (node->no <= speciestree->nodes[i]->end))
//					cout << "ERROR!!!" << endl;
//				if((speciestree->nodes[i]->no >= node->begin) && (speciestree->nodes[i]->no <= node->end))
//					cout << "ERROR!!!" << endl;
				
				
				if (TYPE == 0)   // if not using variable transfer costs
				{
					if (DParray[currnode->idx][speciestree->nodes[i]->idx] < minOut)
					{
						minOut = DParray[currnode->idx][speciestree->nodes[i]->idx];
//						outTrack[currnode->idx][node->idx].optmapping = speciestree->nodes[i];
						count  = sol[currnode->idx][speciestree->nodes[i]->idx];
						trackingQueue.clear();
						tempTrack.numSol = sol[currnode->idx][speciestree->nodes[i]->idx];
						tempTrack.chosenmapping = speciestree->nodes[i];
						trackingQueue.push_back(tempTrack);
						
					}
					else if (DParray[currnode->idx][speciestree->nodes[i]->idx] == minOut)
					{
						count  = count + sol[currnode->idx][speciestree->nodes[i]->idx];
						tempTrack.numSol = sol[currnode->idx][speciestree->nodes[i]->idx];
						tempTrack.chosenmapping = speciestree->nodes[i];
						trackingQueue.push_back(tempTrack);	
					}
				}
				
				if (TYPE == 1)  // if using the simple transfer cost
				{
					
					if (distance >= THRESHOLD)
						transferAdd = ADD;
					else
						transferAdd = 0;
					
					if ((DParray[currnode->idx][speciestree->nodes[i]->idx] + transferAdd) < minOut)
					{
						minOut = (DParray[currnode->idx][speciestree->nodes[i]->idx] + transferAdd);
						count  = sol[currnode->idx][speciestree->nodes[i]->idx];
						trackingQueue.clear();
						tempTrack.numSol = sol[currnode->idx][speciestree->nodes[i]->idx];
						tempTrack.chosenmapping = speciestree->nodes[i];
						trackingQueue.push_back(tempTrack);
					}
					else if ((DParray[currnode->idx][speciestree->nodes[i]->idx] + transferAdd) == minOut)
					{

						count  = count + sol[currnode->idx][speciestree->nodes[i]->idx];
						tempTrack.numSol = sol[currnode->idx][speciestree->nodes[i]->idx];
						tempTrack.chosenmapping = speciestree->nodes[i];
						trackingQueue.push_back(tempTrack);
					}
					
				}
				
				
				if (TYPE == 2)   // if using the more general transfer cost
				{
					
					if (distance >= THRESHOLD)
						transferAdd = (distance/THRESHOLD) * ADD;
					else
						transferAdd = 0;
					
					if ((DParray[currnode->idx][speciestree->nodes[i]->idx] + transferAdd) < minOut)
					{
						minOut = (DParray[currnode->idx][speciestree->nodes[i]->idx] + transferAdd);
						count  = sol[currnode->idx][speciestree->nodes[i]->idx];
						trackingQueue.clear();
						tempTrack.numSol = sol[currnode->idx][speciestree->nodes[i]->idx];
						tempTrack.chosenmapping = speciestree->nodes[i];
						trackingQueue.push_back(tempTrack);
					}
					if ((DParray[currnode->idx][speciestree->nodes[i]->idx] + transferAdd) == minOut)
					{
						count  = count + sol[currnode->idx][speciestree->nodes[i]->idx];
						tempTrack.numSol = sol[currnode->idx][speciestree->nodes[i]->idx];
						tempTrack.chosenmapping = speciestree->nodes[i];
						trackingQueue.push_back(tempTrack);

					}
				}

				
			}

		}
 		
//cout << "value 4 = " << count << endl;
		// get the random number between 1 and count.
		Index = getRandom(count);
//		Index = 1;	
//cout << "end 4" << endl;
		// pick the case to be chosen based on the random number.
		boost::multiprecision::cpp_int tempsum = 0;
		int x = 0;
		for (x = 0; x < trackingQueue.size(); x++)
		{
			if ((Index > tempsum) && (Index <= (tempsum + trackingQueue[x].numSol)))
				break;
			tempsum = tempsum + trackingQueue[x].numSol;
		}
		
		
		outTrack[currnode->idx][node->idx].optmapping = trackingQueue[x].chosenmapping;
	}
	
	
	
	
	
	// output the reconciliation cost and call function for outputting actual reconciliation
	void OutputReconciliationCost(ostream &os, GeneNodeRooted *currnode, bool summary)
	{
		
		int minReconCost = INT_MAX / 10;
		int numSpecies_nodes = speciestree->nodes.size();
		bestMapIdx.clear();
		totalOptSol = 0;
		
		
		for (int i = 0; i < numSpecies_nodes; i++)
		{
			if (DParray[currnode->idx][i] < minReconCost)
			{
				minReconCost = DParray[currnode->idx][i];
				bestMapIdx.clear();
				totalOptSol = sol[currnode->idx][i];
				tempBestMap.bestmapindex = i;
				tempBestMap.bestmapnum = sol[currnode->idx][i];
				bestMapIdx.push_back(tempBestMap);
			}
			else if (DParray[currnode->idx][i] == minReconCost)
			{
				totalOptSol = totalOptSol + sol[currnode->idx][i];
				tempBestMap.bestmapindex = i;
				tempBestMap.bestmapnum = sol[currnode->idx][i];
				bestMapIdx.push_back(tempBestMap);
			}
				
				
		}
		
		
		boost::multiprecision::cpp_int Index;
		Index = getRandom(totalOptSol);

		boost::multiprecision::cpp_int tempsum = 0;
		int x = 0;
		for (x = 0; x < bestMapIdx.size(); x++)
		{
			if ((Index > tempsum) && (Index <= (tempsum + bestMapIdx[x].bestmapnum)))
				break;
			tempsum = tempsum + bestMapIdx[x].bestmapnum;
		}
		
		
		if (!summary)
		{
			os << endl << "Reconciliation:";
			// call function to output the actual reconciliation 
			outputBackTrack(os, ((NamedSpeciesNode*) speciestree->nodes[bestMapIdx[x].bestmapindex]), currnode);
			
			os << endl << endl << "The minimum reconciliation cost is: " << minReconCost;
			os << " (Duplications: " << numDup << ", Transfers: " << numTrans << ", Losses: " << (minReconCost - (DCOST*numDup) - (TCOST*numTrans))/LCOST << ")" << endl;
			os << "Total number of optimal solutions: " << totalOptSol << endl;
//			os << LLONG_MAX << endl;
			os << "Total number of candidates for gene birth: " << bestMapIdx.size() << endl;
		}
		else   // otherwise collect the summary statistics 
		{
			SummaryBackTrack(os, ((NamedSpeciesNode*) speciestree->nodes[bestMapIdx[x].bestmapindex]), currnode);
			totalCost = totalCost + minReconCost;
			totalDup = totalDup + numDup;
			totalTrans = totalTrans + numTrans;
			totalLoss = totalLoss + (minReconCost - (DCOST*numDup) - (TCOST*numTrans))/LCOST;
			
		}
			
				

	}
	
	
	// perform the backtracking to output the reconciliation... in post-order.	
	void outputBackTrack(ostream &os, NamedSpeciesNode *CurrSpeciesNode, GeneNodeRooted *currnode)
	{
		if (currnode->isLeaf()) 
		{
//			os << endl << ((NamedGeneNodeRooted*)currnode)->getName() << ": " << "Leaf Node, Mapping --> " << CurrSpeciesNode->getName(); 
			os << endl << ((NamedGeneNodeRooted*)currnode)->getName() << ": " << "Leaf Node";  
			
			return;
		}
		outputBackTrack(os, ((NamedSpeciesNode*) backTrack[currnode->idx][CurrSpeciesNode->idx].mapping0), currnode->child(0));
		outputBackTrack(os, ((NamedSpeciesNode*) backTrack[currnode->idx][CurrSpeciesNode->idx].mapping1), currnode->child(1));
		
		if (!currnode->isLeaf())
		{
			if (backTrack[currnode->idx][CurrSpeciesNode->idx].event == 0)
				os << endl << "m" << currnode->idx + 1 << " = LCA[" << ((NamedGeneNodeRooted*)currnode->leftMost)->getName() << ", " << ((NamedGeneNodeRooted*)currnode->rightMost)->getName() << "]: " << "Speciation, Mapping --> " << CurrSpeciesNode->getName();  
			if (backTrack[currnode->idx][CurrSpeciesNode->idx].event == 1)
			{
				os << endl << "m" << currnode->idx + 1 << " = LCA[" << ((NamedGeneNodeRooted*)currnode->leftMost)->getName() << ", " << ((NamedGeneNodeRooted*)currnode->rightMost)->getName() << "]: " << "Duplication, Mapping --> " << CurrSpeciesNode->getName();  
				numDup++;
			}
			if (backTrack[currnode->idx][CurrSpeciesNode->idx].event == 2)
			{
				os << endl << "m" << currnode->idx + 1 << " = LCA[" << ((NamedGeneNodeRooted*)currnode->leftMost)->getName() << ", " << ((NamedGeneNodeRooted*)currnode->rightMost)->getName() << "]: " << "Transfer, Mapping --> " << CurrSpeciesNode->getName() << ", Recipient --> " << ((NamedSpeciesNode*)backTrack[currnode->idx][CurrSpeciesNode->idx].recipient)->getName();    
				numTrans++;
			}
		}
	}		
	
	
	
	// perform the backtracking to summarize the reconciliation... in post-order.	
	void SummaryBackTrack(ostream &os, NamedSpeciesNode *CurrSpeciesNode, GeneNodeRooted *currnode)
	{
		if (currnode->isLeaf()) 
		{
			//			os << endl << ((NamedGeneNodeRooted*)currnode)->getName() << ": " << "Leaf Node, Mapping --> " << CurrSpeciesNode->getName(); 
//			os << endl << ((NamedGeneNodeRooted*)currnode)->getName() << ": " << "Leaf Node";  
			
			return;
		}
		SummaryBackTrack(os, ((NamedSpeciesNode*) backTrack[currnode->idx][CurrSpeciesNode->idx].mapping0), currnode->child(0));
		SummaryBackTrack(os, ((NamedSpeciesNode*) backTrack[currnode->idx][CurrSpeciesNode->idx].mapping1), currnode->child(1));
		
		if (!currnode->isLeaf())
		{
			if (backTrack[currnode->idx][CurrSpeciesNode->idx].event == 0)
			{
//				os << endl << "LCA[" << ((NamedGeneNodeRooted*)currnode->leftMost)->getName() << ", " << ((NamedGeneNodeRooted*)currnode->rightMost)->getName() << "]: " << "Speciation, Mapping --> " << CurrSpeciesNode->getName();  
			}
			if (backTrack[currnode->idx][CurrSpeciesNode->idx].event == 1)
			{
//				os << endl << "LCA[" << ((NamedGeneNodeRooted*)currnode->leftMost)->getName() << ", " << ((NamedGeneNodeRooted*)currnode->rightMost)->getName() << "]: " << "Duplication, Mapping --> " << CurrSpeciesNode->getName();  
				numDup++;
			}
			if (backTrack[currnode->idx][CurrSpeciesNode->idx].event == 2)
			{
//				os << endl << "LCA[" << ((NamedGeneNodeRooted*)currnode->leftMost)->getName() << ", " << ((NamedGeneNodeRooted*)currnode->rightMost)->getName() << "]: " << "Transfer, Mapping --> " << CurrSpeciesNode->getName() << ", Recipient --> " << ((NamedSpeciesNode*)backTrack[currnode->idx][CurrSpeciesNode->idx].recipient)->getName();    
				numTrans++;
			}
		}
	}
	
		
	
/*	
	// return the reconciliation cost (for use with evaluating different rootings for unrooted gene trees)
	int returnReconciliationCost(GeneNodeRooted *currnode)
	{
		
		int minReconCost = INT_MAX / 10;
		int numSpecies_nodes = speciestree->nodes.size();
		
		for (int i = 0; i < numSpecies_nodes; i++)
		{
			if (DParray[currnode->idx][i] < minReconCost)
			{
				minReconCost = DParray[currnode->idx][i];
			}
		}
		
		return minReconCost;	

	}
	
*/
	
	void run(ostream &os, const ReRoot reroot, bool summary) 
	{
		bool rerooting = (reroot == ALL);
		createLeafMapping();
		speciestree->establishOrder();
		speciestree->assignIndex();
		
		
		// initialize the counters that store the total cost and number of events across all gene trees for the summary statistics.
		totalCost = 0;
		totalDup = 0;
		totalTrans = 0;
		totalLoss = 0;

		// initialize boosts random number generator
		gen.seed(rand());


		
//		set<int> genetreecount;
//		genetreecount.clear();
		
		msgout << "Computing...\n";


		doNodePreOrder(speciestree->root, 1);
		doBranchPreOrder(speciestree->root);
		
//SetSpeciesTreePOidx(speciestree->root, 0);
		
		// Call the nested post-order procedure to solve the DTL-problem
		for(int i2=0; i2<genetree_rooted.size(); i2++) 
		{
			GeneTreeRooted &tree = *genetree_rooted[i2];
			tree.assignIndex();
			
			// initialize the numDup and numTrans counters (to store the final number of duplications and transfers in the reconciliation)
			numDup = 0;
			numTrans = 0;
			
			// Initialize the dynamic programming table and in and out counters to infinity
			
			int numGene_nodes = tree.nodes.size();
			int numSpecies_nodes = speciestree->nodes.size();
			
			vector<int> row;
			row.clear();
			vector<boost::multiprecision::cpp_int> solrow;
			solrow.clear();
			
			vector<data> datarow;
			datarow.clear();
			data temp1;
			vector<datainout> datainoutrow;
			datainoutrow.clear();
			datainout temp2;
			
			
			for (int i = 0; i< numGene_nodes; i++)
			{		
				
				for (int j = 0; j< numSpecies_nodes;j++)
				{
					row.push_back(INT_MAX / 10);
					solrow.push_back(0);
					temp1.mapping0 = NULL;
					temp1.mapping1 = NULL;
					temp1.recipient = NULL;
					temp1.event = -1;
					datarow.push_back(temp1);
					
					temp2.optmapping = NULL;
					datainoutrow.push_back(temp2);
					
				}
				
				DParray.push_back(row);
				sol.push_back(solrow);
				
				backTrack.push_back(datarow);
				inTrack.push_back(datainoutrow);
				outTrack.push_back(datainoutrow);

				
				row.clear();
				solrow.clear();
				datarow.clear();
				datainoutrow.clear();
			}
			
			
			
			
			// base case of the dynamic programming table (incorporate leaf to leaf mapping)	
			vector<GeneNodeRooted*> &g = tree.nodes;
			SpeciesNode *temp_node;

		

			for(vector<GeneNodeRooted*>::iterator itr=g.begin(); itr!=g.end(); itr++) 
			{
				GeneNodeRooted *genenode = *itr;
				if (genenode->isLeaf())
				{
					SpeciesNode *&mapping = genenode->getMapping();
					DParray[genenode->idx][mapping->idx] = 0;
					
					sol[genenode->idx][mapping->idx] = 1;
					
					temp_node = genenode->getMapping();
					int count = 0;
					while (1)
					{
//						in[genenode->idx][temp_node->idx] = (count * LCOST);

						inTrack[genenode->idx][temp_node->idx].optmapping = genenode->getMapping();
						
						if (temp_node == speciestree->root)
							break;
						temp_node = temp_node->parent();
						count++;
						
					}
					
				}
			}
		
			
			if (tree.rooting == true)    // if the gene tree is rooted
			{
			
				// execute the main nested post order traversal
				PostOrderGeneTree(tree.root);
			
				
				if (!summary)
				{
					// Output the species tree and gene tree
					os << endl << endl << " ------------ Reconciliation for Gene Tree " << i2+1 << " (rooted) -------------" << endl;
					os << "Species Tree: " << endl; 
					speciestree->tree2newick(os);
					os << endl << endl;
					os << "Gene Tree: " << endl;
					tree.genetree2newick(os);
					os << endl;
				}
//	cout << tree << endl;
			
			
				// decorate the nodes of the gene tree with their leftmost and rightmost descendants.
				leftRightDescendants(tree.root);
			
				// output the final reconciliation cost
				OutputReconciliationCost(os,tree.root, summary);
				

			}

			
			
			
/*			
			
			
			
cout << endl<< "DParray " << DParray.size() << " " << DParray[0].size() << endl;
for (int i =0 ; i< DParray.size(); i++)
{
for (int j=0; j < DParray[0].size();j++)
{
cout << DParray[i][j] << " ";
}
cout << endl;			
}
			

			
			
			
cout << endl << "IN" << endl;			
//cout << endl << in.size() << " " << in[0].size() << endl;
for (int i =0 ; i< in.size(); i++)
{
for (int j=0; j < in[0].size();j++)
{
cout << in[i][j] << " ";
}
cout << endl;			
}
		
			

cout << endl << "IN TRACK" << endl;			
//cout << endl << in.size() << " " << in[0].size() << endl;
for (int i =0 ; i< inTrack.size(); i++)
{
for (int j=0; j < inTrack[0].size();j++)
{
if(	inTrack[i][j].optmapping == NULL)
cout << "NULL ";
else
cout << ((NamedSpeciesNode*)inTrack[i][j].optmapping)->getName() << " ";
}
cout << endl;			
}			
			
			
			


		
			
			
cout << endl << "OUT" << endl;			
//cout << endl << in.size() << " " << in[0].size() << endl;
for (int i =0 ; i< out.size(); i++)
{
for (int j=0; j < out[0].size();j++)
{
cout << out[i][j] << " ";
}
cout << endl;			
}
			


cout << endl << "OUT TRACK" << endl;			
//cout << endl << in.size() << " " << in[0].size() << endl;
for (int i =0 ; i< outTrack.size(); i++)
{
for (int j=0; j < outTrack[0].size();j++)
{
if(	outTrack[i][j].optmapping == NULL)
cout << "NULL ";
else
cout << ((NamedSpeciesNode*)outTrack[i][j].optmapping)->getName() << " ";
}
cout << endl;			
}
			
			
			
*/			
			
			
			
			
		
			// clean up the arrays
			for (int i =0 ; i< DParray.size(); i++)
			{
				DParray[i].clear();
				sol[i].clear();
				
				backTrack[i].clear();
				inTrack[i].clear();
				outTrack[i].clear();

				
			}
			DParray.clear();
			sol.clear();

			
			backTrack.clear();
			inTrack.clear();
			outTrack.clear();

			
		}
		
		
		
		if(summary)    // output the summary if needed
		{
			os << endl << endl << "Total reconciliation cost: " << totalCost;
			os << " (Duplications: " << totalDup << ", Transfers: " << totalTrans << ", Losses: " << (totalCost - (DCOST*totalDup) - (TCOST*totalTrans))/LCOST << ")" << endl;
			os << "Total number of optimal solutions: " << totalOptSol << endl;
			os << "Total number of candidates for gene birth: " << bestMapIdx.size() << endl;

		}
		

//	speciestree->tree2newick(os);
		
//	cout << endl << *speciestree << endl;
				
			
	}


	
	
	
	
	
};


#endif
