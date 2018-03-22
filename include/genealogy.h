/** @file genealogy.h
*
*   Copyright (C) 2010 Samuel Neuenschwander <samuel.neuenschwander@unil.ch>
*   Copyright (C) 2018 Frederic Michaud <frederic.michaud@unil.ch>

*   quantiNemo2:
*   quantiNemo2 is an individual-based, genetically explicit stochastic
*   simulation program. It was developed to investigate the effects of
*   selection, mutation, recombination, and drift on quantitative traits
*   with varying architectures in structured populations connected by
*   migration and located in a heterogeneous habitat.
*
*   quantiNemo2 is built on the evolutionary and population genetics
*   programming framework NEMO (Guillaume and Rougemont, 2006, Bioinformatics).
*
*
*   Licensing:
*   This file is part of quantiNemo2.
*
*   quantiNemo2 is free software: you can redistribute it and/or modify
*   it under the terms of the GNU General Public License as published by
*   the Free Software Foundation, either version 3 of the License, or
*   (at your option) any later version.
*
*   quantiNemo2 is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU General Public License for more details.
*
*   You should have received a copy of the GNU General Public License
*   along with quantiNemo2.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef genealogyH
#define genealogyH

#include "types.h"
//#include <vector>
//#include <assert.h>

#include <map>
#include <vector>

using namespace std;

class TLocus;
class TDeme;
class TMetapop;
// ------------------------------------------------------------------------------
/** node of the coalescence tree:
 *  - one ancestor
 *  - 2 leaves
 */

class TTNode {
public:
	unsigned long time;                                         // 0 is now and "time" goes backward in time
	TTNode *desc1;
	TTNode *desc2;
	TTNode *ancestor;
	unsigned int nbMutations;
    
    static unsigned int counter;

	vector<pair<unsigned int, unsigned int> > disp_vec; 	 	// store the IDs of each patch: disp_vec< time / demeID >
																// first entry is the patch and time of coalescence!!!
	vector<pair<unsigned int, unsigned int> > disp_vec_temp;	// this vector is used when recombining (this vector will exchange disp_vec
	vector<pair<unsigned int, unsigned int> >::iterator disp_vec_pos;	// current position in the disp_vec vector
	vector<pair<unsigned int, unsigned int> >::iterator disp_vec_next;	// next position in the disp_vec vector

	unsigned long recombTime;	// if number it is the time when a recombination event will happen on the ancestor branch

	unsigned char sequence;

	unsigned int ID_Node; 			// id node starting at 0

	void mutate(TMetapop* pop, const TTNode& ancestor, TLocus** genome);   // recursively function
	void write_tree(ostream& os, const unsigned int& tree_scale);   // recursively function
	void get_nodeTimes(multimap<unsigned int, TTNode*>&);   // recursively function
	void get_recombinationPos(unsigned long& curSize, multimap<unsigned long, TTNode*>& vecNodes,
                              vector<unsigned long>& vecRecomb);
	TTNode* get_recombinationPosII(unsigned long& curSize, unsigned long recombTime, unsigned long& gen);
	void compute_treeSize(unsigned long & size);
	unsigned int get_currentDeme(unsigned long gen);
    unsigned int get_currentDeme_afterMigration(unsigned long gen);
	unsigned int get_currentDeme(unsigned long gen, vector<pair<unsigned int, unsigned int> >::iterator& pos);
	unsigned int get_currentDeme_temp(unsigned long gen);

	TTNode() {
		++counter;
        time = 0;
		desc1 = NULL;
		desc2 = NULL;
		ancestor = NULL;
		sequence = (unsigned char) 0;
		nbMutations=0;
		ID_Node=my_NAN;
		recombTime = my_NAN;
		// disp_vec_pos=my_NAN;
	}

	TTNode(const TTNode& tt) {	// copy constructor
		++counter;
        time = tt.time;
		desc1 = tt.desc1;
		desc2 = tt.desc2;
		ancestor = tt.ancestor;
		sequence = tt.sequence;
		nbMutations = tt.nbMutations;
		ID_Node = tt.ID_Node;
		disp_vec_pos = tt.disp_vec_pos;
		disp_vec_next = tt.disp_vec_next;
		disp_vec = tt.disp_vec;
		disp_vec_temp = tt.disp_vec_temp;
		recombTime = tt.recombTime;
	}

	bool compare_disp_vec(const pair<unsigned int, unsigned int>& i, const unsigned long& time){
		return (i.first < (unsigned int)time);
	}


	~TTNode() {--counter;}
};

#endif  // genealogy_h sentry.
