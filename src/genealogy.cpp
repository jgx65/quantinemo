/** @file genealogy.cpp
 *
 *   Copyright (C) 2010 Samuel Neuenschwander <samuel.neuenschwander@unil.ch>
 *   Copyright (C) 2018 Frederic Michaud <frederic.michaud@unil.ch>
 *
 *   quantiNemo:
 *   quantiNemo is an individual-based, genetically explicit stochastic
 *   simulation program. It was developed to investigate the effects of
 *   selection, mutation, recombination, and drift on quantitative traits
 *   with varying architectures in structured populations connected by
 *   migration and located in a heterogeneous habitat.
 *
 *   quantiNemo is built on the evolutionary and population genetics
 *   programming framework NEMO (Guillaume and Rougemont, 2006, Bioinformatics).
 *
 *
 *   Licensing:
 *   This file is part of quantiNemo.
 *
 *   quantiNemo is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   quantiNemo is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with quantiNemo.  If not, see <http://www.gnu.org/licenses/>.
 */


#include "genealogy.h"
//#include "tarray.h"
#include "tmetapop.h"
//#include "coal_deme.h"
//#include <algorithm>
using namespace std;

unsigned int TTNode::counter;

//------------------------------------------------------------------------------
/** recursive function to mutate through the tree */
void
TTNode::mutate(TMetapop* pop, const TTNode& ancestor, TLocus** genome)
{
    sequence=ancestor.sequence;
    unsigned int nbMut;
    nbMutations = nbMut = pop->rand().Poisson((*genome)->get_mut_rate()*(ancestor.time-time));
    for(; nbMut>0; --nbMut){                // for each mutation
        (*genome)->mutate_now(&sequence);
    }
    
    if(desc1) desc1->mutate(pop, *this, genome);
    if(desc2) desc2->mutate(pop, *this, genome);
    return;
}

//------------------------------------------------------------------------------
/** recursive function to find the recombination positions through the tree
 * if several recombinations occur on a branch only the lowest one (oldest) has to be considered
 * curSize is the cumulative length of investigated branches
 * vecNodes: map contains the nodes below a recombination event vecNodes<generationOfRecomb, Node>
 * vecRecomb: vector contains the sorted recombination positions (which were drawn randomly)
 *            disp_vec2 is used for the new branch (the part before the recombination is copied to it)
 * */
void
TTNode::get_recombinationPos(unsigned long& curSize, multimap<unsigned long, TTNode*>& vecNodes,
                             vector<unsigned long>& vecRecomb)
{
    if(vecRecomb.empty()) return;	// we have all recombination positions => stop
    
    // check if on this branch there is a recombination event
    curSize += ancestor->time-time;		// branch size above node
    if(vecRecomb.front()<curSize){		// if a recombination event occurs
        // get the corrected recombination time
        recombTime = vecRecomb.front() - curSize + ancestor->time;   	// correct to get it in generations (time)
        vecNodes.insert(pair<unsigned long, TTNode*>(recombTime, this));										// add the node to the vector
        
        // remove the recombination event and all subsequent recombination events on this branch as they are hidden by the first one (oldest one)
        while(vecRecomb.front()<curSize && !vecRecomb.empty()){
            vecRecomb.erase(vecRecomb.begin());
        }
        if(vecRecomb.empty()) return;					// all recombination positions have been found => stop here
    }
    
    // continue recursively searching for the recombination positions
    assert(!vecRecomb.empty());
    if(desc1) desc1->get_recombinationPos(curSize, vecNodes, vecRecomb);
    if(desc2) desc2->get_recombinationPos(curSize, vecNodes, vecRecomb);
    return;
}

//------------------------------------------------------------------------------
/** recursive function to find the recombination positions through the tree
 * when successful the recombining node is returned and the recombTime is set,
 * optherwise NULL is returned and the recombTime is set to my_NAN
 * curSize: cumulative length of the already investigated branches
 * recombTime: the cumulative time of the recombination (was drawn randomly)
 * gen: the generation the recombination event occurs (gen is NaN if not reached)
 * the return value is the node below the recombination event and NULL if not yet reached
 * */
TTNode*
TTNode::get_recombinationPosII(unsigned long& curSize, unsigned long recombTime, unsigned long& gen)
{
    // search recursively
    TTNode* recombNode;
    if(desc1){
        recombNode = desc1->get_recombinationPosII(curSize, recombTime, gen);
        if(recombNode) return recombNode;
    }
    if(desc2){
        recombNode = desc2->get_recombinationPosII(curSize, recombTime, gen);
        if(recombNode) return recombNode;
    }
    
    // check if on this branch there is a recombination event
    curSize += ancestor->time-time;		// branch size above node
    if(recombTime<curSize){		// if a recombination event occurs
        gen = recombTime - curSize + ancestor->time;   	// correct to get it in generations (time)
        return this;                                    // return this node
    }
    return NULL;
}

//------------------------------------------------------------------------------
/** recursive function to get the total length of the tree */
void
TTNode::compute_treeSize(unsigned long & size)
{
    assert(ancestor);
    size += ancestor->time-time;
    if(desc1) desc1->compute_treeSize(size);
    if(desc2) desc2->compute_treeSize(size);
    return;
}

//------------------------------------------------------------------------------
/** recursive function to mutate through the tree */
void
TTNode::write_tree(ostream& os, const unsigned int& tree_scale)
{
    if(desc1){                // if it is a branch
        assert(desc2);
        os << "(";
        desc1->write_tree(os, tree_scale);
        os << ", ";
        desc2->write_tree(os, tree_scale);
        os << ")";
    }
    else os << ID_Node+1;     // if it is a leave
    
    // write the scale
    if(!ancestor) return;                                       // no time scale for the last node
    if(tree_scale==1)  os << ": " << ancestor->time-time;       // scaled by the time of coalescence
    else               os << ": " << nbMutations;               // scaled by the number of mutations
}


//------------------------------------------------------------------------------
/** fill recursively the multimap with the nodes and their time */
void
TTNode::get_nodeTimes(multimap<unsigned int, TTNode*>& nodeMaps)
{
    if(desc1) desc1->get_nodeTimes(nodeMaps);
    if(desc2) desc2->get_nodeTimes(nodeMaps);
    nodeMaps.insert(pair<unsigned int, TTNode*>(time, this));
}

//------------------------------------------------------------------------------
/** return the demeID in which the node is at time "gen"
 *  could well be optimized... TODO
 * if multiple time entries the first one is returned (the one from coalescence)
 */
unsigned int
TTNode::get_currentDeme(unsigned long gen)
{
    // check if the time is within this branch
    assert((unsigned int) time <= gen);
    assert(!ancestor || (unsigned int) ancestor->time >= gen);
    
    //	cout << "\nNode " << ID_Node+1 << " (size: " << disp_vec.size() <<")" << flush;
    vector<pair<unsigned int, unsigned int> >::iterator end = disp_vec.end();;
    disp_vec_pos = disp_vec.begin();
    while(disp_vec_pos!=end && disp_vec_pos->first<gen) {
        //		cout << "\n  Gen " << disp_vec_pos->first << "; Deme " << disp_vec_pos->second << flush;
        ++disp_vec_pos;
    }
    if(disp_vec_pos!=end && disp_vec_pos->first==gen) ++disp_vec_pos; // to avoid if the node migrated just after coaelscence (2 identical time entries)
    
    disp_vec_next = disp_vec_pos;
    --disp_vec_pos;
    
    /*	if(disp_vec_pos == disp_vec.end()) cout << "\n==> Current pos: None" << flush;
     else cout << "\n==> Current pos (" << gen << "): " << disp_vec_pos->first << flush;
     if(disp_vec_next == disp_vec.end()) cout << "\n=> Next pos: None" << flush;
     else cout << "\n=> Next pos (" << gen << "): " << disp_vec_next->first << flush;
     */
    
    assert(disp_vec_pos->first <= gen);
    assert(disp_vec_next== disp_vec.end() || disp_vec_next->first >= gen);
    
    // return the corresponding demeID
    return disp_vec_pos->second;
}

//------------------------------------------------------------------------------
/** return the demeID in which the node is at time "gen"
 *  could well be optimized... TODO
 * if multiple time entries the last one is returned (the one from migration)
 */
unsigned int
TTNode::get_currentDeme_afterMigration(unsigned long gen)
{
    // check if the time is within this branch
    assert((unsigned int) time <= gen);
    assert(!ancestor || (unsigned int) ancestor->time >= gen);
    
    vector<pair<unsigned int, unsigned int> >::iterator end = disp_vec.end();;
    disp_vec_pos = disp_vec.begin();
    while(disp_vec_pos!=end && disp_vec_pos->first<=gen) {
        ++disp_vec_pos;
    }
    disp_vec_next = disp_vec_pos;
    --disp_vec_pos;
    
    assert(disp_vec_pos->first <= gen);
    assert(disp_vec_next== disp_vec.end() || disp_vec_next->first >= gen);
    
    // return the corresponding demeID
    return disp_vec_pos->second;
}

//------------------------------------------------------------------------------
/** return the demeID in which the node is at time "gen" if it is recombining (disp_vec_temp)
 * currently just just within assert()
 *  could well be optimized... TODO
 */
unsigned int
TTNode::get_currentDeme_temp(unsigned long gen)
{
    // check if the time is within this branch
    assert((unsigned int) time <= gen);
    //	assert(!ancestor || (unsigned int) ancestor->time > gen);   // the recombining node can overpass its ancestor
    
    //	cout << "\nNode " << ID_Node+1 << " (size: " << disp_vec.size() <<")" << flush;
    vector<pair<unsigned int, unsigned int> >::iterator end = disp_vec_temp.end();;
    vector<pair<unsigned int, unsigned int> >::iterator cur = disp_vec_temp.begin();
    while(cur!=end && cur->first<=gen) {
        //		cout << "\n  Gen " << disp_vec_pos->first << "; Deme " << disp_vec_pos->second << flush;
        ++cur;
    }
    
    --cur;
    
    /*	if(disp_vec_pos == disp_vec.end()) cout << "\n==> Current pos: None" << flush;
     else cout << "\n==> Current pos (" << gen << "): " << disp_vec_pos->first << flush;
     if(disp_vec_next == disp_vec.end()) cout << "\n=> Next pos: None" << flush;
     else cout << "\n=> Next pos (" << gen << "): " << disp_vec_next->first << flush;
     */
    
    assert(cur->first <= gen);
    assert(cur+1== disp_vec_temp.end() || (cur+1)->first > gen);
    
    // return the corresponding demeID
    return cur->second;
}

//------------------------------------------------------------------------------
/** return the demeID in which the node is at time "gen"
 * when the position is known
 */
unsigned int
TTNode::get_currentDeme(unsigned long gen, vector<pair<unsigned int, unsigned int> >::iterator& pos)
{
    // check if the time is within this branch
    if((unsigned int)  time > gen) return my_NAN;
    if((unsigned int)ancestor->time < gen) return my_NAN;
    assert(pos != disp_vec.end());
    
    // set the iterator
    disp_vec_pos=pos;
    disp_vec_next = pos;
    ++disp_vec_next;
    
    // return the corresponding demeID
    return disp_vec_pos->second;
}






