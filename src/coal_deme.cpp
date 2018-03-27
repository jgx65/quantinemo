/** @file coal_deme.cpp
 *
 *   Copyright (C) 2008 Samuel Neuenschwander <samuel.neuenschwander@unil.ch>
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


#include "coal_deme.h"
#include "lce_coalescence.h"
//#include "treplicate.h"
#include <algorithm>
//#include <map>
//#include "stathandler.cpp"
//#include "version.h"
#include "lce_disperse.h"
//#include <limits.h>

#include <iostream>
using namespace std;

// ---------------------------------------------------------------------------
bool
dbCoalescence::check_db()
{
	map<PATCH_ID, dbPop>::iterator curPop, endPop;
	vector<pair<PATCH_ID, POP_SIZE> >::iterator curFrom, endFrom;
	POP_SIZE size, curPopSize;
	for(curPop=_dbPatches.begin(), endPop=_dbPatches.end(); curPop!=endPop; ++curPop){              // to
		curPopSize = curPop->second.get_popSize();

		// the "pop" deme has to have a population size
		if(!curPopSize){
			cout << "\nERROR: deme " << curPop->first << " should be populated!";
			return false;
		}

		// the sum of immigrants has to be smaler than the pop size
		for(size=0, curFrom=curPop->second.get_immigrants().begin(),
				endFrom=curPop->second.get_immigrants().end();
				curFrom!=endFrom; ++curFrom){  // from
			size += curFrom->second;
		}
		if(size > curPopSize){
			cout << "\nERROR: pop size (" << curPopSize << ") of deme "
					<< curPop->first << " is smaller than the number of immigrants ("<<size<<")!";
			return false;
		}
	}

	return true;
}


// ----------------------------------------------------------------------------------------
/** set the sample size (input is in diploids, Lineages is in haploids) */
void
TDeme::set_sample_size(const unsigned int& size, TTNode** sampleNodes, unsigned int& nbNodes, unsigned int demeID)
{
	_lineages = 2*size; // from diploid individuals to genes

	ChainNodeList.clear();
	TTNode* curNode;
	for(unsigned int i=0; i<_lineages; ++i, ++nbNodes){    // generate the corresponding nodes
		curNode = new TTNode();
		curNode->disp_vec.push_back(pair<unsigned int, unsigned int>(0, demeID));
		curNode->ID_Node = nbNodes;
		ChainNodeList.push_back(curNode);         // add it to the tree vector
		sampleNodes[nbNodes] = curNode;
	}
}

// ----------------------------------------------------------------------------------------
/** after migration the number of lineages is updated
 * returns the number of lineages present
 */
unsigned int
TDeme::post_migration_lineages()
{
	_lineages = (unsigned int)ChainNodeList.size();
	return _lineages;
}

// ----------------------------------------------------------------------------------------
/** performs migration of several lineages at once (to the same sink!)
 * the number of lineages will be later adjusted (second round)
 * the target deme is added to _demesLineags if it gets its first lineage
 *
 * compared to the function below the deme status is directly updated!!!!
 */
void
TDeme::migrate(const PATCH_ID& sinkID, unsigned int nbMigrants, LCE_Coalescence& coal)
{
	migrate(sinkID, nbMigrants, coal, coal.get_demesLineages());
	TDeme* sink = coal.get_demes(sinkID);
	if(!sink->post_migration_lineages()){     // remove the deme from _demesLineages if no lineages: should happen very seldom!!!
		vector<TDeme* >::reverse_iterator curDeme, endDeme;
		for(curDeme=coal.get_demesLineages().rbegin(), endDeme=coal.get_demesLineages().rend(); curDeme!=endDeme; ++curDeme){
			if(*curDeme==sink){
				coal.get_demesLineages().erase((++curDeme).base());
				break;
			}
		}
	}
}
// ----------------------------------------------------------------------------------------
/** performs migration of several lineages at once (to the same sink!)
 * the number of lineages will be later adjusted (second round)
 * the target deme is added to _demesLineags if it gets its first lineage
 * the current deme is not yet removed from _demesLineags if it looses its last lineage:
 * -> this is done in the second round
 */
void
TDeme::migrate(const PATCH_ID& sinkID, unsigned int nbMigrants, LCE_Coalescence& coal,
		vector<TDeme*>& newDemes)
{
	// get the corresponding deme
	TDeme* sink = coal.get_demes(sinkID);

	// check if the sink deme is the first time "colonized" by a lineage
	if(sink->ChainNodeList.empty() && !sink->get_lineages()){     // it can be that all lineages of the deme emigrated and at the end got soem new lineages
		newDemes.push_back(sink); 							     // add it to the traced container
	}

    // perform migration lineage by lineage
	unsigned int randNode;
	for(; nbMigrants; --nbMigrants, --_lineages){
		randNode = _popPtr->rand().Uniform(_lineages);            // get the lineage randomly
		sink->ChainNodeList.push_back(ChainNodeList[randNode]);    // add the lineage to the new deme
		ChainNodeList.erase(ChainNodeList.begin()+randNode);      // remove the lineage from the current deme
	}
}

// ----------------------------------------------------------------------------------------
/** performs migration of several lineages at once (to the same sink!)
 * the number of lineages will be later adjusted (second round)
 * the target deme is added to _demesLineags if it gets its first lineage
 * the current deme is not yet removed from _demesLineags if it looses its last lineage:
 * -> this is done in the second round
 */
void
TDeme::migrate(const PATCH_ID& sinkID, unsigned int nbMigrants, LCE_Coalescence& coal,
		map<PATCH_ID, TDeme*>& demes, map<PATCH_ID, TDeme*>& newDemes)
{
	// find the corresponding sink and if not yet present create it
	//	TDeme* curSink;
	//	map<PATCH_ID, TDeme*>::iterator pos = coal.get_demes().find(sinkID);
	//	if(pos != coal.get_demes().end()) curSink = pos->second;
	//	TDeme * curDeme = new TDeme(sinkID);
	//	newDemes.insert(pair<PATCH_ID, TDeme*>(sinkID, curDeme));

	TDeme* curSink = coal.get_demes(sinkID, newDemes);

	// perform migration lineage by lineage
	unsigned int randNode;
	for(; nbMigrants; --nbMigrants, --_lineages){
		randNode = _popPtr->rand().Uniform(_lineages);            		// get the lineage randomly
		migrate_recomb(curSink, randNode, coal);						// perform migration
    }
}


// ----------------------------------------------------------------------------------------
/** migrate a single gene from this patch to the sink */
void
TDeme::migrate(TDeme* sink, unsigned int geneID, LCE_Coalescence& coal)
{
	sink->ChainNodeList.push_back(ChainNodeList[geneID]);    // add the lineage to the new deme
	ChainNodeList.erase(ChainNodeList.begin()+geneID);      // remove the lineage from the current deme
}

// ----------------------------------------------------------------------------------------
/** migrate a single gene from this patch to the sink and store also the movement
 * do not yet update the number of lineages
 * */
void
TDeme::migrate_recomb(TDeme* sink , unsigned int geneID, LCE_Coalescence& coal)
{
	TTNode* curNode = ChainNodeList[geneID];
	sink->ChainNodeList.push_back(curNode);        				// add the lineage to the new deme
	ChainNodeList.erase(ChainNodeList.begin()+geneID);       // remove the lineage from the current deme
	curNode->disp_vec.push_back(pair<unsigned int, unsigned int>(coal.get_curGen(),sink->get_id()));

#ifdef _DEBUG
    cout << "\n  **** migration node " << curNode->ID_Node+1 << "' (deme "
    << _id+1 << " => " << sink->get_id()+1 << "; gen " << coal.get_curGen() << ")" << flush;
#endif
}

// ------------------------------------------------------------------------------
/** Implements the coalescence process within population */
unsigned int
TDeme::singleCoalescentEvent(const unsigned long& time, LCE_Coalescence& coal)
{
	assert(get_lineages());
	if(_lineages < 2) return 0;         // If there is only one lineage left, do not even try to coalesce it...
	assert(_deme_size);

	// Coalesce two lineages with probability proba
	if(_popPtr->rand().Uniform() < get_probCoal()){
		// Pick the first lineage at random
		unsigned int pick = _popPtr->rand().Uniform(_lineages);   // get the index of the first gene
		TTNode* desc1 = ChainNodeList[pick]; 	                 // get the first node
		ChainNodeList.erase(ChainNodeList.begin()+pick);       // remove this node from the vector

		// Pick the second lineage at random
		pick = _popPtr->rand().Uniform(_lineages-1);              // get the index of the second gene
		TTNode* desc2 = ChainNodeList[pick]; 	                 // get the second node
		ChainNodeList.erase(ChainNodeList.begin()+pick);       // remove this node from the vector

		coalesce(desc1, desc2, time, coal);                    // execute the coalescence

		coal.evacuate_lineages(this, false);                  // evacuate lineages if lin>N
		return 1;                                              // a single coalescent even happened
	}

	coal.evacuate_lineages(this, false);                    	// evacuate lineages if lin>N
	return 0;                                                	// no coalescent event happened
}

// ------------------------------------------------------------------------------
/** True coalescence process (i.e. no approximation), where parents are  randomly drawn for each gene,
 * thus there is the possibility of multiple coalescence events per generation
 * Disadvantage: much slower!
 */
unsigned int
TDeme::multipleCoalescentEvents(const unsigned long& time, LCE_Coalescence& coal)
{
	assert(get_lineages());
	assert(_deme_size);
	if(_lineages < 2) return 0;         // If there is only one lineage left, do not even try to coalesce it...

	// For every Node in the current generation
	// draw 1 ancestor among the _deme_size ancestors in the previous generation
	map<unsigned int, vector<unsigned int> > vec;    // 1: random number; 2: index of the corresponding nodes
	for(unsigned int i = 0; i < _lineages; ++i){
		vec[_popPtr->rand().Uniform(2*_deme_size)].push_back(i);
	}

	// Coalesce all nodes sharing the same ancestor
	unsigned int nbCoal = 0;
	vector<unsigned int> delNodes;			// we need this vector to know which nodes have to be deleted later (not very nice)
	map<unsigned int, vector<unsigned int> >::iterator curPos, endPos;
	for(curPos=vec.begin(), endPos=vec.end(); curPos!=endPos; ++curPos){
		if(curPos->second.size()>1){		// if several nodes are present
			// coalesce nodes with same parent, and attach this new ancestor at end of the vector
			// do not yet delete the coalesced nodes, since this would mix up the id of the nodes:
			// deletion is done afterwards
			nbCoal += curPos->second.size()-1;                      // number of coalescent events
			delNodes.insert(delNodes.end(), curPos->second.begin(), curPos->second.end()); // add the nodes to merge/delete

			TTNode* first = ChainNodeList[curPos->second.front()];  // take the first node as starting point
			for(unsigned int i=0, nb=(unsigned int)curPos->second.size()-1; i<nb; ++i){
				if(i) delNodes.push_back((unsigned int)ChainNodeList.size()-1);       // the ancestor will re-coalesce
				first = coalesce(first, ChainNodeList[curPos->second.back()], time, coal);
				curPos->second.pop_back();
			}
		}
	}

	// remove all the nodes which have coalesced in some way
	if(nbCoal){
		sort(delNodes.rbegin(), delNodes.rend());          // sort in reverse order: first the last node has to be deleted
		for(unsigned int i=0; i<delNodes.size(); ++i){
			ChainNodeList.erase(ChainNodeList.begin() + delNodes[i]);
		}
	}

	assert(coal.evacuate_lineages(this, false));

	return nbCoal;                            // return the number of coalescent events
}

// ------------------------------------------------------------------------------
/** Implements the coalescence process within population
 * Taking into account the recombination between loci
 * */
unsigned int
TDeme::mixedCoalescentEvents(const unsigned long& time, LCE_Coalescence& coal)
{
	assert(get_lineages());
	if(get_probCoal() < coal.get_model_threshold()){
		return singleCoalescentEvent(time, coal);
	}
	return multipleCoalescentEvents(time, coal);
}//
// ----------------------------------------------------------------------------------------
// ChainNodeList_remove
// ----------------------------------------------------------------------------------------
/** remove the given node from ChainNodeList
 * return false if node is not present and true if it was present and has been removed
 */
bool
TDeme::ChainNodeList_remove(TTNode* node)
{
	vector<TTNode*>::iterator cur, end;
	for(cur=ChainNodeList.begin(), end=ChainNodeList.end(); cur!=end; ++cur){
		if(*cur!=node) continue;
		ChainNodeList.erase(cur);
		return true;;
	}
	return false;
}

// ------------------------------------------------------------------------------
/** Implements the coalescence process within a population
 * and return a pointer to the ancestral node
 */
TTNode*
TDeme::coalesce(TTNode* desc1, TTNode* desc2, const unsigned long& time, LCE_Coalescence& coal)
{
	assert(!desc1->ancestor);
	assert(!desc2->ancestor);

	// generate the new node
	TTNode* ancestor = new TTNode;
	ancestor->ID_Node = coal.get_set_nbNodes()++;						// = coal.geneTree.size();
	ChainNodeList.push_back(ancestor);                   // add the node to the current nodes

	--_lineages;

	// Update node information
	ancestor->desc1 = desc1;
	ancestor->desc2 = desc2;
	ancestor->time = time;
	coal.set_mrca(ancestor);

	desc1->ancestor = ancestor;
	desc2->ancestor = ancestor;
    ancestor->ancestor = NULL;
    
    coal.add2treeSize((time-desc1->time)+(time-desc2->time));
    
#ifdef _DEBUG
    cout << "\nnewNode: " << ancestor->ID_Node << " (time " << time << "; "
    << desc1->ID_Node << " & " << desc2->ID_Node << ")" << endl;
#endif

	return ancestor;
}


// ------------------------------------------------------------------------------
/* it can happen that some lineages have to be evacuated if there are more lineages
 * in the deme than genes present. In this case the supernumerous lineages have to
 * evacuated to neighboring demes.
 * max is the population size in diploid numbers
 */
bool
TDeme::evacuate_lineages_island(LCE_Coalescence* pCoal, bool all)
{
	if(_lineages <= 2*_deme_size) return true;   // if there is an excess of lineages, evacuate them
	unsigned int mig_lineages = _lineages - 2*_deme_size; // the number of lineages to migrate
	warning("Coalescence: deme %i has an excess of %i lineages (pop size: %i, gen: %i): -> evacuating lineages!",_id+1,mig_lineages, _deme_size, pCoal->get_curGen());

	// draw for each lineage a potential neighbor (we don't have to pay attention to migrate to ourself, since this deme is full)
	unsigned int pos; 		// index for security reasons
	unsigned int nbNeighbours = (unsigned int)pCoal->get_demes().size();
	map<PATCH_ID, TDeme*>::iterator curPos;
	while(mig_lineages){
		pos = _popPtr->rand().Uniform(nbNeighbours);   // draw randomly a neighbor
		curPos = pCoal->get_demes().begin();
		while(pos){                                 // loop to this neighbor (not nice, but needed since it is a map)
			--pos;
			++curPos;
		}
		if(curPos->second->get_lineages() < 2*curPos->second->get_deme_size()){
			migrate(curPos->first, 1, *pCoal);     // migrate a single lineage
			--mig_lineages;
		}
	}
	return false;
}

// ------------------------------------------------------------------------------template<class NODE>
bool
TDeme::evacuate_lineages_disp_matrix(LCE_Coalescence* pCoal, bool all)
{
	error("Not yet implemetned!");
	/*
	if(_lineages <= 2*_deme_size) return true;   // if there is an excess of lineages, evacuate them
	unsigned int mig_lineages = _lineages - 2*_deme_size; // the number of lineages to migrate
	warning("Coalescence: deme %i has an excess of %i lineages (pop size: %i, gen: %i): -> evacuating lineages!",_id+1,mig_lineages, _deme_size);

	// find the last matrix with the last immigration to this deme
	unsigned int curGen = pCoal->get_nbGen()-pCoal->get_curGen()-1;   // check the current migration matrix
	map<PATCH_ID, dbPop>::iterator curPop;
	do{
		curPop = pCoal->get_dbCoalescence(curGen).get_patches().find(_id);
		if(curPop != pCoal->get_dbCoalescence(curGen).get_patches().end()) break; // we found a migation matix
		++curGen; // check the pervious generation (backward in time)
	}
	while(curGen < pCoal->get_nbGen());

	message("Coalescence: evacuating %i genes with migration matrix %i time difference!\n",
				 mig_lineages, pCoal->get_nbGen() - pCoal->get_curGen() - 1 - curGen);

	if(curGen==pCoal->get_nbGen()) return true; // the lineages in the deme did never immigrate

	// we have now a valid immigration matrix: evacuate the lineages with this matrix
	pCoal->migration(this, curPop->second.get_immigrants(), _deme_size, pCoal->get_demesLineages());
	 */
	return false;
}

// ------------------------------------------------------------------------------
bool
TDeme::evacuate_lineages_stepping_stone(LCE_Coalescence* pCoal, bool all)
{
	if(all) _deme_size = 0;
	if(_lineages <= 2*_deme_size) return true;   // if there is an excess of lineages, evacuate them
	unsigned int mig_lineages = _lineages - 2*_deme_size; // the number of lineages to migrate
	warning("Coalescence: deme %i has an excess of %i lineages (pop size: %i, gen: %i): -> evacuating lineages!",_id+1,mig_lineages, _deme_size, pCoal->get_curGen());
	if(!evacuate_stepping_stone(mig_lineages, 3, pCoal)){
		if(!pCoal->evacuate_stepping_stone(this, _id, mig_lineages, 3)){
			unsigned int newSize = _deme_size + (mig_lineages+1)/2;
			warning("Coalescence: Although all genes of the metapopulation are traced, there is an excess of  %i: lineages (gen: %i): -> pop size adjustement (pop %i: %i -> %i)!",mig_lineages,pCoal->get_curGen(),_id+1, _deme_size, newSize);
			_deme_size = newSize;
		}
	}
	return (_lineages==0);
}

// ------------------------------------------------------------------------------
/** recursively a deme is searched which can take the supernumerous lineages
 * false is returned if this is not possible, since all demes are full
 */
bool
TDeme::evacuate_stepping_stone(unsigned int lineages, unsigned int steps, LCE_Coalescence* pCoal)
{
	if(!lineages) return true;                                           // stop if there are no more lineages to spread
	vector<unsigned int> vNeighbours = pCoal->get_pop_ptr()->get_pDisperse_LCE()->get_neighbours(_id);  // get the neighbouring demes of this deme
	unsigned int i, size = (unsigned int)vNeighbours.size();                         // number of neighbours
	unsigned int* aNeighbours = new unsigned int[size];
	ARRAY::vector2array(vNeighbours, aNeighbours, size);
	_popPtr->rand().randomize(aNeighbours, size);                            // randomize the order of the neighbours

	// spread the lineages on the neighbours if they have free places
	map<PATCH_ID, TDeme*>::iterator cur, end = pCoal->get_demes().end();
	TDeme* curNeighbour;
	bool haveNeighbours = false;
	for(i=0; i<size && lineages; ++i){
		cur = pCoal->get_demes().find(aNeighbours[i]);
		if(cur == end) continue; // stop here if this neighbour is not present, i.e. not populated
		haveNeighbours = true;
		curNeighbour = cur->second;
		while(lineages && curNeighbour->get_lineages() <= 2*curNeighbour->get_deme_size()){
			migrate(curNeighbour->get_id(), 1, *pCoal);     // migrate a single lineage
			--lineages;
		}
	}

	// continue the spreading recursively from the neighbouring demes (for each lineage separately)
	if(steps && haveNeighbours){
		unsigned int max=10*size;
		for(i=0; lineages && max; i = (i==size-1) ? 0 : i+1, --max){  // increment by one and if it is the last one restart at 0
			cur = pCoal->get_demes().find(aNeighbours[i]);
			if(cur == end) continue;                      // if this neighbour is not present, i.e. not populated
			cur->second->evacuate_stepping_stone(1, steps-1, pCoal);    // a single lineage
			--lineages;
		}
	}
	delete[] aNeighbours;
	return (lineages==0);      // return true, if all lineages were spread and false if this was not possible
}


// ------------------------------------------------------------------------------
/** return true if present node is present in the chainNodeList
 */
bool
TDeme::isNodeInChainNodeList(TTNode* pNode)
{
    vector<TTNode*>::iterator cur, end;
    for(cur=ChainNodeList.begin(), end=ChainNodeList.end(); cur!=end; ++cur){
        if(*cur==pNode) return true;
    }
    return false;
}

// ------------------------------------------------------------------------------
/** return the iterator to the given node and one to end of gthe vector if node is not present
 */
vector<TTNode*>::iterator
TDeme::getNodePosInChainNodeList(TTNode* pNode)
{
    vector<TTNode*>::iterator cur, end;
    for(cur=ChainNodeList.begin(), end=ChainNodeList.end(); cur!=end; ++cur){
        if(*cur==pNode) break;
    }
    return cur;
}
