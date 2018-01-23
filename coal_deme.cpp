/** @file coal_deme.cpp
 *
 *   Copyright (C) 2008 Samuel Neuenschwander <samuel.neuenschwander@unil.ch>
 *
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
}

//// ------------------------------------------------------------------------------
///** Implements the coalescence process within population
// * for recombination (to build the first tree)
// * */
//unsigned int
//TDeme::singleCoalescentEventII(const unsigned long& time, LCE_CoalescenceRecomb& coal)
//{
//	assert(get_lineages());
//	if(_lineages < 2) return 0;         // If there is only one lineage left, do not even try to coalesce it...
//	assert(_deme_size);
//
//	// Coalesce two lineages with probability proba
//	if(_popPtr->rand().Uniform() < get_probCoal()){
//		// Pick the first lineage at random
//		unsigned int pick = _popPtr->rand().Uniform(_lineages);   // get the index of the first gene
//		TTNode* desc1 = ChainNodeList[pick]; 	             	// get the first node
//		ChainNodeList.erase(ChainNodeList.begin()+pick);       // remove this node from the vector
//
//		// Pick the second lineage at random
//		pick = _popPtr->rand().Uniform(_lineages-1);              	// get the index of the second gene
//		TTNode* desc2 = ChainNodeList[pick]; 	            	// get the second node
//		ChainNodeList.erase(ChainNodeList.begin()+pick);       	// remove this node from the vector
//
//		coalesce_recomb(desc1, desc2, time, coal, _id);        	// execute the coalescence
//
//		coal.evacuate_lineages(this, false);                  	// evacuate lineages if lin>N
//		return 1;                                              	// a single coalescent even happened
//	}
//
//	coal.evacuate_lineages(this, false);                    	// evacuate lineages if lin>N
//	return 0;                                                	// no coalescent even happened
//}
//
//// ------------------------------------------------------------------------------
///** True coalescence process (i.e. no approximation), where parents are  randomly drawn for each gene,
// * thus there is the possibility of multiple coalescence events per generation
// * Disadvantage: much slower!
// * for recombination (to build the first tree)
// */
//unsigned int
//TDeme::multipleCoalescentEventsII(const unsigned long& time, LCE_CoalescenceRecomb& coal)
//{
//	assert(get_lineages());
//	assert(_deme_size);
//	if(_lineages < 2) return 0;         // If there is only one lineage left, do not even try to coalesce it...
//
//	// For every Node in the current generation
//	// draw 1 ancestor among the _deme_size ancestors in the previous generation
//	map<unsigned int, vector<unsigned int> > vec;    // 1: random number; 2: index of the corresponding nodes
//	for(unsigned int i = 0; i < _lineages; ++i){
//		vec[_popPtr->rand().Uniform(2*_deme_size)].push_back(i);
//	}
//
//	// Coalesce all nodes sharing the same ancestor
//	unsigned int nbCoal = 0;
//	vector<unsigned int> delNodes;			// we need this vector to know which nodes have to be deleted later (not very nice)
//	map<unsigned int, vector<unsigned int> >::iterator curPos, endPos;
//	for(curPos=vec.begin(), endPos=vec.end(); curPos!=endPos; ++curPos){
//		if(curPos->second.size()>1){		// if several nodes are present
//			// coalesce nodes with same parent, and attach this new ancestor at end of the vector
//			// do not yet delete the coalesced nodes, since this would mix up the id of the nodes:
//			// deletion is done afterwards
//			nbCoal += curPos->second.size()-1;                      // number of coalescent events
//			delNodes.insert(delNodes.end(), curPos->second.begin(), curPos->second.end()); // add the nodes to merge/delete
//
//			TTNode* first = ChainNodeList[curPos->second.front()];  // take the first node as starting point
//			for(unsigned int i=0, nb=(unsigned int)curPos->second.size()-1; i<nb; ++i){
//				if(i) delNodes.push_back((unsigned int)ChainNodeList.size()-1);       // the ancestor will re-coalesce
//				first = coalesce_recomb(first, ChainNodeList[curPos->second.back()], time, coal, _id);
//				curPos->second.pop_back();
//			}
//		}
//	}
//
//	// remove all the nodes which have coalesced in some way
//	if(nbCoal){
//		sort(delNodes.rbegin(), delNodes.rend());          // sort in reverse order: first the last node has to be deleted
//		for(unsigned int i=0; i<delNodes.size(); ++i){
//			ChainNodeList.erase(ChainNodeList.begin() + delNodes[i]);
//		}
//	}
//
//	assert(coal.evacuate_lineages(this, false));
//
//	return nbCoal;                            // return the number of coalescent events
//}
//
//// ------------------------------------------------------------------------------
///** Implements the coalescence process within population
// * Taking into account the recombination between loci
// * for recombination (to build the first tree)
// * */
//unsigned int
//TDeme::mixedCoalescentEventsII(const unsigned long& time, LCE_CoalescenceRecomb& coal)
//{
//	assert(get_lineages());
//	if(get_probCoal() < coal.get_model_threshold()){
//		return singleCoalescentEventII(time, coal);
//	}
//	return multipleCoalescentEventsII(time, coal);
//}
//
//// ------------------------------------------------------------------------------
///** Implements the coalescence process within population
// * for recombination (for subsequent trees)
// */
//unsigned int
//TDeme::singleCoalescentEventIII(const unsigned long& time, LCE_CoalescenceRecomb& coal, vector<TTNode*>& recombNodes)
//{
//	assert(_deme_size);
//    
//	// Check if one of the active lineages will coalesce
//	unsigned int sizeRecomb = (unsigned int)recombNodes.size();
//	unsigned int curNbLineages = sizeRecomb + _lineages;
//	double proba = ((double)(curNbLineages)*(curNbLineages-1))/(4*_deme_size); // deme size is in diploid individuals
//	proba *= (double)sizeRecomb/curNbLineages;	// one of the recombining lineages has to be drawn...
//	if(_popPtr->rand().Uniform() < proba){
//		// Pick a recombining lineage at random
//		unsigned int pick2, pick1 = _popPtr->rand().Uniform(sizeRecomb);   // get the index of the first gene
//		TTNode* desc1 = recombNodes[pick1]; 	                    // get the first node
//        
//		// Pick the second old lineage at random
//		// Note: it cannot be another recombining lineage!!!!
//		TTNode* desc2;
//		for(unsigned int i=0; i<10000; ++i){ // 1000 for security reasons
//			pick2 = _popPtr->rand().Uniform(curNbLineages-1);    		// get the index of the second gene
//			if(pick2<_lineages){									// get the second node (it is a present lineage)
//				desc2 = ChainNodeList[pick2];
//				if(desc2->recombTime==RECOMB && desc2!=desc1) continue;	// may not recombine to other recombining trees
//				
//                // detach the nodes
//                recombNodes.erase(recombNodes.begin()+pick1);       // remove the first node from the vector
//                ChainNodeList_remove(desc1);
//                ChainNodeList.erase(ChainNodeList.begin()+pick2);	// remove the second node from the vector
//                
//				coalesce_recomb1(desc1, desc2, time, coal);         // execute the coalescence
//                return 1;
//			}
//			else {													// get the second node (it is an recombining lineage)
//				pick2 -= _lineages;									// correct the pick number
//				desc2 = recombNodes[pick2];
//               
//                // detach the nodes
//                recombNodes.erase(recombNodes.begin()+pick1);       // remove the first node from the vector
//                ChainNodeList_remove(desc1);
//				recombNodes.erase(recombNodes.begin()+pick2);		// remove the second node from the vector
//				ChainNodeList_remove(desc2);
//                
//				desc2 = coalesce_recomb1(desc1, desc2, time, coal); // execute the coalescence
//				recombNodes.push_back(desc2);						// it is not finished, put it back to the recombNodes.
//				return 1;                                           // a single coalescent even happened
//			}
//		}
//	}
//	return 0;                                                       // no coalescent even happened
//}
//
//// ------------------------------------------------------------------------------
///** True coalescence process (i.e. no approximation), where parents are  randomly drawn for each gene,
// * thus there is the possibility of multiple coalescence events per generation
// * Disadvantage: much slower!
// * for recombination (for subsequent trees)
// */
//unsigned int
//TDeme::multipleCoalescentEventsIII(const unsigned long& time, LCE_CoalescenceRecomb& coal, vector<TTNode*>& recombNodes)
//{
//	assert(_deme_size);
//	if(_lineages < 2) return 0;         // If there is only one lineage left, do not even try to coalesce it...
//    
//	// For every Node in the current generation
//	// draw 1 ancestor among the _deme_size ancestors in the previous generation
//	unsigned int sizeRecomb = (unsigned int)recombNodes.size();
//	map<unsigned int, vector<unsigned int> > vec;    // 1: random number; 2: index of the corresponding nodes
//	for(unsigned int i = 0; i < sizeRecomb+_lineages; ++i){
//		vec[_popPtr->rand().Uniform(2*_deme_size)].push_back(i);
//	}
//    
//	// Coalesce all nodes sharing the same ancestor
//	unsigned int nbCoal = 0;
//	vector<unsigned int> delNodes;			// we need this vector to know which nodes have to be deleted later (not very nice)
//	map<unsigned int, vector<unsigned int> >::iterator curPos, endPos;
//	for(curPos=vec.begin(), endPos=vec.end(); curPos!=endPos; ++curPos){
//		if(curPos->second.size()==1) continue; 		// nothing to coalesce
//        
//		// check if a recombining lineage recombines (the first indexes)
//		vector<unsigned int> toRecomb;
//		vector<unsigned int> oldLineage;
//		for(unsigned int i=0, nb=(unsigned int)curPos->second.size(); i<nb; ++i){
//			if(curPos->second[i] < sizeRecomb) toRecomb.push_back(i);
//			else oldLineage.push_back(i-sizeRecomb);
//		}
//        
//		if(toRecomb.empty()) continue;	// recombining lineages do not recombine
//        
//		//small problem: the old nodes are not allowed to recombine: recombine the recombining lineages with
//		TTNode *desc1, *desc2, *ancestor;
//		vector<unsigned int>::iterator curVec, endVec;
//		unsigned int pick1, pick2;
//        
//		if(oldLineage.empty()){		// there are no old lineages (recombine to the first one)
//			// get the first recombining lineage
//			curVec=toRecomb.begin();								// first node
//			pick1 = *curVec;   										// get the index of the first gene
//			desc1 = recombNodes[pick1]; 	                   		// get the first node
//			recombNodes.erase(recombNodes.begin()+pick1);         	// remove this node from the vector
//			ChainNodeList_remove(desc1);
//            
//			// for each other recombining node
//			for(++curVec, endVec=toRecomb.end(); curVec!=endVec; ++curVec){
//				// get second recombining lineage
//				pick2 = *curVec;
//				desc2 = recombNodes[pick2];
//				recombNodes.erase(recombNodes.begin()+pick2);		// remove this node from the vector
//				ChainNodeList_remove(desc2);
//                
//				// coalesce
//				ancestor = coalesce_recomb1(desc1, desc2, time, coal); // execute the coalescence
//				recombNodes.push_back(ancestor);					// it is not finished, put it back to the recombNodes.
//				++nbCoal;
//			}
//		}
//		else{						// there are old lineages: for each recombining lineage draw randomly an old one
//			unsigned int oldLineages_size = (unsigned int)oldLineage.size();
//			for(curVec=toRecomb.begin(), endVec=toRecomb.end(); curVec!=endVec; ++curVec){
//				// Get the recombining lineage
//				pick1 = *curVec;   										// get the index of the first gene
//				desc1 = recombNodes[pick1]; 	                   		// get the first node
//                
//				// Pick the second old lineage at random
//				// Note: it cannot be another recombining lineage!!!!
//				unsigned int i;
//				for(i=0; i<1000; ++i){
//					pick2 = _popPtr->rand().Uniform(oldLineages_size);    		// get the index of the second gene
//					desc2 = ChainNodeList[pick2];
//					if(desc2->recombTime!=RECOMB || desc2==desc1){	// may not recombine to other recombining trees
//                        // we have them: coalesce them (caution with the index!)
//                        ChainNodeList.erase(ChainNodeList.begin()+pick2);		// remove the 2. node from the vector
//                        if(desc1!=desc2) ChainNodeList_remove(desc1);			// remove the 1. node from the vector
//                        recombNodes.erase(recombNodes.begin()+pick1);         	// remove the 1. node from the vector
//                        coalesce_recomb1(desc1, desc2, time, coal);             // execute the coalescence
//                        ++nbCoal;
//                        break;
//                    }
//                }
//            }
//            
//		}
//	}
//	return nbCoal;                            // return the number of coalescent events
//}
//
//// ------------------------------------------------------------------------------
///** Implements the coalescence process within population
// * Taking into account the recombination between loci
// * for recombination (for subsequent trees)
// */
//unsigned int
//TDeme::mixedCoalescentEventsIII(const unsigned long& time, LCE_CoalescenceRecomb& coal, vector<TTNode*>& recombNodes)
//{
//	assert(get_lineages() || !recombNodes.empty());
//	if(get_probCoal((unsigned int)recombNodes.size()) < coal.get_model_threshold()){
//		return singleCoalescentEventIII(time, coal, recombNodes);
//	}
//	return multipleCoalescentEventsIII(time, coal, recombNodes);
//}
//
//// ------------------------------------------------------------------------------
///** Implements the coalescence process within population
// * for recombination (for subsequent trees)
// */
//unsigned int
//TDeme::singleCoalescentEventIV(const unsigned long& time, LCE_CoalescenceRecomb& coal, TTNode*& recombNode)
//{
//	assert(_deme_size);
//    
//	// Check if one of the active lineages will coalesce
//	unsigned int sizeRecomb = 1;
//	unsigned int curNbLineages = sizeRecomb + _lineages;
//	double proba = ((double)(curNbLineages)*(curNbLineages-1))/(4*_deme_size); // deme size is in diploid individuals
//	proba *= (double)sizeRecomb/curNbLineages;	// one of the recombining lineages has to be drawn...
//	if(_popPtr->rand().Uniform() < proba){
//        // Pick the second old lineage at random
//		// Note: it cannot be another recombining lineage!!!!
//		TTNode* desc2;
//        unsigned int pick2;
//		for(unsigned int i=0; i<10000; ++i){ // 1000 for security reasons
//			pick2 = _popPtr->rand().Uniform(_lineages);    		// get the index of the second gene
//            desc2 = ChainNodeList[pick2];
//            if(desc2->recombTime==RECOMB && desc2!=recombNode) continue;	// may not recombine to other recombining trees
//            
//            // detach the nodes
//            ChainNodeList.erase(ChainNodeList.begin()+pick2);           // remove the second node from the vector
//            if(recombNode!=desc2) ChainNodeList_remove(recombNode);     // remove the first node
//            
//            coalesce_recomb1(recombNode, desc2, time, coal);         // execute the coalescence
//            recombNode=NULL;
//            return 1;
//        }
//    }
//    return 0;                                                       // no coalescent even happened
//}
//
//// ------------------------------------------------------------------------------
///** True coalescence process (i.e. no approximation), where parents are randomly drawn for each gene,
// * thus there is the possibility of multiple coalescence events per generation
// * Disadvantage: much slower!
// * for recombination (for subsequent trees)
// */
//unsigned int
//TDeme::multipleCoalescentEventsIV(const unsigned long& time, LCE_CoalescenceRecomb& coal, TTNode*& recombNode)
//{
//	assert(_deme_size);
//	if(_lineages < 2) return 0;         // If there is only one lineage left, do not even try to coalesce it...
//    
//	// For every Node in the current generation
//	// draw 1 ancestor among the _deme_size ancestors in the previous generation
//	unsigned int sizeRecomb = 1;
//	map<unsigned int, vector<unsigned int> > vec;    // 1: random number; 2: index of the corresponding nodes
//	for(unsigned int i = 0; i < sizeRecomb+_lineages; ++i){
//		vec[_popPtr->rand().Uniform(2*_deme_size)].push_back(i);
//	}
//    
//	// Coalesce all nodes sharing the same ancestor
//	vector<unsigned int> delNodes;			// we need this vector to know which nodes have to be deleted later (not very nice)
//	map<unsigned int, vector<unsigned int> >::iterator curPos, endPos;
//	for(curPos=vec.begin(), endPos=vec.end(); curPos!=endPos; ++curPos){
//		if(curPos->second.size()==1) continue; 		// nothing to coalesce
//        
//		// check if a recombining lineage recombines (the first indexes)
//		vector<unsigned int> toRecomb;
//		vector<unsigned int> oldLineage;
//		for(unsigned int i=0, nb=(unsigned int)curPos->second.size(); i<nb; ++i){
//			if(curPos->second[i] < sizeRecomb) toRecomb.push_back(i);
//			else oldLineage.push_back(i-sizeRecomb);
//		}
//
//		if(toRecomb.empty()) continue;	// recombining lineage does not recombine
//		assert(!oldLineage.empty());
//        
//		// the node coalescese: Pick the second old lineage at random
//        unsigned int pick2 = _popPtr->rand().Uniform((unsigned int)oldLineage.size());    		// get the index of the second gene
//        TTNode *desc2 = ChainNodeList[pick2];                   // get the second node
//        ChainNodeList.erase(ChainNodeList.begin()+pick2);		// remove the 2. node from the vector
//        if(recombNode!=desc2) ChainNodeList_remove(recombNode);	// remove the 1. node from the vector
//        coalesce_recomb1(recombNode, desc2, time, coal);        // execute the coalescence
//        recombNode = NULL;
//        return 1;                            // return the number of coalescent events
//	}
//	return 0;                                // no coaelscence occured
//}
//
//// ------------------------------------------------------------------------------
///** Implements the coalescence process within population
// * Taking into account the recombination between loci
// * for recombination (for subsequent trees)
// */
//unsigned int
//TDeme::mixedCoalescentEventsIV(const unsigned long& time, LCE_CoalescenceRecomb& coal, TTNode*& recombNode)
//{
//	assert(get_lineages() || recombNode);
//	if(get_probCoal(1) < coal.get_model_threshold()){
//		return singleCoalescentEventIV(time, coal, recombNode);
//	}
//	return multipleCoalescentEventsIV(time, coal, recombNode);
//}
//
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

//// ------------------------------------------------------------------------------
///** As function coalesce(), however also the positions in the patches over time are stored
// */
//TTNode*
//TDeme::coalesce_recomb(TTNode* desc1, TTNode* desc2, const unsigned long& time, LCE_CoalescenceRecomb& coal, unsigned int demeID)
//{
//	TTNode* ancestor = coalesce(desc1, desc2, time, coal);
//	ancestor->disp_vec.push_back(pair<unsigned int, unsigned int>(time, demeID));
//	return ancestor;
//}
//
///// ------------------------------------------------------------------------------
///** node1 is always a recombining branch
// *  node2 is a present branch
// *  the total number of nodes does not change:
// *  if a new node is created an old one is removed, since one of the children is removed,
// *  we reuse the removed node for the new node
// */
//TTNode*
//TDeme::coalesce_recomb1(TTNode* node1, TTNode* node2, const unsigned long& time, LCE_CoalescenceRecomb& coal)
//{
//#ifdef _DEBUG
//	cout << "\n********* Recombination at generation " << time <<
//			" of nodes " << node1->ID_Node+1 << "' and " <<
//			node2->ID_Node+1 << flush;
//#endif
//
//	TTNode *other, *obsolete;
//	--_lineages;
//    assert(node1->recombTime == RECOMB);
//    node1->recombTime = my_NAN;
//    
//	///////////////////////////////////////////////////////////////////////////
//	// 1. recombining branch recombines back to itself: update disp_vec, all the rest does not change
//	if(node1==node2){
//        obsolete = recombine_back(node1, node2, time, coal);
//		node1->get_currentDeme(time); // update iterators for disp_vec
//#ifdef _DEBUG
//		cout << " back into themselves (ancestor does not change)" << flush;
//        coal.test_tree();
//#endif
//		return obsolete;
//	}
//    
//	///////////////////////////////////////////////////////////////////////////
//	// 2. recombining branch to another lineage
//
//    // get the involved nodes
//    obsolete = node1->ancestor;
//    other = obsolete->desc1==node1 ? obsolete->desc2 : obsolete->desc1;
//    
//    // adapt the recombining disp_vec node
//    node1->disp_vec = node1->disp_vec_temp;
//    node1->disp_vec_temp.clear();
//	coal.add2treeSize(time - obsolete->time);
//    
//    // remove obsolete node: adapt other child
//    vector<pair<unsigned int, unsigned int> >::iterator cur, end;
//    if(!obsolete->ancestor){                // the obsolete node is the mrca: now "other" is the mrca
//        if(node2 == other){    // if the mrca is the ancestor of both nodes: readjust the mrca
//            coal.add2treeSize(time - obsolete->time);   // the second part has also to be removed
//            obsolete->time = time;
//
//            // split disp_vec of node2
//            cur = node2->disp_vec.begin();
//            end = node2->disp_vec.end();;
//            while(cur!=end && cur->first<(unsigned int)time) {
//                ++cur;
//            }
//            //if(cur->first==(unsigned int)time) ++cur;
//            node2->disp_vec.erase(cur, end);                                            // remove upper part
//            assert(node2->disp_vec.back().first < node2->ancestor->time);
//            
//            // adjust disp_vec of the mrca
//            obsolete->disp_vec.clear();
//            obsolete->disp_vec.push_back(pair<unsigned int, unsigned int>(time, _id));  // node entry
//            
//            coal.set_mrca(obsolete);
//            return obsolete;
//        }
//        
//        coal.add2treeSize(other->time - obsolete->time);
//        other->disp_vec.erase(other->disp_vec.begin()+1, other->disp_vec.end());
//        coal.set_mrca(other);
//    }
//    else{
//        if(obsolete->ancestor->desc1 == obsolete) obsolete->ancestor->desc1 = other;
//        else obsolete->ancestor->desc2 = other;
//        cur = obsolete->disp_vec.begin();
//        if(cur->second == other->disp_vec.back().second) ++cur;
//        other->disp_vec.insert(other->disp_vec.end(), cur, obsolete->disp_vec.end());
//        if(node2==obsolete) node2 = other;
//    }
//    other->ancestor = obsolete->ancestor;
//    
//    // insert new node
//    obsolete->ancestor = node2->ancestor;
//    obsolete->desc1 = node1;
//    obsolete->desc2 = node2;
//    obsolete->time = time;
//    if(obsolete->ancestor->desc1 == node2) obsolete->ancestor->desc1 = obsolete;
//    else {
//        assert(obsolete->ancestor->desc2 == node2);
//        obsolete->ancestor->desc2 = obsolete;
//    }
//    node2->ancestor = obsolete;
//    node1->ancestor = obsolete;
//    
//	// split disp_vec of node2
//    cur = node2->disp_vec.begin();
//    end = node2->disp_vec.end();;
//	while(cur!=end && cur->first<(unsigned int)time) {
//		++cur;
//	}
//    obsolete->disp_vec.clear();
//    obsolete->disp_vec.push_back(pair<unsigned int, unsigned int>(time, _id));  // node entry
//    if(cur->second==_id) obsolete->disp_vec.insert(obsolete->disp_vec.end(), cur+1, end);              // upper part
//    else obsolete->disp_vec.insert(obsolete->disp_vec.end(), cur, end);
//    node2->disp_vec.erase(cur, end);                                            // remove upper part
//    assert(node2->disp_vec.back().first < node2->ancestor->time);
//    
//    // if the new node is the MRCA
//    if(!obsolete->ancestor){
//        
//    }
//    
//#ifdef _DEBUG
//	cout << " into node " << obsolete->ID_Node+1 << flush;
//#endif
//    
//    return obsolete;
////    
/*
    
	if(cur!=end && cur->second==node2->disp_vec_temp.back().second) ++cur; // remove one more if the next change is obsolete (the same)
    node1->disp_vec.erase(node1->disp_vec.begin(), cur); // remove the beginning
    assert(node1->disp_vec.empty() || node1->disp_vec.front().first>time);    // the first part is removed: the vector start has now to be larger than "time"
    
	node1->disp_vec.insert(node1->disp_vec.begin(), node1->disp_vec_temp.begin(), node1->disp_vec_temp.end()); // insert the new part of the branch
	node1->disp_vec_temp.clear();	// erase the temp branch;
	return node1->ancestor;
    
    

    for(cur=node2->disp_vec.begin(), end=node2->disp_vec.end(); cur!=end; ++cur){
        if((unsigned long)cur->first < (unsigned long)time) continue;
        
        if(cur->first != (unsigned int)time){	// the disp_vec has to start with the deme position of coalescence
            temp1->disp_vec.push_back(pair<unsigned int, unsigned int>(time, _id));
        }
        temp1->disp_vec.insert(temp1->disp_vec.end(), cur, end); 	// copy the part after the new node to the ancestor
        node2->disp_vec.erase(cur, end);							// remove the same part from desc2
        break;
    }
    if(cur==end){	// if there was no entry
        assert(temp1->disp_vec.empty());
        temp1->disp_vec.push_back(pair<unsigned int, unsigned int>(time, _id));
    }

    
    
    

	// remove obsolete node, update disp_vec
	temp1 = recombining_branch(node1, time, coal);	// always a recombining branch
    
    
    // caution if the coalescence happens on the obsolete node: move it to other
    if(temp1==node2) node2 = temp1->ancestor->desc1==temp1 ? temp1->ancestor->desc2 : temp1->ancestor->desc1;
    
    // if the oboslete node is the MRCA and is not the ancestor of node2
    if(!temp1->ancestor && node2->ancestor){
        // the other child is the MRCA
        temp2 = temp1->desc1==node1 ? temp1->desc2 : temp1->desc1;
        long diffTime = temp2->time - coal.get_tmrca();                        // get the change in tMRCA (must be negative: smaller tMRCA now)
        coal.add2treeSize(diffTime);                                    // since the MRCA is moved the size of the tree has to be adjusted
        coal.set_mrca(temp2);
        
        // the MRCA has changed
        PATCH_ID mrca_patch = temp2->disp_vec.front().second;
        temp2->disp_vec.clear();
        temp2->disp_vec.push_back(pair<unsigned int, unsigned int>(temp2->time, mrca_patch));
    }

	// set up the new ancestor (reuse the removed node temp1)
	if(!temp1) temp1 = new TTNode;							// if not present create it
	temp1->desc1 = node1;
	temp1->desc2 = node2;
	temp1->time = time;
	temp1->ancestor = node2->ancestor;						// it can be that the ancestor is NULL
    
	// set child of ancestor
	if(node2->ancestor){
		if(node2->ancestor->desc1 == node2) node2->ancestor->desc1 = temp1;
		else node2->ancestor->desc2 = temp1;
        
		///////////////////////////////////////////////////////////////////////////
        // 2: node2 is a present branch: split disp_vec between node2 - temp1 - node2->ancestor
        temp1->disp_vec.clear();		// delete the present vector
        vector<pair<unsigned int, unsigned int> >::iterator cur, end;
        for(cur=node2->disp_vec.begin(), end=node2->disp_vec.end(); cur!=end; ++cur){
            if((unsigned long)cur->first < (unsigned long)time) continue;
            
            if(cur->first != (unsigned int)time){	// the disp_vec has to start with the deme position of coalescence
                temp1->disp_vec.push_back(pair<unsigned int, unsigned int>(time, _id));
            }
            temp1->disp_vec.insert(temp1->disp_vec.end(), cur, end); 	// copy the part after the new node to the ancestor
            node2->disp_vec.erase(cur, end);							// remove the same part from desc2
            break;
        }
        if(cur==end){	// if there was no entry
            assert(temp1->disp_vec.empty());
            temp1->disp_vec.push_back(pair<unsigned int, unsigned int>(time, _id));
        }
    }
    else{       // the MRCA is changing
        long diffTime = time - coal.get_tmrca();                        // get the change in tMRCA (has to be negative!!!: smaller tMRCA now)
        assert(diffTime<=0);                                             // the dispersal vector of the childern has to be adapted
        coal.add2treeSize(diffTime);                                    // since the MRCA is moved the size of the tree has to be adjusted

		temp1->disp_vec.clear();
        temp1->disp_vec.push_back(pair<unsigned int, unsigned int>(time, _id));
        coal.set_mrca(temp1);                                           // the MRCA has changed
        
        assert(node1->disp_vec.back().first<time);                      // the recombing node should be ok
        vector<pair<unsigned int, unsigned int> >::iterator cur, end;   // the enot recoming node has to be adapted for the shorter branch
        for(cur=node2->disp_vec.begin(), end=node2->disp_vec.end(); cur!=end; ++cur){
            if((unsigned long)cur->first < (unsigned long)time) continue;
            node2->disp_vec.erase(cur, end);							// remove the last part
            break;
        }
        assert(node2->disp_vec.back().first<time);
#ifdef _DEBUG
		message("\n  ***** The time to the MRCA has changed to %i (node %i)\n", time, coal.get_mrca()->ID_Node+1);
#endif
    }
    
	// set the ancestor for both nodes
	node1->ancestor = temp1;
	node2->ancestor = temp1;
    
	///////////////////////////////////////////////////////////////////////////
	// 3. node2 is also a recombining branch
	if(node2->recombTime == RECOMB)
    {
        assert(node2->recombTime == RECOMB);
        node2->recombTime = my_NAN;
        
		temp2 = recombining_branch(node2, time, coal);
		if(temp2) {		// if present delete it: we need currently just one new node
			delete temp2;
			temp2=NULL;
		}
        
		///////////////////////////////////////////////////////////////////////////
		temp1->get_currentDeme(time);   // check that the current deme is present
	}
    
    // if there is a recombination event on node2 at a later time add it to the new node
    else if(node2->recombTime!=my_NAN)
    {
        assert(node2->recombTime > time);
        temp1->recombTime = node2->recombTime;
        node2->recombTime = my_NAN;
        
        // adapt the vector _vecNodes (replace the obsolete node by other)
        map<unsigned long, TTNode*>::iterator cur = coal._vecNodes_pos;
        for(; cur!= coal._vecNodes_end; ++cur){
            if(cur->second != node2) continue;
            cur->second = temp1;
            break;
        }
    }
    
	ChainNodeList.push_back(temp1);

#ifdef _DEBUG
	cout << " into node " << temp1->ID_Node+1 << flush;
#endif

	return temp1;
*/
//}


//// ------------------------------------------------------------------------------
///** perform recombination if the recombining branch recombines with the same branch
// * (tree does not change, but disp_vec)
// */
//TTNode*
//TDeme::recombine_back(TTNode* node1, TTNode* node2, const unsigned long& time, LCE_CoalescenceRecomb& coal)
//{
//	assert(node1==node2);
//
//	// if it recombines back to the same branch (tree does not change, but disp_vec)
//	// just insert the new dispersal vector
//    vector<pair<unsigned int, unsigned int> >::iterator cur = node1->disp_vec.begin();
//    vector<pair<unsigned int, unsigned int> >::iterator end = node1->disp_vec.end();;
//	while(cur!=end && cur->first<(unsigned int)time) {
//		++cur;
//	}
//	if(cur!=end && cur->second==node1->disp_vec_temp.back().second) ++cur; // remove one more if the next change is obsolete (the same)
//    node1->disp_vec.erase(node1->disp_vec.begin(), cur); // remove the beginning
//    assert(node1->disp_vec.empty() || node1->disp_vec.front().first>=time);    // the first part is removed: the vector start has now to be larger than "time"
//
//	node1->disp_vec.insert(node1->disp_vec.begin(), node1->disp_vec_temp.begin(), node1->disp_vec_temp.end()); // insert the new part of the branch
//	node1->disp_vec_temp.clear();	// erase the temp branch;
//	return node1->ancestor;
//}

//// ------------------------------------------------------------------------------
///** if it is a recombining branch two things have to be done:
// * 1. move disp_vec_temp to disp_vec
// * 2. remove potential obsolete old ancestor (and return it)
// * returns the obsolete node (not anymore used) if available
// */
//TTNode*
//TDeme::recombining_branch(TTNode* node, const unsigned long& time, LCE_CoalescenceRecomb& coal)
//{
//	if(node->disp_vec_temp.empty()) return NULL; 	// not a recombining branch
//
//	// 1. move disp_vec_temp to disp_vec
//	node->disp_vec = node->disp_vec_temp;
//	node->disp_vec_temp.clear();
//	if(!node->ancestor) return NULL;  				// no ancestor => no obsolete node
//
//	// 2. remove obsolete old ancestor
//	TTNode* obsolete = node->ancestor;
//	TTNode* other =  (obsolete->desc1==node) ? obsolete->desc2 : obsolete->desc1; // which child is the other one
//	other->ancestor = obsolete->ancestor;
//	if(obsolete->ancestor){
//		if(obsolete->ancestor->desc1==obsolete) obsolete->ancestor->desc1 = other;
//		else obsolete->ancestor->desc2 = other;
//	}
//
//	// if there is a recombination event on the ancestor branch of the obsolete node add it to other
//	if(obsolete->recombTime!=my_NAN){
//		if(other->recombTime==my_NAN || other->recombTime > obsolete->recombTime){
//			other->recombTime = obsolete->recombTime;
//
//			// adapt the vector _vecNodes (replace the obsolete node by other)
//			multimap<unsigned long, TTNode*>::iterator cur = coal._vecNodes_pos;
//			for(; cur!= coal._vecNodes_end; ++cur){
//				if(cur->second != obsolete) continue;
//				cur->second = other;
//				break;
//			}
//		}
//		obsolete->recombTime = my_NAN;
//	}
//
//	// if the obsolete node is not the MRCA update the disp_vec (remove obsolet entries...)
//	if(other->ancestor){
//        // clean up the merge if the entry is meaningless
//        vector<pair<unsigned int, unsigned int> >::iterator start = obsolete->disp_vec.begin();
//        if(other->disp_vec.back().second== obsolete->disp_vec.front().second) ++start; // ignore the first entry
//        other->disp_vec.insert(other->disp_vec.end(), start, obsolete->disp_vec.end());
//    }
//	return obsolete;
//}
//
///// ------------------------------------------------------------------------------
///** desc1 is always a recombining branch
// *  desc2 is a present branch OR a recombining branch (defined by "recomb2")
// *  the total number of nodes does not change:
// *  if a new node is created an old one is removed, since one of the children is removed,
// *  we reuse the removed node for the new node
// */
//TTNode*
//TDeme::coalesce_recomb2(TTNode* node1, TTNode* node2, const unsigned long& time, LCE_CoalescenceRecomb& coal)
//{
//#ifdef _DEBUG
//	cout << "\n********* Recombination at generation " << time <<
//    " of nodes " << node1->ID_Node+1 << "' and " <<
//    node2->ID_Node+1 << flush;
//#endif
//    
//	TTNode *temp1, *temp2;
//	--_lineages;
//    assert(node1->recombTime == RECOMB);
//    node1->recombTime = my_NAN;
//    
//	///////////////////////////////////////////////////////////////////////////
//	// 1. recombining branch recombines back to itself: update disp_vec, all the rest does not change
//	if(node1==node2) {
//        temp1 = recombine_back(node1, node2, time, coal);
//		node1->get_currentDeme(time); // update iterators for disp_vec
//#ifdef _DEBUG
//		cout << " back into themselves (ancestor does not change)" << flush;
//#endif
//        coal.test_tree();
//		return temp1;
//	}
//    
//	///////////////////////////////////////////////////////////////////////////
//	// recombining branch to another lineage
//	coal.add2treeSize(time - node1->ancestor->time);
//    
//	// remove obsolete node, update disp_vec
//	temp1 = recombining_branch(node1, time, coal);	// always a recombining branch
//    
//    
//    // caution if the coalescence happens on the obsolete node: move it to other
//    if(temp1==node2) node2 = temp1->ancestor->desc1==temp1 ? temp1->ancestor->desc2 : temp1->ancestor->desc1;
//    
//    // if the oboslete node is the MRCA and is not the ancestor of node2
//    if(!temp1->ancestor && node2->ancestor){
//        // the other child is the MRCA
//        temp2 = temp1->desc1==node1 ? temp1->desc2 : temp1->desc1;
//        long diffTime = temp2->time - coal.get_tmrca();                        // get the change in tMRCA (must be negative: smaller tMRCA now)
//        coal.add2treeSize(diffTime);                                    // since the MRCA is moved the size of the tree has to be adjusted
//        coal.set_mrca(temp2);
//        
//        // the MRCA has changed
//        PATCH_ID mrca_patch = temp2->disp_vec.front().second;
//        temp2->disp_vec.clear();
//        temp2->disp_vec.push_back(pair<unsigned int, unsigned int>(temp2->time, mrca_patch));
//    }
//    
//	// set up the new ancestor (reuse the removed node temp1)
//	if(!temp1) temp1 = new TTNode;							// if not present create it
//	temp1->desc1 = node1;
//	temp1->desc2 = node2;
//	temp1->time = time;
//	temp1->ancestor = node2->ancestor;						// it can be that the ancestor is NULL
//    
//	// set child of ancestor
//	if(node2->ancestor){
//		if(node2->ancestor->desc1 == node2) node2->ancestor->desc1 = temp1;
//		else node2->ancestor->desc2 = temp1;
//        
//		///////////////////////////////////////////////////////////////////////////
//        // 2: node2 is a present branch: split disp_vec between node2 - temp1 - node2->ancestor
//        temp1->disp_vec.clear();		// delete the present vector
//        vector<pair<unsigned int, unsigned int> >::iterator cur, end;
//        for(cur=node2->disp_vec.begin(), end=node2->disp_vec.end(); cur!=end; ++cur){
//            if((unsigned long)cur->first < (unsigned long)time) continue;
//            
//            if(cur->first != (unsigned int)time){	// the disp_vec has to start with the deme position of coalescence
//                temp1->disp_vec.push_back(pair<unsigned int, unsigned int>(time, _id));
//            }
//            temp1->disp_vec.insert(temp1->disp_vec.end(), cur, end); 	// copy the part after the new node to the ancestor
//            node2->disp_vec.erase(cur, end);							// remove the same part from desc2
//            break;
//        }
//        if(cur==end){	// if there was no entry
//            assert(temp1->disp_vec.empty());
//            temp1->disp_vec.push_back(pair<unsigned int, unsigned int>(time, _id));
//        }
//    }
//    else{       // the MRCA is changing
//        long diffTime = time - coal.get_tmrca();                        // get the change in tMRCA (has to be negative!!!: smaller tMRCA now)
//        assert(diffTime<=0);                                             // the dispersal vector of the childern has to be adapted
//        coal.add2treeSize(diffTime);                                    // since the MRCA is moved the size of the tree has to be adjusted
//        
//		temp1->disp_vec.clear();
//        temp1->disp_vec.push_back(pair<unsigned int, unsigned int>(time, _id));
//        coal.set_mrca(temp1);                                           // the MRCA has changed
//        
//        assert(node1->disp_vec.back().first<time);                      // the recombing node should be ok
//        vector<pair<unsigned int, unsigned int> >::iterator cur, end;   // the enot recoming node has to be adapted for the shorter branch
//        for(cur=node2->disp_vec.begin(), end=node2->disp_vec.end(); cur!=end; ++cur){
//            if((unsigned long)cur->first < (unsigned long)time) continue;
//            node2->disp_vec.erase(cur, end);							// remove the last part
//            break;
//        }
//        assert(node2->disp_vec.back().first<time);
//#ifdef _DEBUG
//		message("\n  ***** The time to the MRCA has changed to %i (node %i)\n", time, coal.get_mrca()->ID_Node+1);
//#endif
//    }
//    
//	// set the ancestor for both nodes
//	node1->ancestor = temp1;
//	node2->ancestor = temp1;
//    
//	///////////////////////////////////////////////////////////////////////////
//	// 3. node2 is also a recombining branch
//	if(node2->recombTime == RECOMB)
//    {
//        assert(node2->recombTime == RECOMB);
//        node2->recombTime = my_NAN;
//        
//		temp2 = recombining_branch(node2, time, coal);
//		if(temp2) {		// if present delete it: we need currently just one new node
//			delete temp2;
//			temp2=NULL;
//		}
//        
//		///////////////////////////////////////////////////////////////////////////
//		temp1->get_currentDeme(time);   // check that the current deme is present
//	}
//    
//    // if there is a recombination event on node2 at a later time add it to the new node
//    else if(node2->recombTime!=my_NAN)
//    {
//        assert(node2->recombTime > time);
//        temp1->recombTime = node2->recombTime;
//        node2->recombTime = my_NAN;
//        
//        // adapt the vector _vecNodes (replace the obsolete node by other)
//		multimap<unsigned long, TTNode*>::iterator cur = coal._vecNodes_pos;
//        for(; cur!= coal._vecNodes_end; ++cur){
//            if(cur->second != node2) continue;
//            cur->second = temp1;
//            break;
//        }
//    }
//    
//	ChainNodeList.push_back(temp1);
//    
//#ifdef _DEBUG
//	cout << " into node " << temp1->ID_Node+1 << flush;
//#endif
//    
//	return temp1;
//}


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
