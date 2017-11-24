/** @file indfactory.cpp
*
*   Copyright (C) 2006 Frederic Guillaume    <guillaum@zoology.ubc.ca>
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

#include "indfactory.h"
#include "tmetapop.h"

// ----------------------------------------------------------------------------------------
// makePrototype
// ----------------------------------------------------------------------------------------
/** generates the proto traits and thus a prototype individual */
void IndFactory::makePrototype(map< string,TTraitProto* >& TTlist, TMetapop* pMetapop){
  #ifdef _DEBUG
   message("IndFactory::makePrototype\n");
  #endif
	_popPtr = pMetapop;
	_protoGenome->_popPtr = _popPtr;

	//store the traits list for future use:
	_protoTraits = TTlist;

	//then add the traits:
	assert(_TraitsIndex.empty());

	int i = 0;
	map< string, TTraitProto* >::iterator trait;
	for(trait = TTlist.begin(); trait != TTlist.end(); ++trait, ++i) {
		#ifdef _DEBUG
		 message("IndFactory::makePrototype::addTrait: %s\n",trait->first.c_str());
		#endif

		trait->second->init(pMetapop);      					// initialize the prototype traits
		trait->second->set_absolute_index(i);                   // set the absolute trait index
		_TraitsIndex[trait->first] = i;
		_protoTraitsVector.push_back(trait->second);
	}

	_protoGenome->ini_all();

	// generate a prototype individual
	_protoIndividual.clearTraits();                            	// reset the prototype individual
	_protoIndividual.genome.set_protoGenome(_protoGenome);
	_protoIndividual.genome.create_sequence();
	for(i=0, trait = TTlist.begin(); trait != TTlist.end(); ++trait, ++i) {
		_protoIndividual.addTrait(trait->second->hatch(),i);    // create the prototype individual
	}
}

// ----------------------------------------------------------------------------------------
// resetPrototype
// ----------------------------------------------------------------------------------------
/** find the first trati fo the given type */
TTraitProto* IndFactory::getFirstPrototype(string type)
{
    map<string,TTraitProto* >::iterator pos=_protoTraits.find(type);        // search for e.g. "quanti"
    if(pos!=_protoTraits.end()) return pos->second;
    pos=_protoTraits.find(type+"_1");                                       // search for e.g. "quanti_1"
    if(pos!=_protoTraits.end()) return pos->second;
    return NULL;
}

// ----------------------------------------------------------------------------------------
// resetPrototype
// ----------------------------------------------------------------------------------------
void IndFactory::resetPrototype()
{
  #ifdef _DEBUG
   message("IndFactory::resetPrototype\n");
  #endif

	//first, reset the ID counters of the individuals in the patch container:
	vector<Patch*>::iterator curPop = _popPtr->get_vPatch().begin();
	vector<Patch*>::iterator endPop = _popPtr->get_vPatch().end();
	for(; curPop!=endPop; ++curPop){
		(*curPop)->reset_ID_individual();
	}

	vector<TTraitProto* >::iterator trait;
	for(trait = _protoTraitsVector.begin(); trait != _protoTraitsVector.end(); ++trait) {
		(*trait)->reset();
	}

	_protoGenome->ini_genetic_map(); 	// check if the genetic map has to be reset
}

// ----------------------------------------------------------------------------------------
/** returns a vector containing the absolute indexes (across all type of traits)
 * of the linked traits. The position in the vector represents the relative index
 * of the given type of trait.
 */
vector<int> IndFactory::getTraitIndex (string type)
{
    vector<int> index;
    map< string, int >::iterator cur_trait = _TraitsIndex.begin();
    for(; cur_trait != _TraitsIndex.end(); ++cur_trait){
        if(cur_trait->first.find(type) != string::npos){
            index.push_back(cur_trait->second);
        }
    }
    return index;
}

// ----------------------------------------------------------------------------------------
// makeNewIndividual
// ----------------------------------------------------------------------------------------
TIndividual*
IndFactory::makeNewIndividual(TIndividual* mother, TIndividual* father, sex_t sex, Patch* homepatch)
{
  TIndividual* newind;

  if(RecyclingPOOL.empty()) {       //create new Individual
		try{
			newind = _protoIndividual.clone();
		}
		catch(...){
    	error("Out of memory!\n");
		}
		newind->init();                 //allocate memory for the traits' sequences:
	}
  else {                            //recycle an Individual from the POOL
		newind = RecyclingPOOL.back();  //get the individual from the container
		RecyclingPOOL.pop_back();       //remove its place in the vector
		newind->reset();
		newind->reset_counters();       //to reset the matings and fecundity counters
	}

  newind->setSex(sex);
  newind->setNatalPatch(homepatch);
  newind->setID(homepatch->get_next_IDofIndividual());
  newind->setCurrentPatch(homepatch);
  newind->setMother(mother);
	newind->setFather(father);
  return newind;
}

// ----------------------------------------------------------------------------------------
// copyIndividual
// ----------------------------------------------------------------------------------------
TIndividual*
IndFactory::copyIndividual(TIndividual* oldInd)
{
    TIndividual* newind;
    
    if(RecyclingPOOL.empty()) {       //create new Individual
        try{
            newind = new TIndividual(_protoIndividual);
        }
        catch(...){
            error("Out of memory!\n");
        }
        newind->init();                 //allocate memory for the traits' sequences:
    }
    else {                            //recycle an Individual from the POOL
        newind = RecyclingPOOL.back();  //get the individual from the container
        RecyclingPOOL.pop_back();       //remove its place in the vector
    }
    *newind = *oldInd;
    
    if(*newind==*oldInd) cout << "same indvidual" << flush;
    return newind;
}

// ----------------------------------------------------------------------------------------
// makeOffsprg
// ----------------------------------------------------------------------------------------
TIndividual* IndFactory::makeOffsprg(TIndividual* mother, TIndividual* father, sex_t sex, Patch* homepatch)
{
  TIndividual* NewOffsprg = makeNewIndividual(mother,father,sex,homepatch);

  NewOffsprg->create(mother,father);  // inherit and mutate the traits:
	if(father){                         // normal inheritance
		NewOffsprg->setIsSelfed(mother == father);
		bool cat = (mother->getNatalPatch() == father->getNatalPatch());
		mother->DidHaveABaby(cat);
		father->DidHaveABaby(cat);
	}
	else{                               // cloning (cloning is considered as selfing)
		NewOffsprg->setIsSelfed(true);
		mother->DidHaveABaby(true);
	}
	homepatch->add(sex, OFFSx, NewOffsprg);   // add individual to the patch
  return NewOffsprg;
}

