/** @file lce_regulation.cpp
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


#include "lce_regulation.h"
#include "tselection.h"


/* _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/ */

// ******** LCE_Regulation ********

// ----------------------------------------------------------------------------------------
// LCE_Regulation::init
// ----------------------------------------------------------------------------------------
bool LCE_Regulation::init(Metapop* popPtr) {
	LCE::init(popPtr);
    
	TSelection* pSel = _popPtr->get_pSelection(); // if selection acts
	if (pSel && ((_age == OFFSx && pSel->get_selection_position() == 2)
                 // does selection acts at offsrping regulation
                 || (_age == ADLTx && pSel->get_selection_position() == 3))) {
		// does selection acts at adult regulation
        
		switch(pSel->get_selection_level()) {
            case 0:
                regulation = &LCE_Regulation::regulation_fitness_patch;
                break; // patch level
            case 1:
                regulation = &LCE_Regulation::regulation_fitness_metapop; // metapop level
                _popPtr->set_total_carrying_capacity();
                break;
            case 2:
                regulation = &LCE_Regulation::regulation_fitness_hard;
                break; // hard selection
		}
        
	}
	else if ((unsigned int)this->get_parameter_value
             ("regulation_model_" + _ageStr)) { // neutral regulation
		if (_popPtr->isCoalescence())
			regulation = &LCE_Regulation::regulation_coalescence;
		else
			regulation = &LCE_Regulation::regulation_neutral;
	}
	else { // no regulation
		_paramSet->set_isSet(false);
		return false;
	}
	return true;
}

// ----------------------------------------------------------------------------------------
// LCE_Regulation::execute
// ----------------------------------------------------------------------------------------
void LCE_Regulation::execute() {
#ifdef _DEBUG
	string reg;
	if (regulation == &LCE_Regulation::regulation_fitness_patch) reg = "patch selection";
	else if (regulation == &LCE_Regulation::regulation_fitness_metapop) reg = "metapop selection";
	else if (regulation == &LCE_Regulation::regulation_fitness_hard) reg = "hard selection";
	else reg = "neutral";
	message("  LCE_Regulation (%s, %s)... ", _ageStr.c_str(), reg.c_str());
#endif
    
	(this->*regulation)();
    
#ifdef _DEBUG
	message("done! (Patches: %i; Individuals(F/M): [%i/%i; %i/%i])\n",
            _popPtr->get_nbFullPatch(), _popPtr->size(FEM, OFFSx),
            _popPtr->size(MAL, OFFSx), _popPtr->size(FEM, ADLTx),
            _popPtr->size(MAL, ADLTx));
    assert(_popPtr->individual_container_ok());
#endif
}

// ----------------------------------------------------------------------------------------
// LCE_Regulation::execute
// ----------------------------------------------------------------------------------------
void LCE_Regulation::regulation_neutral() {
	unsigned int N, K;
	Patch* current_patch;
    
#ifdef _DEBUG
	_ex_cnt = _col_cnt = _ph_cnt = 0;
#endif
    
	vector<Patch*>::iterator curPop = _popPtr->get_vFullPatch().begin();
	vector<Patch*>::iterator endPop = _popPtr->get_vFullPatch().end();
	for (; curPop != endPop; ++curPop) {
		current_patch = *curPop;
        
		// females
		K = current_patch->get_KFem();
		N = current_patch->size(FEM, _age);
		if (N > K) {
			if (N < K / 2)
				drawSuccessfullIndividuals(current_patch, K, FEM);
			else
				drawUnSuccessfullIndividuals(current_patch, K, FEM); // if(N>K/2)
		}
        
		// males
		K = current_patch->get_KMal();
		N = current_patch->size(MAL, _age);
		if (N > K) {
			if (N < K / 2)
				drawSuccessfullIndividuals(current_patch, K, MAL);
			else
				drawUnSuccessfullIndividuals(current_patch, K, MAL); // if(N>K/2)
		}
        
#ifdef _DEBUG
		_ex_cnt += current_patch->get_isExtinct();
		_col_cnt += (current_patch->get_isExtinct()
                     ? current_patch->nbKolonisers : 0);
#endif
	} // end for patches
}

// ----------------------------------------------------------------------------------------
// LCE_Regulation::execute
// ----------------------------------------------------------------------------------------
void LCE_Regulation::regulation_coalescence() {
	unsigned int N, K;
	Patch* current_patch;
    
#ifdef _DEBUG
	_ex_cnt = _col_cnt = _ph_cnt = 0;
#endif
    
	vector<Patch*>::iterator curPop = _popPtr->get_vFullPatch().begin();
	vector<Patch*>::iterator endPop = _popPtr->get_vFullPatch().end();
	for (; curPop != endPop; ++curPop) {
		current_patch = *curPop;
        
		// females
		K = current_patch->get_KFem();
		N = current_patch->size(FEM, _age);
		if (N > K)
			current_patch->set_size(FEM, _age, K);
        
#ifdef _DEBUG
		_ex_cnt += current_patch->get_isExtinct();
		_col_cnt += (current_patch->get_isExtinct()
                     ? current_patch->nbKolonisers : 0);
#endif
	} // end for patches
}

// ----------------------------------------------------------------------------------------
// LCE_Regulation::drawSuccessfulIndividuals
// ----------------------------------------------------------------------------------------
/** regulates randomly the size of the pop, by randomly drawing survivors
 * (this needs a new array of the succeeders which will then be exchanged)
 */
void LCE_Regulation::drawSuccessfullIndividuals(Patch* curPatch,
                                                const unsigned int& K, const sex_t& SEX) {
	unsigned int ind, cur_size = curPatch->size(SEX, _age);
    
	vector<Individual*>vec;
	vec.assign(K, NULL);
    
	for (unsigned int i = 0; i < K; ++i, --cur_size) {
		ind = get_pop_ptr()->rand().Uniform(cur_size);
        
		assert(!(curPatch->get_isExtinct() && curPatch->get(SEX, _age,
                                                            ind)->getCurrentPatch()->get_ID() == curPatch->get_ID()));
        
		// this individual survived, congratulations!
		vec[i] = curPatch->get(SEX, _age, ind);
		curPatch->remove(SEX, _age, ind);
	}
    
	// reassign the container
	vector<Individual*>&container = curPatch->get_containers()[SEX][_age];
	delete &container; // delete the remaining individuals
	container = vec; // assign the new array to the container
	curPatch->get_sizes()[SEX][_age] = K;
    
	assert(curPatch->size(SEX, _age) == K);
}

// ----------------------------------------------------------------------------------------
// LCE_Regulation::drawUnSuccessfulIndividuals
// ----------------------------------------------------------------------------------------
/** regulates randomly the pop size, by removing randomly supernumerous individuals */
void LCE_Regulation::drawUnSuccessfullIndividuals(Patch* curPatch,
                                                  const unsigned int& K, const sex_t& SEX) {
	unsigned int ind, nbInd;
    
	// remove randomly supernumerous individuals
	for (nbInd = curPatch->size(SEX, _age); nbInd > K; --nbInd) {
        assert(nbInd == curPatch->size(SEX, _age));
		ind = get_pop_ptr()->rand().Uniform(nbInd);
		curPatch->recycle(SEX, _age, ind);
	}
	assert(curPatch->size(SEX, _age) == K);
}

// ----------------------------------------------------------------------------------------
