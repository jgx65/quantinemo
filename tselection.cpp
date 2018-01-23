/** @file tselection.cpp
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

//---------------------------------------------------------------------------
#include "tselection.h"
#include "tselectiontype.h"
#include "metapop.h"
#include "stathandler.cpp"
#include "ttquanti.h"


// -----------------------------------------------------------------------------
// destructor
// -----------------------------------------------------------------------------
TSelection::~TSelection ( )
{
	if(_selTrait){
		for(unsigned int i=0; i<_vTraitsSize; ++i){
			if(_selTrait[i]) delete _selTrait[i];
		}
		delete[] _selTrait;
	}

	if(_aPheno)   delete[] _aPheno;
    if(_selTrait_fitnessDependent) delete[] _selTrait_fitnessDependent;
    
    _fem_sex_allocation_value = my_NAN;
    

	for(int i=0; i<2; ++i){
		if(_fit[i])   delete _fit[i];
	}
}

//-----------------------------------------------------------------------------
// set_phenotype
//-----------------------------------------------------------------------------
void
TSelection::set_phenotype(Individual* ind)
{
	for(unsigned int t = 0; t<_vTraitsSize; ++t){
		(_selTrait[t]->*(_selTrait[t]->func_ptr_set_phenotype))(ind);
	}
}


//-----------------------------------------------------------------------------
// set_phenotype
//-----------------------------------------------------------------------------
void
TSelection::set_phenotype(Individual* ind, unsigned int qtraitID)
{
	(_selTrait[qtraitID]->*(_selTrait[qtraitID]->func_ptr_set_phenotype))(ind);
}

//-----------------------------------------------------------------------------
// set_phenotype
//-----------------------------------------------------------------------------
/** set the phenotype of all individuals of the given age */
void
TSelection::set_phenotype_of_all_individuals(age_idx AGE)
{
	vector<Individual*>::iterator curInd, endInd;
	vector<Patch*>::iterator curPop = _popPtr->get_vFullPatch().begin();
	vector<Patch*>::iterator endPop = _popPtr->get_vFullPatch().end();
	for(;curPop!=endPop; ++curPop){                                 // for each patch
		// females
		vector<Individual*>& curFem = (*curPop)->get_sampled_inds(FEM, AGE);
		for(curInd=curFem.begin(), endInd=curFem.end(); curInd!=endInd; ++curInd){
			set_phenotype(*curInd);  // set the phenotype
		}

		// males
		vector<Individual*>& curMal = (*curPop)->get_sampled_inds(MAL, AGE);
		for(curInd=curMal.begin(), endInd=curMal.end(); curInd!=endInd; ++curInd){
			set_phenotype(*curInd);  // set the phenotype
		}
	}
}

//-----------------------------------------------------------------------------
// LCE_Breed_fitness::
//-----------------------------------------------------------------------------
/** computes the fitness of the individuals of the patch and age. returns an array of the fitness of all individuals of the specified patch, sex and age
 * First the phenotypes will be computed and set.
 * patch:   individuals of this patch will be considered
 * age:     only this age class will be considered
 * sort:    if the selection is based on something like fittest or less fittest,
 *          it is necessary to sort the fitness array. Therefore a two-dimensional
 *          array has to be passed which informs about the sort
 *          (0: no (default); 1: upwards; 2: downwards)?
 */
void
TSelection::set_fitness(Patch* patch, age_idx age, int* sort, int subset)
{
    // does frequency dependend selection act?
    if(_selTrait_fitnessDependentSize){
        for(unsigned int t=0; t<_selTrait_fitnessDependentSize ; ++t){
            TTQuantiSH* pStats = dynamic_cast< TTQuantiSH*>(_selTrait_fitnessDependent[t]->_stats);
            pStats->get_locusGenotypeFreqs_ofPatch_allInds(patch, age, t, _selTrait_fitnessDependent[t]->get_locusFreqs());
        }
    }
    
    // compute the fitnesses
    set_fitness(patch, FEM, age);
    set_fitness(patch, MAL, age);
    
    // add female sex allocation if needed
    (this->*female_sex_allocation_func_ptr)(patch, age, sort, subset);
    
    // sort and make cumulative
	if(sort){
        _fit[FEM]->sort(sort[FEM], subset);
        _fit[MAL]->sort(sort[MAL], subset);
    }
    else{
        _fit[FEM]->sort(10, subset);
        _fit[MAL]->sort(10, subset);
    }
}

//-----------------------------------------------------------------------------
// set_fitness::
//-----------------------------------------------------------------------------
/** returns an array of the fitness of all individuals of the specified patch, sex and age
 * First the phenotypes will be computed and set.
 * patch:   individuals of this patch will be considered
 * sex:     only this sex will be considered
 * age:     only this age class will be considered
 */
void
TSelection::set_fitness(Patch* patch, sex_t sex, age_idx age)
{
	unsigned int i, newSize;
	Individual* ind;

	newSize = patch->size(sex, age);

	// resize TPatchFitness to correct number of individuals
	if(!_fit[sex]) _fit[sex] = new TPatchFitness(_popPtr, newSize);
	else _fit[sex]->resize(newSize);

	// if the patch is empty for this sex, return
	if(!newSize) return;

	// get the selection pressure of the patch
	for(i=0; i<_vTraitsSize; ++i){
		(_selTrait[i]->*(_selTrait[i]->_get_selection_pressure_tot_func_ptr))(patch, sex);
	}


	// get the current arrays
	double*      aFit = _fit[sex]->_aFit;
	Individual** aInd = _fit[sex]->_aInd;

	// for each individual
	for(i=0; i<newSize; ++i){
		aInd[i] = ind = patch->get(sex, age, i);            // get the individuals
		set_phenotype(ind);              					// set the phenotypes of this individual at each trait
		aFit[i] = _get_fitness_multiplicative(ind);         // compute overall fitness
		ind->setFitness(aFit[i]);                           // set the fitness to the individual
	}
}

//-----------------------------------------------------------------------------
// sort_fitness
//-----------------------------------------------------------------------------
/** sort the array if necessary and make it cumulative */
void
TSelection::sort_fitness(sex_t SEX, int how, int subset)
{
	assert(_fit[SEX]->_sort == 10);
	_fit[SEX]->sort(how, subset);
}

//-----------------------------------------------------------------------------
// LCE_Breed_fitness::
//-----------------------------------------------------------------------------
/** fitness are multiplicative */
double
TSelection::_get_fitness_multiplicative(Individual* ind)
{
	double cur_w, w = 1;
	for(unsigned int t=0; t<_vTraitsSize; ++t){
		assert(_selTrait[t]->get_fitness() != my_NAN);
		cur_w = _selTrait[t]->get_fitness();
		ind->setTraitFitness(t, cur_w);
		w *= cur_w;
	}
	return w;
}

//-----------------------------------------------------------------------------
// TSelection::female_sex_allocation
//-----------------------------------------------------------------------------
/** if the female sex allocation of hermaphrodites may evolve 
 * copy the individualarray to the males (all individuals are present in both "sexes")
 * get the female allocation trait and its current value for each indvidual
 * multiple the sex allocation to male and female overall fitness
 * re-sort the fitness arrays
 */
void
TSelection::female_sex_allocation_fix(Patch* patch, age_idx age, int* sort, int subset)
{
    assert(_fem_sex_allocation_value>=0);
    assert(_fem_sex_allocation_value<=1);
    
    // get the "females" (they are set)
    double*      femFit = _fit[FEM]->_aFit;     // all females without sex allocation
    Individual** femInd = _fit[FEM]->_aInd;     // all females without sex allocation
    unsigned int nbFem = patch->size(FEM, age);
    
    // get the "males" (they are not yet set)
    if(!_fit[MAL]) _fit[MAL] = new TPatchFitness(_popPtr, nbFem);
    else _fit[MAL]->resize(nbFem);
    double*      malFit = _fit[MAL]->_aFit;
    Individual** malInd = _fit[MAL]->_aInd;
    
    // for each individualget the overall fitness and correct it for the sex allocation
    for(unsigned int i=0; i<nbFem; ++i){
        malInd[i] = femInd[i];                                      // female and male array are at the moment identical
        malFit[i] = (1-_fem_sex_allocation_value)*femFit[i];        // compute male allocation proportion and multiply it with the overall fitness
        femFit[i] *= _fem_sex_allocation_value;                     // mutiply overall fitness with female sex allocation
    }
    
    // re-sort arrays
    _fit[FEM]->sort(sort[FEM], subset);                        	// sort
    _fit[MAL]->sort(sort[MAL], subset);                        	// sort
}

//-----------------------------------------------------------------------------
// TSelection::female_sex_allocation
//-----------------------------------------------------------------------------
/** if the female sex allocation of hermaphrodites may evolve
 * copy the individualarray to the males (all individuals are present in both "sexes")
 * get the female allocation trait and its current value for each indvidual
 * multiple the sex allocation to male and female overall fitness
 * re-sort the fitness arrays
 */
void
TSelection::female_sex_allocation_G(Patch* patch, age_idx age, int* sort, int subset)
{
    assert(_fem_sex_allocation_value != my_NAN);
    
    // get the "females" (they are set)
    double*      femFit = _fit[FEM]->_aFit;     // all females without sex allocation
    Individual** femInd = _fit[FEM]->_aInd;     // all females without sex allocation
    unsigned int nbFem = patch->size(FEM, age);
    
    // get the "males" (they are not yet set)
    if(!_fit[MAL]) _fit[MAL] = new TPatchFitness(_popPtr, nbFem);
    else _fit[MAL]->resize(nbFem);
    double*      malFit = _fit[MAL]->_aFit;
    Individual** malInd = _fit[MAL]->_aInd;
    
    // for each individualget the overall fitness and correct it for the sex allocation
    double fem_alloc;
    for(unsigned int i=0; i<nbFem; ++i){
        // get the female allocation proportion and check for its validity
        fem_alloc = femInd[i]->getTraitGenotype(_fem_sex_allocation_value);
        if(fem_alloc < 0) error("female_sex_allocation is below 0. It should be between 0 and 1.");
        if(fem_alloc > 1) error("female_sex_allocation is above 1. It should be between 0 and 1.");
        
        malInd[i] = femInd[i];                      // female and male array are at the moment identical
        malFit[i] = (1-fem_alloc)*femFit[i];        // compute male allocation proportion and multiply it with the overall fitness
        femFit[i] *= fem_alloc;                     // mutiply overall fitness with female sex allocation
    }
    
    // re-sort arrays
    _fit[FEM]->sort(sort[FEM], subset);                        	// sort
    _fit[MAL]->sort(sort[MAL], subset);                        	// sort
}

//-----------------------------------------------------------------------------
// TSelection::female_sex_allocation
//-----------------------------------------------------------------------------
/** if the female sex allocation of hermaphrodites may evolve
 * copy the individualarray to the males (all individuals are present in both "sexes")
 * get the female allocation trait and its current value for each indvidual
 * multiple the sex allocation to male and female overall fitness
 * re-sort the fitness arrays
 */
void
TSelection::female_sex_allocation_equation(Patch* patch, age_idx age, int* sort, int subset)
{
    assert(_female_sex_allocation);
    
    // get the "females" (they are set)
    double*      femFit = _fit[FEM]->_aFit;     // all females without sex allocation
    Individual** femInd = _fit[FEM]->_aInd;     // all females without sex allocation
    unsigned int nbFem = patch->size(FEM, age);
    
    // get the "males" (they are not yet set)
    if(!_fit[MAL]) _fit[MAL] = new TPatchFitness(_popPtr, nbFem);
    else _fit[MAL]->resize(nbFem);
    double*      malFit = _fit[MAL]->_aFit;
    Individual** malInd = _fit[MAL]->_aInd;
    
    // for each individualget the overall fitness and correct it for the sex allocation
    double fem_alloc;
    for(unsigned int i=0; i<nbFem; ++i){
        // get the female allocation proportion and check for its validity
        fem_alloc = _female_sex_allocation->getValue(femInd[i]);
        if(fem_alloc < 0) error("female_sex_allocation is below 0. It should be between 0 and 1.");
        if(fem_alloc > 1) error("female_sex_allocation is above 1. It should be between 0 and 1.");
        
        malInd[i] = femInd[i];                      // female and male array are at the moment identical
        malFit[i] = (1-fem_alloc)*femFit[i];        // compute male allocation proportion and multiply it with the overall fitness
        femFit[i] = fem_alloc*femFit[i];            // mutiply overall fitness with female sex allocation
    }
    
    // re-sort arrays
    _fit[FEM]->sort(sort[FEM], subset);                        	// sort
    _fit[MAL]->sort(sort[MAL], subset);                        	// sort
}

//-----------------------------------------------------------------------------
// TSelection::female_sex_allocation
//-----------------------------------------------------------------------------
/** if the female sex allocation of hermaphrodites may evolve
 * copy the individualarray to the males (all individuals are present in both "sexes")
 * get the female allocation trait and its current value for each indvidual
 * multiple the sex allocation to male and female overall fitness
 * re-sort the fitness arrays
 */
void
TSelection::female_sex_allocation_Z(Patch* patch, age_idx age, int* sort, int subset)
{
    assert(_fem_sex_allocation_value != my_NAN);
    
    // get the "females" (they are set)
    double*      femFit = _fit[FEM]->_aFit;     // all females without sex allocation
    Individual** femInd = _fit[FEM]->_aInd;     // all females without sex allocation
    unsigned int nbFem = patch->size(FEM, age);
    
    // get the "males" (they are not yet set)
    if(!_fit[MAL]) _fit[MAL] = new TPatchFitness(_popPtr, nbFem);
    else _fit[MAL]->resize(nbFem);
    double*      malFit = _fit[MAL]->_aFit;
    Individual** malInd = _fit[MAL]->_aInd;
    
    // for each individualget the overall fitness and correct it for the sex allocation
    double fem_alloc;
    for(unsigned int i=0; i<nbFem; ++i){
        // get the female allocation proportion and check for its validity
        fem_alloc = femInd[i]->getTraitPhenotype(_fem_sex_allocation_value);
        if(fem_alloc < 0) error("female_sex_allocation is below 0. It should be between 0 and 1.");
        if(fem_alloc > 1) error("female_sex_allocation is above 1. It should be between 0 and 1.");
        
        malInd[i] = femInd[i];                      // female and male array are at the moment identical
        malFit[i] = (1-fem_alloc)*femFit[i];        // compute male allocation proportion and multiply it with the overall fitness
        femFit[i] *= fem_alloc;                     // mutiply overall fitness with female sex allocation
    }
    
    // re-sort arrays
    _fit[FEM]->sort(sort[FEM], subset);                        	// sort
    _fit[MAL]->sort(sort[MAL], subset);                        	// sort
}

//-----------------------------------------------------------------------------
// TSelection::init
//-----------------------------------------------------------------------------
/** the function returns true if selection is used and false if not */
bool
TSelection::init(Metapop* popPtr)
{
#ifdef _DEBUG
	message("TSelection::init\n");
#endif
    
	// initialize parameters /////////////////////////////////////////////////////
	for(int i=0; i<2; ++i){
		_fit[i] = NULL;
	}
	_selTrait = NULL;
	_aPheno   = NULL;
	_popPtr = popPtr;
    _selTrait_fitnessDependent=NULL;
    _fem_sex_allocation_value=my_NAN;

	// check if selection is used ////////////////////////////////////////////////
	_selection_position   = (unsigned int)_popPtr->get_paramset()->getValue("selection_position");
	if(_selection_position == 4) return false;                 // 4: no selection at all

	// check if quantitative traits are specified
	_vTraits = _popPtr->getTraitIndex("quanti");  // get the indexes to the quanti traits
	_vTraitsSize = (unsigned int)_vTraits.size();
	if(!_vTraitsSize) return false;               // no quantitative traits => no selection at all

	// what specifies the mean fitness of a population ///////////////////////////
	switch((unsigned int)_popPtr->get_paramset()->getValue("patch_mean_fitness")){
		case 0: getMeanFitness_func_ptr = &TSelection::getMeanFitnessBoth; break; // fem + mal (default)
		case 1: getMeanFitness_func_ptr = &TSelection::getMeanFitnessFem;  break; // just fem
		case 2: getMeanFitness_func_ptr = &TSelection::getMeanFitnessMal;  break; // just mal
	}

	// selection acts, prepare the object ////////////////////////////////////////
	// get the linked traits
	_aPheno      = new double[_vTraitsSize];
	_nbPops      = _popPtr->getPatchNbr();

	_optima_sd        = sqrt(_popPtr->get_paramset()->getValue("patch_stab_sel_optima_var"));
	_intensity_sd     = sqrt(_popPtr->get_paramset()->getValue("patch_stab_sel_intensity_var"));

	_min_sd           = sqrt(_popPtr->get_paramset()->getValue("patch_dir_sel_min_var"));
	_max_sd           = sqrt(_popPtr->get_paramset()->getValue("patch_dir_sel_max_var"));
	_growth_rate_sd   = sqrt(_popPtr->get_paramset()->getValue("patch_dir_sel_growth_rate_var"));
	_max_growth_sd    = sqrt(_popPtr->get_paramset()->getValue("patch_dir_sel_max_growth_var"));
	_symmetry_sd      = sqrt(_popPtr->get_paramset()->getValue("patch_dir_sel_symmetry_var"));

	_selCoef_sd       = sqrt(_popPtr->get_paramset()->getValue("patch_coef_sel_var"));

	// create the TSelectionTrait objects
	_selTrait = new TSelectionTrait*[_vTraitsSize];
	unsigned int selKind = 0;                 // bit: 1: neutral, 2: stabilizing, 4: directional
	for(unsigned int t = 0; t < _vTraitsSize; ++t) {
		// find the kind of selection of the trait
		assert((dynamic_cast <TTQuantiProto*> (&_popPtr->getTraitPrototype(_vTraits[t]))));
		switch((dynamic_cast <TTQuantiProto*> (&_popPtr->getTraitPrototype(_vTraits[t])))->get_selection_model()){
			case 0: // neutral:     1. digit
				_selTrait[t] = new TSelectionNeutral(this, t);
				selKind |= 1;
				break;
			case 1: // stabilizing: 2. digit
				_selTrait[t] = new TSelectionStabilizing(this, t);
				selKind |= 2;
				break;
			case 2: // directional: 3. digit
				_selTrait[t] = new TSelectionDirectional(this, t);
				selKind |= 4;
				break;
			case 3: // fitness landscape
				_selTrait[t] = new TSelectionFitnessLandscape(this, t);
				selKind |= 8;
				break;
			case 4: // selection coefficient:     1. digit
				_selTrait[t] = new TSelectionSelectionCoefficient(this, t);
				selKind |= 16;
				break;
			}
	}


	// set the nb linked traits in the pops
	for(int p=0; p< _nbPops; ++p){
		_popPtr->get_vPatch(p)->set_LinkedTraits(_vTraitsSize, _selTrait, (bool)_popPtr->get_sexInitRatio());
	}

	// environmental variance
	set_ve_mean();
	set_ve_var();

	// what must be read?
	if(selKind & 2){      // stabilizing selection
		_popPtr->set_patch_parameter(_vTraitsSize, "patch_stab_sel_optima",    "optima",    &Patch::set_localOptima);    // default 0
		_popPtr->set_patch_parameter(_vTraitsSize, "patch_stab_sel_intensity", "intensity", &Patch::set_localIntensity); // default 1
	}
	if(selKind & 4){      // directional selection
		_popPtr->set_patch_parameter(_vTraitsSize, "patch_dir_sel_growth_rate", "growth_rate", &Patch::set_localGrowthRate);  // default 1
		_popPtr->set_patch_parameter(_vTraitsSize, "patch_dir_sel_max_growth",  "max_growth",  &Patch::set_localMaxGrowth);   // default 1
		_popPtr->set_patch_parameter(_vTraitsSize, "patch_dir_sel_symmetry",    "symmetry",    &Patch::set_localSymmetry);    // default 0.5
		_popPtr->set_patch_parameter(_vTraitsSize, "patch_dir_sel_min",         "min",         &Patch::set_localMin);         // default 0
		_popPtr->set_patch_parameter(_vTraitsSize, "patch_dir_sel_max",         "max",         &Patch::set_localMax);         // default 1
	}
	if(selKind & 8){      // fitness landscape (phenotype has to be set after fitness due to sorting)
		_popPtr->set_patch_parameter_array(_vTraitsSize, "patch_fitness_landscape",   "fitness_landscape", &Patch::set_fitnessLandscape_fitness);
		_popPtr->set_patch_parameter_array(_vTraitsSize, "patch_phenotype_landscape", "phenotype_landscape", &Patch::set_fitnessLandscape_phenotype);
	}
	if(selKind & 16){      // selection coefficient
		_popPtr->set_patch_parameter(_vTraitsSize, "patch_coef_sel_AA",  "selection_coefficient_AA",&Patch::set_localSelCoefAA);  // default 1
		_popPtr->set_patch_parameter(_vTraitsSize, "patch_coef_sel",    "selection_coefficient",    &Patch::set_localSelCoef);    // default 0
	}

	_selection_level      = (unsigned int)_popPtr->get_paramset()->getValue("selection_level");
	set_selection_level_coef();
    set_frequency_dependend_selection();
	return true;
}

//-----------------------------------------------------------------------------
// TSelection::init
//-----------------------------------------------------------------------------
/** the function returns true if selection is used and false if not */
bool
TSelection::init2(Metapop* popPtr)
{
#ifdef _DEBUG
	message("TSelection::init2\n");
#endif

    // initialize parameters /////////////////////////////////////////////////////
	for(int i=0; i<2; ++i){
		_fit[i] = NULL;
	}
	_selTrait = NULL;
	_aPheno   = NULL;
	_popPtr = popPtr;
    _selTrait_fitnessDependent=NULL;
    _fem_sex_allocation_value=my_NAN;

	// check if selection is used ////////////////////////////////////////////////
	_selection_position   = (unsigned int)_popPtr->get_paramset()->getValue("selection_position");
	if(_selection_position == 4) return false;                 // 4: no selection at all
    
    // check if frequency dependend selection is used //////////////////////////
    

	// check if quantitative traits are specified
	_vTraits = _popPtr->getTraitIndex("quanti");  // get the indexes to the quanti traits
	_vTraitsSize = (unsigned int)_vTraits.size();
	if(!_vTraitsSize) return false;               // no quantitative traits => no selection at all

	// what specifies the mean fitness of a population ///////////////////////////
	switch((unsigned int)_popPtr->get_paramset()->getValue("patch_mean_fitness")){
		case 0: getMeanFitness_func_ptr = &TSelection::getMeanFitnessBoth; break; // fem + mal (default)
		case 1: getMeanFitness_func_ptr = &TSelection::getMeanFitnessFem;  break; // just fem
		case 2: getMeanFitness_func_ptr = &TSelection::getMeanFitnessMal;  break; // just mal
	}


	// selection acts, prepare the object ////////////////////////////////////////
	// get the linked traits
	_aPheno      = new double[_vTraitsSize];
	_nbPops      = _popPtr->getPatchNbr();

	_optima_sd        = sqrt(_popPtr->get_paramset()->getValue("patch_stab_sel_optima_var"));
	_intensity_sd     = sqrt(_popPtr->get_paramset()->getValue("patch_stab_sel_intensity_var"));

	_min_sd           = sqrt(_popPtr->get_paramset()->getValue("patch_dir_sel_min_var"));
	_max_sd           = sqrt(_popPtr->get_paramset()->getValue("patch_dir_sel_max_var"));
	_growth_rate_sd   = sqrt(_popPtr->get_paramset()->getValue("patch_dir_sel_growth_rate_var"));
	_max_growth_sd    = sqrt(_popPtr->get_paramset()->getValue("patch_dir_sel_max_growth_var"));
	_symmetry_sd      = sqrt(_popPtr->get_paramset()->getValue("patch_dir_sel_symmetry_var"));

	_selCoef_sd       = sqrt(_popPtr->get_paramset()->getValue("patch_coef_sel_var"));

	// create the TSelectionTrait objects
	_selTrait = new TSelectionTrait*[_vTraitsSize];
	string trait;
	TTQuantiProto* pQuanti;
	for(unsigned int t = 0; t < _vTraitsSize; ++t) {
		assert((dynamic_cast <TTQuantiProto*> (&_popPtr->getTraitPrototype(_vTraits[t]))));
		pQuanti = (dynamic_cast <TTQuantiProto*> (&_popPtr->getTraitPrototype(_vTraits[t])));
		switch(pQuanti->get_selection_model()){          // find the kind of selection of the trait
			case 0: _selTrait[t] = new TSelectionNeutral(this, t);         	break;  // neutral:           1. digit
			case 1: _selTrait[t] = new TSelectionStabilizing(this, t);      break;  // stabilizing:       2. digits
			case 2: _selTrait[t] = new TSelectionDirectional(this, t);      break;  // directional:       3. digits
			case 3: _selTrait[t] = new TSelectionFitnessLandscape(this, t); break;  // fitness landscape  x. digits
			case 4: _selTrait[t] = new TSelectionSelectionCoefficient(this, t); 	break;  // selCoef:     1. digit
		}
	}

	// set the linked traits in the patches
	for(int p=0; p< _nbPops; ++p){
		_popPtr->get_vPatch(p)->set_LinkedTraits(_vTraitsSize, _selTrait, (bool)_popPtr->get_sexInitRatio());
	}

	// set the selection pressure (only now since several traits may be under the same selection preessure, and this step has to be executed only onces)
	for(unsigned int t = 0; t < _vTraitsSize; ++t) {
		assert((dynamic_cast <TTQuantiProto*> (&_popPtr->getTraitPrototype(_vTraits[t]))));
		pQuanti = (dynamic_cast <TTQuantiProto*> (&_popPtr->getTraitPrototype(_vTraits[t])));
		trait = pQuanti->get_trait_indexStr_();

		switch(pQuanti->get_selection_model()){          // find the kind of selection of the trait
			case 0: // neutral:     1. digit
				break;
			case 1: // stabilizing: 2. digit
				_popPtr->set_patch_parameter_ofTrait(pQuanti, t, trait, "quanti_stab_sel_optima",    "optima",    &Patch::set_localOptima);    // default 0
				_popPtr->set_patch_parameter_ofTrait(pQuanti, t, trait, "quanti_stab_sel_intensity", "intensity", &Patch::set_localIntensity); // default 1
				break;
			case 2: // directional: 3. digit
				_popPtr->set_patch_parameter_ofTrait(pQuanti, t, trait, "quanti_dir_sel_growth_rate", "growth_rate", &Patch::set_localGrowthRate);  // default 1
				_popPtr->set_patch_parameter_ofTrait(pQuanti, t, trait, "quanti_dir_sel_max_growth",  "max_growth",  &Patch::set_localMaxGrowth);   // default 1
				_popPtr->set_patch_parameter_ofTrait(pQuanti, t, trait, "quanti_dir_sel_symmetry",    "symmetry",    &Patch::set_localSymmetry);    // default 0.5
				_popPtr->set_patch_parameter_ofTrait(pQuanti, t, trait, "quanti_dir_sel_min",         "min",         &Patch::set_localMin);         // default 0
				_popPtr->set_patch_parameter_ofTrait(pQuanti, t, trait, "quanti_dir_sel_max",         "max",         &Patch::set_localMax);         // default 1
				break;
			case 3: // fitness landscape
				_popPtr->set_patch_parameter_array_ofTrait(pQuanti, t, trait, "quanti_fitness_landscape",   "fitness_landscape", &Patch::set_fitnessLandscape_fitness);
				_popPtr->set_patch_parameter_array_ofTrait(pQuanti, t, trait, "quanti_phenotype_landscape", "phenotype_landscape", &Patch::set_fitnessLandscape_phenotype);
				break;
			case 4: // selection coefficient: 1. digit
				_popPtr->set_patch_parameter_ofTrait(pQuanti, t, trait, "quanti_coef_sel_AA", "selection_coefficient_AA", &Patch::set_localSelCoefAA);  // default 1
				_popPtr->set_patch_parameter_ofTrait(pQuanti, t, trait, "quanti_coef_sel",    "selection_coefficient",    &Patch::set_localSelCoef);    // default 0
				break;
		}
	}
	// environmental variance
	set_ve_mean();
	set_ve_var();

	_selection_level      = (unsigned int)_popPtr->get_paramset()->getValue("selection_level");
	set_selection_level_coef();
    set_frequency_dependend_selection();
	return true;
}

//-----------------------------------------------------------------------------
// TSelection::set_frequency_dependend_selection
//-----------------------------------------------------------------------------
/** this funciton initializes the frequency depended selection 
 */
void
TSelection::set_frequency_dependend_selection()
{
	string trait;
	TTQuantiProto* pQuanti;
    vector<TTQuantiProto*> tempTrait;
	for(unsigned int t = 0; t < _vTraitsSize; ++t) {
		assert((dynamic_cast <TTQuantiProto*> (&_popPtr->getTraitPrototype(_vTraits[t]))));
		pQuanti = (dynamic_cast <TTQuantiProto*> (&_popPtr->getTraitPrototype(_vTraits[t])));
        if(pQuanti->get_fitnessFactor_freqDepend()) tempTrait.push_back(pQuanti);
	}
    _selTrait_fitnessDependent=NULL;
    ARRAY::vector2array(tempTrait, _selTrait_fitnessDependent, _selTrait_fitnessDependentSize);
}

//-----------------------------------------------------------------------------
// TSelection::set_fem_sex_allocation
//-----------------------------------------------------------------------------
/** this funciton initializes the sex allocation settings
 */
void
TSelection::set_fem_sex_allocation(unsigned int model, double value)
{
    switch(model){
        default:
        case 1: female_sex_allocation_func_ptr = &TSelection::female_sex_allocation_none; break;
        case 2: female_sex_allocation_func_ptr = &TSelection::female_sex_allocation_fix; break;
        case 3: female_sex_allocation_func_ptr = &TSelection::female_sex_allocation_G; break;
        case 4: female_sex_allocation_func_ptr = &TSelection::female_sex_allocation_Z; break;
    }
    _fem_sex_allocation_value=value;
}

// ----------------------------------------------------------------------------------------
// set_selection_level_coef
// ----------------------------------------------------------------------------------------
// gets the selection_level_coeficient and sets according to it the getNbOff_SoftHardSelection function pointer */
void
TSelection::set_selection_level_coef()
{
	_selection_level_coef = _popPtr->get_paramset()->getValue("selection_level_coef");
	if(_selection_level==0) _selection_level_coef = 0;    // this is soft selection
	if(_selection_level_coef == 0)      get_SoftHardSelection_func_ptr = &TSelection::getSoft;
	else if(_selection_level_coef == 1) get_SoftHardSelection_func_ptr = &TSelection::getHard;
	else                                get_SoftHardSelection_func_ptr = &TSelection::getSoftHard;
}

// ----------------------------------------------------------------------------------------
// reset_selectionTypes
// ----------------------------------------------------------------------------------------
// reset the selection types, i.e. needed when the selection settings vary between generations
void
TSelection::reset_selectionTypes()
{
	_optima_sd        = sqrt(_popPtr->get_paramset()->getValue("patch_stab_sel_optima_var"));
	_intensity_sd     = sqrt(_popPtr->get_paramset()->getValue("patch_stab_sel_intensity_var"));

	_min_sd           = sqrt(_popPtr->get_paramset()->getValue("patch_dir_sel_min_var"));
	_max_sd           = sqrt(_popPtr->get_paramset()->getValue("patch_dir_sel_max_var"));
	_growth_rate_sd   = sqrt(_popPtr->get_paramset()->getValue("patch_dir_sel_growth_rate_var"));
	_max_growth_sd    = sqrt(_popPtr->get_paramset()->getValue("patch_dir_sel_max_growth_var"));
	_symmetry_sd      = sqrt(_popPtr->get_paramset()->getValue("patch_dir_sel_symmetry_var"));

	_selCoef_sd       = sqrt(_popPtr->get_paramset()->getValue("patch_coef_sel_var"));

	for(unsigned int t = 0; t < _vTraitsSize; ++t){
		_selTrait[t]->init();
	}
}

// ----------------------------------------------------------------------------------------
// set_ve_mean
// ----------------------------------------------------------------------------------------
// set the environmental mean if used and the corresponding function pointer
void
TSelection::set_ve_mean()
{
	_Ve_mean_set = (_popPtr->get_parameter("patch_ve_mean")->isSet() ||_popPtr->get_parameter("patch_ve_mean_fem")->isSet());
	if(_Ve_mean_set){
		_popPtr->set_patch_parameter(_vTraitsSize, "patch_ve_mean", "mean_environmental_effect", &Patch::set_localMeanVe);
	}
	for(unsigned int t = 0; t < _vTraitsSize; ++t){    // set the function pointers
		_selTrait[t]->set_ve_mean_func_ptr();
	}
}

// ----------------------------------------------------------------------------------------
// set_ve_var
// ----------------------------------------------------------------------------------------
// set the environmental variance
// here only the input is passed to the patches,
// Ve itself can only be computed when pops are populated if h2 or H2 is defined
// therefore this has to be done in the function execute_before_each_replicate()
// at this time the corresponding function pointer will also be set
void
TSelection::set_ve_var()
{
	_Ve_var_set = (_popPtr->get_parameter("patch_ve_var")->isSet() || _popPtr->get_parameter("patch_ve_var_fem")->isSet());
	if(_Ve_var_set){            // Ve is set by the parameter patch_ve_var (has precedence)
		_popPtr->set_patch_parameter(_vTraitsSize, "patch_ve_var", "heritability", &Patch::set_localh2Ve);

		// heritability cannot change over time for environmental model 1 and 3
		if(_popPtr->get_parameter("patch_ve_var")->isTemporalParam()){
			for(unsigned int t = 0; t < _vTraitsSize; ++t) {
				if(_selTrait[t]->get_Ve_model() == 1 || _selTrait[t]->get_Ve_model() == 3){
					error("The heritability cannot change over time if the parameter 'quanti_environmental_model is set to 1 or 3!\n");
				}
			}
		}
	}
	else{                    	// Ve is set by the parameter quanti_heritability
		for(unsigned int t = 0; t < _vTraitsSize; ++t){
			_selTrait[t]->set_quantiHeritability();
		}
	}
}

// ----------------------------------------------------------------------------------------
// reset_Ve
// ----------------------------------------------------------------------------------------
/** compute the environmental variance for each patch and trait
 * all is set by default to true, i.e. all Ve are set
 * if all=False, the environment has to be set at each generation (model 3)
 */
void
TSelection::reset_Ve(bool all){
	// set the Ve in each patch and trait
	for(unsigned int t = 0; t < _vTraitsSize; ++t){
		if(all || _selTrait[t]->get_Ve_model() == 2 || _selTrait[t]->get_Ve_model() == 4){
			_selTrait[t]->set_Ve();
		}
	}
}

// ----------------------------------------------------------------------------------------
// executeBeforeEachReplicate::
// ----------------------------------------------------------------------------------------
/** compute the environmental variance if used */
void
TSelection::executeBeforeEachReplicate(const unsigned int& rep){
	reset_Ve(true); // it has to be done here, since the pops have to be populated

}

// ----------------------------------------------------------------------------------------
// temporal_change::
// ----------------------------------------------------------------------------------------
/** compute the environmental variance if used */
void
TSelection::temporal_change(const unsigned int& gen){
	reset_Ve(false); // check if the Ve has to be recomputed
}

// ----------------------------------------------------------------------------------------
//
// ----------------------------------------------------------------------------------------








