/** @file TSelectionTrait.cpp
 *
 *   Copyright (C) 2008 Samuel Neuenschwander <samuel.neuenschwander@unil.ch>
 *   Copyright (C) 2018 Frederic Michaud <frederic.a.michaud@gmail.com>

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

#include "tselectiontrait.h"
#include "tselectiontype.h"
#include "tselection.h"
#include "ttquanti.h"
#include "tmetapop.h"
#include "lce_breed.h"
// ----------------------------------------------------------------------------------------
// constructor
// ----------------------------------------------------------------------------------------
TSelectionTrait::TSelectionTrait() : func_ptr_get_meanVe(0), func_ptr_get_sdVe(0)
{
	_selection_pressure[FEM] = NULL;
	_selection_pressure[MAL] = NULL;
	_selection_sd            = NULL;
	_get_selection_pressure_func_ptr = NULL;
	_nb_selection_params = 0;
}

// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
void
TSelectionTrait::init(TSelection* s, const int& trait)
{
	_pSel         = s;
	_popPtr       = _pSel->get_popPtr();
	_quantiIndex  = trait;
	_traitIndex   = _pSel->get_vTraits(trait);
	_pQuantiProto = dynamic_cast <TTraitQuantiProto*> (&_popPtr->getTraitPrototype(_traitIndex));

	// get parameters from quantiProto
	_Ve_prop      = _pQuantiProto->get_parameter("quanti_environmental_proportion"+_pQuantiProto->get_trait_indexStr_())->get_value();
	_Ve_model     = (int)_pQuantiProto->get_parameter("quanti_environmental_model"+_pQuantiProto->get_trait_indexStr_())->get_value();

	// is there any fitness factor set? if fitnessFactor=0 is lethal, then the fitness factor is set just after breeding
	// in the later case the fitness factor has to be set explicitly for the first generation.
	if(_pQuantiProto->fitnessFactor_used() && !_popPtr->get_pBreed_LCE()->get_fitness_factor_zero_isLethal()){
		func_ptr_set_phenotype = &TSelectionDirectional::set_phenotype_and_fitness_factor;
	}
    else func_ptr_set_phenotype = &TSelectionDirectional::set_phenotype;

	// check if the environment shapes the phenotype
	if(!(_pQuantiProto->get_parameter_isSet("quanti_heritability"+_pQuantiProto->get_trait_indexStr_())
			|| _pQuantiProto->get_parameter_isSet("quanti_heritability_fem"+_pQuantiProto->get_trait_indexStr_())
			|| _pQuantiProto->get_parameter_isSet("quanti_heritability_mal"+_pQuantiProto->get_trait_indexStr_())
			|| _popPtr->get_parameter_isSet("patch_ve_mean")
			|| _popPtr->get_parameter_isSet("patch_ve_mean_fem")
			|| _popPtr->get_parameter_isSet("patch_ve_mean_mal")
			|| _popPtr->get_parameter_isSet("patch_ve_var")
			|| _popPtr->get_parameter_isSet("patch_ve_var_fem")
			|| _popPtr->get_parameter_isSet("patch_ve_var_mal"))
			|| _pQuantiProto->get_selection_model()==4){  // selection coefficient

		if(func_ptr_set_phenotype==&TSelectionDirectional::set_phenotype_and_fitness_factor)
			func_ptr_set_phenotype=&TSelectionDirectional::set_phenotype_and_fitness_factor_noVe;
		else if(func_ptr_set_phenotype==&TSelectionDirectional::set_phenotype)
			func_ptr_set_phenotype=&TSelectionDirectional::set_phenotype_noVe;
	}

	if(_pQuantiProto->get_selection_model()==4 && _Ve_model!=0){
		error("Parameter 'quanti_environmental_model' has to be set to 0 when using selection coefficient!");
	}
    
    // is there any fitness factor set?
	if(_pQuantiProto->fitnessFactor_used()) _func_ptr_get_fitness = &TSelectionTrait::get_fitnessFactor;
	else                                    _func_ptr_get_fitness = &TSelectionTrait::get_fitness_none;

}



// ----------------------------------------------------------------------------------------
// set_ve_mean_func_ptr
//-----------------------------------------------------------------------------------------
void
TSelectionTrait::set_ve_mean_func_ptr()
{
	// is there a constant effect of the environment?
	if(_pSel->get_Ve_mean_set()){
		if(_Ve_prop == 1)      func_ptr_get_meanVe = &TSelectionTrait::_get_meanVe_current;
		else if(_Ve_prop == 0) func_ptr_get_meanVe = &TSelectionTrait::_get_meanVe_natal;
		else                   func_ptr_get_meanVe = &TSelectionTrait::_get_meanVe_mix;
	}
	else                     func_ptr_get_meanVe = &TSelectionTrait::_get_null;   // there is no constant effect
}

// ----------------------------------------------------------------------------------------
// _get_sel_pressure
// ----------------------------------------------------------------------------------------
void
TSelectionTrait::set_quantiHeritability(){
	if(_pSel->get_Ve_var_set()) return;

	Param* pParam = _pQuantiProto->get_parameter("quanti_heritability"+_pQuantiProto->get_trait_indexStr_());
	unsigned int nbPatch = _popPtr->getPatchNbr();
	if(pParam->is_matrix()){
		TMatrix* m = pParam->get_matrix();
		unsigned int count_patch = m->get_dims(NULL);

		// check the number of h2 per number of patches
		if(count_patch>nbPatch) warning("There are more heritabilities specified than patches! Only a part of the heritabilities is considered!\n");
		else if(nbPatch % count_patch) warning("The number of heritabilities is not an entire subset of the number of patches!\n");

		// set the h2 for each patch
		for(unsigned int p = 0; p < nbPatch; ++p) {
			_popPtr->get_vPatch(p)->set_h2(p % count_patch, _quantiIndex);
		}
		delete m;
	}
	else{         // a single value
		double val = pParam->get_value();
		for(unsigned int p = 0; p < nbPatch; ++p) {
			_popPtr->get_vPatch(p)->set_h2(val, _quantiIndex);
		}
	}
}

// ----------------------------------------------------------------------------------------
// _get_sel_pressure
// ----------------------------------------------------------------------------------------
double
TSelectionTrait::_get_selection_pressure_var(double value, const int& i){
	return _popPtr->rand().Normal(value, _selection_sd[i]);
}

double
TSelectionTrait::_get_selection_pressure_const(double value, const int& i){
	return value;
}

//-----------------------------------------------------------------------------
// _get_Ve
//-----------------------------------------------------------------------------
/** returns the environmental variance according to the selected proportion
 * between the Ve of the natal patch and the one of the mating patch.
 * _Ve_prop is the proportion of the current patch
 */
//-----------------------------------------------------------------------------
double TSelectionTrait::_get_meanVe_natal(){
	return _curInd->getNatalPatch()->get_meanVe(_curInd->getSex(),_quantiIndex);
}

double TSelectionTrait::_get_meanVe_current(){
	return _curInd->getCurrentPatch()->get_meanVe(_curInd->getSex(), _quantiIndex);
}

double TSelectionTrait::_get_meanVe_mix(){
	return ((1-_Ve_prop) * _curInd->getNatalPatch()->get_meanVe(_curInd->getSex(), _quantiIndex)      // Ve of natal patch
			+ _Ve_prop   * _curInd->getCurrentPatch()->get_meanVe(_curInd->getSex(), _quantiIndex));  // Ve of current patch
}

//-----------------------------------------------------------------------------
double TSelectionTrait::_get_sdVe_natal(){
	return _curInd->getNatalPatch()->get_sdVe(_curInd->getSex(), _quantiIndex);
}

double TSelectionTrait::_get_sdVe_current(){
	return _curInd->getCurrentPatch()->get_sdVe(_curInd->getSex(), _quantiIndex);
}

double TSelectionTrait::_get_sdVe_mix(){
	return ((1-_Ve_prop) * _curInd->getNatalPatch()->get_sdVe(_curInd->getSex(), _quantiIndex)      // Ve of natal patch
			+ _Ve_prop   * _curInd->getCurrentPatch()->get_sdVe(_curInd->getSex(), _quantiIndex));  // Ve of current patch
}

//-----------------------------------------------------------------------------
// set_phenotype_and_fitness_factor
//-----------------------------------------------------------------------------
/** sets the phenotype of each trait of the individuum taking into account Va, Vd, Vep and Ve */
void
TSelectionTrait::set_phenotype_and_fitness_factor(TIndividual* ind)
{
	//_curInd = ind;      // is done in set_phenotype
	set_phenotype(ind);
	_curInd->setFitnessFactor(_traitIndex);
}

//-----------------------------------------------------------------------------
// set_phenotype_and_fitness_factor
//-----------------------------------------------------------------------------
/** sets the phenotype of each trait of the individuum without any Ve */
void
TSelectionTrait::set_phenotype_and_fitness_factor_noVe(TIndividual* ind)
{
	//_curInd = ind;      // is done in set_phenotype
	set_phenotype_noVe(ind);
	_curInd->setFitnessFactor(_traitIndex);
}

//-----------------------------------------------------------------------------
// set_phenotype
//-----------------------------------------------------------------------------
/** sets the phenotype of each trait of the individual taking into account Va, Vd, Vep and Ve */
void
TSelectionTrait::set_phenotype(TIndividual* ind)
{
	_curInd = ind;

	if(_curInd->getTraitPhenotype(_traitIndex) == my_NAN){	// if phenotype is not yet set
		// v=Normal(0, sqrt(Ve))
		double sd = (this->*func_ptr_get_sdVe)();
		if(sd) sd *= _popPtr->rand().Normal();
		_curInd->setTraitValue(_traitIndex, sd+(this->*func_ptr_get_meanVe)());  // P = G + Normal()
	}
}

//-----------------------------------------------------------------------------
// set_phenotype_no_ve
//-----------------------------------------------------------------------------
/** the phenotype is the genotypic value (no environmental variance)
 */
void
TSelectionTrait::set_phenotype_noVe(TIndividual* ind)
{
	_curInd = ind;

	if(_curInd->getTraitPhenotype(_traitIndex) == my_NAN){	// if phenotype is not yet set
		_curInd->setTraitValue(_traitIndex, 0);               // no Ve
	}
}

//-----------------------------------------------------------------------------
// set_phenotype_and_fitness_factor
//-----------------------------------------------------------------------------
/** sets the phenotype of each trait of the individual taking into account Va, Vd, Vep and Ve */
void
TSelectionTrait::set_fitness_factor(TIndividual* ind)
{
	_curInd = ind;
	_curInd->setFitnessFactor(_traitIndex);
}

//-----------------------------------------------------------------------------
// set_phenotype_and_fitness_factor
//-----------------------------------------------------------------------------
/** sets the phenotype of each trait of the individual taking into account Va, Vd, Vep and Ve */
double
TSelectionTrait::set_get_fitness_factor(TIndividual* ind)
{
	_curInd = ind;
	return _curInd->set_getFitnessFactor(_traitIndex);
}

//-----------------------------------------------------------------------------
// get_selection_pressure_tot_var
//-----------------------------------------------------------------------------
/** at least a selection parameter changes between generations, thus all parameters
 * have to be copied to the _selection_pressure array (this array has thus also to be created and deleted!)
 */
void
TSelectionTrait::_get_selection_pressure_tot_var(TPatch* patch, sex_t SEX)
{
	_curPatch = patch;
	_curSex = SEX;
	for(unsigned int i = 0; i< _nb_selection_params; ++i){
		_selection_pressure[_curSex][i]
		                             = (this->*_get_selection_pressure_func_ptr[i])(patch->get_localSelection(SEX,_quantiIndex,i),i);
	}
}

//-----------------------------------------------------------------------------
// get_selection_pressure_tot_const
//-----------------------------------------------------------------------------
/** the selection pressure is constant between generations, thus set simply the
 * pointer to the corresponding patch array (_selection_pressure must not be created nor deleted!)
 */
void
TSelectionTrait::_get_selection_pressure_tot_const(TPatch* patch, sex_t SEX)
{
	_curPatch = patch;
	_curSex = SEX;
	_selection_pressure[_curSex] = patch->get_localSelection(SEX, _quantiIndex);
}

//-----------------------------------------------------------------------------
// get_phenotype
//-----------------------------------------------------------------------------
double
TSelectionTrait::get_phenotype(){
	return _curInd->getTraitValue(_traitIndex);
}

//-----------------------------------------------------------------------------
// get_genotype
//-----------------------------------------------------------------------------
double
TSelectionTrait::get_genotype(){
    return _curInd->getTraitGenotype(_traitIndex);
}


//-----------------------------------------------------------------------------
// get_fitnessFactor
//-----------------------------------------------------------------------------
double
TSelectionTrait::get_fitnessFactor_trait(){
	return _curInd->getFitnessFactor(_traitIndex);
}

//-----------------------------------------------------------------------------
// set_Ve
//-----------------------------------------------------------------------------
/** set the environmental variance for each patch
 * and set the pointer accordingly
 */
void
TSelectionTrait::set_Ve(){
	unsigned int Ve_used = 0;
	for(unsigned int p = 0; p < _popPtr->getPatchNbr(); ++p){
		Ve_used |= _popPtr->get_vPatch(p)->set_Ve(_Ve_model, _quantiIndex, FEM);
		if(_popPtr->get_sexInitRatio()) Ve_used |= _popPtr->get_vPatch(p)->set_Ve(_Ve_model, _quantiIndex, MAL);
	}
	if(Ve_used){
		if(_Ve_prop == 1)      func_ptr_get_sdVe = &TSelectionTrait::_get_sdVe_current;
		else if(_Ve_prop == 0) func_ptr_get_sdVe = &TSelectionTrait::_get_sdVe_natal;
		else                   func_ptr_get_sdVe = &TSelectionTrait::_get_sdVe_mix;
	}
	else                     func_ptr_get_sdVe = &TSelectionTrait::_get_null;
}




