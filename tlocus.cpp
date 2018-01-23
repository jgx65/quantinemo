/** @file tlocus.cpp
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


#include "tlocus.h"
#include "metapop.h"
#include "ttrait.h"
//---------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------
// set_ini_sequence_of_locus
// ----------------------------------------------------------------------------------------
/** set the sequence initialization model of locus l */
void
TLocus::set_ini_sequence_model(const unsigned int& model)
{
	switch(model){
		case 0: _ini_seq_model_func_ptr = &TLocus::_ini_sequence_uniform;    break; // uniformly random
		case 1: _ini_seq_model_func_ptr = &TLocus::_ini_sequence_monomorph;  break; // a single allele
		case 2: _ini_seq_model_func_ptr = &TLocus::_ini_sequence_dist;       break; // following a distribution
	}
}

// ----------------------------------------------------------------------------------------
/** set the sequence initialization model of locus l */
void
TLocus::set_ini_sequence_model(const ini_model_t& model)
{
	switch(model){
		case INI_UNIF: _ini_seq_model_func_ptr = &TLocus::_ini_sequence_uniform;    break; // uniformly random
		case INI_MONO: _ini_seq_model_func_ptr = &TLocus::_ini_sequence_monomorph;  break; // a single allele
		case INI_DIST: _ini_seq_model_func_ptr = &TLocus::_ini_sequence_dist;       break; // following a distribution
	}
}

// ----------------------------------------------------------------------------------------
// constructor
// ----------------------------------------------------------------------------------------
TLocus::TLocus(){
	reset();
}

// ----------------------------------------------------------------------------------------
// destructor
// ----------------------------------------------------------------------------------------
TLocus::~TLocus(){
	if(_allelic_effect) delete[] _allelic_effect;
	if(_mutationFreq) delete[] _mutationFreq;
	if(_ini_frequencies) delete[] _ini_frequencies;
}

// ----------------------------------------------------------------------------------------
// functions to initialize a sequence
// ----------------------------------------------------------------------------------------
/** initialization of the sequences based on the _initAlleleFreq array (settings file) */
void
TLocus::_ini_sequence_dist(unsigned char* seq, Patch* patch, unsigned int size)
{
	double* freq_dist = _pTrait->get_ini_alleleFreq(patch, _locus_id_trait);
	for(unsigned int a = 0; a < size; ++a) {
		seq[a] = (unsigned char)_popPtr->rand().AfterDistribution(freq_dist, _nb_allele);
	}
}

/** initialization of the sequences: alleles are randomly distributed -> maximal variance */
void
TLocus::_ini_sequence_uniform(unsigned char* seq, Patch* patch, unsigned int size)
{
	for(unsigned int a = 0; a < size; ++a) {
		seq[a] = (unsigned char)_popPtr->rand().Uniform(_nb_allele);
	}
}

/** initialization of the sequences: populations are fixed for the middle allele -> minimal variance
 * if number of alleles is even, then the lower allele is used, ie.e nb_allele=2 => 0*/
void
TLocus::_ini_sequence_monomorph(unsigned char* seq, Patch* patch, unsigned int size)
{
	for(unsigned int a = 0; a < size; ++a) {
		seq[a] = (unsigned char)(_nb_allele-1)/2;
	}
}

// ----------------------------------------------------------------------------------------
// functions to initialize a sequence
// ----------------------------------------------------------------------------------------
void
TLocus::set_mutationModel(const mut_model_t& model)
{
	_mut_model = model;
	switch(model){
		case KAM: _mut_model_func_ptr = &TLocus::_mutate_KAM;  break; // KAM
		case SSM: _mut_model_func_ptr = &TLocus::_mutate_SSM;  break; // SSM
		case RMM: _mut_model_func_ptr = &TLocus::_mutate_RMM;  break; // RMM
		case IMM: _mut_model_func_ptr = &TLocus::_mutate_IMM;  break; // IMM
		case NO:  _mut_model_func_ptr = &TLocus::_mutate_none; break; // none
	}

}


// ----------------------------------------------------------------------------------------
// mutate a single locus
// ----------------------------------------------------------------------------------------
/** check if a mutation occurs at this locus */
void
TLocus::mutate(unsigned char* seq)
{
	unsigned int NbMut;
	for(NbMut = _popPtr->rand().Poisson(ploidy*_mut_rate) ; NbMut != 0; --NbMut) {
		(this->*_mut_model_func_ptr)(&seq[_popPtr->rand().Uniform((unsigned int)ploidy)]);
	}
}

// ----------------------------------------------------------------------------------------
/** a mutation occurs, make it */
void
TLocus::mutate_now(unsigned char* curPos)
{
	(this->*_mut_model_func_ptr)(curPos);
}

// ----------------------------------------------------------------------------------------
// mutate_KAM
// ----------------------------------------------------------------------------------------
void
TLocus::_mutate_KAM(unsigned char* curPos)
{
	//assign an arbitrary allele value:
	unsigned char mut;

	do{
		mut = (unsigned char) _popPtr->rand().Uniform(_nb_allele);
	} while (mut == *curPos);  // it has to change
	*curPos = mut;
}

// ----------------------------------------------------------------------------------------
// mutate_SSM_all
// ----------------------------------------------------------------------------------------
/** SMM mutation model with a common mutation model and rate */
void
TLocus::_mutate_SSM(unsigned char* curPos)
{
	//alleles values are from 0 to NtrlAll - 1 !!!
	if(_popPtr->rand().Bool() && (unsigned int)*curPos < _nb_allele-1){   // Direction
		*curPos += 1; //one step to the right
	}
	else if(*curPos > 0){ // !direction || all==_nb_allele
		*curPos -= 1; //one step to the left
	}
	else{ //!direction && all == 0
		*curPos += 1; // on step to the right if the current allele is 0
	}
}

// ----------------------------------------------------------------------------------------
// mutate_RMM
// ----------------------------------------------------------------------------------------
/** At each mutation a new allelic effect is drawn proportional to their frequencies
	* RMM: rekursive mutation model
	*/
void
TLocus::_mutate_RMM(unsigned char* curPos)
{
	unsigned char mutPos;
	unsigned int i=0;

	//assign an arbitrary allele value:
	do{
		mutPos = (unsigned char) (_popPtr->rand().AfterDistribution(_pTrait->get_mutationFreq(_locus_id_trait), _nb_allele));
		++i;                                 // used to asure that the loop ends
    } while (mutPos == *curPos && i<TGenomeProto::MUTATION_TRIAL);  // it has to change (except if after 1e4 trials it has not yet changed)
	*curPos = mutPos;
}

// ----------------------------------------------------------------------------------------
// mutate_IMM
// ----------------------------------------------------------------------------------------
/** Incremental mutation: at each mutation an allelic effect is drawn proportional to their frequencies
	* This drawn allelic affect is added to the present one.
	* If the new allele is out of range, nothing is made (following fl)
	*/
void
TLocus::_mutate_IMM(unsigned char* curPos)
{
	unsigned int mutStep, newPos;
	mutStep = _popPtr->rand().AfterDistribution(_pTrait->get_mutationFreq(_locus_id_trait), _nb_allele) - (_nb_allele-1)/2;  // _nb_allele/2 is the middle index
	newPos = mutStep+(unsigned int)(*curPos);
	if(newPos<_nb_allele){   // do nothing if the new allele is out of range
		*curPos += (unsigned char)mutStep;
	}
}

// ----------------------------------------------------------------------------------------
/** a mutation occurs, make it (the random deviate is passed)
	* true is returned if mutation was successful and false if not
	*/
bool
TLocus::mutate_now(unsigned char* curPos, const double& ran)
{
	return (this->*_mut_model2_func_ptr)(curPos, ran);
}

// ----------------------------------------------------------------------------------------
// mutate_RMM
// ----------------------------------------------------------------------------------------
/** At each mutation a new allelic effect is drawn proportional to their frequencies
	* RMM: rekursive mutation model
	*/
bool
TLocus::_mutate_RMM(unsigned char* curPos, const double& ran)
{
	unsigned char mutPos;
	unsigned int i=0;

	//assign an arbitrary allele value:
	do{
		mutPos = (unsigned char) (_popPtr->rand().AfterDistribution(_pTrait->get_mutationFreq(_locus_id_trait), _nb_allele));
		++i;                                 // used to asure that the loop ends
		} while (mutPos == *curPos && i<1e4);  // it has to change (except if after 1e4 trials it has not yet changed)
	*curPos = mutPos;
	return true;
}

// ----------------------------------------------------------------------------------------
// mutate_IMM
// ----------------------------------------------------------------------------------------
/** Incremental mutation: at each mutation an allelic effect is drawn proportional to their frequencies
	* This drawn allelic affect is added to the present one.
	* If the new allele is out of range, nothing is made (following fl)
	*/
bool
TLocus::_mutate_IMM(unsigned char* curPos, const double& ran)
{
	unsigned int mutStep, newPos;
	mutStep = _popPtr->rand().AfterDistribution(_pTrait->get_mutationFreq(_locus_id_trait), _nb_allele) - (_nb_allele-1)/2;  // _nb_allele/2 is the middle index
	newPos = mutStep+(unsigned int)(*curPos);
	if(newPos<_nb_allele){   // do nothing if the new allele is out of range
		*curPos += (unsigned char)mutStep;
	}
	return true;
}

// ----------------------------------------------------------------------------------------
// set_locus_id_tot
// ----------------------------------------------------------------------------------------
/** sets the locus index of the genome at the locus level and at the trait level */
void
TLocus::set_locus_id_tot(const unsigned int& n){
	_locus_id_tot = n;
	_pTrait->add_locus_index(_locus_id_trait, _locus_id_tot);
}

// ----------------------------------------------------------------------------------------
// reset
// ----------------------------------------------------------------------------------------
void
TLocus::reset(){
	_pTrait = NULL;
    _popPtr = NULL;
	_allelic_effect = NULL;
	_mutationFreq = NULL;
	_mut_model_func_ptr = NULL;
	_ini_frequencies = NULL;
	_ini_seq_model_func_ptr = NULL;
	_MNR_param = NULL;
	_locusPos = NULL;
}

// ----------------------------------------------------------------------------------------
// set_traitPointer
// ----------------------------------------------------------------------------------------
void
TLocus::set_traitPointer(TTraitProto* p)
{
    _pTrait = p;
    _popPtr = p->get_popPtr();
}

