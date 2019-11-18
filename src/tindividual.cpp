/** @file individual.cpp
*
*   Copyright (C) 2006 Frederic Guillaume    <guillaum@zoology.ubc.ca>
*   Copyright (C) 2018 Frederic Michaud <frederic.a.michaud@gmail.com>

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


#include "tindividual.h"
#include "patch.h"


// ----------------------------------------------------------------------------------------
// constructor
// ----------------------------------------------------------------------------------------
TIndividual::TIndividual ( )
: _sex(MAL),_mother(NULL),_father(NULL),
	_natalPatch(NULL),_currentPatch(NULL),
	_isSelfed(false),_fecundity(0),_trait_nb(0)
{
    _matings[0] = _matings[1] = 0;
	_realizedFecundity[0] = _realizedFecundity[1] = 0;
}

// ----------------------------------------------------------------------------------------
// copy constructor
// ----------------------------------------------------------------------------------------
/** copy constructor used to clone an individual */
TIndividual::TIndividual(const TIndividual& ind)
: _sex(MAL),_mother(NULL),_father(NULL),
	_natalPatch(NULL),_currentPatch(NULL),
	_isSelfed(false),_fecundity(0),_trait_nb(0)
{
    _matings[0] = _matings[1] = 0;
	_realizedFecundity[0] = _realizedFecundity[1] = 0;
	genome.set_protoGenome(ind.genome.get_protoGenome());
	genome.clone(ind.genome);
}

// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
TIndividual * TIndividual::init ()
{
	_isSelfed = false;
  _mother = NULL;
	_father = NULL;
	_natalPatch = NULL;
	_currentPatch = NULL;
	_matings[0] = _matings[1] = 0;
  _realizedFecundity[0] = _realizedFecundity[1] = 0;
  _fitness = my_NAN;

	assert(_trait_nb == Traits.size());

  for(unsigned int i = 0; i < _trait_nb; i++){
		Traits[i]->ini(this);
  }

  return this;
}

// ----------------------------------------------------------------------------------------
// clear
// ----------------------------------------------------------------------------------------
/**Clears the traits container.**/
void
TIndividual::clearTraits ( )
{
	genome.clear();

	for(unsigned int i = 0; i < _trait_nb; ++i){
		delete Traits[i];
	}
	Traits.clear();
	_trait_nb = 0;
}

// ----------------------------------------------------------------------------------------
// reset
// ----------------------------------------------------------------------------------------
void TIndividual::reset ()
{
  _id.clear();
  _sex = MAL;
  _isSelfed = false;
  _mother = NULL;
	_father = NULL;
	_motherID.clear();
	_fatherID.clear();
  _natalPatch = NULL;
  _currentPatch = NULL;
  _matings[0] = _matings[1] = 0;
  _realizedFecundity[0] = _realizedFecundity[1] = 0;
  _fitness = my_NAN;

  if(_trait_nb != Traits.size()){
    error("TIndividual::reset: trait counter and table size differ, resetting\n");
    _trait_nb = (unsigned int)Traits.size();
  }

  reset_traits();
}

// ----------------------------------------------------------------------------------------
// reset_traits
// ----------------------------------------------------------------------------------------
void TIndividual::reset_traits ()
{
  for(unsigned int i = 0; i < _trait_nb; i++){
    Traits[i]->reset();
  }
}

// ----------------------------------------------------------------------------------------
// switch_sex
// ----------------------------------------------------------------------------------------
/** change the sex of the individual.
  * The individual itself does not know its age, not its current index in the container ...
  */
void TIndividual::switch_sex (age_idx AGE, const int& i){
  if(this->_sex == MAL){    // it is a male
      _sex = FEM;                           // change sex
      assert(_currentPatch->get(MAL, AGE, i)==this);
      _currentPatch->add(FEM, AGE, this);   // add individual to the female container
      _currentPatch->remove(MAL, AGE, i);   // remove individual from the male container
  }
  else{                     // it is a female
      _sex = MAL;
      assert(_currentPatch->get(FEM, AGE, i)==this);
      _currentPatch->add(MAL, AGE, this);   // add individual to the male container
      _currentPatch->remove(FEM, AGE, i);   // remove individual from the female container
  }
}

// ----------------------------------------------------------------------------------------
// show_up
// ----------------------------------------------------------------------------------------
void TIndividual::show_up ()
{
  message("\nIndividual ID: %s\n\
          sex: %i\n\
       mother: %s\n\
       father: %s\n\
  natal patch: %i\n\
current patch: %i\n\
        traits values: \n",_id.c_str(),_sex,_motherID.c_str(),_fatherID.c_str(),_natalPatch->get_ID(),_currentPatch->get_ID());

  for(unsigned int i = 0; i < _trait_nb; i++){
    Traits[i]->show_up();
  }
}
// ----------------------------------------------------------------------------------------
// clone
// ----------------------------------------------------------------------------------------
TIndividual* TIndividual::clone ()
{
	TIndividual* myClone = new TIndividual(*this);

  for(unsigned int i = 0; i < _trait_nb; i++){
		myClone->addTrait(Traits[i]->clone(), i);
  }

	return myClone;
}

// ----------------------------------------------------------------------------------------
// inherit
// ----------------------------------------------------------------------------------------
/**Calls the inheritance procedure of all the traits present in the individual.
  * @param mother the mother
  * @param father the father
  **/
void
TIndividual::inherit (TIndividual* mother, TIndividual* father){
	genome.inherit(mother, father);
}

// ----------------------------------------------------------------------------------------
// operator=
// ----------------------------------------------------------------------------------------
TIndividual& TIndividual::operator=(const TIndividual& i)
{
  if(this != &i) {

    if(Traits.size() != i.Traits.size()) error("TIndividual::operator=:not same number of traits in left and right sides of assignment\n");
    if(_trait_nb != i._trait_nb) {
      error("TIndividual::operator=:trait counters differ, restting\n");
      _trait_nb = i._trait_nb;
    }
    _sex                  = i._sex;
    _mother               = i._mother;
    _father               = i._father;
    _natalPatch           = i._natalPatch;
    _currentPatch         = i._currentPatch;
    _isSelfed             = i._isSelfed;
    _fecundity            = i._fecundity;
    _matings[0]           = i._matings[0];
    _matings[1]           = i._matings[1];
    _realizedFecundity[0] = i._realizedFecundity[0];
    _realizedFecundity[1] = i._realizedFecundity[1];

    for(unsigned int t = 0; t < _trait_nb; t++){
      if(Traits[t]->pTraitProto->get_type_index() != i.Traits[t]->pTraitProto->get_type_index()){
        error("TIndividual::operator=: not same kinds of traits in left and right sides of assignment\n");
      }
      (*Traits[t]) = (*i.Traits[t]);
    }
  }
  return *this;
}
// ----------------------------------------------------------------------------------------
// operator==
// ----------------------------------------------------------------------------------------
bool TIndividual::operator==(const TIndividual& i)
{
  if(this != &i) {
    if(Traits.size() != i.Traits.size()) return false;

    for(unsigned int t = 0; t < Traits.size(); t++){
      if((*Traits[t]) != (*i.Traits[t])) return false;
    }
  }
  return true;
}

// ----------------------------------------------------------------------------------------
// operator!=
// ----------------------------------------------------------------------------------------
bool TIndividual::operator!=(const TIndividual& i)
{
  if(!((*this) == i)) return true;
  return false;
}

