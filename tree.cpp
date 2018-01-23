/** @file tree.cpp
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

#include "tree.h"
#include "types.h"

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                                ******** Tree ********

// ----------------------------------------------------------------------------------------
// Tree
// ----------------------------------------------------------------------------------------
template <class T>
Tree<T>::Tree(unsigned int nbloc, unsigned int nball)
{
  _nb_locus = nbloc;
  _nb_all = nball;
  _nb_branches = nball*(nball+1) / 2;

  _root.set(nbloc-1, 0, _nb_branches);

  _mapper = new unsigned int* [nball];
	for(unsigned int i = 0; i < nball; ++i){
		_mapper[i]    = new unsigned int [nball];
	}
	_un_mapper = new T*[_nb_branches];
	for(unsigned int i = 0; i < _nb_branches; ++i){
    _un_mapper[i] = new T[2];  // ploidy
  }

  unsigned int cntr = 0;
	for(unsigned int i = 0; i < nball; ++i) {
		_un_mapper[cntr][0] = (T)i;
		_un_mapper[cntr][1] = (T)i;
		_mapper[i][i] = cntr++;
		for(unsigned int j = i + 1; j < nball; ++j){
			_un_mapper[cntr][0] = (T)i;
			_un_mapper[cntr][1] = (T)j;
			_mapper[i][j] = _mapper[j][i] = cntr++;
		}
	}

	_coord = new unsigned int[nbloc];
}

// ----------------------------------------------------------------------------------------
// ~Tree
// ----------------------------------------------------------------------------------------
template <class T>
Tree<T>::~Tree()
{
  if(_coord) delete [] _coord;
  if(_mapper){
    for(unsigned int i = 0; i < _nb_all; ++i){
	    delete [] _mapper[i];
    }
    delete [] _mapper;
  }
  if(_un_mapper){
    for(unsigned int i = 0; i < _nb_branches; ++i){
	    delete [] _un_mapper[i];
    }
    delete [] _un_mapper;
  }
}

// ----------------------------------------------------------------------------------------
// get_value
// ----------------------------------------------------------------------------------------
template <class T>
double Tree<T>::get_value(T** genotype)
{
  /**The genotype is converted in a tree coordinate, each locus receives a value corresponding to
    its genotype's value in the mapper. This coordinate is used to scan the tree to find the phenotype leaf.*/
  for(unsigned int i = 0; i < _nb_locus; i++){
		_coord[i] = _mapper[ (unsigned int)genotype[i][0] ][ (unsigned int)genotype[i][1] ];
	}

  return _root.get_value(_coord, _nb_locus-1);
}

// ----------------------------------------------------------------------------------------
// get_next
// ----------------------------------------------------------------------------------------
/** returns the first avaiable genotype and its a value (the genotype is altered accordignly).
  * If there is no one, it returns my_NAN
  */
template <class T>
double Tree<T>::get_first(T** genotype)
{
  // translate the genotype
  for(unsigned int i = 0; i < _nb_locus; i++){
	  _coord[i] = 0;       // first one
  }

  // find the next value
  double value = _root.get_next(_coord, _nb_locus-1);

  // recreate the genotype
  for(unsigned int i = 0; i < _nb_locus; i++){
    genotype[i][0] = _un_mapper[_coord[i]][0];
    genotype[i][1] = _un_mapper[_coord[i]][1];
  }

  return value;
}

// ----------------------------------------------------------------------------------------
// get_next
// ----------------------------------------------------------------------------------------
/** returns the next avaiable genotype and its value (the genotype is altered accordignly).
  * If there is no one, it returns my_NAN
  */
template <class T>
double Tree<T>::get_next(T** genotype)
{
  // translate the genotype
  for(unsigned int i = 0; i < _nb_locus; ++i){
	  _coord[i] = _mapper[ (unsigned int)genotype[i][0] ][ (unsigned int)genotype[i][1] ];
  }

  // set it to the next possible genotype
  int i = (int) _nb_locus-1;
  ++_coord[i];
  while(_coord[i] >= _nb_branches){
    _coord[i]=0;
    --i;
    if(i<0){
      return my_NAN;   // if it was the last element
    }
    ++_coord[i];
  }

  // find the next value
  double value = _root.get_next(_coord, _nb_locus-1);

  // recreate the genotype
  for(unsigned int i = 0; i < _nb_locus; i++){
    genotype[i][0] = _un_mapper[_coord[i]][0];
    genotype[i][1] = _un_mapper[_coord[i]][1];
  }

  return value;
}

// ----------------------------------------------------------------------------------------
// set_value
// ----------------------------------------------------------------------------------------
template <class T>
void Tree<T>::set_value(T** genotype, double value)
{
  /**The genotype is converted in a tree coordinate, each locus receives a value corresponding to
    its genotype's value in the mapper. This coordinate is used to scan the tree to find the phenotype leaf.*/
  for(unsigned int i = 0; i < _nb_locus; i++){
	  _coord[i] = _mapper[ (unsigned int)genotype[i][0] ][ (unsigned int)genotype[i][1] ];
  }

  return _root.set_value(_coord, _nb_locus-1, value);
}

