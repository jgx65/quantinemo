/** @file tree.cpp
*
*   Copyright (C) 2006 Frederic Guillaume    <guillaum@zoology.ubc.ca>
*   Copyright (C) 2008 Samuel Neuenschwander <samuel.neuenschwander@unil.ch>
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

#include "genotype_value_map.h"
#include "types.h"

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                          ******** GenotypeValueMap ********

// ----------------------------------------------------------------------------------------
// GenotypeValueMap
// ----------------------------------------------------------------------------------------
GenotypeValueMap::GenotypeValueMap(unsigned int nbloc, unsigned int nball)
{
  _nb_locus = nbloc;
  _nb_all = nball;
  _nb_branches = nball*(nball+1) / 2;

  _root.set(nbloc-1, 0, _nb_branches);

  // _pair_to_coord[a1][a2]: unordered locus genotype -> compact coordinate (branch index)
  _pair_to_coord = new unsigned int* [nball];
  for(unsigned int i = 0; i < nball; ++i){
    _pair_to_coord[i] = new unsigned int [nball];
  }

  unsigned int cntr = 0;
  for(unsigned int i = 0; i < nball; ++i) {
    _pair_to_coord[i][i] = cntr++;
    for(unsigned int j = i + 1; j < nball; ++j){
      _pair_to_coord[i][j] = _pair_to_coord[j][i] = cntr++;
    }
  }

  _coord = new unsigned int[nbloc];
}

// ----------------------------------------------------------------------------------------
// ~GenotypeValueMap
// ----------------------------------------------------------------------------------------
GenotypeValueMap::~GenotypeValueMap()
{
  if(_coord) delete [] _coord;
  if(_pair_to_coord){
    for(unsigned int i = 0; i < _nb_all; ++i){
	    delete [] _pair_to_coord[i];
    }
    delete [] _pair_to_coord;
  }
}

// ----------------------------------------------------------------------------------------
// get_value
// ----------------------------------------------------------------------------------------
double GenotypeValueMap::get_value(AlleleContainer& genotype, const unsigned int* genome_locus)
{
  // compress each locus genotype to its coordinate, then walk the trie to the leaf
  for(unsigned int i = 0; i < _nb_locus; i++){
		_coord[i] = _pair_to_coord[ (unsigned int)genotype.allele(genome_locus?genome_locus[i]:i, 0) ][ (unsigned int)genotype.allele(genome_locus?genome_locus[i]:i, 1) ];
	}

  return _root.get_value(_coord, _nb_locus-1);
}

// ----------------------------------------------------------------------------------------
// set_value
// ----------------------------------------------------------------------------------------
void GenotypeValueMap::set_value(AlleleContainer& genotype, const unsigned int* genome_locus, double value)
{
  for(unsigned int i = 0; i < _nb_locus; i++){
	  _coord[i] = _pair_to_coord[ (unsigned int)genotype.allele(genome_locus?genome_locus[i]:i, 0) ][ (unsigned int)genotype.allele(genome_locus?genome_locus[i]:i, 1) ];
  }

  return _root.set_value(_coord, _nb_locus-1, value);
}

