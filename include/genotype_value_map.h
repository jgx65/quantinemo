/** @file genotype_value_map.h
*
*   Copyright (C) 2006 Frederic Guillaume    <guillaum@zoology.ubc.ca>
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

#ifndef genotypeValueMapH
#define genotypeValueMapH

//#include <iostream>

#include "node.h"
#include "types.h"   // AlleleContainer
//#include "functions.h"
//using namespace std;


/** Sparse, lazily-grown trie that stores a value (an epistatic contribution, a genome-wide
  * fitness factor, ...) for a whole diploid multi-locus genotype. The order of the two alleles
  * at a locus does not matter, so each locus genotype is first compressed to a compact
  * coordinate. Values are authoritative (drawn/read once, then stored), not a recomputable cache.
  */
class GenotypeValueMap {

private:
  /**The depth of the tree, = the number of loci of the trait.*/
  unsigned int _nb_locus;
  /**The number of branches per node, = the number of possible genotypes at a locus.*/
  unsigned int _nb_branches;
  /**The number of allelic states of the trait.*/
  unsigned int _nb_all;
  /**The root node of the trie.*/
  Node _root;
  /**nb_all x nb_all matrix mapping an (unordered) locus genotype to its compact coordinate.*/
  unsigned int** _pair_to_coord;
  /**Scratch of length _nb_locus: the current genotype's coordinate at each locus.*/
  unsigned int * _coord;

public:

  GenotypeValueMap(unsigned int nbloc, unsigned int nball);
  ~GenotypeValueMap();
  /**Value stored for `genotype` (read through `trait_to_genome_locus`, NULL = identity), or my_NAN.*/
  double get_value(AlleleContainer& genotype, const unsigned int* trait_to_genome_locus);
  void   set_value(AlleleContainer& genotype, const unsigned int* trait_to_genome_locus, double value);

};

#endif //genotypeValueMapH

