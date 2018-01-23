/** @file node.h
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

#ifndef nodeH
#define nodeH

//#include <iostream>
//using namespace std;


/**Abstraction of the node of a Tree, connecting branches and leaves.*/
class Node {

private:
  /**The position of the Node in the Tree.*/
  unsigned int _myrank;
  /**The number of branches leaving this Node.*/
  unsigned int _nb_branches;
  /**The table of branches pointers.*/
  Node  ** _branches;
  /**The table of phenotypic values, used only if the Node is a leaf.*/
  double * _values;

public:

  Node  ( ) { }

  /**@param depth the depth of the Tree, number of ranks (locus) in the Tree
    @param rank the rank of this Node, if rank == depth, this is a leaf
    @param size the number of branches and values to store in the Node.*/
  Node  (unsigned int depth, unsigned int rank, unsigned int size);

  ~Node ();

  /**Set the node params, branches are NULL by default, they are created only if needed.
    @param depth the depth of the Tree, number of ranks (locus) in the Tree
    @param rank the rank of this Node, if rank == depth, this is a leaf
    @param size the number of branches and values to store in the Node.*/
  void set   (unsigned int depth, unsigned int rank, unsigned int size);

  /**Scans the Tree to find the phenotype value of the mapped genotype.
    If a node at a given branch is missing, it is created.
    @param coord the coordinate values of the genotype
    @param depth the depth of the tree, used to determine if we are at a leaf
  */
  double get_value(unsigned int *coord, unsigned int depth);
  double get_next(unsigned int *coord, unsigned int depth);
  void   set_value(unsigned int *coord, unsigned int depth, double value);
};
#endif //NODE_H

