/** @file node.cpp
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

#include "node.h"
#include "types.h"
//#include "functions.h"

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                                ******** Node ********

// ----------------------------------------------------------------------------------------
// Node
// ----------------------------------------------------------------------------------------
Node::Node(unsigned int depth, unsigned int rank, unsigned int size)
{
  set(depth, rank, size);
}

// ----------------------------------------------------------------------------------------
// set
// ----------------------------------------------------------------------------------------
void Node::set(unsigned int depth, unsigned int rank, unsigned int size)
{
  _myrank = rank;
  _nb_branches = size;

  if(rank == depth) {      // is a terminal leaf
  	_branches = NULL;
	  _values = new double [size];
	  for(unsigned int i = 0; i < size; ++i){
	    _values[i] = my_NAN;
    }
  }
  else {                  // is a node
	  _values = NULL;
	  _branches = new Node* [size];
	  for(unsigned int i = 0; i < size; ++i){
	    _branches[i] = NULL;
    }
  }
}

// ----------------------------------------------------------------------------------------
// ~Node
// ----------------------------------------------------------------------------------------
Node::~Node()
{
  if(_values) delete[]_values;
  if(_branches) {
	  for(unsigned int i = 0; i < _nb_branches; ++i){
	    if(_branches[i]) delete _branches[i];
    }
	  delete [] _branches;
  }
}

// ----------------------------------------------------------------------------------------
// get_value
// ----------------------------------------------------------------------------------------
double Node::get_value(unsigned int *coord, unsigned int depth)
{
  // if it is the leave
  if(depth == _myrank){
	  return _values[coord[_myrank]];
  }

  //if the node does not exist
  if(!_branches[ coord[_myrank]]) return my_NAN;

  // if the node exists
  return _branches[coord[_myrank]]->get_value(coord, depth);
}

// ----------------------------------------------------------------------------------------
// get_next
// ----------------------------------------------------------------------------------------
double Node::get_next(unsigned int *coord, unsigned int depth)
{
  unsigned int pos = coord[_myrank];

  // if it is the leave
  if(depth == _myrank){
    for(; pos<_nb_branches; ++pos){
      if(_values[pos] != my_NAN){
        coord[_myrank] = pos;
        return _values[pos];
      }
    }
    coord[_myrank]=0;
	  return my_NAN;
  }

  // find the next node
  double value;
  for(; pos<_nb_branches; ++pos){
    if(_branches[pos]){
      value = _branches[pos]->get_next(coord, depth);
      if(value != my_NAN){
        coord[_myrank] = pos;
        return value;
      }
    }
  }
  coord[_myrank]=0;
  return my_NAN;
}

// ----------------------------------------------------------------------------------------
// set_value
// ----------------------------------------------------------------------------------------
void Node::set_value(unsigned int *coord, unsigned int depth, double value)
{
  // if it is the leave
  if(depth == _myrank){
	  _values[coord[_myrank]] = value;
    return;
  }
  // if the node does not exist create it
  if(!_branches[ coord[_myrank]]){
	  _branches[coord[_myrank]] = new Node(depth, _myrank+1, _nb_branches);
  }

  // go to the node
  return _branches[coord[_myrank]]->set_value(coord, depth, value);
}
