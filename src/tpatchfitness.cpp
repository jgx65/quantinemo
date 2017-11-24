/** @file tpatchfitness.cpp
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

#include "tpatchfitness.h"
#include "tmetapop.h"
//---------------------------------------------------------------------------
// constructor
// -----------------------------------------------------------------------------
void
TPatchFitness::init(TMetapop* ptr, int size)
{
    assert(ptr);
    _pop=ptr;
    if(size){
        _nbInd = _capacity = size;
        _aInd  = new TIndividual*[_nbInd];
        _aFit  = new double[_nbInd];
    }
    else{
        _nbInd = _capacity = 0;
        _aInd  = NULL;
        _aFit  = NULL;
    }
    _nbSubset = 0;
    _sort  = 10;
}

// -----------------------------------------------------------------------------
// constructor
// -----------------------------------------------------------------------------
void
TPatchFitness::resize(const unsigned int& size)
{
    if(size > _capacity){
        if(_aInd) delete[] _aInd;
        if(_aFit) delete[] _aFit;
        _aInd  = new TIndividual*[size];
        _aFit  = new double[size];
        
    }
    _nbInd = size;
}

// -----------------------------------------------------------------------------
// constructor
// -----------------------------------------------------------------------------
double
TPatchFitness::getMeanFitness(){
    if(!_nbInd)return 0;
    return getSumFitness()/_nbInd;
}

double
TPatchFitness::getSumFitness(){
    if(!_nbInd)return 0;
    
    switch(_sort){
            // cumulative
        case -3:
        case -2:
        case -1:
        case 0 :  return _aFit[_nbInd-1];
            
            // reverse cumulative
        case 1 :
        case 2 :
        case 3 :{
            double  sum = 1.0/_aFit[0];	// first indvidual
            for(unsigned int i=1; i<_nbInd; ++i){
                sum += 1.0/(_aFit[i]-_aFit[i-1]);
            }
            return sum;
        }
            
            // array is not cumulative
        case 10: 	return ARRAY::sum(_aFit, _nbInd);
    }
    error("TPatchFitness::getSumFitness not correctly implemented!\n");
    return 1;
}

// -----------------------------------------------------------------------------
// _sort
// -----------------------------------------------------------------------------
/** the function is called when the first time a individual is drawn (_sort == -1)
	* sort:   10: the fitnesses were never used before (don't make a cumulative array
	*         -3: decremental sort (highest fitness first) subset random most fittest
	*         -2: decremental sort (highest fitness first) subset fix most fittest
	*         -1: decremental sort (highest fitness first) -> normal cumulative
	*          0: no sort                                  -> normal cumulative
	*          1: incremental sort (lowest fitness first)  -> reverse cumulative
	*          2: incremental sort (lowest fitness first) subset fix less fittest
	*          3: incremental sort (lowest fitness first) subset random less fittest
	* nbInd:   0: all individuals are used (default)
	*          else : a subset has to be created
	*/
void
TPatchFitness::sort(int sort, unsigned int nbInd){
    _sort = sort;
    
    // check the subset size
    if(sort==-3 || sort==-2){       // decrement sort subset
        assert(nbInd>0);              // nbInd must be specified
        if(nbInd>=_nbInd)_sort = -1;  // if the subset is bigger then the number of individuals no subset is needed
    }
    else if(sort==2 || sort==3){    // increment sort subset
        assert(nbInd>0);              // nbInd must be specified
        if(nbInd>=_nbInd)_sort = 1;   // if the subset is bigger then the number of individuals no subset is needed
    }
    
    switch(sort){
            // decrement (most fittest):
        case -3: // random subset
            randomize_order_mostFit(_aFit, _aInd, _nbInd, nbInd);
            ARRAY::cumulative(_aFit, _nbInd);
            break;
        case -2: // fixed subset
            ARRAY::quicksortBoth(_aFit, _nbInd, _aInd, false);  // fittest first
            ARRAY::cumulative(_aFit, _nbInd);
            break;
        case -1: // no subset
            ARRAY::cumulative(_aFit, _nbInd);
            break;
            
            // neutral (fitness not considered, for testing)
        case  0:
            randomize_order_noFit(_aFit, _aInd, _nbInd);      // the order has to be randomized
            ARRAY::cumulative(_aFit, _nbInd);
            break;
            
            // increment (less fittest)
        case  1: // no subset
            ARRAY::reverse_cumulative(_aFit, _nbInd);
            break;
        case  2: // fixed subset
            ARRAY::quicksortBoth(_aFit, _nbInd, _aInd, true);   // less fittest first
            ARRAY::reverse_cumulative(_aFit, _nbInd);
            break;
        case  3: // random subset
            randomize_order_lessFit(_aFit, _aInd, _nbInd, nbInd);
            ARRAY::reverse_cumulative(_aFit, _nbInd);
            break;
    }
}

// ----------------------------------------------------------------------------------------
// set_subset_random::
// ----------------------------------------------------------------------------------------
/** Randomize the order of the indivduals */
void
TPatchFitness::randomize_order_noFit(double* aFit, TIndividual** aInd, unsigned int& size, unsigned int nbSubset){
    // we start randomly drawing individuals for the first position
    // this fiest individual is then swapted with the drawn individual
    unsigned int pos;
    if(!size || size<nbSubset) nbSubset=size;
    for(unsigned int i=0; i<nbSubset; ++i){
        pos = _pop->rand().Uniform(size-i);  // get the random individual form the curent position (included) to the last position
        if(pos != i){                 			 // if it is not the same element swap them
            ARRAY::swap(aFit, pos, i+pos);     // swap the fitness
            ARRAY::swap(aInd, pos, i+pos);     // swap the individual
        }
    }
}

// ----------------------------------------------------------------------------------------
// set_subset_random::
// ----------------------------------------------------------------------------------------
/** Randomize the order of the indivduals */
void
TPatchFitness::randomize_order_mostFit(double* aFit, TIndividual** aInd, unsigned int& size, unsigned int nbSubset){
    // we start randomly drawing individuals for the first position
    // this first individual is then swapted with the drawn individual
    unsigned int pos;
    if(!size || size<nbSubset) _nbSubset = size;
    else                       _nbSubset = nbSubset;
    
    for(unsigned int i=0; i<nbSubset; ++i){
        pos = _pop->rand().AfterDistribution(aFit+i, size-i);  // get the random individual form the current position (included) to the last position
        if(pos != i){                                   			 // if it is not the same element swap them
            ARRAY::swap(aFit, pos, i+pos);                       // swap the fitness
            ARRAY::swap(aInd, pos, i+pos);                       // swap the individual
        }
    }
}

// ----------------------------------------------------------------------------------------
// set_subset_random::
// ----------------------------------------------------------------------------------------
/** Randomize the order of the indivduals */
void
TPatchFitness::randomize_order_lessFit(double* aFit, TIndividual** aInd, unsigned int& size, unsigned int nbSubset){
    // we start randomly drawing individuals for the first position
    // this first individual is then swapted with the drawn individual
    unsigned int pos;
    if(!size || size<nbSubset) _nbSubset = size;
    else                       _nbSubset = nbSubset;
    for(unsigned int i=0; i<_nbSubset; ++i){
        pos = _pop->rand().AfterDistribution(aFit+i, size-i, -1); // get the random individual form the current position (included) to the last position
        if(pos != i){                                   					// if it is not the same element swap them
            ARRAY::swap(aFit, pos, i+pos);                       		// swap the fitness
            ARRAY::swap(aInd, pos, i+pos);                       		// swap the individual
        }
    }
}

// ----------------------------------------------------------------------------------------
// remove::
// ----------------------------------------------------------------------------------------
/* remove the specified individual from the list and adjust the rest
 * this can only be done if a subset is not randomly drawn
 */
void
TPatchFitness::remove(unsigned int i){
    assert(i<_nbInd);
    assert((_sort>=-2 && _sort<=2) || _sort==10);
    
    // remove the individual
    if(_sort == 10){	// fitness array is not cumulative
        ++i;
        for(; i<_nbInd; ++i){
            _aFit[i-1] = _aFit[i];
            _aInd[i-1] = _aInd[i];
        }
    }
    else {          	// fitness array is cumulative
        double diff = i ? _aFit[i]-_aFit[i-1] : _aFit[0];	// get the fitness of the individual
        ++i;
        for(; i<_nbInd; ++i){
            _aFit[i-1] = _aFit[i]-diff;
            _aInd[i-1] = _aInd[i];
        }
    }
    
    --_nbInd;	// decrement the pop size by 1
}

//////////////////////////////////////////////////////////////////////////////
// functions to get an individual depending on its fitness
// random not taking into account the fitness
TIndividual* TPatchFitness::get_RAND_noFit(){                                          // of all individuals
    return _aInd[_pop->rand().Uniform(_nbInd)];
}
TIndividual* TPatchFitness::get_RAND_noFit_index(unsigned int& index){                                          // of all individuals
    index = _pop->rand().Uniform(_nbInd);
    return _aInd[index];
}
TIndividual* TPatchFitness::get_RAND_noFit_subset(const int& nbInd){                   // of a subset of nbInd individuals
    return _aInd[_pop->rand().Uniform(_nbInd)];
}

// most fittest
TIndividual* TPatchFitness::get_mostFit(){                                             // get the fittest
    assert(_sort==-2);
    return _aInd[0];
}
TIndividual* TPatchFitness::get_RAND_mostFit(){                                        // get randomly the fittest
    assert(_sort<0);
    return _aInd[_pop->rand().AfterDistribution(_aFit, _nbInd)];
}
TIndividual* TPatchFitness::get_RAND_mostFit_index(unsigned int& index){               // "returns" the index of the individual
    assert(_sort<0);
    index = _pop->rand().AfterDistribution(_aFit, _nbInd);
    return _aInd[index];
}
TIndividual* TPatchFitness::get_RAND_mostFit_of_mostFit(const unsigned int& nb){       // fixed subset
    assert(_sort==-2);
    if(nb>_nbInd) return _aInd[_pop->rand().AfterDistribution(_aFit, _nbInd)];
    return               _aInd[_pop->rand().AfterDistribution(_aFit, nb)];
}
TIndividual* TPatchFitness::get_RAND_mostFit_of_RAND_mostFit(const unsigned int& nb){  // random subset
    assert(_sort==-3);
    if(nb>_nbInd) return _aInd[_pop->rand().AfterDistribution(_aFit, _nbInd)];
    assert(nb==_nbSubset);
    return _aInd[_pop->rand().AfterDistribution(_aFit, nb)];
}

// less fittest
TIndividual* TPatchFitness::get_lessFit(){                                             // get the less fittest
    assert(_sort==2);
    return _aInd[0];
}
TIndividual* TPatchFitness::get_RAND_lessFit(){                                        // get randomly the less fittest
    assert(_sort>0);
    return _aInd[_pop->rand().AfterDistribution(_aFit, _nbInd)];
}
TIndividual* TPatchFitness::get_RAND_lessFit_index(unsigned int& index){               // "returns" the index of the individual
    assert(_sort>0);
    index = _pop->rand().AfterDistribution(_aFit, _nbInd);
    return _aInd[index];
}
TIndividual* TPatchFitness::get_RAND_lessFit_of_lessFit(const unsigned int& nb){      // fixed subset
    assert(_sort==2);
    if(nb>_nbInd) return _aInd[_pop->rand().AfterDistribution(_aFit, _nbInd)];
    return               _aInd[_pop->rand().AfterDistribution(_aFit, nb)];
}
TIndividual* TPatchFitness::get_RAND_lessFit_of_RAND_lessFit(const unsigned int& nb){  // random subset
    assert(_sort==3 && nb==_nbSubset);
    if(nb>_nbInd) return _aInd[_pop->rand().AfterDistribution(_aFit, _nbInd)];
    assert(nb==_nbSubset);
    return _aInd[_pop->rand().AfterDistribution(_aFit, _nbInd)];
}


