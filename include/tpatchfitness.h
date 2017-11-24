/** @file tpatchfitness.h
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
//---------------------------------------------------------------------------

#ifndef tpatchfitnessH
#define tpatchfitnessH
//---------------------------------------------------------------------------

class Individual;
class TSelection;
class TMetapop;


class TPatchFitness{
    friend class TSelection;
private:
    Individual**  _aInd;      // an array of pointers to the individulas (don't delete the individuals here!)
    double*       _aFit;
    unsigned int  _nbInd;     // USED size of the array (i.e. number of individuals)
    unsigned int  _nbSubset;   // if a subset is used
    unsigned int  _capacity;  // size of the array
    int           _sort;      // 10: the fitnesses were never used before (not sorted nor made cumulative)
    // -3: decremental sort (highest fitness first) subset random most fittest
    // -2: decremental sort (highest fitness first) subset fix most fittest
    // -1: decremental sort (highest fitness first) -> normal cumulative
    //  0: no sort                                  -> normal cumulative
    //  1: incremental sort (lowest fitness first)  -> reverse cumulative
    //  2: incremental sort (lowest fitness first) subset fix less fittest
    //  3: incremental sort (lowest fitness first) subset random less fittest
    
    TMetapop* _pop;
    
    void init(TMetapop* ptr, int size);
    void resize(const unsigned int& size);
    
    // constructors
    TPatchFitness(){init(0, 0);}
    TPatchFitness(TMetapop* ptr, const unsigned int& size){init(ptr, size);}
    
    // destructor
    ~TPatchFitness(){
        if(_aInd) delete[] _aInd;
        if(_aFit) delete[] _aFit;
    }
    
   
    //////////////////////////////////////////////////////////////////////////////
    // functions to get an individual depending on its fitness
    // random not taking into account the fitness
    Individual* get_RAND_noFit();                                          // of all individuals
    Individual* get_RAND_noFit_index(unsigned int& index);                                          // of all individuals
    Individual* get_RAND_noFit_subset(const int& nbInd);                   // of a subset of nbInd individuals
    
    // most fittest
    Individual* get_mostFit();                                             // get the fittest
    Individual* get_RAND_mostFit();                                        // get randomly the fittest
    Individual* get_RAND_mostFit_index(unsigned int& index);               // "returns" the index of the individual
    Individual* get_RAND_mostFit_of_mostFit(const unsigned int& nb);       // fixed subset
    Individual* get_RAND_mostFit_of_RAND_mostFit(const unsigned int& nb);  // random subset

    // less fittest
    Individual* get_lessFit();                                            // get the less fittest
    Individual* get_RAND_lessFit();                                        // get randomly the less fittest
    Individual* get_RAND_lessFit_index(unsigned int& index);               // "returns" the index of the individual
    Individual* get_RAND_lessFit_of_lessFit(const unsigned int& nb);      // fixed subset
    Individual* get_RAND_lessFit_of_RAND_lessFit(const unsigned int& nb);  // random subset
    
    void sort(int sort, unsigned int nbInd=0);
    void randomize_order_noFit  (double* aFit, Individual** aInd, unsigned int& size, unsigned int nbSubset=0);
    void randomize_order_mostFit(double* aFit, Individual** aInd, unsigned int& size, unsigned int nbSubset=0);
    void randomize_order_lessFit(double* aFit, Individual** aInd, unsigned int& size, unsigned int nbSubset=0);
    
    
    double getMeanFitness();
    double getSumFitness();
    
    void remove(unsigned int i);
    
    Individual** get_aInd() {return _aInd;}
    double*      get_aFit() {return _aFit;}
};



#endif
