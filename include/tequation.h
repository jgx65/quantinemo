/** @file tequation.h
 *
 *   Copyright (C) 2008 Samuel Neuenschwander <samuel.neuenschwander@unil.ch>
 *   Copyright (C) 2018 Frederic Michaud <frederic.michaud@unil.ch>

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


#ifndef tequationH
#define tequationH

#include <vector>
#include "tindividual.h"
using namespace std;

class TEquation{
public:
    // containers
    vector<double (TEquation::*)(double val1, double val2)> arithmetic_vec;        // has one less argument
    vector<double (TEquation::*)(TIndividual* ind, double val)> values_vec;    // G: genotpyic value, Z: phenotype, C: constant
    vector<double> index_vec;          // the trait index or a constant
    unsigned int nb_values;         // number of values (values.size() & quantity.size()
    
    // iterators
    vector<double (TEquation::*)(double val1, double val2)>::iterator cur_arithmetic, end_arithmetic;
    vector<double (TEquation::*)(TIndividual* ind, double val)>::iterator cur_value;
    vector<double>::iterator cur_index;
    
    TEquation(string input, vector<int> qtrait){
        read_input(input, qtrait);
    }

    void read_input(string input, vector<int> qtrait);
    double getValue(TIndividual* ind);
    void add_index(int index, string input, vector<int> qtrait);
    
    // sign
    double (TEquation::*arithmetic)(double val1, double val2);
    double plus (double val1, double val2){return val1 + val2;}
    double minus(double val1, double val2){return val1 - val2;}
    double mult (double val1, double val2){return val1 * val2;}
    double div  (double val1, double val2){return val1 / val2;}
   
    // source
    double (TEquation::*getValuePos)(TIndividual* ind, double val);
    double getPheno(TIndividual* ind, double val){return ind->getTraitPhenotype((unsigned int)val);}
    double getGeno (TIndividual* ind, double val){return ind->getTraitGenotype((unsigned int)val);}
    double getConst(TIndividual* ind, double val){return val;}
    
};

#endif /* defined(__quantiNemo2__tequation__) */
