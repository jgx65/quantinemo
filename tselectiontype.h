/** @file tselectiontype.h
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

#ifndef tselectiontypeH
#define tselectiontypeH
//---------------------------------------------------------------------------
#include "tselectiontrait.h"
#include "tselection.h"


class TSelection;

class TSelectionNeutral: public TSelectionTrait
{
public:
    double get_fitness_none();
    
    TSelectionNeutral(TSelection* s, const int& t){
        TSelectionTrait::init(s, t);
        init();
    }
    ~TSelectionNeutral(){}
    void init();
};

class TSelectionStabilizing: public TSelectionTrait
{
public:
    double get_fitness_none();
    
    TSelectionStabilizing(TSelection* s, const int& t){
        TSelectionTrait::init(s, t);
        init();
    }
    ~TSelectionStabilizing(){}
    void init();
};

class TSelectionDirectional: public TSelectionTrait
{
public:
    double get_fitness_none();
    
    TSelectionDirectional(TSelection* s, const int& t){
        TSelectionTrait::init(s, t);
        init();
    }
    ~TSelectionDirectional(){}
    void init();
};

class TSelectionFitnessLandscape: public TSelectionTrait
{
public:
    double get_fitness_none();
    
    TSelectionFitnessLandscape(TSelection* s, const int& t){
        TSelectionTrait::init(s, t);
        init();
    }
    ~TSelectionFitnessLandscape(){}
    void init();
};

class TSelectionSelectionCoefficient: public TSelectionTrait
{
public:
    double get_fitness_none();
    
    TSelectionSelectionCoefficient(TSelection* s, const int& t){
        TSelectionTrait::init(s, t);
        init();
    }
    ~TSelectionSelectionCoefficient(){}
    void init();
};


#endif
