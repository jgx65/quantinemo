/** @file stathandlerbase.cpp
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

#include "stathandlerbase.h"
#include "tmetapop.h"
using namespace std;



// ----------------------------------------------------------------------------------------
// functions to get patches
// ----------------------------------------------------------------------------------------
// all patch (empty or populated, i corresponds to the id of the patch
vector<Patch*>& StatHandlerBase::get_vPatch() {return _popPtr->get_vPatch();}
Patch*          StatHandlerBase::get_vPatch(const unsigned int& i) {return _popPtr->get_vPatch(i);}
unsigned int    StatHandlerBase::get_nbPatch() {return _popPtr->get_nbPatch();}

// sampled AND currently populated patches (i does not correspond to anything)
unsigned int    StatHandlerBase::get_nbTotSamplePatch() {return _popPtr->get_nbTotSamplePatch();}
vector<Patch*>& StatHandlerBase::get_vSamplePatch() {return _popPtr->get_vSamplePatch();}


// ----------------------------------------------------------------------------------------
// fucntions to get the number of populated patches for the specified sex and/or age
// ----------------------------------------------------------------------------------------
unsigned int
StatHandlerBase::get_nbSamplePatch(const sex_t& SEX)
{
    unsigned int size=0;
    vector<Patch*>::iterator curPop = get_vSamplePatch().begin();
    vector<Patch*>::iterator endPop = get_vSamplePatch().end();
    for(; curPop != endPop; ++curPop){
        if((*curPop)->size(SEX, ALL)) ++size;
    }
    return size;
}

unsigned int
StatHandlerBase::get_nbSamplePatch(const sex_t& SEX, const age_idx& AGE)
{
    unsigned int size=0;
    vector<Patch*>::iterator curPop = get_vSamplePatch().begin();
    vector<Patch*>::iterator endPop = get_vSamplePatch().end();
    for(; curPop != endPop; ++curPop){
        if((*curPop)->size(SEX, AGE)) ++size;
    }
    return size;
}

unsigned int
StatHandlerBase::get_nbSamplePatch(const sex_t& SEX, const age_t& AGE)
{
    unsigned int size=0;
    vector<Patch*>::iterator curPop = get_vSamplePatch().begin();
    vector<Patch*>::iterator endPop = get_vSamplePatch().end();
    for(; curPop != endPop; ++curPop){
        if((*curPop)->size(SEX, AGE)) ++size;
    }
    return size;
}

unsigned int
StatHandlerBase::get_nbSamplePatch(const age_idx& AGE)
{
    unsigned int size=0;
    vector<Patch*>::iterator curPop = get_vSamplePatch().begin();
    vector<Patch*>::iterator endPop = get_vSamplePatch().end();
    for(; curPop != endPop; ++curPop){
        if((*curPop)->size(AGE)) ++size;
    }
    return size;
}

unsigned int
StatHandlerBase::get_nbSamplePatch(const age_t& AGE)
{
    unsigned int size=0;
    vector<Patch*>::iterator curPop = get_vSamplePatch().begin();
    vector<Patch*>::iterator endPop = get_vSamplePatch().end();
    for(; curPop != endPop; ++curPop){
        if((*curPop)->size( AGE)) ++size;
    }
    return size;
}

// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
bool StatHandlerBase::init()
{
    _pStat_db = _popPtr->get_pStat_db();
    
    return true;
}

//-----------------------------------------------------------------------------------------
// set_age_class
// ----------------------------------------------------------------------------------------
age_t StatHandlerBase::get_age_class()
{
    age_t age = NONE;
    list< StatRecBase* >::iterator pos, end;
    for(pos=_stats.begin(), end=_stats.end(); pos!=end; ++pos) {
        age |= (*pos)->getAge();
    }
    return age;
}

// ----------------------------------------------------------------------------------------
// update
// ----------------------------------------------------------------------------------------
void StatHandlerBase::update()
{
    execute();
}

// ----------------------------------------------------------------------------------------
// update
// ----------------------------------------------------------------------------------------
void StatHandlerBase::update(ostream& FH)
{
    execute(FH);
}
// ----------------------------------------------------------------------------------------
// reset
// ----------------------------------------------------------------------------------------
void StatHandlerBase::reset ( )
{
    for(list< StatRecBase* >::iterator IT = _stats.begin(); IT != _stats.end(); ++IT) {
        delete (*IT);
    }
    
    _stats.clear();
    clear();
}

// ----------------------------------------------------------------------------------------
// get_current_nbSamplePatch
// ----------------------------------------------------------------------------------------
unsigned int
StatHandlerBase::get_current_nbSamplePatch()
{
    assert(_popPtr);
    return _popPtr->_current_nbSamplePatch;
}

// ----------------------------------------------------------------------------------------
// get_current_replicate
// ----------------------------------------------------------------------------------------
unsigned int
StatHandlerBase::get_current_replicate()
{
    assert(_popPtr);
    return _popPtr->_current_replicate;
}

// ----------------------------------------------------------------------------------------
// get_last_nbSamplePatch
// ----------------------------------------------------------------------------------------
unsigned int
StatHandlerBase::get_last_nbSamplePatch()
{
    assert(_popPtr);
    return _popPtr->_last_nbSamplePatch;
}

// ----------------------------------------------------------------------------------------
// get_current_index_stat_db
// ----------------------------------------------------------------------------------------
unsigned int
StatHandlerBase::get_current_index_stat_db()
{
    assert(_popPtr);
    return _popPtr->_current_index_stat_db;
}

// ----------------------------------------------------------------------------------------
// get_current_generation
// ----------------------------------------------------------------------------------------
unsigned int
StatHandlerBase::get_current_generation()
{
    assert(_popPtr);
    return _popPtr->getCurrentGeneration();
}

// ----------------------------------------------------------------------------------------
// getStatRecIndex
// ----------------------------------------------------------------------------------------
/*
unsigned int StatHandlerBase::getStatRecIndex (unsigned int i)
{
    if(_stats.empty()) return my_NAN;
    assert(i <= (*_stats.begin())->getIndexSize());
    
    //find the first statRecorder that has computed some stats
    for(STAT_IT el=_stats.begin(); el != _stats.end(); ++el){
        if((*el)->getIndex(i) != my_NAN) return (*el)->getIndex(i);
    }
    return my_NAN;
}
*/
