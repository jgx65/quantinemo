/** @file stat_rec_base.cpp
 *
 *   Copyright (C) 2006 Frederic Guillaume    <guillaum@zoology.ubc.ca>
 *   Copyright (C) 2008 Samuel Neuenschwander <samuel.neuenschwander@unil.ch>
 *   Copyright (C) 2018 Frederic Michaud <frederic.michaud@unil.ch>
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


#include "stat_rec_base.h"
#include "tmetapop.h"
#include "tstat_db.h"
#include "lce_misc.h"
#include "lce_coalescence_base.h"
#include "ttneutral.h"
#include "ttquanti.h"

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/
//
//                               ****** StatRecBase ******
// ----------------------------------------------------------------------------------------
// StatRecBase::~StatRecBase()
// ----------------------------------------------------------------------------------------
StatRecBase::~StatRecBase()
{
    //  if(_val){
    //    delete [] _val;
    //  }
}

// ----------------------------------------------------------------------------------------
// StatRecBase::init()
// ----------------------------------------------------------------------------------------
void StatRecBase::init()
{
    _title = "";
    _name = "";
    _age = ALL;
    _ordering = FLAT;
    _arg = 0;
}

// ----------------------------------------------------------------------------------------
// StatRecBase::set()
// ----------------------------------------------------------------------------------------
void StatRecBase::set(string T, string N, st_order Order, age_t AGE, unsigned int ARG)
{
    _title = T;
    _name = N;
    _ordering = Order;
    _age = AGE;
    _arg = ARG;
}


// ----------------------------------------------------------------------------------------
// StatRecorder::init_stat_rec
// ----------------------------------------------------------------------------------------
/** initialize StatRecorder and link it to the Recorder object of main db */
template<class S>
void StatRecorder<S>::init_stat_rec(TMetapop* pop)
{
    if(pop->_current_replicate==my_NAN) _recAll = _pStat_db->add_statRecBase(this); // initialize the stat
    else{
        _recAll = *pop->_current_statRecBase++;
        assert(_recAll->_name == getName());
        
        assert(!_stat);
        ARRAY::create_1D(_stat, _pStat_db->get_tot_occurrence(), (double)my_NAN);
        _recAll->set_stat_pointer(_stat, pop->_current_replicate);
    }
}


// ----------------------------------------------------------------------------------------
// StatRecorder::StatRecorder()
// ----------------------------------------------------------------------------------------
template <class S>
StatRecorder<S>::StatRecorder(TStat_db* db): _getStat(0),    _getStatBool(0), _getStatUI(0),
_getStatAGE(0), _getStatBoolAGE(0), _getStatUIAGE(0), _stat(0), _recAll(0)
{
    _pStat_db = db;
    
    //init the params, set all table cells to "NaN"
    init();
}

// ----------------------------------------------------------------------------------------
// StatRecorder::setStats
// ----------------------------------------------------------------------------------------
/** store te stats in the db */
template<class S>
void StatRecorder<S>::setStats(age_t AGE, int crnt_gen, int rpl_cntr, S* StatHandler, TMetapop* pop)
{
    double statValue=0;
    age_t age = getAge();
    
    // test if the stat has to be computed for this age
    if(age & AGE) {
#ifdef _DEBUG
        message(" %s\n",getName().c_str());
#endif
        
        //get the value:
        if(_getStat)           		statValue = (StatHandler->*_getStat)       ();
        else if(_getStatBool)  		statValue = (StatHandler->*_getStatBool)   ((bool)getArg());
        else if(_getStatUI)    		statValue = (StatHandler->*_getStatUI)     (getArg());
        else if(_getStatAGE)      	statValue = (StatHandler->*_getStatAGE)    (age);
        else if(_getStatBoolAGE)  	statValue = (StatHandler->*_getStatBoolAGE)((bool)getArg(), age);
        else if(_getStatUIAGE)    	statValue = (StatHandler->*_getStatUIAGE)  (getArg(), age);
        else error("StatRecorder::setVal2db: no _getStat fct ptr !!\n");
        
        assert(_stat);
        
        _stat[pop->_current_index_stat_db] = statValue;
    }
}

// ----------------------------------------------------------------------------------------
// needed since compiler has to know which classes are templated.
// ----------------------------------------------------------------------------------------
template class StatRecorder<LCE_CoalescenceSH>;
template class StatRecorder<LCE_StatSH>;
template class StatRecorder<MetapopSH>;
template class StatRecorder<TTNeutralSH>;
template class StatRecorder<TTQuantiSH>;
