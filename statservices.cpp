/** @file statservices.cpp
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


#include "statservices.h"
#include <sstream>
#include <list>
#include "stathandler.h"
#include "metapop.h"
#include "functions.h"
#include "tsimulation.h"
#include "treplicate.h"
#include "tstat_db.h"
#include "param.h"

//-----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
bool StatServices::init()
{
    _pStat_db = _popPtr->get_pSimulation()->stats;
    if(!_pStat_db->get_tot_occurrence()) return false;

    if(!init_paramArgs() & !init_statArgs()) return false;       // no stats are computed (both have to be computed: & and not &&!!!)
    
    // write the heading of the stat file if the stats are outputted regularerly
    if(_pStat_db->get_save_choice() == 1){    // every replicate: the output is directly written to the file (no db)
#ifdef _DEBUG
        message("  LCE_StatFH::FHwrite (%s)\n",_pStat_db->get_file_name_stats().c_str());
#endif
        assert(!file_stat);
        file_stat = new ofstream(_pStat_db->get_file_name_stats().c_str(),ios::out);
        if(!(*file_stat)) error("Could not open stat output file \"%s\"!\n",_pStat_db->get_file_name_stats().c_str());
        
        _pStat_db->printStatHeaders(*file_stat, FLAT | PARAM);     	//print the stat names
    }
    
    // compute the age classes used to compute the stats
    set_age_class();
    
    return true;
}


//-----------------------------------------------------------------------------------------
// init_statArgs
// ----------------------------------------------------------------------------------------
bool
StatServices::init_statArgs()
{
    if(_statArg.empty()) return false;
    list< StatHandlerBase* >::iterator HIT;
    
#ifdef _DEBUG
    message("Loaded stats:\n");
    message("  init StatHandler (");
    for(HIT = _children.begin(); HIT != _children.end(); ++HIT) {
        message(" %s", (*HIT)->getName().c_str());
    }
    message(")\n");
#endif
    
    //first init the stat handlers -> will define the dims of the stat tables
    for(HIT = _children.begin(); HIT != _children.end(); ++HIT) {
        (*HIT)->set_popPtr(_popPtr);
        (*HIT)->init();
    }
    
    // for each stat of the input file
    TMatrixVar<string> m(_statArg);
    unsigned int i, j, s1, s2;
    string statName;
    bool is_set, global_set=false;
    for(i=0, s1=m.getNbRows(); i<s1; ++i){
        for(j=0, s2=(unsigned int)m.get(i)->size(); j<s2; ++j){
            is_set = false;
            statName = m.get(i,j);      // get the next stat key word
            
#ifdef _DEBUG
            message("  %s (", statName.c_str());  // find the corresponding stat handler
#endif
            
            for(HIT = _children.begin(); HIT != _children.end(); ++HIT) {
                if(!(*HIT)->setStatRecorders(statName)) continue;
                is_set = global_set = true;
            }
            
#ifdef _DEBUG
            message(")\n");
#endif
            
            if(!is_set) {
                warning("The string \"%s\" is not a valid statistic option and is therfore not considered!\n", statName.c_str());
            }
        }
    }
    return global_set;
}

//-----------------------------------------------------------------------------------------
// init_paramArgs
// ----------------------------------------------------------------------------------------
bool
StatServices::init_paramArgs()
{
    if(_paramArg.empty()) return false;
    
#ifdef _DEBUG
    message("Loaded params:\n");
#endif

    // for each param of the input file
    assert(_paramList.empty());
    list<ParamSet*>* lParams = &_popPtr->get_pReplicate()->get_allParams();
    list< ParamSet* >::iterator HIT;
    Param* curParam;
    unsigned int i, j, s1, s2;
    string paramName;
    bool global_set=false;
    TMatrixVar<string> m(_paramArg);
    for(i=0, s1=m.getNbRows(); i<s1; ++i){
        for(j=0, s2=(unsigned int)m.get(i)->size(); j<s2; ++j){
            paramName = m.get(i,j);      // get the next param key word
            
#ifdef _DEBUG
            message("  %s", paramName.c_str());
#endif
       
            for(HIT = lParams->begin(); HIT != lParams->end(); ++HIT) {
                curParam = (*HIT)->find_param(paramName);
                if(!curParam) continue;
                _paramList.push_back(curParam);
                global_set = true;
                break;
            }
                        
            if(HIT==lParams->end()) {
                warning("The string \"%s\" is not a valid parameter option for the output and is therfore not considered!\n", paramName.c_str());
            }
        }
    }
    
    // now we have the stats : generate now the corresponding arrays to store the arguments
    StatRecBaseAll* recBase;
    assert(!_paramStats);
    ARRAY::create_1D<string*>(_paramStats, (unsigned int)_paramList.size(), NULL); // array for each stat
    list<Param*>::iterator curList=_paramList.begin(), endList=_paramList.end();
    for(unsigned int i=0; curList!=endList; ++curList, ++i){
        
        if(_popPtr->_current_replicate==my_NAN) _pStat_db->add_paramRecBase(*curList);
        else{
            recBase = *(_popPtr->_current_statRecBase++);
            assert(recBase->_name == (*curList)->get_name());
            
            assert(!_paramStats[i]);
            ARRAY::create_1D(_paramStats[i], _pStat_db->get_tot_occurrence());
            recBase->set_param_pointer(_paramStats[i], _popPtr->_current_replicate);
        }
    }
    
    return global_set;
}

//-----------------------------------------------------------------------------------------
// get_nb_stats
// ----------------------------------------------------------------------------------------
unsigned int
StatServices::get_nb_stats (){
    unsigned int nb=0;
    list< StatHandlerBase* >::iterator pos, end;
    for(pos=_children.begin(), end=_children.end(); pos!=end; ++pos) {
        nb +=(*pos)->getNbRecorders();
    }
    return nb;
}

//-----------------------------------------------------------------------------------------
// load
// ----------------------------------------------------------------------------------------
void StatServices::load ( SimComponent* sc ){
    sc->loadStatServices(this);
}

//-----------------------------------------------------------------------------------------
// attach
// ----------------------------------------------------------------------------------------
void StatServices::attach ( StatHandlerBase* SH){
    Service::attach(SH);
    _children.push_back(SH);
    SH->set_service(this);
}


//-----------------------------------------------------------------------------------------
// set_age_class
// ----------------------------------------------------------------------------------------
void StatServices::set_age_class()
{
    _age = NONE;
    list< StatHandlerBase* >::iterator pos, end;
    for(pos=_children.begin(), end=_children.end(); pos!=end; ++pos) {
        _age |= (*pos)->get_age_class();
    }
}

//-----------------------------------------------------------------------------------------
// notify
// ----------------------------------------------------------------------------------------
void StatServices::notify ()
{
    if(file_stat){          // if stats are directly written to file (file remains open during the entire simulation)
        //	ofstream FH(_file_stats.c_str(), ios::app);
        //	if (!FH) error("Could not open stat output file \"%s\"!\n", _file_stats.c_str());
       // _popPtr->printGenRep2File(*file_stat); // already done with the params
        Service::notify(*file_stat);
        *file_stat << endl;       // flush the stream...
        //	FH.close();
    }
    else Service::notify();              // if the stats are stored in db
}


//-----------------------------------------------------------------------------------------
// notify
// ----------------------------------------------------------------------------------------
void StatServices::update_states ()
{
    _popPtr->set_sampledInds(_age);
    _popPtr->update_patch_states();
}

//-----------------------------------------------------------------------------------------
// notify
// ----------------------------------------------------------------------------------------
void StatServices::notify_params ()
{
    if(_paramList.empty()) return;
    
    if(file_stat){          // if stats are directly written to file (file remians open during the entire simulation)
        //	ofstream FH(_file_stats.c_str(), ios::app);
        //	if (!FH) error("Could not open stat output file \"%s\"!\n", _file_stats.c_str());
        _popPtr->printGenRep2File(*file_stat);
        
        
        assert(_paramStats);
        list<Param*>::iterator curParam=_paramList.begin(), endParam=_paramList.end();
        for(; curParam!=endParam; ++curParam){
            *file_stat<<"\t";
            (*file_stat).width(12);
            *file_stat << (*curParam)->get_arg();
        }
        //	FH.close();
    }
    else {        // if the stats are stored in db
        list<Param*>::iterator curParam=_paramList.begin(), endParam=_paramList.end();
        for(unsigned int i=0; curParam!=endParam; ++curParam, ++i){
            assert(*curParam);
            string tt = (*curParam)->get_arg();
            _paramStats[i][get_pop_ptr()->_current_index_stat_db] = (*curParam)->get_arg();
        }
    }
}

//-----------------------------------------------------------------------------------------
// reset
// ----------------------------------------------------------------------------------------
void StatServices::reset (){
    Service::reset_observers();
    _children.clear(); // StatHandlerBase themself are deleted by the corresponding traits
}



