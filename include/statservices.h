/** @file statservices.h
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

#ifndef statservicesH
#define statservicesH

#include "service.h"




class TMetapop;
class TStat_db;
class Param;

class StatHandlerBase;
/**The Service class used to manage the StatHandler objects.*/
class StatServices : public Service {
    
private:
    
    TMetapop* _popPtr;
    TStat_db* _pStat_db;
    list< StatHandlerBase* > _children;
    string _statArg;
    string _paramArg;
    
    // parameters to store the arguments to output
    list<Param*> _paramList;      // contains a pointer to the param object
    string** _paramStats;         //  _paramStats[stat][generation] same order as in _paramList (delete only outer array)

    
    age_t         _age;
    
    string        _NaN;         // placeholder for NaN or not available
    
    ofstream* file_stat;
    
public:
    
    typedef list< StatHandlerBase* >::const_iterator stat_it;
    
    StatServices ( ) : _popPtr(0), _paramStats(0), _pStat_db(0), _age(0), file_stat(0){ }
    
    virtual ~StatServices ( ) {
        if(file_stat){file_stat->close(); delete file_stat;}
        if(_paramStats) {delete[] _paramStats;}
    }
    
    virtual bool init ( );
    virtual bool init_statArgs ( );
    virtual bool init_paramArgs ( );
    
    TMetapop* get_pop_ptr          ( )      {return _popPtr;}
    
    void set_pop_ptr              (TMetapop* pop) {_popPtr=pop;}
    
    void set_statArg  (string& str) {_statArg = str;}
    void set_paramArg (string& str) {_paramArg = str;}
    
    void set_age_class();

    unsigned int get_nb_stats();
    unsigned int get_nb_params(){return (unsigned int)_paramList.size();}
    
    stat_it getFirst () {return _children.begin();}
    
    stat_it getLast () {return _children.end();}
    
    //   unsigned int getStatRecIndex  (unsigned int i);
    /*
     *
     */
    virtual void update_states ();      // update the sampled inds
    virtual void notify ();             // output the stats
    virtual void notify_params ();      // output the params
    
    /**
     *  tell the SimComponent to load its stat handlers
     *  @param sc the SimComponent
     */
    virtual void load ( SimComponent* sc );
    
    /**
     *  attach the StatHandler to the current list (_children) of the StatServices
     *  @param SH the StatHandler
     */
    virtual void attach ( StatHandlerBase* SH);
    
    /**
     *  clear the list of StatHandler
     */
    virtual void reset ( );
};


#endif //STATSERVICES_H

