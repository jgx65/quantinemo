/** @file stathandler.h
 *
 *   Copyright (C) 2006 Frederic Guillaume    <guillaum@zoology.ubc.ca>
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

#ifndef stathandlerbaseH
#define stathandlerbaseH

#include <list>
#include <map>
#include "handler.h"
#include "statservices.h"
using namespace std;

class TMetapop;
class TPatch;
class TStat_db;

class StatRecBase;
/**Base class of the StatHandler class, implements the Handler interface.
 * This class stores the Handler state and the list of links to the StatRecBase.
 * This list is duplicated in the StatHandler as a list of the StatRecorder templates.
 * It allows the StatService to access the stat recorders without knowledge of the
 * actual StatHandler and StatRecorder template type.
 */
class StatHandlerBase : public Handler {
    
private:
    
    list<StatRecBase*> _stats;           	/**Link to the StatRecorder list elements in the StatHandler derived class.*/
    TStat_db* _pStat_db;
    map<unsigned int, unsigned int> *_index;           // don't delete it
    StatServices*      _service;          /**Link to the StatService.*/
    unsigned int       _GenerationOccurrence;  /**Occurence of the stats recording as set by the user (see parameter "stat_log_time").*/
    
protected:
    
public:
    StatHandlerBase( ) {}
    
    virtual           ~StatHandlerBase     ( ) {  }
    
    ///@name Accessors
    ///@{
    TMetapop*          get_pop_ptr_popHandler( )            {return _popPtr;}
    void              set_service      (StatServices* srv) {_service = srv;}
    StatServices*     get_service      ( )                 {return _service;}
    unsigned int      getNbRecorders   ( )                 {return (unsigned int)_stats.size();}
    
    list<StatRecBase*>& getStats       ( )                 {return _stats;}
    virtual void      add              (StatRecBase* rec)  {_stats.push_back(rec);}
    
    unsigned int get_current_nbSamplePatch();
    unsigned int get_last_nbSamplePatch();
    unsigned int get_current_index_stat_db();
    unsigned int get_current_generation();
    unsigned int get_current_replicate();

    ///@}
    /**Empties the _stats list and calls clear() (defined in the derived class).*/
    virtual void      reset            ( );
    
    // all sampled patch (empty or populated)
    vector<TPatch*>&   get_vPatch();
    TPatch*            get_vPatch(const unsigned int& i);
    unsigned int      get_nbPatch();
    unsigned int      get_nbTotSamplePatch();
    
    // functions to get the total populated number of patches for the specified sex and/or age class
    unsigned int      get_nbSamplePatch(const sex_t& SEX);
    unsigned int      get_nbSamplePatch(const sex_t& SEX, const age_idx& AGE);
    unsigned int      get_nbSamplePatch(const sex_t& SEX, const age_t& AGE);
    unsigned int      get_nbSamplePatch(const age_idx& AGE);
    unsigned int      get_nbSamplePatch(const age_t& AGE);
    
    
    // just sampled and populated patches
    vector<TPatch*>&   get_vSamplePatch();
    unsigned int      get_nbSamplePatch();
    
    ///@name Stat recorder interface
    ///@{
    //	unsigned int      getStatRecIndex  (unsigned int i);
    
    age_t             get_age_class();
    
    void              setMean       ( );
    
    ///@}
    ///@name Handler implementation
    ///@{
    virtual bool      init             ( );
    virtual void      update           ( );
    virtual void      update           (ostream& FH);
    ///@}
    ///@name StatHandler interface declaration
    ///@{
    virtual void      execute          ( ){};
    virtual void      execute          (ostream& FH){}
    
    virtual bool      setStatRecorders (const string& token) = 0;
    
    virtual string getName() = 0;
    
    virtual void      clear            ( ) = 0;
    
    void set_popPtr(TMetapop* pop){_popPtr=pop;}
    ///@}
};

#endif //STATHANDLERBASE_H

