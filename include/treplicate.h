/** @file treplicate.h
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

#ifndef treplicateH
#define treplicateH


#include "basicsimulation.h"

using namespace std;

class TSimulation;
class TStat_db;

class TReplicate: public SimBuilder{
private:
    map<string, string> _params;
    
    unsigned int _isCoalescence; // 0: individual base; 1: coal with speed optimization; 2: coal with meory optimization
    
    TSimulation* _pSimulation;  // don't delete
    TStat_db*    _pStat_db;     // don't delete
    
    TMetapop* _thePop;
    
    /** time parameters */
    char      _startTime[20], _endTime[20];
    clock_t   _meanElapsedTimeRepl;
    
    unsigned int _current_thread;       // the thrread the replciate is running on (to pick the random engine)
    unsigned int _nbThreads;            // the nubmer of threads for this replciate (not used yet)

    
public:
    TReplicate(TSimulation* ptr, map<string, string> params, map<string, string> keys,
               unsigned int replicate, unsigned int curThread, unsigned int nbThreads);
    TReplicate(TSimulation* ptr);
    ~TReplicate();
    
    RAND* rand;     // get the actual random generator from TSimulation (don't delete it)
    
    
    unsigned int _current_replicate;    // the current replciate.
    unsigned int _generations;          // Number of generations to run (used for db).

    
    unsigned int get_nbTraits(string name, map<string, string>& inputParams);
    unsigned int isCoalescence(map<string, string>& inputParams);
    void init_paramset();
    
    void loadDefaultTemplates(map<string, string>& inputParams);
    
    unsigned int get_isCoalescence(){return _isCoalescence;}
    

    void generate_LCE_individual();
    void generate_LCE_coalescence();
    void generate_traits(map<string, string>& inputParams);
    
    bool run_replicate(map<string, string>& params, map<string, string>& keys,
                       unsigned int currentReplicate);
    bool run_replicate_ind(map<string, string>& params, map<string, string>& keys,
                       unsigned int currentReplicate);
    bool run_replicate_coal(map<string, string>& params, map<string, string>& keys,
                       unsigned int currentReplicate);
    
    unsigned int get_current_replicate() {return _current_replicate;}
    
    FileServices* _FileServices;
    StatServices* _StatServices;
    list<ParamSet*>* _allParams;
    
    void print_info(TMetapop* thePop);
    bool setup(map<string, string>& simparams, map<string, string>& simkeaywords, TMetapop* thePop);
    void register_services(TSimComponent* cmpt);
    bool build_pop(TMetapop* thePop);
    
    void register_all(TMetapop* thePop);
    void build_stat_recorders(TMetapop* thePop);
    
    clock_t print_start_replicate(unsigned int currentReplicate);
    
    TSimulation* get_pSimulation(){return _pSimulation;}
    TStat_db*    get_pStat_db(){return _pStat_db;}
    
    bool test_replicate_and_setUpStats(map<string, string>& params, map<string, string>& keays);
    void print_help(ostream& os, unsigned int wide1=20, unsigned int wide2=60,
                    char fill='.', unsigned int importance=5, string arg="");
    clock_t  getMeanElapsedTimeRepl    ( ) {return _meanElapsedTimeRepl;}
    
    TMetapop* get_popPtr(){return _thePop;}
    unsigned int get_nbThreads(){return _nbThreads;}

};

#endif /* defined(treplicateH) */
