/** @file lce_misc.h
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

#ifndef lce_miscH
#define lce_miscH


#include "lifecycleevent.h"
#include "filehandler.h"
#include "tstat_db.h"

class TStat_db;

// LCE_FileServicesNotifier
//
/**Event used to notify all file handlers to update their state through the FileServices::notify() interface.*/
class LCE_FileServicesNotifier: public LCE{
    
private:
    FileServices* _service;
    
    
public:
    LCE_FileServicesNotifier(int rank = my_NAN)
    : LCE("save_files", "output of files", "", rank), _service(0){
    }
    
    virtual ~LCE_FileServicesNotifier( ) { }
    
    virtual void execute ();
    
    virtual LCE_FileServicesNotifier* clone ( ) {return new LCE_FileServicesNotifier();}
    
    //SimComponent overrides:
    virtual void  loadFileServices ( FileServices* loader ) {
        _service = loader;
    }
    virtual void  loadStatServices (StatServices* loader ) {}
    virtual age_t removeAgeClass ( ) {return 0;}
    virtual age_t addAgeClass ( ) {return 0;}
    virtual age_t requiredAgeClass () {return 0;}
    
};

// LCE_StatFH
//
/**FileHandler of the LCE_StatServiceNotifier class, writes the recorded stats to txt files.*/
class LCE_StatFH : public FileHandler {
public:
    StatServices* _statService;
    
public:
    
    LCE_StatFH () : _statService(0) {set_extension(".txt");}
    
    ~LCE_StatFH() { };
    
    virtual string getName() {return "LCE_SataFH";}
    
    void set_statService(StatServices* srv) {
        _statService = srv;
        
        TStat_db* pStats = _popPtr->get_pStat_db();
        
        // set the filenames
        string base = get_path() + _popPtr->getBaseFileName();
        pStats->set_file_name_legend(base + "_legend.txt");
        pStats->set_file_name_stats(base + "_stats.txt");
        pStats->set_file_name_mean(base + "_mean.txt");
        pStats->set_file_name_median(base + "_median.txt");
        pStats->set_file_name_var(base + "_var.txt");
    }
    
    virtual void FHwrite(){}
    
    void set_save_choice(const int& i){_popPtr->get_pStat_db()->set_save_choice(i);}
    
    
    };

// LCE_StatSH
//
/**StatHandler of the LCE_StatServiceNotifier class, adds a default StatRecorder to the recorders list (alive.rpl).*/
class LCE_StatSH : public StatHandler<LCE_StatSH> {
    
public:
    LCE_StatSH () {}
    ~LCE_StatSH () {}
    
    virtual string getName() {return "LCE_StatSH";}
    virtual bool init ( ){
        StatHandler<LCE_StatSH>::init();
#ifdef _DEBUG
        message("  default (");
#endif
        
        if (_popPtr->get_pStat_db()->get_save_choice() != 1){
            add("Nb alive replicates","alive.rpl",GEN,ALL,0,&LCE_StatSH::get_isAlive,0,0,0);
        }
        
#ifdef _DEBUG
        message(")\n");
#endif
        return true;
    }
    
    virtual bool setStatRecorders(const string& token) {return false;}
    
    double get_isAlive ();
};

// LCE_StatServiceNotifier
//
/**Initiates the StatServices' parameters (log time) when registering, calls StatServices::notify() when executing.
 * Registers the file handler used to save the stats to the '.txt' and '_mean.txt' files.*/
class LCE_StatServiceNotifier: public LCE
{
private:
    StatServices* _service;
    
    unsigned int _save_choiceTemp;
    
    TStat_db* _pStat_db;
    
    string _stat_arg;       // list of stats to print in the stat file
    string _param_arg;      // list of parameters to print in the stat file
    
    LCE_StatFH _fileHandler;
    LCE_StatSH _statHandler;
    
public:
    
    LCE_StatServiceNotifier(int rank = my_NAN);
    
    virtual ~LCE_StatServiceNotifier ( ) { }
    
    virtual string getName() {return "LCE_StatServiceNotifier";}
    
    virtual void  execute ();
    
    virtual LCE_StatServiceNotifier* clone ( ) {return new LCE_StatServiceNotifier();}
    
    virtual bool init (Metapop* popPtr);
    
    //SimComponent overrides:
    virtual void loadFileServices ( FileServices* loader ) {loader->attach(&_fileHandler);}
    virtual void loadStatServices ( StatServices* loader );
    virtual age_t removeAgeClass ( ) {return 0;}
    virtual age_t addAgeClass ( ) {return 0;}
    virtual age_t requiredAgeClass () {return 0;}
    
    virtual void executeBeforeEachGeneration(const unsigned int& gen){}
    void temporal_change(const unsigned int& gen);
    
    vector<unsigned int> get_logtime_occurences(Param* pParam, unsigned int totGen);

    
};



// LCE_Aging
//
/**Removes all adults from the patches adult containers.
 * This is the only LCE that actually removes the adults from the patches. It is thus
 * mandatory to use it once in the life cycle, preferably after breeding!**/
class LCE_Aging: public LCE
{
public:
    
    LCE_Aging(int rank = my_NAN) : LCE("aging","aging","",rank) {}
    
    virtual ~LCE_Aging( ) { }
    
    virtual void execute ();
    
    virtual LCE_Aging* clone ( ) {return new LCE_Aging();}
    
    //implementations:
    virtual void loadFileServices ( FileServices* loader ) {}
    virtual void loadStatServices ( StatServices* loader ) {}
    virtual age_t removeAgeClass ( ) {return ADULTS;}
    virtual age_t addAgeClass ( ) {return 0;}
    virtual age_t requiredAgeClass () {return ADULTS | OFFSPRG;}
};


// LCE_store_popSizes
//
/** This life cycle allows to store the curent population sizes in the db for the coalecence **/
class LCE_store_popSizes: public LCE
{
public:
    LCE_store_popSizes(int rank = my_NAN) : LCE("store_pop_size","population size storage (LCE coalescence)","",rank) {}
    
    virtual ~LCE_store_popSizes( ) { }
    
    virtual bool init (Metapop* popPtr);
    
    virtual void execute ();
    
    virtual LCE_store_popSizes* clone ( ) {return new LCE_store_popSizes();}
    
    //SimComponent overrides:
    virtual void loadFileServices ( FileServices* loader ) {}
    virtual void loadStatServices ( StatServices* loader ) {}
    virtual age_t removeAgeClass ( ) {return 0;}
    virtual age_t addAgeClass ( ) {return 0;}
    virtual age_t requiredAgeClass () {return ADULTS | OFFSPRG;}
};



#endif //LCEMISC_H
