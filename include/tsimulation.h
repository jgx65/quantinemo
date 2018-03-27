/** @file tsimulation.h
 *
 *   Copyright (C) 2008 Samuel Neuenschwander <samuel.neuenschwander@unil.ch>
 *
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

#ifndef tsimulationH
#define tsimulationH

#include "basicsimulation.h"
using namespace std;


/** TSimulation contains the stats db and launches each replicate separately */


class TSimulation: public SimBuilder{
public:
    map<string, string> _params;
    
    string _postexec_script;
    string _preexec_script;
    
    string _iniFile;
    
    unsigned int _threads;  // number of threads for this simulation
    
    /** should existing files be overwriten without asking? */
    bool      _overwriteFiles;
    int       _logfile_type;       // 0: as input, 1: minimal, 2: maximal
    string    _logfile;
    
    string _exe_directory;      // directory of the exe
    string _iniFile_directory;  // folder of the ini file
    string _working_directory;  // folder of any output, or relative to this
    string _basename;           // parameter "filename"
    string _simfolder;          // folder name contining all the simulation outputs
    string _simfolder_short;
    
    unsigned int _current_sim;      // the current simulation
    unsigned int _nbSims;           // the total number of simulations
    
    RAND* randEngines;          // an array of random generators for each thread
    bool  _fixedSeed;           // was the random generator seed fixed by the user?
    vector<string> _vSeeds;     // the seed for each random number engine
    bool  _randomPerReplicate;  // true: random numbers in input file are drawn for each raplicate separately
                                // false: random numbers in input file are identical for all replicates
    
    /** time parameters */
    char      _startTimeSimStr[20], _endTimeSimStr[20];
    time_t    _startTimeSim, _endTimeSim;
    string    _elapsedTimeSim;
    double*   _elapsedTimeRepl; // stores for each replciate the computation time
    
    unsigned int* _genLength;   // stores for each replciate the number of simualted generations (before exinction)
    
    TReplicate* _testRepl;

public:
    TSimulation();
    ~TSimulation();
    
    TStat_db* stats;        // only available if stats are computed
    
    list<ParamSet*>* _allParams;
    
    unsigned int _replicates;           // Number of replicates to iterate.
    unsigned int _generations;          // Number of generations to run (used for db).
    unsigned int _nb_stats;             // Number of stats (used for stat db)
    
    
    // random generator
    void init_randGenerator();
    bool get_fixedSeed() {return _fixedSeed;}

    
    void set_elapsedTimeRepl(clock_t time, unsigned int pos){_elapsedTimeRepl[pos]=time;}
    void set_genLength(unsigned int gen, unsigned int pos){_genLength[pos]=gen;}
    
    void run_sim(map<string, string> params, map<string, string> keys,
                 unsigned int i, unsigned int nbSims, string working_directory,
                 unsigned int nbThreads);
    void print_help(ostream& os, unsigned int wide1=20, unsigned int wide2=60,
                    char fill='.', unsigned int importance=5, string arg="");
    unsigned int get_nbTraits(string name, map<string, string>& inputParams);
    unsigned int isCoalescence(map<string, string>& inputParams);
    bool init_fileServices();
    void init_paramset();
    
    // log file
    void printLog();
    void printLogHeader();
    void printLogEnd();
    
    unsigned int get_replicates() {return _replicates;}
    
    void init_stat_db(map<string, string> params);
    
    clock_t print_start_replicate(unsigned int currentReplicate);
    

    bool getOverwriteFiles()             {return _overwriteFiles;}
    void setOverwriteFiles(bool b)       {_overwriteFiles=b;}
    void set_logfile_type(int val)       {_logfile_type=val;}
    void set_logfile(string s)           {_logfile=s;}

    void print_log_params(list< ParamSet* >&  params);
    void save_simparams(list< ParamSet* >&  params, ostream& FH);
    string get_fileName(const string& name);
    
    string get_exe_directory()      {return _exe_directory;}
    string get_iniFile_directory()  {return _iniFile_directory;}
    string get_working_directory()  {return _working_directory;}
    string get_basename()           {return _basename;}
    string get_simfolder()          {return _simfolder;}
    string get_simfolder_short()    {return _simfolder_short;}

    void set_exe_directory(string s);
    void set_iniFile_directory(string s);
    void set_working_directory(string s);
    void set_basename(string name) {_basename = name;}
    void set_simfolder(string name);
    
    unsigned int get_threads(){return _threads;}
    
    unsigned int get_current_sim() {return _current_sim;}
    
    bool ini_stats(map<string, string> params, map<string, string> keys,
                   unsigned int i, unsigned int nbSims,
                   string exe_directory, unsigned int nbThreads);

};
void run_replicates(TSimulation* ptr, map<string, string> params,
                    map<string, string> keys, unsigned int from,
                    unsigned int to, unsigned int curThread, unsigned int nbThreads);

#endif /* defined(__quantiNemo__tsimulation__) */
