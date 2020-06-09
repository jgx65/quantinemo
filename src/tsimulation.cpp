/** @file tsimulation.cpp
 *
 *   Copyright (C) 2008 Samuel Neuenschwander <samuel.neuenschwander@unil.ch>
 *   Copyright (C) 2018 Frederic Michaud <frederic.a.michaud@gmail.com>

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

#include "tsimulation.h"
#include "tstat_db.h"
#include "filehandler.h"
#include "treplicate.h"
#include "version.h"
#include <errno.h>
#include <limits>
#include <algorithm>


using namespace std;


// static members
//----------------------------------------------------------------------------------------
// TSimulation
// ----------------------------------------------------------------------------------------
TSimulation::TSimulation(): _replicates(0), stats(0), _threads(0), _elapsedTimeRepl(0), _genLength(0), _testRepl(0), randEngines(0)
{
#ifdef _DEBUG
    message(" TSimulation::TSimulation()\n");
#endif
    init_paramset();
}

//----------------------------------------------------------------------------------------
// ~TSimulation
// ----------------------------------------------------------------------------------------
TSimulation::~TSimulation()
{
#ifdef _DEBUG
    message(" TSimulation::~TSimulation()...\n");
#endif
    if(stats) delete stats;
    delete[] _elapsedTimeRepl;
    delete[] _genLength;
    if(_testRepl) delete _testRepl;
    if(randEngines) delete[] randEngines;
#ifdef _DEBUG
    message(" TSimulation::~TSimulation() done.\n");
#endif
}

//----------------------------------------------------------------------------------------
// init_paramset
// ----------------------------------------------------------------------------------------
void
TSimulation::init_paramset()
{
    _paramSet.set_name("simulation");
    _paramSet.set_name_full("simulation");
    _paramSet.set_isRequired(true);
    _paramSet.set_pParamManager(this);
    
    _paramSet.add_param("_settings_file",STR,false,my_NAN,my_NAN,"quantiNemo.ini",false,
                        "Settings file name with path.",5);
    
    _paramSet.add_param("filename",STR,false,my_NAN,my_NAN,"simulation",false,
                        "Base filename for all outputs, if they are not specified individually.",1);
    
    _paramSet.add_param("folder",STR,false,my_NAN,my_NAN,"simulation_yyyy-mm-dd_hh-mm-ss",false,
                        "Output folder of the output. If not specified the " \
                        "folder will have a date and time stamp. \"\" may be " \
                        "used to store the output directly in the working directory.",0);
    
    _paramSet.add_param("working_directory",STR,false,my_NAN,my_NAN,"0",false,
                        "The working directory will be the:\n" \
                        "  0: current working directory\n" \
                        "  1: settings file directory\n" \
                        "  2: executable directory\n" \
                        "  string: explicitly defined directory",3);
    
    _paramSet.add_param("logfile",STR,false,my_NAN,my_NAN,"quantiNemo.log",false,
                        "Filename (including extension) of the simulation log. " \
                        " It is stored next to the executable.",5);
    
    _paramSet.add_param("overwrite",INT2,false,0,1,"0",false,
                        "How to tread present output:\n"
                        "  0: no (the user will be asked)" \
                        "  1: yes (present output is overwritten without asking)",1);
    
    _paramSet.add_param("preexec_script",STR,false,my_NAN,my_NAN,"",false,
                        "Name of script to execute before the simulation.",5);

    _paramSet.add_param("postexec_script",STR,false,my_NAN,my_NAN,"",false,
                        "Name of script to execute after the simulation.",5);
    
    _paramSet.add_param("seed",INT_MAT,false,my_NAN,my_NAN,"-1",false,
                        "Seed to initialize the random generator. If not " \
                        "specified the time is used.",2);
    
    _paramSet.add_param("logfile_type",INT2,false,0,2,"0",false,
                        "How extensive should the logfile be:\n" \
                        "  0: as input (log contains the same parameters as " \
                        "the settings file, plus the parameter seed)\n" \
                        "  1: minimal (log contains a minimal set of parameters: " \
                        "parameters with default values are removed)\n" \
                        "  2: maximal (log contains all available parameters with details)",2);
    
    _paramSet.add_param("all_combinations",INT2,false,0,1,"1",false,
                        "How to treat multiple sequential parameters:\n" \
                        "  0: order (simulations with the different orders will " \
                        "be simulated)\n" \
                        "  1: full (all possible combinations will be simulated)",5);
    
    _paramSet.add_param("replicates",INT2,false,0,my_NAN,"1",false,
                        "Number of replicates to perform per simulation.",0);
    
    _paramSet.add_param("stat_NaN", STR, false, my_NAN, my_NAN, my_NANstr,false,
                        "String to use for incalculable statistics.",5);
    
    
    _paramSet.add_param("random_per_replicate", INT2, false, 0, 1, "0",false,
                        "How to treat random macros with multiple replicates:\n" \
                        "  0: same (all replicates share the same random number)\n" \
                        "  1: individual (all replicates have an individual random number)"
                        ,5);
    
    _paramSet.add_param("help", INT2, false, my_NAN, my_NAN, "",false,
                        "Show this nice help (--help or -h).",5);
    
    _paramSet.add_param("version", INT2, false, my_NAN, my_NAN, "",false,
                        "Show the program version (--version or -v).",5);
    
}

//----------------------------------------------------------------------------------------
// help
// ----------------------------------------------------------------------------------------
/** run each simulation indvidually
 * loop here over the replciates
 * the stats have to be stored here
 * metapop has to be created/deleted within the replicate loop!
 */
void TSimulation::print_help(ostream& os, unsigned int wide1, unsigned int wide2,
                             char fill, unsigned int importance, string arg)
{
    TReplicate rep(this);
    rep.print_help(os, wide1, wide2, fill, importance, arg);
 }

//----------------------------------------------------------------------------------------
// run_sim
// ----------------------------------------------------------------------------------------
/** run each simulation indvidually
 * loop here over the replciates
 * the stats have to be stored here
 * metapop has to be created/deleted within the replicate loop!
 */
void TSimulation::run_sim(map<string, string> params,
                          map<string, string> keys,
                          unsigned int i, unsigned int nbSims,
                          string exe_directory, unsigned int nbThreads)
{
    _current_sim = i;
    _nbSims = nbSims;
    _current_sim = i;
    _exe_directory=exe_directory;
	_threads=nbThreads;
    
    if(!_preexec_script.empty()) runScript(get_fileName(_preexec_script), get_fileName(_iniFile));
    
    //build the list of parameters from the record:
    // get all parameters
    build_allParams();
    if (!build_currentParams(params, keys, NULL, true)) error("SimRunner::run:couldn't build current params\n");
    init_randGenerator();
    
    time (&_startTimeSim);
    strftime(_startTimeSimStr, 20, "%d-%m-%Y %H:%M:%S", localtime (&_startTimeSim));

    // test all parameters and initialize the stat_db
    if(!ini_stats(params, keys, i, _nbSims, exe_directory, nbThreads)) return;
    
    // supress the same warning messages of the replicates
    verbose_warning=0;  
        
    if (_nbSims > 1) message("\n\n----- SIMULATION %i/%i -----\n", _current_sim+1, _nbSims);
    
    // set_seed(params);   // used here since some macros need the random generator
    
    // get the replicate number
    assert(!_elapsedTimeRepl);
    assert(!_genLength);
    ARRAY::create_1D(_elapsedTimeRepl, _replicates,(double) my_NAN);
    ARRAY::create_1D(_genLength, _replicates,(unsigned int) my_NAN);
    

    
    //-------------------------------------------------------------------------------------
    message("\n\nSIMULATION");
    
    
    
    /////////////////////////////////////////////////////////////////////////////
    // for each replicate
    run_replicates(this, params, keys, 0, _replicates, 0, _threads);

    
    /////////////////////////////////////////////////////////////////////////////
    
    if(stats) stats->FHwrite();
    
    
    //---------------------------------------------------------------------------------------
    if(!_postexec_script.empty()) runScript(get_fileName(_postexec_script), get_fileName(_iniFile));

    time (&_endTimeSim);
    strftime(_endTimeSimStr, 20, "%d-%m-%Y %H:%M:%S", localtime (&_endTimeSim));

    _elapsedTimeSim = getElapsedTime(_endTimeSim, _startTimeSim);
    
    
    if (nbSims == 1) message("\n\n----- SIMULATION done (CPU time: %ss) -----\n", _elapsedTimeSim.c_str());
    else message("\n\n----- SIMULATION %i/%i done (CPU time: %ss) -----\n", _current_sim+1, _nbSims, _elapsedTimeSim.c_str());
    
    printLog();
    printLogEnd(); // write the duration of the simulation to the log file
}

// ----------------------------------------------------------------------------------------
// loadDefaultTemplates
// ----------------------------------------------------------------------------------------
/** run the number of replicates per thread */
void
run_replicates(TSimulation* ptr, map<string, string> params, map<string, string> keys,
               unsigned int from, unsigned int to, unsigned int curThread,
               unsigned int nbThreads)
{
#ifdef _DEBUG
    unsigned int fromConst=from;
    message("\nTSimulation::run_replicates (multithreaded): replicates %i to %i ...\n", from+1, to);
#endif
    assert(to<=ptr->_replicates);
    for(; from<to; ++from) {
        TReplicate(ptr, params, keys, from, curThread, nbThreads);
    }
#ifdef _DEBUG
    message("\nTSimulation::run_replicates (multithreaded): replicates %i to %i ... done.\n", fromConst+1, to);
#endif
}

// ----------------------------------------------------------------------------------------
// SimRunner::printLogHeader()
// ----------------------------------------------------------------------------------------
void TSimulation::printLogHeader()
{
    // if file is present, just skip the step
    ifstream IF(_logfile.c_str());
    if (IF) {
        IF.close();
        return;
    }
    
    ofstream FH(_logfile.c_str());
    if (!FH) error("could not create simulation logfile \"%s\"\n", _logfile.c_str());
    
    FH << "--- q u a n t i N E M O 2 ---\n" << "    LOGFILE\n\n\n";
    FH
    << "| basename                                | simfolder                               |      start time     |      stop time      | duration   |"
    << " repl done | mean rpl   | mean gen | version                 | hostname             | output files \n";
    
    FH.close();
}

//--------------------------------------------------------------------------------------
// SimRunner::printLog()
//-------------------------------------------------------------------------------------
void TSimulation::printLog()
{
    ofstream FH(_logfile.c_str(), ios::app);
    
    if (!FH.is_open()) {
        error("could not open simulation logfile \"%s\"\n", _logfile.c_str());
        return;
    }
    
    string simfolder = get_basename();
    
    FH << "| ";
    FH.width(40);
    FH.setf(ios::left, ios::adjustfield);
    if (simfolder.empty()) FH << "-";
    else FH << simfolder;
    
    FH << "| ";
    FH.width(40);
    FH.setf(ios::left, ios::adjustfield);
    FH << get_simfolder();
    
    FH << "| " << _startTimeSimStr << " | " << _endTimeSimStr << " | ";
    
    FH.width(10);
    FH.setf(ios::right, ios::adjustfield);
    FH << _elapsedTimeSim << " | ";
    
    FH.width(9);
    FH << _replicates << " | ";
    
    FH.width(10);
    FH << getElapsedTime(ARRAY::mean(_elapsedTimeRepl, _replicates)) << " | ";
    
    FH.width(8);
    FH << ARRAY::mean(_genLength, _replicates) << " | ";
    
    FH << "[" << VERSION_DATE << "; " << VERSION_TIME << "]";
    
    FH << " | ";
    FH.width(20);
    FH.setf(ios::left, ios::adjustfield);
    char* host;
    if ((host = getenv("HOST")) != NULL) FH << host << " |";
    else if ((host = getenv("HOSTNAME")) != NULL) FH << host << " |";
    else FH << "-" << " |";
    
    list< FileHandler* >::iterator cur, end = _testRepl->_FileServices->get_children().end();
    for (cur = _testRepl->_FileServices->get_children().begin(); cur!=end; ++cur){
        string rrr =(*cur)->get_extension();
        string ggg = (*cur)->get_short_path();
        FH << " \"" << (*cur)->get_extension() << "\":" << (*cur)->get_short_path();
    }
    FH << "\n";
    
    FH.close();
}

// ----------------------------------------------------------------------------------------
/** get the number of traits of the given trait type */
unsigned int
TSimulation::get_nbTraits(string name, map<string, string>& inputParams)
{
    // if the number of loci is not specified no trait will be inizialized
    map<string, string>::iterator pos = inputParams.find(name + "_loci"); // general param
    if (pos == inputParams.end()) pos = inputParams.find(name + "_loci_1"); // or at least for the frist trait
    if (pos == inputParams.end()) return 0; // the number of loci is not specified
    if (pos == inputParams.end()) return 0; // the number of loci is not specified
    try {
        if (!strTo<unsigned int>(pos->second)) return 0;
    }
    catch(...) {
        error("The parameter '%s_loci' (%s) must be an integer!\n",
              name.c_str(), pos->second.c_str());
    }
    
    // get the number of traits
    pos = inputParams.find(name + "_nb_trait");
    if (pos == inputParams.end()) return 1; // the number of traits is not specified
    unsigned int nb;
    try {
        nb = strTo<unsigned int>(pos->second);
    }
    catch(...) {
        error(
              "The parameter '%s_nb_trait' must have a positive number as argument!\n",
              name.c_str());
    }
    return nb;                          // return the specified number of traits
}


// ----------------------------------------------------------------------------------------
/** is it a normal or a coalescent simulation */
unsigned int
TSimulation::isCoalescence(map<string, string>& inputParams)
{
    map<string, string>::iterator pos = inputParams.find("coalescence"); // general param
    if (pos == inputParams.end()) return false; // the parameter is not specified
    try {
        return strTo<bool>(pos->second);
    }
    catch(...) {
		error("The parameter 'coalescence' requires an integer as argument (%s)!\n",
              pos->second.c_str());
    }
    return false;
}

//----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
bool
TSimulation::init_fileServices()
{
    string base_name =_paramSet.getArg("filename");
    if (_nbSims > 1) base_name += "-" + to_string(_current_sim);
    set_basename(base_name);
	string wd = _paramSet.getArg("working_directory");
    if (isNumber<unsigned int>(wd)) {
        switch (strTo<unsigned int>(wd)) {
            case 0: // working directory is the current working directory
                set_working_directory(get_cwd());
                break;
            case 1: // working directory is the directory containing the ini file
                set_working_directory(get_iniFile_directory());
                break;
            case 2: // working directory is the directory containing the executable
                set_working_directory(get_exe_directory());
                break;
        }
	}
    else set_working_directory(wd); // working directory is explicitly defined
    string name;
    // set the simulation folder name
    if (!this->_paramSet.isSet("folder")) {
        // generate the default folder name with current date
        char datetime[20];
        time_t t = time(NULL);
        strftime(datetime, 20, "%Y-%m-%d_%H-%M-%S", localtime(&t));
        name = STRING::get_file_basename(this->_paramSet.getArg("_settings_file")) + "_";
        //string name = "simulation_";
        name += datetime;
    }
    else name = this->_paramSet.getArg("folder"); // folder name is set
    set_simfolder(name);
    set_logfile(get_working_directory() + _paramSet.getArg("logfile"));
    set_logfile_type((int) _paramSet.getValue("logfile_type"));
    setOverwriteFiles((bool) this->_paramSet.getValue("overwrite"));
    
    printLogHeader();
    
    _postexec_script = _paramSet.getArg("postexec_script");
    _preexec_script = _paramSet.getArg("preexec_script");
    
    return true;
}

// ----------------------------------------------------------------------------------------
// run_replicate_thread
// ----------------------------------------------------------------------------------------
clock_t
TSimulation::print_start_replicate(unsigned int currentReplicate)
{
#ifdef _DEBUG
#ifdef _SHOW_MEMORY
    message("\n\r**** replicate %i/%i [%s] ****   RAM: %f MB\n",currentReplicate+1,_replicates
            ,getElapsedTime(0).c_str(), process_mem_usage());
#else
    message("\n\r**** replicate %i/%i [%s] ****\n",currentReplicate+1,_replicates
            ,getElapsedTime(0).c_str());
#endif
#else
#ifdef _SHOW_MEMORY
    message("\n\r    replicate %i/%i [%s]   RAM: %f MB",currentReplicate+1,_replicates
            ,getElapsedTime(0).c_str(), process_mem_usage());
#else
    message("\n\r    replicate %i/%i [%s]",currentReplicate+1,_replicates
            ,getElapsedTime(0).c_str());
#endif
#endif
    fflush(stdout);
    return clock();
}

//----------------------------------------------------------------------------------------
// print_info
// ----------------------------------------------------------------------------------------
void
TSimulation::set_exe_directory(string s)
{
    if(s.empty())                       _exe_directory = ""; 	// nothing passed: to prevent problems
    else if(s[0] != '/' && s[1] != ':') _exe_directory = "";	// no absolute path passed, just the name
    else{                                                     	// absolute path passed: remove the file name from it
        assert(s.find_last_of("/\\") != string::npos);
        _exe_directory = s.substr(0, s.find_last_of("/\\")+1);
    }
    
#ifdef _DEBUG
    message("Exe directory: '%s'\n", _exe_directory.c_str());
#endif
}

//----------------------------------------------------------------------------------------
// print_info
// ----------------------------------------------------------------------------------------
void
TSimulation::set_iniFile_directory(string s)
{
    if(!s.empty()){
        if(s[s.length()-1]=='/' || s[s.length()-1]=='\\') _iniFile_directory = s;
        else _iniFile_directory = s + SEP;
    }
#ifdef _DEBUG
    message("\nIniFile directory: '%s'", _iniFile_directory.c_str());
#endif
}

//----------------------------------------------------------------------------------------
// print_info
// ----------------------------------------------------------------------------------------
void
TSimulation::set_working_directory(string s)
{
    if(!s.empty()){
        if(s[s.length()-1]=='/' || s[s.length()-1]=='\\') _working_directory = s;
        else _working_directory = s + SEP;
    }
#ifdef _DEBUG
    message("\nWorking directory: '%s'\n", _working_directory.c_str());
#endif
}

// ----------------------------------------------------------------------------------------
// save_simparams
// ----------------------------------------------------------------------------------------
/**Sets the folder name of the simulation*/
void
TSimulation::set_simfolder (string name)
{
    _simfolder_short = name;
    _simfolder = get_working_directory() + name;
    if(!_simfolder.empty() && _simfolder[_simfolder.length()-1] != SEP) _simfolder += SEP;
    check_path(_simfolder);
}

// ----------------------------------------------------------------------------------------
// save_simparams
// ----------------------------------------------------------------------------------------
void
TSimulation::printLogEnd()
{
    string logfile = _simfolder + _basename + ".log";
    ofstream FH(logfile.c_str(),ios::app);
    if(!FH) error("FileServices::init: could not output sim parameters to\"%s\"\n",logfile.c_str());
    else {
        FH<<"\n#############################################"
        <<  "\n#  simulation ended on " << _endTimeSimStr << "  #"
        <<  "\n#                                           #"
        <<  "\n#  duration of the simulation was " << _elapsedTimeSim << "  #"
        <<  "\n#############################################\n";
        
        FH.close();
    }
}

// ----------------------------------------------------------------------------------------
// save_simparams
// ----------------------------------------------------------------------------------------
void
TSimulation::print_log_params(list< ParamSet* >&  params)
{
    // write log file
    string logfile = _simfolder + _basename + ".log";
    check_path(_simfolder);
    ofstream FH(logfile.c_str(),ios::out);
    if(!FH){
        error("FileServices::init: could not open log file '%s'!\n",logfile.c_str());
    }
    else {
        save_simparams(params, FH);
        FH.close();
    }
}



// ----------------------------------------------------------------------------------------
// save_simparams
// ----------------------------------------------------------------------------------------
void
TSimulation::save_simparams(list< ParamSet* >&  params, ostream& FH)
{
    FH << "#############################################"
    <<  "\n#                quantiNemo 2               #"
    <<  "\n#   version "<<RELEASE<<"."<<REVISION<<"."<<MINOR_VERSION<<TEMP_VERSION
    <<" ["<<VERSION_DATE << "; " << VERSION_TIME <<"]   #"
#ifdef VERSIONGIT
        <<  "\n#       Commit N° : " << VERSIONGIT <<  "           #"
#endif
    << "\n#                                           #"
    << "\n#              Simulation log               #"
    << "\n# simulation started on " << _startTimeSimStr << " #\n";
    
    switch(_logfile_type){
            case 0:  // as input file
            FH << "#              (log as input)               #"
            <<  "\n#############################################\n\n";
            for(list< ParamSet* >::iterator Pit = params.begin(); Pit != params.end(); ++Pit) {
                (*Pit)->print(FH, 1);
            }
            break;
            case 1: // minimal
            FH << "#               (minimal log)               #"
            <<  "\n#############################################\n\n";
            for(list< ParamSet* >::iterator Pit = params.begin(); Pit != params.end(); ++Pit) {
                (*Pit)->print_minimal(FH, 1);
            }
            break;
            case 2: // maximal
            FH << "#               (maximal log)               #"
            <<  "\n#############################################\n\n";
            for(list< ParamSet* >::iterator Pit = params.begin(); Pit != params.end(); ++Pit) {
                FH << "\n";
                (*Pit)->print_maximal(FH, 1);
            }
            break;
    }
    FH<<endl;
}

//----------------------------------------------------------------------------------------
// init_randGenerator
// ----------------------------------------------------------------------------------------
/** Initialize a random generator for each thread (set seed)
 */
void
TSimulation::init_randGenerator()
{
    vector<unsigned int> seed;
    Param* pParam = _paramSet.get_param("seed");
    
    // check if the seed is set
    string text = pParam->get_arg();
    if (text != "-1") {               // if default seed (-1) was given stop
        // seed is given by input
        _fixedSeed = true;
        if (text[0] == '{') {       		// seed is given as an array
            seed = STRING::strMatrix2vector<unsigned int>(text);
        }
        else {                     	// seed is given as a single value
            seed.push_back(strTo<unsigned int>(text));
        }
#ifdef _DEBUG
        message("Random generator seed taken from the setting file: %s\n", text.c_str());
#endif
    }
    else {             				// if initialized by time
        _fixedSeed = false;
        seed.push_back(time(NULL) % std::numeric_limits<unsigned int>::max()); // seconds since January 1, 1970.
        seed.push_back(clock() % std::numeric_limits<unsigned int>::max()); // number of clock ticks elapsed since the program started.
        pParam->set_isSet(true);
#ifdef _DEBUG
        message("Random generator is initialized with a time seed: %i %i\n", seed[0], seed[1]);
#endif
    }
    
    // create the random engines for each thread
    assert(!randEngines);
    
    randEngines = new RAND[_threads];
    vector<unsigned int>::iterator cur, end;
    string curSeed;
    for(unsigned int i=0; i<_threads; ++i){
        switch(i){
            case 0: break;                          // use the given (original) seed
            case 1: seed.push_back(i+1); break;     // add a supplement seed: the thread
            default: seed.back() = i+1; break;      // replace the supplement seed by the index
        }
        
        curSeed="{";
        for(cur=seed.begin(), end=seed.end(); cur!=end; ++cur){
            if(cur!=seed.begin()) curSeed += " ";
            curSeed += toStr<unsigned int>(*cur);
        }
        curSeed += "}";
       // pParam->set_arg(curSeed);           // set the seed used for this replicate
        _vSeeds.push_back(curSeed);
        randEngines[i].set_seed(seed);      // initialize the random Engine
    }
}


// ----------------------------------------------------------------------------------------
// get_fileName
// ----------------------------------------------------------------------------------------
/** returns the filename as it is if the path is absolute or adds the program directory path
 * if relative
 */
string TSimulation::get_fileName(const string& name)
{
    if (name[0] == '/' || name[1] == ':') return name;    // absolute path
    return get_iniFile_directory() + name;    // relative path
}

// ----------------------------------------------------------------------------------------
// ini_stats
// ----------------------------------------------------------------------------------------
/* runs quantiNemo as test to check if all varaibles are well set and initializes the stat_db */
bool
TSimulation::ini_stats(map<string, string> params, map<string, string> keys,
                       unsigned int i, unsigned int nbSims, string exe_directory,
                       unsigned int nbThreads)
{
#ifdef _DEBUG
    message("\n\n**************************************************");
    message("\n***** Check parameters and set up stats ...");
#endif
    _nbSims = nbSims;
    _current_sim = i;
    _exe_directory=exe_directory;
	_threads=nbThreads;
    
    set_iniFile_directory(STRING::get_file_path(params["_settings_file"]));
    
    init_fileServices();
    
    // get the replicate number
    _replicates = (unsigned int)_paramSet.getValue("replicates");
    _randomPerReplicate = (bool)_paramSet.getValue("random_per_replicate");
    

    stats = new TStat_db(params, _replicates);
    
    // generate metapop in order to asses all stats
    _testRepl = new TReplicate(this);
    _allParams = &(_testRepl->get_allParams());
    _testRepl->rand = randEngines;
    _testRepl->test_replicate_and_setUpStats(params, keys);
    
    //init the file services:
    //->files check, false means the user wants to skip this simulation.
    //->save simparameters in log files
    if(!_testRepl->_FileServices->check_file_overwrite()) return false;

    list< ParamSet* >& _ParamSet(_testRepl->_FileServices->get_params());
    for(list< ParamSet* >::iterator Pit = _ParamSet.begin(); Pit != _ParamSet.end(); ++Pit) {
        if((*Pit)->_name == "simulation") (*Pit)->set_param("seed",_vSeeds[0], _testRepl);

    }
    
    print_log_params( _ParamSet);
    
    
#ifdef _DEBUG
    message("\n***** Check parameters and set up stats ... end");
    message("\n**************************************************\n");
#endif
    return true;
}


