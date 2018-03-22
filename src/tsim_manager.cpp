/** @file tsim_manager.cpp
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


#include "tsim_manager.h"
#include "version.h"
#include "tsimulation.h"
#include "functions.h"



//#include <iostream>
using namespace std;


// ----------------------------------------------------------------------------------------
// TSimManager
// ----------------------------------------------------------------------------------------
TSimManager::TSimManager(int ARGC, char **ARGV)
{
    _exe_directory = ARGV[0];
    readInputFiles(ARGC, ARGV, _simSettings, _keyParams);
    _nbSims = (unsigned int)_simSettings.size();
    _nbThreadsTot = get_nbThreads();
    
    
    
#ifdef _DEBUG
    // ouptut all simulation settings
    vector<map<string, string> >::iterator curSim, endSim=_simSettings.end();
    vector<map<string, string> >::iterator curKey, endKey=_keyParams.end();
    map<string, string>::iterator curMap, endMap;
    unsigned int i;
    for (i=0, curSim=_simSettings.begin(), curKey=_keyParams.begin();
         curSim!=endSim; ++curSim, ++curKey, ++i){
        assert(curKey!=endKey);
        message("\n\n********** Simulation %i/%i **********", i+1, _nbSims);
        for(curMap=curKey->begin(), endMap=curKey->end(); curMap!=endMap; ++curMap){
            cout << "\n  set " << curMap->first << " [" << curMap->second << "]"<< flush;
        }
        for(curMap=curSim->begin(), endMap=curSim->end(); curMap!=endMap; ++curMap){
            cout << "\n  " << curMap->first << " [" << curMap->second << "]"<< flush;
        }
    }
#endif
}

//----------------------------------------------------------------------------------------
// run the simulations
// ----------------------------------------------------------------------------------------
/* run simulations sequentially or parallel when multithreaded
 _nbThreadsTot: total number of avaialble threads (specified)
 _nbSims:       total number of differnet simulations (specified)
 nbThread:      number of threads to use
 step:          number of simulations per thread
 restSims:      number of threads getting one more simulation
 nbThreatPerSim: number of threads per simulation
 nbThreatPerSimRest: number of sims getting one more thread
 */
void
TSimManager::print_help(ostream& os, unsigned int wide1, unsigned int wide2,
                                   char fill, unsigned int importance, string arg)
{
#ifdef _DEBUG
    message("\n\n**************************************************");
    message("\n***** Help ...");
#endif

    TSimulation theSim;
    theSim.print_help(os, wide1, wide2, fill, importance, arg);

#ifdef _DEBUG
    message("\n***** Help ... done");
#endif
}

//----------------------------------------------------------------------------------------
// run the simulations
// ----------------------------------------------------------------------------------------
/* run simulations sequentially or parallel when multithreaded
 _nbThreadsTot: total number of avaialble threads (specified)
 _nbSims:       total number of differnet simulations (specified)
 nbThread:      number of threads to use
 step:          number of simulations per thread
 restSims:      number of threads getting one more simulation
 nbThreatPerSim: number of threads per simulation
 nbThreatPerSimRest: number of sims getting one more thread
 */
void
TSimManager::run()
{
#ifdef _DEBUG
    message("\n\n**************************************************");
    message("\n***** Run simulations ...");
#endif
    
    run_sims_withinThread(this, 0, _nbSims, _nbThreadsTot);

    
#ifdef _DEBUG
    message("\n***** Run simulations ... done");
#endif
}

//----------------------------------------------------------------------------------------
// run_sim
// ----------------------------------------------------------------------------------------
/* run the specified simulations sequentially */
void
run_sims_withinThread(TSimManager* ptr, unsigned int from, unsigned int to, unsigned int nbThread)
{
    for (; from<to; ++from){
        TSimulation theSim;
        theSim.run_sim(ptr->_simSettings[from], ptr->_keyParams[from], from, ptr->_nbSims, ptr->_exe_directory, nbThread);
    }
}

//----------------------------------------------------------------------------------------
// get_nbThreads
// ----------------------------------------------------------------------------------------
/** multithreading is possible a the simulation and replicate level. Get the
 number of available threads. If the threads are differently defined in the 
 differentsettings files, output a warning and take the highest number 
 */
unsigned int
TSimManager::get_nbThreads()
{
    return 1;
}


//----------------------------------------------------------------------------------------
// read command line arguments and settings file
// ----------------------------------------------------------------------------------------
/** the param ARGV[0] was changed such as to be the folder name of the executable (working directory) */
void
TSimManager::readInputFiles(int ARGC, char **ARGV,
                            vector<map<string, string> >& paramRecord,
                            vector<map<string, string> >& keysRecord)
{
#ifdef _DEBUG
    message("\n\n**************************************************");
    message("\n***** Read parameters and settings files ...");
#endif
    map<string, string> args;   // all params passed as arguments
    vector<string> inputfiles;  // all file names (params not declared as ("--")
    
    // read arguments passed to the program
    if (ARGC == 1) {  // if no argument is passed use the default file name
        // is there a default settings file?
        if (STRING::file_exists("quantinemo2.ini")) inputfiles.push_back("quantinemo2.ini");
        else { // if not present ask for a name
            string name;
            message("\nPlease enter the settings file name (or 'help'): ");
            cin >> name;
            inputfiles.push_back(name);
        }
    }
    else {  // get all arguments
        // get all parameters passed as arguments (e.g. --generations 5000)
        // spaces within {}, (), "" are allowed, however not other spaces as for batch mode would be necessary
        
        // concatenate all words to a string with spaces
        string INstr;
        for (int i = 1; i < ARGC; ++i) {
            INstr += (string) ARGV[i] + " ";
        }
        
        // get the parameters and their arguments or the setting file names
        unsigned int curLine = 0;
        bool otherParams;
        char c;
        istringstream IN(INstr);
        while (!IN.eof() && IN.good()) {
            IN >> ws;                                  // remove foregoing space
            if (IN.peek() == '-') {                        // is it a parameter?
                // get parameter name
                IN.get();                                    // remove first '-'
                if (IN.peek() == '-') {               // parameter "--" (string)
                    IN.get();                               // remove second '-'
                    INstr = STRING::find_first_of(IN, "({", true); // get parameter name
                }
                else {                          // parameter "-" (single char)
                    while (IN.get(c)) { // if char parameters are aligned (-abc)
                        INstr = c;
                        if (isspace(IN.peek())) break;
                        args[INstr] = string(""); // no argument (if concatenated)
                    }
                }
                
                // get potential argument
                IN >> ws;
                if (IN.peek() == '-') args[INstr] = string(""); // no argument
                else if (INstr == "settings") inputfiles.push_back(INstr); // settings file name
                else args[INstr] = readArguments(IN, curLine, '\0', otherParams); // get the argument
                IN >> ws;
            }
            else {                                        // it is a file name
                IN >> INstr >> ws;
                inputfiles.push_back(INstr);
            }
        }
    }
    
#ifdef _DEBUG
    // output all passed arguments
    map<string, string>::iterator cur1, end1=args.end();
    message("\nParameters passed as argument:\n");
    for(cur1=args.begin();cur1!=end1; ++cur1) {
        message("   %s: %s\n", cur1->first.c_str(), cur1->second.c_str());
    }
#endif
    
    // check parameters for special commands (h, v,...)
    program_args(args, inputfiles);
    
    // read the settings files (first all parameters passed as arguments have to be read)
    vector<string>::iterator cur = inputfiles.begin(), end=inputfiles.end();
    for (; cur != end; ++cur) {
        map<string, vector<string> > params;
        map<string, vector<string> > keywords;
        
        // get all inputs (file and console)
        if (!read_settings_file(*cur, params, keywords)){         // get parameters from settings file
            error("Could not read the settings file '%s'!\nCurrent working directory: %s",
                  cur->c_str(), get_cwd().c_str());
        }
        add_programArgs(args, params);  // add the parameters passed as argument to the program
        
        // generate the simulation records
        map<string, vector<string> >::iterator posAll=params.find("all_combinations");
        if(posAll==params.end() || posAll->second.front()=="1"){
            build_records_full(paramRecord, keysRecord, params, keywords);     // build the different sim records (all combinations)
        }
        else build_records_reduced(paramRecord, keysRecord, params, keywords); // build the different sim records (reduced combiantions)
    }
    
#ifdef _DEBUG
    message("\n***** Read parameters and settings files ... done");
#endif
}


//----------------------------------------------------------------------------------------
//  printVersion
// ----------------------------------------------------------------------------------------
/** print the current version of quantiNemo and usefull information about the author/the commit */
void TSimManager::printVersion()
{
       string date = VERSION_DATE;

    message("Program:    quantiNemo 2 (quantitative genetics simulator)");
    message("\n\nVersion:    %i.%i.%i%s [%s; %s]", RELEASE,
            REVISION, MINOR_VERSION, TEMP_VERSION, VERSION_DATE,
            VERSION_TIME, date.substr(date.rfind(' ') + 1).c_str());
#ifdef VERSIONGIT
    message("\nCommit N°:  %s",VERSIONGIT);
#endif
    message("\n\nAuthors:    Samuel Neuenschwander (samuel.neuenschwander@unil.ch) &");
    message(  "\n            Frederic Michaud (frederic.michaud@unil.ch)");
    message(  "\n            Jerome Goudet (jerome.goudet@unil.ch)");
    message(  "\n            Department of Ecology and Evolution");
    message(  "\n            University of Lausanne, Switzerland\n");
    
}

//----------------------------------------------------------------------------------------
//  program_args
// ----------------------------------------------------------------------------------------
/** check the program args */
void
TSimManager::program_args(map<string, string>& args, vector<string>& inputfiles)
{
    //print the header
    printVersion();
    // If only the version is required, stop here
    map<string, string>::iterator cur, end = args.end();
    if (args.find("v") != end || args.find("version") != end) {
        throw 1111;
    }
    


    
    // is the help desired
    if (args.find("h") != end || args.find("help") != end || (!inputfiles.empty() && inputfiles.front()=="help")) {
        message(  "\nUsage:      quantiNemo2               (prompt for settings file)" \
                "\n            quantiNemo2 settings.ini  (parameters defined in settings file)" \
                "\n            quantiNemo2 args          (parameters passed as arguments)");
        
        message("\n\nManual:     An extensive manual can be downloaded from" \
                "\n            http://www.unil.ch/popgen/softwares/quantinemo/quantinemo_files/quantiNemo.pdf\n");
        
        // find argument if present
        string arg;
        cur=args.find("h");
        if(cur!= end) arg = cur->second;
        else{
            cur=args.find("help");
            if(cur!= end) arg = cur->second;
        }
        
        print_help(cout, 30, 60, '.', 5, arg);
        
        throw 1111;
    }
}

//------------------------------------------------------------------------------
// read_settings_file
// ----------------------------------------------------------------------------------------
/** the format of the parameters:
 * map<patch_capacity, vector({1000 1000}; rep(1000,2))>
 * map key: parameter name
 * 1. pos of vector: full argument (without any macro)
 * 2. pos of vector: original argument with macro        (optional)
 * return true if the file was read and false if there was a problem
 */
bool
TSimManager::read_settings_file(string name, map<string, vector<string> >& parsedParams,
                                map<string, vector<string> >& keywords)
{
    // open the settings file
    ifstream IN(name.c_str());
    if (!IN) return false;
    
    message("\nReading settings file '%s' ...", name.c_str());
    
    unsigned int linecnt = 0;
    string key;                     // string to store parameter name
    string args;                    // output string to collect parameter arguments
    string keyword;                 // the keyword specified as "set LAST 100"
    stringstream info;              // information on the success of the reading process
    vector<string> argvect;         // vector to store sequential arguments
    bool otherParams;               // true if there are more arguments ot be read for the same parameter
    parsedParams.clear();           // remove previous parameters if present
    string comment;                 // string for temp output or when error happens
    
    // add the settings file name
    set_iniFile_directory(STRING::get_file_path(name));
    add_parsedParam("_settings_file", name, parsedParams);
    
    //--------------------------------------------------------------------------------
    try {
        //read the file parameter by parameter, respectively line by line
        while (IN.good() && !IN.eof()) {
            argvect.clear();
            otherParams = true;
            keyword.clear();
            
            // remove the foregoing space of the line
            if (!STRING::removeCommentAndSpace(IN, linecnt, 0, '\n')) continue;
            
            //read the parameter name:
            key = STRING::find_first_of(IN, "{(", true);
            
            // remove space between parameter name and argument (if no argument jump the key)
            if (!STRING::removeCommentAndSpace(IN, linecnt, 0, '\n')) continue;
            
            comment = "\n  line " + toStr(linecnt) + ": " + key;
            info << comment;
            
#ifdef _DEBUG
            message("%s", comment.c_str());
            cout << flush;
#endif
            
            // get the arguments
            try {
                while (otherParams) {	// for each argument
                    args = readArguments(IN, linecnt, '\n', otherParams);
                    
                    // is the argument passed in an external file?
                    if (args[0] == '$') {
                        ifstream EXTERN(get_fileName(args.substr(1)).c_str());
                        //ifstream EXTERN(args.substr(1).c_str());
                        if (!EXTERN)
                            error("External file '%s' could not be found!\n",
                                  //FileHandler::get_fileName(args.substr(1)).c_str());
                                  args.substr(1).c_str());
                        unsigned int curLine = 0; // use a temp line counter for the external file
                        args = readArguments(EXTERN, curLine, '\0', otherParams);
                        EXTERN.close();
                    }
                    
                    if(key=="set" && keyword.empty()){ // a keyword is specified
                         keyword = args;
                        otherParams=true;
                    }
                    else if (!args.empty()) argvect.push_back(args);
                }
            }
            catch(const char* err) {
                string t = info.str();
                error(
                      "Could not read the settings file:%s\n***ERROR*** Parameter '%s' in line %i could not be read: %s\n",
                      t.c_str(), key.c_str(), linecnt, err);
            }
            
            if(key=="set") comment = " " + keyword + ":";
            else           comment = ":";
            unsigned int nb = (unsigned int)argvect.size();
            for (unsigned int i = 0; i < nb; ++i) {
                if (i) comment += " |";
                comment += " " + argvect[i];
            }
            comment += " (" + toStr(nb) + " args)";
            info << comment;
            
#ifdef _DEBUG
            message("%s", comment.c_str());
            cout << flush;
#endif
            
            if(key=="set") add_parsedKeys(keyword, argvect, keywords);
            else           add_parsedParam(key, argvect, parsedParams);
        }    //__END__WHILE__
    }
    catch(const char* err) {
        throw -1;     // that the message is not re-thrown
    }
    catch(...) {
        comment = info.str();
        message("%s\n", comment.c_str());
        error("Settings file could not be read!\n");
    }
    
    IN.close();
    
    message("\nReading settings file '%s' done (%u parameters)\n", name.c_str(),
            parsedParams.size());
    
    return true;
}

//------------------------------------------------------------------------------
/** read the arguments character by character until the "end" character,
 * which has to be a ws.
 */
string
TSimManager::readArguments(istream & IN, unsigned int& linecnt, char end,
                     bool& otherParams)
{
    string args, temp;
    char c;
    
    // remove all forgoing space
    IN >> ws;
    
    // character by character
    while (IN.get(c)) {
        if (isspace(c)) {
            IN.putback(c);
            if (!STRING::removeCommentAndSpace(IN, linecnt, 0, end)) { // tab
                otherParams = false; // that was the last argument
            }
            return args;
        }
        else {
            switch (c) {
                default: // read
                    args += c;
                    break;
                    
                case '\\': // line continuation
                    while (IN.get(c) && IN.good() && !IN.eof() && !isEOL(c, IN)) {
                    }   // remove rest of line
                    linecnt++;
                    break;
                    
                case '{':
                    IN.putback(c);
                    args += STRING::readBrackets2String(IN, linecnt, '}');
                    break;
                    
                case '(':
                    IN.putback(c);
                    args += STRING::readBrackets2String(IN, linecnt, ')');
                    break;
                    
                case '[':
                    IN.putback(c);
                    args += STRING::readBrackets2String(IN, linecnt, ']');
                    break;
                    
                case '\"':
                    IN.putback(c);
                    args += STRING::readUntilCharacter(IN, linecnt, '\"', true,
                                                       true);
                    break;
                    
            }
        }
    } //end while read args
    
    // remove trailing space
    unsigned int i = (unsigned int)args.size();
    while (i && isspace(args[i - 1])) {
        --i;
    }
    if (i != args.size()) args.erase(i);
    
    if (IN.eof()) otherParams = false;
    return args;
}

//------------------------------------------------------------------------------
/** add the parameters directly passed as arguemtn to the program to teh parameter
 * list and overwrite previous settings
 * read the arguments character by character until the "end" character,
 * which has to be a ws.
 */
void
TSimManager::add_programArgs(map<string, string>& args, map<string, vector<string> >& parsedParams)
{
    map<string, string>::iterator cur, end = args.end();
    for (cur = args.begin(); cur != end; ++cur) {
        add_parsedParam(cur->first, cur->second, parsedParams);
    }
}

//------------------------------------------------------------------------------
/** adds the new parameter with its argument to the map.
 * if the parameter has previously been set a warning is plotted
 * this function is used if parameters are passed as arguments to the executable
 */
void
TSimManager::add_parsedParam(string param, string argvect, map<string, vector<string> >& parsedParams)
{
    if (parsedParams.find(param) != parsedParams.end())
        warning("Parameter '%s' has previously been set!\n", param.c_str());
    parsedParams[param].clear(); // the previous porams are overwritten and not just added!
    parsedParams[param].push_back(argvect);
}

//------------------------------------------------------------------------------
/** adds the new parameter with its argument to the map.
 * if the parameter has previously been set a warning is plotted
 */
void
TSimManager::add_parsedParam(string& param, const vector<string>& argvect,
                             map<string, vector<string> >& parsedParams)
{
    if (parsedParams.find(param) != parsedParams.end())
        warning("Parameter '%s' has previously been set!\n", param.c_str());
    parsedParams[param] = argvect;
}

//------------------------------------------------------------------------------
/** adds the new keayword with its argument to the map.
 * if the parameter has previously been set a warning is plotted
 */
void
TSimManager::add_parsedKeys(string keyword,
                               const vector<string>& argvect,
                               map<string, vector<string> >& parsedKeywords)
{
    if (parsedKeywords.find(keyword) != parsedKeywords.end())
        warning("Keyword '%s' has previously been set!\n", keyword.c_str());
    parsedKeywords[keyword] = argvect;
}

//----------------------------------------------------------------------------------------
// print_info
// ----------------------------------------------------------------------------------------
void
TSimManager::set_iniFile_directory(string s)
{
    if(!s.empty()){
        if(s[s.length()-1]=='/' || s[s.length()-1]=='\\') _ini_file_directory = s;
        else _ini_file_directory = s + SEP;
    }
}

// ----------------------------------------------------------------------------------------
// get_fileName
// ----------------------------------------------------------------------------------------
/** returns the filename as it is if the path is absolute or adds the program directory path
 * if relative
 */
string
TSimManager::get_fileName(const string& name)
{
    if (name[0] == '/' || name[1] == ':') return name;    // absolute path
    return _ini_file_directory + name;    // relative path
}

//----------------------------------------------------------------------------------------
// build_records
// ----------------------------------------------------------------------------------------
/** generates the individual simulation parameter
 * if multiple sequential parameters are used, then all combiantions are run
 * keywordRecord has same length as paramRecord
 */
void
TSimManager::build_records_full(vector<map<string, string> >& paramRecord,
                                vector<map<string, string> >& keysRecord,
                                map<string, vector<string> >& initParams,
                                map<string, vector<string> >& initKeys)
{
    //print2file();     // de-comment it if all parameters should be dumped to file
    
    map<string,string> params;
    unsigned int RecNb = 1, SeqParam;
    vector<unsigned int> sequence;    //stores the number of args of each sequence parameters
    vector<unsigned int> sequenceKey; //as sequence, but for keywords
    vector<string> currSeqArg;
    
    string arg, NAME = STRING::get_file_basename(initParams["_settings_file"][0]);
    
    //then process the InitParams map to find sequence parameters:
    map< string,vector<string> >::iterator Pit;
    for(Pit = initParams.begin(); Pit != initParams.end(); ++Pit) {
        if(Pit->second.size() > 1) {
            sequence.push_back((unsigned int)Pit->second.size()); //fetch the number of the current sequence parameter:
            RecNb *= Pit->second.size();                          //increase the total number of simulations records:
        }
    }
    for(Pit = initKeys.begin(); Pit != initKeys.end(); ++Pit) {
        if(Pit->second.size() > 1) {
            sequenceKey.push_back((unsigned int)Pit->second.size()); //fetch the number of the current sequence parameter:
            RecNb *= Pit->second.size();                          //increase the total number of simulations records:
        }
    }
    
    if(RecNb > 100) warning("Your settings file has %u sequential simulations specified! Is that your aim, or is there a problem in your settings file?\n", RecNb);
    
    // for each sequence
    for(unsigned int i=0; i<RecNb; ++i) {
        //now build the simulation records with the right params!
        //the map 'param' will get all the params used for one simulation
        //it is then added to the list of simulations' parameters map
        
        SeqParam = 0;//used as index of the 'sequence' vector
        
        // for each input parameter
        build_records_fullII(params, initParams, sequence, RecNb, SeqParam, i);
        if(RecNb>1) params["filename"] = setFilename(NAME,i+1,RecNb,currSeqArg);   // sequence simulation
        paramRecord.push_back(params);

        // for each input keyword
        build_records_fullII(params, initKeys, sequenceKey, RecNb, SeqParam, i);
        keysRecord.push_back(params);
        
    } // end of each sequence of simulation
}

//----------------------------------------------------------------------------------------
// build_records
// ----------------------------------------------------------------------------------------
/** generates the individual simulation parameter sets
 * if multiple sequential parameters are used,the higest number of elements is searched
 * and then all other arguements are expanded as a vector to meet the length
 */
void
TSimManager::build_records_fullII(map<string,string>& params,
                                  map<string, vector<string> >& initMap,
                                  vector<unsigned int> sequence,
                                  unsigned int RecNb, unsigned int& SeqParam,
                                  unsigned int i)
{
    string NAME, arg;
    unsigned int ArgNb, BlockSize;
    vector<string> currSeqArg;
    params.clear();
    
    // for each input parameter
    map< string,vector<string> >::iterator Pit;
    for(Pit = initMap.begin(); Pit != initMap.end(); ++Pit) {
        // if it is the filename
        if (Pit->first == "filename") {
            if(RecNb>1) NAME = Pit->second[0];       // sequence simulations
            else  params["filename"] = Pit->second[0];
            continue;
        }
        
        // for all other parameters
        //get the number of arguments for the current parameter:
        switch(Pit->second.size()){
            case 0:
                //current param has no argument (bool type param) we give it value 1 (true)
                params[Pit->first] = "1";
                break;
                
            case 1:
                //the current param has one argument
                params[Pit->first] = Pit->second[0];
                break;
                
            default:
                //the current param is a sequence param (ArgNb > 1)
                //increase the index of the sequence parameter
                SeqParam++;
                //then compute the right argument to give to the current simulation record:
                BlockSize = RecNb;
                ArgNb = (unsigned int)Pit->second.size();
                
                //message("\n%s (%i): %s",Pit->first.c_str(), Pit->second.size(), Pit->second[ (i/BlockSize) % ArgNb ].c_str());
                
                for(unsigned int j=0;j < SeqParam;++j){
                    BlockSize /= sequence[j];
                }
                arg = Pit->second[ (i/BlockSize) % ArgNb ];
                if(arg[0]=='{' || arg[0]=='('){         // if matrix or temporal -> argument is too big, replaced by the number
                    currSeqArg.push_back(toStr(1+(i/BlockSize) % ArgNb, ArgNb));
                }
                else currSeqArg.push_back(arg);        // add argument
                params[Pit->first] = Pit->second[ (i/BlockSize) % ArgNb ];
                break;
        }
    }
}

//----------------------------------------------------------------------------------------
// build_records
// ----------------------------------------------------------------------------------------
/** generates the individual simulation parameter sets
 * if multiple sequential parameters are used,the higest number of elements is searched
 * and then all other arguements are expanded as a vector to meet the length
 */
void
TSimManager::build_records_reduced(vector<map<string, string> >& paramRecord,
                                   vector<map<string, string> >& keysRecord,
                                   map<string, vector<string> >& initParams,
                                   map<string, vector<string> >& initKeys)
{
    //print2file();     // de-comment it if all parameters should be dumped to file
    
    map<string,string> params;
    unsigned int RecNb = 1;
    vector<unsigned int> sequence;    //stores the number of args of each sequence parameters
    vector<unsigned int> sequenceKey; //same as "sequence" but for keywords
    vector<string> currSeqArg;
    
    string arg, NAME = STRING::get_file_basename(initParams["_settings_file"][0]);
    
    //then process the InitParams map to find sequence parameters:
    map< string,vector<string> >::iterator Pit;
    for(Pit = initParams.begin(); Pit != initParams.end(); ++Pit) {
        if(Pit->second.size() > 1) {
            sequence.push_back((unsigned int)Pit->second.size());           // fetch the number of the current sequence parameter:
            if(RecNb<Pit->second.size()) RecNb=(unsigned int)Pit->second.size();   //increase the total number of simulations records:
        }
    }
    
    if(RecNb > 100) warning("Your settings file has %u sequential simulations specified! Is that your aim, or is there a problem in your settings file?\n", RecNb);
    
    // for each sequence
    for(unsigned int i=0; i<RecNb; ++i) {
        //now build the simulation records with the right params!
        //the map 'param' will get all the params used for one simulation
        //it is then added to the list of simulations' parameters map
        
        // for each input parameter
        build_records_reducedII(params, initParams, sequence, RecNb, i);
        if(RecNb>1) params["filename"] = setFilename(NAME,i+1,RecNb,currSeqArg);   // sequence simulation
        paramRecord.push_back(params);
        
        // for each input keyword
        build_records_reducedII(params, initKeys, sequenceKey, RecNb, i);
        keysRecord.push_back(params);
    } // end of each sequence of simulation
}

//----------------------------------------------------------------------------------------
// build_records
// ----------------------------------------------------------------------------------------
/** generates the individual simulation parameter sets
 * if multiple sequential parameters are used,the higest number of elements is searched
 * and then all other arguements are expanded as a vector to meet the length
 */
void
TSimManager::build_records_reducedII(map<string,string>& params,
                                  map<string, vector<string> >& initMap,
                                  vector<unsigned int> sequence,
                                  unsigned int RecNb, unsigned int i)
{
    string NAME, arg;
    unsigned int ArgNb;
    vector<string> currSeqArg;
    params.clear();
    
    // for each input parameter
    map< string,vector<string> >::iterator Pit;
    for(Pit = initMap.begin(); Pit != initMap.end(); ++Pit) {
        // if it is the filename
        if (Pit->first == "filename") {
            if(RecNb>1) NAME = Pit->second[0];       // sequence simulations
            else  params["filename"] = Pit->second[0];
            continue;
        }
        
        // for all other parameters
        //get the number of arguments for the current parameter:
        switch(Pit->second.size()){
            case 0:
                //current param has no argument (bool type param) we give it value 1 (true)
                params[Pit->first] = "1";
                break;
                
            case 1:
                //the current param has one argument
                params[Pit->first] = Pit->second[0];
                break;
                
            default:
                // expand any parameter to the max size
                ArgNb = (unsigned int)Pit->second.size();
                arg = Pit->second[i % ArgNb];
                params[Pit->first] = arg;
                
                if(arg[0]=='{' || arg[0]=='('){         // if matrix or temporal -> argument is too big, replaced by the number
                    currSeqArg.push_back(toStr(1+i%ArgNb, ArgNb));
                }
                else currSeqArg.push_back(arg);        // add argument
                break;
        }
    }
}

// ----------------------------------------------------------------------------------------
// setFilename
// ----------------------------------------------------------------------------------------
string
TSimManager::setFilename(string& fstring,const unsigned int& sim, const unsigned& nbSims, vector<string>& args)
{
    string out, tail;
    int cur, digits;
    unsigned int next;
    unsigned int index, i, size = (unsigned int)args.size();
    
    // control array that all parameters were used, i.e. that the file name is unique!
    int* usedParams = new int[size];
    for(i=0; i<size; ++i){
        usedParams[i] = 0;    // set it to false
    }
    
    // parse the name
    cur  = 0;
    next = (unsigned int)fstring.find_first_of('%');
    while(next != (unsigned int)string::npos){
        out += fstring.substr(cur, next-cur);
        tail = fstring.substr(next+1, string::npos);
        index = (unsigned int) strtol(tail.c_str(),NULL,10);          // get the param number
        if(index > args.size()) error("too many arguments in filename!\n");
        digits = getNbDigits(index);
        
        out += args[index-1];
        usedParams[index-1] = 1;    // this param was used
        cur = next + digits + 1;
        next = (unsigned int)fstring.find('%',cur);  // get the next position
    }
    out += fstring.substr(cur, next-cur);
    
    // control if all params were used
    for(i=0; i<size; ++i){
        if(!usedParams[i]) break;
    }
    delete[] usedParams;
    
    if(i != size) out += "-" + toStr(sim, nbSims);    // add the simulation number if the name is not unique
    
    return out;
}


