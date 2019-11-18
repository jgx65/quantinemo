/** @file tsim_manager.h
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


#ifndef tsim_managerH
#define tsim_managerH

#include <map>
#include <vector>
#include <string>
#include <iostream>

using namespace std;

class TSimManager{
public:
    vector<map<string, string> > _simSettings; // simSettings[sim][param][arg]
    vector<map<string, string> > _keyParams;    // _keyParams[sim][keay][arg]
    unsigned int _nbThreadsTot;              // total number of threatds as specified by the user

    unsigned int _nbSims;
    string _exe_directory;
    
    void readInputFiles(int ARGC, char **ARGV,
                        vector<map<string, string> >& paramRecord,
                        vector<map<string, string> >& keysRecord);
    
    void program_args(map<string, string>& args, vector<string>& inputfiles);
    
    bool read_settings_file(string name, map<string, vector<string> >& parsedParams,
                            map<string, vector<string> >& keywords);
    string readArguments(istream & IN, unsigned int& linecnt, char end,
                         bool& otherParams);
    
    void add_programArgs(map<string, string>& args, map<string, vector<string> >& parsedParams);
    
    void add_parsedKeys(string keyword,
                        const vector<string>& argvect,
                        map<string, vector<string> >& parsedKeys);
    void add_parsedParam(string param, string argvect,
                         map<string, vector<string> >& parsedParams);
    void add_parsedParam(string& param, const vector<string>& argvect,
                         map<string, vector<string> >& parsedParams);
    
    void build_records_full(vector<map<string, string> >& paramRecord,
                            vector<map<string, string> >& keyRecord,
                            map<string, vector<string> >& initParams,
                            map<string, vector<string> >& initKeys);
    void build_records_reduced(vector<map<string, string> >& paramRecord,
                               vector<map<string, string> >& keysRecord,
                               map<string, vector<string> >& initParams,
                               map<string, vector<string> >& initKeys);
    
    void build_records_fullII(map<string,string>& params,
                              map<string, vector<string> >& initMap,
                              vector<unsigned int> sequence,
                              unsigned int RecNb, unsigned int& SeqParam,
                              unsigned int i);
    void build_records_reducedII(map<string,string>& params,
                              map<string, vector<string> >& initMap,
                              vector<unsigned int> sequence,
                              unsigned int RecNb, unsigned int i);
    
    string setFilename(string& fstring,const unsigned int& sim,
                       const unsigned& nbSims, vector<string>& args);
    
    unsigned int get_nbThreads();
    void printVersion();
    string _ini_file_directory;

public:
    TSimManager(int ARGC, char **ARGV);
    ~TSimManager(){}
    void run();
    void print_help(ostream& os, unsigned int wide1=20, unsigned int wide2=60,
                    char fill='.', unsigned int importance=5, string arg="");
    
    void set_iniFile_directory(string s);
    string get_fileName(const string& name);
    
    
//    void run_sims(unsigned int from, unsigned int to, unsigned int nbThread);
};

void
run_sims_withinThread(TSimManager* ptr, unsigned int from, unsigned int to, unsigned int nbThread);

#endif /* defined(tsim_managerH) */
