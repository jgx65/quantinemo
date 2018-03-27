/** @file fileservices.h
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

#ifndef fileservicesH
#define fileservicesH

#include <list>
#include "service.h"
#include "param.h"

class TMetapop;

class FileHandler;

/**A class to manage the files associated with each components of the simulation.
 
 Implements the Observer design pattern (is the concrete subject), stores the base filename of the simulation and
 updates the replicate filenames. It also performs files checking and saves the simulation parameters on init.
 */

class FileServices : public Service {
    
private:
    /**a pointer to the current Metapop*/
    TMetapop*    _popPtr;
    
    /**the list of the FileHandler's registered by the SimComponent*/
    list< FileHandler* > _children;
    
    
    /**the list of the simulation parameters*/
    list< ParamSet* > _params;
    
    /** should existing files be overwriten without asking? */
    bool      _overwriteFiles;
    int       _logfile_type;       // 0: as input, 1: minimal, 2: maximal
    
    
public:
    
    FileServices ( ) : _popPtr(0){ }
    
    virtual ~FileServices ( ) { }
    
    virtual bool init ( ) {return false;}
    
    /**Checks if files with _basename already exist and save the simulation parameters in log files.
     * @param params a ref to the list of the current parameters of the simulation
     * @return true if the files check is ok
     * @return false if the user wants to skip this simulation
     */
    bool check_file_overwrite ();
    void set_params(list<ParamSet* >& params){_params=params;}
    /**
     * @return the pointer to current Metapop
     */
    virtual TMetapop*   get_pop_ptr ( )      {return _popPtr;}
    
    /**Sets the Metapop reference.*/
    virtual void set_pop_ptr (TMetapop* pop) {_popPtr=pop;}
    
    /**Saves the current simulation parameters in log files.
     * @param params a ref to the list of the current parameters of the simulation
     * @param FH the handle to the log file
     */
    void save_simparams(list< ParamSet* >&  params, ostream& FH);
    
    /**
     * @return the list of the current parameters of the simulation
     */
    list< ParamSet* >& get_params() {return _params;};
    
    list< FileHandler* >& get_children() {return _children;}
    
    /**
     * @return the current replicate file name
     */
    
    /**Tells the SimComponent to load its file handlers.
     *  @param sc the SimComponent
     */
    virtual void load ( TSimComponent* sc );
    
    /**Attaches the FileHandler to the current list (_children) of the FileServices.
     *  @param FH the FileHandler
     */
    virtual void attach ( FileHandler* FH );
    
    /**Clears the list of FileHandlers. */
    virtual void reset ( ) {
        //   for(list< FileHandler* >::iterator HIT  = _children.begin(); HIT != _children.end(); ++HIT) {
        //     (*HIT)->reset();
        //   }
        Service::reset_observers();
        _children.clear();
    }
};
#endif //FILESERVICES_H

