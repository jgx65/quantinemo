/** @file fileservices.cpp
 *
 *   Copyright (C) 2006 Frederic Guillaume    <guillaum@zoology.ubc.ca>
 *   Copyright (C) 2008 Samuel Neuenschwander <samuel.neuenschwander@unil.ch>
 *   Copyright (C) 2018 Frederic Michaud <frederic.michaud@unil.ch>
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


#include "fileservices.h"
#include "simcomponent.h"
#include "filehandler.h"
#include "tmetapop.h"

using namespace std;


// ----------------------------------------------------------------------------------------
// attach
// ----------------------------------------------------------------------------------------
void FileServices::attach ( FileHandler* FH )
{
    Service::attach(FH);
    _children.push_back(FH);
    FH->set_service(this);
    _popPtr->set_service(this);
}
// ----------------------------------------------------------------------------------------
// load
// ----------------------------------------------------------------------------------------
void FileServices::load ( TSimComponent* sc ) {
    sc->loadFileServices(this);
}

// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
bool FileServices::check_file_overwrite ()
{
    char yn;
    bool ok;
    list< FileHandler* >::iterator HIT;
    string name;
     // while the filename is not properly set
    do{
        ok = true;
        for(HIT = _children.begin(); HIT != _children.end(); ++HIT) {
            ok &= (*HIT)->init();
        }
        
        if(!ok && !_overwriteFiles) {
            message(" Do you want to overwrite all the files that use it ? (y/n/s(kip)): \n");
            cin>>yn;
            
            switch(yn) {
                case 'y':
                    ok = true;
                    break;
                case 's':
                    message(" Skipping this simulation\n");
                    return false;
                default: {
                    message(" Please give a new output filename: ");
                    cin>>name;
                }
            }
        }
    } while(!ok && !_overwriteFiles);
    
    return true;
}






