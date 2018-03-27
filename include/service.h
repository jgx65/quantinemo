/** @file service.h
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

#ifndef serviceH
#define serviceH

#include <list>
#include "handler.h"
#include "functions.h"
using namespace std;


class TSimComponent;

/**Interface for the simulation services (files and stats).
 * Implements the observer pattern. Notify the observers (Handler) to update their state.
 * Contains the observer list. Provides interface to attach the observers to the service.
 **/
class Service {
private:
    /**The list of observers. */
    list<Handler*> _observers;
    
public:
    
    Service( ) { }
    
    virtual ~Service( ) {  _observers.clear();}
    
    /**Inits internals. */
    virtual bool init ( ) = 0;
    
    /**Notifies all observers to update their state. */
    virtual void notify() {
#ifdef _DEBUG
        message(" %i ...\n",_observers.size());
#endif
        
        unsigned int i=0;
        for(list< Handler* >::iterator HIT = _observers.begin(); HIT != _observers.end(); ++HIT, ++i) {
            (*HIT)->update();
        }
    }
    
    /**Notifies all observers to update their state. */
    virtual void notify (ostream& FH) {
#ifdef _DEBUG
        message(" %i ...\n",_observers.size());
#endif
        
        for(list< Handler* >::iterator HIT = _observers.begin(); HIT != _observers.end(); ++HIT) {
            (*HIT)->update(FH);
        }
    }
        
    /**Interface to used by a simulation component to load its obervers onto a service provider.*/
    virtual void load ( TSimComponent* sc ) = 0;
    
    /**Adds an observer to the list. */
    void attach ( Handler* h ) {
        _observers.push_back(h);
    }
    
    /**Clears the observers list. */
    virtual void reset_observers ( ) {
        _observers.clear();
    }
    
};

#endif //SERVICE_H






