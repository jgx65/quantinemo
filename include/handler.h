/** @file handler.h
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

#ifndef handlerH
#define handlerH

#include <fstream>
using namespace std;

class TMetapop;

/**Service handler (an observer).
 * Interface for the FileHandler and StatHandler classes.
 * Implements the observer design pattern
 */
class Handler {
    
public:
    
    Handler(): _popPtr(0){}
    
    /**Inits state. */
    virtual bool init() = 0;
    
    /**Updates the handler state.  */
    virtual void update(ostream& FH){}
    virtual void update(){}
    
    /**Updates the current populated patches (used only to compute stats).  */
    virtual void update_patch_states(){};
    
    //virtual char getName() = 0;
    
    /**Link to the current population, set through the link to the StatService.*/
    TMetapop* _popPtr;

    
    virtual ~Handler(){}
};
#endif //HANDLER_H
