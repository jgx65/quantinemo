/** @file lifecycleevent.h
 *
 *   Copyright (C) 2006 Frederic Guillaume    <guillaum@zoology.ubc.ca>
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

#ifndef lifecycleH
#define lifecycleH

#include "simcomponent.h"
#include "tmetapop.h"

/**Base class of the Life Cycle Events, declares the LCE interface.
 * A Life Cycle Event (LCE) is a kind of population operator that modifies the population state. It declares the execute() method
 * that is called during the life cycle loop (\c Metapop::Generation_LOOP()). No interface is given to manage the periodicity of the event, this
 * task is left to the derived class. Each LCE has a link to the current Metapop instance which is set by a call to \c LCE::init()
 * during the life cycle setup (\c Metapop::setLifeCycle()).
 * Each LCE has a name by which it is called from the user's defined parameter file to specify when in the life cycle that particular event
 * has to be executed. The position in the life cycle is the LCE's rank (i.e. the value given in input to the parameter "name"). As only one
 * rank value is allowed per LCE, it is executed only once in the life cycle, this is a limitation that might be removed in future
 * versions if need be!
 * This interface offers the possibility to link the LCE with one particular trait so that it is easily and efficiently accessed
 * throughout the simulation. Links (the index of the trait in the individual's trait table) are set by a call to \c IndFactory::getTraitIndex()
 * accessed through the Metapop pointer at initialization.
 
 * \b Note: If more than one trait is used by your LCE, think about overloading the init() procedure to allow the
 * registration of several trait links within your derived LCE.
 
 * \b Note: the LCEs \a must be inited \a after a call to IndFactory::makePrototype() has been issued so that the trait links can be set!
 */
class LCE: public TSimComponent {
private:
    /**The param name to be read in the init file.*/
    string _event_name;
    
protected:
    /**The ptr to the current Metapop.*/
    TMetapop* _popPtr;
    
    /**The name of the linked trait.*/
    string        _linkedTraitType;
    
    /** default rank of the LCE */
    int           _default_rank;
    
    /**The index in the individual's trait table of the linked trait. A value of -1 means the link is broken.*/
    vector<int>   _linkedTraitIndex;
    unsigned int  _linkedTraitSize;  // size of the _linkedTraitIndex vector
    
public:
    /**Cstor.
     @param name the name of the LCE as it must appear in the parameter input file
     @param trait_link the name of the linked trait used by this LCE
     @callgraph
     */
    LCE (string name, string name_long, string trait_link, int rank)
    : _popPtr(0), _linkedTraitType(trait_link), _default_rank(rank){
        
        set_event_name(name, name_long);
    }
    
    virtual ~LCE ( ) { }
    
    /**Sets the pointer to the current Metapop and the trait link if applicable.
     * DEV: Don't forget to explicitly call this function (\c LCE::init()) when you overload it!
     * @param popPtr the pointer to the current instance of \c Metapop
     * @callgraph
     */
    virtual bool init(TMetapop* popPtr){
        _popPtr = popPtr;
        if(!_linkedTraitType.empty()) {
            _linkedTraitIndex = popPtr->getTraitIndex(_linkedTraitType.c_str());
            _linkedTraitSize = (int) _linkedTraitIndex.size();
            if(_linkedTraitIndex.empty()){
                error("Trait '%s' not present in prototypes, cannot include life cycle event '%s'!\n",_linkedTraitType.c_str(),_event_name.c_str());
            }
        }
        else _linkedTraitSize = 0;
        return true;
    }
    
    /**Set the name of the event (name of the ParamSet) and add the corresponding parameter to the set.   */
    virtual void set_event_name (string name, string name_long) {
        _event_name = name;
        set_paramset(name, name_long, false, _popPtr);
        add_parameter(name,INT2,false,my_NAN,my_NAN,toStr(_default_rank));
    }
    
    /**Accessor to the LCE's name.*/
    virtual string& get_event_name ( )              {return _event_name;}
    
    /**Accessor to the LCE rank in the life cycle.*/
    virtual int get_rank ( ) {return (int)get_parameter_value(_event_name);}
    
    /**Accessors for the population pointer
     * @param popPtr The pointer to the current Metapop
     */
    virtual void set_pop_ptr (TMetapop* popPtr) {_popPtr=popPtr;}
    virtual TMetapop*   get_pop_ptr ( ) {return _popPtr;}
    virtual vector<int>* get_linkedTraitIndex()   {return &_linkedTraitIndex;}
    
    ///@name LCE interface
    ///@{
    /**Execute the event on the pop. */
    virtual void  execute () = 0;
    
    /**Cloning interface. */
    virtual LCE*  clone () = 0;
    
    /**Removes the required age-class(es) flag from the current Metapop age-class flags.*/
    virtual age_t removeAgeClass () = 0;
    
    /**Adds the required age-class(es) flag to the current Metapop age-class flags.*/
    virtual age_t addAgeClass () = 0;
    
    /**Specifies what age-classes are required by the LCE to execute.*/
    virtual age_t requiredAgeClass () = 0;
    ///@}
    
    virtual void executeBeforeEachReplicate(const unsigned int& rep){}
    virtual void executeAfterEachReplicate(const unsigned int& rep){}
    
    virtual void temporal_change(const unsigned int& gen){}
    virtual void executeBeforeEachGeneration(const unsigned int& gen){}
    virtual void executeAfterEachGeneration(const unsigned int& gen){}
    
    
    
    
};

#endif //LIFECYCLEEVENT_H

