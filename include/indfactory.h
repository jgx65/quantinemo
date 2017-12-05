/** @file indfactory.h
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

#ifndef indfactoryH
#define indfactoryH

#include <map>
#include "tindividual.h"
#include "types.h"
#include "ttrait.h"
#include <vector>


class TSelection;

/**Factory of Individual, stores the individual prototype and the trait prototypes, manages the individual garbage collector.
 * Provides methods to generate new individuals within the Metapop. Each new individual is created by cloning a prototype itself
 * created at the simulation setup. New individuals are decorated with the appropriate traits as set by the trait prototypes and
 * receives a unique ID (unique within a simulation).
 **/
class IndFactory {
protected:
    /**Map of the trait prototypes.*/
    map< string,TTraitProto* > _protoTraits;
    
    TGenomeProto* _protoGenome; 		// pointer to the genetic map (has to be deleted)
    
    TMetapop* _popPtr;
    
    /**Table containing the absolute index across all types of traits*/
    map< string, int > _TraitsIndex;
    
    /** vector containing a pointer to the TraitPrototype (allows to access the prototypes directly) */
    vector<TTraitProto*> _protoTraitsVector;
    
    
    /**Garbage collector for unused Individual's.*/
    vector<TIndividual*> RecyclingPOOL;
    
public:
    
    /**The individuals prototype used to create any new individual in a simulation.*/
    TIndividual _protoIndividual;
    
    IndFactory ( ):_popPtr(0) {
        _protoGenome = new TGenomeProto();
    }
    
    virtual ~IndFactory ( ) {
        for(unsigned int i=0; i<RecyclingPOOL.size(); ++i){
            delete RecyclingPOOL[i];
        }
        RecyclingPOOL.clear();
        _protoIndividual.genome.clear();		// has to be deon before _protoGenoem is deleted
        delete _protoGenome;
    }
    
    /**Put an individual in the recycling pool.*/
    void recycle(TIndividual* ind) {
        assert(ind);
        RecyclingPOOL.push_back(ind);
    }
    
    /**Creates the individuals prototype from the selected trait prototypes.
     Resets the individual's ID counter to 0 and sets the traits index table.
     @callgraph
     @param TTlist the list of the current trait prototype selected from the current simulation parameters.
     **/
    void                    makePrototype               (map< string,TTraitProto* >& TTlist, TMetapop* pMetapop);
    
    /**Resets the individuals prototype between replicates from the selected trait prototypes.
     Resets the individual's ID counter to 0 and sets the traits index table.
     @callgraph
     @param TTlist the list of the current trait prototype selected from the current simulation parameters.
     **/
    void                    resetPrototype               ( );
    
    
    /**Creates a blank individual which has to be decorated.
     * ID is set and new traits are allocated but no genetic data is created. Sex has to be set too.
     * @callgraph
     **/
    TIndividual*             getNewIndividual() {return makeNewIndividual(NULL,NULL,MAL,NULL);}
    
    /**Creates an individual with pointers to parents, sex and home ID set but no genetic data.
     * No inheritance or mutations on the trait sequences are done.
     * @callgraph
     * @param mother ptr to the mother
     * @param father ptr to the father
     * @param sex gender of the individual
     * @param homepatch ID of the Patch where this individual is born, usually the current position in the Patch array
     **/
    TIndividual*    makeNewIndividual           (TIndividual* mother, TIndividual* father, sex_t sex, TPatch* homepatch);
    TIndividual*    copyIndividual(TIndividual* oldInd);

    
    /**Completely creates an individual with inheritance and mutations on all traits.
     * @callgraph
     * @param mother ptr to the mother
     * @param father ptr to the father
     * @param sex gender of the individual
     * @param homepatch ID of the Patch where this individual is born, usually the current position in the Patch array
     **/
    TIndividual*     makeOffsprg (TIndividual* mother, TIndividual* father, sex_t sex, TPatch* homepatch);
    
    /**Individual prototype accessor.*/
    const TIndividual*       getIndividualProtoype       ( )   {return &_protoIndividual;}
    
    /**Accessor to the list of TTraitProto's.*/
    map< string,TTraitProto* >& getTraitPrototypes  ( )  {return _protoTraits;}
    TTraitProto* getFirstPrototype(string type);
    
    TTraitProto& getTraitPrototype(int i) {
        assert(i>=0 && i<(int)_protoTraitsVector.size());
        return *_protoTraitsVector[i];
    }
    
    unsigned int getTraitPrototypeSize() {return (unsigned int)_protoTraitsVector.size();}
    
    /**Gives the index of trait with \a type.
     @param type the type of the trait*/
    vector<int> getTraitIndex (string type);
};

#endif

