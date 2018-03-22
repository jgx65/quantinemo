/** @file lce_regulation.h
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


#ifndef lce_regulationH
#define lce_regulationH

#include "lifecycleevent.h"


// LCE_Regulation
//
/** Regulates the population after dispersal.
 * Moves individuals from post-dispersal to adult containers according to Patch capacities.
 * Sets the age flag of the population to ADULTS.
 **/
class LCE_Regulation: public LCE{
protected:
	void (LCE_Regulation::* regulation)();     // pointer to the correct regulation function
    
	string  _ageStr;
	age_idx _age;
    
#ifdef _DEBUG
    unsigned int _ex_cnt, _col_cnt, _ph_cnt;
#endif
    
public:
    
	LCE_Regulation(age_idx age = ADLTx, int rank = my_NAN) :
    LCE("regulation_"+ (string)(age==ADLTx ? "adults":"offspring"),
         "population size regulation of "+ (string)(age==ADLTx ? "adults":"offspring"),"",rank) {
		_age = age;
		_ageStr = age==ADLTx ? "adults":"offspring";
		add_parameter("regulation_model_"+_ageStr,INT2,false,0,1,"0", false,
                            "Population size regulation of "+_ageStr+":\n" \
                            "  0: no regulation\n" \
                            "  1: random regulation (population size is down-regulated to carrying capacity)",3);
	}
    
	virtual LCE_Regulation* clone () {return new LCE_Regulation();}
    
    virtual ~LCE_Regulation( ) { }
    
	virtual void execute ();
    virtual bool init(TMetapop* popPtr);
	virtual void regulation_neutral();
	virtual void regulation_coalescence();
	virtual void regulation_fitness_patch()  {_popPtr->regulate_selection_fitness_patch(_age);}
	virtual void regulation_fitness_metapop(){_popPtr->regulate_selection_fitness_metapop(_age);}
	virtual void regulation_fitness_hard()   {_popPtr->regulate_selection_fitness_hard(_age);}
	virtual void drawSuccessfullIndividuals(TPatch* curPatch, const unsigned int& K, const sex_t& SEX);
    virtual void drawUnSuccessfullIndividuals(TPatch* curPatch, const unsigned int& K, const sex_t& SEX);
    
    //SimComponent overrides:
    virtual void loadFileServices ( FileServices* loader ) {}
	virtual void loadStatServices ( StatServices* loader ) {}
	virtual age_t removeAgeClass ( ) {return 0;}
	virtual age_t addAgeClass ( ) {return 0;}
	virtual age_t requiredAgeClass () {return _age==ADLTx ? ADULTS : OFFSPRG;}
};
#endif /* defined(__quantiNemo2__lce_regulation__) */
