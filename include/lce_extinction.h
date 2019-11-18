/** @file lce_extinction.h
*
*   Copyright (C) 2006 Frederic Guillaume    <guillaum@zoology.ubc.ca>
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
*///---------------------------------------------------------------------------

#ifndef lce_extinctionH
#define lce_extinctionH
//---------------------------------------------------------------------------

#include "lifecycleevent.h"


// LCE_Extinction
//
/**Randomly extincts patches according to the extinction rate parameter.
 * Sets the patches extinction flag accordingly.**/
class LCE_Extinction: public LCE
{
protected:
  /**Patch extinction probability.*/
	double*       _Xtion_rate;            // array with extinction rates
	unsigned int  _Xtion_rate_size;       // number of such values (size of array)
	double*       _survival_rate[2];      // array with survival rate (<1: ratio, >=1 absolute) per sex
	unsigned int  _survival_rate_size[2]; // size of the array _survival_rate per sex

public:

	LCE_Extinction(int rank = my_NAN);

	virtual ~LCE_Extinction( ) {
		if(_Xtion_rate) delete[] _Xtion_rate;
		if(_survival_rate[FEM]) delete[] _survival_rate[FEM];
		if(_survival_rate[MAL]) delete[] _survival_rate[MAL];
	}

	virtual bool init (TMetapop* popPtr);

	virtual void execute (){(this->*extinction_func_ptr)();}

	// extinction rate
	bool set_extinction_rate();
	void (LCE_Extinction::* extinction_func_ptr) ();
	void extinction_variable_total ();              // patch specific and entire extinction
	void extinction_homogenous_high_total ();       // identical among patches (faster) and entire extinction
	void extinction_homogenous_low_total ();        // identical among patches (faster) and entire extinction
	void extinction_variable_partial ();            // patch specific and partial extinction
	void extinction_homogenous_high_partial ();     // identical among patches (faster) and partial extinction
	void extinction_homogenous_low_partial_1sex (); // identical among patches (faster) and partial extinction
	void extinction_homogenous_low_partial_2sex (); // identical among patches (faster) and partial extinction

	// survival rate
	bool set_survival_rate();
	bool set_survival_rate(Param* param, sex_t SEX, double correction = 1);
	unsigned int (LCE_Extinction::* survivors_func_ptr[2]) (TPatch* curPatch, const sex_t & SEX);
	unsigned int survivors_relative_var (TPatch* curPatch, const sex_t & SEX);   // survivors are only relativly defined
	unsigned int survivors_absolute_var (TPatch* curPatch, const sex_t & SEX);   // survivors are only absolutivly defined
	unsigned int survivors_relative_const (TPatch* curPatch, const sex_t & SEX); // survivors are only relativly defined
	unsigned int survivors_absolute_const (TPatch* curPatch, const sex_t & SEX); // survivors are only absolutivly defined
	unsigned int survivors_mixed (TPatch* curPatch, const sex_t & SEX);      // survivors are relativly and absolutivly defined
	unsigned int survivors_none (TPatch* curPatch, const sex_t & SEX);       // all individuals  of this sex die
	unsigned int survivors_all (TPatch* curPatch, const sex_t & SEX);        // all survive of this sex (nothing to do)

	unsigned int survivors_relative_var_coal (TPatch* curPatch, const sex_t & SEX);   // survivors are only relativly defined
	unsigned int survivors_absolute_var_coal (TPatch* curPatch, const sex_t & SEX);   // survivors are only absolutivly defined
	unsigned int survivors_relative_const_coal (TPatch* curPatch, const sex_t & SEX); // survivors are only relativly defined
	unsigned int survivors_absolute_const_coal (TPatch* curPatch, const sex_t & SEX); // survivors are only absolutivly defined
	unsigned int survivors_mixed_coal (TPatch* curPatch, const sex_t & SEX);      // survivors are relativly and absolutivly defined

    void executeBeforeEachGeneration(const unsigned int& gen){}
    void temporal_change(const unsigned int& gen);

	virtual LCE_Extinction* clone ( ) {return new LCE_Extinction();}

  //SimComponent overrides:
  virtual void loadFileServices ( FileServices* loader ) {}
  virtual void loadStatServices ( StatServices* loader ) {}
  virtual age_t removeAgeClass ( ) {return 0;}
  virtual age_t addAgeClass ( ) {return 0;}
  virtual age_t requiredAgeClass () {return ADULTS | OFFSPRG;}
};
#endif
