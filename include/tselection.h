/** @file tselection.h
 *
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
//---------------------------------------------------------------------------

#ifndef tselectionH
#define tselectionH

#include "simcomponent.h"
#include "tpatchfitness.h"
#include "tequation.h"

class TTraitQuantiProto;
class LCE_Breed;
class LCE;
class TMetapop;
class TPatch;
class TIndividual;
class FileServices;
class StatServices;
class TSelectionTrait;


class TSelection: public TSimComponent {
private:
	// object pointers
	TMetapop*            _popPtr;
//	TTraitQuantiProto**     _pQuantiProto;


	TPatchFitness*      _fit[2];          // pointer to obejcts (male/female) containg the fitness, and individuals

	double*             _aPheno;          // temp array to stoe the phenotype of all traits
    
    // frequency dependent selection
	TTraitQuantiProto**     _selTrait_fitnessDependent; // array for each trait with frequency dependent selection
    unsigned int        _selTrait_fitnessDependentSize;
    
    // sex allocation trait
    double     _fem_sex_allocation_value;  // relative trait index (among qtraits) defining the female sex allocation proportion or the proportion
    

	// linked traits
	vector<int>         _vTraits;         // vector of indexes to the qTraits
	unsigned int        _vTraitsSize;     // number of qTraits

	// number of populations (used for the destructor as metapop is already destroit)
	int                 _nbPops;

	TSelectionTrait**   _selTrait;        // array for each selective trait

//	double              _Ve_prop;         // Proportion of the used Ve (0: only natal patch; 1: only mating patch(default))
//	int                 _Ve_model;        // environmental model: 0: Ve directly, 1: h2; 2: h2 at every gen; 3: H2; 4: H" at every gen

	double _optima_sd;
	double _intensity_sd;
	double _growth_rate_sd;
	double _max_growth_sd;
	double _symmetry_sd;
	double _min_sd;
	double _max_sd;
	double _selCoef_sd;

	unsigned int  _selection_level;    // 0: soft selection: fitness relative to patch
	// 1: metapop selection (soft to metapop): fitness relative to metapopulation
	// 2: hard selection (soft to hard): absolute fitness
	// 3: hard selection (metapop to hard): absolute fitness

	double _selection_level_coef;      // 0: soft selection; 1: hard selection

	unsigned int  _selection_position; // when does selection acts?:
	//      0: reproductive success (fitness of adults, nbOffs in relation to adults)
	//      1: reproductive success (fitness of offsprings, nbOffs in relation to adults))
	//      2: pre-dispersal regulation (offspring)
	//      3: post-dispersal regulation (adults)
	// 		 4: no selection at all (or when no quantitative traits are specified)

	bool  _Ve_mean_set;     // is there a constant environmental factor?
	bool  _Ve_var_set;      // is there an environmental effect set for each patch separately?
	double _get_fitness_multiplicative(TIndividual* ind);

public:
	// constructor
	TSelection(){}
	TSelection(TMetapop* popPtr){init(popPtr);}

	bool init(TMetapop* popPtr);
	bool init2(TMetapop* popPtr);
	void reset_selectionTypes();

	// destructor
	~TSelection ( );

	unsigned int get_selection_level()             {return _selection_level;}
	void         set_selection_level(const unsigned int& i) {_selection_level = i;}
	unsigned int get_selection_position()          {return _selection_position;}

	//////////////////////////////////////////////////////////////////////////////
	// functions to get an individual depending on its fitness
	// random not taking into account the fitness
	inline TIndividual* get_RAND_noFit(sex_t sex){
		return _fit[sex]->get_RAND_noFit();
	}
	inline TIndividual* get_RAND_noFit_index(sex_t sex, unsigned int& index){
		return _fit[sex]->get_RAND_noFit_index(index);
	}
	inline TIndividual* get_RAND_noFit(sex_t sex, unsigned int nb){
		return _fit[sex]->get_RAND_noFit_subset(nb);
	}

	// most fittest
	inline TIndividual* get_mostFit(sex_t sex){
		return _fit[sex]->get_mostFit();
	}
	inline TIndividual* get_RAND_mostFit(sex_t sex){
		return _fit[sex]->get_RAND_mostFit();
	}
	inline TIndividual* get_RAND_mostFit_index(sex_t sex, unsigned int& index){     // changes the index
		return _fit[sex]->get_RAND_mostFit_index(index);
	}
	inline TIndividual* get_RAND_mostFit_of_mostFit(sex_t sex, const unsigned int& nb){
		return _fit[sex]->get_RAND_mostFit_of_mostFit(nb);
	}
	inline TIndividual* get_RAND_mostFit_of_RAND_mostFit(sex_t sex, const unsigned int& nb){
		return _fit[sex]->get_RAND_mostFit_of_RAND_mostFit(nb);
	}

	// less fittest
	inline TIndividual* get_lessFit(sex_t sex){
		return _fit[sex]->get_lessFit();
	}
	inline TIndividual* get_RAND_lessFit(sex_t sex){
		return _fit[sex]->get_RAND_lessFit();
	}
	inline TIndividual* get_RAND_lessFit_index(sex_t sex, unsigned int& index){
		return _fit[sex]->get_RAND_lessFit_index(index);
	}
	inline TIndividual* get_RAND_lessFit_of_lessFit(sex_t sex, const unsigned int& nb){
		return _fit[sex]->get_RAND_lessFit_of_lessFit(nb);
	}
	inline TIndividual* get_RAND_lessFit_of_RAND_lessFit(sex_t sex, const unsigned int& nb){
		return _fit[sex]->get_RAND_lessFit_of_RAND_lessFit(nb);
	}

	//////////////////////////////////////////////////////////////////////////////
	// returns the object containg the fitness and the individuals
	// (the array has to be deleted afterwards if it is detached)
	TPatchFitness* get_FitnessObject(sex_t sex, bool detachObject = false){
		assert(_fit[sex]);
		TPatchFitness* temp = _fit[sex];
		if(detachObject) _fit[sex] = NULL;
		return temp;
	}

	// adds any fitness arrays, which can then be used. It will be deleted later!
	void set_FitnessObject(sex_t sex, TPatchFitness* obj){
		if(_fit[sex]) delete _fit[sex];
		_fit[sex] = obj;
	}

	//////////////////////////////////////////////////////////////////////////////
	void set_selection_level_coef();
    void set_frequency_dependend_selection();
	inline unsigned int get_SoftHardSelection(const unsigned int& soft, const unsigned int& hard)
	{return (this->*get_SoftHardSelection_func_ptr)(soft, hard);}
	unsigned int (TSelection::*get_SoftHardSelection_func_ptr)(const unsigned int& soft, const unsigned int& hard);

	inline unsigned int getSoft(const unsigned int& soft, const unsigned int& hard){return soft;} // soft
	inline unsigned int getHard(const unsigned int& soft, const unsigned int& hard){return hard;} // hard
	inline unsigned int getSoftHard(const unsigned int& soft, const unsigned int& hard){          // intermediate
		assert(_selection_level_coef>0 && _selection_level_coef<1);
		return my_round(_selection_level_coef*hard + (1-_selection_level_coef)*soft);
	}

	inline double getMeanFitness(sex_t sex){   // of a single sex
		assert(_fit[sex]);
		return _fit[sex]->getMeanFitness();
	}

	inline double getMeanFitness(){return (this->*getMeanFitness_func_ptr)();}
	double (TSelection::*getMeanFitness_func_ptr)();
	inline double getMeanFitnessFem(){return getMeanFitness(FEM);}
	inline double getMeanFitnessMal(){return getMeanFitness(MAL);}
	inline double getMeanFitnessBoth(){            // of both sexes
		assert(_fit[FEM && _fit[MAL]]);
		return (getSumFitness(MAL) + getSumFitness(FEM))/(get_nbInd(MAL) + get_nbInd(FEM));
	}

	inline double getSumFitness(sex_t sex){
		if(!_fit[sex]) return 0;
		return _fit[sex]->getSumFitness();
	}

	void   set_selection_level_coef(double d){_selection_level_coef = d;}
	double get_selection_level_coef(){return _selection_level_coef;}

	TIndividual**  get_aInd(sex_t sex)  {return _fit[sex]->_aInd;}
	double*       get_aFit(sex_t sex)  {return _fit[sex]->_aFit;}
	unsigned int  get_nbInd(sex_t sex) {return _fit[sex]->_nbInd;}
	void          remove(sex_t sex, unsigned int i){_fit[sex]->remove(i);}

	// function setting the fitness for individuals of an entire patch
	void set_fitness(TPatch*, age_idx, int* sort=NULL, int subset=0);
	void set_fitness(TPatch*, sex_t, age_idx); // only for a single sex

	void sort_fitness(sex_t SEX, int how, int subset=0);

    void    set_phenotype(TIndividual* ind, unsigned int qtraitID);
	void    set_phenotype(TIndividual* ind);
	void    set_phenotype_of_all_individuals(age_idx AGE);

	void executeBeforeEachReplicate(const unsigned int& rep);
	void executeBeforeEachGeneration(const unsigned int& gen){}
	void temporal_change(const unsigned int& gen);
    
    // female sex allocation
    TEquation*  _female_sex_allocation;
    void (TSelection::*female_sex_allocation_func_ptr)(TPatch* patch, age_idx age, int* sort, int subset);
    void female_sex_allocation_equation(TPatch* patch, age_idx age, int* sort, int subset);
    void female_sex_allocation_none(TPatch* patch, age_idx age, int* sort, int subset){};

    void female_sex_allocation_fix(TPatch* patch, age_idx age, int* sort, int subset);
    void female_sex_allocation_G(TPatch* patch, age_idx age, int* sort, int subset);
    void female_sex_allocation_Z(TPatch* patch, age_idx age, int* sort, int subset);
    void set_fem_sex_allocation(unsigned int model, double value);
    void set_fem_sex_allocation(TEquation* p){
        _female_sex_allocation=p;
        female_sex_allocation_func_ptr = p ? &TSelection::female_sex_allocation_equation : &TSelection::female_sex_allocation_none;
    }

	void      reset_Ve(bool all=true);
	void      set_ve_mean();
	void      set_ve_var();

	TMetapop*  get_popPtr()                        {return _popPtr;}
	int       get_vTraits(const int& i) const     {return _vTraits[i];}
	vector<int> get_vTraits() const               {return _vTraits;}
	TSelectionTrait* get_selTrait(const int& i)   {return _selTrait[i];}
	double    get_optima_sd()                     {return _optima_sd;}
	double    get_intensity_sd()                  {return _intensity_sd;}
	double    get_min_sd()                        {return _min_sd;}
	double    get_max_sd()                        {return _max_sd;}
	double    get_growth_rate_sd()                {return _growth_rate_sd;}
	double    get_max_growth_sd()                 {return _max_growth_sd;}
	double    get_symmetry_sd()                   {return _symmetry_sd;}
	double    get_selCoef_sd()                    {return _selCoef_sd;}
	bool      get_Ve_mean_set()                   {return _Ve_mean_set;}
	bool      get_Ve_var_set()                    {return _Ve_var_set;}

	virtual void loadFileServices ( FileServices* loader ) {}
	virtual void loadStatServices ( StatServices* loader ) {}
};

#endif
