/** @file lce_breed.h
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
 */

#ifndef lce_breedH
#define lce_breedH


#include "lifecycleevent.h"
#include "tequation.h"


class TSelection;

// Class LCE_Breed
//
/**Class for the breeding (and mating) life cycle events.
   This class registers the \a mating_system, \a mating_proportion, and \a mean_fecundity
   parameters to the \c ParamSet. Sets the mating function variable to the right mating function.

   Implementation of the basic breeding and mating procedures, does not link to any trait.
   Individuals mate according to the mating system chosen. The number of offspring by female is driven
   from a Poisson distribution with mean equal to the \a mean_fecundity parameter value. The mated adults
   are not removed from the population. The offspring containers are filled with the new generation. Note that
   they are first emptied if still containing offspring individuals.

   The population's age is set to \c ALL. The adults mating and realized fecundity counters
   are updated.
 **/
class LCE_Breed : public LCE  {

public:
	LCE_Breed(int rank = my_NAN);

	virtual ~LCE_Breed   ( ) {
		if(_aMatingPairs[MAL]) delete[] _aMatingPairs[MAL];
		if(_aMatingPairs[FEM]) delete[] _aMatingPairs[FEM];
	}

protected:
	void (LCE_Breed::* breed)();

	void create_mating_pairs(TPatch* cur_patch);
	void createOffspring(TPatch* cur_patch, unsigned int nbDaugthers, unsigned int nbSons);

	TIndividual**  _aMatingPairs[2];      // 0: the male, 1: the female (used for the monogamy mating system
	unsigned int  _aMatingPairs_size;    // the size of aMatingPairs

	int           _mating_system;         // 0: random mating (hermaphrodite, selfing by chance allowed (1/N))(default)
	// 1: selfing       (hermaphrodite, selfing rate depends on parameter mating_proportion)
	// 2: cloning       (hermaphrodite, sexual reproduction depends on parameter mating_proportion)
	// 3: promiscuity / random mating (two sexes)
	// 4: polygyny
	// 5: monogamy
	// 6: no mating/reproduction occurs (used to just compute stats from ini genotype file)
	double        _threshold;
	unsigned int  _fitness_factor_zero_isLethal; // 0: no; 1: yes

	unsigned int  _mating_males;
	double        _mating_proportion;
	double        _mean_fecundity;
	double        _growth_rate;       // a growth rate of 0 is stable!!!!
	double        _sex_ratio;         // saved as males/(males+females); input as males/females
    
    TEquation*    _fem_sex_allocation;
    unsigned int  _fem_sex_allocation_quanti_trait;  // fealtive quanti trait index (if not used my_NAN)
    unsigned int  _fem_sex_allocation_all_trait;     // absolute trait index (if not used my_NAN)
    double        _fem_sex_allocation_prop;          // proportion (if not used my_NAN)

	unsigned int  _nbOffspring_model;
    // 0: carrying capacity:	          nbOffs = K (default)
	// 1: keep number                     nbOffs = nbAdult
	// 3: fecundity:                      nbOffs = ranPossion(nbFemales*meanFecundity)
	// 3: fecundity simple:               nbOffs = nbFemales*meanFecundity
	// 3: fecundity limited:              nbOffs = ranPossion(nbFemales*meanFecundity), but max K
	// 3: fecundity simple & limited:     nbOffs = nbFemales*meanFecundity, but max K
	// 4: logistic regulated: 	          nbOffs = logisticGrowth(nbAdults, carrying capacity)
	// 5: logistic regulated (stochastic):nbOffs = ranPoisson(logisticGrowth(nbAdults, K))

	// temporaire varaibles
	double _nbIndividuals[2]; // current number of females and males (real number due to female sex allocation)
	// temporare array for each sex when selection occurs (male: 0; female: 1)
	int           _sort[2]; // how to use the fitness: 0: get random fittest;
	//                         1: get random fittest of a subset of fittest;
	//                         2: get random less fittest;
	//                         3: get random less fittest of a subset of less fittest;

	TIndividual* (LCE_Breed::* getMother_func_ptr)(TPatch*, unsigned int&, sex_t sex);///< A pointer to a mating function for females
	TIndividual* (LCE_Breed::* getFather_func_ptr)(TPatch*, unsigned int&, sex_t sex);///< A pointer to a mating function for males


	void preMating(TPatch* cur_patch);

	/** function to set the number of females and males and checks if the conditions are met for mating */
    bool (LCE_Breed::* isMatingPossible_func_ptr)(TPatch*);
    bool isMatingPossible_1_sex(TPatch* cur_patch);
    bool isMatingPossible_2_sex(TPatch* cur_patch);
    bool isMatingPossible_sex_allocation_fix(TPatch* cur_patch);
    bool isMatingPossible_sex_allocation_G(TPatch* cur_patch);
    bool isMatingPossible_sex_allocation_Z(TPatch* cur_patch);
    bool isMatingPossible_sex_allocation_equation(TPatch* cur_patch);
    
	/** function to compute the number of sons/daugthers depending on the total number of childs
	 * nbSons and nbDaughters are changed
	 * nbBaby is the total number of new offspring
	 * nbMAL and nbFEM are the number of females and males currently in the patch
	 *     (nbMAL and nbFEM are used to compute the current sex ratio)
	 */
	void (LCE_Breed::* setSexRatio_func_ptr)(const unsigned int&, unsigned int&, unsigned int&, const unsigned int&, const unsigned int&);
    void setSexRatio_Selfing(const unsigned int& nbBaby, unsigned int& nbSons, unsigned int& nbDaugthers, const unsigned int& nbMAL, const unsigned int& nbFEM);
        void setSexRatio_NoSelfing(const unsigned int& nbBaby, unsigned int& nbSons, unsigned int& nbDaugthers, const unsigned int& nbMAL, const unsigned int& nbFEM);
        void setSexRatio_KeepSexRatio(const unsigned int& nbBaby, unsigned int& nbSons, unsigned int& nbDaugthers, const unsigned int& nbMAL, const unsigned int& nbFEM);

	sex_t (LCE_Breed::* getRandomSex_func_ptr)(const unsigned int&, const unsigned int&);
    sex_t getRandomSex_Selfing(const unsigned int& nbMAL, const unsigned int& nbFEM);
    sex_t getRandomSex_NoSelfing(const unsigned int& nbMAL, const unsigned int& nbFEM);
    sex_t getRandomSex_KeepSexRatio(const unsigned int& nbMAL, const unsigned int& nbFEM);

	sex_t _maleSex;    // normaly MAL, but when only one sex is used this is FEM

	TSelection* _pSelection; // do not delete this object here (it belongs to the metapop)

public:

    // function to compute the number of daugthers and sons (returns the total number of offspring)
    unsigned int (LCE_Breed::* setNbOffspring_func_ptr) (double, double, unsigned int);
    unsigned int setNbOffspring_KeepNb(double nbMAL, double nbFEM, unsigned int K){
        return nbMAL+nbFEM;
    }
    unsigned int setNbOffspring_CarryCapacity(double nbMAL, double nbFEM, unsigned int K){
        return K;
    }
    unsigned int setNbOffspring_Logistic(double nbMAL, double nbFEM, unsigned int K){
        // return logisticGrowth(_growth_rate, K, nbMAL+nbFEM);
        return my_round(beverton_hold(_growth_rate, K, nbMAL+nbFEM));
    }
    unsigned int setNbOffspring_RandLogistic(double nbMAL, double nbFEM, unsigned int K){
        // return SimRunner::r.Poisson(logisticGrowth(_growth_rate, K, nbMAL+nbFEM));
        return get_pop_ptr()->rand().Poisson(beverton_hold(_growth_rate, K, nbMAL+nbFEM));
    }
    unsigned int setNbOffspring_Fecundity(double nbMAL, double nbFEM, unsigned int K){
        return get_pop_ptr()->rand().Poisson(nbFEM*_mean_fecundity);
    }
    unsigned int setNbOffspring_FecunditySimple(double nbMAL, double nbFEM, unsigned int K){
        return my_round(nbFEM*_mean_fecundity);
    }
    unsigned int setNbOffspring_FecundityBinomial(double nbMAL, double nbFEM, unsigned int K){
        double valDbl = nbFEM*_mean_fecundity;
        unsigned int valInt = (unsigned int)valDbl;
        return valInt + get_pop_ptr()->rand().Binomial(valDbl - valInt, 1);
    }
    unsigned int setNbOffspring_FecundityLimited(double nbMAL, double nbFEM, unsigned int K){
        unsigned int size = setNbOffspring_Fecundity(nbMAL, nbFEM, K);
        return (size>K ? K : size);
    }
    unsigned int setNbOffspring_FecunditySimpleLimited(double nbMAL, double nbFEM, unsigned int K){
        unsigned int size = setNbOffspring_FecunditySimple(nbMAL, nbFEM, K);
        return (size>K ? K : size);
    }
    unsigned int setNbOffspring_FecundityBinomialLimited(double nbMAL, double nbFEM, unsigned int K){
        unsigned int size = setNbOffspring_FecundityBinomial(nbMAL, nbFEM, K);
        return (size>K ? K : size);
    }

	/**Link to the mating function, used to get the father from the mother in a Patch
    @param thePatch Patch instance of the current where the father has to be fetched
    @param mother index of the mother in the current Patch, used in the \a polyginy an \a monoginy mating systems
    @return the pointer to the father chosen following the mating scheme chosen
	 **/
	virtual TIndividual* getMotherPtr (TPatch* thePatch, unsigned int& index)
	{ return (this->*getMother_func_ptr)(thePatch, index, FEM);  }
	virtual TIndividual* getFatherPtr (TPatch* thePatch, unsigned int& index)
	{ return (this->*getFather_func_ptr)(thePatch, index, _maleSex);  }

	///@name getter
	///@{

	string getMatingSystem_str () {
		switch(_mating_system){
		case 0: return "random mating (hermaphrodite)";
		case 1: return "selfing (hermaphrodite)";
		case 2: return "cloning (hermaphrodite)";
		case 3: return "random mating (promiscuity)";
		case 4: return "polygyny";
		case 5: return "monogamy";
		case 6: return "no mating/reproduction";
		}
		return "";
	}

	///@}

	///@name Mating functions
	///@{
	// get the individual of the index ///////////////////////////////////////////
	TIndividual* Index_MatingFunc                        (TPatch*, unsigned int&, sex_t); // get the individual of the index
	TIndividual* NULL_pointer                            (TPatch*, unsigned int&, sex_t); // used for cloning

	// random mating /////////////////////////////////////////////////////////////
	TIndividual* Random_MatingFunc                       (TPatch*, unsigned int&, sex_t);
	TIndividual* Random_Index_MatingFunc                 (TPatch*, unsigned int&, sex_t); // return the index of the selected individual

	TIndividual* Random_S_MatingFunc                     (TPatch*, unsigned int&, sex_t);
	TIndividual* Random_Index_S_MatingFunc               (TPatch*, unsigned int&, sex_t);

	// full polygyny /////////////////////////////////////////////////////////////
	TIndividual* fullPolygyny_oneMale_MatingFunc         (TPatch*, unsigned int&, sex_t);
	TIndividual* fullPolygyny_manyMales_MatingFunc       (TPatch*, unsigned int&, sex_t);

	TIndividual* fullPolygyny_oneMale_S_MatingFunc       (TPatch*, unsigned int&, sex_t); // the most fittest (not random)
	TIndividual* fullPolygyny_manyMales_S_MatingFunc     (TPatch*, unsigned int&, sex_t); // get the x most fittest (not random)

	TIndividual* fullPolygyny_oneMale_S_MatingFunc2      (TPatch*, unsigned int&, sex_t); // get the randomly chosen most fittest
	TIndividual* fullPolygyny_manyMales_S_MatingFunc2    (TPatch*, unsigned int&, sex_t); // get the randomly chosen x most fittest

	// partial polygyny //////////////////////////////////////////////////////////
	TIndividual* partialPolygyny_oneMale_MatingFunc      (TPatch*, unsigned int&, sex_t);
	TIndividual* partialPolygyny_manyMales_MatingFunc    (TPatch*, unsigned int&, sex_t);

	TIndividual* partialPolygyny_oneMale_S_MatingFunc    (TPatch*, unsigned int&, sex_t); // the most fittest (not random)
	TIndividual* partialPolygyny_manyMales_S_MatingFunc  (TPatch*, unsigned int&, sex_t); // get the x most fittest (not random)

	TIndividual* partialPolygyny_oneMale_S_MatingFunc2   (TPatch*, unsigned int&, sex_t); // get the randomly chosen most fittest
	TIndividual* partialPolygyny_manyMales_S_MatingFunc2 (TPatch*, unsigned int&, sex_t); // get the randomly chosen x most fittest

	// monogamy //////////////////////////////////////////////////////////////////
	TIndividual* Monogyny_MatingFunc                     (TPatch*, unsigned int&, sex_t); // mating pairs have to be first fixed
	TIndividual* Monogyny_S_MatingFunc                   (TPatch*, unsigned int&, sex_t); // mating pairs have to be first fixed

	// one sex ///////////////////////////////////////////////////////////////////
	TIndividual* oneSex_notSameIndex_MatingFunc          (TPatch*, unsigned int&, sex_t); // get random individual, but not the same index
	TIndividual* partialSelfing_MatingFunc               (TPatch*, unsigned int&, sex_t); // partial sefling, the rest is random mating

	TIndividual* oneSex_notSameIndex_S_MatingFunc        (TPatch*, unsigned int&, sex_t); // get random individual, but not the same index
	TIndividual* partialSelfing_S_MatingFunc             (TPatch*, unsigned int&, sex_t); // partial sefling, the rest is random mating

	TIndividual* partialCloning_MatingFunc               (TPatch*, unsigned int&, sex_t); // partial sefling, the rest is random mating
	TIndividual* partialCloning_S_MatingFunc             (TPatch*, unsigned int&, sex_t); // partial sefling, the rest is random mating

	///@}

	///@name Implementations
	///@{
	bool init(TMetapop* popPtr);
    void set_nbOffspring_model
    ();
	age_t removeAgeClass ( ) {return 0;}
	age_t addAgeClass ( ) {return OFFSPRG;}
	age_t requiredAgeClass () {return ADULTS;}
	///@}

public:
	/** selection acts at the offspring stage:
	 * 1. a number of offspring is generated depending on the fecundity of the females (parameter "mean_fecundity")
	 * 2. the number of surviving offspring is determined (parameter "mating_nb_offspring_model")
	 * 3. downregulation of his number is made depending on the fitness of the offspring: the fittest have a higher chance to survive
	 */
	void breed_offspring_soft2metapop( );
	void breed_offspring_soft2hard   ( );
	void breed_offspring_soft2hard2   ( );
	void breed_offspring_metapop2hard( );
	void breed_none( ){return;}                   // the life cycle is jumped (used to just compute stats)

	/** selection acts at the adults age: the higher the fitness is the more offspring an individual has
	 * 1. the nubmer of total offspring is determined (parameter "mating_nb_offspring_model")
	 * 2. for each offspring parents are randomly drawn depending on their fitness: the higher the parents fitness is the
	 *    the more offspring they get
	 */
	void breed_neutral               ( );    // no selection acts
	void breed_soft2metapop          ( );
	void breed_soft2hard             ( );
	void breed_metapop2hard          ( );

	unsigned int get_fitness_factor_zero_isLethal(){return _fitness_factor_zero_isLethal;}
	void remove_inds_zero_fitnessFactor(age_idx AGE);
	void reset_sex_after_phenotype(age_idx AGE);

	///@name Implementations
	///@{

	virtual void  execute ();

	virtual LCE_Breed* clone ( ) {return new LCE_Breed();}

	virtual void loadFileServices ( FileServices* loader ) {}
	virtual void loadStatServices ( StatServices* loader ) {}

	virtual void executeBeforeEachReplicate(const unsigned int& rep);
	virtual void executeBeforeEachGeneration(const unsigned int& gen){}
    
    void temporal_change(const unsigned int& gen);
	///@}

};


class LCE_BreedCoalescence: public LCE_Breed {
public:
	LCE_BreedCoalescence(int rank = my_NAN) : LCE_Breed(rank){};
	~LCE_BreedCoalescence(){}
	void execute();
	bool init(TMetapop* popPtr);
	bool breed_coal();
};

#endif //LCEBREED_H

