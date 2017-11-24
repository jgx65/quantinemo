/** @file patch.h
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
//---------------------------------------------------------------------------

#ifndef patchH
#define patchH

#include <list>
#include <vector>
#include <map>
#include <time.h>
#include "types.h"
#include "indfactory.h"
#include "metapop_sh.h"
//---------------------------------------------------------------------------
//CLASS PATCH

/**Second class in the metapopulation design structure, between the Metapop and Individual classes.
 * The Patch class is an abstraction of a sub-population or patch concept (also called a deme) in a metapopulation context.
 * It contains the individual containers for the different age classes. Three main age classes are currently implemented,
 * the offspring, post-dispersal and adult classes (see the age_t enum) which are all subdivided into male and female individuals.
 * These containers are accessed using the interface defined here or through the Metapop class interface.
 * The different LCE will use these interfaces to handle the individuals. They are also responsible to change the age flag of the population (see Metapop).
 *
 * The individuals are accessed using their age index value and not the age flag value (e.g., using the Patch::get() method).
 * These indexes are defined in the age_idx enum (see type.h) and differ from the age class flag values. For instance the adults'
 * containers have the index 2 (ADLTx = 2) whereas their age flag is 4 (ADULTS = 4). This might be confusing but it saves a lot
 * of checks at runtime! It also allows to give a flag containing several class bits set instead of a unique index value when needed
 * (see the Patch::size() and Metapop::size() suite of functions).
 */

class TSelectionTrait;
class Patch
{
	unsigned int _ID;        // id of the patch (starting at 0)
	string       _IDstr;     // same as _ID (but starting at 1), but stored as a string   num2str(_ID+1)
	unsigned int _nbPatch;   // total number of patches
	TMetapop* _popPtr;

public:
	inline void         set_size(sex_t SEX, age_idx AGE, const unsigned int& s){_sizes[SEX][AGE] = s;}
	inline void         add_size(sex_t SEX, age_idx AGE, const unsigned int& s){_sizes[SEX][AGE] += s;}
	inline unsigned int get_size(sex_t SEX, age_idx AGE) {return _sizes[SEX][AGE];}
	void                addImmigrant(Patch* p, const unsigned int& s);
	void                addImmigrant(const unsigned int& p, const unsigned int& s);

protected:

	unsigned long _ID_individual;            /** counter for the individual id (starts at 0 for each simulation) */
	unsigned int _K, _Ksex[2];               /** Carrying capacities */
	unsigned int _N_ini, _N_ini_sex[2];      /** Sex specific initial population sizes */
	double _N_sample, _N_sample_sex[2];      /** Sample size of the patch, if <1 then proportion of the pop size, if >1 then number of samples */
	unsigned int _sampleID;                  /** if the patch is sampled AND CURRETLY populated it contains the  sample id (index in the vector _vFullPatch)
	 * if sampled but CURRENTLY not populated it is 55555
	 * if not sampled at all it is NaN
	 */

	double _friction, _friction_sex[2];       /** Friction of the patch (i.e. disp_final = _friction*disp_rate) */

	// selection pressure  selection[sex][trait][param]
	double*** _localSelection;         // if stabilizing selection:     _localSelection[sex][trait][param]
	//   0. pos: optima       (default: 0)
	//   1. pos: intensity    (default: 1)
	// if directional selection:
	//   0. pos: min          (min fitness,                 default: 0)
	//   1. pos: max          (max fitness,                 default: 1)
	//   2. pos: max growth   (phenotype of maximal growth, default: 1)
	//   3. pos: growth rate  (steepness of the slope,      default: 1)
	//   4. pos: symmetry     (symmetry of the slope,       default: 1)
	// if fitness landscape  (tot_size = 2*size+1)
	//   0. pos: number of data points (size)
	//   1 -> size: phenotype landscape
	//   size+1 -> 2*size: fitness landscape

	TSelectionTrait** _pSelectionType; // pointer to the selection trait neutral, stabilizing, directional, or fitness landscape

	double*   _meanVe[2];              // mean environmental Factor  _meanVe[sex][trait]
	double*   _sdVe[2];                // environmental sd computed not input   _sdVe[sex][trait]
	double*   _h2[2];                  // Ve or h2: stored input     _sdVe[sex][trait]
	unsigned int _nbLinkedTraits;

	bool _isExtinct;                    /** Extinction flag */
	vector<Individual*>** _containers;  /** Individuals containers, sex X age ([sex][age]) not used for coalescence */
	unsigned int** _sizes;              /** Containers size counters, sex X age ([sex][age]  also used for coalescence */
	vector<Individual*>** _sampled_inds;/** Individuals to sample for stats and output (_sampled_inds[sex][age]) */

public:
	//counters:
	unsigned int nbEmigrant, nbImmigrant, nbKolonisers;  // not used in coalescence!!!

//	//methods:
//	Patch() : _ID(0), _ID_individual(0), _K(0), _N_ini(0),
//			_localSelection(0), _pSelectionType(0), _nbLinkedTraits(0), _isExtinct(1), _popPtr(0)
//	{
//		init_containers();
//	}
//
//	// copy constructor (however gives just a empty initial Patch)
//	Patch(const Patch& p) : _ID(0), _ID_individual(0), _K(0), _N_ini(0),
//			_localSelection(0), _pSelectionType(0), _nbLinkedTraits(0), _isExtinct(1), _popPtr(0)
//	{
//		init_containers();
//	}
    
    Patch(TMetapop* p, unsigned int i);


	~Patch();

	Patch*        init                      (unsigned int id);
	Patch*        init_coal                 (unsigned int id);
	void          init_containers           ( );
	inline void   set_ID_individual         (unsigned int i) {_ID_individual = i;}   // used when genotypes are passed
	void          set_PopSizes_ini          (unsigned int nbfem, unsigned int nbmal);
	void          set_friction              (double nbfem, double nbmal);
	void          set_PopSizes_ini_carrying_capacity();



	///@name Setters
	///@{
	void          set_isExtinct             (bool status);
	void          set_localMeanVe		        (double *array, sex_t SEX);
	void          set_localh2Ve		          (double *array, sex_t SEX);
	void          set_localOptima		        (double *array, sex_t SEX);
	void          set_localIntensity		    (double *array, sex_t SEX);
	void          set_localSelCoefAA 		    (double *array, sex_t SEX);
	void          set_localSelCoef   		    (double *array, sex_t SEX);
	void          set_localGrowthRate		    (double *array, sex_t SEX);
	void          set_localMaxGrowth		    (double *array, sex_t SEX);
	void          set_localSymmetry 		    (double *array, sex_t SEX);
	void          set_localMin 		          (double *array, sex_t SEX);
	void          set_localMax      		    (double *array, sex_t SEX);
	void          set_localOptima		        (unsigned int t, double array, sex_t SEX);
	void          set_localIntensity		    (unsigned int t, double array, sex_t SEX);
	void          set_localSelCoefAA        (unsigned int t, double array, sex_t SEX);
	void          set_localSelCoef          (unsigned int t, double array, sex_t SEX);
	void          set_localGrowthRate		    (unsigned int t, double array, sex_t SEX);
	void          set_localMaxGrowth		    (unsigned int t, double array, sex_t SEX);
	void          set_localSymmetry 		    (unsigned int t, double array, sex_t SEX);
	void          set_localMin 		          (unsigned int t, double array, sex_t SEX);
	void          set_localMax      		    (unsigned int t, double array, sex_t SEX);
	void          set_fitnessLandscape_phenotype(unsigned int trait, double *array, unsigned int size, sex_t SEX);
	void          set_fitnessLandscape_fitness  (unsigned int trait, double *array, unsigned int size, sex_t SEX);
	bool          set_Ve                    (const int& model, const int& trait, sex_t sex);
	void          set_h2                    (const double& v, const int& i){   // input of the h2
		if(_h2[0]) _h2[0][i] = v;   // check if there are males present
		_h2[1][i] = v;              // females are always present
	}
	void   set_nbPatch               (const unsigned int& n){_nbPatch=n;}
	void   set_pMetapop              (TMetapop* p)           {_popPtr=p;}

	template <typename T>
	void set_localParameter(T value, void (Patch::*pt2Func)(T)){
		(this->*pt2Func)(value);
	}

	template <typename T>
	void set_localParameter(T value, sex_t sex, void (Patch::*pt2Func)(T, sex_t)){
		(this->*pt2Func)(value, sex);
	}

	template <typename T>
	void set_localParameter(T fem, T mal, T gen, void (Patch::*pt2Func)(T, T, T)){
		(this->*pt2Func)(fem, mal, gen);
	}

	template <typename T>
	void set_localParameter_ofTrait(unsigned int t, T value, void (Patch::*pt2Func)(unsigned int, T)){
		(this->*pt2Func)(t, value);
	}

	template <typename T>
	void set_localParameter_ofTrait(unsigned int t, T value, sex_t sex, void (Patch::*pt2Func)(unsigned int, T, sex_t)){
		(this->*pt2Func)(t, value, sex);
	}

	template <typename T>
	void set_localParameter_ofTrait(unsigned int t, T fem, T mal, T gen, void (Patch::*pt2Func)(unsigned int, T, T, T)){
		(this->*pt2Func)(t, fem, mal, gen);
	}

	void set_localParameter(double* array, sex_t sex, void (Patch::*pt2Func)(double*, sex_t));
	void set_localParameter_matrix(unsigned int t, double* array, unsigned int size, sex_t sex,
			void (Patch::*pt2Func)(unsigned int,  double*, unsigned int, sex_t));
	void set_localParameter_matrix_ofTrait(unsigned int t, double* array, unsigned int size, sex_t sex,
			void (Patch::*pt2Func)(unsigned int, double*, unsigned int, sex_t));
	template<typename T>
	void setSexSpecificPatchParam_sexRatio(T (Patch::*getGeneral)(),
			void (Patch::*setSexSpecific)(T, sex_t), double sexRatio){
		T val = (this->*getGeneral)();
		T malVal = (T)(val * sexRatio);
		(this->*setSexSpecific)(malVal, MAL);
		(this->*setSexSpecific)(val-malVal, FEM);
	}
	template<typename T>
	void setGeneralPatchParam_fem(T (Patch::*getSexSpecific)(sex_t), void (Patch::*setGeneral)(T)){
		(this->*setGeneral)((this->*getSexSpecific)(FEM));
	}
	template<typename T>
	void setGeneralPatchParam_sum(T (Patch::*getSexSpecific)(sex_t), void (Patch::*setGeneral)(T)){
		(this->*setGeneral)((this->*getSexSpecific)(FEM) + (this->*getSexSpecific)(MAL));
	}

	unsigned int  set_sampledInds          (sex_t SEX, age_idx AGE); 	// for the current generation
	bool          set_sampledInds          (sex_t SEX, age_t AGE);   	// for the current generation
	unsigned int  get_sampleID(){return _sampleID;}
	void          set_sampleID(unsigned int i){_sampleID=i;}

	void          reset_ID_individual       ()               {_ID_individual=0;}


	///@}
	///@name Getters
	///@{
	inline unsigned int  get_ID             ()               {return _ID;}
	string        get_next_IDofIndividual   ()               {return toStr(_ID_individual++)+"_"+_IDstr;}  // "234_1": individual 234 of patch 1
	unsigned int  get_K                     ()               {return _K;}
	unsigned int  get_K                (sex_t SEX)           {return _Ksex[SEX];}
	unsigned int  get_KFem                  ()               {return _Ksex[FEM];}
	unsigned int  get_KMal                  ()               {return _Ksex[MAL];}
	unsigned int  get_N_ini                 ()               {return _N_ini;}
	unsigned int  get_N_ini            (sex_t SEX)           {return _N_ini_sex[SEX];}
	double        get_N_sample              ()               {return _N_sample;}
	double        get_N_sample         (sex_t SEX)           {return _N_sample_sex[SEX];}
	unsigned int  get_curN_sample      (sex_t SEX, age_idx AGE); 	// for the current generation

	double        get_density   (const sex_t& SEX, const age_t& AGE){return (double)size(SEX, AGE)/_Ksex[SEX];}
	double        get_density   (const age_t& AGE)                  {return (double)size(AGE)/_K;}
	double        get_density   (const sex_t& SEX, const age_idx& AGE){return (double)size(SEX, AGE)/_Ksex[SEX];}
	double        get_density   (const age_idx& AGE)                {return (double)size(AGE)/_K;}

	unsigned int**        &get_sizes        ()               {return _sizes;}
	vector <Individual*>**&get_containers   ()               {return _containers;}
	vector <Individual*>*& get_containers   (const sex_t& SEX) {return _containers[SEX];}
	vector <Individual*>&  get_containers   (const sex_t& SEX, const age_idx& AGE){return _containers[SEX][AGE];}
	vector <Individual*>&  get_all_inds     (const sex_t& SEX, const age_idx& AGE){return _containers[SEX][AGE];}

	// current samples
	vector<Individual*>& get_sampled_inds(const sex_t& SEX, const age_idx& AGE){return _sampled_inds[SEX][AGE];}
	unsigned int  sampleSize(const sex_t SEX, const age_idx& AGE);
	unsigned int  sampleSize(const sex_t SEX, const age_t& AGE);
	unsigned int  sampleSize(const age_idx& AGE){return sampleSize(FEM, AGE) + sampleSize(MAL,AGE);}
	unsigned int  sampleSize(const age_t& AGE){return sampleSize(FEM, AGE) + sampleSize(MAL,AGE);}

	double        get_friction              ( )              {return _friction;}
	double        get_friction              (sex_t SEX)      {return _friction_sex[SEX];}

	bool          get_isExtinct             ()               {return _isExtinct;}
	double**      get_localSelection        (sex_t SEX)      {return _localSelection[SEX];}
	double*       get_localSelection        (sex_t s, int t) {return _localSelection[s][t];}
	double        get_localSelection        (sex_t s, int t, int p) {return _localSelection[s][t][p];}
	double*       get_meanVe                (sex_t s)        {return _meanVe[s];}
	double        get_meanVe                (sex_t s,int pos){return _meanVe[s][pos];}
	double*       get_sdVe                  (sex_t s)        {return _sdVe[s];}
	double        get_sdVe                  (sex_t s,int pos){return _sdVe[s][pos];}
	double*       get_h2                    (sex_t s)        {return _h2[s];}
	double        get_h2                    (sex_t s,int pos){return _h2[s][pos];}
	int           get_nbLinkedTraits        ()               {return _nbLinkedTraits;}
	bool          isEmpty                   ()               {return (size(ALL) == 0);}
	unsigned int  getAdultsNumber           ()               {return size(ADLTx);}
	unsigned int  getOffspringNumber        ()               {return size(OFFSx);}
	unsigned int  getIndividualNumber       ()               {return size(ALL);}
	double        getPatchOffsprgSatLevel   ()               {return (double)size(OFFSx)/_K;}
	double        getPatchAdultsSatLevel    ()               {return (double)size(ADLTx)/_K;}
    TSelectionTrait** get_pSelectionType    ()               {return _pSelectionType;}

	unsigned int  size       (sex_t SEX, age_t AGE);
	unsigned int  size       (age_t AGE=ALL);
	unsigned int  size       (sex_t SEX, age_idx AGE);
	unsigned int  size       (age_idx AGE);
	Individual*   get        (sex_t SEX, age_idx AGE, unsigned int at);
	void          set        (sex_t SEX, age_idx AGE, unsigned int at, Individual* ind);
	void          add        (const sex_t& SEX, const age_idx& AGE, Individual* ind);
	void          assign     (sex_t SEX, age_idx AGE, unsigned int n);
	void          remove     (sex_t SEX, age_idx AGE, unsigned int at);
	void          remove     (sex_t SEX, age_idx AGE, Individual* ind);
	void          recycle    (sex_t SEX, age_idx AGE, unsigned int at);
	void          recycle    (sex_t SEX, age_idx AGE, Individual* ind);
	void          clear      ();
	void          clear      (age_idx AGE);
	void          clear      (sex_t SEX, age_idx AGE);
    bool          individual_container_ok();

	void 	        (Patch::*func_ptr_swap)(sex_t SEX, age_idx from, age_idx to);
	void          swap       (sex_t SEX, age_idx from, age_idx to);
	void          swap_ind  (sex_t SEX, age_idx from, age_idx to);
	void          swap_coal  (sex_t SEX, age_idx from, age_idx to);
	void          swap       (age_idx from, age_idx to);

	void          move       (sex_t SEX, age_idx from, age_idx to, unsigned int at);
	void          move       (sex_t SEX, age_idx from, age_idx to);
	void          move_ind   (sex_t SEX, age_idx from, age_idx to);
	void          move_coal  (sex_t SEX, age_idx from, age_idx to);
	void 	        (Patch::*func_ptr_move)(sex_t SEX, age_idx from, age_idx to);
	void          move       (age_idx from, age_idx to);

	void          survive_randomly_inds_relative(sex_t SEX, age_idx AGE, double ratio);
	void          survive_randomly_inds_absolute(sex_t SEX, age_idx AGE, unsigned int size);
	void          remove_randomly_inds_relative(sex_t SEX, age_idx AGE, double ratio);
	void          remove_randomly_inds_absolute(sex_t SEX, age_idx AGE, unsigned int size);

	void          aging();

	void          set_func_pointer();
	void          set_func_pointer_coal();

	void          (Patch::*func_ptr_flush)(sex_t SEX, age_idx AGE);
	void          flush      ();
	void          flush      (age_t AGE);
	void          flush      (age_idx AGE);
	void          flush      (sex_t SEX, age_idx AGE);
	void          flush_ind  (sex_t SEX, age_idx AGE);
	void          flush_coal (sex_t SEX, age_idx AGE){_sizes[SEX][AGE]=0;}

	void reset_counters();

	void set_LinkedTraits(const unsigned int& size, TSelectionTrait** type, bool bothSexes);
	void reset_LinkedTraits();

	/** return an array with all phenotypes of the specific trait, sex and age */
	double* getPhenotypes(sex_t SEX, age_idx AGE, const int& trait);

	/** pos dize regulation based on their fitness */
	void regulate_selection_fitness(const unsigned int& K, TSelection* pSel, sex_t SEX, age_idx AGE);
	void regulation_selection_draw_survivors(TSelection* pSel, const unsigned int& K, sex_t SEX, age_idx AGE);
	void regulation_selection_draw_loosers(TSelection* pSel, const unsigned int& K, sex_t SEX, age_idx AGE);
	void regulation_selection_hard(age_idx AGE);

	/**Fills the patch containers corresponding to the age flags passed, for both sexes.*/
	unsigned int setNewGeneration(age_t AGE);
	unsigned int setNewGeneration(age_idx AGE);
	unsigned int setNewGeneration_coal(age_t AGE);

	void set_N_sample(double Nfem, double Nmal, double N){
		_N_sample_sex[FEM] = Nfem; _N_sample_sex[MAL] = Nmal; _N_sample = N;
	}
	void set_N_sample(double size, sex_t sex){_N_sample_sex[sex] = size;}
	void set_N_sample(double size){_N_sample = size;}

	void set_K(unsigned int Kfem, unsigned int Kmal, unsigned int K){
		_Ksex[FEM] = Kfem; _Ksex[MAL] = Kmal; _K = K;
	}
	void set_K(unsigned int K, sex_t sex){_Ksex[sex] = K;}
	void set_K(unsigned int K){_K = K;}

	void set_N_ini(unsigned int Nfem, unsigned int Nmal, unsigned int N){
		_N_ini_sex[FEM] = Nfem; _N_ini_sex[MAL] = Nmal; _N_ini = N;
	}
	void set_N_ini(unsigned int N, sex_t sex){_N_ini_sex[sex] = N;}
	void set_N_ini(unsigned int N){_N_ini = N;}

	void set_friction(double Kfem, double Kmal, double K){
		_friction_sex[FEM] = Kfem; _friction_sex[MAL] = Kmal; _friction = K;
	}
	void set_friction(double K, sex_t sex){_friction_sex[sex] = K;}
	void set_friction(double K){_friction = K;}
};
#endif
