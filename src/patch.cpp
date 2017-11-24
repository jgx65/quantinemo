 /** @file patch.cpp
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

#include "patch.h"
#include "tmetapop.h"
//#include "random.h"
//#include "functions.h"
#include "tselectiontrait.h"
#include "ttquanti.h"
#include "tselection.h"
#include "stathandler.cpp"  // used by the icpc compiler...

using namespace std;

// ----------------------------------------------------------------------------------------
// reset_LinkedTraitsParams
// ----------------------------------------------------------------------------------------
/** creates the space for parameters linked to the quantitative traits */
void
Patch::reset_LinkedTraits()
{
	if(_pSelectionType) {delete[] _pSelectionType; _pSelectionType = NULL;}
	ARRAY::delete_3D(_localSelection, 2, _nbLinkedTraits);
	_nbLinkedTraits = 0;

	for(int i=0; i<2; ++i){      // for each sex
		if(_meanVe[i]){delete[] _meanVe[i]; _meanVe[i]= NULL;}
		if(_sdVe[i])  {delete[] _sdVe[i];   _sdVe[i]= NULL;}
		if(_h2[i])    {delete[] _h2[i];     _h2[i]= NULL;}
	}
}

// ----------------------------------------------------------------------------------------
// set_LinkedTraits
// ----------------------------------------------------------------------------------------
/** creates the space for parameters linked to the quantitative traits */
void
Patch::set_LinkedTraits(const unsigned int& nbTraits, TSelectionTrait** type, bool bothSexes)
{
	reset_LinkedTraits();

	_nbLinkedTraits = nbTraits;
	unsigned int size;
	ARRAY::create_2D(_localSelection, 2, _nbLinkedTraits);
	_pSelectionType = new TSelectionTrait*[_nbLinkedTraits];
	for(unsigned int i=0; i<_nbLinkedTraits; ++i){
		_pSelectionType[i] = type[i];
		size = type[i]->get_nb_selection_params();
		if(size){
			_localSelection[FEM][i] = new double[type[i]->get_nb_selection_params()];
			if(bothSexes) _localSelection[MAL][i] = new double[type[i]->get_nb_selection_params()];
			else          _localSelection[MAL][i] = NULL;
		}
		else{
			_localSelection[FEM][i] = NULL;
			_localSelection[MAL][i] = NULL;
		}
	}

	if(bothSexes){
		_meanVe[MAL] = new double [_nbLinkedTraits];
		_sdVe[MAL]   = new double [_nbLinkedTraits];
		_h2[MAL]     = new double [_nbLinkedTraits];
	}
	_meanVe[FEM] = new double [_nbLinkedTraits];
	_sdVe[FEM]   = new double [_nbLinkedTraits];
	_h2[FEM]     = new double [_nbLinkedTraits];
}

// ----------------------------------------------------------------------------------------
// set_Ve
// ----------------------------------------------------------------------------------------
/** transfers the input _h2 to _sdVe depending on the selected model
 * computes the Ve for the specified trait and both sexes
 * returns true if Ve is larger than zero for males and/or females
 */
bool
Patch::set_Ve(const int& model, const int& trait, sex_t sex)
{
	// for each quantitative trait
	double h2 = _h2[sex][trait];

	// is Ve directly set?
	if(model==0){
		if(h2<0) error("The environmental variance has to be positive (%f)!\n", h2);
		_sdVe[sex][trait] = sqrt(h2);    // input is the variance, here we need the sd!!!
	}

	// if heritability is 1 => Ve=0
	else if(h2==1){ // fully heritable
		_sdVe[sex][trait] = 0;
	}

	// Ve is set by the heritability
	else{
		if(h2<=0 || h2>1) error("The heritability is out of range %f!\n", h2); // check the range of h2

		// check if the pop is populated in case where Ve is only set at the beginning of the sim
		if((model==1 || model==3) && !size(ADLTx)){
			error("The environmental variance cannot be computed since all patches have to be populated at the start of a simulation (true for quanti_environmental_model set to 1 or 3)!\n");
		}

		// narrow or broad-sense heritability?
		double var, mean;
		TTQuantiSH* pStats = _pSelectionType[trait]->get_pQuantiProto()->_stats;
		if(model==1 || model==2){ // narrow-sense heritability
			unsigned int nb_locus  = _pSelectionType[trait]->get_pQuantiProto()->get_nb_locus();
			map<unsigned char, double>* freqs = new map<unsigned char, double>[nb_locus];
			map<unsigned char, double>* freqs2 = new map<unsigned char, double>[nb_locus];  // just a placeholder for the global allele frequency
			pStats->set_alleleFreq_ofPatch_allInds(this, ADLTx, freqs, freqs2);
			pStats->get_Va_ofPatch(this, ADLTx, mean, var, freqs);  // compute Va
			if(var == my_NAN){      // if it could not be computed -> use the method for random mating to compute Va
				warning("Ve of patch %i could not be computed since Va could not be estimated -> using Va computed for random mating to set Ve (see manual parameter 'quanti_va_model')!\n",_ID+1);
				pStats->get_Va_ofPatch_random_mating(this, ADLTx, mean, var, freqs);
			}
			delete[] freqs;
		}
		else{   // broad-sense heritability
			pStats->setMeanAndVar_Vg_ofPatch_allInds(this, ADLTx, mean, var); // compute Vg
		}

		_sdVe[sex][trait] = var * (1-h2)/h2;    // set the Ve
	}

	return (_sdVe[sex][trait] != 0);  // return true if sdVe is not zero
}


// ----------------------------------------------------------------------------------------
// set_PopSizes_ini_carrying_capacity
// ----------------------------------------------------------------------------------------
void
Patch::set_PopSizes_ini_carrying_capacity()
{
	_N_ini_sex[FEM] = _Ksex[FEM];
	_N_ini_sex[MAL] = _Ksex[MAL];
	_N_ini    = _K;
	//  _isExtinct = (_N_ini==0);
}

// ----------------------------------------------------------------------------------------
// get_curN_sample
// ----------------------------------------------------------------------------------------
/** returns the number of individuals to sample in this population for the current generation */
unsigned int
Patch::get_curN_sample(sex_t SEX, age_idx AGE)
{
	unsigned int N = size(SEX, AGE);
	if(_N_sample_sex[SEX]==my_NAN || _N_sample_sex[SEX]>=N) return N;      // not specified or exceeds the number: all inds sampled
	if(_N_sample_sex[SEX] >= 1)  return (unsigned int) _N_sample_sex[SEX]; // absolute sampling
	return (unsigned int) my_round(_N_sample_sex[SEX]*N);                  // relative sampling
}

// ----------------------------------------------------------------------------------------
// set_sampledInds
// ----------------------------------------------------------------------------------------
/** at a given generation the individuals to sample are selected randomly
 * AGE: contains the age classes to sample now
 * returns yes if sampling revealed some individuals and no if no inds are sampled
 */
bool
Patch::set_sampledInds(sex_t SEX, age_t AGE)
{
	// perform the sampling for the given age classes
	unsigned int sampleSize=0, a, mask;
	for(a = 0, mask=1; a < NB_AGE_CLASSES; ++a, mask <<= 1) {
		if(mask & AGE) sampleSize += set_sampledInds(SEX, static_cast<age_idx>(a));  // for each sampled age class
	}

	return (sampleSize!=0);
}

// ----------------------------------------------------------------------------------------
// set_sampledInds
// ----------------------------------------------------------------------------------------
/** at a given generation the individuals to sample are selected randomly
 * the number of individuals to sample is returned
 */
unsigned int
Patch::set_sampledInds(sex_t SEX, age_idx AGE)
{
	vector<TIndividual*>& vSample = _sampled_inds[SEX][AGE];
	//assert(vSample.empty());
	vSample.clear();
	unsigned int curSize = size(SEX, AGE);                         // get the current patch size
	if(!curSize) return 0;                                         // if patch is empty: nothing to sample

	unsigned int curSample = get_curN_sample(SEX, AGE);            // get the number of inds to sample
	vSample.insert(vSample.end(), _containers[SEX][AGE].begin(), _containers[SEX][AGE].begin()+curSize);
	if(curSample<curSize) {
		vSample = _popPtr->rand().sample(vSample, curSample); // sample without replacement curSample individuals
	}
	return curSample;
}

// ----------------------------------------------------------------------------------------
// sampleSize
// ----------------------------------------------------------------------------------------
/** get the current sample size
 */
unsigned int
Patch::sampleSize(const sex_t SEX, const age_t& AGE)
{
	unsigned int a, mask, size = 0;
	for(a = 0, mask=1; a < NB_AGE_CLASSES; ++a, mask <<= 1) {
		if(mask & AGE) size += sampleSize(SEX, static_cast<age_idx>(a));  // for each sampled age class
	}
	return size;
}

// ----------------------------------------------------------------------------------------
// sampleSize
// ----------------------------------------------------------------------------------------
/** get the current sample size
 * check that the sampled individuls are set for this generation
 */
unsigned int
Patch::sampleSize(const sex_t SEX, const age_idx& AGE)
{
	assert((_popPtr->get_isSampled()[0] == _popPtr->getCurrentGeneration())
			&&(_popPtr->get_isSampled()[1] == _popPtr->getCurrentReplicate()));
	//	 &&(AGE & _popPtr->get_isSampled()[2]));
	return (unsigned int)_sampled_inds[SEX][AGE].size();
}

// ----------------------------------------------------------------------------------------
// set_isExtinct
// ----------------------------------------------------------------------------------------
/** set if a patch is extinced or not */
void
Patch::set_isExtinct(bool status)
{
	_isExtinct = status;

	if(status==true){          // sampling is zero in this case (no individuals)
		set_sampledInds(FEM, ALL);
		set_sampledInds(MAL, ALL);
	}
}

// ----------------------------------------------------------------------------------------
// init_containers
// ----------------------------------------------------------------------------------------
void
Patch::init_containers()
{
	_containers   = ARRAY::new_2D<vector<TIndividual*> >(2, NB_AGE_CLASSES);
	_sampled_inds = ARRAY::new_2D<vector<TIndividual*> >(2, NB_AGE_CLASSES);
	_sizes        = ARRAY::new_2D<unsigned int>(2, NB_AGE_CLASSES, (unsigned int)0);

	_sampleID = my_NAN;
	_N_sample = my_NAN;
	_friction = 1;
	for(unsigned int i = 0; i < 2; i++) {  // for each sex
		_meanVe[i]= NULL;
		_sdVe[i]  = NULL;
		_h2[i]    = NULL;
		_Ksex[i]  = 0;
		_N_ini_sex[i] = 0;
		_N_sample_sex[i] = my_NAN;
		_friction_sex[i] = 1;
	}
}

// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
Patch*
Patch::init(unsigned int id)
{
	_ID = id;
	_IDstr = toStr(id+1);
	reset_counters();

	_N_sample = _N_sample_sex[FEM] = _N_sample_sex[MAL] = my_NAN;

	for(unsigned int i=0; i < NB_AGE_CLASSES; i++) {
		_containers[MAL][i].assign( _Ksex[MAL], NULL );
		_containers[FEM][i].assign( _Ksex[FEM], NULL );
		_sizes[MAL][i] = 0;
		_sizes[FEM][i] = 0;
	}

	return this;
}

// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
Patch*
Patch::init_coal(unsigned int id)
{
	_ID = id;
	_IDstr = toStr(id+1);
	reset_counters();

	_N_sample = _N_sample_sex[FEM] = _N_sample_sex[MAL] = my_NAN;

	// a single age class, but for security reasons, initialize all
	for(unsigned int i=0; i < NB_AGE_CLASSES; i++) {
		_sizes[MAL][i] = 0;
		_sizes[FEM][i] = 0;
	}

	return this;
}

// ----------------------------------------------------------------------------------------
// set_func_pointer
// ----------------------------------------------------------------------------------------
void
Patch::set_func_pointer()
{
	func_ptr_move  = &Patch::move_ind;
	func_ptr_flush = &Patch::flush_ind;
	func_ptr_swap  = &Patch::swap_ind;
}

// ----------------------------------------------------------------------------------------
// set_func_pointer
// ----------------------------------------------------------------------------------------
void
Patch::set_func_pointer_coal()
{
	func_ptr_move  = &Patch::move_coal;
	func_ptr_flush = &Patch::flush_coal;
	func_ptr_swap  = &Patch::swap_coal;
}

// ----------------------------------------------------------------------------------------
// environmental variance
// ----------------------------------------------------------------------------------------
void Patch::set_localMeanVe (double *array, sex_t SEX)
{
	for(unsigned int i = 0; i<_nbLinkedTraits; i++) {
		_meanVe[SEX][i] = array[i];
	}
}

void Patch::set_localh2Ve (double *array, sex_t SEX)
{
	for(unsigned int i = 0; i<_nbLinkedTraits; i++) {
		_h2[SEX][i] = array[i];
	}
}

// ----------------------------------------------------------------------------------------
// stabilizing selection pressure
// ----------------------------------------------------------------------------------------
void Patch::set_localOptima (double *array, sex_t SEX)
{
	for(unsigned int i = 0; i<_nbLinkedTraits; i++) {
		if(_pSelectionType[i]->get_nb_selection_params() == 2){
			_localSelection[SEX][i][0] = array[i];     // 1. pos
		}
	}
}

void Patch::set_localIntensity (double *array, sex_t SEX)
{
	for(unsigned int i = 0; i<_nbLinkedTraits; i++) {
		if(_pSelectionType[i]->get_nb_selection_params() == 2){
			_localSelection[SEX][i][1] = array[i];     // 2. pos
		}
	}
}

// ----------------------------------------------------------------------------------------
// selecrtion coefficient selection pressure
// ----------------------------------------------------------------------------------------
void Patch::set_localSelCoefAA (double *array, sex_t SEX)
{
	for(unsigned int i = 0; i<_nbLinkedTraits; i++) {
		if(_pSelectionType[i]->get_nb_selection_params() == 2){
			_localSelection[SEX][i][0] = array[i];     // 1. pos
		}
	}
}

void Patch::set_localSelCoef (double *array, sex_t SEX)
{
	for(unsigned int i = 0; i<_nbLinkedTraits; i++) {
		if(_pSelectionType[i]->get_nb_selection_params() == 2){
			_localSelection[SEX][i][1] = array[i];     // 2. pos
		}
	}
}

// ----------------------------------------------------------------------------------------
// directional selection pressure
// ----------------------------------------------------------------------------------------
void Patch::set_localMin (double *array, sex_t SEX)
{
	for(unsigned int i = 0; i<_nbLinkedTraits; i++) {
		if(_pSelectionType[i]->get_nb_selection_params() == 5){
			_localSelection[SEX][i][0] = array[i];     // 1. pos
		}
	}
}

void Patch::set_localMax (double *array, sex_t SEX)
{
	for(unsigned int i = 0; i<_nbLinkedTraits; i++) {
		if(_pSelectionType[i]->get_nb_selection_params() == 5){
			_localSelection[SEX][i][1] = array[i];     // 1. pos
		}
	}
}

void Patch::set_localMaxGrowth (double *array, sex_t SEX)
{
	for(unsigned int i = 0; i<_nbLinkedTraits; i++) {
		if(_pSelectionType[i]->get_nb_selection_params() == 5){
			_localSelection[SEX][i][2] = array[i];     // 2. pos
		}
	}
}

void Patch::set_localGrowthRate (double *array, sex_t SEX)
{
	for(unsigned int i = 0; i<_nbLinkedTraits; i++) {
		if(_pSelectionType[i]->get_nb_selection_params() == 5){
			_localSelection[SEX][i][3] = array[i];     // 3. pos
		}
	}
}

void Patch::set_localSymmetry (double *array, sex_t SEX)
{
	for(unsigned int i = 0; i<_nbLinkedTraits; i++) {
		if(_pSelectionType[i]->get_nb_selection_params() == 5){
			_localSelection[SEX][i][4] = array[i];     // 4. pos
		}
	}
}

// ----------------------------------------------------------------------------------------
// value for a specific trait
// ----------------------------------------------------------------------------------------
void Patch::set_localOptima		 (unsigned int t, double val, sex_t SEX){
	assert(_pSelectionType[t]->get_nb_selection_params() == 2);
	_localSelection[SEX][t][0]=val;
}
void Patch::set_localIntensity (unsigned int t, double val, sex_t SEX){
	assert(_pSelectionType[t]->get_nb_selection_params() == 2);
	_localSelection[SEX][t][1]=val;
}
void Patch::set_localSelCoefAA (unsigned int t, double val, sex_t SEX){
	assert(_pSelectionType[t]->get_nb_selection_params() == 2);
	_localSelection[SEX][t][0]=val;
}
void Patch::set_localSelCoef (unsigned int t, double val, sex_t SEX){
	assert(_pSelectionType[t]->get_nb_selection_params() == 2);
	_localSelection[SEX][t][1]=val;
}
void Patch::set_localMin 		   (unsigned int t, double val, sex_t SEX){
	assert(_pSelectionType[t]->get_nb_selection_params() == 5);
	_localSelection[SEX][t][0]=val;
}
void Patch::set_localMax       (unsigned int t, double val, sex_t SEX){
	assert(_pSelectionType[t]->get_nb_selection_params() == 5);
	_localSelection[SEX][t][1]=val;
}
void Patch::set_localMaxGrowth (unsigned int t, double val, sex_t SEX){
	assert(_pSelectionType[t]->get_nb_selection_params() == 5);
	_localSelection[SEX][t][2]=val;
}
void Patch::set_localGrowthRate(unsigned int t, double val, sex_t SEX){
	assert(_pSelectionType[t]->get_nb_selection_params() == 5);
	_localSelection[SEX][t][3]=val;
}
void Patch::set_localSymmetry  (unsigned int t, double val, sex_t SEX){
	assert(_pSelectionType[t]->get_nb_selection_params() == 5);
	_localSelection[SEX][t][4]=val;
}

// ----------------------------------------------------------------------------------------
// fitness landscape selection pressure
// ----------------------------------------------------------------------------------------
/** fitness landscape must be set before phenotype landscape!
 * the size of the array is 2*nbDatapoints+1;
 * array[0]:                number of datapoints
 * array[1 -> size]:        fitness landscape
 * array[size+1 -> 2*size]: phenotype landscape
 */
void Patch::set_fitnessLandscape_fitness (unsigned int trait, double *array, unsigned int size, sex_t SEX)
{
	if(!_localSelection[SEX][trait] || _localSelection[SEX][trait][0]!=size){
		delete[] _localSelection[SEX][trait];
		_localSelection[SEX][trait] = new double[2*size+1];
		_localSelection[SEX][trait][0] = size;
	}
	for(unsigned int j=0; j<size; ++j){
		_localSelection[SEX][trait][j+size+1] = array[j];
	}
}

/** fitness landscape must be set before phenotype landscape! */
void Patch::set_fitnessLandscape_phenotype (unsigned int trait,  double *array, unsigned int size, sex_t SEX)
{
	assert(_localSelection[SEX][trait]);
	assert(_localSelection[SEX][trait][0]);
	if(_localSelection[SEX][trait][0] != size){
		error("Fitness landscape: the array sizes of phenotype (%u) and fitness landscapes (%u) do not match!",
				size, (unsigned int)_localSelection[SEX][trait][0]);
	}
	for(unsigned int j=0; j<size; ++j){
		_localSelection[SEX][trait][j+1] = array[j];
	}

	// sort the array following the phenotypes
	ARRAY::quicksortBoth(_localSelection[SEX][trait]+1, size, _localSelection[SEX][trait]+1+size);
}

// ----------------------------------------------------------------------------------------
// set_localParameter
// ----------------------------------------------------------------------------------------
void Patch::set_localParameter(double* array, sex_t sex,
		void (Patch::*pt2Func)(double*, sex_t))
{
	(this->*pt2Func)(array, sex);
}

// ----------------------------------------------------------------------------------------
// set_localParameter
// ----------------------------------------------------------------------------------------
void Patch::set_localParameter_matrix(unsigned int trait, double* array, unsigned int size, sex_t sex,
		void (Patch::*pt2Func)(unsigned int, double*, unsigned int, sex_t))
{
	(this->*pt2Func)(trait, array, size, sex);
}

// ----------------------------------------------------------------------------------------
// set_localParameter
// ----------------------------------------------------------------------------------------
void Patch::set_localParameter_matrix_ofTrait(unsigned int t, double* array, unsigned int size, sex_t sex,
		void (Patch::*pt2Func)(unsigned int, double*, unsigned int, sex_t))
{
	(this->*pt2Func)(t, array, size, sex);
}

// ----------------------------------------------------------------------------------------
// reset_counters
// ----------------------------------------------------------------------------------------
void Patch::reset_counters()
{
	nbEmigrant = 0;
	nbImmigrant = 0;
	nbKolonisers = 0;
}
// ----------------------------------------------------------------------------------------
// setNewGeneration
// ----------------------------------------------------------------------------------------
/** returns the new popualtion size */
unsigned int
Patch::setNewGeneration(age_t AGE)
{
	unsigned int i, mask = 1, size=0;
	for(i = 0; i < NB_AGE_CLASSES; i++, mask<<=1) {
		if(mask & AGE) size += setNewGeneration(static_cast<age_idx>(i));
	}
	return size;
}

// ----------------------------------------------------------------------------------------
// setNewGeneration
// ----------------------------------------------------------------------------------------
/** returns the popualtion size */
unsigned int
Patch::setNewGeneration_coal(age_t AGE)
{
	unsigned int mask = 1, size=0;
	for(unsigned int i = 0; i < NB_AGE_CLASSES; i++, mask<<=1) {
		if(mask & AGE){
			size += _sizes[FEM][static_cast<age_idx>(i)] = _N_ini_sex[FEM];
		}
	}
	return size;
}

// ----------------------------------------------------------------------------------------
// setNewGeneration
// ----------------------------------------------------------------------------------------
unsigned int
Patch::setNewGeneration(age_idx AGE)
{
	TIndividual *new_ind;
	assert(!size(FEM, AGE));
	for(unsigned int i = 0; i < _N_ini_sex[FEM]; i++) {
		new_ind = _popPtr->makeNewIndividual(0,0,FEM,this);
		new_ind->create();
		add(FEM, AGE, new_ind);
	}

	//males: same as for the females....
	assert(!size(MAL, AGE));
	for(unsigned int i = 0; i < _N_ini_sex[MAL]; i++) {
		new_ind = _popPtr->makeNewIndividual(0,0,MAL,this);
		new_ind->create();
		add(MAL, AGE, new_ind);
	}
	return (_N_ini_sex[FEM] + _N_ini_sex[MAL]);
}

// ----------------------------------------------------------------------------------------
// Patch
// ----------------------------------------------------------------------------------------
Patch::Patch(TMetapop* p, unsigned int i) : _ID(0), _ID_individual(0), _K(0), _N_ini(0),
_localSelection(0), _pSelectionType(0), _nbLinkedTraits(0), _isExtinct(1), _popPtr(0)
{
    _popPtr = p;
    _nbPatch = _popPtr->getPatchNbr();
    init_containers();
    init(i);
}

// ----------------------------------------------------------------------------------------
// ~Patch
// ----------------------------------------------------------------------------------------
Patch::~Patch()
{
	//#ifdef _DEBUG
	//  message("Patch::~Patch\n");
	//#endif
	flush();

	ARRAY::delete_2D(_containers, 2);
	ARRAY::delete_2D(_sizes, 2);
	ARRAY::delete_2D(_sampled_inds, 2);

	reset_LinkedTraits();
}

//---------------------------------------------------------------------------
///@}
//---------------------------------------------------------------------------
void
Patch::addImmigrant(const unsigned int& p, const unsigned int& s)
{
	_sizes[FEM][ADLTx] += s;
	_popPtr->add_immigrants(_ID, p, s);     // (to, from, nbImmigrant)
}

//---------------------------------------------------------------------------
///@}
//---------------------------------------------------------------------------
void
Patch::addImmigrant(Patch* p, const unsigned int& s)
{
	_sizes[FEM][ADLTx] += s;
	_popPtr->add_immigrants(_ID, p->_ID, s);     // (to, from, nbImmigrant)
}

//---------------------------------------------------------------------------
///@}
/**Returns the size of the container for the appropriate sex and age classes present in the age flag.
	@param SEX the sex class
	@param AGE the flag value of the age class
 */
unsigned int
Patch::size(sex_t SEX, age_t AGE){
	unsigned int mask = 1, s = 0, i;
	for(i = 0; i < NB_AGE_CLASSES; i++, mask <<= 1) {
		if(mask & AGE) s += _sizes[SEX][i];
	}
	return s;
}

//---------------------------------------------------------------------------
/**Returns the size of the container of the appropriate age class(es) for both sexes.
	@param AGE the flag value of the age class
 */
unsigned int
Patch::size(age_t AGE){
	return size(MAL,AGE) + size(FEM,AGE);
}

//---------------------------------------------------------------------------
/**Returns the size of the container for the appropriate sex and age class.
	@param SEX the sex class
	@param AGE the index of the age class
 */
unsigned int
Patch::size(sex_t SEX, age_idx AGE){
	return _sizes[SEX][AGE];
}

//---------------------------------------------------------------------------
/**Returns the size of the container for the appropriate age class for both sexes.
	@param AGE the index of the age class
 */
unsigned int
Patch::size(age_idx AGE){
    return _sizes[0][AGE] + _sizes[1][AGE];
}

//---------------------------------------------------------------------------
/** returns true if individual container is in agreement with size container
 */
bool
Patch::individual_container_ok(){
     for(unsigned int s, a=0; a<2; ++a){
        for(s=0; s<2; ++s){
            if(_sizes[s][a] != (unsigned int)_containers[s][a].size()){
                cout << "ERROR individual_container_ok(): patchID " << _IDstr
                << " _sizes: " << _sizes[s][a]
                << "; _container.size(): " <<_containers[s][a].size()
                << " (sex " << (s==0 ? "MAL" : "FEM")
                << "; AGE: " << (a==0 ? "OFFS" : "ADLT") << ")!" << endl;
                return false;
            }
        }
    }
    
    return true;
}

//---------------------------------------------------------------------------
/**Returns a pointer to the individual sitting at the index passed.
	@param SEX the sex class of the individual
	@param AGE the index of the age class
	@param at the index of the individual in the container
 */
TIndividual*
Patch::get(sex_t SEX, age_idx AGE, unsigned int at){
	return  _containers[SEX][AGE][at];
}

//---------------------------------------------------------------------------
/**Removes all individual pointers of the appropriate sex and age class and flush them into the recycling pool.
	\b Note: not memory operation is performed, the total amount of memory allocated is
	left untouched.
	@param SEX the sex class of the individual
	@param AGE the index of the age class
	@param pop the pointer to the metapop for access to the recycling pool
 */
void
Patch::flush_ind(sex_t SEX, age_idx AGE)
{
	for(unsigned int i = 0; i < _sizes[SEX][AGE]; ++i) {
		_popPtr->recycle(_containers[SEX][AGE][i]);    // remove the individual
	}
    _containers[SEX][AGE].clear();
	_sizes[SEX][AGE] = 0;
    assert(_containers[SEX][AGE].size() == _sizes[SEX][AGE]);
}

//---------------------------------------------------------------------------
/**Removes all individual pointers of the appropriate sex and age class and flush them into the recycling pool.
	@param AGE an unsigned int containing the flags of the age classes to flush
	@param pop the pointer to the metapop for access to the recycling pool
	@see flush()
 */
void
Patch::flush(sex_t SEX, age_idx AGE)
{
	(this->*func_ptr_flush)(SEX, AGE);
}

//---------------------------------------------------------------------------
/**Removes all individual pointers of the appropriate sex and age class and flush them into the recycling pool.
	@param AGE an unsigned int containing the flags of the age classes to flush
	@param pop the pointer to the metapop for access to the recycling pool
	@see flush()
 */
void
Patch::flush(age_idx AGE)
{
	(this->*func_ptr_flush)(FEM, AGE);
	(this->*func_ptr_flush)(MAL, AGE);
}

//---------------------------------------------------------------------------
void
Patch::flush(age_t AGE)
{
	for(unsigned int mask=1, i=0; i<NB_AGE_CLASSES; i++, mask <<= 1) {
		if(mask & AGE) {
			(this->*func_ptr_flush)(MAL, static_cast<age_idx>(i));
			(this->*func_ptr_flush)(FEM, static_cast<age_idx>(i));
		}
	}
}

//---------------------------------------------------------------------------
/**Removes all individual pointers of all sex and age classes and flush them into the recycling pool.*/
void
Patch::flush(){
	for(unsigned int i = 0; i < NB_AGE_CLASSES; i++)
	{
		flush(static_cast<age_idx>(i));
	}
}

//---------------------------------------------------------------------------
/**Modifies the appropriate container with value of the pointer given.
	@param SEX the sex class of the individual
	@param AGE the index of the age class
	@param at the index of the individual in the container
	@param ind the pointer to the individual
 */
void
Patch::set(sex_t SEX, age_idx AGE, unsigned int at, TIndividual* ind)
{
	cout << "\nPatch::set(): not yet changed!\n" << endl;
	//_containers[SEX][AGE][at] = ind;
}

//---------------------------------------------------------------------------
/**Adds an individual to the appropriate container, increments its size, eventually resizing it.
	@param SEX the sex class of the individual
	@param AGE the index of the age class
	@param ind the pointer to the individual
 */
void
Patch::add(const sex_t& SEX, const age_idx& AGE, TIndividual* idx)
{
	assert(_sizes[SEX][AGE]<=_containers[SEX][AGE].size());
	if(_sizes[SEX][AGE] < _containers[SEX][AGE].size()) _containers[SEX][AGE][_sizes[SEX][AGE]] = idx;
	else     _containers[SEX][AGE].push_back(idx);
	++_sizes[SEX][AGE];
    assert((unsigned int)_containers[SEX][AGE].size()==_sizes[SEX][AGE]);
}

//---------------------------------------------------------------------------
/**Assigns a new container of given size for the sex and age class passed, sets all values to NULL.*/
void
Patch::assign(sex_t SEX, age_idx AGE, unsigned int n)
{
	_containers[SEX][AGE].assign(n,0);
	_sizes[SEX][AGE] = 0;
    assert((unsigned int)_containers[SEX][AGE].size()==_sizes[SEX][AGE]);
}

//---------------------------------------------------------------------------
/**Removes the individual sitting at the given index in the appropriate container,
 * but does not delete the individual.
	@param SEX the sex class of the individual
	@param AGE the index of the age class
	@param at the index of the individual in the container
 */
void
Patch::remove(sex_t SEX, age_idx AGE, unsigned int at)
{
	assert(at<_sizes[SEX][AGE]);
	_containers[SEX][AGE][at] = _containers[SEX][AGE][_sizes[SEX][AGE] - 1];  // swap the last individual to the position at 'at'
    _containers[SEX][AGE].pop_back();
	_sizes[SEX][AGE]--;
    assert((unsigned int)_containers[SEX][AGE].size()==_sizes[SEX][AGE]);
}

//---------------------------------------------------------------------------
/**Removes the individual in the appropriate container (function is not very efficient).
	@param SEX the sex class of the individual
	@param AGE the index of the age class
	@param ind pointer of an individual to remove
 */
void
Patch::remove(sex_t SEX, age_idx AGE, TIndividual* ind)
{
	for(unsigned int i=0; i<_sizes[SEX][AGE]; ++i){
		if(ind == _containers[SEX][AGE][i]){
			remove(SEX, AGE, i);
			return;
		}
	}
	assert(1!=1); 	// function should never pass here
}

//---------------------------------------------------------------------------
/**Removes and deletes the individual sitting at the given index in the appropriate container.
	@param SEX the sex class of the individual
	@param AGE the index of the age class
	@param at the index of the individual in the container
 */
void
Patch::recycle(sex_t SEX, age_idx AGE, unsigned int at)
{
	assert(at<_sizes[SEX][AGE]);
 	_popPtr->recycle(_containers[SEX][AGE][at]);    // remove the individual
	_containers[SEX][AGE][at] = _containers[SEX][AGE][_sizes[SEX][AGE] - 1];  // swap the last individual to the position at 'at'
    _containers[SEX][AGE].pop_back();
	_sizes[SEX][AGE]--;
    assert(_sizes[SEX][AGE] == _containers[SEX][AGE].size());
}

//---------------------------------------------------------------------------
/**Removes and deletes the individual in the appropriate container (function is not very efficient).
	@param SEX the sex class of the individual
	@param AGE the index of the age class
	@param ind pointer of an individual to remove
 */
void
Patch::recycle(sex_t SEX, age_idx AGE, TIndividual* ind)
{
	for(unsigned int i=0; i<_sizes[SEX][AGE]; ++i){
		if(ind == _containers[SEX][AGE][i]){
			recycle(SEX, AGE, i);
			break;
		}
	}
}

//---------------------------------------------------------------------------
// MOVE: does not overwrite the individuals in the to container
//---------------------------------------------------------------------------
/**Moves all individual from an age class to an other one.
	\b Note: both containers are transformed by this operation. The 'from'
					 container size is set to zero while the 'to' container size
					 is increased by the number of the 'from' container.
	@param from the original age class of the individual
	@param to the destination age class of the individual
 */
void
Patch::move(age_idx from, age_idx to)
{
	(this->*func_ptr_move)(FEM, from, to);
	(this->*func_ptr_move)(MAL, from, to);
}

//---------------------------------------------------------------------------
/**Moves all individual from an age class to an other one.
	\b Note: both containers are transformed by this operation. The 'from'
					 container size is set to zero while the 'to' container size
					 is increased by the number of the 'from' container.
	@param from the original age class of the individual
	@param to the destination age class of the individual
 */
void
Patch::move(sex_t SEX, age_idx from, age_idx to)
{
	(this->*func_ptr_move)(SEX, from, to);
}

//---------------------------------------------------------------------------
/**Moves all individual from an age class to another one.
	\b Note: both containers are transformed by this operation. The 'from'
					 container size is set to zero while the 'to' container size
					 is increased by the number of the 'from' container.
	@param SEX the sex class of the individual
	@param from the original age class of the individual
	@param to the destination age class of the individual
 */
void
Patch::move_ind(sex_t SEX, age_idx from, age_idx to)
{
	for(int i=0, size=_sizes[SEX][from]; i<size; ++i){
		move(SEX, from, to, 0);
	}
}

//---------------------------------------------------------------------------
/**Moves all individual from an age class to another one.
	\b Note: both containers are transformed by this operation. The 'from'
					 container size is set to zero while the 'to' container size
					 is increased by the number of the 'from' container.
	@param SEX the sex class of the individual
	@param from the original age class of the individual
	@param to the destination age class of the individual
 */
void
Patch::move_coal(sex_t SEX, age_idx from, age_idx to)
{
	_sizes[SEX][to] += _sizes[SEX][from];
	_sizes[SEX][from] = 0;
}

//---------------------------------------------------------------------------
/**Moves an individual from an age class to an other one.
	\b Note: both containers are transformed by this operation. The 'from'
					 container size is reduced by one while the 'to' container size
					 is increased by one.
	@param SEX the sex class of the individual
	@param from the original age class of the individual
	@param to the destination age class of the individual
	@param at the index of the individual in the container
 */
void
Patch::move(sex_t SEX, age_idx from, age_idx to, unsigned int at)
{
	add(SEX, to, _containers[SEX][from][at]);
	remove(SEX, from, at);
}

//---------------------------------------------------------------------------
/** this function removes randomly all until 1-n% individuals  of the given age and sex remain */
void
Patch::survive_randomly_inds_relative(sex_t SEX, age_idx AGE, double ratio)
{
	assert(ratio>=0 && ratio<=1);
	unsigned int size = _sizes[SEX][AGE];

	remove_randomly_inds_absolute(SEX, AGE, size - my_round(size*ratio));
}

//---------------------------------------------------------------------------
/** this function removes randomly all until n individuals  of the given age and sex remain */
void
Patch::survive_randomly_inds_absolute(sex_t SEX, age_idx AGE, unsigned int size)
{
	assert(size==0 || size>=1);
	if(_sizes[SEX][AGE]>size) remove_randomly_inds_absolute(SEX, AGE, _sizes[SEX][AGE] - size);
}

//---------------------------------------------------------------------------
/** this function removes randomly n% individuals  of the given age and sex */
void
Patch::remove_randomly_inds_relative(sex_t SEX, age_idx AGE, double ratio)
{
	assert(ratio>=0 && ratio<=1);
	remove_randomly_inds_absolute(SEX, AGE, my_round(_sizes[SEX][AGE]*ratio));
}

//---------------------------------------------------------------------------
/** this function removes randomly n% individuals  of the given age and sex */
void
Patch::remove_randomly_inds_absolute(sex_t SEX, age_idx AGE, unsigned int size)
{
	assert(size==0 || size>=1);
	vector<unsigned int> vIndex = _popPtr->rand().sample<unsigned int>(0, _sizes[SEX][AGE], size);  // get random indexes
	sort(vIndex.rbegin(), vIndex.rend());           // start form the back
	for(unsigned int i=0; i<size; ++i){
		assert(vIndex[i]<_sizes[SEX][AGE]);
		recycle(SEX, AGE, vIndex[i]);                 // that is a looser
	}
}

//---------------------------------------------------------------------------
/** everything moves by a generation. The oldest age class will be flushed
 * and the youngest age class will be empty
 */
void
Patch::aging()
{
	for(unsigned int s=0; s<2; ++s){ // for each sex
		_containers[s][ADLTx].clear();
		_sizes[s][ADLTx] = 0;
	}

	/*
	unsigned int a, s;
	for(s=0; s<2; ++s){ // for each sex
		// oldest age class: remove them (can be ignored here)
		a=NB_AGE_CLASSES-1;
		vector<unsigned int>& temp = _containers[s][a];
		temp.clear();

		// for each other age class: move them one up
		for(; a>0; --a){
			_containers[s][a] = _containers[s][a-1];      // is a vector
			_sizes[s][a]      = _sizes[s][a-1];           // is a number
		}

		// set the youngest age class to zero, i.e. setting just the size to zero is enough
		_sizes[s][0] = 0;
		_containers[s][0] = temp;
	}
	 */
}

//---------------------------------------------------------------------------
// SWAP: overwrites the individuals in the to container, thus the to container has to be empty
//---------------------------------------------------------------------------
/**Copies all elements in the 'from' age-class container to the 'to' age-class container.
	The previous elements of the 'to' container are overwritten, it is thus worth considering using flush()
	before swaping to avoid memory leaks! The size of the 'from' container is set to 0.
	@param SEX the sex class of the individual
	@param from the original age class of the individual
	@param to the destination age class of the individual
 */
void
Patch::swap(age_idx from, age_idx to)
{
	(this->*func_ptr_swap)(FEM, from, to);
	(this->*func_ptr_swap)(MAL, from, to);
}

//---------------------------------------------------------------------------
// SWAP: overwrites the individuals in the to container, thus the to container has to be empty
//---------------------------------------------------------------------------
/**Copies all elements in the 'from' age-class container to the 'to' age-class container.
	The previous elements of the 'to' container are overwritten, it is thus worth considering using flush()
	before swaping to avoid memory leaks! The size of the 'from' container is set to 0.
	@param SEX the sex class of the individual
	@param from the original age class of the individual
	@param to the destination age class of the individual
 */
void
Patch::swap(sex_t SEX, age_idx from, age_idx to)
{
	(this->*func_ptr_swap)(SEX, from, to);
}

//---------------------------------------------------------------------------
/**Copies all elements in the 'from' age-class container to the 'to' age-class container of the same sex.
	The previous elements of the 'to' container are overwritten, it is thus worth considering using flush()
	before swaping to avoid memory leaks! The size of the 'from' container is set to 0.
	@param SEX the sex class of the individual
	@param from the original age class of the individual
	@param to the destination age class of the individual
 */
void
Patch::swap_ind(sex_t SEX, age_idx from, age_idx to)
{
	assert(!_sizes[SEX][to]);

	//store values of the "to" container temporarily
	vector<TIndividual*> temp_container = _containers[SEX][to];

	// reassign the "from" container
	_containers[SEX][to]    = _containers[SEX][from];
	_sizes[SEX][to]         = _sizes[SEX][from];

	// emtpy the "to" container
	_containers[SEX][from]  = temp_container;
	_sizes[SEX][from]       = 0;
}

//---------------------------------------------------------------------------
/**Copies all elements in the 'from' age-class container to the 'to' age-class container of the same sex.
	The previous elements of the 'to' container are overwritten, it is thus worth considering using flush()
	before swaping to avoid memory leaks! The size of the 'from' container is set to 0.
	@param SEX the sex class of the individual
	@param from the original age class of the individual
	@param to the destination age class of the individual
 */
void
Patch::swap_coal(sex_t SEX, age_idx from, age_idx to)
{
	assert(!_sizes[SEX][to]);
	_sizes[SEX][to]         = _sizes[SEX][from];                   // reassign the "from" container
	_sizes[SEX][from]       = 0;                                   // emtpy the "to" container
}

//---------------------------------------------------------------------------
/**Sets the size of the appropriate container to zero.
	\b Note: no memory operation is performed, the capacity of the container is thus not affected.
					The individual pointers are not flushed to the recycling pool. It is thus a good idea
					to consider using Patch::flush to be sure no pointers remained in the container.
	@see flush()
	@param SEX the sex class
	@param AGE the index of the age class
 */
void
Patch::clear(sex_t SEX, age_idx AGE)
{
	_sizes[SEX][AGE] = 0;
}

//---------------------------------------------------------------------------
void
Patch::clear(age_idx AGE)
{
	clear(MAL, AGE);
	clear(FEM, AGE);
}

//---------------------------------------------------------------------------
void
Patch::clear()
{
	for(unsigned int i = 0; i < NB_AGE_CLASSES; i++){
		clear(static_cast<age_idx>(i));
	}
}

//---------------------------------------------------------------------------
/** return an array with all phenotypes of the specific trait.
 * The array has to be deleted after use.
 */
double*
Patch::getPhenotypes(sex_t SEX, age_idx AGE, const int& trait)
{
	unsigned int i, size=_sizes[SEX][AGE];
	double* pheno = new double[size];
	for (i = 0; i < size; ++i) {
		pheno[i] = _containers[SEX][AGE][i]->getTraitPhenotype(trait);
	}
	return pheno;
}

//---------------------------------------------------------------------------
/** downregulates the pop size to the given K if necessary. The survivors are drawn
 * randomly in relation to their fitness
 */
void
Patch::regulate_selection_fitness(const unsigned int& K, TSelection* pSel, sex_t SEX, age_idx AGE)
{
	// if K is zero: remove all individuals
	if(!K){
		flush(SEX, AGE);
		return;
	}

	unsigned int nbInd = pSel->get_nbInd(SEX);
	if(nbInd>K){
		if((nbInd-K)>K/2) regulation_selection_draw_survivors(pSel, K, SEX, AGE);
		else              regulation_selection_draw_loosers(pSel, nbInd-K, SEX, AGE);
	}
}

//---------------------------------------------------------------------------
void
Patch::regulation_selection_draw_survivors(TSelection* pSelection, const unsigned int& K, sex_t SEX, age_idx AGE)
{
	assert(K && K<=size(SEX, AGE));
	// make the array cumulative to draw the MOST fittest
	pSelection->sort_fitness(SEX, -3, K);
	TIndividual** aInd = pSelection->get_aInd(SEX);

	// remove the loosers (second part of the array)
	for(unsigned int i = K, N = size(SEX, AGE); i<N; ++i){
		recycle(SEX, AGE, aInd[i]);     // remove the "unsorted" individuals
	}
}

//---------------------------------------------------------------------------
void
Patch::regulation_selection_draw_loosers(TSelection* pSelection, const unsigned int& K, sex_t SEX, age_idx AGE)
{
	assert(K && K<=size(SEX, AGE));
	// make the array cumulative to draw the MOST fittest
	pSelection->sort_fitness(SEX, 3, K);
	TIndividual** aInd = pSelection->get_aInd(SEX);

	// remove the loosers (first sorted part of the array)
	for(unsigned int i = 0; i<K; ++i){
		recycle(SEX, AGE, aInd[i]);     // remove the "sorted" individuals
	}
}

//---------------------------------------------------------------------------
/** viability selection for the specified sex and age */
void
Patch::regulation_selection_hard(age_idx AGE)
{
	_popPtr->get_pSelection()->set_fitness(this, AGE);    // compute the fitnesses, but don't sort or make yet the array cumulative
	unsigned int SEX, nbInd;
	TIndividual** aInd;
	double*      aFit;

	for(SEX=0; SEX<2; ++SEX){                     // for each sex
		aInd = _popPtr->get_pSelection()->get_aInd((sex_t)SEX);
		aFit = _popPtr->get_pSelection()->get_aFit((sex_t)SEX);
		nbInd = size((sex_t)SEX, AGE);

		for(unsigned int i = 0; i<nbInd; ++i){      // for each individual
			if(aFit[i] < _popPtr->rand().Uniform()) recycle((sex_t)SEX, AGE, aInd[i]);     // does the individual survive?
		}
	}
}

//---------------------------------------------------------------------------




