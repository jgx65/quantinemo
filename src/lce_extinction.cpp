/** @file LCE_extinction.cpp
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

#include "lce_extinction.h"

#include <algorithm>

using namespace std;

// -----------------------------------------------------------------------------
// constructor
// -----------------------------------------------------------------------------
LCE_Extinction::LCE_Extinction(int rank) :
LCE("extinction", "population extinction", "", rank),
_Xtion_rate(0), _Xtion_rate_size(0) {
	for (unsigned int i = 0; i < 2; ++i) {
		_survival_rate[i] = NULL;
		_survival_rate_size[i] = 0;
	}

    add_parameter("extinction_rate", DBL_MAT, true, 0, 1, "0", true,
                        "Extinction probability of a patch at each generation", 0);
    
	add_parameter("extinction_rate_survival", DBL_MAT, false, 0, my_NAN,
                        "0", true,
                        "How serious is an extinction:\n"\
                        "  integer: absolute number of individuals will survive.\n" \
                        "  decimal: the given population ratio will survive.\n" \
                        "The survivors are drawn randomly.", 2);
    
	add_parameter("extinction_rate_survival_fem", DBL_MAT, false, 0,
                        my_NAN, "0", true,
                        "How serious is an extinction for females:\n"\
                        "  integer: absolute number of females will survive.\n" \
                        "  decimal: the given female ratio will survive.\n" \
                        "The survivors are drawn randomly.", 2);
    
	add_parameter("extinction_rate_survival_mal", DBL_MAT, false, 0,
                  my_NAN, "0", true,
                        "How serious is an extinction for males:\n"\
                        "  integer: absolute number of males will survive.\n" \
                        "  decimal: the given male ratio will survive.\n" \
                        "The survivors are drawn randomly.", 2);
}

// -----------------------------------------------------------------------------
// LCE_Extinction::extinction_homogenous_low
// -----------------------------------------------------------------------------
/** LOW extinction rate is common among patches: draw first the number of total extinction
	* before drawing the patch to get extinct (more efficient since normaly extinctions are rare)
	* exticntions may be partial
	*/
void
LCE_Extinction::extinction_homogenous_low_partial_1sex()
{
	#ifdef _DEBUG
	 message("  LCE_Extinction (homogenous low, partial) ... ");
	#endif

	unsigned int nbPatch, nbExt, i, rand;
	Patch *current_patch;

	// get the number of extinctions
	nbPatch = _popPtr->get_nbFullPatch();
	nbExt = (unsigned int)get_pop_ptr()->rand().Binomial(*_Xtion_rate, nbPatch);
    vector<Patch*>::iterator curPos, endPos; // not really used
    vector<unsigned int> PatchGoingExtinct = get_pop_ptr()->rand().Uniforms(nbPatch, nbExt);
	if (nbExt) {
		for (i = 0; i < nbExt;) {
			rand = get_pop_ptr()->rand().Uniform(nbPatch);                               // get a patch randomly
			current_patch = _popPtr->get_vFullPatch()[PatchGoingExtinct[i]];
			if(!(this->*survivors_func_ptr[FEM])(current_patch, FEM)){          // extinction female
				assert(!current_patch->size());
				curPos = _popPtr->get_vFullPatch().begin() + rand;              // curPos must be a non-const variable...
                _popPtr->new_emptyPatch(curPos, endPos);                        // remove if now empty
				--nbPatch;                                                      // decrement the number of patches
			}
			++i;        // increment the number of extinctions
		}
	}

	#ifdef _DEBUG
	 message("done! (Extinct patches: %i; Patches: %i; Individuals(F/M): [%i/%i; %i/%i])\n",
			nbExt, _popPtr->get_nbFullPatch(),
			_popPtr->size(FEM, OFFSx),  _popPtr->size(MAL, OFFSx),
			_popPtr->size(FEM, ADLTx),  _popPtr->size(MAL, ADLTx));
	#endif
}

// -----------------------------------------------------------------------------
// LCE_Extinction::extinction_homogenous_low
// -----------------------------------------------------------------------------
/** LOW extinction rate is common among patches: draw first the number of total extinction
	* before drawing the patch to get extinct (more efficient since normaly extinctions are rare)
	* exticntions may be partial
	* the patches mnay only be delted later, since not all patches get entirely  extinced
	* which are hit by an exticntion
	*/
void
LCE_Extinction::extinction_homogenous_low_partial_2sex()
{
	#ifdef _DEBUG
	 message("  LCE_Extinction (homogenous low, partial) ... ");
	#endif

	unsigned int nbPatch, nbExt, i, j, rand;
	Patch *current_patch;

	// get the number of extinctions
	nbPatch = _popPtr->get_nbFullPatch();
	nbExt = (unsigned int)get_pop_ptr()->rand().Binomial(*_Xtion_rate, nbPatch);
	if (nbExt) {
		unsigned int* extinctions = new unsigned int[nbExt]; // to store the patches of extinction in order to prevent to extinct twice a patch
		vector<unsigned int> toDelete;
		for (i = 0; i < nbExt;) {
			// get a patch randomly
			rand = get_pop_ptr()->rand().Uniform(nbPatch);
			for (j = 0; j < i; ++j) {               // check if this patch was already drawn
				if (extinctions[j] == rand) break;
			}
			if (j < i) continue; 		// the patch was already choosen
            
			// perform extinction
			extinctions[i] = rand; // add the patch to the temp array
			current_patch = _popPtr->get_vFullPatch()[rand];
			if(!((this->*survivors_func_ptr[FEM])(current_patch, FEM)            // extinction female
				 + (this->*survivors_func_ptr[MAL])(current_patch, MAL))){         // extinction male
				assert(!current_patch->size());
				toDelete.push_back(rand);                                          // patch is empty and has to be delted later
			}
			++i;        // increment the number of extinctions
		}

		// remove the empty patches
        vector<Patch*>::iterator curPos, endPos;
        vector<unsigned int>::iterator curDel, endDel;
		sort(toDelete.begin(), toDelete.end(), std::greater<unsigned int>());   // start from behind, otherwise the index is wrong
		for(curDel=toDelete.begin(), endDel=toDelete.end(); curDel!=endDel; ++curDel){
			assert(!_popPtr->get_vFullPatch()[toDelete.back()]->size());
			curPos = _popPtr->get_vFullPatch().begin() + toDelete.back();   // curPos must be a non-const variable...
            _popPtr->new_emptyPatch(curPos, endPos);                        // remove if now empty
        }

		delete[] extinctions;
	}

	#ifdef _DEBUG
    message("done! (%i extinct patches; Patches: %i; Individuals(F/M): [%i/%i; %i/%i])\n",
			nbExt, _popPtr->get_nbFullPatch(),
			_popPtr->size(FEM, OFFSx),  _popPtr->size(MAL, OFFSx),
			_popPtr->size(FEM, ADLTx),  _popPtr->size(MAL, ADLTx));
	#endif
}

// -----------------------------------------------------------------------------
// LCE_Extinction::extinction_homogenous_low
// -----------------------------------------------------------------------------
/** LOW extinction rate is common among patches: draw first the number of total extinction
	* and entire extinction
	* before drawing the patch to get extinct (more efficient since normaly extinctions are rare)
	*/
void
LCE_Extinction::extinction_homogenous_low_total()
{
	#ifdef _DEBUG
	 message("  LCE_Extinction (homogenous low, total) ... ");
	#endif

	// get the number of extinctions
	unsigned int nbPatch = _popPtr->get_nbFullPatch();
   	unsigned int nbExt = (unsigned int)get_pop_ptr()->rand().Binomial(*_Xtion_rate, nbPatch);
    vector<unsigned int> PatchGoingExtinct = get_pop_ptr()->rand().Uniforms(nbPatch, nbExt);
    if (nbExt) {
        vector<Patch*>::iterator curPos, endPos = _popPtr->get_vFullPatch().end();
		for (unsigned int i = 0; i < nbExt; ++i, --nbPatch) {
            curPos = _popPtr->get_vFullPatch().begin() + PatchGoingExtinct[i]; // curPos must be a non-const variable...
            (*curPos)->flush();
            _popPtr->new_emptyPatch(curPos, endPos);                        // remove if now empty
        }
	}

	#ifdef _DEBUG
    message("done! (Extinct patches: %i; Patches: %i; Individuals(F/M): [%i/%i; %i/%i])\n",
			nbExt, _popPtr->get_nbFullPatch(),
			_popPtr->size(FEM, OFFSx),  _popPtr->size(MAL, OFFSx),
			_popPtr->size(FEM, ADLTx),  _popPtr->size(MAL, ADLTx));
	#endif
}

// -----------------------------------------------------------------------------
// LCE_Extinction::extinction_homogenous_high
// -----------------------------------------------------------------------------
/** HIGH extinction rate is common among patches and extinction may be partial:
	* more efficient to check each popualted patch seaprately for extinction
	*/
void
LCE_Extinction::extinction_homogenous_high_partial()
{
	#ifdef _DEBUG
	 message("  LCE_Extinction (homogenous high, partial) ... ");
	 unsigned int nbExt = 0;
	#endif

	// for each popualted patch
	vector<Patch*>::iterator curPop = _popPtr->get_vFullPatch().begin();
	vector<Patch*>::iterator endPop = _popPtr->get_vFullPatch().end();
	for (; curPop != endPop;) {
		if (get_pop_ptr()->rand().Uniform() >= *_Xtion_rate) {
            ++curPop;
            continue;
        }
        
		if(!((this->*survivors_func_ptr[FEM])(*curPop, FEM)          // extinction female
			 + (this->*survivors_func_ptr[MAL])(*curPop, MAL))){       // extinction male
			assert(!(*curPop)->size());
            _popPtr->new_emptyPatch(curPop, endPop);   // remove if now empty (curPop and endPop are adjusted)
		}

		#ifdef _DEBUG
		 ++nbExt;
		#endif
	}

	#ifdef _DEBUG
    message("done! (Extinct patches: %i; Patches: %i; Individuals(F/M): [%i/%i; %i/%i])\n",
			nbExt, _popPtr->get_nbFullPatch(),
			_popPtr->size(FEM, OFFSx),  _popPtr->size(MAL, OFFSx),
			_popPtr->size(FEM, ADLTx),  _popPtr->size(MAL, ADLTx));
	#endif
}

// -----------------------------------------------------------------------------
// LCE_Extinction::extinction_homogenous_high
// -----------------------------------------------------------------------------
/** HIGH extinction rate is common among patches and extinction may be partial:
	* more efficient to check each popualted patch seaprately for extinction
	* extinction are total
	*/
void
LCE_Extinction::extinction_homogenous_high_total()
{
	#ifdef _DEBUG
	 message("  LCE_Extinction (homogenous high, total) ... ");
	 unsigned int nbExt = 0;
	#endif

	// for each popualted patch
	vector<Patch*>::iterator curPop = _popPtr->get_vFullPatch().begin();
	vector<Patch*>::iterator endPop = _popPtr->get_vFullPatch().end();
	for (; curPop != endPop;) {
		if (get_pop_ptr()->rand().Uniform() >= *_Xtion_rate){
            ++curPop;
            continue;
        }

		(*curPop)->flush();
		_popPtr->new_emptyPatch(curPop, endPop);   // remove if now empty (curPop and endPop are adjusted)

		#ifdef _DEBUG
		 ++nbExt;
		#endif
	}

	#ifdef _DEBUG
    message("done! (Extinct patches: %i; Patches: %i; Individuals(F/M): [%i/%i; %i/%i])\n",
			nbExt, _popPtr->get_nbFullPatch(),
			_popPtr->size(FEM, OFFSx),  _popPtr->size(MAL, OFFSx),
			_popPtr->size(FEM, ADLTx),  _popPtr->size(MAL, ADLTx));
	#endif
}

// -----------------------------------------------------------------------------
// LCE_Extinction::extinction_variable
// -----------------------------------------------------------------------------
/** extinction rate varies among patches and extinction may be partial
	* each popuated patch has to be checked for extinction
	*/
void
LCE_Extinction::extinction_variable_partial()
{
	#ifdef _DEBUG
	 message("  LCE_Extinction (variable, partial)... ");
	 unsigned int nbExt = 0;
	#endif

	// for each populated patch (going backwards avoids deletion problems)
	vector<Patch*>::iterator curPop = _popPtr->get_vFullPatch().begin();
	vector<Patch*>::iterator endPop = _popPtr->get_vFullPatch().end();
	for (; curPop != endPop;) {
		if (get_pop_ptr()->rand().Uniform() >= _Xtion_rate[(*curPop)->get_ID() % _Xtion_rate_size]){
            ++curPop;
            continue;
        }

		if(!((this->*survivors_func_ptr[FEM])(*curPop, FEM)          // extinction female
			 + (this->*survivors_func_ptr[MAL])(*curPop, MAL))){       // extinction male
			assert(!(*curPop)->size());
			_popPtr->new_emptyPatch(curPop, endPop);   // remove if now empty (curPop and endPop are adjusted)
		}

		#ifdef _DEBUG
		 ++nbExt;
		#endif
	}

	#ifdef _DEBUG
    message("done! (Extinct patches: %i; Patches: %i; Individuals(F/M): [%i/%i; %i/%i])\n",
			nbExt, _popPtr->get_nbFullPatch(),
			_popPtr->size(FEM, OFFSx),  _popPtr->size(MAL, OFFSx),
			_popPtr->size(FEM, ADLTx),  _popPtr->size(MAL, ADLTx));
	#endif
}

// -----------------------------------------------------------------------------
// LCE_Extinction::extinction_variable
// -----------------------------------------------------------------------------
/** extinction rate varies among patches, and there is no partial extinction
  * each popuated patch has to be checked for extinction
	*/
	//#include <fstream>
	using namespace std;
void
LCE_Extinction::extinction_variable_total()
{
	#ifdef _DEBUG
	 message("  LCE_Extinction (variable total)... ");
	 unsigned int nbExt = 0;
	#endif

	// for each populated patch (going backwards avoids deletion problems)
	vector<Patch*>::iterator curPop = _popPtr->get_vFullPatch().begin();
	vector<Patch*>::iterator endPop = _popPtr->get_vFullPatch().end();
	for (; curPop != endPop;) {
		if (get_pop_ptr()->rand().Uniform() >= _Xtion_rate[(*curPop)->get_ID() % _Xtion_rate_size]){
            ++curPop;
            continue;
        }
            
		(*curPop)->flush();
		_popPtr->new_emptyPatch(curPop, endPop);   // remove if now empty (curPop and endPop are adjusted)

		#ifdef _DEBUG
		 ++nbExt;
		#endif
	}

	#ifdef _DEBUG
    message("done! (Extinct patches: %i; Patches: %i; Individuals(F/M): [%i/%i; %i/%i])\n",
			nbExt, _popPtr->get_nbFullPatch(),
			_popPtr->size(FEM, OFFSx),  _popPtr->size(MAL, OFFSx),
			_popPtr->size(FEM, ADLTx),  _popPtr->size(MAL, ADLTx));
	#endif
}

// -----------------------------------------------------------------------------
// LCE_Extinction::survivors
// -----------------------------------------------------------------------------
/** if extinction hits the patch what to do */
unsigned int         // survivors are only absolutely defined (identical among patches)
LCE_Extinction::survivors_absolute_const(Patch* curPatch, const sex_t & SEX)
{
	assert(_survival_rate[SEX][0]==0 || _survival_rate[SEX][0]>=1);
	curPatch->survive_randomly_inds_absolute(SEX, ADLTx, _survival_rate[SEX][0]);
	return _survival_rate[SEX][0];
}

unsigned int         // survivors are only absolutely defined (patch sepecific settings)
LCE_Extinction::survivors_absolute_var(Patch* curPatch, const sex_t & SEX)
{
	unsigned int size = _survival_rate[SEX][curPatch->get_ID() % _survival_rate_size[SEX]];
	assert(size==0 || size>=1);
	curPatch->survive_randomly_inds_absolute(SEX, ADLTx, size);
	return size;
}

unsigned int         // survivors are only relativly defined (identical among patches)
LCE_Extinction::survivors_relative_const(Patch* curPatch, const sex_t & SEX)
{
	assert(_survival_rate[SEX][0]>=0 && _survival_rate[SEX][0]<=1);
	curPatch->survive_randomly_inds_relative(SEX, ADLTx, _survival_rate[SEX][0]);
	return curPatch->size(SEX, ADLTx);
}

unsigned int         // survivors are only relativly defined (patch specific settings)
LCE_Extinction::survivors_relative_var(Patch* curPatch, const sex_t & SEX)
{
	double value = _survival_rate[SEX][curPatch->get_ID() % _survival_rate_size[SEX]];
	assert(value>=0 && value<=1);
	curPatch->survive_randomly_inds_relative(SEX, ADLTx, value);
	return curPatch->size(SEX, ADLTx);
}

unsigned int        // survivors are specified relativly and absolutely defined (patch specific settings)
LCE_Extinction::survivors_mixed(Patch* curPatch, const sex_t & SEX)
{
	double value = _survival_rate[SEX][curPatch->get_ID() % _survival_rate_size[SEX]];
	if(value>=1 || value==0) curPatch->survive_randomly_inds_absolute(SEX, ADLTx, value);
	else                     curPatch->survive_randomly_inds_relative(SEX, ADLTx, value);
	return curPatch->size(SEX, ADLTx);
}

unsigned int        // all individuals  of this sex die (identical among patches)
LCE_Extinction::survivors_none(Patch* curPatch, const sex_t & SEX)
{
	curPatch->flush(SEX, ADLTx);
	return 0;
}

unsigned int        // all survive of this sex (nothing to do) (identical among patches)
LCE_Extinction::survivors_all(Patch* curPatch, const sex_t & SEX)
{
	return curPatch->size(SEX, ADLTx);
}

// -----------------------------------------------------------------------------
// LCE_Extinction::init
// -----------------------------------------------------------------------------
/** returns true if at any stage extinction may occur */
bool LCE_Extinction::init(Metapop* popPtr)
{
	LCE::init(popPtr);
	if (!set_extinction_rate()) {
		_paramSet->set_isSet(false);       // extinction rate is not set
		return false;
	}
	set_survival_rate();
	return true;
}

// -----------------------------------------------------------------------------
// LCE_Extinction::set_extinction_rate
// -----------------------------------------------------------------------------
/** returns false if the parameter is not set or the default value is set */
bool LCE_Extinction::set_extinction_rate()
{
	if (_Xtion_rate) {delete[]_Xtion_rate; _Xtion_rate = NULL; _Xtion_rate_size=0;}

	Param* param = this->get_parameter("extinction_rate");
	if (!param->is_matrix()) { // if it is a single value or not used
		double val = param->get_value();
		if (!val){                                   // currently no extinction may occur
			if(param->isTemporalParam()) return true;  // maybe later extinctions may occur
			else                         return false; // never an extinction may occur
		}

		_Xtion_rate = new double[1];
		_Xtion_rate_size = 1;
		_Xtion_rate[0] = val;
		if(val>0.5) extinction_func_ptr = &LCE_Extinction::extinction_homogenous_high_total;
		else        extinction_func_ptr = &LCE_Extinction::extinction_homogenous_low_total;
		return true;
	}

	// if it is a matrix
	TMatrix* m = param->get_matrix();
	double* array = m->get();
	unsigned int count_patch = m->length();

	// check the number of values per number of patches
	unsigned int nbPatch = _popPtr->get_nbPatch();
	if (count_patch > nbPatch) warning("There are more extinction rates specified than patches! Only a part of the extinction rates  is considered!\n");
	else if (nbPatch % count_patch) 	warning("The number of extinction rates is not a entire subset of the number of patches!\n");

	_Xtion_rate = new double[count_patch];
	_Xtion_rate_size = count_patch;
	ARRAY::copy(_Xtion_rate, array, _Xtion_rate_size);
	delete m;

	extinction_func_ptr = &LCE_Extinction::extinction_variable_total;
	return true;
}

// -----------------------------------------------------------------------------
// LCE_Extinction::set_survival_rate
// -----------------------------------------------------------------------------
/** the sruvival rate is the proportion/absolute number of individuals  which
 * survive if an extincation happens.
 * Default is when the value is 0, i.e. nobody survives. In this case the array
 * _survival_rate is not used.
 */
bool LCE_Extinction::set_survival_rate()
{
	// sex specific settings
	bool set = false;
	Param* paramF = this->get_parameter("extinction_rate_survival_fem");
	Param* paramM = this->get_parameter("extinction_rate_survival_mal");
	Param* param = this->get_parameter("extinction_rate_survival");
	if (paramF->isSet()) {
		set = set_survival_rate(paramF, FEM);

		if (_popPtr->get_sexInitRatio()) {
			if (paramM->isSet()) set = (set_survival_rate(paramM, MAL) || set);
			else error("Only one sex specific extinction_rate_survival is specified: both are required!\n");
		}
		else survivors_func_ptr[MAL] = &LCE_Extinction::survivors_all;
	}
	else if (_popPtr->get_sexInitRatio() && paramM->isSet()){
		error("Only one sex specific extinction_rate_survival is specified: both are required!\n");
	}

	// general settings
	else if (param->isSet()) {
		double sexRatio = _popPtr->get_sexInitRatio();
		if(sexRatio){            // female and male
			set = (set_survival_rate(param, MAL, sexRatio)
					|| set_survival_rate(param, MAL, sexRatio-(1e-6))); // bricolage: 1e-6 forces 0.5 to be ceiled instead of floored
		}
		else{
			survivors_func_ptr[MAL] = &LCE_Extinction::survivors_all;
			set = set_survival_rate(param, FEM);                   // hemaphrodites
		}
	}
	else{           // nothing is set
	   survivors_func_ptr[FEM] = &LCE_Extinction::survivors_none;
	   survivors_func_ptr[MAL] = &LCE_Extinction::survivors_none;
	}

	if(!set) return false;    // partial/sex specifiuc settings are not available


	// partial/sex specific extinction
	if(	extinction_func_ptr == &LCE_Extinction::extinction_variable_total){
		extinction_func_ptr = &LCE_Extinction::extinction_variable_partial;
	}
	else if(extinction_func_ptr == &LCE_Extinction::extinction_homogenous_high_total){
		extinction_func_ptr = &LCE_Extinction::extinction_homogenous_high_partial;
	}
	else{  // extinction_func_ptr == &LCE_Extinction::extinction_homogenous_low_total
		// if only  sex is available there is a function which is much more efficient
		if(this->_popPtr->get_sexInitRatio()) extinction_func_ptr = &LCE_Extinction::extinction_homogenous_low_partial_2sex;
		else                                  extinction_func_ptr = &LCE_Extinction::extinction_homogenous_low_partial_1sex;
	}

	return true;
}

// -----------------------------------------------------------------------------
// LCE_Extinction::set_survival_rate
// -----------------------------------------------------------------------------
/** the survival rate is the proportion/absolute number of individuals  which
 * survive if an extincation happens.
 * Default is when the value is 0, i.e. nobody survives. In this case the array
 * _survival_rate is not used.
 * function returnes false if nobody (at any time) is removed, i.e. survival rate is 1
 */
bool LCE_Extinction::set_survival_rate(Param* param, sex_t SEX, double correction)
 {
	if (!param->is_matrix()) { // if it is a single value or not used
		double val = param->get_value();
		if (val == my_NAN && !param->isTemporalParam()) { // there is never an extinction
			if (_survival_rate[SEX]) { delete[]_survival_rate[SEX]; _survival_rate[SEX] = NULL;}
			_survival_rate_size[SEX] = 0;
			survivors_func_ptr[SEX] = &LCE_Extinction::survivors_all;
			return false;     // newer an extinction will happen
		}
		if(val==0){         // individuals  are entirely removed
			survivors_func_ptr[SEX] = &LCE_Extinction::survivors_none;
			return true;
		}
		if (_survival_rate_size[SEX] != 1) {  // set the value (a single value)
			if (_survival_rate[SEX]) delete[]_survival_rate[SEX];
			_survival_rate[SEX] = new double[1];
			_survival_rate_size[SEX] = 1;
		}
		if(val >= 1){              // absolute number
			survivors_func_ptr[SEX] = &LCE_Extinction::survivors_absolute_const;
			_survival_rate[SEX][0] = my_round(val*correction);
		}
		else{                     // relative number
			survivors_func_ptr[SEX] = &LCE_Extinction::survivors_relative_const;
			_survival_rate[SEX][0] = val;
    }
	}
	else{                            // if it is a matrix
		TMatrix* m = param->get_matrix();
		double* array = m->get();
		unsigned int count_patch = m->length();

		// check the number of values per number of patches
		unsigned int nbPatch = _popPtr->get_nbPatch();
		if (count_patch > nbPatch) warning("There are more survivor rates specified than patches! Only a part of the survivor rates  is considered!\n");
		else if (nbPatch % count_patch) 	warning("The number of survivor rates is not a entire subset of the number of patches!\n");

		if (_survival_rate_size[SEX] != count_patch) {
			if (_survival_rate[SEX]) delete[]_survival_rate[SEX];
			_survival_rate[SEX] = new double[count_patch];
			_survival_rate_size[SEX] = count_patch;
		}
		bool absolute = false, relative = false;
		double val;
		for (unsigned int p = 0; p < count_patch; ++p) {
			_survival_rate[SEX][p] = val = array[p];  // get the current value
			if(val!=0){
				if(val<1) relative = true;                            // relative
				else if(val!=my_NAN){                                 // absolute
					_survival_rate[SEX][p] = my_round(val*correction);  // correction just for absolute values
					absolute = true;
        }
			}
		}
		delete m;

		if (absolute) {
			if (relative) survivors_func_ptr[SEX] = &LCE_Extinction::survivors_mixed;
			else          survivors_func_ptr[SEX] = &LCE_Extinction::survivors_absolute_var;
		}
		else if (relative) survivors_func_ptr[SEX] = &LCE_Extinction::survivors_relative_var;
		else if (param->isTemporalParam()) survivors_func_ptr[SEX] = &LCE_Extinction::survivors_all; // there is no extinction at the moment
		else {              // there is never an extinction
			if (_survival_rate[SEX]) {delete[]_survival_rate[SEX]; _survival_rate[SEX] = NULL;}
			_survival_rate_size[SEX] = 0;
			return false;
		}
	}

	// check if it is a coalescence simulation
	if(_popPtr->getCoalescence()){
		if(     survivors_func_ptr[SEX] == &LCE_Extinction::survivors_relative_var)   survivors_func_ptr[SEX] = &LCE_Extinction::survivors_relative_var_coal;   // survivors are only relativly defined
		else if(survivors_func_ptr[SEX] == &LCE_Extinction::survivors_absolute_var)   survivors_func_ptr[SEX] = &LCE_Extinction::survivors_absolute_var_coal; // survivors are only absolutivly defined
		else if(survivors_func_ptr[SEX] == &LCE_Extinction::survivors_relative_const) survivors_func_ptr[SEX] = &LCE_Extinction::survivors_relative_const_coal; // survivors are only relativly defined
		else if(survivors_func_ptr[SEX] == &LCE_Extinction::survivors_absolute_const) survivors_func_ptr[SEX] = &LCE_Extinction::survivors_absolute_const_coal; // survivors are only absolutivly defined
		else if(survivors_func_ptr[SEX] == &LCE_Extinction::survivors_mixed)          survivors_func_ptr[SEX] = &LCE_Extinction::survivors_mixed_coal;      // survivors are relativly and absolutivly defined
	}
	return true;
}

// -----------------------------------------------------------------------------
void LCE_Extinction::temporal_change(const unsigned int& gen)
{
	// temporal parameters
    if(_paramSet->getTemporalParams(gen) && _paramSet->updateTemporalParams(gen)){
        init(_popPtr);
	}
}

// -----------------------------------------------------------------------------
// LCE_Extinction::survivors
// -----------------------------------------------------------------------------
/** if extinction hits the patch what to do */
unsigned int         // survivors are only absolutely defined (identical among patches)
LCE_Extinction::survivors_absolute_const_coal(Patch* curPatch, const sex_t & SEX)
{
	unsigned int size = _survival_rate[SEX][0];
	assert(size==0 || size>=1);
	if(size >= curPatch->size(SEX, ADLTx)) return curPatch->size(SEX, ADLTx); // not enough individuals  to cut
	curPatch->set_size(SEX, ADLTx, size);
	return size;
}

unsigned int         // survivors are only absolutely defined (patch sepecific settings)
LCE_Extinction::survivors_absolute_var_coal(Patch* curPatch, const sex_t & SEX)
{
	unsigned int size = _survival_rate[SEX][curPatch->get_ID() % _survival_rate_size[SEX]];
	assert(size==0 || size>=1);
	if(size >= curPatch->size(SEX, ADLTx)) return curPatch->size(SEX, ADLTx); // not enough individuals  to cut
	curPatch->set_size(SEX, ADLTx, size);
	return size;
}

unsigned int         // survivors are only relativly defined (identical among patches)
LCE_Extinction::survivors_relative_const_coal(Patch* curPatch, const sex_t & SEX)
{
	assert(_survival_rate[SEX][0]>=0 && _survival_rate[SEX][0]<=1);
	unsigned int size = curPatch->size(SEX, ADLTx);
	size = my_round(size*_survival_rate[SEX][0]);
	curPatch->set_size(SEX, ADLTx, size);
	return size;
}

unsigned int         // survivors are only relativly defined (patch specific settings)
LCE_Extinction::survivors_relative_var_coal(Patch* curPatch, const sex_t & SEX)
{
	double value = _survival_rate[SEX][curPatch->get_ID() % _survival_rate_size[SEX]];
	assert(value>=0 && value<=1);
	unsigned int size = curPatch->size(SEX, ADLTx);
	size = my_round(size*value);
	curPatch->set_size(SEX, ADLTx, size);
	return size;
}

unsigned int        // survivors are specified relativly and absolutely defined (patch specific settings)
LCE_Extinction::survivors_mixed_coal(Patch* curPatch, const sex_t & SEX)
{
	double value = _survival_rate[SEX][curPatch->get_ID() % _survival_rate_size[SEX]];
	if(value>=1 || value==0) return survivors_absolute_var_coal(curPatch, SEX);
	return survivors_relative_var_coal(curPatch, SEX);
}


