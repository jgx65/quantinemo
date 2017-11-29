/** @file LCEdisperse.cpp
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

#include "lce_disperse.h"
using namespace std;

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                             ******** LCE_Disperse ********
// ----------------------------------------------------------------------------------------
LCE_Disperse::LCE_Disperse(int rank) :
LCE("dispersal", "dispersal", "", rank), _disp_model(0), _border_model(0),
_lattice_range(0), _disp_propagule_prob(0), _x_size(0), _y_size(0) {
    
	for (int i = 0; i < 2; ++i) {
		_dispMatrix[i] = NULL;
		_migr_rate[i] = 0;
		_tot_emigRate[i] = 0;
		_disp_factor[i] = 0;
		_disp_long_range_coef[i] = 0;
		_rel_abs_disp_rate[i] = 0;
	}
    
    add_parameter("dispersal_model", INT2, false, 0, 4, "0", false,
                  "Dispersal models:\n" \
                  "  0: migration-pool island model\n" \
                  "  1: propagule-pool island model\n" \
                  "  2: 1D stepping-stone model\n" \
                  "  3: 2D stepping-stone model\n", 1);
    
    add_parameter("dispersal_border_model", INT2, false, 0, 2, "0", false,
                  "Border models for stepping-stone models:\n"
                  "  0: circle/torus\n" \
                  "  1: reflecting boundaries\n" \
                  "  2: absorbing boundaries.", 1);
    
    add_parameter("dispersal_lattice_range", INT2, false, 0, 1, "0", false,
                  "Dispersal range for 2D stepping-stone model:\n" \
                  "  0: 4 neighbors (horizontal & vertical)\n" \
                  "  1: 8 neighbors (horizontal & vertical & diagonal).", 3);
    
	add_parameter("dispersal_lattice_dims", MAT, false, my_NAN, my_NAN,"", false,
                  "The dimension of the 2D steppings-stone area. " \
                  "Matrix with two values is expected: {nb_rows nb_cols}.", 3);
    
	add_parameter("dispersal_propagule_prob", DBL, false, 0, 1, "1",
                  true,
                  "Specifies the probability that an emigrant migrates to the " \
                  "propagule-assigned patch. Only used for propagule-island migration model.", 3);
    
    
    add_parameter("dispersal_rate", DBL_MAT, false, 0, my_NAN, "0", true,
                  "The emigration rate.", 0);
    
	add_parameter("dispersal_rate_fem", DBL_MAT, false, 0, my_NAN, "0",
                  true,
                  "The emigration rate of females.", 0);
    
	add_parameter("dispersal_rate_mal", DBL_MAT, false, 0, my_NAN, "0",
                  true,
                  "The emigration rate of males.", 0);
    
    
	// density dependent dispersal rate
    add_parameter("dispersal_rate_model", INT2, false, 0, 2, "0", false,
                  "What model we use for the dispersal rate:\n "
                  "  0: Flat rate\n" \
                  "  1: Rate depend on density with an effective rate of migration_rate when the density is one \n" \
                  "  2: Rate depend on density ans is smoothly adjusted so that after emigration, density is lower than one", 1);
    
    
    add_parameter("dispersal_rate_model_fem", INT2, false, 0, 1, "0", false,
                  "If we want a density dependent dispersal for femal (1 = yes)", 4);
  
        
    add_parameter("dispersal_rate_model_mal", INT2, false, 0, 1, "0", false,
                  "If we want a density dependent dispersal for male (1 = yes)", 4);
    
    
    
    add_parameter("dispersal_direction", INT2, false, 0, 1,  "0", false,
                  "Direction of dispersal:\n" \
                  "  0: emigration\n" \
                  "  1: immigration", 5);
    
}

// ----------------------------------------------------------------------------------------
LCE_Disperse::~LCE_Disperse() {
	for (unsigned int s = 0; s < 2; ++s) {
		if (_dispMatrix[s]) delete _dispMatrix[s];
		if (_disp_factor[s]) delete[] _disp_factor[s];
	}
}

// ----------------------------------------------------------------------------------------
// LCE_Disperse::setDispMatrix
// ----------------------------------------------------------------------------------------
/** returns false if no migration */
bool LCE_Disperse::setDispMatrix(TMatrix* mat) {
	bool used = checkDispMatrix(mat, FEM);
	_rel_abs_disp_rate[MAL] = _rel_abs_disp_rate[FEM];
    
	if (_dispMatrix[MAL]) delete _dispMatrix[MAL];
	if (_dispMatrix[FEM]) delete _dispMatrix[FEM];
    
	_dispMatrix[MAL] = new TMatrix(*mat);
	_dispMatrix[FEM] = mat;               // reuse the input matrix
	return used;
}

// ----------------------------------------------------------------------------------------
// LCE_Disperse
// ----------------------------------------------------------------------------------------
/** This function is set, when there is no migration or if there is only a single population.
 * All offsprings are simply transfered to the adult container
 */
void LCE_Disperse::migrate_zeroMigration() {
	// for each patch
	vector<Patch*>::iterator curPop = _popPtr->get_vFullPatch().begin();
	vector<Patch*>::iterator endPop = _popPtr->get_vFullPatch().end();
	for (; curPop != endPop; ++curPop) {
		(*curPop)->swap(OFFSx, ADLTx); // "move" would be correct, but swap is much faster and also ok for here
	}
}

// ----------------------------------------------------------------------------------------
// LCE_Disperse
// ----------------------------------------------------------------------------------------
/** This function is set, when there is no migration or if there is only a single population.
 * All offsprings are simply transfered to the adult container
 */
void LCE_DisperseCoalescence::migrate_zeroMigration() {
	// for each patch
	vector<Patch*>::iterator curPop = _popPtr->get_vFullPatch().begin();
	vector<Patch*>::iterator endPop = _popPtr->get_vFullPatch().end();
	for (; curPop != endPop; ++curPop) {
		(*curPop)->swap_coal(FEM, OFFSx, ADLTx); // "move" would be correct, but swap is much faster and also ok for here
	}
}

// ----------------------------------------------------------------------------------------
// LCE_Disperse::setDispMatrix
// ----------------------------------------------------------------------------------------
/** returns false if no migration */
bool LCE_Disperse::setDispMatrix(sex_t sex, TMatrix* mat) {
	bool used = checkDispMatrix(mat, sex);
	if (_dispMatrix[sex]) delete _dispMatrix[sex];
	_dispMatrix[sex] = mat;
	return used;
}

// ----------------------------------------------------------------------------------------
// LCE_Disperse::checkDispMatrix
// ----------------------------------------------------------------------------------------
/** checks the size of the matrix and the sum of migr rates and returns the proportion of rersidents */
bool LCE_Disperse::checkDispMatrix(TMatrix* mat, sex_t SEX) {
	assert(mat);
    
	double sum, val;
	unsigned int i, j;
	mat->get_dims(_x_size, _y_size);
    
	if (_x_size != _y_size || _x_size != _nb_patch) {
		error(
              "The dimension of the dispersal matrix (%ix%i) does not match the number of patches (%i)!\n",
              _x_size, _y_size, _nb_patch);
	}
    
	bool absolute = false;         // are absolute and/or relative numbers used?
	bool relative = false;
	for (i = 0; i < _x_size; ++i) {
		sum = 0;                      // sum of emigration and immigration rates
		for (j = 0; j < _y_size; ++j) {
			val = mat->get(i, j);
			if (i == j) {	// residents
				if(val<0) error("Parameter dispersal_rate cannot have negative numbers in the diagonal (residents)!");
				continue;
			}
			if (val > 1) {
				absolute = true;
				sum += mat->get(i, j);
			}  // [1;-[    absolute
			else if (val > 0) {
				relative = true;
				sum += mat->get(i, j);
			}  // ]0;1[    relative
			else if (val <= -1) {
				absolute = true;
				sum -= mat->get(i, j);
			}  // ]-;-1]   absolute
			else if (val < 0) {
				relative = true;
				sum -= mat->get(i, j);
			}  // ]-1;0[   relative
		}
        
		if (absolute)
			mat->set(i, i, sum);  // set the number of EMIGRANTS to the diagonal
		else if (relative) {
			if (abs(1 - sum - mat->get(i, i)) > 1e-4) { // (sum != 1): floating numbers can never be the same
				warning(
						"The elements of row %i of the dispersal matrix do not sum up to 1! Residents fraction adjusted\n",
						i + 1);
			}
			mat->set(i, i, sum); // set the proportion of EMIGRANTS to the diagonal
		}
		else mat->set(i, i, 1); // nobody emigrants
	}
    
	if (absolute && relative)
		error(
              "Parameter dispersal_rate cannot have relative AND absolute numbers in a matrix!");
	if (relative)
		_rel_abs_disp_rate[SEX] = 0;   // relative number of dispersers
	else if (absolute)
		_rel_abs_disp_rate[SEX] = 1;   // absolute number of dispersers
	else return false;                               // no migration
	return true;
}

// ----------------------------------------------------------------------------------------
// LCE_Disperse
// ----------------------------------------------------------------------------------------
void LCE_Disperse::execute() {
#ifdef _DEBUG
	message("  LCE_Disperse ... ");
#endif
    
	preDispersal();
	(this->*migration_func_ptr)();
	(this->*func_ptr_postDispersal)();
    
#ifdef _DEBUG
	message("done! (Patches: %i; Individuals(F/M): [%i/%i; %i/%i])\n",
			_popPtr->get_nbFullPatch(),
			_popPtr->size(FEM, OFFSx), _popPtr->size(MAL, OFFSx),
			_popPtr->size(FEM, ADLTx), _popPtr->size(MAL, ADLTx));
    assert(_popPtr->individual_container_ok());
#endif
}

// ----------------------------------------------------------------------------------------
// LCE_Disperse
// ----------------------------------------------------------------------------------------
void LCE_DisperseCoalescence::execute() {
#ifdef _DEBUG
	message("  LCE_DisperseCoalescence ... ");
#endif
    
	//(this->*func_ptr_preDispersal)();
	(this->*migration_func_ptr)();
	(this->*func_ptr_postDispersal)();
    
#ifdef _DEBUG
	message("done! (Patches: %i; Individuals(F/M): [%i/%i; %i/%i])\n",
			_popPtr->get_nbFullPatch(),
			_popPtr->size(FEM, OFFSx), _popPtr->size(MAL, OFFSx),
			_popPtr->size(FEM, ADLTx), _popPtr->size(MAL, ADLTx));
#endif
}

// ----------------------------------------------------------------------------------------
// LCE_Disperse::preDispersal
// ----------------------------------------------------------------------------------------
/** prepares the dispersal, resets counters,... */
void LCE_Disperse::preDispersal() {
	assert(!_popPtr->size(ADLTx));           // adult containers should be empty
	Patch* curPatch;
	vector<Patch*>::iterator curPop = _popPtr->get_vFullPatch().begin();
	vector<Patch*>::iterator endPop = _popPtr->get_vFullPatch().end();
	for (; curPop != endPop; ++curPop) {
		curPatch = (*curPop);
		curPatch->nbImmigrant = 0;              // reset counters
		curPatch->nbEmigrant = 0;              // reset counters
		curPatch->nbKolonisers = my_NAN;
	}
}

// ----------------------------------------------------------------------------------------
// LCE_Disperse::postDispersal
// ----------------------------------------------------------------------------------------
/** set the demographic information about colonization and migration of the current patch
 * nbEmigrant:  is not incremented here but in the function _sendEmigrants
 * nbImmigrant: is not incremented here but in the function _sendEmigrants
 * if a patch was newly colonized it was already added to the container _vFullPatch
 * if a patch is newly deserted it has to be removed HERE from the cotainer _vFullPatch !!!
 */
void LCE_Disperse::postDispersal_noSample_withFull()
{
	// for all patches which were populated check if they are now empty
	vector<Patch*>::iterator curPop, endPop = _popPtr->get_vFullPatch().end();
	for (curPop = _popPtr->get_vFullPatch().begin(); curPop != endPop;) { // for each populated patch
		if ((*curPop)->size(ALL)) ++curPop;
		else {                                                 // patch is empty
			assert((*curPop)->nbEmigrant);
			_popPtr->new_emptyPatch(curPop, endPop); // remove the patch from the populated patches
		}
	}
    
	// for all NEWLY colonized patches
	endPop = _popPtr->get_vTempPatch().end();
	for (curPop = _popPtr->get_vTempPatch().begin(); curPop != endPop;
         ++curPop) {  // for each populated patch
		assert(
               (*curPop)->size(ADLTx) && (*curPop)->size(ADLTx)==(*curPop)->nbImmigrant && !(*curPop)->nbEmigrant);
		(*curPop)->nbKolonisers = (*curPop)->nbImmigrant;
	}
    
	// merge the two containers
	_popPtr->add_tempPatch_noSample_withFull();
    
}

// ----------------------------------------------------------------------------------------
// LCE_Disperse::postDispersal
// ----------------------------------------------------------------------------------------
/** set the demographic information about colonization and migration of the current patch
 * nbEmigrant:  is not incremented here but in the function _sendEmigrants
 * nbImmigrant: is not incremented here but in the function _sendEmigrants
 * if a patch was newly colonized it was already added to the container _vFullPatch
 * if a patch is newly deserted it has to be removed HERE from the cotainer _vFullPatch !!!
 */
void LCE_DisperseCoalescence::postDispersal_noSample_withFull() {
	// for all patches which were populated check if they are now empty
	vector<Patch*>::iterator curPop, endPop = _popPtr->get_vFullPatch().end();
	for (curPop = _popPtr->get_vFullPatch().begin(); curPop != endPop;) { // for each populated patch
		if ((*curPop)->size(ALL)) ++curPop;
		else _popPtr->new_emptyPatch(curPop, endPop); // patch is empty: remove it
	}
    
	// merge the two containers
	_popPtr->add_tempPatch();
    
}

// ----------------------------------------------------------------------------------------
// LCE_Disperse::postDispersal
// ----------------------------------------------------------------------------------------
/** set the demographic information about colonization and migration of the current patch
 * nbEmigrant:  is not incremented here but in the function _sendEmigrants
 * nbImmigrant: is not incremented here but in the function _sendEmigrants
 * if a patch was newly colonized it was already added to the container _vFullPatch
 * if a patch is newly deserted it has to be removed HERE from the cotainer _vFullPatch !!!
 */
void LCE_Disperse::postDispersal_withSample_withFull() {
	// for all patches which were populated check if they are now empty
	vector<Patch*>::iterator curPop, endPop = _popPtr->get_vFullPatch().end();
	for (curPop = _popPtr->get_vFullPatch().begin(); curPop != endPop;) { // for each populated patch
		if ((*curPop)->size(ALL)) ++curPop;
		else {                                                 // patch is empty
			assert((*curPop)->nbEmigrant);
			_popPtr->new_emptyPatch(curPop, endPop);   // curPop and endPop are adjusted
		}
	}
    
	// for all NEWLY colonized patches
	endPop = _popPtr->get_vTempPatch().end();
	for (curPop = _popPtr->get_vTempPatch().begin(); curPop != endPop;
         ++curPop) {  // for each populated patch
		assert(
               (*curPop)->size(ADLTx) && (*curPop)->size(ADLTx)==(*curPop)->nbImmigrant && !(*curPop)->nbEmigrant);
		(*curPop)->nbKolonisers = (*curPop)->nbImmigrant;
	}
    
	// merge the two containers
	_popPtr->add_tempPatch();
    
}

// ----------------------------------------------------------------------------------------
// LCE_Disperse::postDispersal
// ----------------------------------------------------------------------------------------
/** set the demographic information about colonization and migration of the current patch
 * nbEmigrant:  is not incremented here but in the function _sendEmigrants
 * nbImmigrant: is not incremented here but in the function _sendEmigrants
 * if a patch was newly colonized it was already added to the container _vFullPatch
 * if a patch is newly deserted it has to be removed HERE from the cotainer _vFullPatch !!!
 */
void LCE_Disperse::postDispersal_withSample_noFull() {
	// for all patches which were populated check if they are now empty
	vector<Patch*>::iterator curPop, endPop = _popPtr->get_vFullPatch().end();
	for (curPop = _popPtr->get_vFullPatch().begin(); curPop != endPop;) { // for each populated patch
		if ((*curPop)->size(ALL)) ++curPop;                             // populated
		else _popPtr->new_emptyPatch(curPop, endPop);   // curPop and endPop are adjusted
	}
}

// ----------------------------------------------------------------------------------------
// LCE_Disperse::postDispersal
// ----------------------------------------------------------------------------------------
/** set the demographic information about colonization and migration of the current patch
 * nbEmigrant:  is not incremented here but in the function _sendEmigrants
 * nbImmigrant: is not incremented here but in the function _sendEmigrants
 * if a patch was newly colonized it was already added to the container _vFullPatch
 * if a patch is newly deserted it has to be removed HERE from the cotainer _vFullPatch !!!
 */
void LCE_Disperse::postDispersal_noSample_noFull() {
	// for all patches which were populated check if they are now empty
	Patch* curPatch;
	vector<Patch*>::iterator curPop, endPop = _popPtr->get_vPatch().end();
	for (curPop = _popPtr->get_vPatch().begin(); curPop != endPop; ++curPop) { // for each patch
		curPatch = *curPop;
		if (!curPatch->size(ALL)) {  // if the patch is empty
			if (curPatch->nbEmigrant) curPatch->set_isExtinct(true); // newly empty
		}
		else if (!curPatch->nbEmigrant
                 && curPatch->size(ADLTx) == curPatch->nbImmigrant) { // newly colonized
			curPatch->nbKolonisers = curPatch->nbImmigrant;
			assert(!curPatch->get_isExtinct());         // has already been done
		}
	}
}

// ----------------------------------------------------------------------------------------
// LCE_Disperse::postDispersal
// ----------------------------------------------------------------------------------------
/** set the demographic information about colonization and migration of the current patch
 * nbEmigrant:  is not incremented here but in the function _sendEmigrants
 * nbImmigrant: is not incremented here but in the function _sendEmigrants
 * if a patch was newly colonized it was already added to the container _vFullPatch
 * if a patch is newly deserted it has to be removed HERE from the cotainer _vFullPatch !!!
 */
void LCE_DisperseCoalescence::postDispersal_noSample_noFull() {
	// for all patches which were populated check if they are now empty
	vector<Patch*>::iterator curPop, endPop = _popPtr->get_vPatch().end();
	for (curPop = _popPtr->get_vPatch().begin(); curPop != endPop; ++curPop) { // for each patch
		assert(((*curPop)->size(ALL) && !(*curPop)->get_isExtinct()) // this has already been done
               ||(!(*curPop)->size(ALL)));// needed that all is within assert()
		(*curPop)->set_isExtinct(!(*curPop)->size(ALL)); // if the patch is empty
	}
}

// ----------------------------------------------------------------------------------------
// LCE_Disperse
// ----------------------------------------------------------------------------------------
/** if the migration rates are defined by a matrix */
void LCE_Disperse::migrate_matrix() {
	Patch *curPatch;
	unsigned int nbInd[2];        // nbInd is set by _computeTotEmigrants()
	double factor[2], sum_m[2], migr[2]; // factor is set by get_migr_factor() & sum_m is set by _computeTotEmigrants()
	unsigned int home, target;
    
	// for each populated patch
	vector<Patch*>::iterator curPop = _popPtr->get_vFullPatch().begin();
	vector<Patch*>::iterator endPop = _popPtr->get_vFullPatch().end();
	for (; curPop != endPop; ++curPop) {
		curPatch = *curPop;
		home = curPatch->get_ID();
		_tot_emigRate[FEM] = _dispMatrix[FEM]->get(home, home);
		_tot_emigRate[MAL] = _dispMatrix[MAL]->get(home, home);
		get_migr_factor(curPatch, factor);
        
		if (_computeTotEmigrants(curPatch, nbInd, _tot_emigRate, sum_m,
                                 factor)) { // total number of emigrants
			// migration to each other patch (also empty ones)
			for (target = 0; target < _popPtr->get_nbPatch(); ++target) {
				if (home == target) continue; // if it is the same patch continue
				migr[FEM] = _dispMatrix[FEM]->get(home, target);
				migr[MAL] = _dispMatrix[MAL]->get(home, target);
				if (!_sendEmigrants(curPatch, home, target, nbInd, migr, sum_m,
                                    factor)) break;   // send the emigrants
			}
		} // end for each target
        
		curPatch->move(OFFSx, ADLTx); // all not migrated individuals remain in the original patch
	}                   //end for each home patch
}

// ----------------------------------------------------------------------------------------
// LCE_Disperse
// ----------------------------------------------------------------------------------------
/** if the migration rates are defined by a matrix */
void LCE_DisperseCoalescence::migrate_matrix() {
	Patch *curPatch;
	unsigned int nbEmigr;        // nbInd is set by _computeTotEmigrants()
	double factor, sum_m; // factor is set by get_migr_factor() & sum_m is set by _computeTotEmigrants()
	unsigned int home, target;
    
	// for each populated patch
	vector<Patch*>::iterator curPop = _popPtr->get_vFullPatch().begin();
	vector<Patch*>::iterator endPop = _popPtr->get_vFullPatch().end();
	for (; curPop != endPop; ++curPop) {
		curPatch = *curPop;
		home = curPatch->get_ID();
		_tot_emigRate[FEM] = _dispMatrix[FEM]->get(home, home);
		factor = (this->*get_migr_factor_funcPtr[FEM])(curPatch, FEM); // adjust factor
        
		if (_computeTotEmigrants(curPatch, nbEmigr, _tot_emigRate[FEM], sum_m,
                                 factor)) {   // total number of emigrants
			// migration to each other patch (also empty ones)
			for (target = 0; target < _nb_patch; ++target) {
				if (home == target) continue; // if it is the same patch continue
				_migr_rate[FEM] = _dispMatrix[FEM]->get(home, target);
				if (!_sendEmigrants(curPatch, home, target, nbEmigr,
                                    _migr_rate[FEM], sum_m, factor)) break;
			}
		} // end for each target
        
		//curPatch->move(OFFSx, ADLTx);                   // all not migrated individuals remain in the original patch
	} //end for each home patch
}

// ----------------------------------------------------------------------------------------
/** this function computes the total number of emigrants
 * returns true if there are emigrants and false if not
 * the following parameters are set by this function:
 *   - nbEmigrants
 *   - sum_m
 */
bool
LCE_Disperse::_computeTotEmigrants(Patch* curPatch, unsigned int* nbEmigrants,
                                   double* migrTotRate, double* sum_m, double* factor)
{
	nbEmigrants[FEM] = _computeTotEmigrants(curPatch, migrTotRate[FEM], sum_m[FEM], factor[FEM], FEM);
	nbEmigrants[MAL] = _computeTotEmigrants(curPatch, migrTotRate[MAL], sum_m[MAL], factor[MAL], MAL);
	return (nbEmigrants[FEM] || nbEmigrants[MAL]);
}

// ----------------------------------------------------------------------------------------
/** this function computes and returns the total number of emigrants
 * sum_m is adjusted/set
 */
unsigned int
LCE_Disperse::_computeTotEmigrants(Patch* curPatch, const double& migrTotRate,
                                   double& sum_m, const double& factor, const sex_t& SEX)
{
	unsigned int popSize = curPatch->size(SEX, OFFSx); // get current population size
	if (!popSize) {               // pop size is 0: no emigrants
		return 0;
	}
    
	double m = migrTotRate * factor;            // compute total emigration rate
    if(_paramSet->getValue("dispersal_rate_model")==2) m = factor;   
	sum_m = m;                        // this is the sum of all emigration rates
	if (!m) return 0;
    
	// compute number of emigrants
	switch (_rel_abs_disp_rate[SEX]) {
		case 0: return (unsigned int) get_pop_ptr()->rand().Binomial(m, popSize); // emigration rate defined
		case 1: return ((unsigned int) m < popSize) ? (unsigned int) m : popSize; // number of emigrants defined
	}
	return popSize; // should never happen, but it avoids a compilation message
}

// ----------------------------------------------------------------------------------------
/** this function computes and returns the total number of emigrants
 * sum_m is adjusted/set
 */
bool LCE_DisperseCoalescence::_computeTotEmigrants(Patch* curPatch,
                                                   unsigned int& totEmigr, const double& migrTotRate, double& sum_m,
                                                   const double& factor) {
	unsigned int popSize = curPatch->size(FEM, OFFSx); // get current population size
	if (!popSize) {               // pop size is 0: no emigrants
		totEmigr = 0;
		return 0;
	}
	double m = migrTotRate * factor;            // compute total emigration rate
 
    if(_paramSet->getValue("dispersal_rate_model")==2) m = factor;
	
	sum_m = m;                        // this is the sum of all emigration rates
	if (!m) return 0;
    
	// compute number of emigrants
	switch (_rel_abs_disp_rate[FEM]) {
		case 0:
			totEmigr = (unsigned int) get_pop_ptr()->rand().Binomial(m, popSize);
			break; // emigration rate defined
		case 1:
			totEmigr =
            ((unsigned int) m < popSize) ? (unsigned int) m : popSize;
			break; // number of emigrants defined
	}
    
	// the ones which do not emigrate become adult in this patch (only in coalescence possible like that)
	curPatch->add_size(FEM, ADLTx, popSize - totEmigr);
	curPatch->set_size(FEM, OFFSx, 0);                // rest the juveniles to 0
    
	return totEmigr;
}

// ----------------------------------------------------------------------------------------
/** function returns the number of emigrants for the given migration rate.
 * if m>=1 then the migration rate is the absolute number of emigrants!
 * if m<1  then m is the migration rate
 * if m is a migration rate, then sum_m has to be adjusted as in this case a multinomial distribution is used
 * totEmigr is adjusted
 */
unsigned int LCE_Disperse::get_nbEmigrant(const double& m, double& sum_m,
                                          unsigned int& totEmigr) {
	unsigned int nbEmigr;
	if (m >= 1.0)
		nbEmigr = ((unsigned int) m < totEmigr) ? (unsigned int) m : totEmigr; // number of emigrants defined
	else {                                            // emigration rate defined
		nbEmigr = get_pop_ptr()->rand().Binomial(m / sum_m, totEmigr); // adjust m (multinomial distribution)
		sum_m -= m;                      // adapt the sum of the migration rates
	}
    
	totEmigr -= nbEmigr;                                 // adjust the totEmigr!
	return nbEmigr;
}

// ----------------------------------------------------------------------------------------
/** function to send emigrants from the home to the target patch for both sexes
 * if the migration rate is negative, the emigrating individuals are removed (absorbed)
 * returns true if the patch is still populated and false if empty
 */
bool LCE_Disperse::_sendEmigrants(Patch* curPatch, const unsigned int& home,
                                  const unsigned int& target, unsigned int* totEmigr,
                                  double* migrRate, double* sum_m, double* factor)
{
	assert(curPatch->get_ID()==home);
	_sendEmigrants(curPatch, home, target, totEmigr[FEM], migrRate[FEM],
                   sum_m[FEM], factor[FEM], FEM);
	_sendEmigrants(curPatch, home, target, totEmigr[MAL], migrRate[MAL],
                   sum_m[MAL], factor[FEM], MAL);
	return (totEmigr[FEM] || totEmigr[MAL]);
}

// ----------------------------------------------------------------------------------------
/** function to send emigrants form the home to the target patch for a specific sex
 * if the migration rate is negative, the emigrating individuals are removed (absorbed)
 * parameters which are changed:
 *     - totEmigr
 *     - sum_m
 * the number of not yet send emigrants is returned
 */
unsigned int LCE_Disperse::_sendEmigrants(Patch* curPatch,
                                          const unsigned int& home, const unsigned int& target,
                                          unsigned int& totEmigr, const double& migrRate, double& sum_m,
                                          const double& factor, const sex_t& SEX)
{
	if (!totEmigr) return totEmigr;
	assert(sum_m > 0);
    
	// compute the corrected migration rate (multinomial distribution)
	double m = migrRate * factor;         // compute current migration rate
    if(_paramSet->getValue("dispersal_rate_model")==2) m = factor;	
	
	if (!m) return totEmigr;               // check if migration occurs
    
	// perform migration
	unsigned int nb_emigrants;
	if (m > 0) {           // normal migration
		nb_emigrants = get_nbEmigrant(m, sum_m, totEmigr);
		if (!nb_emigrants) return totEmigr;             // no emigrants : return
		_sendEmigrant(curPatch, home, target, nb_emigrants, SEX,
                      curPatch->size(SEX, OFFSx));
	}
	else {              // emigrants are removed (absorbing boundaries)
		m *= -1;                                             // make m possitive
		nb_emigrants = get_nbEmigrant(m, sum_m, totEmigr);
		if (!nb_emigrants) return totEmigr;             // no emigrants : return
		_absorbEmigrant(curPatch, nb_emigrants, SEX,
                        curPatch->size(SEX, OFFSx));  // absorbe the individuals
	}
    
	return totEmigr;                            // remaining number of emigrants
}

// ----------------------------------------------------------------------------------------
/** function to send emigrants form the home to the target patch for a specific sex
 * if the migration rate is negative, the emigrating individuals are removed (absorbed)
 */
unsigned int LCE_DisperseCoalescence::_sendEmigrants(Patch* curPatch,
                                                     const unsigned int& home, const unsigned int& target,
                                                     unsigned int& totEmigr, const double& migrRate, double& sum_m,
                                                     const double& factor) {
	if (!totEmigr) return totEmigr;
	assert(curPatch->get_ID()==home);
	assert(sum_m > 0);
    
	// compute the corrected migration rate (multinomial distribution)
	double m = migrRate * factor;         // compute current migration rate
	if(_paramSet->getValue("dispersal_rate_model")==2) m = factor;	
	if (!m) return totEmigr;               // check if migration occurs
    
	// perform migration
	unsigned int nb_emigrants;
	if (m > 0) {           // normal migration
		nb_emigrants = get_nbEmigrant(m, sum_m, totEmigr);
		if (!nb_emigrants) return totEmigr;             // no emigrants : return
		_sendEmigrant(curPatch, home, target, nb_emigrants); // send the emigrants
	}
	else {              // emigrants are removed (absorbing boundaries)
		m *= -1;                                             // make m possitive
		nb_emigrants = get_nbEmigrant(m, sum_m, totEmigr);
		if (!nb_emigrants) return totEmigr;             // no emigrants : return
	}
	return totEmigr;                            // remaining number of emigrants
}

// ----------------------------------------------------------------------------------------
/** function to send a single emigrant from the home to the target patch for a specific sex
 * the pop size is computed here
 */
void LCE_Disperse::_sendEmigrant(Patch* curPatch, const unsigned int& home,
                                 const unsigned int& target, const unsigned int& nbEmigr,
                                 const sex_t& SEX) {
	_sendEmigrant(curPatch, home, target, nbEmigr, SEX,
                  curPatch->size(SEX, OFFSx));
}

// ----------------------------------------------------------------------------------------
/** function to send a single emigrant from the home to the target patch for a specific sex
 * the pop size is given here
 */
void LCE_Disperse::_sendEmigrant(Patch* curPatch, const unsigned int& home,
                                 const unsigned int& target, const unsigned int& nbEmigr,
                                 const sex_t& SEX, const unsigned int& popSize) {
	assert(popSize == curPatch->size(SEX, OFFSx));
	curPatch->nbEmigrant += nbEmigr;        // increment the number of emigrants
	Patch* targetPatch = _popPtr->get_vPatch(target);
	targetPatch->nbImmigrant += nbEmigr; // increment the number of immigrants of the target patch
	for (unsigned int i = 0; i < nbEmigr; ++i) {
		_popPtr->move(SEX, OFFSx, home, ADLTx, target, get_pop_ptr()->rand().Uniform(popSize - i));
	}
	if (targetPatch->get_isExtinct()) _popPtr->new_fullPatch(targetPatch); // if the patch is newly colonized
}

// ----------------------------------------------------------------------------------------
/** function to send a single emigrant from the home to the target patch for a specific sex
 * the pop size is given here
 */
void LCE_DisperseCoalescence::_sendEmigrant(Patch* curPatch,
                                            const unsigned int& home, const unsigned int& target,
                                            const unsigned int& nbEmigr) {
	Patch* targetPatch = _popPtr->get_vPatch(target);
	targetPatch->addImmigrant(curPatch, nbEmigr);      // that are the emigrants
	if (targetPatch->get_isExtinct()) _popPtr->new_fullPatch(targetPatch); // if the patch is newly colonized
}

// ----------------------------------------------------------------------------------------
/** function to send a single emigrant from the home to the target patch for a specific sex
 * the pop size is given here
 * popSize is not used, but needed for the function pointers
 */
void LCE_DisperseCoalescence::_sendEmigrant(Patch* curPatch,
                                            const unsigned int& home, const unsigned int& target,
                                            const unsigned int& nbEmigr, const sex_t& SEX,
                                            const unsigned int& popSize) {
	assert(SEX==FEM);
	_sendEmigrant(curPatch, home, target, nbEmigr, SEX, popSize);
}

// ----------------------------------------------------------------------------------------
/** function to absorbe a single emigrant (emigration outside of the world) for a specific sex
 * the pop size is computed here
 */
void LCE_Disperse::_absorbEmigrant(Patch* curPatch, const unsigned int& nbEmigr,
                                   const sex_t& SEX) {
	_absorbEmigrant(curPatch, nbEmigr, SEX, curPatch->size(SEX, OFFSx));
}

// ----------------------------------------------------------------------------------------
/** function to absorbe a single emigrant (emigration outside of the world) for a specific sex
 * the pop size is given here
 */
void LCE_Disperse::_absorbEmigrant(Patch* curPatch, const unsigned int& nbEmigr,
                                   const sex_t& SEX, const unsigned int& popSize) {
	assert(popSize == curPatch->size(SEX, OFFSx));
	curPatch->nbEmigrant += nbEmigr;        // increment the number of emigrants
	for (unsigned int i = 0; i < nbEmigr; ++i) {
		curPatch->recycle(SEX, OFFSx, get_pop_ptr()->rand().Uniform(popSize - i));
	}
}

// ----------------------------------------------------------------------------------------
// LCE_Disperse
// ----------------------------------------------------------------------------------------
/** island migration model with a female and a male migration rate
 * migration rate to any patch is m/(_nb_patch-1)
 */
void LCE_Disperse::migrate_island() {
	Patch *curPatch;
	unsigned int home, target;
	unsigned int nbInd[2];          // nbInd is set by _computeTotEmigrants()
    double factor[2];               // factor is set by get_migr_factor()
    double sum_m[2];                // sum_m is set by _computeTotEmigrants()
    
	// for each populated patch
	vector<Patch*>::iterator curPop, endPop = _popPtr->get_vFullPatch().end();
	for (curPop = _popPtr->get_vFullPatch().begin(); curPop != endPop; ++curPop) {
		curPatch = *curPop;
		home = curPatch->get_ID();
		get_migr_factor(curPatch, factor);   // adjust factor
        
		if (_computeTotEmigrants(curPatch, nbInd, _tot_emigRate, sum_m, factor)) { // total number of emigrants
			// migration to each other patch
			for (target = 0; target < _nb_patch; ++target) {
				if (home == target) continue; // if it is the same patch continue
				if (!_sendEmigrants(curPatch, home, target, nbInd, _migr_rate,
                                    sum_m, factor)) break;
			}
		}
		curPatch->move(OFFSx, ADLTx); // all not migrated individuals remain in the original patch
	}           //end for each home patch
}

// ----------------------------------------------------------------------------------------
/** immigration island model with a female and a male migration rate
 * migration rate from any patch is m/(_nb_patch-1)
 * migration with replacement!!!
 */
void LCE_Disperse::immigrate_island()
{
    unsigned int nbInd[2];        // nbInd is set by _computeTotEmigrants()
    double factor[2], sum_m[2]; // factor is set by get_migr_factor() & sum_m is set by _computeTotEmigrants()
    
    // get the metapopulation size for each sex
    unsigned int neighbourSizes[2];
    neighbourSizes[FEM] = neighbourSizes[MAL] = 0;
    vector<Patch*>::iterator curPop, endPop = _popPtr->get_vFullPatch().end();
    for (curPop = _popPtr->get_vFullPatch().begin(); curPop != endPop; ++curPop) {
        neighbourSizes[FEM] += (*curPop)->size(FEM, OFFSx);
        neighbourSizes[MAL] += (*curPop)->size(MAL, OFFSx);
    }
    
    // perform migration
    for (curPop = _popPtr->get_vFullPatch().begin(); curPop != endPop; ++curPop) {
        get_migr_factor(*curPop, factor);   // adjust factor
        _computeTotEmigrants(*curPop, nbInd, _tot_emigRate, sum_m, factor);            // total number of emigrants (or immigrants in this case)
        assert(nbInd[FEM]<=(*curPop)->size(FEM, OFFSx));
        
        immigrate(*curPop, _popPtr->get_vFullPatch(), neighbourSizes[FEM], nbInd[FEM], FEM);
        immigrate(*curPop, _popPtr->get_vFullPatch(), neighbourSizes[MAL], nbInd[MAL], MAL);
    }
    
    // empty the offspring containers (migration with replacement!)
    _popPtr->flush(OFFSx);
}

// ----------------------------------------------------------------------------------------
/** Immigration to the given patch:
 * - homePop: the patch in focus
 * - vNeighbours: vector of all neighbouring patches
 * - neihgbourSizes: vector of population sizes of the given sex in same order as vNeighbours
 * - totSizes: sum of neighbourSize (size of all neighbours AND current patch of the given sex)
 * - nbMigr: total number of immigrants
 * - SEX: the sex in question
 */
void LCE_Disperse::immigrate(Patch* homePop, vector<Patch*> vNeighbours, unsigned int totSize, unsigned int nbMigr, sex_t SEX, age_idx fromAge, age_idx toAge)
{
    unsigned int homeSize = homePop->size(SEX, fromAge);  // population size of focal patch
    unsigned int curMigr, curSize;
    unsigned int neighbourSize = totSize - homeSize;      // population size of all neighbours
    double cur_m;
    
    //cout << "\n\nPopulation " << homePop->get_ID()+1 << "\tsize: " << homeSize << "\tNbMigr: " << nbMigr << "\tneighbourSize: " << neighbourSize << endl;
    
    // if immigration is not possible
    if(!neighbourSize && nbMigr){
        warning("Immigration: no neighour populations => immigration is not possbile => skipped!\n");
        nbMigr = 0;
    }
    
    // resident ones
    homePop->nbImmigrant += nbMigr;                     // increment the number of immigrants
    if(homeSize) _popPtr->copyMove_random_withReplacement(SEX, fromAge, homePop, toAge, homePop, homeSize-nbMigr);
    
    // for each neighbour
    vector<Patch*>::iterator curPop, endPop;
    for (curPop = vNeighbours.begin(), endPop = vNeighbours.end(); nbMigr && curPop != endPop; ++curPop) {
        if(*curPop==homePop) continue; // local patch, continue
        
        // multinomial random deviates
        curSize = (*curPop)->size(SEX, fromAge);        // get population size
        cur_m = (double)curSize/neighbourSize;          // compute migration rate
        neighbourSize -= curSize;                       // correct the remianing neighbour size
        curMigr = get_pop_ptr()->rand().Binomial(cur_m, nbMigr); // draw randomly the number of migrants
        if(!curMigr) continue;
        
        (*curPop)->nbEmigrant += curMigr; // increment the number of immigrants of the target patch
        //cout << "\tneighbour " << (*curPop)->get_ID()+1 << "\tsize: " << curSize << "\tnbMigr: " << nbMigr << "\tneighbourSize: " << neighbourSize+curSize << "\tcur_m: " << cur_m << "\tcurMigr: " << curMigr << endl;
        _popPtr->copyMove_random_withReplacement(SEX, fromAge, *curPop, toAge, homePop, curMigr);
        
        nbMigr -= curMigr;
    }
    
    assert(!nbMigr);        // all individuals have to have migrated!!!
}

// ----------------------------------------------------------------------------------------
/** island migration model with a female and a male migration rate
 * migration rate to any patch is m/(_nb_patch-1)
 */
void LCE_DisperseCoalescence::migrate_island() {
	Patch *curPatch;
	unsigned int home, target;
	unsigned int nbEmigr;   // nbEmigr is set by _computeTotEmigrants()
	double factor, sum_m; // factor is set by get_migr_factor() & sum_m is set by _computeTotEmigrants()
    
	// for each populated patch
	vector<Patch*>::iterator curPop, endPop = _popPtr->get_vFullPatch().end();
	for (curPop = _popPtr->get_vFullPatch().begin(); curPop != endPop;
         ++curPop) {
		curPatch = *curPop;
		home = curPatch->get_ID();
		factor = (this->*get_migr_factor_funcPtr[FEM])(curPatch, FEM); // adjust factor
        
		if (_computeTotEmigrants(curPatch, nbEmigr, _tot_emigRate[FEM], sum_m,
                                 factor)) {   // total number of emigrants
			// migration to each other patch
			for (target = 0; target < _nb_patch; ++target) {
				if (home == target) continue; // if it is the same patch continue
				if (!_sendEmigrants(curPatch, home, target, nbEmigr,
                                    _migr_rate[FEM], sum_m, factor)) break;
			}
		}
		//curPatch->move(OFFSx,ADLTx);  // already done in _computeTotEmigrants_coal
	}   //end for each home patch
}

// ----------------------------------------------------------------------------------------
// LCE_Disperse
// ----------------------------------------------------------------------------------------
/** island migration model with a female and a male migration rate
 * migration rate to any patch is m/(_nb_patch-1)
 */
void LCE_Disperse::migrate_island_propagule() {
	Patch *curPatch;
	unsigned int home, target, propagule;
	unsigned int nbInd[2];        // nbInd is set by _computeTotEmigrants()
	double factor[2], sum_m[2]; // factor is set by get_migr_factor() & sum_m is set by _computeTotEmigrants()
    
	// for each patch
	vector<Patch*>::iterator curPop, endPop = _popPtr->get_vFullPatch().end();
	for (curPop = _popPtr->get_vFullPatch().begin(); curPop != endPop;
         ++curPop) {
		curPatch = *curPop;
		home = curPatch->get_ID();
		get_migr_factor(curPatch, factor);   // adjust factor
        
		if (_computeTotEmigrants(curPatch, nbInd, _tot_emigRate, sum_m,
                                 factor)) {            // total number of emigrants
			// migration to the propagule patch
			do {                // draw randomly the propagule patch
				propagule = get_pop_ptr()->rand().Uniform(_nb_patch);
			} while (home == propagule);
			if (_sendEmigrants(curPatch, home, propagule, nbInd,
                               _migr_rate_propagule, sum_m, factor)) {
                
				// migration to each other patch
				for (target = 0; target < _nb_patch; ++target) {
					if (home == target) continue; // if it is the same patch continue
					if (target == propagule) continue; // this have already been done
					if (!_sendEmigrants(curPatch, home, target, nbInd,
                                        _migr_rate, sum_m, factor)) break; // send emigrants
				}
			}
		}
		curPatch->move(OFFSx, ADLTx); // all not migrated individuals remain in the original patch
	}          //end for each home patch
}

// ----------------------------------------------------------------------------------------
// LCE_Disperse
// ----------------------------------------------------------------------------------------
/** island migration model with a female and a male migration rate
 * migration rate to any patch is m/(_nb_patch-1)
 */
void LCE_DisperseCoalescence::migrate_island_propagule() {
	Patch *curPatch;
	unsigned int home, target, propagule;
	unsigned int nbEmigr;   // nbEmigr is set by _computeTotEmigrants()
	double factor, sum_m; // factor is set by get_migr_factor() & sum_m is set by _computeTotEmigrants()
    
	// for each patch
	vector<Patch*>::iterator curPop, endPop = _popPtr->get_vFullPatch().end();
	for (curPop = _popPtr->get_vFullPatch().begin(); curPop != endPop;
         ++curPop) {
		curPatch = *curPop;
		home = curPatch->get_ID();
		factor = (this->*get_migr_factor_funcPtr[FEM])(curPatch, FEM); // adjust factor
        
		if (_computeTotEmigrants(curPatch, nbEmigr, _tot_emigRate[FEM], sum_m,
                                 factor)) {            // total number of emigrants
			// migration to the propagule patch
			do {                // draw randomly the propagule patch
				propagule = get_pop_ptr()->rand().Uniform(_nb_patch);
			} while (home == propagule);
			if (_sendEmigrants(curPatch, home, propagule, nbEmigr,
                               _migr_rate_propagule[FEM], sum_m, factor)) {
                
				// migration to each other patch
				for (target = 0; target < _nb_patch; ++target) {
					if (home == target) continue; // if it is the same patch continue
					if (target == propagule) continue; // this have already been done
					if (!_sendEmigrants(curPatch, home, target, nbEmigr,
                                        _migr_rate[FEM], sum_m, factor)) break; // send emigrants
				}
			}
		}
		//curPatch->move(OFFSx,ADLTx);  // already done in _computeTotEmigrants_coal
	}  //end for each home patch
}

// ----------------------------------------------------------------------------------------
// migrate_2_neighbours
// ----------------------------------------------------------------------------------------
/** migration to 2 neighburs (1D-SS model) */
void LCE_Disperse::migrate_2_neighbours(Patch* curPatch,
                                        const unsigned int& home, const unsigned int& n1,
                                        const unsigned int& n2, double* m1, double* m2, double* tot_m) {
	assert(curPatch->get_ID()==home);
	unsigned int nbInd[2];        // nbInd is set by _computeTotEmigrants()
	double factor[2], sum_m[2]; // factor is set by get_migr_factor() & sum_m is set by _computeTotEmigrants()
	get_migr_factor(curPatch, factor);
    
	if (_computeTotEmigrants(curPatch, nbInd, tot_m, sum_m, factor)) // total number of emigrants
		if (_sendEmigrants(curPatch, home, n1, nbInd, m1, sum_m, factor)) // migrate left
			if (_sendEmigrants(curPatch, home, n2, nbInd, m2, sum_m, factor)) {
			};   // migrate right
	curPatch->move(OFFSx, ADLTx); // not migrated individuals remain in the original patch
	assert(!nbInd[FEM] && !nbInd[MAL]);
}

// ----------------------------------------------------------------------------------------
// migrate_2_neighbours
// ----------------------------------------------------------------------------------------
/** migration to 2 neighburs (1D-SS model) */
void LCE_DisperseCoalescence::migrate_2_neighbours(Patch* curPatch,
                                                   const unsigned int& home, const unsigned int& n1,
                                                   const unsigned int& n2, double* m1, double* m2, double* tot_m) {
	assert(curPatch->get_ID()==home);
	unsigned int nbEmigr;   // nbEmigr is set by _computeTotEmigrants()
	double factor, sum_m; // factor is set by get_migr_factor() & sum_m is set by _computeTotEmigrants()
	factor = (this->*get_migr_factor_funcPtr[FEM])(curPatch, FEM); // female factor
    
	if (_computeTotEmigrants(curPatch, nbEmigr, tot_m[FEM], sum_m, factor)) // total number of emigrants
		if (_sendEmigrants(curPatch, home, n1, nbEmigr, m1[FEM], sum_m, factor)) // migrate left
			if (_sendEmigrants(curPatch, home, n2, nbEmigr, m2[FEM], sum_m,
                               factor)) {
			};   // migrate right
	//curPatch->move(OFFSx,ADLTx);  // already done in _computeTotEmigrants_coal
	assert(!nbEmigr);
}

// ----------------------------------------------------------------------------------------
// migrate_4_neighbours
// ----------------------------------------------------------------------------------------
/** migration to 4 neighbors (2D-SS model with 4 neighbors) */
void LCE_Disperse::migrate_4_neighbours(Patch* curPatch,
                                        const unsigned int& home, const unsigned int& n1,
                                        const unsigned int& n2, const unsigned int& n3, const unsigned int& n4,
                                        double* m1, double* m2, double* m3, double* m4, double* tot_m) {
	assert(curPatch->get_ID()==home);
	unsigned int nbInd[2];        // nbInd is set by _computeTotEmigrants()
	double factor[2], sum_m[2]; // factor is set by get_migr_factor() & sum_m is set by _computeTotEmigrants()
	get_migr_factor(curPatch, factor);
    
	if (_computeTotEmigrants(curPatch, nbInd, tot_m, sum_m, factor)) // total number of emigrants
		if (_sendEmigrants(curPatch, home, n1, nbInd, m1, sum_m, factor)) // migrate N
			if (_sendEmigrants(curPatch, home, n2, nbInd, m2, sum_m, factor)) // migrate E
				if (_sendEmigrants(curPatch, home, n3, nbInd, m3, sum_m,
                                   factor))       // migrate S
					if (_sendEmigrants(curPatch, home, n4, nbInd, m4, sum_m,
                                       factor)) {
					}   // migrate W
	curPatch->move(OFFSx, ADLTx); // not migrated individuals remain in the original patch
	assert(!nbInd[FEM] && !nbInd[MAL]);
}

// ----------------------------------------------------------------------------------------
// migrate_4_neighbours
// ----------------------------------------------------------------------------------------
/** migration to 4 neighbors (2D-SS model with 4 neighbors) */
void LCE_DisperseCoalescence::migrate_4_neighbours(Patch* curPatch,
                                                   const unsigned int& home, const unsigned int& n1,
                                                   const unsigned int& n2, const unsigned int& n3, const unsigned int& n4,
                                                   double* m1, double* m2, double* m3, double* m4, double* tot_m) {
	assert(curPatch->get_ID()==home);
	unsigned int nbInd;        // nbInd is set by _computeTotEmigrants()
	double factor, sum_m; // factor is set by get_migr_factor() & sum_m is set by _computeTotEmigrants()
	factor = (this->*get_migr_factor_funcPtr[FEM])(curPatch, FEM); // female factor
    
	if (_computeTotEmigrants(curPatch, nbInd, tot_m[FEM], sum_m, factor)) // total number of emigrants
		if (_sendEmigrants(curPatch, home, n1, nbInd, m1[FEM], sum_m, factor)) // migrate N
			if (_sendEmigrants(curPatch, home, n2, nbInd, m2[FEM], sum_m,
                               factor))         // migrate E
				if (_sendEmigrants(curPatch, home, n3, nbInd, m3[FEM], sum_m,
                                   factor))       // migrate S
					if (_sendEmigrants(curPatch, home, n4, nbInd, m4[FEM],
                                       sum_m, factor)) {
					}   // migrate W
	//curPatch->move(OFFSx,ADLTx);  // already done in _computeTotEmigrants_coal
	assert(!nbInd);
}

// ----------------------------------------------------------------------------------------
// migrate_8_neighbours
// ----------------------------------------------------------------------------------------
/** migration to 8 neighbors (2D-SS model with 8 neighbors) */
void LCE_Disperse::migrate_8_neighbours(Patch* curPatch,
                                        const unsigned int& home, const unsigned int& n1,
                                        const unsigned int& n2, const unsigned int& n3, const unsigned int& n4,
                                        const unsigned int& n5, const unsigned int& n6, const unsigned int& n7,
                                        const unsigned int& n8, double* m1, double* m2, double* m3, double* m4,
                                        double* m5, double* m6, double* m7, double* m8, double* tot_m) {
	assert(curPatch->get_ID()==home);
	unsigned int nbInd[2];        // nbInd is set by _computeTotEmigrants()
	double factor[2], sum_m[2]; // factor is set by get_migr_factor() & sum_m is set by _computeTotEmigrants()
	get_migr_factor(curPatch, factor);
    
	if (_computeTotEmigrants(curPatch, nbInd, tot_m, sum_m, factor)) // total number of emigrants
		if (_sendEmigrants(curPatch, home, n1, nbInd, m1, sum_m, factor)) // migrate NW
			if (_sendEmigrants(curPatch, home, n2, nbInd, m2, sum_m, factor)) // migrate N
				if (_sendEmigrants(curPatch, home, n3, nbInd, m3, sum_m,
                                   factor))             // migrate NE
					if (_sendEmigrants(curPatch, home, n4, nbInd, m4, sum_m,
                                       factor))           // migrate E
						if (_sendEmigrants(curPatch, home, n5, nbInd, m5, sum_m,
                                           factor))         // migrate SE
							if (_sendEmigrants(curPatch, home, n6, nbInd, m6,
                                               sum_m, factor))       // migrate S
								if (_sendEmigrants(curPatch, home, n7, nbInd,
                                                   m7, sum_m, factor))     // migrate SE
									if (_sendEmigrants(curPatch, home, n8,
                                                       nbInd, m8, sum_m, factor)) {
									} // migrate W
	curPatch->move(OFFSx, ADLTx); // not migrated individuals remain in the original patch
	assert(!nbInd[FEM] && !nbInd[MAL]);
}

// ----------------------------------------------------------------------------------------
// migrate_8_neighbours
// ----------------------------------------------------------------------------------------
/** migration to 8 neighburs (2D-SS model with 8 neighbours) */
void LCE_DisperseCoalescence::migrate_8_neighbours(Patch* curPatch,
                                                   const unsigned int& home, const unsigned int& n1,
                                                   const unsigned int& n2, const unsigned int& n3, const unsigned int& n4,
                                                   const unsigned int& n5, const unsigned int& n6, const unsigned int& n7,
                                                   const unsigned int& n8, double* m1, double* m2, double* m3, double* m4,
                                                   double* m5, double* m6, double* m7, double* m8, double* tot_m) {
	assert(curPatch->get_ID()==home);
	unsigned int nbInd;        // nbInd is set by _computeTotEmigrants()
	double factor, sum_m; // factor is set by get_migr_factor() & sum_m is set by _computeTotEmigrants()
	factor = (this->*get_migr_factor_funcPtr[FEM])(curPatch, FEM); // female factor
    
	if (_computeTotEmigrants(curPatch, nbInd, tot_m[FEM], sum_m, factor)) // total number of emigrants
		if (_sendEmigrants(curPatch, home, n1, nbInd, m1[FEM], sum_m, factor)) // migrate NW
			if (_sendEmigrants(curPatch, home, n2, nbInd, m2[FEM], sum_m,
                               factor))               // migrate N
				if (_sendEmigrants(curPatch, home, n3, nbInd, m3[FEM], sum_m,
                                   factor))             // migrate NE
					if (_sendEmigrants(curPatch, home, n4, nbInd, m4[FEM],
                                       sum_m, factor))           // migrate E
						if (_sendEmigrants(curPatch, home, n5, nbInd, m5[FEM],
                                           sum_m, factor))         // migrate SE
							if (_sendEmigrants(curPatch, home, n6, nbInd,
                                               m6[FEM], sum_m, factor))       // migrate S
								if (_sendEmigrants(curPatch, home, n7, nbInd,
                                                   m7[FEM], sum_m, factor))  // migrate SE
									if (_sendEmigrants(curPatch, home, n8,
                                                       nbInd, m8[FEM], sum_m, factor)) {
									} // migrate W
	//curPatch->move(OFFSx,ADLTx);  // already done in _computeTotEmigrants_coal
	assert(!nbInd);
}

// ----------------------------------------------------------------------------------------
// LCE_Disperse
// ----------------------------------------------------------------------------------------
/** select model of disersal depending on the number of populated patches */
void LCE_Disperse::migrate_1D_ss() {
	if (_popPtr->get_fullPatch_ratio() < 0.8)
		migrate_1D_ss_full();
	else migrate_1D_ss_all();
}

// ----------------------------------------------------------------------------------------
// LCE_Disperse
// ----------------------------------------------------------------------------------------
/** all patches are iterated */
void LCE_Disperse::migrate_1D_ss_all() {
	unsigned int id = 0;
	vector<Patch*>::iterator curPop = _popPtr->get_vPatch().begin();
    
	// first patch (edge patch)
	migrate_2_neighbours(*curPop, id, _nb_patch - 1, 1, _migr_rateOut,
                         _migr_rateIn, _tot_emigRate);
    
	// for each patch in the middle: 1-_nb_patch-2
	for (++curPop, ++id; id < _nb_patch - 1; ++curPop, ++id) {
		migrate_2_neighbours(*curPop, id, id - 1, id + 1, _migr_rate,
                             _migr_rate, _tot_emigRate);
	}
    
	// edge patches (last patch)
	migrate_2_neighbours(*curPop, id, id - 1, 0, _migr_rateIn, _migr_rateOut,
                         _tot_emigRate);
}

// ----------------------------------------------------------------------------------------
// LCE_Disperse
// ----------------------------------------------------------------------------------------
/** only the populated patches are iterated */
void LCE_Disperse::migrate_1D_ss_full() {
	unsigned int id;
	vector<Patch*>::iterator curPop = _popPtr->get_vFullPatch().begin();
	vector<Patch*>::iterator endPop = _popPtr->get_vFullPatch().end();
	for (; curPop != endPop; ++curPop) {
		id = (*curPop)->get_ID();
		if (id == 0)
			migrate_2_neighbours(*curPop, id, _nb_patch - 1, 1, _migr_rateOut,
                                 _migr_rateIn, _tot_emigRate);
		else if (id == _nb_patch - 1)
			migrate_2_neighbours(*curPop, id, id - 1, 0, _migr_rateIn,
                                 _migr_rateOut, _tot_emigRate);
		else migrate_2_neighbours(*curPop, id, id - 1, id + 1, _migr_rate,
                                  _migr_rate, _tot_emigRate);
	}
}

// ----------------------------------------------------------------------------------------
// LCE_Disperse
// ----------------------------------------------------------------------------------------
/** select model of dispersal depending on the number of populated patches */
void LCE_Disperse::migrate_2D_ss_4Neighbour() {
	if (_popPtr->get_fullPatch_ratio() < 0.8)
		migrate_2D_ss_4Neighbour_full();
	else migrate_2D_ss_4Neighbour_all();
}

// ----------------------------------------------------------------------------------------
// LCE_Disperse
// ----------------------------------------------------------------------------------------
/** all patches are iterated */
void LCE_Disperse::migrate_2D_ss_4Neighbour_all() {
	unsigned int id = 0, x, y;
	// unsigned int xx, yy;
	vector<Patch*>::iterator curPop = _popPtr->get_vPatch().begin();
    
	// upper edge ////////////////////////////////////////////////////////////////
	// left corner       order of neighbours: N,            E,                S,                W
	//ID2Coord(id, xx, yy);
	//assert(xx==0 && yy==0);
	migrate_4_neighbours(*curPop, id, _nb_patch - _x_size, 1, _x_size,
                         _x_size - 1, _migr_rateOut, _migr_rateCorner, _migr_rateCorner,
                         _migr_rateOut, _tot_emigRate);
    
	// upper edge
	for (x = 1, ++curPop, ++id; x <= _x_size - 2; ++x, ++curPop, ++id) {
		//ID2Coord(id, xx, yy);
		//assert(yy==0);
		migrate_4_neighbours(*curPop, id, _nb_patch - _x_size + id, id + 1,
                             id + _x_size, id - 1, _migr_rateOut, _migr_rateIn, _migr_rateIn,
                             _migr_rateIn, _tot_emigRate);
	}
    
	// right corner
	//ID2Coord(id, xx, yy);
	//assert(xx==_x_size-1 && yy==0);
	migrate_4_neighbours(*curPop, id, _nb_patch - 1, 0, id + _x_size, id - 1,
                         _migr_rateOut, _migr_rateOut, _migr_rateCorner, _migr_rateCorner,
                         _tot_emigRate);
    
	// centre ////////////////////////////////////////////////////////////////////
	for (y = 1, ++curPop, ++id; y <= _y_size - 2; ++y, ++curPop, ++id) {
		// left edge
		//ID2Coord(id, xx, yy);
		//assert(xx==0);
		migrate_4_neighbours(*curPop, id, id - _x_size, id + 1, id + _x_size,
                             id + _x_size - 1, _migr_rateIn, _migr_rateIn, _migr_rateIn,
                             _migr_rateOut, _tot_emigRate);
        
		// centre
		for (x = 1, ++curPop, ++id; x <= _x_size - 2; ++x, ++curPop, ++id) {
			//ID2Coord(id, xx, yy);
			//assert(xx>0 && xx<_x_size-1 && yy>0 && yy<_y_size-1);
			migrate_4_neighbours(*curPop, id, id - _x_size, id + 1,
                                 id + _x_size, id - 1, _migr_rate, _migr_rate, _migr_rate,
                                 _migr_rate, _tot_emigRate);
		}
        
		// right edge
		//ID2Coord(id, xx, yy);
		//assert(xx==_x_size-1);
		migrate_4_neighbours(*curPop, id, id - _x_size, id - _x_size + 1,
                             id + _x_size, id - 1, _migr_rateIn, _migr_rateOut, _migr_rateIn,
                             _migr_rateIn, _tot_emigRate);
	}
    
	// lower edge ////////////////////////////////////////////////////////////////
	// left corner       order of neighbours: N,            E,                S,                W
	//ID2Coord(id, xx, yy);
	//assert(xx==0);
	migrate_4_neighbours(*curPop, id, id - _x_size, id + 1, 0, _nb_patch - 1,
                         _migr_rateCorner, _migr_rateCorner, _migr_rateOut, _migr_rateOut,
                         _tot_emigRate);
    
	// lower edge
	for (x = 1, ++curPop, ++id; x <= _x_size - 2; ++x, ++curPop, ++id) {
		//ID2Coord(id, xx, yy);
		//assert(yy==_y_size-1);
		migrate_4_neighbours(*curPop, id, id - _x_size, id + 1,
                             _x_size - (_nb_patch - id), id - 1, _migr_rateIn, _migr_rateIn,
                             _migr_rateOut, _migr_rateIn, _tot_emigRate);
	}
    
	// right corner
	//ID2Coord(id, xx, yy);
	//assert(xx==_x_size-1 && yy==_y_size-1);
	migrate_4_neighbours(*curPop, id, id - _x_size, _nb_patch - _x_size,
                         _x_size - 1, id - 1, _migr_rateCorner, _migr_rateOut, _migr_rateOut,
                         _migr_rateCorner, _tot_emigRate);
}

// ----------------------------------------------------------------------------------------
// LCE_Disperse
// ----------------------------------------------------------------------------------------
/** only the populated patches are iterated */
void LCE_Disperse::migrate_2D_ss_4Neighbour_full() {
	unsigned int id, cur_x, cur_y, x = _x_size, y = _y_size;
	vector<Patch*>::iterator curPop, endPop = _popPtr->get_vFullPatch().end();
	for (curPop = _popPtr->get_vFullPatch().begin(); curPop != endPop;
         ++curPop) {   // for each populated patch
		id = (*curPop)->get_ID();
        
		ID2Coord(id, cur_x, cur_y);
		assert((*curPop)->size(OFFSx));
        
		// upper edge //////////////////////////////////////////////////////////////
		if (cur_y == 0) {
			// left corner
			if (cur_x == 0)
				migrate_4_neighbours(*curPop, id, _nb_patch - x, 1, x, x - 1,
                                     _migr_rateOut, _migr_rateCorner, _migr_rateCorner,
                                     _migr_rateOut, _tot_emigRate);
			// uper edge
			else if (cur_x < x - 1)
				migrate_4_neighbours(*curPop, id, _nb_patch - x + id, id + 1,
                                     id + x, id - 1, _migr_rateOut, _migr_rateIn,
                                     _migr_rateIn, _migr_rateIn, _tot_emigRate);
			// right corner
			else migrate_4_neighbours(*curPop, id, _nb_patch - 1, 0, id + x,
                                      id - 1, _migr_rateOut, _migr_rateOut, _migr_rateCorner,
                                      _migr_rateCorner, _tot_emigRate);
		}
        
		// centre //////////////////////////////////////////////////////////////////
		else if (cur_y < y - 1) {
			// left edge
			if (cur_x == 0)
				migrate_4_neighbours(*curPop, id, id - x, id + 1, id + x,
                                     id + x - 1, _migr_rateIn, _migr_rateIn, _migr_rateIn,
                                     _migr_rateOut, _tot_emigRate);
			// centre
			else if (cur_x < x - 1)
				migrate_4_neighbours(*curPop, id, id - x, id + 1, id + x,
                                     id - 1, _migr_rate, _migr_rate, _migr_rate, _migr_rate,
                                     _tot_emigRate);
			// right edge
			else migrate_4_neighbours(*curPop, id, id - x, id - x + 1, id + x,
                                      id - 1, _migr_rateIn, _migr_rateOut, _migr_rateIn,
                                      _migr_rateIn, _tot_emigRate);
		}
        
		// lower edge //////////////////////////////////////////////////////////////
		else {
			// left corner
			if (cur_x == 0)
				migrate_4_neighbours(*curPop, id, id - x, id + 1, 0,
                                     _nb_patch - 1, _migr_rateCorner, _migr_rateCorner,
                                     _migr_rateOut, _migr_rateOut, _tot_emigRate);
			// lower edge
			else if (cur_x < x - 1)
				migrate_4_neighbours(*curPop, id, id - x, id + 1,
                                     id - _nb_patch + x, id - 1, _migr_rateIn, _migr_rateIn,
                                     _migr_rateOut, _migr_rateIn, _tot_emigRate);
			// right corner
			else migrate_4_neighbours(*curPop, id, id - x, _nb_patch - x, x - 1,
                                      id - 1, _migr_rateCorner, _migr_rateOut, _migr_rateOut,
                                      _migr_rateCorner, _tot_emigRate);
		}
	}
}

// ----------------------------------------------------------------------------------------
// LCE_Disperse
// ----------------------------------------------------------------------------------------
/** select model of dispersal depending on the number of populated patches */
void LCE_Disperse::migrate_2D_ss_8Neighbour() {
	if (_popPtr->get_fullPatch_ratio() < 0.8)
		migrate_2D_ss_8Neighbour_full();
	else migrate_2D_ss_8Neighbour_all();
}

// ----------------------------------------------------------------------------------------
// LCE_Disperse
// ----------------------------------------------------------------------------------------
/** all patches are iterated */
void LCE_Disperse::migrate_2D_ss_8Neighbour_all() {
	unsigned int id = 0, x, y;
	// unsigned int xx, yy;
	vector<Patch*>::iterator curPop = _popPtr->get_vPatch().begin();
    
	// upper edge ////////////////////////////////////////////////////////////////
	// left corner       order of neighbours: NW,  N,   NE,    E,   SE,   S,   SW,   W
	//ID2Coord(id, xx, yy);
	//assert(yy==0 && xx==0);
	migrate_8_neighbours(*curPop, id, _nb_patch - 1, _nb_patch - _x_size,
                         _nb_patch - _x_size + 1, 1, _x_size + 1, _x_size, 2 * _x_size - 1,
                         _x_size - 1, _migr_rateOut, _migr_rateOut, _migr_rateOut,
                         _migr_rateCorner, _migr_rateCorner, _migr_rateCorner, _migr_rateOut,
                         _migr_rateOut, _tot_emigRate);
    
	// upper edge
	for (x = 1, ++curPop, ++id; x <= _x_size - 2; ++x, ++curPop, ++id) {
		//ID2Coord(id, xx, yy);
		//assert(yy==0);
		migrate_8_neighbours(*curPop, id, _nb_patch - _x_size + id - 1,
                             _nb_patch - _x_size + id, _nb_patch - _x_size + id + 1, id + 1,
                             id + _x_size + 1, id + _x_size, id + _x_size - 1, id - 1,
                             _migr_rateOut, _migr_rateOut, _migr_rateOut, _migr_rateIn,
                             _migr_rateIn, _migr_rateIn, _migr_rateIn, _migr_rateIn,
                             _tot_emigRate);
	}
    
	// right corner
	//ID2Coord(id, xx, yy);
	//assert(yy==0 && xx==_x_size-1);
	migrate_8_neighbours(*curPop, id, _nb_patch - 2, _nb_patch - 1,
                         _nb_patch - _x_size, 0, _x_size, id + _x_size, id + _x_size - 1,
                         id - 1, _migr_rateOut, _migr_rateOut, _migr_rateOut, _migr_rateOut,
                         _migr_rateOut, _migr_rateCorner, _migr_rateCorner, _migr_rateCorner,
                         _tot_emigRate);
    
	// centre ////////////////////////////////////////////////////////////////////
	for (y = 1, ++curPop, ++id; y <= _y_size - 2; ++y, ++curPop, ++id) {
		// left edge
		//ID2Coord(id, xx, yy);
		//assert(yy>0 && yy<_y_size-1 && xx==0);
		migrate_8_neighbours(*curPop, id, id - 1, id - _x_size,
                             id - _x_size + 1, id + 1, id + _x_size + 1, id + _x_size,
                             id + 2 * _x_size - 1, id + _x_size - 1, _migr_rateOut,
                             _migr_rateIn, _migr_rateIn, _migr_rateIn, _migr_rateIn,
                             _migr_rateIn, _migr_rateOut, _migr_rateOut, _tot_emigRate);
        
		// centre
		for (x = 1, ++curPop, ++id; x <= _x_size - 2; ++x, ++curPop, ++id) {
			//ID2Coord(id, xx, yy);
			//assert(yy>0 && yy<_y_size && xx<_x_size && xx>0);
			migrate_8_neighbours(*curPop, id, id - 1, id - _x_size,
                                 id - _x_size + 1, id + 1, id + _x_size + 1, id + _x_size,
                                 id + _x_size - 1, id - x - 1, _migr_rate, _migr_rate,
                                 _migr_rate, _migr_rate, _migr_rate, _migr_rate, _migr_rate,
                                 _migr_rate, _tot_emigRate);
		}
        
		// right edge
		//ID2Coord(id, xx, yy);
		//assert(yy>0 && yy<_y_size && xx==_x_size-1);
		migrate_8_neighbours(*curPop, id, id - _x_size - 1, id - _x_size,
                             id - 2 * _x_size + 1, id - _x_size + 1, id + 1, id + _x_size,
                             id + _x_size - 1, id - 1, _migr_rateIn, _migr_rateIn,
                             _migr_rateOut, _migr_rateOut, _migr_rateOut, _migr_rateIn,
                             _migr_rateIn, _migr_rateIn, _tot_emigRate);
	}
    
	// lower edge ////////////////////////////////////////////////////////////////
	// left corner
	//ID2Coord(id, xx, yy);
	//assert(yy==_y_size-1 && xx==0);
	migrate_8_neighbours(*curPop, id, id - 1, id - _x_size, id - _x_size + 1,
                         id + 1, 1, 0, _x_size - 1, _nb_patch - 1, _migr_rateOut,
                         _migr_rateCorner, _migr_rateCorner, _migr_rateCorner, _migr_rateOut,
                         _migr_rateOut, _migr_rateOut, _migr_rateOut, _tot_emigRate);
    
	// lower edge
	for (x = 1, ++curPop, ++id; x <= _x_size - 2; ++x, ++curPop, ++id) {
		//ID2Coord(id, xx, yy);
		//assert(yy==_y_size-1 && xx<_x_size && xx>0);
		migrate_8_neighbours(*curPop, id, id - _x_size - 1, id - _x_size,
                             id - _x_size + 1, id + 1, _x_size - _nb_patch + id + 1,
                             _x_size - _nb_patch + id, _x_size - _nb_patch + id - 1, id - 1,
                             _migr_rateIn, _migr_rateIn, _migr_rateIn, _migr_rateIn,
                             _migr_rateOut, _migr_rateOut, _migr_rateOut, _migr_rateIn,
                             _tot_emigRate);
	}
    
	// right corner
	//ID2Coord(id, xx, yy);
	//assert(yy==_y_size-1 && xx==_x_size-1);
	migrate_8_neighbours(*curPop, id, id - _x_size - 1, id - _x_size,
                         id - 2 * _x_size + 1, _nb_patch - _x_size, 0, _x_size - 1,
                         _x_size - 2, id - 1, _migr_rateCorner, _migr_rateCorner,
                         _migr_rateOut, _migr_rateOut, _migr_rateOut, _migr_rateOut,
                         _migr_rateOut, _migr_rateCorner, _tot_emigRate);
}

// ----------------------------------------------------------------------------------------
// LCE_Disperse
// ----------------------------------------------------------------------------------------
/** only the populated patches are iterated */
void LCE_Disperse::migrate_2D_ss_8Neighbour_full() {
	unsigned int id, cur_x, cur_y, x = _x_size, y = _y_size;
	vector<Patch*>::iterator curPop = _popPtr->get_vFullPatch().begin();
	vector<Patch*>::iterator endPop = _popPtr->get_vFullPatch().end();
	for (; curPop != endPop; ++curPop) {   // for each popualted patch
		id = (*curPop)->get_ID();
		ID2Coord(id, cur_x, cur_y);
		assert((*curPop)->size(ALL));
        
		// upper edge //////////////////////////////////////////////////////////////
		if (cur_y == 0) {
			// left corner
			if (cur_x == 0)
				migrate_8_neighbours(*curPop, id, _nb_patch - 1, _nb_patch - x,
                                     _nb_patch - x + 1, 1, x + 1, x, 2 * x - 1, x - 1,
                                     _migr_rateOut, _migr_rateOut, _migr_rateOut,
                                     _migr_rateCorner, _migr_rateCorner, _migr_rateCorner,
                                     _migr_rateOut, _migr_rateOut, _tot_emigRate);
			// uper edge
			else if (cur_x < x - 1)
				migrate_8_neighbours(*curPop, id, _nb_patch - x + id - 1,
                                     _nb_patch - x + id, _nb_patch - x + id + 1, id + 1,
                                     id + x + 1, id + x, id + x - 1, id - 1, _migr_rateOut,
                                     _migr_rateOut, _migr_rateOut, _migr_rateIn,
                                     _migr_rateIn, _migr_rateIn, _migr_rateIn, _migr_rateIn,
                                     _tot_emigRate);
			// right corner
			else migrate_8_neighbours(*curPop, id, _nb_patch - 2, _nb_patch - 1,
                                      _nb_patch - x, 0, x, id + x, id + x - 1, id - 1,
                                      _migr_rateOut, _migr_rateOut, _migr_rateOut, _migr_rateOut,
                                      _migr_rateOut, _migr_rateCorner, _migr_rateCorner,
                                      _migr_rateCorner, _tot_emigRate);
		}
        
		// centre //////////////////////////////////////////////////////////////////
		else if (cur_y < y - 1) {
			// left edge
			if (cur_x == 0)
				migrate_8_neighbours(*curPop, id, id - 1, id - x,
                                     id - x + 1, id + 1, id + x + 1, id + x, id + 2 * x - 1,
                                     id + x - 1, _migr_rateOut, _migr_rateIn, _migr_rateIn,
                                     _migr_rateIn, _migr_rateIn, _migr_rateIn, _migr_rateOut,
                                     _migr_rateOut, _tot_emigRate);
			// centre
			else if (cur_x < x - 1)
				migrate_8_neighbours(*curPop, id, id - 1, id - x, id - x + 1,
                                     id + 1, id + x + 1, id + x, id + x - 1, id - 1,
                                     _migr_rate, _migr_rate, _migr_rate, _migr_rate,
                                     _migr_rate, _migr_rate, _migr_rate, _migr_rate,
                                     _tot_emigRate);
			// right edge
			else migrate_8_neighbours(*curPop, id, id - x - 1, id - x,
                                      id - 2 * x + 1, id - x + 1, id + 1, id + x, id + x - 1,
                                      id - 1, _migr_rateIn, _migr_rateIn, _migr_rateOut,
                                      _migr_rateOut, _migr_rateOut, _migr_rateIn, _migr_rateIn,
                                      _migr_rateIn, _tot_emigRate);
		}
        
		// lower edge //////////////////////////////////////////////////////////////
		else {
			// left corner
			if (cur_x == 0)
				migrate_8_neighbours(*curPop, id, id - 1, id - x, id - x + 1,
                                     id + 1, 1, 0, x - 1, _nb_patch - 1, _migr_rateOut,
                                     _migr_rateCorner, _migr_rateCorner, _migr_rateCorner,
                                     _migr_rateOut, _migr_rateOut, _migr_rateOut,
                                     _migr_rateOut, _tot_emigRate);
			// lower edge
			else if (cur_x < x - 1)
				migrate_8_neighbours(*curPop, id, id - x - 1, id - x,
                                     id - x + 1, id + 1, id - _nb_patch + x + 1,
                                     id - _nb_patch + x, id - _nb_patch + x - 1, id - 1,
                                     _migr_rateIn, _migr_rateIn, _migr_rateIn, _migr_rateIn,
                                     _migr_rateOut, _migr_rateOut, _migr_rateOut,
                                     _migr_rateIn, _tot_emigRate);
			// right corner
			else migrate_8_neighbours(*curPop, id, id - x - 1, id - x,
                                      id - 2 * x + 1, _nb_patch - x, 0, x - 1, x - 2, id - 1,
                                      _migr_rateCorner, _migr_rateCorner, _migr_rateOut,
                                      _migr_rateOut, _migr_rateOut, _migr_rateOut, _migr_rateOut,
                                      _migr_rateCorner, _tot_emigRate);
		}
	}
}

// ----------------------------------------------------------------------------------------
// _long_range_geometric
// ----------------------------------------------------------------------------------------
/** this function allows to simulate long range dispersal following a geometric function.
 * the edges are threated as abosrbing, i.e. all emigrants crossing an edge will disapear
 */
void LCE_Disperse::migrate_2D_ss_geometric() {
	migrate_2D_ss_geometric(FEM);
	migrate_2D_ss_geometric(MAL);
}

// ----------------------------------------------------------------------------------------
// _long_range_geometric
// ----------------------------------------------------------------------------------------
/** this function allows to simulate long range dispersal following a geometric function.
 * the edges are threated as abosrbing, i.e. all emigrants crossing an edge will disapear
 */
void LCE_Disperse::migrate_2D_ss_geometric(const sex_t& SEX) {
	if (!_migr_rate[SEX]) return;
	assert(
           _disp_long_range_coef[SEX]>=1e-10 && _disp_long_range_coef[SEX]<=1.0);
    
	double dist, angle;
	unsigned int popSize, nbEmigrants, x, y, id;
	int newX, newY;          // can be negative!!!
	vector<Patch*>::iterator curPop, endPop = _popPtr->get_vFullPatch().end();
	for (curPop = _popPtr->get_vFullPatch().begin(); curPop != endPop;
         ++curPop) {  // for each populated patch
		popSize = (*curPop)->size(SEX, OFFSx);
		if (!popSize) continue;                                   // empty patch
		nbEmigrants = get_pop_ptr()->rand().Binomial(_migr_rate[SEX], popSize);
		if (!nbEmigrants) continue;                              // no emigrants
		id = (*curPop)->get_ID();
		ID2Coord(id, x, y);
        
		for (; nbEmigrants; --nbEmigrants, --popSize) {     // for each emigrant
			// draw randomly a target patch
			dist = 1 + get_pop_ptr()->rand().Geometric(_disp_long_range_coef[SEX]);
			angle = 2 * PI * get_pop_ptr()->rand().Uniform();
            
			// get the coordinates of the target patch
			if (angle <= (PI / 2.0)) {                            // 1. quadrant
				newX = x + my_round(sin(angle) * dist);
				newY = y + my_round(cos(angle) * dist);
			}
			else if (angle <= PI) {                             // 2. quadrant
				angle -= (PI / 2.0);
				newX = x + my_round(cos(angle) * dist);
				newY = y - my_round(sin(angle) * dist);
			}
			else if (angle <= (3.0 * PI / 2.0)) {                 // 3. quadrant
				angle -= (PI);
				newX = x - my_round(sin(angle) * dist);
				newY = y - my_round(cos(angle) * dist);
			}
			else {                                       		// 4. quadrant
				angle -= (3.0 * PI / 2.0);
				newX = x - my_round(cos(angle) * dist);
				newY = y + my_round(sin(angle) * dist);
			}
            
			// migrate to this patch
			if ((newX < 0 || (unsigned int) newX >= _x_size)
                || (newY < 0 || (unsigned int) newY >= _y_size))
				_absorbEmigrant(*curPop, 1, SEX, popSize); // outside of the range
			else _sendEmigrant(*curPop, id, Coord2ID(newX, newY), 1, SEX,
                               popSize);                // a successful candidate
		}
	}
}

// ----------------------------------------------------------------------------------------
// LCE_Disperse::init
// ----------------------------------------------------------------------------------------
bool LCE_Disperse::init(TMetapop* popPtr) {
	LCE::init(popPtr);
    
	// if only one patch is used there is simply no migration
	_nb_patch = popPtr->getPatchNbr();
	if (_nb_patch == 1) {
		migration_func_ptr = &LCE_Disperse::migrate_zeroMigration;
		return true;
	}
    
	// get parameters
	_disp_model = (int) _paramSet->getValue("dispersal_model");
	_border_model = (int) _paramSet->getValue("dispersal_border_model");
	_lattice_range = (int) _paramSet->getValue("dispersal_lattice_range");
    
	//reset matrixes
	if (_dispMatrix[FEM]) {
		delete _dispMatrix[FEM];
		_dispMatrix[FEM] = NULL;
	}
	if (_dispMatrix[MAL]) {
		delete _dispMatrix[MAL];
		_dispMatrix[MAL] = NULL;
	}
    
	// dispersal is set by matrixes
	if ((_paramSet->isSet("dispersal_rate")
         || _paramSet->isSet("dispersal_rate_fem")
         || _paramSet->isSet("dispersal_rate_mal"))
        && (_paramSet->isMatrix("dispersal_rate")
            || _paramSet->isMatrix("dispersal_rate_fem")
            || _paramSet->isMatrix("dispersal_rate_mal"))) {
            setDispersalMatrix();
        }
    
	// dispersal is set by a single dispersal rate (or default)
	else {
		if (_disp_model != 0 && _nb_patch == 2) { // if there are only 2 patches use the island model
			_disp_model = 0;   // it results always in an island migration model
			warning(
					"With 2 populations all dispersal models are equivalent to the island model: The island model is used!\n");
		}
		if (_disp_model >= 3) _get_lattice_dims(); // if it is a 2D stepping stone model
		setDispersalRate();
	}
    
	_setDispersalFactor();
    _setDispersal_direction();

    
	if (!_popPtr->get_sexInitRatio()) _migr_rate[MAL] = 0; // if hermaphrodites are simulated
    
	return true;
}

// ----------------------------------------------------------------------------------------
// LCE_Disperse::setDispersalFactor
// ----------------------------------------------------------------------------------------
/** set the genralized logistic function parameters
 * check if the function is really needed, i.e. if there is really a change
 */
void LCE_Disperse::_setDispersalFactor() {
	_setDispersalFactor(FEM);
	_setDispersalFactor(MAL);
}

// ----------------------------------------------------------------------------------------
// LCE_Disperse::setDispersalFactor
// ----------------------------------------------------------------------------------------
void LCE_Disperse::_setDispersalFactor(const sex_t& SEX) {
    
	// get the values
	if (_disp_factor[SEX]) delete[] _disp_factor[SEX];
	_disp_factor[SEX] = new double[5];
    
	string curSex = SEX == FEM ? "_fem" : "_mal";
    
	if (_paramSet->get_param("dispersal_rate_model" + curSex)->isSet())
		_disp_factor[SEX][0] = _paramSet->getValue("dispersal_rate_model" + curSex);
	else _disp_factor[SEX][0] = _paramSet->getValue("dispersal_rate_model");

	if(_disp_factor[SEX][0] == 1){
		get_migr_factor_funcPtr[SEX] =
        &LCE_Disperse::get_migr_factor_density_dependant;
	}
	else if(_disp_factor[SEX][0] == 0){
        get_migr_factor_funcPtr[SEX] =
        &LCE_Disperse::get_migr_factor_one;
	}
	else{
        get_migr_factor_funcPtr[SEX] =
        &LCE_Disperse::get_migr_factor_saturation;
	}

}


// ----------------------------------------------------------------------------------------
// LCE_Disperse::_setDispersal_direction
// ----------------------------------------------------------------------------------------
/* emigration or immigration? */
void LCE_Disperse::_setDispersal_direction()
{
    if(_paramSet->getValue("dispersal_direction")==0) return;
    if(migration_func_ptr == &LCE_Disperse::migrate_zeroMigration) return;
    
    if(migration_func_ptr == &LCE_Disperse::migrate_island){
        migration_func_ptr = &LCE_Disperse::immigrate_island;
        warning("Immigration as dispersal is an experimental feature. Use by your own risk!\n");
    }
    else error("Immigration is currently only implemented for the island model!\n");
}

// ----------------------------------------------------------------------------------------
// LCE_Disperse::setDispersalMatrix
// ----------------------------------------------------------------------------------------
/** get the dimesnions of the 2D stepping stone matrix */
void LCE_Disperse::_get_lattice_dims() {
	if (_paramSet->isSet("dispersal_lattice_dims")) {
		TMatrix* m = _paramSet->getMatrix("dispersal_lattice_dims");
        
		// check if the dimensions of the matrix are correct
		if (m->get_dims(NULL) != 2)
			error(
                  "The parameter disersal_lattice_dims should have a matrix with two values!\n");
        
		// get the dimension of the lattice
		_x_size = (unsigned int) m->get(0, 0);
		_y_size = (unsigned int) m->get(0, 1);
		if (_x_size * _y_size != _nb_patch) {
			error(
                  "Parameter disersal_lattice_dims: The dimension of the lattice (%ix%i) does not mach the number of patches (%i)!\n",
                  _x_size, _y_size, _nb_patch);
		}
		delete m;
	}
	else { // if the dimensions are not set, we assume that x=y
		_x_size = _y_size = (unsigned int) sqrt((double) _nb_patch);
        
		if (_x_size * _y_size != _nb_patch) {
			error(
                  "Parameter disersal_lattice_dims: The dimension of the lattice (%ix%i) does not mach the number of patches (%i)!\n",
                  _x_size, _y_size, _nb_patch);
		}
	}
    
	// check if it is really a 2D and not a 1D lattice
	if (_x_size == 1 || _y_size == 1) {
		_disp_model = 2;            // use a 1D stepping stone model
	}
}

// ----------------------------------------------------------------------------------------
// LCE_Disperse::setDispersalMatrix
// ----------------------------------------------------------------------------------------
/** dispersal is set by a dispersal matrix
 * sex specific matrixes have precedence
 * if only one sex specific matrix is set, an error is drawn
 */
void LCE_Disperse::setDispersalMatrix() {
	bool used=false;
	if (_paramSet->isSet("dispersal_rate_fem")) {
		if (_paramSet->isSet("dispersal_rate_mal")) {
			used = setDispMatrix(FEM,
                                 _paramSet->getMatrix("dispersal_rate_fem"))
            | setDispMatrix(MAL,
							_paramSet->getMatrix("dispersal_rate_mal"));
		}
		else error(
                   "Only one sex specific migration matrix is specified: both are required!\n");
	}
    
	// general dispersal matrix
	else used = setDispMatrix(_paramSet->getMatrix("dispersal_rate"));
    
	// function pointer to the migration function
	if (used)
		migration_func_ptr = &LCE_Disperse::migrate_matrix;
	else migration_func_ptr = &LCE_Disperse::migrate_zeroMigration;
}

// ----------------------------------------------------------------------------------------
// LCE_Disperse::setDispersalRate
// ----------------------------------------------------------------------------------------
/** dispersal is set by a dispersal rate
 * sex specific rates have precedence
 * if only one sex specific rate is set, an error is drawn
 */
void LCE_Disperse::setDispersalRate() {
	if (_paramSet->isSet("dispersal_rate_fem")) {
		if (_paramSet->isSet("dispersal_rate_mal")) {
			_migr_rate[FEM] = _paramSet->getValue("dispersal_rate_fem");
			_migr_rate[MAL] = _paramSet->getValue("dispersal_rate_mal");
		}
		else error("Only one sex specific migration rate is specified: both are required!\n");
	}
    
	// general dispersal rate
	else _migr_rate[FEM] = _migr_rate[MAL] = _paramSet->getValue("dispersal_rate");
    
    
    // if there is no migration
    int _disp_rate_model = _paramSet->getValue("dispersal_rate_model");
    if (!_migr_rate[MAL] && !_migr_rate[FEM] && _disp_rate_model < 2) {
        migration_func_ptr = &LCE_Disperse::migrate_zeroMigration;
        return;
    }
    
	// check if a dispersal rate or an absolute number of dispersers is defined
	if (_migr_rate[FEM] > 1)      _rel_abs_disp_rate[FEM] = 1;   // absolute
	else if (_migr_rate[FEM] > -1) _rel_abs_disp_rate[FEM] = 0;   // relative
	else _rel_abs_disp_rate[FEM] = 1;   // absolute
    
	if (_migr_rate[MAL] > 1)      _rel_abs_disp_rate[MAL] = 1;   // absolute
	else if (_migr_rate[MAL] > -1) _rel_abs_disp_rate[MAL] = 0;   // relative
	else _rel_abs_disp_rate[MAL] = 1;   // absolute
    
	// proportion of individuals remaining in the patch
	_tot_emigRate[FEM] = _migr_rate[FEM];
	_tot_emigRate[MAL] = _migr_rate[MAL];
    
	switch (_disp_model) {
		case 0: // island migration model
			_migr_rate[FEM] /= _nb_patch - 1;
			_migr_rate[MAL] /= _nb_patch - 1;
			migration_func_ptr = &LCE_Disperse::migrate_island;
			break;
            
		case 1: // island migration model with propagule pool
			_disp_propagule_prob = _paramSet->getValue(
                                                       "dispersal_propagule_prob");
			_migr_rate_propagule[FEM] = _migr_rate[FEM] * _disp_propagule_prob;
			_migr_rate_propagule[MAL] = _migr_rate[MAL] * _disp_propagule_prob;
			_migr_rate[FEM] *= (1.0 - _disp_propagule_prob) / (_nb_patch - 2);
			_migr_rate[MAL] *= (1.0 - _disp_propagule_prob) / (_nb_patch - 2);
			migration_func_ptr = &LCE_Disperse::migrate_island_propagule;
			break;
            
		case 2: // 1D steppings stone model
			int factorIn, factorOut;
			migration_func_ptr = &LCE_Disperse::migrate_1D_ss;
            
			// edge effect
			switch (_border_model) {
				default:
				case 0:
					factorIn = 2;
					factorOut = 2;
					break;        // circle
				case 1:
					factorIn = 1;
					factorOut = 0;
					break;				// reflecting
				case 2:
					factorIn = 2;
					factorOut = -2;
					break;// absorbing (negative number means removing individuals)
			}
            
			_migr_rateIn[FEM] = _migr_rate[FEM] / factorIn;
			_migr_rateIn[MAL] = _migr_rate[MAL] / factorIn;
			_migr_rateOut[FEM] = factorOut ? _migr_rate[FEM] / factorOut : 0;
			_migr_rateOut[MAL] = factorOut ? _migr_rate[MAL] / factorOut : 0;
			_migr_rate[FEM] /= 2;
			_migr_rate[MAL] /= 2;
			break;
            
		case 3: // 2D stepping stone
			int factorIn4, factorOut4, factorInCorner4;
			int factorIn8, factorOut8, factorInCorner8;
            
			// edge effect
			switch (_border_model) {
				default:
				case 0: // torus
					factorIn4 = 4;
					factorIn8 = 8;
					factorOut4 = 4;
					factorOut8 = 8;
					factorInCorner4 = 4;
					factorInCorner8 = 8;
					break;
				case 1: // reflecting
					factorIn4 = 3;
					factorIn8 = 5;
					factorOut4 = 0;
					factorOut8 = 0;
					factorInCorner4 = 2;
					factorInCorner8 = 3;
					break;
				case 2: // absorbing
					factorIn4 = 4;
					factorIn8 = 8;
					factorOut4 = -4;
					factorOut8 = -8; // negative number means removing individuals
					factorInCorner4 = 4;
					factorInCorner8 = 8;
					break;
			}
            
			// number of neighbors
			if (_lattice_range == 0) {         // 4 Neighbors
				migration_func_ptr = &LCE_Disperse::migrate_2D_ss_4Neighbour;
				_migr_rateIn[FEM] = _migr_rate[FEM] / factorIn4;
				_migr_rateIn[MAL] = _migr_rate[MAL] / factorIn4;
				_migr_rateOut[FEM] =
                factorOut4 ? _migr_rate[FEM] / factorOut4 : 0;
				_migr_rateOut[MAL] =
                factorOut4 ? _migr_rate[MAL] / factorOut4 : 0;
				_migr_rateCorner[FEM] = _migr_rate[FEM] / factorInCorner4;
				_migr_rateCorner[MAL] = _migr_rate[MAL] / factorInCorner4;
				_migr_rate[FEM] /= 4;
				_migr_rate[MAL] /= 4;
			}
			else {                          // 8 Neighbors
				migration_func_ptr = &LCE_Disperse::migrate_2D_ss_8Neighbour;
				_migr_rateIn[FEM] = _migr_rate[FEM] / factorIn8;
				_migr_rateIn[MAL] = _migr_rate[MAL] / factorIn8;
				_migr_rateOut[FEM] =
                factorOut8 ? _migr_rate[FEM] / factorOut8 : 0;
				_migr_rateOut[MAL] =
                factorOut8 ? _migr_rate[MAL] / factorOut8 : 0;
				_migr_rateCorner[FEM] = _migr_rate[FEM] / factorInCorner8;
				_migr_rateCorner[MAL] = _migr_rate[MAL] / factorInCorner8;
				_migr_rate[FEM] /= 8;
				_migr_rate[MAL] /= 8;
			}
			break;
            
		case 4: // 2D ss geometric
			if (_paramSet->isSet("dispersal_long_range_coef_fem")) {
				if (_paramSet->isSet("dispersal_long_range_coef_mal")) {
					_disp_long_range_coef[FEM] = _paramSet->getValue(
                                                                     "dispersal_long_range_coef_fem");
					_disp_long_range_coef[MAL] = _paramSet->getValue(
                                                                     "dispersal_long_range_coef_mal");
				}
				else error(
                           "Only one sex specific dispersal coefficient is specified: both are required!\n");
			}
			else _disp_long_range_coef[FEM] = _disp_long_range_coef[MAL] =
                _paramSet->getValue("dispersal_long_range_coef");
            
			if ((_disp_long_range_coef[FEM] < 1e-10
                 || _disp_long_range_coef[FEM] > 1.0)
                || (_popPtr->get_sexInitRatio()
                    && (_disp_long_range_coef[MAL] < 1e-10
                        || _disp_long_range_coef[MAL] > 1.0)))
				throw("Long range dispersal: The parameter 'dispersal_long_range_coef' is out of range!\n");
			migration_func_ptr = &LCE_Disperse::migrate_2D_ss_geometric;
			break;
            
		default:
			error("\nDispersal model '%i' not available!\n", _disp_model);
			break;
	}
}

//-----------------------------------------------------------------------------
// temporal_change
//-----------------------------------------------------------------------------
void LCE_Disperse::temporal_change(const unsigned int& gen)
{
	// temporal parameters
	if (_paramSet->getTemporalParams(gen) && _paramSet->updateTemporalParams(gen)){
        init(_popPtr);
    }
}

//-----------------------------------------------------------------------------
// get_neighbours
//-----------------------------------------------------------------------------
/** returns a vector with all neighbors as id of the patch */
vector<unsigned int> LCE_Disperse::get_neighbours(
                                                  const unsigned int& curPatch) {
	vector<unsigned int> vec;
	switch (_disp_model) {
		case 0:   // island
		case 1:   // island propagule
		{
			for (unsigned int i = 0; i < _nb_patch; ++i) {
				if (i != curPatch) vec.push_back(i);
			}
		}
			break;
		case 2:   // 1D SS
		{
			// left neighbor
			if (curPatch != 0)
				vec.push_back(curPatch - 1);
			else if (_border_model == 0) vec.push_back(_nb_patch - 1); // if circle
            
			// right neighbor
			if (curPatch != _nb_patch - 1)
				vec.push_back(curPatch + 1);
			else if (_border_model == 0) vec.push_back(0);          // if circle
		}
			break;
		case 3:   // 2D SS
		{
			unsigned int x_max = _x_size - 1, y_max = _y_size - 1;
			unsigned int x, y;
			ID2Coord(curPatch, x, y);
            
			// right neighbor
			if (x == x_max) {                                   // if right edge
				if (_border_model == 0) vec.push_back(Coord2ID(0, y));
			}
			else vec.push_back(Coord2ID(x + 1, y));
            
			// left neighbor
			if (x == 0) {                                        // if left edge
				if (_border_model == 0) vec.push_back(Coord2ID(x_max, y));
			}
			else vec.push_back(Coord2ID(x - 1, y));
            
			// upper neighbor
			if (y == 0) {                                       // if upper edge
				if (_border_model == 0) vec.push_back(Coord2ID(x, y_max));
			}
			else vec.push_back(Coord2ID(x, y - 1));
            
			// lower neighbor
			if (y == y_max) {                                   // if lower edge
				if (_border_model == 0) vec.push_back(Coord2ID(x, 0));
			}
			else vec.push_back(Coord2ID(x, y + 1));
            
			// if 8 neighbors the diagonals have to be added
			if (_lattice_range == 1) {
				// right neighbors (upper and lower right)
				if (x == x_max) {                				// if right edge
					if (_border_model == 0) {        // neighbours only if torus
						// upper right neighbour
						if (y == 0)
							vec.push_back(Coord2ID(0, y_max)); // if upper corner
						else vec.push_back(Coord2ID(0, y - 1)); // upper right neighbour
						// lower right neighbour
						if (y == y_max)
							vec.push_back(Coord2ID(0, 0));    // if lower corner
						else vec.push_back(Coord2ID(0, y + 1)); // upper right neighbour
					}
				}
				else {                                      // if not right edge
					// upper right neighbour
					if (y == 0) {                               // if upper edge
						if (_border_model == 0)
							vec.push_back(Coord2ID(x + 1, y_max));
					}
					else vec.push_back(Coord2ID(x + 1, y - 1)); // upper right neighbour
					// lower right neighbour
					if (y == y_max) {
						if (_border_model == 0)
							vec.push_back(Coord2ID(x + 1, 0)); // if lower corner
						else vec.push_back(Coord2ID(x + 1, y + 1)); // upper right neighbour
					}
				}
                
				// left neighbours (upper and lower left)
				if (x == 0) {                					 // if left edge
					if (_border_model == 0) {        // neighbours only if torus
						// upper right neighbour
						if (y == 0)
							vec.push_back(Coord2ID(x_max, y_max)); // if upper corner
						else vec.push_back(Coord2ID(x_max, y - 1)); // upper left neighbour
						// lower right neighbour
						if (y == y_max)
							vec.push_back(Coord2ID(x_max, 0)); // if lower corner
						else vec.push_back(Coord2ID(x_max, y + 1)); // upper left neighbour
					}
				}
				else {                                       // if not left edge
					// upper left neighbour
					if (y == 0) {                               // if upper edge
						if (_border_model == 0)
							vec.push_back(Coord2ID(x - 1, y_max));
					}
					else vec.push_back(Coord2ID(x - 1, y - 1)); // upper right neighbour
					// lower right neighbour
					if (y == y_max) {
						if (_border_model == 0)
							vec.push_back(Coord2ID(x - 1, 0)); // if lower corner
						else vec.push_back(Coord2ID(x - 1, y + 1)); // upper right neighbour
					}
				}
			}
		}
			break;
	}
    
	return vec;
}

//-----------------------------------------------------------------------------
// get_neighbours_ptr
//-----------------------------------------------------------------------------
/** returns a vector with all neighbors as pointers to the patch*/
vector<Patch*> LCE_Disperse::get_neighbours(Patch* curPatch) {
	vector<unsigned int> idx = get_neighbours(curPatch->get_ID());
    
	vector<Patch*> vec;
	vector<unsigned int>::iterator cur, end;
	for (cur = idx.begin(), end = idx.end(); cur != end; ++cur) {
		vec.push_back(_popPtr->get_vPatch(*cur));
	}
    
	return vec;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------------------------
// LCE_Disperse::execute
// ----------------------------------------------------------------------------------------
string LCE_Disperse::get_disp_model_str() {
	if (migration_func_ptr == &LCE_Disperse::migrate_zeroMigration)
		return "no migration";
	if (migration_func_ptr == &LCE_Disperse::migrate_matrix) return "matrix";
	switch (_disp_model) {
		case 0:
			return "island";
		case 1:
			return "propagule island";
		case 2:
			return "1D stepping stone";
		case 3:
			return "2D stepping stone";
            
	}
	return "not specified";
}

