/** @file tselectionneutral.cpp
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


#include "tselectiontype.h"
#include "ttquanti.h"
#include "tmetapop.h"

// ----------------------------------------------------------------------------------------
// ini_neutral_selection
// ----------------------------------------------------------------------------------------
/** set the selection optima and intensity of the patches */
void
TSelectionNeutral::init()
{
	_nb_selection_params = 0;

	// selection prussure does not change between generations
	_get_selection_pressure_func_ptr = NULL;
	_get_selection_pressure_tot_func_ptr = &TSelectionTrait::_get_selection_pressure_tot_const;
}

// ----------------------------------------------------------------------------------------
// ini_stabilizing_selection
// ----------------------------------------------------------------------------------------
/** set the selection optima and intensity of the patches */
void
TSelectionStabilizing::init()
{
	_nb_selection_params = 2;

	// get the selection variance values
	if(!_selection_sd) _selection_sd = new double[_nb_selection_params];
	_selection_sd[0] = _pSel->get_optima_sd();           // does the optimum varies from generation to generation?
	_selection_sd[1] = _pSel->get_intensity_sd();        // does the intensity varies from generation to generation?

	// check if at least one param has fluctuations between generations (i.e. is not zero)
	unsigned int i;
	for(i=0; i<_nb_selection_params; ++i){
		if(_selection_sd[i]) break;
	}

	// check if all are const, i.e. do not change between generations
	if(i != _nb_selection_params){
		if(!_get_selection_pressure_func_ptr){
			_get_selection_pressure_func_ptr = new _func_ptr[_nb_selection_params];
			_selection_pressure[FEM]         = new double[_nb_selection_params];
			_selection_pressure[MAL]         = new double[_nb_selection_params];
		}

		_get_selection_pressure_tot_func_ptr = &TSelectionTrait::_get_selection_pressure_tot_var;

		// set the correct function pointer
		for(i=0; i<_nb_selection_params; ++i){
			if(_selection_sd[i]) _get_selection_pressure_func_ptr[i] = &TSelectionTrait::_get_selection_pressure_var;
			else                 _get_selection_pressure_func_ptr[i] = &TSelectionTrait::_get_selection_pressure_const;
		}
	}
	else{           // no changes between generations
		if(_get_selection_pressure_func_ptr){
			delete[] _get_selection_pressure_func_ptr; _get_selection_pressure_func_ptr = NULL;
			delete[] _selection_pressure[FEM];         _selection_pressure[FEM]=NULL;
			delete[] _selection_pressure[MAL];         _selection_pressure[MAL]=NULL;
		}
		delete[] _selection_sd; _selection_sd=NULL;
		_get_selection_pressure_tot_func_ptr = &TSelectionTrait::_get_selection_pressure_tot_const;
	}
}

// ----------------------------------------------------------------------------------------
// ini_directional_selection
// ----------------------------------------------------------------------------------------
void
TSelectionDirectional::init()
{
	_nb_selection_params = 5;

	// get the selection variance values
	if(!_selection_sd) _selection_sd = new double[_nb_selection_params];
	_selection_sd[0] = _pSel->get_min_sd();           // does the lower asymptote varies from generation to generation?
	_selection_sd[1] = _pSel->get_max_sd();           // does the upper asymptote varies from generation to generation?
	_selection_sd[2] = _pSel->get_max_growth_sd();    // does the maximal growth varies from generation to generation?
	_selection_sd[3] = _pSel->get_growth_rate_sd(); 	// does the growth rate varies from generation to generation?
	_selection_sd[4] = _pSel->get_symmetry_sd();      // does the symmetry varies from generation to generation?

	// check if at least one param has fluctuations between generations (i.e. is not zero)
	unsigned int i;
	for(i=0; i<_nb_selection_params; ++i){
		if(_selection_sd[i]) break;
	}

	// check if all are const, i.e. do not change between generations
	if(i != _nb_selection_params){
		if(!_get_selection_pressure_func_ptr){
			_get_selection_pressure_func_ptr = new _func_ptr[_nb_selection_params];
			_selection_pressure[FEM]         = new double[_nb_selection_params];
			_selection_pressure[MAL]         = new double[_nb_selection_params];
		}

		_get_selection_pressure_tot_func_ptr = &TSelectionTrait::_get_selection_pressure_tot_var;

		// set the correct function pointer
		for(i=0; i<_nb_selection_params; ++i){
			if(_selection_sd[i]) _get_selection_pressure_func_ptr[i] = &TSelectionTrait::_get_selection_pressure_var;
			else                 _get_selection_pressure_func_ptr[i] = &TSelectionTrait::_get_selection_pressure_const;
		}
	}
	else{           // no changes between generations
		if(_get_selection_pressure_func_ptr){
			delete[] _get_selection_pressure_func_ptr; _get_selection_pressure_func_ptr = NULL;
			delete[] _selection_pressure[FEM];         _selection_pressure[FEM]=NULL;
			delete[] _selection_pressure[MAL];         _selection_pressure[MAL]=NULL;
		}
		delete[] _selection_sd; _selection_sd=NULL;
		_get_selection_pressure_tot_func_ptr = &TSelectionTrait::_get_selection_pressure_tot_const;
	}
}

//-----------------------------------------------------------------------------
// _get_fitness_neutral
//-----------------------------------------------------------------------------
/** if no selection acts we return a fitness of 1 (function should never be called) */
double
TSelectionNeutral::get_fitness_none()
{
 return 1;
}

//-----------------------------------------------------------------------------
// get_fitness_stabilizing:
//-----------------------------------------------------------------------------
/** standard stabilizing function
  * w = exp(-[(vp-optima)/omega]^2 / 2))
  * intensity is the sd
	* special is the weight when multiple traits are used
  */
double
TSelectionStabilizing::get_fitness_none()
{
  double value = (get_phenotype()-_selection_pressure[_curSex][0])/_selection_pressure[_curSex][1];
  value *= value;
  return exp(value/(-2));
}

//-----------------------------------------------------------------------------
// get_fitness_directionl
//-----------------------------------------------------------------------------
/** directional selection with a general logistic curve.
  *  A: the lower asymptote;
  *   C: the upper asymptote;
  *   M: the time of maximum growth;
  *   B: the growth rate;
  *   T: affects near which asymptote maximum growth occurs.
  */
double
TSelectionDirectional::get_fitness_none()
{
	return generalLogisticCurve(get_phenotype(),
			_selection_pressure[_curSex][0],  // min
			_selection_pressure[_curSex][1],  // max
			_selection_pressure[_curSex][2],  // max growth
			_selection_pressure[_curSex][3],  // growth rate
			_selection_pressure[_curSex][4]); // symmetry
}

// ----------------------------------------------------------------------------------------
// TSelectionFitnessLandscape
// ----------------------------------------------------------------------------------------
void
TSelectionFitnessLandscape::init()
{
	_nb_selection_params = 0;

	// selection prussure does not change between generations
	_get_selection_pressure_func_ptr = NULL;
	_get_selection_pressure_tot_func_ptr = &TSelectionTrait::_get_selection_pressure_tot_const;
}

//-----------------------------------------------------------------------------
// TSelectionFitnessLandscape
//-----------------------------------------------------------------------------
/** selection pressure is defined by a fitness landscape
	* _selection_pressure[SEX][0]:               number of datapoints
	* _selection_pressure[SEX][1 -> size]:       phenotype array;
	* _selection_pressure[SEX][size+1 ->2*size]: phenotype array;
	*/
double
TSelectionFitnessLandscape::get_fitness_none()
{
	unsigned int size  = (unsigned int)_selection_pressure[_curSex][0];
	double* aPhenotype = _selection_pressure[_curSex] + 1 ;        // sorted phenotpyes
	double* aFitness   = _selection_pressure[_curSex] + 1 + size;  // corresponding fitnesses
	double pheno = get_phenotype();
	unsigned int pos = ARRAY::searchBest(aPhenotype, size, pheno); // find the same or next bigger phenotype

	// out of range?
	if(pos==my_NAN){
		if(pheno<aPhenotype[0]) return aFitness[0];       					 // phenotype too small
		else                    return aFitness[size-1];             // phenotype too large
	}


	double upper_pheno = aPhenotype[pos];
	if(pheno==upper_pheno) return aFitness[pos];                   // precise fit?

	// linear interpolation  (pos==0 is not possible...)
	assert(pheno < upper_pheno);
	double lower_pheno = aPhenotype[pos-1];
	return aFitness[pos-1]+(pheno-lower_pheno)/(upper_pheno-lower_pheno)*(aFitness[pos]-aFitness[pos-1]);
}

// ----------------------------------------------------------------------------------------
// ini_stabilizing_selection
// ----------------------------------------------------------------------------------------
/** set the selection coefficient of the patches */
void
TSelectionSelectionCoefficient::init()
{
	_nb_selection_params = 2;

	// get the selection variance values
	if(!_selection_sd) _selection_sd = new double[_nb_selection_params];
	_selection_sd[0] = _pSel->get_selCoef_sd();           // does the fitenss for AA varies from generation to generation?
	_selection_sd[1] = _pSel->get_selCoef_sd();           // does the selection coeffcient varies from generation to generation?

	// check if at least one param has fluctuations between generations (i.e. is not zero)
	unsigned int i;
	for(i=0; i<_nb_selection_params; ++i){
		if(_selection_sd[i]) break;
	}

	// check if all are const, i.e. do not change between generations
	if(i != _nb_selection_params){
		if(!_get_selection_pressure_func_ptr){
			_get_selection_pressure_func_ptr = new _func_ptr[_nb_selection_params];
			_selection_pressure[FEM]         = new double[_nb_selection_params];
			_selection_pressure[MAL]         = new double[_nb_selection_params];
		}

		_get_selection_pressure_tot_func_ptr = &TSelectionTrait::_get_selection_pressure_tot_var;

		// set the correct function pointer
		for(i=0; i<_nb_selection_params; ++i){
			if(_selection_sd[i]) _get_selection_pressure_func_ptr[i] = &TSelectionTrait::_get_selection_pressure_var;
			else                 _get_selection_pressure_func_ptr[i] = &TSelectionTrait::_get_selection_pressure_const;
		}
	}
	else{           // no changes between generations
		if(_get_selection_pressure_func_ptr){
			delete[] _get_selection_pressure_func_ptr; _get_selection_pressure_func_ptr = NULL;
			delete[] _selection_pressure[FEM];         _selection_pressure[FEM]=NULL;
			delete[] _selection_pressure[MAL];         _selection_pressure[MAL]=NULL;
		}
		delete[] _selection_sd; _selection_sd=NULL;
		_get_selection_pressure_tot_func_ptr = &TSelectionTrait::_get_selection_pressure_tot_const;
	}
}

//-----------------------------------------------------------------------------
// get_fitness_selCoef:
//-----------------------------------------------------------------------------
/** standard selection coefficient function
	* AA 1          X
	* Aa 1-hs       X-hs
	* aa 1-s        X-s
	*/
double
TSelectionSelectionCoefficient::get_fitness_none()
{
	//The genotype is computed using the formula 2[(1−h)a_1 +h*a_2] with a_1 and a_2 either 0 or 1, so the only missing piece is s
    double fitness=_selection_pressure[_curSex][0] - get_genotype()*_selection_pressure[_curSex][1]/2;

	// is there environmental variance
	double sd = (this->*func_ptr_get_sdVe)();
	if(sd) fitness += sd*_popPtr->rand().Normal();

	if(fitness<0)return 0;       // if fitness is negative return 0
	return fitness;
}

//-----------------------------------------------------------------------------







