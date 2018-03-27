/** @file simcomponent.h
 *
 *   Copyright (C) 2006 Frederic Guillaume    <guillaum@zoology.ubc.ca>
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

#ifndef simcomponentH
#define simcomponentH

#include "fileservices.h"
#include "statservices.h"
#include "param.h"

/**Interface to all basic components of a simulation (traits, life cycle events, pop, etc.\ ).
 * Implements the interface to handle the files and stats handler and the parameters container.
 * Contains the components parameters.
 */

class TSimComponent {
protected:
	/**The parameters container. */
	ParamSet* _paramSet;  // has to be deleted

public:

	TSimComponent() :
			_paramSet(0) {
	}
	//**Dstor. Deletes the parameter container. */
	virtual ~TSimComponent() {
		if (_paramSet) delete _paramSet;
	}

	///@name Services interface
	///@{
	/**Loads the component's FileHandler onto the FileServices.
	 @param loader the file service
	 **/
	virtual void loadFileServices(FileServices* loader) = 0;

	/**Loads the component's StatHandler onto the StatServices.
	 @param loader the stat service
	 **/
	virtual void loadStatServices(StatServices* loader) = 0;
	///@}

	/**Sets a new ParamSet and name it.
	 @param name the name of the parameter container
	 **/
	///@name Parameters handling
	///@{
	virtual void set_paramset(string name, string name_long, bool isRequired, TMetapop* ptr) {
		if (_paramSet) delete _paramSet;
		_paramSet = new ParamSet(name, name_long, isRequired, ptr);
	}

	/**ParamSet accessor. **/
	virtual ParamSet* get_paramset() {
		return _paramSet;
	}

	/**Interface to add a parameter to the set.
	 * @param param The parameter to add to the set
	 */
	virtual void add_parameter(Param* param) {
		_paramSet->add_param(param);
	}

	/**Interface to add a parameter to the set.
	 * @param Name The param string name as read in the init file
	 * @param Type The type of the argument (DBL=double, INT=integer, BOOL=boolean, STR=string, etc., see type.h)
	 * @param isRequired True if the parameter is mandatory
	 * @param low_bnd The lower value the argument can take
	 * @param up_bnd The upper value the argument can take
	 * @param default value
	 * @param temp can the parameter have temporal changes?
	 */
	virtual void add_parameter(string Name, param_t Type, bool isRequired,
			double low_bnd, double up_bnd, string def, bool temp = false,
                               string text="", unsigned int importance=5) {
		_paramSet->add_param(Name, Type, isRequired, low_bnd, up_bnd, def,
				temp, text, importance);
	}

	/**Param getter. Asks the ParamSet to return a pointer to Param \a name.
	 @see ParamSet::get_param
	 */
	virtual Param* get_parameter(string name) {
		return _paramSet->get_param(name);
	}
	virtual double get_parameter_value(string name) {
		return _paramSet->getValue(name);
	}
	virtual bool get_parameter_isSet(string name) {
		return _paramSet->isSet(name);
	}
	virtual bool get_parameter_isMatrix(string name) {
		return _paramSet->isMatrix(name);
	}
	virtual TMatrix* get_parameter_matrix(string name) {
		return _paramSet->getMatrix(name);
	}
	///}

};

#endif  //SIMCOMPONENT_H
