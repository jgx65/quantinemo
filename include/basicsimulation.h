/** @file basicsimulation.h
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

#ifndef basicsimulationH
#define basicsimulationH

#include <list>
#include <map>

#include "ttrait.h"
#include "lifecycleevent.h"
using namespace std;

/*/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/
/** Class to manage the simulation components and prototypes.
 *  This class stores and provides accessors to the simulation components. It also stores
 *  the trait prototype and life cycle event template instances.
 **/
/*\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\*/

class ComponentManager {

private:

public:
	ComponentManager() {
#ifdef _DEBUG
		message(" ComponentManager()\n");
#endif
	}
	~ComponentManager() {
		reset_sim_components();
	}
	/**Clears and builds the list of all components from the lists of trait prototype templates and life cycle event templates. **/
	void build_component_list();

	/**Push a component at the back of the component list. **/
	void add_component(TSimComponent* cmpt) {
		_components.push_back(cmpt);
	}

	/**Add a trait prototype to the template and component lists. **/
	void add_trait_template(TTraitProto* trait) {
		_TTrait_Templates.push_back(trait);
		_components.push_back(trait);
	}

	/**Reset all trait prototype of the template and component lists. **/
	void reset_sim_components() {
       list<TTraitProto*>::iterator pos2 = _TTrait_Templates.begin();
		for (; pos2 != _TTrait_Templates.end(); ++pos2) {
 			if (*pos2) delete *pos2;
		}
		_TTrait_Templates.clear();

        list<LCE*>::iterator pos = _LCE_Templates.begin();
        for (; pos != _LCE_Templates.end(); ++pos) {
            if (*pos) delete *pos;
        }
        _LCE_Templates.clear();

		_components.clear();
	}

	/**Add a life cycle event to the template and component lists. **/
	void add_LCE_template(LCE* event) {
		_LCE_Templates.push_back(event);
		_components.push_back(event);
	}

	/**Search for component with "name" in the trait prototype list.
	 * @return NULL if component with "name" not found, a pointer to that component otherwise
	 **/
	TTraitProto* get_trait_template(string& name);

	/**Search for component with "name" in the life cycle events list.
	 * @return NULL if component with "name" not found, a pointer to that component otherwise
	 **/
	LCE* get_LCE_template(string& name);

protected:

	/**List of all the simulation components. **/
	list<TSimComponent*> _components;			// don't delete elements

	/**List of all trait prototypes of the simulation, a subset of _components list. */
	list<TTraitProto*> _TTrait_Templates;		// delete elements

	/**List of all the life-cycle events of the simulation, a subset of _components list. */
	list<LCE*> _LCE_Templates;					// delete elements
};

/*/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/
/**Class to manage the sets of parameters of the simulation components.
 *  This class performs parameters setting and checking for the whole set of the simulation
 *  components. Provides access to derived classes to the complete list of parameter sets.
 *  Also sets the list of simulations parameters in case of sequential parameters found in
 *  input.
 *  It stores and builds the simulation parameters set.
 *
 *  @see ParamsParser and derived classes. These are the input parameters providers.
 **/
/*\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\*/
class ParamManager: public ComponentManager {
public:
	/**Cstor. Builds the simulation PramaSet. **/
	ParamManager();

	~ParamManager(); // {/*message("ParamManager::~ParamManager\n");*/}

	/**Adds a ParamSet to the list of the parameter sets of the simulation. **/
	void add_paramset(ParamSet* paramset) {
        paramset->set_pParamManager(this);
        _allParams.push_back(paramset);
	}

	/**Looks for paramset with "name" in the list of parameter sets.
	 * @return NULL if parameter set "name" not found
	 **/
	ParamSet* get_paramset(string& name);

	ParamSet& get_paramset() {
		return _paramSet;
	}  // the current paramset

	/**Clears and fills the _allParams list with the ParamSet's of the simulation components. **/
	void build_allParams();
    void clear_allParams();

	/**Accessor of the whole list of the parameter sets.
	 *@return the reference to the list of ParamSet
	 **/
	list<ParamSet*>& get_allParams() {
		return _allParams;
	}

	/**Sets the parameters of the simulation with the argument strings given in input.
	 *Scans the _allParams list to set the parameters present in the input map simparams. Each ParamSet checks
	 *internally for the presence of a Param with the given name string and sets its value with the given argument, if present.

	 \b Note: all ParamSet owning a Param with the same name will use the same argument string. The input map is not a
	 * multimap, each param name is present only once.
	 *
	 *@param simparams a map containing the parameter names and their argument string
	 *@return the status of the \c param_consistency_check() function
	 *@callgraph
	 **/
    bool set_parameters(map<string, string>& simparams, TReplicate* pRep, bool quiet);
    void set_keywords(map<string, string>& simkeys, TReplicate* , TMetapop* pPop,  bool quiet);
//	void check_parameter_names(map<string, vector<string> >& inputParam);
//	void check_parameter_names(map<string, string>& inputParam);

	/**Checks if all the mandatory parameters are set so that the simulation can be launched.
	 *@return TRUE if all \c ParamSet::check_consistency() returned true
	 *@callgraph
	 **/
	bool param_consistency_check();

	/**Builds the list of simulation parameters from the parsed input file(s). @callgraph**/
	void build_records(map<string, vector<string> >& initParams);

	/**Accessor to the simulations parameter list. **/
	list<map<string, string> >& get_simRecords() {
		return _simRecords;
	}

	/**Accessor to the first element in the simulations parameter list. **/
	map<string, string>& get_firstRecord() {
		return (*_simRecords.begin());
	}

	/**Accessor to the size of the simulations parameter list, i.e. the number of simulations to perform. **/
	int get_nbSims() {
		return (int)_simRecords.size();
	}

	void reset_params();
    
    ParamSet* get_paramSetKeys(){return _paramSetKeys;}

    /**The ParamSet param set of the simulation. **/
 //   ParamSet _paramSet;

    void print_help(ostream& os, unsigned int wide1=20, unsigned int wide2=60,
                    char fill='.', unsigned int importance=5, string arg="");
    void print_help(ParamSet* pParamSet, ostream& os, unsigned int wide1=20, unsigned int wide2=60,
                    char fill='.', unsigned int importance=5, string arg="");
    

protected:
	/**A list of all the parameter sets of all the simulation components loaded in the _component list of the ComponentManager. **/
	list<ParamSet*> _allParams;

	/**A map of the parameters and their arguments of the current (running) simulation. **/
    map<string, string> _inputParams;  //multimap?

	/**Sets of parameters of all the simulations to perform. **/
	list<map<string, string> > _simRecords;
    
    map<string, string> _allKeywords; // keywords of the simulation

    /**The ParamSet param set of the simulation. **/
    ParamSet                     _paramSet;
    
    /** keys used in this simulation */
    map<string, string> _inputKeys;
    ParamSet*           _paramSetKeys;

private:

	string setFilename(string& fstring, const unsigned int& sim,
			const unsigned& nbSims, vector<string>& args);
	string lowercase(string& input);

	void print2file();	// print all parameter names to file

};


/*/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/
/**Provides methods to build the user's selected set of life cycle events and traits from the parameters.
 *  This class implements methods to build the lists of selected traits and life cycle events from the user's defined parameters.
 *  Each simulation component that has its ParamSet in the "set" state is elligible to be part of the current simulation.
 *  Accessors to these components are provided. This class does however not provide a runnable simulation object.
 **/
/*\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\*/
class SimBuilder: public ParamManager {
public:
	SimBuilder() {
#ifdef _DEBUG
		message(" SimBuilder::SimBuilder\n");
#endif
	}

	/**copy cstor.*/
	SimBuilder(const SimBuilder& SB);

	~SimBuilder() {
#ifdef _DEBUG
		message(" SimBuilder::~SimBuilder\n");
#endif
	}
	/**Builds the list of parameters from user's defined input parameters.
	 * @param simparams Hashtable of the parsed input parameters
	 * @return TRUE if parameters are consitently set
	 **/
	bool build_currentParams(map<string, string>& simparams,
                             map<string, string>& simkeywords,
                             TReplicate* pRep, bool quiet=false);

	/**Selects the trait prototypes that have their parameters set.
	 * @return Hashtable of the trait prototypes
	 **/
	map<string, TTraitProto*>& build_currentTraits();

	/**Selects the life cycle events that have their parameters set.
	 * @return Hashtable of the current life cycle events.
	 **/
	map<int, LCE*>& build_currentLifeCycle();

	/**Accessor to the list of current trait prototypes.
	 * @param type the trait type
	 * @return ptr to the prototype
	 **/
	TTraitProto* get_current_trait(string type);

	/**Accessor to the list of current LCEs.
	 * @param name the name of the LCE
	 * @return ptr to the LCE
	 **/
	LCE* get_current_event(string& name);

	/**Accessor to the list of the selected parameter sets.
	 * @return list of ParamSet
	 **/
	list<ParamSet*>& get_currentParams() {
		return _currentParams;
	}

	bool checkLCEconsistency();

protected:
	/**List of the selected simulation components from the user defined input parameters. **/
	list<ParamSet*> _currentParams;

	/**List of the selected trait prototypes from the user defined input parameters. **/
	map<string, TTraitProto*> _currentTraits;

	/**List of the selected life cycle events from the user defined input parameters. **/
	map<int, LCE*> _currentLifeCycle;

};
#endif

