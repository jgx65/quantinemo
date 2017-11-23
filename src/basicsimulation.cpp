/** @file basicsimulation.cpp
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


#include "basicsimulation.h"

//----------------------------------------------------------------------------------------
// build_component_list
// ----------------------------------------------------------------------------------------
void ComponentManager::build_component_list ()
{
	assert(_components.empty());
	_components.clear();
	list< TTraitProto* >::iterator TT_it = _TTrait_Templates.begin();
	for(; TT_it != _TTrait_Templates.end(); ++TT_it) {
		_components.push_back( (*TT_it) );
	}

	list< LCE* >::iterator LCE_it = _LCE_Templates.begin();
	for(;LCE_it != _LCE_Templates.end(); ++LCE_it) {
		_components.push_back( (*LCE_it) );
	}
}
//----------------------------------------------------------------------------------------
// get_trait_template
// ----------------------------------------------------------------------------------------
TTraitProto* ComponentManager::get_trait_template (string& name)
{
	list< TTraitProto* >::iterator TT = _TTrait_Templates.begin();
	for(; TT != _TTrait_Templates.end(); ++TT) {
		if( (*TT)->get_paramset()->getName() == name) return (*TT);
	}
	return NULL;
}
//----------------------------------------------------------------------------------------
// get_LCE_template
// ----------------------------------------------------------------------------------------
LCE* ComponentManager::get_LCE_template (string& name)
{
	list< LCE* >::iterator LCE = _LCE_Templates.begin();
	for(; LCE != _LCE_Templates.end(); ++LCE) {
		if( (*LCE)->get_event_name() == name) return (*LCE);
	}
	return NULL;
}
/*/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                             ******ParamManager*******

//----------------------------------------------------------------------------------------
// cstor
// ---------------------------------------------------------------------------------------
ParamManager::ParamManager(): _paramSetKeys(0)
{
}

// ---------------------------------------------------------------------------------------
/** destructor */
ParamManager::~ParamManager ( )
{
#ifdef _DEBUG
	message(" ParamManager::~ParamManager (%i ParamSets)\n", _allParams.size());
#endif
    clear_allParams();
    if(_paramSetKeys) delete _paramSetKeys;
}

//----------------------------------------------------------------------------------------
// clear_allParams
// ----------------------------------------------------------------------------------------
void ParamManager::clear_allParams()
{
    list<ParamSet*>::iterator cur=_allParams.begin(), end=_allParams.end();
    ++cur;
    for (; cur != end; cur++) {
     //   cout << "\n*******" << (*cur)->getName() << "..." << flush;
        //delete (*cur);
     //   cout << "done" << endl;
    }
    _allParams.clear();
}

//----------------------------------------------------------------------------------------
// build_allParams
// ----------------------------------------------------------------------------------------
void ParamManager::build_allParams ()
{
#ifdef _DEBUG
	message("ParamManager::build_allParams()\n");
#endif
    
    assert(_allParams.empty());
    _allParams.push_back(&_paramSet);
    
    list< SimComponent* >::iterator cmpt = this->_components.begin(), end = this->_components.end();
    for(; cmpt != end; ++cmpt) {
        _allParams.push_back( (*cmpt)->get_paramset() );
    }
}

//----------------------------------------------------------------------------------------
// reset_temporalParams
// ----------------------------------------------------------------------------------------
void ParamManager::reset_params ()
{
	list<ParamSet*>::iterator pos = _allParams.begin();
	for(; pos != _allParams.end(); ++pos){
		(*pos)->reset();
	}
}

//----------------------------------------------------------------------------------------
// param_consistency_check
// ----------------------------------------------------------------------------------------
bool ParamManager::param_consistency_check ()
{
	bool check = true;

	list<ParamSet*>::iterator current_paramset = _allParams.begin();
	for(;current_paramset != _allParams.end(); ++current_paramset){
		if(!(*current_paramset)->check_consistency()){
			error("ParamManager::param_consistency_check::consistency not satisfied for \"%s\"\n",
					(*current_paramset)->getName().c_str());
			check = false;
		}
	}
	return check;
}


//----------------------------------------------------------------------------------------
// set_parameters
// ----------------------------------------------------------------------------------------
/**  attribute the input to the parameters */
void
ParamManager::set_keywords (map< string,string >& simkeys, TReplicate* pRep,
                                 Metapop* pPop, bool quiet=false)
{
    _inputKeys=simkeys;
    assert(!_paramSetKeys);
    if(simkeys.empty()) return;
    
    _paramSetKeys = new ParamSet("keywords", "keywords", false, pPop);
    _paramSetKeys->set_pParamManager(this);
    map< string,string >::iterator cur=simkeys.begin(), end=simkeys.end();
    for(; cur!=end; ++cur){
        _paramSetKeys->add_param(cur->first, STR, false, my_NAN, my_NAN, "", true);
        _paramSetKeys->set_param(cur->first,cur->second, pRep);
    }
    _allParams.push_front(_paramSetKeys);
}

//----------------------------------------------------------------------------------------
// set_parameters
// ----------------------------------------------------------------------------------------
/**  attribute the input to the parameters */
bool
ParamManager::set_parameters (map< string,string >& simparams, TReplicate* pRep,
                                   bool quiet=false)
{
	_inputParams = simparams;
	list<ParamSet*>::iterator cur_paramset;
	map<string, Param*>::iterator cur_param, end_param;
	map< string,string >::iterator input;
	bool set;
	vector< string > notSetParams;
	vector< string >::iterator iterVec;
	string generalName, token, paramName, tail;
	int cur_number, number;
	unsigned int length;
	string::size_type pos;
    
    // set the ParamManager pointer (not very nice to do it here...
    for(cur_paramset = _allParams.begin(); cur_paramset != _allParams.end(); ++cur_paramset){
        (*cur_paramset)->set_pParamManager(this);
    }


	// 1: set all parameters for the parameters given in the settings file ("precise" parameters)
	//   - general parameters without numbers and parameters with a 1 are identical (eg quanti_nb_all = quanti_nb_all_1)
	for(input = _inputParams.begin(); input != _inputParams.end(); ++input) {
        string nn = input->first;
        if(input->second == "NOT_SET") continue;
        set = false;
        generalName = input->first;
        // find the corresponding parameter set
        for(cur_paramset = _allParams.begin(); cur_paramset != _allParams.end(); ++cur_paramset){
            set |= (*cur_paramset)->set_param(generalName,input->second, pRep); //is true if at least one paramset could set the param
        }
        if(set) continue;
        
        // if the parameter has an index 1 then remove it and retry
        if(generalName[generalName.size()-1]=='1' && (generalName[generalName.size()-2]=='_')){
            generalName.erase(generalName.size()-2);
            for(cur_paramset = _allParams.begin(); cur_paramset != _allParams.end(); ++cur_paramset){
                set |= (*cur_paramset)->set_param(generalName, input->second, pRep); //is true if at least one paramset could set the param
            }
        }
        if(set) continue;
        
        // if the parameter was not set add it to not set parameters
        notSetParams.push_back(input->first);    // has to be the full name, therefore input->first and not generalName!!!
	}

	// 2: check for each multiple trait if a parameter with a lower index is given:
	//   multiple traits possibilities (example parameter: quanti_nb_all_5):
	//   - for each trait a specific parameter is passed in the settings file (eg "quanti_nb_all_5")
	//   - if such a parameter is not set it is verified if a parameter with a lower number is set (eg "quanti_nb_all_3")
	//   - if such a parameter is not set it is verified if a general parameter is set (eg "quanti_nb_all"),
	for(cur_paramset = _allParams.begin(); cur_paramset != _allParams.end(); ++cur_paramset){
		// check if this parameterSET is a multiple trait
		token = (*cur_paramset)->getName();
		pos = token.rfind("_");
		if(pos == string::npos) continue;                              // "_" did not exist
		try{
			if(!isNumber<int>(token.substr(pos+1).c_str())) continue;
			number = strTo<int>(token.substr(pos+1).c_str());            // check if the last substring is a number
		}catch(BadConversion& f) {continue;                             // if the substring is not a number (don't know why this is used)
		}catch(...) {continue;}                                        // if the substring is not a number
		length = (unsigned int)token.length()-(unsigned int)pos;
		tail = token.substr(pos);

		// this paramSET is a multiple trait:
		// scroll through all the trait parameters and search the correct input parameter
		cur_param = (*cur_paramset)->getAllParams().begin();
		end_param = (*cur_paramset)->getAllParams().end();
		for(; cur_param != end_param; ++cur_param){         // for each parameter of the paramSet
			if(cur_param->second->isSet()) continue;          // if the parameter is already set continue

			// check if the parameter is grouped (i.e. the parameter has a lower number)
			// or finally no number (general parameter)
			generalName = cur_param->first;
			if(generalName.substr(generalName.size()-length) != tail) continue; // if this parameter of the set is not indexed, stop (// this parameter may not be set for each trait separately)
			generalName = generalName.substr(0, generalName.size()-length);     // get the general name without the ending number for the trait
			for(cur_number=number; cur_number>0; --cur_number){					// retest of its own due to NOT_SET
				string tt = generalName+"_"+toStr(cur_number);
				input = simparams.find(generalName+"_"+toStr(cur_number));
				if(input != simparams.end() ) break;  // is this parameter set?
			}
			if(cur_number == 0) input = simparams.find(generalName);            // not yet found try the general parameter
			if(input == simparams.end()) continue;                            // no parameter set for this parameter: stop here
			if(input->second=="NOT_SET") continue;

			// we have got the corresponding parameter: set the general parameter
			cur_param->second->set(input->second, *cur_paramset, pRep);

			// remove this parameter from the not set vector
			for(iterVec = notSetParams.begin(); iterVec != notSetParams.end(); ++iterVec){
				if(*iterVec== generalName){
					notSetParams.erase(iterVec);
					break;
				}
			}
		}
	}

    // output all not used input parameters
    if(!quiet){
        for(iterVec = notSetParams.begin(); iterVec != notSetParams.end(); ++iterVec){
            warning("ParamManager::could not set param '%s'\n",(*iterVec).c_str());
        }
    }

	return param_consistency_check();
}

//----------------------------------------------------------------------------------------
// get_paramset
// ----------------------------------------------------------------------------------------
ParamSet* ParamManager::get_paramset (string& name)
{
	list<ParamSet*>::iterator pset = _allParams.begin();
	for(; pset != _allParams.end(); ++pset) {
		if( !(*pset)->getName().compare(name)) return (*pset);
	}
	return NULL;
}

// ----------------------------------------------------------------------------------------
// setFilename
// ----------------------------------------------------------------------------------------
string ParamManager::setFilename(string& fstring,const unsigned int& sim,const unsigned& nbSims, vector<string>& args)
{
	string out, tail;
	int cur, digits;
	unsigned int next;
	unsigned int index, i, size = (unsigned int)args.size();

	// control array that all parameters were used, i.e. that the file name is unique!
	int* usedParams = new int[size];
	for(i=0; i<size; ++i){
		usedParams[i] = 0;    // set it to false
	}

	// parse the name
	cur  = 0;
	next = (unsigned int)fstring.find_first_of('%');
	while(next != (unsigned int)string::npos){
		out += fstring.substr(cur, next-cur);
		tail = fstring.substr(next+1, string::npos);
		index = (unsigned int) strtol(tail.c_str(),NULL,10);          // get the param number
		if(index > args.size()) error("too many arguments in filename!\n");
		digits = getNbDigits(index);

		out += args[index-1];
		usedParams[index-1] = 1;    // this param was used
		cur = next + digits + 1;
		next = (unsigned int)fstring.find('%',cur);  // get the next position
	}
	out += fstring.substr(cur, next-cur);

	// control if all params were used
	for(i=0; i<size; ++i){
		if(!usedParams[i]) break;
	}
	delete[] usedParams;

	if(i != size) out += "-" + toStr(sim, nbSims);    // add the simulation number if the name is not unique

	return out;
}

// ----------------------------------------------------------------------------------------
// lowercase
// ----------------------------------------------------------------------------------------
string ParamManager::lowercase(string& input)
{
	for(unsigned int i=0;i<input.size();++i){
		input[i] = tolower(input[i]);
	}
	return input;
}

// ----------------------------------------------------------------------------------------
// print2file
// ----------------------------------------------------------------------------------------
/** print all parameters to file (used for example to set up the syntax highlighting file) */
void
ParamManager::print2file(){
	ofstream FILE("parameters.txt");

	map<string, Param*>::iterator curParam, endParam;
	list<ParamSet*>::iterator curParamSet = _allParams.begin();
	for(; curParamSet != _allParams.end(); ++curParamSet){   // for each parameter set
		curParam = (*curParamSet)->getAllParams().begin();
		endParam = (*curParamSet)->getAllParams().end();
		for(; curParam != endParam; ++curParam){
			FILE << "		<string>"<< curParam->first << "</string>\n";   // for TextWrangler (Mac)
			//FILE << curParam->first << "\n";
		}
	}
}

// ----------------------------------------------------------------------------------------
// print_help
// ----------------------------------------------------------------------------------------
/** output all parameters with their description */
void
ParamManager::print_help(ostream& os, unsigned int wide1, unsigned int wide2, char fill,
                         unsigned int importance, string arg)
{
    list<ParamSet*>::iterator curP = _allParams.begin(), endP = _allParams.end();
    for(; curP!=endP; ++curP){
        (*curP)->print_help(os, wide1, wide2, fill, importance, arg);
    }
}


/*/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                             ******SimBuilder*******

//----------------------------------------------------------------------------------------
// copy cstor
// ---------------------------------------------------------------------------------------
SimBuilder::SimBuilder (const SimBuilder& SB)
{
	error("SimBuilder::SimBuilder(const SimBuilder&) not used and thus not implemented yet!\n");
}

//----------------------------------------------------------------------------------------
// get_current_trait
// ----------------------------------------------------------------------------------------
TTraitProto* SimBuilder::get_current_trait (string type)
{
	map<string, TTraitProto* >::iterator trait = _currentTraits.find(type);
	if(trait != _currentTraits.end()) return trait->second;
	return NULL;
}

//----------------------------------------------------------------------------------------
// get_current_event
// ----------------------------------------------------------------------------------------
LCE* SimBuilder::get_current_event (string& name)
{
	map< int, LCE* >::iterator LCE = _currentLifeCycle.begin();
	for(; LCE != _currentLifeCycle.end(); ++LCE) {
		if( LCE->second->get_event_name() == name) return LCE->second;
	}
	return NULL;
}

//----------------------------------------------------------------------------------------
// build_currentParams
// ----------------------------------------------------------------------------------------
bool SimBuilder::build_currentParams (map<string, string>& simparams,
                                      map<string, string>& simkeys,
                                      TReplicate* pRep, bool quiet)
{
    this->set_keywords(simkeys, pRep, NULL, quiet);
    if(! this->set_parameters(simparams, pRep, quiet) ) return false;

	//empty the current parameters container:
	//_currentParams.clear();
    assert(_currentParams.empty());

	//add all set parameters to the list of the current params:
	list<ParamSet*>::iterator current_paramset = this->_allParams.begin();
	for(; current_paramset != this->_allParams.end(); ++current_paramset){
		if((*current_paramset)->isSet()) {
#ifdef _DEBUG
            message("SimBuilder::build_currentParams:adding paramset %s\n",(*current_paramset)->getName().c_str());
#endif
			_currentParams.push_back(*current_paramset);
		}
	}
	return true;
}

//----------------------------------------------------------------------------------------
// build_currentTraits
// ----------------------------------------------------------------------------------------
map< string,TTraitProto* >& SimBuilder::build_currentTraits()
{
	//build the list of Traits for this simulation:
	assert(_currentTraits.empty());

	list< TTraitProto* >::iterator trait = this->_TTrait_Templates.begin();
	for(; trait != this->_TTrait_Templates.end(); ++trait) {
		if( (*trait)->get_paramset()->isSet() ){
			_currentTraits[(*trait)->get_type_index()] = (*trait);
		}
	}
	return _currentTraits;
}

//----------------------------------------------------------------------------------------
// build_currentLifeCycle
// ----------------------------------------------------------------------------------------
//build the list of the life cycle events for this simulation:
map< int, LCE* >& SimBuilder::build_currentLifeCycle()
{
    assert(_currentLifeCycle.empty());

	list< LCE* >::iterator LCEs = this->_LCE_Templates.begin();
	for(;LCEs != this->_LCE_Templates.end();++LCEs){
		_currentLifeCycle[(int)(*LCEs)->get_rank()] = (*LCEs);
	}
	return _currentLifeCycle;
}


//----------------------------------------------------------------------------------------
// checkLCEconsistency
// ----------------------------------------------------------------------------------------
/** check if the sequence of LCE makes sense */
bool
SimBuilder::checkLCEconsistency()
{
	// check if a LCE is really used
	map< int, LCE* >::iterator iterLCE = _currentLifeCycle.begin();
	for(; iterLCE != _currentLifeCycle.end();){
		string name = iterLCE->second->get_event_name();
		if(!iterLCE->second->get_paramset()->isSet()){
			_currentLifeCycle.erase(iterLCE++);
		}
		else  ++iterLCE;
	}

	//find the age class required to start the life cycle with:
	age_t availableAge = NONE;
	iterLCE = _currentLifeCycle.begin();
	for(; iterLCE != _currentLifeCycle.end() && availableAge == NONE; ++iterLCE){
		availableAge = iterLCE->second->requiredAgeClass();
	}

	// check if the LCE makes sense
	LCE* curLCE, *lastLCE = NULL;
	bool stop = false;
	iterLCE = _currentLifeCycle.begin(); // start with the first LCE
	while(!stop){
		// check if the last and the first LCE correspond
		if(iterLCE == _currentLifeCycle.end()){
			stop = true;
			iterLCE = _currentLifeCycle.begin();
		}

		// get the LCE
		curLCE = iterLCE->second;

		// if a age class is required
		if(curLCE->requiredAgeClass()){
			// check if the required age class is available or if no age classes are needed
			if(!(availableAge & curLCE->requiredAgeClass())){
				string name = curLCE->get_event_name();
				error("LCE order does not make sense (%s -> %s)\n",
						lastLCE->get_event_name().c_str(), curLCE->get_event_name().c_str());
			}
			lastLCE = curLCE;
		}

		// update the current age class
		availableAge &= ~curLCE->removeAgeClass();
		availableAge |= curLCE->addAgeClass();
		++iterLCE;
	}
	return true;
}
