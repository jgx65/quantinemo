/** @file param.cpp
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

#include "param.h"

#include "basicsimulation.h"
#include "tmetapop.h"
#include <sstream>
#include <iostream>
#include <string>
#include "treplicate.h"
#include "tsimulation.h"
using namespace std;



// ----------------------------------------------------------------------------------------
// Param
// ----------------------------------------------------------------------------------------
Param::Param (string& Name,param_t Type,bool mandatory,double low_bnd,
		double up_bnd, string def, bool temp, string text, unsigned int importance)
: _name(Name),_arg(def),_default_arg(def),_type(Type),_isSet(0),_isRequired(mandatory),
_parSet(NULL), _temporalParamAllowed(temp), _isRandom(false), _description(text),
_importance(importance)
{
	_bounds[0] = low_bnd;
	_bounds[1] = up_bnd;

	_macros_func_ptr[0] = &Param::rep;           	_macros[0] = "rep";
	_macros_func_ptr[1] = &Param::seq;           	_macros[1] = "seq";
	_macros_func_ptr[2] = &Param::runif;         	_macros[2] = "runif";
	_macros_func_ptr[3] = &Param::rnorm;         	_macros[3] = "rnorm";
	_macros_func_ptr[4] = &Param::rgamma;        	_macros[4] = "rgamma";
	_macros_func_ptr[5] = &Param::rbeta;         	_macros[5] = "rbeta";
	_macros_func_ptr[6] = &Param::rpois;         	_macros[6] = "rpois";
	_macros_func_ptr[7] = &Param::rlnorm;        	_macros[7] = "rlnorm";
	_macros_func_ptr[8] = &Param::rbinom;        	_macros[8] = "rbinom";
	_macros_func_ptr[9] = &Param::rsample;       	_macros[9] = "rsample";
	_macros_func_ptr[10]= &Param::seq2D;         	_macros[10]= "seq2D";
	_macros_func_ptr[11]= &Param::seq2Db;        	_macros[11]= "seq2Db";
	_macros_func_ptr[12]= &Param::product;       	_macros[12]= "product";
	_macros_func_ptr[13]= &Param::round;          _macros[13]= "round";
	_macros_func_ptr[14]= &Param::ceil;           _macros[14]= "ceil";
	_macros_func_ptr[15]= &Param::floor;          _macros[15]= "floor";
	_macros_func_ptr[16]= &Param::trunc;          _macros[16]= "trunc";
	_macros_func_ptr[17]= &Param::logistic;       _macros[17]= "logistic";
}

// ----------------------------------------------------------------------------------------
// copy constructor
// ----------------------------------------------------------------------------------------
Param::Param(const Param& i)
{
    _name = i._name;
    _arg =i._arg;
    _tot_arg = i._tot_arg;
    _default_arg = i._default_arg;
    _type = i._type;
    _isSet = i._isSet;
    _isRequired = i._isRequired;
    _bounds[0] = i._bounds[0];
    _bounds[1] = i._bounds[1];
    _parSet = i._parSet;
    _temporalParamAllowed = i._temporalParamAllowed;
    _temporalArgs = i._temporalArgs;
    _isRandom = i._isRandom;
    _tot_arg_no_macro = i._tot_arg_no_macro;
    _description = i._description;
    _importance = i._importance;
    
    for (unsigned int s=0; s<NB_MACROS; ++s){
        _macros[s] = i._macros[s];
        _macros_func_ptr[s] = i._macros_func_ptr[s];
    }

}

// ----------------------------------------------------------------------------------------
// update_arg
// ----------------------------------------------------------------------------------------
/** the argument is updated by the new temporal argument */
void
Param::update_arg(const int& gen)
{
	assert(_temporalArgs.find(gen) != _temporalArgs.end());
    set_new_arg(_temporalArgs.find(gen)->second);
}

// ----------------------------------------------------------------------------------------
// set_arg_for_gen
// ----------------------------------------------------------------------------------------
/** if temporal parameter the argument is updated for the given generation
 * if arg ist not set for this generation then search backwards
 */
void
Param::set_arg_for_gen(int gen)
{
	if(_temporalArgs.empty()) return;            // if it is not temporal
	map<int, string>::iterator pos;
	do{
		pos = _temporalArgs.find(gen);             // check if the generation is set
		--gen;
	}while(pos == _temporalArgs.end());           /// stop if found (latest at generation 1)

	_arg = pos->second;
}

// ----------------------------------------------------------------------------------------
// reset
// ----------------------------------------------------------------------------------------
void
Param::reset ()
{
	_arg = _default_arg;
	_isSet = false;
	_temporalArgs.clear();
}

// ----------------------------------------------------------------------------------------
// set
// ----------------------------------------------------------------------------------------
/** control if the value is within the bounds */
bool
Param::check_arg (string arg)
{
    //integer or decimal:
    if(_type == INT2 || _type == DBL || _type == INT_MAT || _type == DBL_MAT) {
        //!!a matrix may also be specified!! in that case, we don't check values here
        if(arg[0] != '{') {
            double val;
            try{
                val = strTo<double>(arg);
            }
            catch(...){
                error("The argument of parameter '%s' should be a number (arg: '%s')!\n", _name.c_str(), arg.c_str());
            }
            
           // if the value is outside the bounds
            if((val < _bounds[0] && _bounds[0]!= my_NAN) ||
               (val > _bounds[1] && _bounds[1]!= my_NAN)) {
                if(_type == INT2 || _type == INT_MAT){
                    if(_bounds[0]== my_NAN){
                        error("The value of parameter '%s' is outside of the bounds (value: %i; bounds ] - %i])!\n",
                              _name.c_str(), (int)val, (int)_bounds[1]);
                    }
                    else if(_bounds[1]== my_NAN){
                        error("The value of parameter '%s' is outside of the bounds (value: %i; bounds [%i - [)!\n",
                              _name.c_str(), (int)val, (int)_bounds[0]);
                    }
                    else{
                        error("The value of parameter '%s' is outside of the bounds (value: %i; bounds [%i - %i])!\n",
                              _name.c_str(), (int)val, (int)_bounds[0], (int)_bounds[1]);
                    }
                }
                else{
                    if(_bounds[0]== my_NAN){
                        error("The value of parameter '%s' is outside of the bounds (value: %f; bounds ] - %f])!\n",
                              _name.c_str(), val, _bounds[1]);
                    }
                    else if(_bounds[1]== my_NAN){
                        error("The value of parameter '%s' is outside of the bounds (value: %f; bounds [%f - [)!\n",
                              _name.c_str(), val, _bounds[0]);
                    }
                    else{
                        error("The value of parameter '%s' is outside of the bounds (value: %f; bounds [%f - %f])!\n",
                              _name.c_str(), val, _bounds[0], _bounds[1]);
                    }
                }
                return false;
            }
        }
        else if(_type == INT2 || _type == DBL){
            error("The value of parameter '%s' should be a single number, not a matrix!\n",
                  _name.c_str());
        }
    }
    else if(_type == MAT || _type == MAT_VAR || _type == INT_MAT || _type == DBL_MAT) {
        if(arg[0] != '{') {
            error("The argument of parameter '%s' should be a matrix!\n", _name.c_str());
            return false;
        }
    }
    return true;
}

// ----------------------------------------------------------------------------------------
// set
// ----------------------------------------------------------------------------------------
/** first setting of the parameter argument */
bool
Param::set(string arg, ParamSet* parSet, TReplicate* pRep)
{
    _tot_arg = arg;
	_parSet  = parSet;
    _pReplicate = pRep;
    
    set_up_main_random_param(arg);
    
    // temporal argument?
    arg = set_temporal_args(arg);
    
    bool ok = set_new_arg(arg);
    
    set_up_main_random_param(arg);

    return ok;
}

// ----------------------------------------------------------------------------------------
// set_up_main_radnom_param
// ----------------------------------------------------------------------------------------
/** if it is a main parameter (_current_generation==my_NAN and is:
 * - random numbers are identical among replicates (parameter _randomPerReplicate)
 * - is temporal
 * store the common rundom numbers in _common_random_number
 */
void
Param::set_up_main_random_param(string arg)
{
    if(!_isRandom) return;
    assert(_pReplicate);
    if(_pReplicate->_current_replicate!=my_NAN) return;
    if(_pReplicate->get_pSimulation()->_randomPerReplicate) return;
    if(_temporalArgs.empty()) return;
    
    
    // ok, we have to store the args for the different times.
    map<int, string>::iterator cur=_temporalArgs.begin(), end=_temporalArgs.end();
    for(; cur!=end; ++cur){
        _temporalArgs[cur->first] = get_new_arg(cur->second);
    }

}

// ----------------------------------------------------------------------------------------
// set
// ----------------------------------------------------------------------------------------
/** set the current arg at any time (used for temporal arguments)
 */
string
Param::get_new_arg(string arg)
{
    // if the string is enclosed with "" remove them
    if(arg[0] == '\"'){
        if(arg[arg.length()-1] != '\"') error("Parameter %s has a problem with quotes!\n", _name.c_str());
        arg = arg.substr(1, arg.length()-2);
    }
    
    // keywords?
    arg = replaceKeywords(arg);
    
    // macros?
    return replaceMacros(arg);
}

// ----------------------------------------------------------------------------------------
// set
// ----------------------------------------------------------------------------------------
/** set the current arg at any time (used for temporal arguments)
 */
bool
Param::set_new_arg(string arg)
{
    _arg = get_new_arg(arg);
    
    // check current argument
    if(!check_arg(_arg)) return false;
    
    _isSet = true;
    return true;
}

// ----------------------------------------------------------------------------------------
// replaceKeywords
// ----------------------------------------------------------------------------------------
/** replaces all occurences of keywords by the argument of the keyword
 */
string
Param::replaceKeywords(string arg)
{
    assert(_parSet),
    assert(_parSet->get_pParamManager());
    if(!_parSet->get_pParamManager()->get_paramSetKeys()) return arg;   // there are no kewords
    
    map<string, Param*>* pKeys=&_parSet->get_pParamManager()->get_paramSetKeys()->getAllParams();
    map<string, Param*>::iterator curKey;
    
    // start from the back and it has to be int and not unsigned int!
    std::size_t start=arg.length(), last, startReplace, lengthReplace;
    string word;
    char c;// = arg[start];
    while(start>0){
        // get each word separately
        last=start;
        do{
            c = arg[--start];
        }while(start>0 && !isspace(c) && c!='(' && c!=')' && c!='{'  && c!='}' && c!='[' && c!=']');
    
        if(start){  // if it is not the first word
            startReplace = start+1;
            lengthReplace = last-start-1;
        }
        else{
            startReplace = start;
            lengthReplace = last-start;
        }
        
        word = arg.substr(startReplace, lengthReplace);
        curKey = pKeys->find(word);
        if(curKey==pKeys->end() || word==_name) continue; // check if keyword, avoid circularity
        
        
        // it is a keyword: replace the keyword by its argument
        arg.replace(startReplace, lengthReplace, curKey->second->get_arg());
        
        // if the keyword is temporal: add the changes also to this param
        if(!curKey->second->get_temporal_args()->empty()){
            if(!_temporalParamAllowed){
                error("Temporal parameter: Parameter %s can not change over time!\n", _name.c_str());
            }
            
            // add the new temporal changes to it ("keeping the same argument")
            map<int, string>* key_map = curKey->second->get_temporal_args();
            map<int, string>::iterator pos, curT=key_map->begin(), endT=key_map->end();
            for(; curT!=endT; ++curT){
                if(_temporalArgs.find(curT->first)!=_temporalArgs.end()) continue;
                
                // is tha thte first temporal parameter?
                if(_temporalArgs.empty()) {
                    _temporalArgs.insert(pair<int, string>(curT->first, _tot_arg));
                }
                else {
                    pos = _temporalArgs.upper_bound(curT->first),
                    --pos;
                    _temporalArgs[curT->first] = pos->second;
                }
            }
            
            // add the generation times to the list in the paramSet
            for(map<int, string>::iterator pos = _temporalArgs.begin(); pos != _temporalArgs.end(); ++pos){
                _parSet->set_temporal_param(pos->first, this);
            }
        }
        
    }
    return arg;
}

// ----------------------------------------------------------------------------------------
// replaceKeywords
// ----------------------------------------------------------------------------------------
/** replaces all occurences of keywords by the argument of the keyword
 *
string
Param::replaceKeywords(string arg)
{
    assert(_parSet),
    assert(_parSet->get_pParamManager());
    if(!_parSet->get_pParamManager()->get_paramSetKeys()) return arg;   // there are no kewords
    
    map<string, Param*>* pKeys=&_parSet->get_pParamManager()->get_paramSetKeys()->getAllParams();
    map<string, Param*>::iterator curKey;
    std::size_t start=0, last;
    while(1){
        start = arg.find("$", start);           // find next occurence
        if(start==string::npos) break;          // if not present: stop
        
        // there seems to be a keyword: get the keyword
        last=start;                             // search end of keyword
        while(isalpha(arg[++last])){}
        curKey = pKeys->find(arg.substr(start+1, last-start-1));
        if(curKey==pKeys->end()) continue;      // the $ does not belong to a key word
        
        // it is a keyword: replace the keyword by its argument
        arg.replace(start, last-start, curKey->second->get_arg());
        
        // if the keyword is temporal: add the changes also to this param
        if(!curKey->second->get_temporal_args()->empty()){
            if(!_temporalParamAllowed){
                error("Temporal parameter: Parameter %s can not change over time!\n", _name.c_str());
            }
            
            // add the new temporal changes to it ("keeping the same argument")
            map<int, string>* key_map = curKey->second->get_temporal_args();
            map<int, string>::iterator pos, curT=key_map->begin(), endT=key_map->end();
            for(; curT!=endT; ++curT){
                if(_temporalArgs.find(curT->first)!=_temporalArgs.end()) continue;
                
                // is tha thte first temporal parameter?
                if(_temporalArgs.empty()) {
                    _temporalArgs.insert(pair<int, string>(curT->first, _tot_arg));
                }
                else {
                    pos = _temporalArgs.upper_bound(curT->first),
                    --pos;
                    _temporalArgs[curT->first] = pos->second;
                }
            }
            
            // add the generation times to the list in the paramSet
            for(map<int, string>::iterator pos = _temporalArgs.begin(); pos != _temporalArgs.end(); ++pos){
                _parSet->set_temporal_param(pos->first, this);
            }
        }
        
    }
    return arg;
}
*/
// ----------------------------------------------------------------------------------------
// replaceMacros
// ----------------------------------------------------------------------------------------
/** replaces all macros
 * if macros were present the new value is stored in _tot_arg_no_macro,
 * otherwise _tot_arg_no_macro=""
 */
string
Param::replaceMacros(string arg)
{
    if(!_pReplicate) return arg;
    
    _tot_arg_no_macro = replace_macros(arg); // _tot_arg_no_macro==NULL if arg has no macros
    if(_tot_arg_no_macro.empty()) return arg;       // if no macros are involved
    
    if(!_isRandom) return _tot_arg_no_macro;
    if(_pReplicate->get_pSimulation()->_randomPerReplicate) return _tot_arg_no_macro;
    
    unsigned int curRep = _pReplicate->_current_replicate;
    if(curRep==my_NAN) return _tot_arg_no_macro;
    
    // random numbers are identical among replicates
    // get the main paramter (_current_replciate==my_NAN)
    Param* curParam;
    list<ParamSet*>* pMainList = _pReplicate->get_pSimulation()->_allParams;
    list<ParamSet*>::iterator cur, end;
    for(cur=pMainList->begin(), end=pMainList->end(); cur!=end; ++cur){
        curParam = (*cur)->find_param(_name);
        if(curParam) break;
    }
    assert(cur!=end);
    
    // if it is not temproal just overake _arg
    if(_temporalArgs.empty() || !_pReplicate->get_popPtr()->getCurrentGeneration()){
        _tot_arg_no_macro = curParam->_tot_arg_no_macro;
        return _tot_arg_no_macro;
    }
    
    // it is a temproal parameter: search the correct value
    assert(curParam->get_temporal_args()->find(_pReplicate->get_popPtr()->getCurrentGeneration())!=curParam->get_temporal_args()->end());
    _tot_arg_no_macro = curParam->get_temporal_args()->find(_pReplicate->get_popPtr()->getCurrentGeneration())->second;
    return _tot_arg_no_macro;
}

// ----------------------------------------------------------------------------------------
// set
// ----------------------------------------------------------------------------------------
/** read a temporal parameter and store each change in _temporalArgs[gen]==argument. 
 * A temporal parameter is enclosed by () and changes are difined by a pair of time 
 * followed by the argument (separeted by space). The changes themselves are separated 
 * by , or ;
 * Comments have previously been removed.
 * Returned is the current (first) argument.
 */
string
Param::set_temporal_args (const string& arg)
{
    if(arg[0] != '(') return arg;  // it is not a temporal parameter
	
    if(!_temporalParamAllowed){
		error("Temporal parameter: Parameter %s can not change over time!\n", _name.c_str());
	}

	int gen, counter=0;
	unsigned int linecnt=0;
	string curArg;
	istringstream IN(arg);
	char c;
	_temporalArgs.clear();

	IN.get(c);      // remove "("

	// for each temporal setting
	while(IN.good() && !IN.eof() && c != ')'){
		// read the generation time
		if(!STRING::removeCommentAndSpace(IN, linecnt)) error("Reading temporal series!\n");
		IN >> gen;
		IN >> ws;

		// read the argument
		curArg = "";
		while((IN.get(c) && IN.good() && !IN.eof() && c != ';' && c != ',' && c != ')') || counter>0){
            switch (c) {
                default: // read
                    curArg += c;
                    break;
                    
                case '\\': // line continuation
                    while (IN.get(c) && IN.good() && !IN.eof() && !isEOL(c, IN)) {
                    }   // remove rest of line
                    linecnt++;
                    break;
                    
                case '{':
                    IN.putback(c);
                    curArg += STRING::readBrackets2String(IN, linecnt, '}');
                    break;
                    
                case '(':
                    IN.putback(c);
                    curArg += STRING::readBrackets2String(IN, linecnt, ')');
                    break;
                    
                case '[':
                    IN.putback(c);
                    curArg += STRING::readBrackets2String(IN, linecnt, ']');
                    break;
                    
                case '\"':
                    IN.putback(c);
                    curArg += STRING::readUntilCharacter(IN, linecnt, '\"', true,
                                                       true);
                    break;
                    
            }
		}

		if(counter) error("Temporal parameter: The parameter %s has a problem with the brackets '{}'!\n", _name.c_str());

		if(gen < 1) gen=1;                // the first generation time has to be 1!!!

		// store the argument
		if(_temporalArgs.find(gen) != _temporalArgs.end()){
			warning("Temporal parameter: Generation %i of parameter %s already set!\n", gen, _name.c_str());
		}
		if(!check_arg(get_new_arg(curArg))) return curArg;   // this forces quantiNemo to stop
		_temporalArgs[gen] = curArg;      // store the variables
	}

	// the first argument should be the one used to start the simulation (generations are sorted by the map)
	if(_temporalArgs.begin()->first != 1) error("Temporal parameter: The parameter %s is not defined for the first generation!\n", _name.c_str());

	// if the temporal parameter has only a single element it is indeed not a temporal parameter
	if(_temporalArgs.size() == 1){
		_temporalArgs.clear();
	}
	else{
		// add the generation times to the list in the paramSet
		for(map<int, string>::iterator pos = _temporalArgs.begin(); pos != _temporalArgs.end(); ++pos){
			_parSet->set_temporal_param(pos->first, this);
		}
	}

	return _temporalArgs.begin()->second;
}

// ----------------------------------------------------------------------------------------
// get_value
// ----------------------------------------------------------------------------------------
double
Param::get_value ()
{
	assert(!(is_matrix() || _type == STR));
	try{
	   	//cout << endl << _name << "            " << _arg.c_str();
		return strTo<double>(_arg.c_str());
	}catch(...){
		if(_arg=="NaN") return my_NAN;
		throw;
	}
	return my_NAN;	// function has to return something althoug not used here...
}

// ----------------------------------------------------------------------------------------
// get_matrix
// ----------------------------------------------------------------------------------------
TMatrix*
Param::get_matrix ()
{
	if( is_matrix() && _isSet){
		try{
			return new TMatrix(_arg);
		}
		catch(const char* err){
			error("Parameter '%s': %s\n", _name.c_str(), err);
		}
	}
	return NULL;
}

// ----------------------------------------------------------------------------------------
// get_as_matrix: forces to be a matrix
// ----------------------------------------------------------------------------------------
/** if allowStr is set to true, then text is stored in a double matrix
 */
TMatrix*
Param::get_as_matrix (bool allowStr)
{
	assert(_type==MAT || _type==MAT_VAR || _type==INT_MAT || _type==DBL_MAT || _type==STR_MAT); // be sure that this may be a matrix
	try{
		if(_arg[0] == '{') return new TMatrix(_arg, allowStr);
		return new TMatrix("{"+_arg+"}", allowStr);
	}
	catch(const char* err){
		error("Parameter '%s': %s\n", _name.c_str(), err);
	}
	return NULL;
}

// ----------------------------------------------------------------------------------------
// show_up
// ----------------------------------------------------------------------------------------
void
Param::show_up ()
{
    message("\n%s\t%s\t%s\t%i", _name.c_str(), get_type_str().c_str(),
            _default_arg.empty() ? "\"\"":_default_arg.c_str(), _isRequired ? 1:0);
}

// ----------------------------------------------------------------------------------------
// help
// ----------------------------------------------------------------------------------------
void
Param::print_help (ostream& os, unsigned int wide1, unsigned int wide2, char fill,
                   unsigned int importance)
{
    if(_importance != importance) return;
    
    // the parameter name
    os << "\n\n--" << left << setw(wide1-4) << setfill(fill) << _name << setfill(' ');
    
    // the parameter type
    os << "  [" << get_type_str() << "]";
    
    // argument range
    if(_bounds[0]!=my_NAN || _bounds[1]!=my_NAN){
        if(_bounds[0]==my_NAN) os << " ]-;";
        else os << " [" << _bounds[0] << ";";
        if(_bounds[1]==my_NAN) os << "-[";
        else os << _bounds[1] << "]";
    }
    
    // default - required - temporal
    if(!_default_arg.empty() || _temporalParamAllowed || _isRequired){
        os << " (";
        if(_temporalParamAllowed) os << "temporal";
        if(_temporalParamAllowed && (_isRequired || !_default_arg.empty())) os << "/";
        if(_isRequired) os << "required";
        if(_isRequired && !_default_arg.empty()) os << "/";
        if(!_default_arg.empty()) os << "default: " << _default_arg;
        os << ")";
    }
    
    os << "  {" << importance << "}";
    

    // description
    if(!_description.empty()){
        bool stop;
        std::size_t posStart=0, posSpace2=0, posReturn=0, posSpace=0, startWS=0, size=_description.length();
        do{
            posReturn = _description.find_first_of("\n", posStart);
            stop=false;
            startWS=0;
            do{
                // get the end position
                posSpace = posStart+2.5*startWS+wide2;
                if(posSpace >= posReturn){posSpace=posReturn; stop=true;}
                else if(posSpace >= size){posSpace=size; stop=true;}
                else posSpace  = _description.find_last_of(" ", posSpace);
                
                // remove trailing white space
                posSpace2=posSpace;
                while(isspace(_description[posSpace-1])) {posSpace--;}
                
                 //cout << endl << "'" << _description.substr(posStart, posSpace-posStart) << "' " << posStart << " - " << posSpace << endl;
                os<<"\n"<<setw(wide1+(unsigned int)(2.5*startWS))<<""<<_description.substr(posStart, posSpace-posStart)<<flush;
                
                // get the number of forgoing white spaces
                if(!startWS){
                    while(isspace(_description[posStart+startWS])) {++startWS;}
                }

                posStart = posSpace2 + 1;
            }while(!stop);
        }while(posReturn!=std::string::npos);
    }
}

// ----------------------------------------------------------------------------------------
// show_up
// ----------------------------------------------------------------------------------------
string
Param::get_type_str()
{
	switch(_type){
	case DBL:       return "double";
	case INT2:      return "integer";
	case STR:       return "string";
	case MAT:       return "matrix";
	case DIST:      return "distribution";
	case MAT_VAR:   return "variable matrix";
	case INT_MAT:   return "integer/matrix";
	case DBL_MAT:   return "double/matrix";
	case STR_MAT:   return "string/matrix";
	}
	return "";
}

// ----------------------------------------------------------------------------------------
// get_matrix
// ----------------------------------------------------------------------------------------
TMatrixVar<string>*
Param::get_matrixVarStr ()
{
	if( is_matrixVar() && _isSet ){
		TMatrixVar<string>* mat = new TMatrixVar<string>();
		try{
			mat->read(_arg);
			return mat;
		}
		catch(const char* err){
			if(mat) delete mat;
			error("Parameter '%s': %s\n", _name.c_str(), err);
		}
		catch(...){
			error("Parameter '%s': Could not parse the matrix!\n", _name.c_str());
		}
	}
	return NULL;
}

// ----------------------------------------------------------------------------------------
// Random variables
// ----------------------------------------------------------------------------------------
double Param::Uniform(){
    _isRandom=true;
    if(_pReplicate->rand) return _pReplicate->rand->Uniform();
    return _pReplicate->get_pSimulation()->randEngines[0].Uniform();
}
double Param::Uniform(const double& min, const double& max){
    _isRandom=true;
    return (max-min)*getParamSet()->_popPtr->rand().Uniform()+min;
}
double Param::Poisson(const double& mean){
    _isRandom=true;
    return getParamSet()->_popPtr->rand().Poisson(mean);
}
double Param::Normal(const double& mean, const double& sd){
    _isRandom=true;
    return getParamSet()->_popPtr->rand().Normal(mean, sd);
}
double Param::LogNormal(const double& mean, const double& sd){
    _isRandom=true;
    return getParamSet()->_popPtr->rand().LogNormal(mean, sd);
}
double Param::Gamma(const double& a, const double& b){
    _isRandom=true;
    return getParamSet()->_popPtr->rand().Gamma(a, b);
}
double Param::Binomial(const unsigned int& size, const double& prob){
    _isRandom=true;
    return getParamSet()->_popPtr->rand().Binomial(prob, size);
}
double Param::Beta(const double& a, const double& b){
    _isRandom=true;
    return getParamSet()->_popPtr->rand().Beta(a, b, 0, 1);
}
double Param::Beta(const double& a, const double& b, const double& min, const double& max){
    _isRandom=true;
    return getParamSet()->_popPtr->rand().Beta(a, b, min, max);
}

// ----------------------------------------------------------------------------------------
// runif
// ----------------------------------------------------------------------------------------
/** the text is reaplced by random deviates form a uniform distribution
 * runif(nb, min, max)
 * random numbers between [min, max[
 * nb:    number of random deviates
 */
void
Param::runif(string& s, double min, double max)
{
	string::size_type pos=0; // count the number of commas
	unsigned int nbComma=0;
	while(pos!=string::npos){
		pos = s.find(',', ++pos); // first position of the string should not (is not allowed) to be a comma!
		++nbComma;
	}
	--nbComma; // it was counted one too much
	vector<string> vec;
	if(nbComma==0) vec = macro_0args(s, "runif", &Param::Uniform);
	else           vec = macro_2args<double, double>(s, "runif", &Param::Uniform);
	if(min==my_NAN) s = vector2string(vec);
	else{
		unsigned int n = strTo<unsigned int>(s.substr(0, s.find(',')));   // at the start
		vector<double> time_vec = get_seq(min, max, n, true);
		s = vector2string(vec, time_vec);
	}
}

// ----------------------------------------------------------------------------------------
// rpois
// ----------------------------------------------------------------------------------------
/** the text is reaplced by random deviates form a poisson distribution
 * rpois(nb, mean)
 * mean:  mean of the distribution
 * nb:    number of random deviates
 */
void
Param::rpois(string& s, double min, double max)
{
	vector<string> vec = macro_1args<double>(s, "rpois", &Param::Poisson);
	if(min==my_NAN) s = vector2string(vec);
	else{
		unsigned int n = strTo<unsigned int>(s.substr(0, s.find(',')));   // at the start
		vector<double> time_vec = get_seq(min, max, n, true);
		s = vector2string(vec, time_vec);
	}
}

// ----------------------------------------------------------------------------------------
// rnorm
// ----------------------------------------------------------------------------------------

/** the text is reaplced by random deviates from a normal distribution
 * rnorm(nb, mean, sd)
 * nb:    number of random deviates
 * mean:  mean of the distribution (may be a matrix)
 * sd:    standard deviation of the distribution (may be a matrix)
 */
void
Param::rnorm(string& s, double min, double max)
{
	vector<string> vec = macro_2args<double, double>(s, "rnorm", &Param::Normal);
	if(min==my_NAN) s = vector2string(vec);
	else{
		unsigned int n = strTo<unsigned int>(s.substr(0, s.find(',')));   // at the start
		vector<double> time_vec = get_seq(min, max, n, true);
		s = vector2string(vec, time_vec);
	}
}

// ----------------------------------------------------------------------------------------
// rgamma
// ----------------------------------------------------------------------------------------
/** the text is reaplced by random deviates form a gamma distribution
 * rgamma(nb, shape, scale)
 * shape: shape of the gamma distribution (may be a matrix)
 * scale: r.Gamma(shape)/scale            (may be a matrix)
 * nb:    number of random deviates
 */
void
Param::rgamma(string& s, double min, double max)
{
	vector<string> vec = macro_2args<double, double>(s, "rgamma", &Param::Gamma);
	if(min==my_NAN) s = vector2string(vec);
	else{
		unsigned int n = strTo<unsigned int>(s.substr(0, s.find(',')));   // at the start
		vector<double> time_vec = get_seq(min, max, n, true);
		s = vector2string(vec, time_vec);
	}
}

// ----------------------------------------------------------------------------------------
// rbeta
// ----------------------------------------------------------------------------------------
/** the text is replaced by random deviates form a beta distribution
 * rbeta(nb, alpha, beta)
 * rbeta(nb, alpha, beta, min, max)
 * alpha: scale parameter 1
 * beta:  scale parameter 2
 * min:   minimum range
 * max:   maximal range
 * nb:    number of random deviates
 */
void
Param::rbeta(string& s, double min, double max)
{
	string::size_type pos=0; // count the number of commas
	unsigned int nbComma=0;
	while(pos!=string::npos){
		pos = s.find(',', ++pos); // first position of the string should not (is not allowed) to be a comma!
		++nbComma;
	}
	--nbComma; // it was counted one too much
	vector<string> vec;
	if(nbComma==2) vec = macro_2args<double, double>(s, "rbeta", &Param::Beta);
	else           vec = macro_4args<double, double, double, double>(s, "rbeta", &Param::Beta);

	if(min==my_NAN) s = vector2string(vec);
	else{
		unsigned int n = strTo<unsigned int>(s.substr(0, s.find(',')));   // at the start
		vector<double> time_vec = get_seq(min, max, n, true);
		s = vector2string(vec, time_vec);
	}
}

// ----------------------------------------------------------------------------------------
// rlnorm
// ----------------------------------------------------------------------------------------
/** the text is replaced by random deviates form a log normal distribution
 * rlnorm(nb, meanlog, sdlog)
 * meanlog: mean and standard deviation of the distribution on the log scale (may also be a matrix)
 * sdlog:   mean and standard deviation of the distribution on the log scale (may also be a matrix)
 * nb:      number of random deviates
 */
void
Param::rlnorm(string& s, double min, double max)
{
	vector<string> vec = macro_2args<double, double>(s, "rlnorm", &Param::LogNormal);
	if(min==my_NAN) s = vector2string(vec);
	else{
		unsigned int n = strTo<unsigned int>(s.substr(0, s.find(',')));   // at the start
		vector<double> time_vec = get_seq(min, max, n, true);
		s = vector2string(vec, time_vec);
	}
}

// ----------------------------------------------------------------------------------------
// rbinom
// ----------------------------------------------------------------------------------------
/** the text is replaced by random deviates form a binomial distribution
 * rbinom(nb, size, prob)
 * nb:    number of random deviates
 * size:  number of trials
 * prob:  probability of success on each trial
 */
void
Param::rbinom(string& s, double min, double max)
{
	vector<string> vec = macro_2args<unsigned int, double>(s, "rbinom", &Param::Binomial);
	if(min==my_NAN) s = vector2string(vec);
	else{
		unsigned int n = strTo<unsigned int>(s.substr(0, s.find(',')));   // at the start
		vector<double> time_vec = get_seq(min, max, n, true);
		s = vector2string(vec, time_vec);
	}
}

// ----------------------------------------------------------------------------------------
// rsample
// ----------------------------------------------------------------------------------------
/** the text is replaced by randomly sampling inthe given number with or without replacement
 * rsample(nb, replace, val1, val2, ...)
 * replace: 0 for no (in this case nb must be equal or smaller than the number of vals) and 1 for yes
 * valX:    a series of numbers or text separated by a comma to sample
 * nb:      number of random deviates
 */
void
Param::rsample(string& s, double min, double max)
{
	vector<string> vec = rsample_vec(s);
	if(min==my_NAN) s = vector2string(vec);
	else{
		unsigned int n = strTo<unsigned int>(s.substr(0, s.find(',')));   // at the start
		vector<double> time_vec = get_seq(min, max, n, true);
		s = vector2string(vec, time_vec);
	}
}

// ----------------------------------------------------------------------------------------
vector<string>
Param::rsample_vec(const string& t)
{
	istringstream IN;                             // make a new stream
	IN.str(t);                                    // allocate the line to the new stream
	unsigned int n;                               // number of deviates
	unsigned int replace;                         // should sampling be done with our without replacement
	char c;                                       // comma
	vector<string> vec;                           // the text to be sampled
	string valStr;

	// get the number of elements
	IN >> n >> c;
	if(c!=',') error("rsample '%s': The parameters have to be separated by a comma ','!\n", t.c_str());
	if(n<1)       error("rsample '%s': The number of deviates has to be positive!\n", t.c_str());

	// get if drawing is with replacement
	IN >> replace >> c;
	if(c!=',') error("rsample '%s': The parameters have to be separated by a comma ','!\n", t.c_str());
	if(replace!=0 && replace!=1)  error("rsample '%s': The second parameter has to be either 0 (without replacement) or 1 (with replacement)!\n", t.c_str());
	if(IN.eof()) error("Sequence 'seq(%s)': Problems reading the arguments!\n", t.c_str());

	while(!IN.eof()){
		// read the element
		valStr.clear();
		IN >> ws;                                  // remove forgoing white space
		while(IN.get(c) && c!=',' && c!=')'){      // get the text between '(' and ','
			valStr += c;
		}
		rem_trailing_space(valStr);
		if(!valStr.empty()) vec.push_back(valStr); // store the element
	}

	// create the sequence
	vector<string> newVec(getParamSet()->_popPtr->rand().sample(vec, n, (bool) replace));
	return newVec;
}

// ----------------------------------------------------------------------------------------
// seq
// ----------------------------------------------------------------------------------------
/** function seq adapted to the function seq in R:
 * format: seq(from, to, by)              // simple sequence (from and to may be a matrix, but does not really makes sense...)
 * input: the entire format including seq and the parantheses  "seq(1,10,2)"
 * output: the sequence specified  "1 3 5 7 9"
 */
void
Param::seq(string& s, double min, double max)
{
	vector<string> vec = seq_vec(s);
	if(min==my_NAN) s = vector2string(vec);
	else{
		unsigned int n = strTo<unsigned int>(s.substr(s.rfind(',')+1));   // at the end
		vector<double> time_vec = get_seq(min, max, n, true);
		s = vector2string(vec, time_vec);
	}
}

// ----------------------------------------------------------------------------------------
vector<string>
Param::seq_vec(const string& t)
{
	istringstream IN;                             // make a new stream
	IN.str(t);                                    // allocate the line to the new stream
	double from=my_NAN, to=my_NAN;
	TMatrix *mat1=NULL, *mat2=NULL;
	unsigned int n;                               // number of steps
	char c;                                       // comma 1, comma 2
	unsigned int linecnt = 0;                              // not used, but must be present

	// read from
	IN >> ws;
	if(IN.peek()=='{') mat1 = new TMatrix(STRING::readBrackets2String(IN, linecnt, '}'));
	else IN >> from;
	IN >> ws >> c;
	if(c!=',') error("Sequence 'seq(%s)': The parameters have to be separated by a comma ','!\n", t.c_str());
	if(IN.eof()) error("Sequence 'seq(%s)': Problems reading the arguments!\n", t.c_str());

	// read to
	IN >> ws;
	if(IN.peek()=='{') mat2 = new TMatrix(STRING::readBrackets2String(IN, linecnt, '}'));
	else IN >> to;
	IN >> ws >> c;
	if(c!=',')    error("Sequence 'seq(%s)': The parameters have to be separated by a comma ','!\n", t.c_str());
	if(IN.eof())  error("Sequence 'seq(%s)': Problems reading the arguments!\n", t.c_str());

	//read nb step
	IN >> n >> ws;
	if(n<=1)      error("Sequence 'seq(%s)': The step number has to be 2 or higher (present: %i)!\n", t.c_str(), n);
	if(!IN.eof()) error("Sequence 'seq(%s)': Problems reading the arguments!\n", t.c_str());

	// create the sequence
	vector<string> vec;
	vec.reserve(n);
	if(mat1 || mat2){ // at least one parameter is specified by a matrix		if(!from_m) from_m = new TMatrix("{"+toStr(from)+"}");  // get it as a matrix
		if(!mat1)   mat1   = new TMatrix("{"+toStr(from)+"}");    // get it as a matrix
		if(!mat2)   mat2   = new TMatrix("{"+toStr(to)+"}");      // get it as a matrix

		vector<TMatrix> vecMat = get_seq(*mat1, *mat2, n); 	// write the matrix
		for(unsigned int i=0; i<n; ++i){
			vec.push_back(vecMat[i].toStr());
		}
		delete mat1;
		delete mat2;
	}
	else{	// simple sequence
		vector<double> vecVal = get_seq(from, to, n);
		for(unsigned int i=0; i<n; ++i){
			vec.push_back(toStr(vecVal[i]));
		}
	}
	return vec;
}

// ----------------------------------------------------------------------------------------
// seq
// ----------------------------------------------------------------------------------------
/** returns a vector of sequences between from and to and the length by.
 * if asInt is true, then the values are rounded for integers
 */
vector<double>
Param::get_seq(double from, const double& to, const unsigned int& by, bool asInt)
{
	vector<double> vec;
	vec.reserve(by);
	double step = (to-from)/(by-1);
	for(unsigned int i=0; i<by; ++i, from+=step){
		if(asInt) vec.push_back(my_round(from));
		else      vec.push_back(from);
	}
	return vec;
}

// ----------------------------------------------------------------------------------------
// seq
// ----------------------------------------------------------------------------------------
/** returns a vector of sequences between from and to and the length by.
 * if asInt is true, then the values are roudned fro integers
 */
vector<TMatrix>
Param::get_seq(const TMatrix& fromMatrix, const TMatrix& toMatrix, const unsigned int& by, bool asInt)
{
	unsigned int nbCol1 = fromMatrix.getNbCols();
	unsigned int nbCol2 = toMatrix.getNbCols();
	unsigned int nbRow1 = fromMatrix.getNbRows();
	unsigned int nbRow2 = toMatrix.getNbRows();
	unsigned int nbCol = nbCol1>nbCol2 ? nbCol1 : nbCol2;
	unsigned int nbRow = nbRow1>nbRow2 ? nbRow1 : nbRow2;

	// get step size matrix and the starting matrix
	TMatrix steps(nbRow, nbCol);     // contains the step size change
	TMatrix start(nbRow, nbCol);     // contains the starting values
	unsigned int x, y;
	double from, to;
	for(x=0; x<nbRow; ++x){
		for(y=0; y<nbCol; ++y){
			from = fromMatrix.get(x % nbRow1, y % nbCol1);
			to   = toMatrix.get(x % nbRow2, y % nbCol2);
			start.set(x,y,from);
			steps.set(x,y,(to-from)/(by-1));
		}
	}

	// write the vector
	vector<TMatrix> vec;
	vec.reserve(by);
	for(unsigned int i=0; i<by; ++i){
		vec.push_back(TMatrix(nbRow, nbCol));
		for(x=0; x<nbRow; ++x){
			for(y=0; y<nbCol; ++y){
				if(asInt) vec.back().set(x, y, my_round(start.get(x,y) + i*steps.get(x,y)));
				else      vec.back().set(x, y, start.get(x,y) + i*steps.get(x,y));
			}
		}
	}

	return vec;
}

// ----------------------------------------------------------------------------------------
// logistic
// ----------------------------------------------------------------------------------------
/** returns a vector of sequences between from and to following a logistic slope and the length by.
 * if asInt is true, then the values are rounded for integers

 ** that is an implementation of the general logistic curve
 * (Richards, F.J. 1959 A flexible growth function for empirical use. J. Exp. Bot. 10: 290--300.)
 * The curve is defined by 5 parameters.
 * Do to size problems a slope of more than (+/-)1e4 is regarded as instantenious change from min to max

 * double generalLogisticCurve(const double& x,     // time
                                   const double& min,   // the lower asymptote
                                   const double& max,   // the upper asymptote
                                   const double& max_r, // the time of maximum growth
                                   const double& r,     // the growth rate
																	 const double& s)     // affects near which asymptote maximum growth occurs
 */
vector<double>
Param::get_logistic(const double& from, const double& to, const double& r, const double r_max,
		const double& s, const unsigned int& by, bool asInt)
{
	vector<double> vec;
	vec.reserve(by);
	double step = 1.0/(by-1);        // range is 1 (set arbitrary...)
	double val, cur = 0;
	for(unsigned int i=0; i<by; ++i, cur+=step){
		val = generalLogisticCurve(cur,from,to,r_max,r,s);
		if(asInt) vec.push_back(my_round(val));
		else      vec.push_back(val);
	}
	return vec;
}

// ----------------------------------------------------------------------------------------
// logistic
// ----------------------------------------------------------------------------------------
/** returns a vector of sequences between from and to following a logistic slope and the length by.
 * if asInt is true, then the values are roudned for integers
 */
vector<TMatrix>
Param::get_logistic(const TMatrix& fromMatrix, const TMatrix& toMatrix,
		const TMatrix& rMatrix,    const TMatrix& r_maxMatrix,
		const TMatrix& sMatrix, const unsigned int& by, bool asInt)
{
	unsigned int nbCol1 = fromMatrix.getNbCols();
	unsigned int nbCol2 = toMatrix.getNbCols();
	unsigned int nbCol3 = rMatrix.getNbCols();
	unsigned int nbCol4 = r_maxMatrix.getNbCols();
	unsigned int nbCol5 = sMatrix.getNbCols();
	unsigned int nbRow1 = fromMatrix.getNbRows();
	unsigned int nbRow2 = toMatrix.getNbRows();
	unsigned int nbRow3 = rMatrix.getNbRows();
	unsigned int nbRow4 = r_maxMatrix.getNbRows();
	unsigned int nbRow5 = sMatrix.getNbRows();
	unsigned int nbCol = nbCol1>nbCol2 ? nbCol1 : nbCol2;    // find max number of columns
	if(nbCol<nbCol3) nbCol = nbCol3;
	if(nbCol<nbCol4) nbCol = nbCol4;
	if(nbCol<nbCol5) nbCol = nbCol5;
	unsigned int nbRow = nbRow1>nbRow2 ? nbRow1 : nbRow2;    // find max number of rows
	if(nbRow<nbRow3) nbRow = nbRow3;
	if(nbRow<nbRow4) nbRow = nbRow4;
	if(nbRow<nbRow5) nbRow = nbRow5;

	// get step size matrix and the starting matrix
	TMatrix steps(nbRow, nbCol);     // contains the step size change
	unsigned int x, y;
	double val;
	for(x=0; x<nbRow; ++x){
		for(y=0; y<nbCol; ++y){
			val = 1.0/(by-1);                  // range set to 1 (arbitrary)
			steps.set(x,y,val);
		}
	}

	// write the vector
	vector<TMatrix> vec;
	vec.reserve(by);
	for(unsigned int i=0; i<by; ++i){
		vec.push_back(TMatrix(nbRow, nbCol));
		for(x=0; x<nbRow; ++x){
			for(y=0; y<nbCol; ++y){
				val = generalLogisticCurve(i*steps.get(x,y),
						fromMatrix.get(x % nbRow1, y % nbCol1),
						toMatrix.get(x % nbRow2, y % nbCol2),
						rMatrix.get(x % nbRow3, y % nbCol3),
						r_maxMatrix.get(x % nbRow4, y % nbCol4),
						sMatrix.get(x % nbRow5, y % nbCol5));
				if(asInt) vec.back().set(x, y, my_round(val));
				else      vec.back().set(x, y, val);
			}
		}
	}

	return vec;
}

// ----------------------------------------------------------------------------------------
// seq
// ----------------------------------------------------------------------------------------
/** function seq2D allows to generate a sequence for a two dimensional world (and works just in this case)
 * format: seq2D(xlim, ylim, patch1, patch2, value1, value2, by)
 * input: the entire format including seq2D and the parantheses  "seq2D(1,10,2)"
 * output: the sequence specified  "1 3 5 7 9"
 * alpha is the angle of the line connecting the two points, starting at horizontal
 */
void
Param::seq2D(string& s, double min, double max)
{
	vector<string> vec = seq2D_vec(s);
	if(min==my_NAN) s = vector2string(vec);
	else{
		unsigned int n = strTo<unsigned int>(s.substr(s.rfind(',')+1));   // at the end
		vector<double> time_vec = get_seq(min, max, n, true);
		s = vector2string(vec, time_vec);
	}
}

// ----------------------------------------------------------------------------------------
vector<string>
Param::seq2D_vec(const string& t)
{
	istringstream IN;                             // make a new stream
	int xlim, ylim, patchID1, patchID2;
	IN.str(t);             // allocate the line to the new stream
	double val1, val2;
	int by;                                       // number of steps
	char c1;                                       // comma 1, comma 2

	IN >> xlim >> ws >> c1;
	if(c1!=',') error("Sequence 'seq2D(%s)': The parameters have to be separated by a comma ','!\n", t.c_str());
	IN >> ylim >> ws >> c1;
	if(c1!=',') error("Sequence 'seq2D(%s)': The parameters have to be separated by a comma ','!\n", t.c_str());
	IN >> patchID1 >> ws >> c1;
	if(c1!=',') error("Sequence 'seq2D(%s)': The parameters have to be separated by a comma ','!\n", t.c_str());
	IN >> patchID2 >> ws >> c1;
	if(patchID1==patchID2) error("Sequence 'seq2D(%s)': The two points may not be identical!\n", t.c_str());
	if(c1!=',') error("Sequence 'seq2D(%s)': The parameters have to be separated by a comma ','!\n", t.c_str());
	IN >> val1 >> ws >> c1;
	if(c1!=',') error("Sequence 'seq2D(%s)': The parameters have to be separated by a comma ','!\n", t.c_str());
	IN >> val2 >> ws >> c1;
	if(c1!=',') error("Sequence 'seq2D(%s)': The parameters have to be separated by a comma ','!\n", t.c_str());
	IN >> by >> ws;
	if(by<=1) error("Sequence 'seq2D(%s)': The step number has to be 2 or higher (present: %i)!\n", t.c_str(), by);
	if(!IN.eof()) error("Sequence 'seq2D(%s)': Problems reading the arguments!\n", t.c_str());

	// get coordiantes of patch 1
	patchID1 -= 1;    // internal counter starts at 0;
	int patch1y = patchID1/xlim;
	int patch1x = patchID1 - patch1y*xlim;
	if(patch1y>=ylim) error("Sequence 'seq2D(%s)': Patch 1 is out of range!\n", t.c_str());

	// get coordiantes of patch 2
	patchID2 -= 1;    // internal counter starts at 0;
	int patch2y = patchID2/xlim;
	int patch2x = patchID2 - patch2y*xlim;
	if(patch2y>=ylim) error("Sequence 'seq2D(%s)': Patch 2 is out of range!\n", t.c_str());

	// check that patch2 is right of patch1
	if(patch1x > patch2x){        // if patch2 is left of patch1 swap them
		swap(patch1x, patch2x);
		swap(patch1y, patch2y);
		swap(val1, val2);
	}

	// get the angle
	double alpha, beta, stepSize, stepChange;
	double a, b, c;
	int dX, dY, x, y;

	a = patch2x - patch1x;
	assert(a>=0);	// it was swaped if not met
	b = patch2y - patch1y;
	c = sqrt((double)a*a + b*b);
	stepSize = c/(by-1);
	stepChange = (val2-val1)/(by-1);
	if(a) alpha = atan(((double) b)/a);      // slope of the line
	else if(b>0) alpha = PI/2.0;
	else         alpha = -PI/2.0;      // if b<0

	ostringstream out;
	for(y=0; y<ylim; ++y){                   // for each element of the 2D matrix
		for(x=0; x<xlim; ++x){
			dX = x-patch1x;
			dY = y-patch1y;
			if(dX) beta = atan(((double)dY)/dX);
			else if(dY>0) beta = PI/2.0;
			else          beta = -PI/2.0;

			if(x || y) out << ' ';
			if(dX<0) out << (val1 -  std::ceil(sqrt((double)dX*dX + dY*dY) * cos(beta-alpha)/stepSize)*stepChange);
			else     out << (val1 + std::floor(sqrt((double)dX*dX + dY*dY) * cos(beta-alpha)/stepSize)*stepChange);
		}
	}
	vector<string> vec;
	vec.push_back(out.str());
	return vec;
}

// ----------------------------------------------------------------------------------------
// seq
// ----------------------------------------------------------------------------------------
/** function seq2D allows to generate a sequence for a two dimensional world (and works just in this case)
 * format: seq2D(xlim, ylim, patch1, patch2, vaue1, value2, by)
 * input: the entire format including seq2D and the parantheses  "seq2D(1,10,2)"
 * output: the sequence specified  "1 3 5 7 9"
 * alpha is the angle of the line connecting the two points, starting at horizontal
 */
void
Param::seq2Db(string& s, double min, double max)
{
	vector<string> vec = seq2Db_vec(s);
	if(min==my_NAN) s = vector2string(vec);
	else{
		unsigned int n = strTo<unsigned int>(s.substr(s.rfind(',')+1));   // at the end
		vector<double> time_vec = get_seq(min, max, n, true);
		s = vector2string(vec, time_vec);
	}
}

// ----------------------------------------------------------------------------------------
vector<string>
Param::seq2Db_vec(const string& t)
{
	istringstream IN;                             // make a new stream
	int xlim, ylim, patchID1, patchID2;
	IN.str(t);             // allocate the line to the new stream
	double val1, val2;
	int by;                                       // number of steps
	char c1;                                       // comma 1, comma 2

	IN >> xlim >> ws >> c1;
	if(c1!=',') error("Sequence 'seq2Db(%s)': The parameters have to be separated by a comma ','!\n", t.c_str());
	IN >> ylim >> ws >> c1;
	if(c1!=',') error("Sequence 'seq2Db(%s)': The parameters have to be separated by a comma ','!\n", t.c_str());
	IN >> patchID1 >> ws >> c1;
	if(c1!=',') error("Sequence 'seq2Db(%s)': The parameters have to be separated by a comma ','!\n", t.c_str());
	IN >> patchID2 >> ws >> c1;
	if(patchID1==patchID2) error("Sequence 'seq2Db(%s)': The two points may not be identical!\n", t.c_str());
	if(c1!=',') error("Sequence 'seq2Db(%s)': The parameters have to be separated by a comma ','!\n", t.c_str());
	IN >> val1 >> ws >> c1;
	if(c1!=',') error("Sequence 'seq2Db(%s)': The parameters have to be separated by a comma ','!\n", t.c_str());
	IN >> val2 >> ws >> c1;
	if(c1!=',') error("Sequence 'seq2Db(%s)': The parameters have to be separated by a comma ','!\n", t.c_str());
	IN >> by >> ws;
	if(by<=1) error("Sequence 'seq2Db(%s)': The step number has to be 2 or higher (present: %i)!\n", t.c_str(), by);
	if(!IN.eof()) error("Sequence 'seq2Db(%s)': Problems reading the arguments!\n", t.c_str());

	// get coordiantes of patch 1
	patchID1 -= 1;    // internal counter starts at 0;
	int patch1y = patchID1/xlim;
	int patch1x = patchID1 - patch1y*xlim;
	if(patch1y>=ylim) error("Sequence 'seq2D(%s)': Patch 1 is out of range!\n", t.c_str());

	// get coordiantes of patch 2
	patchID2 -= 1;    // internal counter starts at 0;
	int patch2y = patchID2/xlim;
	int patch2x = patchID2 - patch2y*xlim;
	if(patch2y>=ylim) error("Sequence 'seq2D(%s)': Patch 2 is out of range!\n", t.c_str());

	// get the angle
	string text;
	double stepSize, stepChange;
	double dist1, dist2;
	int dX, dY, x, y;

	stepSize = 1.0/(by-1);
	stepChange = (val2-val1)/(by-1);

	ostringstream out;
	for(y=0; y<ylim; ++y){                   // for each element of the 2D matrix
		for(x=0; x<xlim; ++x){
			// distance to patch1
			dX = x-patch1x;
			dY = y-patch1y;
			dist1 = sqrt((double)dX*dX + dY*dY);

			// distance to patch2
			dX = x-patch2x;
			dY = y-patch2y;
			dist2 = sqrt((double)dX*dX + dY*dY);

			if(x || y) out << ' ';
			out << (val1 + std::ceil(dist1/(dist1+dist2)/stepSize)*stepChange);
		}
	}
	vector<string> vec;
	vec.push_back(out.str());
	return vec;
}

// ----------------------------------------------------------------------------------------
// rep
// ----------------------------------------------------------------------------------------
/** function rep adapted to the function rep in R:
 * format: rep(val, nb)
 * input: the entire format including seq and the parantheses  "rep(1,10)"
 * output: the sequence specified  "1 1 1 1 1 1 1 1 1 1"
 * note that everything between the opening bracket '(' and the ',' is repeated
 */
void
Param::rep(string& s, double min, double max)
{
	vector<string> vec = rep_vec(s);
	if(min==my_NAN) s = vector2string(vec);
	else{
		unsigned int n = strTo<unsigned int>(s.substr(s.rfind(',')+1));   // at the end
		vector<double> time_vec = get_seq(min, max, n, true);
		s = vector2string(vec, time_vec);
	}
}

// ----------------------------------------------------------------------------------------
vector<string>
Param::rep_vec(const string& t)
{
	istringstream IN;                             // make a new stream
	IN.str(t);                                    // allocate the line to the new stream
	string valStr;
	int n;                                        // number of steps
	char c;                                       // comma 1, comma 2

	IN >> ws;
	while(IN.get(c) && c!=','){      // get the text between '(' and ','
		valStr += c;
	}
	rem_trailing_space(valStr);

	if(c!=',') error("rep '%s': The parameters have to be separated by a comma ','!\n", t.c_str());

	IN >> n >> ws;                                     // read the number
	if(!IN.eof()) error("rep '%s': Problems reading the arguments!\n", t.c_str());

	// create the sequence
	vector<string> vec;
	vec.reserve(n);
	for(; n; --n){
		vec.push_back(valStr);
	}
	return vec;
}

// ----------------------------------------------------------------------------------------
// logsitic
// ----------------------------------------------------------------------------------------
/** general logistic function
 * (Richards, F.J. 1959 A flexible growth function for empirical use. J. Exp. Bot. 10: 290--300.)
 * Do to size problems a slope of more than (+/-)1e4 is regarded as instantenious change from min to max
 * format: logistic(from, to, growth, max_growth, symmetry, by)
 * input: the entire format including logistic and the parantheses  "logistic(0,1,1,0.5,1,5)"
 * output: the sequence specified  "1 3 5 7 9"
 * from:   minimal value (if growth is positive)
 * to:     maximal value (if growth is positive)
 * growth: growth rate of the slope (if negative "from" and "to" are inversed)
 * max_growth: value with teh maximum growth rate
 * symmetry: a value of 1 implies a symmetric slope
inline double generalLogisticCurve(const double& x,     // time
																	 const double& min,   // the lower asymptote
																	 const double& max,   // the upper asymptote
																	 const double& max_r, // the time of maximum growth
																	 const double& r,     // the growth rate (default=0.5)
																	 const double& s)     // affects near which asymptote maximum growth occurs (default=1)


 */
void
Param::logistic(string& s, double min, double max)
{
	vector<string> vec = logistic_vec(s);
	if(min==my_NAN)s = vector2string(vec);
	else{
		unsigned int n = strTo<unsigned int>(s.substr(s.rfind(',')+1));   // at the end
		vector<double> time_vec = get_seq(min, max, n, true);
		s = vector2string(vec, time_vec);
	}
}

// ----------------------------------------------------------------------------------------
vector<string>
Param::logistic_vec(const string& t)
{
	istringstream IN;                             // make a new stream
	IN.str(t);                                    // allocate the line to the new stream
	double from=my_NAN, to=my_NAN, r=my_NAN, r_max=my_NAN, s=my_NAN;
	TMatrix *mat1=NULL, *mat2=NULL, *mat3=NULL, *mat4=NULL, *mat5=NULL;
	unsigned int n;                               // number of steps
	char c;                                       // comma
	unsigned int linecnt = 0;                              // not used, but must be present

	// read from
	IN >> ws;
	if(IN.peek()=='{') mat1 = new TMatrix(STRING::readBrackets2String(IN, linecnt, '}'));
	else IN >> from;
	IN >> ws >> c;
	if(c!=',') error("Sequence 'logistic(%s)': The parameters have to be separated by a comma ','!\n", t.c_str());
	if(IN.eof()) error("Sequence 'logistic(%s)': Problems reading the arguments!\n", t.c_str());

	// read to
	IN >> ws;
	if(IN.peek()=='{') mat2 = new TMatrix(STRING::readBrackets2String(IN, linecnt, '}'));
	else IN >> to;
	IN >> ws >> c;
	if(c!=',')    error("Sequence 'logistic(%s)': The parameters have to be separated by a comma ','!\n", t.c_str());
	if(IN.eof())  error("Sequence 'logistic(%s)': Problems reading the arguments!\n", t.c_str());

	// read growth rate
	IN >> ws;
	if(IN.peek()=='{') mat3 = new TMatrix(STRING::readBrackets2String(IN, linecnt, '}'));
	else IN >> r;
	IN >> ws >> c;
	if(c!=',')    error("Sequence 'logistic(%s)': The parameters have to be separated by a comma ','!\n", t.c_str());
	if(IN.eof())  error("Sequence 'logistic(%s)': Problems reading the arguments!\n", t.c_str());

	// read max growth
	IN >> ws;
	if(IN.peek()=='{') mat4 = new TMatrix(STRING::readBrackets2String(IN, linecnt, '}'));
	else IN >> r_max;
	IN >> ws >> c;
	if(c!=',')    error("Sequence 'logistic(%s)': The parameters have to be separated by a comma ','!\n", t.c_str());
	if(IN.eof())  error("Sequence 'logistic(%s)': Problems reading the arguments!\n", t.c_str());

	// read symmetry
	IN >> ws;
	if(IN.peek()=='{') mat5 = new TMatrix(STRING::readBrackets2String(IN, linecnt, '}'));
	else IN >> s;
	IN >> ws >> c;
	if(c!=',')    error("Sequence 'logistic(%s)': The parameters have to be separated by a comma ','!\n", t.c_str());
	if(IN.eof())  error("Sequence 'logistic(%s)': Problems reading the arguments!\n", t.c_str());

	//read nb step
	IN >> n >> ws;
	if(n<=1)      error("Sequence 'logistic(%s)': The step number has to be 2 or higher (present: %i)!\n", t.c_str(), n);
	if(!IN.eof()) error("Sequence 'logistic(%s)': Problems reading the arguments!\n", t.c_str());

	// create the sequence
	vector<string> vec;
	vec.reserve(n);
	if(mat1 || mat2 || mat3 || mat4 || mat5){ // at least one parameter is specified by a matrix		if(!from_m) from_m = new TMatrix("{"+toStr(from)+"}");  // get it as a matrix
		if(!mat1)   mat1   = new TMatrix("{"+toStr(from)+"}");    // get it as a matrix
		if(!mat2)   mat2   = new TMatrix("{"+toStr(to)+"}");      // get it as a matrix
		if(!mat3)   mat3   = new TMatrix("{"+toStr(r)+"}");       // get it as a matrix
		if(!mat4)   mat4   = new TMatrix("{"+toStr(r_max)+"}");   // get it as a matrix
		if(!mat5)   mat5   = new TMatrix("{"+toStr(s)+"}");       // get it as a matrix

		vector<TMatrix> vecMat = get_logistic(*mat1, *mat2, *mat3, *mat4, *mat5, n); 	// write the matrix
		for(unsigned int i=0; i<n; ++i){
			vec.push_back(vecMat[i].toStr());
		}
		delete mat1;
		delete mat2;
	}
	else{	// simple sequence
		vector<double> vecVal = get_logistic(from, to, r, r_max, s, n);
		for(unsigned int i=0; i<n; ++i){
			vec.push_back(toStr(vecVal[i]));
		}
	}
	return vec;
}

// ----------------------------------------------------------------------------------------
// product
// ----------------------------------------------------------------------------------------
/** function product computes the product of the two given matrixes or values:
 * format: product(val1, val2), product(val1, mat2)or product(mat1, mat2)
 * output: matrix or single value
 * to be able to reuse the generic function a 1 is added in front of the string
 */
void
Param::product(string& s, double min, double max)
{
	if(min!=my_NAN) error("Macro 'product' cannot be temporal!\n");
	vector<string> vec = macro_2args<double, double>(string("1, ") + s, "product", &Param::Product);
	s = vector2string(vec);
}

// ----------------------------------------------------------------------------------------
// round
// ----------------------------------------------------------------------------------------
/** this function rounds values in the text to entire numbers;
 */
void
Param::round(string& s, double min, double max)
{
	if(min!=my_NAN) error("Macro 'round' cannot be temporal!\n");
	istringstream is(s);
	ostringstream out;
	double val;
	char c;
	while(is.get(c)){
		if(isNumber(c)){               // it is a number
			is.putback(c);
			is >> val;
			out << my_round(val);
		}
		else out << c;                 // it is a space or a text
	}
	s = out.str();
}

// ----------------------------------------------------------------------------------------
// floor
// ----------------------------------------------------------------------------------------
/** this function floor values in the text to entire numbers;
 */
void
Param::floor(string& s, double min, double max)
{
	if(min!=my_NAN) error("Macro 'floor' cannot be temporal!\n");
	istringstream is(s);
	ostringstream out;
	double val;
	char c;
	while(is.get(c)){
		if(isNumber(c) || c=='-'){     // it is a number
			is.putback(c);
			is >> val;
			out << std::floor(val);
		}
		else out << c;                 // it is a space or a text
	}
	s = out.str();
}

// ----------------------------------------------------------------------------------------
// ceil
// ----------------------------------------------------------------------------------------
/** this function ceils values in the text to entire numbers;
 */
void
Param::ceil(string& s, double min, double max)
{
	if(min!=my_NAN) error("Macro 'ceil' cannot be temporal!\n");
	istringstream is(s);
	ostringstream out;
	double val;
	char c;
	while(is.get(c) || c=='-'){
		if(isNumber(c)){               // it is a number
			is.putback(c);
			is >> val;
			out << std::ceil(val);
		}
		else out << c;                 // it is a space or a text
	}
	s = out.str();
}

// ----------------------------------------------------------------------------------------
// trunc
// ----------------------------------------------------------------------------------------
/** all negative values are set to zero;
 */
void
Param::trunc(string& s, double min, double max)
{
	if(min!=my_NAN) error("Macro 'trunc' cannot be temporal!\n");
	istringstream is(s);
	ostringstream out;
	double val;
	char c;
	while(is.get(c)){
		if(isNumber(c)|| c=='-'){               // it is a number
			is.putback(c);
			is >> val;
			if(val<0) out << 0.0;
			else      out << val;
		}
		else out << c;                 // it is a space or a text
	}
	s = out.str();
}

// ----------------------------------------------------------------------------------------
// replace_macros
// ----------------------------------------------------------------------------------------
/** if there were any macros the text with replaced macros is returned. i
 * if no macros are present en empty string is returned
 * if at the end of the key word a 'I' is placed, then the values are rounded to integers
 */
string
Param::replace_macros(string& text)
{
	char c;
	string arg, t;
	string::iterator curPos, start;
    string::size_type pos;
	int b;
	bool hasMacro = false, time;						// flag to set if a macro was present
	unsigned int i, j, m, totSize, size;
	for(i=0; i<text.size(); ++i){
		c = text[i];
		if(c != '('){                   // get the string in front of the '('
			if(isChar(c) || isNumber(c)) arg += c;       // add char to arg
			else          arg.clear();    // reset argument
			continue;
		}

		// check if it is a temporal macro
		size = totSize = (unsigned int)arg.length();
		if(!size) continue;
		if(size>4 && arg.substr(size-4)=="Time"){
			time = true;
			arg.erase(size-4);
			size -= 4;
		}
		else time = false;

		// now arg is the macro name: check if it corresponds to a macro
		for(m=0; m<NB_MACROS; ++m){            // check for each macro
			if(size < _macros[m].size()) continue;
			if(arg.substr(size-_macros[m].size()) == _macros[m]) break;
		}
		if(m==NB_MACROS) continue;           // it was not a macro: stop here

		// it is a macro: get the text between the brackets (it may contain another macro...)
		arg.clear();                 // reset arg
		j=i-totSize;                 // get the starting postion for the replacment
		++i;
		for(b=1; i<text.size(); ++i){
			c = text[i];
			if(c==')') --b;      // closing bracket
			else if(c=='(') ++b; // opening bracket
			if(!b) break;        // stop here we have found the closing bracket
			arg += c;            // add the char to the text
		}

		replace_macros(arg);  // check recursively if there is another macro WITHIN the brackets

		try{
			if(time){ // if temporal macro remove the time from and to at the end of this sequence
				pos = arg.rfind(',');                // get time to
				unsigned int time_to = strTo<unsigned int>(arg.substr(pos+1));
				arg.erase(pos);
				pos = arg.rfind(',');                // get time from
				unsigned int time_from = strTo<unsigned int>(arg.substr(pos+1));
				arg.erase(pos);
				(this->*_macros_func_ptr[m])(arg, time_from, time_to);    // replace this macro
			}
			else (this->*_macros_func_ptr[m])(arg, my_NAN, my_NAN);   // not temporal macro
		}catch(const char* err){
			error("Macro '%s':%s!\n", text.c_str(), err);
			throw 1;
		}catch(...){
			error("Macro '%s' could not be read!\n", text.c_str());
		}

		// finilization
		text.replace(j, i-j+1, arg);
		i = j+(unsigned int)arg.size();
		arg.clear();             // reset the arguemnt, maybe there is a further macro
		hasMacro = true;
	}

	if(hasMacro) return text;
	return "";
}

// ----------------------------------------------------------------------------------------
// vector2string
// ----------------------------------------------------------------------------------------
/** this function rounds values in the text to entire numbers;
 */
string
Param::vector2string(const vector<string>& vec)
{
	ostringstream out;
	string::size_type pos, n = vec.size();
	for(pos=0; pos<n; ++pos){
		if(pos) out << ' ';
		out << vec[pos];
	}
	return out.str();
}

// ----------------------------------------------------------------------------------------
// vector2string
// ----------------------------------------------------------------------------------------
/** this function rounds values in the text to entire numbers;
 */
string
Param::vector2string(const vector<string>& vec, const vector<double>& time_vec)
{
	ostringstream out;
	assert(vec.size()==time_vec.size());
	string::size_type pos, n = vec.size();
	for(pos=0; pos<n; ++pos){
		if(pos) out << ", ";
		out << time_vec[pos] << ' ' << vec[pos];;
	}
	return out.str();
}

// ----------------------------------------------------------------------------------------
// macro_0args
// ----------------------------------------------------------------------------------------
/** the text is reaplced by random deviates from a normal distribution for a macro with no arguments: Example:
 * runif(nb)
 * nb:    number of random deviates
 */
vector<string>
Param::macro_0args(const string& t, const string& name, double (Param::*pt2Func)())
{
	istringstream IN;                             // make a new stream
	IN.str(t);                                    // allocate the line to the new stream
	unsigned int n;                               // number of deviates

	// read the number of deviates
	IN >> n >> ws;
	if(n<1)    error("%s '%s': The number of deviates has to be positive!\n", name.c_str(), t.c_str());
	if(!IN.eof()) error("%s '%s': Problems reading the arguments!\n", name.c_str(), t.c_str());

	// create the sequence
	vector<string> vec;
	vec.reserve(n);
	for(; n; --n){
		vec.push_back(toStr((this->*pt2Func)()));
	}
	return vec;
}

// ----------------------------------------------------------------------------------------
// macro_1args
// ----------------------------------------------------------------------------------------
/** the text is reaplced by random deviates from a normal distribution for a macro with one arguments: Example:
 * rpois(nb, mean)
 * nb:    number of random deviates
 * mean:  mean of the distribution (may be a matrix)
 */
template<typename F>
vector<string>
Param::macro_1args(const string& t, const string& name, double (Param::*pt2Func)(const F&))
{
	istringstream IN;                             // make a new stream
	IN.str(t);                                    // allocate the line to the new stream
	F val1=my_NAN;                                // mean and var if they are numbers
	TMatrix* mat1=NULL;                           // mean and var if they are matrixes
	unsigned int n;                               // number of deviates
	unsigned int linecnt=0;                       // temporal variable which is not really used but demanded
	char c;                                       // comma

	// read the number of deviates
	IN >> n >> c >> ws;
	if(c!=',') error("%s '%s': The parameters have to be separated by a comma ','!\n", name.c_str(), t.c_str());
	if(n<1)    error("%s '%s': The number of deviates has to be positive!\n", name.c_str(), t.c_str());
	if(IN.eof()) error("%s '%s': Problems reading the arguments!\n", name.c_str(), t.c_str());

	// read first element
	if(IN.peek()=='{') mat1 = new TMatrix(STRING::readBrackets2String(IN, linecnt, '}'));
	else IN >> val1;
	IN >> c >> ws;
	if(c!=',') error("%s '%s': The parameters have to be separated by a comma ','!\n", t.c_str());
	if(!IN.eof()) error("%s '%s': Problems reading the arguments!\n", name.c_str(), t.c_str());

	// create the sequence
	vector<string> vec;
	vec.reserve(n);
	if(mat1){ // at least one parameter is specified by a matrix
		unsigned int nbCol = mat1->getNbCols();
		unsigned int nbRow = mat1->getNbRows();

		// get step size matrix and the starting matrix
		TMatrix matrix(nbRow, nbCol);     // contains the step size change
		unsigned int y, x;
		for(; n; --n){
			ostringstream out;
			out << '{';
			for(x=0; x<nbRow; ++x){
				out << '{';
				for(y=0; y<nbCol; ++y){
					if(y) out << ' ';
					out << (this->*pt2Func)((F)mat1->get(x, y));
				}
				out << '}';
			}
			out << '}';
			vec.push_back(out.str());
		}
		delete mat1;
	}
	else{               // both parameters are dfiend by a simple number
		for(; n; --n){
			vec.push_back(toStr((this->*pt2Func)(val1)));
		}
	}
	return vec;
}

// ----------------------------------------------------------------------------------------
// macro_2args
// ----------------------------------------------------------------------------------------
/** the text is reaplced by random deviates from a normal distribution for a macro with two arguments: Example:
 * rnorm(nb, mean, sd)
 * nb:    number of random deviates
 * mean:  mean of the distribution (may be a matrix)
 * sd:    standard deviation of the distribution (may be a matrix)
 */
template<typename F, typename S>
vector<string>
Param::macro_2args(const string& t, const string& name, double (Param::*pt2Func)(const F&, const S&))
{
	istringstream IN;                             // make a new stream
	IN.str(t);                                    // allocate the line to the new stream
	F val1=my_NAN;
	S val2=my_NAN;                                // mean and var if they are numbers
	TMatrix* mat1=NULL, *mat2=NULL;               // mean and var if they are matrixes
	unsigned int n;                               // number of deviates
	unsigned int linecnt=0;                       // temporal variable which is not really used but demanded
	char c;                                       // comma

	// read the number of deviates
	IN >> n >> c >> ws;
	if(c!=',') error("%s '%s': The parameters have to be separated by a comma ','!\n", name.c_str(), t.c_str());
	if(n<1)    error("%s '%s': The number of deviates has to be positive!\n", name.c_str(), t.c_str());
	if(IN.eof()) error("%s '%s': Problems reading the arguments!\n", name.c_str(), t.c_str());

	// read first element
	if(IN.peek()=='{') mat1 = new TMatrix(STRING::readBrackets2String(IN, linecnt, '}'));
	else IN >> val1;
	IN >> c >> ws;
	if(c!=',') error("%s '%s': The parameters have to be separated by a comma ','!\n", name.c_str(), t.c_str());
	if(IN.eof()) error("%s '%s': Problems reading the arguments!\n", name.c_str(), t.c_str());

	// read second element
	if(IN.peek()=='{') mat2 = new TMatrix(STRING::readBrackets2String(IN, linecnt, '}'));
	else IN >> val2;
	IN >> ws;
	if(!IN.eof()) error("%s '%s': Problems reading the arguments!\n", name.c_str(), t.c_str());

	// create the sequence
	vector<string> vec;
	vec.reserve(n);
	if(mat1 || mat2){ // at least one parameter is specified by a matrix
		if(!mat1) mat1 = new TMatrix("{"+toStr(val1)+"}");  // get it as a matrix
		if(!mat2) mat2 = new TMatrix("{"+toStr(val2)+"}");  // get it as a matrix

		unsigned int nbCol1 = mat1->getNbCols();
		unsigned int nbCol2 = mat2->getNbCols();
		unsigned int nbRow1 = mat1->getNbRows();
		unsigned int nbRow2 = mat2->getNbRows();
		unsigned int nbCol  = nbCol1>nbCol2 ? nbCol1 : nbCol2;
		unsigned int nbRow  = nbRow1>nbRow2 ? nbRow1 : nbRow2;

		// get step size matrix and the starting matrix
		TMatrix matrix(nbRow, nbCol);     // contains the step size change
		unsigned int y, x;
		for(; n; --n){
			ostringstream out;
			out << '{';
			for(x=0; x<nbRow; ++x){
				out << '{';
				for(y=0; y<nbCol; ++y){
					if(y) out << ' ';
					out << (this->*pt2Func)((F)mat1->get(x % nbRow1, y % nbCol1), (S)mat2->get(x % nbRow2, y % nbCol2));
				}
				out << '}';
			}
			out << '}';
			vec.push_back(out.str());
		}
		delete mat1;
		delete mat2;
	}
	else{               // both parameters are defined by a simple number
		for(; n; --n){
			vec.push_back(toStr((this->*pt2Func)(val1, val2)));
		}
	}
	return vec;
}

// ----------------------------------------------------------------------------------------
// macro_4args
// ----------------------------------------------------------------------------------------
/** the text is reaplced by random deviates from a normal distribution for a macro with four arguments: Example:
 * rbeta(nb, a, b, min, max)
 * nb:    number of random deviates
 * a:  shape
 * b:    scale
 * min/max: range of the deviates
 *          if they are not set (just 2 argumetns, then then min is set to 0 and max is set to 1
 */
template<typename F, typename S, typename T1, typename T2>
vector<string>
Param::macro_4args(const string& t, const string& name, double (Param::*pt2Func)(const F&, const S&, const T1&, const T2&))
{
	istringstream IN;                             // make a new stream
	IN.str(t);                                    // allocate the line to the new stream
	F val1=my_NAN;
	S val2=my_NAN;
	T1 val3=my_NAN;
	T2 val4=my_NAN;                               // mean and var if they are numbers
	TMatrix* mat1=NULL, *mat2=NULL, *mat3=NULL, *mat4=NULL;               // mean and var if they are matrixes
	unsigned int n;                               // number of deviates
	unsigned int linecnt=0;                       // temporal variable which is not really used but demanded
	char c;                                       // comma

	// read the number of deviates
	IN >> n >> c >> ws;
	if(c!=',') error("%s '%s': The parameters have to be separated by a comma ','!\n", name.c_str(), t.c_str());
	if(n<1)    error("%s '%s': The number of deviates has to be positive!\n", name.c_str(), t.c_str());
	if(IN.eof()) error("%s '%s': Problems reading the arguments!\n", name.c_str(), t.c_str());

	// read first element
	if(IN.peek()=='{') mat1 = new TMatrix(STRING::readBrackets2String(IN, linecnt, '}'));
	else IN >> val1;
	IN >> c >> ws;
	if(c!=',') error("%s '%s': The parameters have to be separated by a comma ','!\n", t.c_str());
	if(IN.eof()) error("%s '%s': Problems reading the arguments!\n", name.c_str(), t.c_str());

	// read second element
	if(IN.peek()=='{') mat2 = new TMatrix(STRING::readBrackets2String(IN, linecnt, '}'));
	else IN >> val2;
	IN >> c >> ws;
	if(c!=',') error("%s '%s': The parameters have to be separated by a comma ','!\n", t.c_str());
	if(IN.eof()) error("%s '%s': Problems reading the arguments!\n", name.c_str(), t.c_str());

	// read third element
	if(IN.peek()=='{') mat3 = new TMatrix(STRING::readBrackets2String(IN, linecnt, '}'));
	else IN >> val3;
	IN >> c >> ws;
	if(c!=',') error("%s '%s': The parameters have to be separated by a comma ','!\n", t.c_str());
	if(IN.eof()) error("%s '%s': Problems reading the arguments!\n", name.c_str(), t.c_str());

	// read fourth element
	if(IN.peek()=='{') mat4 = new TMatrix(STRING::readBrackets2String(IN, linecnt, '}'));
	else IN >> val4;
	IN >> ws;
	if(!IN.eof()) error("%s '%s': Problems reading the arguments!\n", name.c_str(), t.c_str());

	// create the sequence
	vector<string> vec;
	vec.reserve(n);
	if(mat1 || mat2 || mat3 || mat4){ // at least one parameter is specified by a matrix
		if(!mat1) mat1 = new TMatrix("{"+toStr(val1)+"}");  // get it as a matrix
		if(!mat2) mat2 = new TMatrix("{"+toStr(val2)+"}");  // get it as a matrix
		if(!mat3) mat3 = new TMatrix("{"+toStr(val3)+"}");  // get it as a matrix
		if(!mat4) mat4 = new TMatrix("{"+toStr(val4)+"}");  // get it as a matrix

		unsigned int nbCol1 = mat1->getNbCols();
		unsigned int nbCol2 = mat2->getNbCols();
		unsigned int nbCol3 = mat3->getNbCols();
		unsigned int nbCol4 = mat4->getNbCols();
		unsigned int nbRow1 = mat1->getNbRows();
		unsigned int nbRow2 = mat2->getNbRows();
		unsigned int nbRow3 = mat3->getNbRows();
		unsigned int nbRow4 = mat4->getNbRows();
		unsigned int nbCol  = nbCol1>nbCol2 ? nbCol1 : nbCol2;
		if(nbCol<nbCol3) nbCol = nbCol3;
		if(nbCol<nbCol4) nbCol = nbCol4;
		unsigned int nbRow  = nbRow1>nbRow2 ? nbRow1 : nbRow2;
		if(nbRow<nbRow3) nbRow = nbRow3;
		if(nbRow<nbRow4) nbRow = nbRow4;

		// get step size matrix and the starting matrix
		TMatrix matrix(nbRow, nbCol);     // contains the step size change
		unsigned int y, x;
		for(; n; --n){
			ostringstream out;
			out << '{';
			for(x=0; x<nbRow; ++x){
				out << '{';
				for(y=0; y<nbCol; ++y){
					if(y) out << ' ';
					out << (this->*pt2Func)((F)mat1->get(x % nbRow1, y % nbCol1), (S)mat2->get(x % nbRow2, y % nbCol2),
							(T1)mat3->get(x % nbRow3, y % nbCol3), (T2)mat4->get(x % nbRow4, y % nbCol4));
				}
				out << '}';
			}
			out << '}';
			vec.push_back(out.str());
		}
		delete mat1;
		delete mat2;
		delete mat3;
		delete mat4;
	}
	else{               // both parameters are dfiend by a simple number
		for(; n; --n){
			vec.push_back(toStr((this->*pt2Func)(val1, val2, val3, val4)));
		}
	}
	return vec;
}

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                             ******** ParamSet ********

// ----------------------------------------------------------------------------------------
// ParamSet
// ----------------------------------------------------------------------------------------
ParamSet::ParamSet(string name, string name_long, bool isRequired, TMetapop* p) :
_isSet(0), _pParamManager(0)
{
    _popPtr = p;
    _name=name;
    _name_long = name_long;
    _isRequired=isRequired;
#ifdef _DEBUG
    message(" ParamSet::ParamSet(%s)\n", _name.c_str());
#endif
}

// ----------------------------------------------------------------------------------------
// ParamSet
// ----------------------------------------------------------------------------------------
ParamSet::ParamSet( ) : _isSet(0), _isRequired(0), _pParamManager(0), _popPtr(0)
{
#ifdef _DEBUG
    message(" ParamSet::ParamSet(not set)\n", _name.c_str());
#endif
}

// ----------------------------------------------------------------------------------------
// ~ParamSet
// ----------------------------------------------------------------------------------------
ParamSet::~ParamSet ()
{
#ifdef _DEBUG
    message(" ParamSet::~ParamSet(%s)\n", _name.c_str());
#endif
    map<string, Param*>::iterator param = _params.begin();
    for(; param != _params.end(); ++ param) {
        if(param->second) {
            // cout << "removeParam: " << param->first << endl;	// delete Param
            delete param->second;
        }
    }
    _params.clear();
    
    multimap<int, map<string, Param*>* >::iterator pos = _temporalParams.begin();
    for(; pos != _temporalParams.end(); ++pos){
        if(pos->second) delete pos->second;		// delete map<string, Param*>
    }
    _temporalParams.clear();
}

// ----------------------------------------------------------------------------------------
// reset
// ----------------------------------------------------------------------------------------
void
ParamSet::reset ()
{
    _isSet = false;
    
    map<string, Param*>::iterator param =  _params.begin();
    for(; param != _params.end(); ++ param) {
        param->second->reset();
    }
    
    multimap<int, map<string, Param*>* >::iterator pos = _temporalParams.begin();
    for(; pos != _temporalParams.end(); ++pos){
        if(pos->second) delete pos->second;
    }
    _temporalParams.clear();
}

// ----------------------------------------------------------------------------------------
// add_param
// ----------------------------------------------------------------------------------------
void
ParamSet::add_param (Param* param)
{
    assert(_params.find(param->get_name()) == _params.end());
    _params[param->get_name()] = param;
}

// ----------------------------------------------------------------------------------------
// add_param
// ----------------------------------------------------------------------------------------
void
ParamSet::add_param (string Name,param_t Type,bool mandatory,
                     double low_bnd,double up_bnd,string def,bool temp,
                     string text, unsigned int importance)
{
    // control that not another object is overridden if the names are identical
    // cout << "addParam: " << Name << endl;
    assert(_params.find(Name) == _params.end());
    _params[Name] = new Param(Name, Type, mandatory, low_bnd, up_bnd, def, temp,
                              text, importance);
}

// ----------------------------------------------------------------------------------------
// set_param
// ----------------------------------------------------------------------------------------
bool
ParamSet::set_param (const string& Name, const string& Arg, TReplicate* pRep)
{
    map<string, Param*>::iterator param = _params.find(Name);
    if(param != _params.end()) return param->second->set(Arg, this, pRep);
    else return false;
}

// ----------------------------------------------------------------------------------------
// check_param_name
// ----------------------------------------------------------------------------------------
/** returns true if the passed name is a prameter name, otherwise false */
bool
ParamSet::check_param_name (const string& Name)
{
    return (_params.find(Name) != _params.end());
}

// ----------------------------------------------------------------------------------------
// set_temporal_param
// ----------------------------------------------------------------------------------------
void
ParamSet::set_temporal_param(const int& gen, Param* parm)
{
    multimap<int, map<string, Param*>* >::iterator pos = _temporalParams.find(gen);
    if(pos != _temporalParams.end()){
        pos->second->insert(pair<string, Param*>(parm->get_name(), parm));      // add the new parm
    }
    else{        // this generation has not yet been used
        map<string, Param*>* secMap = new map<string, Param*>;                  // make the second map with new
        secMap->insert(pair<string, Param*>(parm->get_name(), parm));           // set second map
        _temporalParams.insert(pair<int, map< string, Param*>* >(gen, secMap)); // set multimap
    }
    
}

// ----------------------------------------------------------------------------------------
// getTemporalParams
// ----------------------------------------------------------------------------------------
/** returns the map if the passed generation time is available, NULL if not */
map<string, Param*>*
ParamSet::getTemporalParams(const int& gen)
{
    multimap<int, map<string, Param*>* >::iterator pos;
    pos = _temporalParams.find(gen);
    
    // if present return the map, otherwise return NULL
    if(pos != _temporalParams.end()) return pos->second;
    return NULL;
}

// ----------------------------------------------------------------------------------------
// updateTemporalParams
// ----------------------------------------------------------------------------------------
/** updates the current value by the temporal value
 * returns the map with the params if a paramter was updated, NULL if no parameter was updated
 */
map<string, Param*>*
ParamSet::updateTemporalParams(const int& gen)
{
    multimap<int, map<string, Param*>* >::iterator pos;
    pos = _temporalParams.find(gen);
    
    // if not present go on
    if(pos == _temporalParams.end()) return NULL;
    
    // if present update the parameters and return true
    map<string, Param*>::iterator map_pos = pos->second->begin();
    for(; map_pos != pos->second->end(); ++map_pos){
        map_pos->second->update_arg(gen);
    }
    return pos->second;
}

// ----------------------------------------------------------------------------------------
// get_param
// ----------------------------------------------------------------------------------------
/** returns a pointer to the param if found and an error if not found */
Param*
ParamSet::get_param (string Name)
{
    map<string, Param*>::iterator param = _params.find(Name);
    assert(param != _params.end());
    return param->second;
}

// ----------------------------------------------------------------------------------------
// find_param
// ----------------------------------------------------------------------------------------
/** returns a pointer to the param if found and NULL if not found */
Param*
ParamSet::find_param (string Name)
{
    map<string, Param*>::iterator param = _params.find(Name);
    if(param == _params.end()) return NULL;
    return param->second;
}
// ----------------------------------------------------------------------------------------
// check_consistency
// ----------------------------------------------------------------------------------------
bool
ParamSet::check_consistency ()
{
    map<string, Param*>::iterator param = _params.begin();
    bool isOK = true;
    // bool touched = false;
    bool firstOK = true;  // show only the first error
    
    for(int i=0; param != _params.end(); ++param, ++i) {
        //check if all required fields have been set properly
        if(!i)firstOK = param->second->isSet();
        if(param->second->isRequired()) {
            isOK &= param->second->isSet();
            if(!param->second->isSet() && firstOK){
                firstOK = false;
                warning("Parameter set '%s' is disabled as parameter '%s' is required but not set!\n",
                        _params.begin()->second->get_name().c_str(), param->second->get_name().c_str());
            }
        }
    }
    _isSet = isOK;
    //return isOk or check if _isRequired in case no params are set (untouched params)
    return ( isOK | (!_isRequired));// & !touched) );
}
// ----------------------------------------------------------------------------------------
// show_up
// ----------------------------------------------------------------------------------------
void
ParamSet::show_up ()
{
    message("\n%s",_name.c_str());
    map<string, Param*>::iterator param = _params.begin();
    while(param != _params.end()) {
        param->second->show_up();
        param++;
    }
}
// ----------------------------------------------------------------------------------------
// print
// ----------------------------------------------------------------------------------------
void
ParamSet::print (ostream& FILE, bool macros)
{
    map<string, Param*>::iterator cur, end= _params.end();
    
    // test if the ParamSet has set params
    for(cur=_params.begin(); cur!=end; ++cur) {
        if(cur->second->isSet()) break;
    }
    if(cur==end) return;    // noting is set: stop here
    
    
    FILE << "\n## " << _name << "\n";
    for(cur=_params.begin(); cur!=end; ++cur) {
        if(!cur->second->isSet()) continue;
        FILE << cur->second->get_name() << " ";
        if(macros) FILE << cur->second->get_tot_arg();
        else       FILE << cur->second->get_tot_arg_no_macro();
        FILE << "\n";
    }
}

// ----------------------------------------------------------------------------------------
// print
// ----------------------------------------------------------------------------------------
void
ParamSet::print_minimal (ostream& FILE, bool macros)
{
    map<string, Param*>::iterator cur, end= _params.end();
    
    // test if the ParamSet has set params
    for(cur=_params.begin(); cur!=end; ++cur) {
        if(cur->second->isSet() && cur->second->get_tot_arg() != cur->second->get_default_arg()) break;
    }
    if(cur==end) return;    // noting is set: stop here
    
    
    FILE << "\n## " << _name << "\n";
    for(cur=_params.begin(); cur!=end; ++cur) {
        if(!cur->second->isSet() || cur->second->get_tot_arg() == cur->second->get_default_arg()) continue;
        FILE << cur->second->get_name() << " ";
        if(macros) FILE << cur->second->get_tot_arg();
        else       FILE << cur->second->get_tot_arg_no_macro();
        FILE << "\n";
    }
}

// ----------------------------------------------------------------------------------------
// print
// ----------------------------------------------------------------------------------------
void
ParamSet::print_maximal (ostream& FILE, bool macros)
{
    string name;
    
    // name of the param set
    FILE << "\n#### " << getName() << " ####\n";
    
    map<string, Param*>::iterator param = _params.begin();
    for(;param != _params.end(); ++param) {
        // parameter + argument
        name = param->second->get_name();
        if(name[0]=='_') continue;
        name += + "  ";
        if(param->second->isSet()){
            if(macros) name += param->second->get_tot_arg();
            else       name += param->second->get_tot_arg_no_macro();
        }
        else         name += param->second->get_default_arg();
        FILE.width(40);
        FILE.setf(ios::left,ios::adjustfield);
        FILE << name;
        
        // argument type
        FILE << "# type: ";
        name = param->second->get_type_str() + ";";
        FILE.width(16);
        FILE << name;
        
        // default value
        if(param->second->get_default_arg().empty()) name = "\"\";";
        else name = param->second->get_default_arg() + ";";
        FILE << " default: ";
        FILE.width(16);
        FILE << name;
        
        if(param->second->get_type()!=STR && param->second->get_type()!=STR_MAT){
            // temporal parameter
            FILE.width(11);
            if(param->second->isTemporalArgumentAllowed()) FILE << " temporal; ";
            else                                           FILE << "";
            
            // range
            double min = param->second->get_bound(0);
            double max = param->second->get_bound(1);
            if(min != my_NAN) FILE << " range: [" << min << "; ";
            else              FILE << " range: ]-; ";
            if(max != my_NAN) FILE << max << "];";
            else              FILE << "-[;";
        }
        
        FILE << "\n";
    }
}

// ----------------------------------------------------------------------------------------
// operator=
// ----------------------------------------------------------------------------------------
ParamSet&
ParamSet::operator=(const ParamSet& i)
{
    if(this != &i) {
        _name = i._name;
        _name_long = i._name_long;
        _isSet = i._isSet;
        _isRequired = i._isRequired;
        _popPtr = i._popPtr;
        
        _temporalParams = i._temporalParams;    // contains pointers, so that is fine
        
        map<string, Param* >::const_iterator cur, end;
        for(cur=i._params.begin(), end=i._params.end(); cur!=end; ++cur){
            _params[cur->first] = new Param(*cur->second);
        }
        
        multimap<int, map<string, Param*>* >::const_iterator curT, endT;
        for(curT=i._temporalParams.begin(), endT=i._temporalParams.end(); curT!=endT; ++curT){
            map<string, Param*>* secMap = new map<string, Param*>;   // make the second map with new
            _temporalParams.insert(pair<int, map< string, Param*>* >(curT->first, secMap)); // set multimap
            
            for(cur=curT->second->begin(), end=curT->second->end(); cur!=end; ++cur){
                assert(_params.find(cur->first) != _params.end());
                secMap->insert(pair<string, Param*>(cur->first, _params[cur->first]));
            }
        }
        
        
    }
    return *this;
}


// ----------------------------------------------------------------------------------------
// print_help
// ----------------------------------------------------------------------------------------
/** output all parameters with their description 
 * print first the more important parameters and then the les simprotant ones until 
 * "importance"
 * if arg is "": print only the paramSet names
 * if arg == the current _name_long: print the parameters of this ParamSet
 * if arg == "all": print all parameters odf all ParamSets
 */
void
ParamSet::print_help(ostream& os, unsigned int wide1, unsigned int wide2, char fill,
                         unsigned int importance, string arg)
{
    if(!_params.empty() && arg==""){
        os << "\n         " << left << setw(22) << _name << "(" << _name_long << ")";
        return;
    }
    else if(_params.empty() ||
            (_params.size()==1 && _params.begin()->first==_name) ||
            (arg!="all" && arg!=_name)) return;
    
    os << "\n\n\n" << _name_long;
    os << "\n" << setw((unsigned int)_name_long.length()) << setfill('=') << "";
    
    for(unsigned int i=0; i<=importance; ++i){
       map<string, Param*>::iterator cur = _params.begin(), end = _params.end();
        for(; cur!=end; ++cur){
            if(cur->second->get_name()==_name) continue;
            cur->second->print_help(os, wide1, wide2, fill, i);
        }
    }
    
}

