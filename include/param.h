/** @file param.h
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

#ifndef paramH
#define paramH

#include <map>
#include "ttree.h"
#include "tmatrix.h"
using namespace std;

const unsigned int NB_MACROS = 18;

/**This structure stores one parameter, its definition and its string argument. **/
class TReplicate;
class ParamSet;
class Param
{
    
private:
    string          _name;
    string          _arg;         // contains the current argument
    string          _tot_arg;     // contains the argument as passed to quantiNemo2
    string          _tot_arg_no_macro;  // as _tot_arg_no_keyword but macros were replaced ("" if no macro)
    string          _default_arg; // contains the default argument
    param_t         _type;
    bool            _isSet;
    bool            _isRequired;
    double          _bounds[2];   // my_NAN stands for not bounded at this side
    
    ParamSet*        _parSet;
    
    // parameters for temporal changes
    bool             _temporalParamAllowed;
    map<int, string> _temporalArgs; // _temporalArgs[generation]=_tot_arg(for the current time)
    
    // parameters for macros
    bool             _isRandom;
    TReplicate*      _pReplicate;
    map<int, string> _common_random_number; // map[generation] = arg
    
    string           _description;
    unsigned int     _importance;       // how important is the parameter to be plotted to the help
                                        // 0: very important
                                        // 1: important
                                        // 2: less improtant
                                        // 20: hacks
    
    
public:
    /**
     * @param Name the name of the parameter as read in the init file
     * @param Type the type of the parameter argument (see types.h), used to convert the argument string into a value
     * @param mandatory specifies whether this parameter is mandatory for the ParamSet owning it
     * @param low_bnd the lower bound of the parameter value
     * @param up_bnd the upper bound
     * @param temp can tis parameter change over time?
     **/
    Param  (string& Name,param_t Type,bool mandatory = 0,
            double low_bnd = 0,double up_bnd = 0, string def_arg = "",
            bool temp = 0, string text="", unsigned int importance=5);
    
    Param(const Param& i);
    
    ~Param () {	}
    
    /**@brief clears the _isSet flag and the argument string**/
    void    reset   ();
    
    //accessors:
    /**@brief Sets the parameter's name**/
    void            set_name             (string value)          {_name = value;}
    
    /**@brief Sets the parameter's argument**/
    void            set_arg              (string value)          {_arg = value;}
    void            set_tot_arg          (string value)          {_tot_arg = value;}
    
    /**@brief Sets the parameter's type (see types.h)**/
    void            set_type             (param_t value)         {_type = value;}
    
    /**@brief Sets the _isSet flag**/
    void            set_isSet            (bool value)            {_isSet = value;}
    
    /** the argument is updated by the new temporal argument */
    void            update_arg   (const int& gen);
    void            set_arg_for_gen(int gen);
    
    // random variables
    double Uniform();
    double Uniform(const double& min, const double& max);
    double Poisson(const double& mean);
    double Normal(const double& mean, const double& sd);
    double LogNormal(const double& mean, const double& sd);
    double Gamma(const double& a, const double& b);
    double Binomial(const unsigned int& size, const double& prob);
    double Beta(const double& a, const double& b);
    double Beta(const double& a, const double& b, const double& min, const double& max);
    double Product(const double& a, const double& b){return a*b;}
    
    void Uniform(string& s, double min=my_NAN, double max=my_NAN);
    void Poisson(string& s, double min=my_NAN, double max=my_NAN);
    void Normal(string& s, double min=my_NAN, double max=my_NAN);
    void LogNormal(string& s, double min=my_NAN, double max=my_NAN);
    void Gamma(string& s, double min=my_NAN, double max=my_NAN);
    void Binomial(string& s, double min=my_NAN, double max=my_NAN);
    void Beta(string& s, double min=my_NAN, double max=my_NAN);
    void Product(string& s, double min=my_NAN, double max=my_NAN);
    
    string vector2string(const vector<string>& vec);
    string vector2string(const vector<string>& vec, const vector<double>& time_vec);
    
    // macro replace functions
    typedef void (Param::*_func_ptr)(string& s, double min, double max);
    _func_ptr _macros_func_ptr[NB_MACROS];
    string _macros[NB_MACROS];
    
    string replace_macros(string& text);
    vector<string> macro_0args(const string& t, const string& name, double (Param::*pt2Func)());
    template<typename F>
    vector<string> macro_1args(const string& t, const string& name, double (Param::*pt2Func)(const F&));
    template<typename F, typename S>
    vector<string> macro_2args(const string& t, const string& name, double (Param::*pt2Func)(const F&, const S&));
    template<typename F, typename S, typename T1, typename T2>
    vector<string> macro_4args(const string& t, const string& name, double (Param::*pt2Func)(const F&, const S&, const T1&, const T2&));
    vector<double> get_seq(double from, const double& to, const unsigned int& by, bool asInt=false);
    vector<TMatrix> get_seq(const TMatrix& from, const TMatrix& to, const unsigned int& by, bool asInt=false);
    vector<double>  get_logistic(const double& from, const double& to, const double& r, const double r_max,
                                 const double& s, const unsigned int& by, bool asInt=false);
    vector<TMatrix> get_logistic(const TMatrix& fromMatrix, const TMatrix& toMatrix,
                                 const TMatrix& rMatrix,    const TMatrix& r_maxMatrix,
                                 const TMatrix& sMatrix, const unsigned int& by, bool asInt=false);
    
    vector<string> rsample_vec(const string& t);
    vector<string> seq_vec(const string& t);
    vector<string> seq2D_vec(const string& t);
    vector<string> seq2Db_vec(const string& t);
    vector<string> rep_vec(const string& t);
    vector<string> logistic_vec(const string& t);
    
    void seq(string& s, double min=my_NAN, double max=my_NAN);
    void seq2D(string& s, double min=my_NAN, double max=my_NAN);
    void seq2Db(string& s, double min=my_NAN, double max=my_NAN);
    void product(string& s, double min=my_NAN, double max=my_NAN);
    void rep(string& s, double min=my_NAN, double max=my_NAN);
    void logistic(string& s, double min=my_NAN, double max=my_NAN);
    void runif(string& s, double min=my_NAN, double max=my_NAN);
    void rnorm(string& s, double min=my_NAN, double max=my_NAN);
    void rgamma(string& s, double min=my_NAN, double max=my_NAN);
    void rbeta(string& s, double min=my_NAN, double max=my_NAN);
    void rpois(string& s, double min=my_NAN, double max=my_NAN);
    void rlnorm(string& s, double min=my_NAN, double max=my_NAN);
    void rbinom(string& s, double min=my_NAN, double max=my_NAN);
    void rsample(string& s, double min=my_NAN, double max=my_NAN);
    void round(string& s, double min=my_NAN, double max=my_NAN);
    void ceil(string& s, double min=my_NAN, double max=my_NAN);
    void floor(string& s, double min=my_NAN, double max=my_NAN);
    void trunc(string& s, double min=my_NAN, double max=my_NAN);
    
    
    string&         get_name             ()                      {return _name;}
    string&         get_arg              ()                      {return _arg;}
    string&         get_default_arg      ()                      {return _default_arg;}
    string&         get_tot_arg          ()                      {return _tot_arg;}
    string&         get_tot_arg_no_macro () {return _tot_arg_no_macro.empty() ?  _tot_arg : _tot_arg_no_macro;}
    map<int,string>* get_temporal_args   ()                      {return &_temporalArgs;}
    param_t         get_type             ()                      {return _type;}
    bool            isSet                ()                      {return _isSet;}
    bool            isRequired           ()                      {return _isRequired;}
    double          get_bound            (unsigned int i)        {return _bounds[i];}
    bool            isTemporalParam      ()                      {return !_temporalArgs.empty();}
    bool            isTemporalArgumentAllowed()                  {return _temporalParamAllowed;}
    
    /**@brief Sets the _isSet flag to true and _arg to arg if the arg is of the right type and whithin the bounds**/
    bool            set                  (string arg, ParamSet* parSet, TReplicate* pRep);
    bool            set_new_arg          (string arg);
    string          get_new_arg          (string arg);
    string          set_temporal_args    (const string& arg);   // set temporal arguments
    void            set_up_main_random_param(string arg);
    bool            check_arg            (string arg);
    
    /**@brief Get the argument value according to its type
     *@return my_NAN if the parameter is not set or not of a the right type
     **/
    double          get_value            ();
    
    /**@brief Checks if the argument is of matrix type**/
    bool            is_matrix            () {return (_type == MAT || (_arg.size() ? _arg[0] == '{' : false));}
    
    
    unsigned int    matrix_depth         () {
        unsigned int i=0, nb=0;
        unsigned int argSize=(unsigned int)_arg.size();
        char c;
        for(;i<argSize; ++i){
            c=_arg[i];
            if(c=='{') {++nb; continue;}
            if(isspace(c)) continue;
            break;
        }
        return nb;
    }
    
    template<class T>
    TTree<unsigned int, T>* get_tree(){
        try{
            if(_arg[0] != '{')_arg = "{" + _arg + "}"; // add brackets around the string if they are missing
            return new TTree<unsigned int, T>(_arg);
        }
        catch(const char* err){
            error("Parameter '%s': %s\n", _name.c_str(), err);
        }
        catch(...){
            error("Parameter '%s': Could not parse the tree!\n", _name.c_str());
        }
        return NULL;
    }
    
    
    
    /**@brief Gets the matrix from the argument string if the parameter is set and of matrix type
     *@return 0 if parameter not set or of the wrong type
     **/
    TMatrix*        get_matrix           ();
    TMatrix*        get_as_matrix        (bool allowStr = false);  // returns a matrix even if a single value was passed
    
    /**@brief Checks if the argument is of matrix type**/
    bool            is_matrixVar() {return (_type == MAT_VAR || (_arg.size() ? _arg[0] == '{' : false) );}
    
    /**@brief Gets the matrix from the argument string if the parameter is set and of matrix type
     *@return 0 if parameter not set or of the wrong type
     **/
    
    TMatrixVar<double>* get_matrixVar ()
    {
        if( is_matrixVar() && _isSet ){
            try{
                return new TMatrixVar<double>(_arg);
            }
            catch(const char* err){
                error("Parameter '%s': %s\n", _name.c_str(), err);
            }
            catch(...){
                error("Parameter '%s': Could not parse the matrix!\n", _name.c_str());
            }
        }
        return NULL;
    }
    
    TMatrixVar<string>*     get_matrixVarStr ();
    
    /**@brief Print state to stdout**/
    void            show_up              ();
    void            print_help           (ostream& os, unsigned int wide1=20,
                                          unsigned int wide2=80, char fill='.',
                                          unsigned int level=0);
    string          get_type_str         ();
    string          replaceKeywords      (string arg);
    string          replaceMacros        (string arg);
    
    ParamSet*       getParamSet          ()                      {return _parSet;}
    
};


/**Parameters container, implemented in each SimComponent. **/
class TMetapop;
class ParamManager;
class ParamSet
{
public:
    string                  _name;
    string                  _name_long;
    bool                    _isSet;
    bool                    _isRequired;
    
    ParamManager*           _pParamManager;
    
    map<string, Param*>     _params;    // Param* has to be deleted
    
    
    // temporal parameter (don't delete the Param* object!)
    multimap<int, map<string, Param*>* > _temporalParams;
    
public:
    TMetapop*                _popPtr;    // used for random numbers
    
    ParamSet                       ( );
    ParamSet                       (string name, string name_long, bool isRequired, TMetapop* p);
    ~ParamSet                      ( );
    
    /**@brief Put the container in the unset state, reset each Param it contains**/
    void            reset                ( );
    
    /**@brief Sets the container's name**/
    void            set_name             (string value)          {_name = value;}
    void            set_name_full        (string value)          {_name_long = value;}
    
    void            set_pParamManager(ParamManager* p){_pParamManager=p;}
    ParamManager*   get_pParamManager(){return _pParamManager;}
    
    /**@brief Sets the _isRequired flag meaning this container is mandatory and must be set in order to run a simulation**/
    void            set_isRequired       (bool value)            {_isRequired = value;}
    
    /**@brief Checks for the status of the required parameters
     *@return TRUE if all required parameters are or if nothing is set and the container is not required
     **/
    bool            check_consistency    ( );
    
    /**@brief print info to stdout**/
    void            show_up              ( );
    
    /**@brief print all set parameters to the outpout file stream**/
    void            print                (ostream& FILE, bool macros=1);
    void            print_minimal        (ostream& FILE, bool macros=1);
    void            print_maximal        (ostream& FILE, bool macros=1);
    
    /**@brief Returns the number of parameters contained**/
    unsigned  int             size                 ( )                     {return (unsigned int)_params.size();}
    
    /**@brief Returns the complete list of parameters**/
    map<string, Param*>& getAllParams   ( )                     {return _params;}
    
    /** returns the map if the passed generation time is available, NULL if not */
    multimap<int, map<string, Param*>* >* getTemporalParams() {return &_temporalParams;}
    map<string, Param*>* getTemporalParams(const int& gen);
    map<string, Param*>* updateTemporalParams(const int& gen);
    
    ///@name Accessors to Param members.
    ///@{
    /**@brief Adds the param argument to the list**/
    void            add_param            (Param* param);
    
    /**@brief Adds a new param specified by arguments to the list.
     *@param Name the name of the parameter
     *@param Type the type of the parameter
     *@param mandatory specifies if this parameter is required and must be set for the container to gain the "set" status
     *@param low_bnd the lower value the parameter can take
     *@param up_bnd the upper value the parameter can take
     *@param default_param the default value
     *@param temp can the parameter change over time?
     **/
    void            add_param            (string Name,param_t Type,bool mandatory,
                                          double low_bnd,double up_bnd,
                                          string default_param, bool temp=false,
                                          string text="", unsigned int importance=5);
    
    /**@brief Look for a param named "Name" and try to set it with the "Arg" argument string.
     *@return TRUE if param Name has been found and set with Arg
     *@return FALSE otherwise
     *@param Name the name of the parameter to find in the list
     *@param Arg the argument string as found in the init params
     **/
    bool            set_param            (const string& Name, const string& Arg, TReplicate* pRep);
    bool            check_param_name     (const string& Name);
    void            set_temporal_param   (const int& gen, Param* parm);  // for temporal chainging parameters
    bool            set_general_parameter(map< string,string >* input);
    
    /**@brief Look for a param "name" in its parameters list.
     *@return NULL if no Param with _name = name exists
     **/
    Param*         get_param            (string name);     // returns an error if not found
    Param*         find_param           (string name);     // returns NULL if not found
    map<string, Param* >* get_params(){return &_params;}
    
    /**@brief Accessor to the status flag. **/
    bool            isSet                ()                      {return _isSet;}
    
    /**@brief Accessor to the parameters status flag. **/
    bool            isSet                (string name)           {return (get_param(name))->isSet();}
    
    /**@brief Accessor to the parameters status flag. **/
    void            set_isSet            (bool choice)           {_isSet = choice;}
    
    /**@brief Accessor to the mandatory flag. **/
    bool            isRequired           ()                      {return _isRequired;}
    
    bool            isTemporalParamSet   ()                      {return !_temporalParams.empty();}
    
    /**@brief Check if the parameter "name" is of matrix type. **/
    bool            isMatrix             (string name)           {return (get_param(name))->is_matrix();}
    unsigned int    matrixDepth          (string name)           {return (get_param(name))->matrix_depth();}
    
    /**@brief Name accessor.**/
    string&         getName              ()                      {return _name;}
    string&         getArg               (string name)           {return (get_param(name))->get_arg();}
    double          getValue             (string name)           {return (get_param(name))->get_value();}
    TMatrix*        getMatrix            (string name)           {return (get_param(name))->get_matrix();}
    TMatrix*        getAsMatrix          (string name)           {return (get_param(name))->get_as_matrix();}
    
    ParamSet& operator=(const ParamSet& i);
    
    void print_help(ostream& os, unsigned int wide1=20, unsigned int wide2=60,
                    char fill='.',  unsigned int importance=5, string arg="");
    
    ///@}
};


////////////////////////////////////////////////////////////////////////////////

#endif
