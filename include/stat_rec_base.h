/** @file stat_rec_base.h
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

#ifndef stat_rec_baseH
#define stat_rec_baseH

#include "functions.h"

using namespace std;

class TStat_db;
class StatRecBaseAll;
class ParamRecBaseAll;
class TMetapop;

/**Base class for the StatRecorder's, stores the stat values in a matrix.
 * The way the stat values are stored in the _val matrix is set by the _ordering flag (see StatRecorder::setValDB()).
 * one stat per object
 * **/
class StatRecBase {
    
protected:
    /**The title of the stat recorder, longer and more explicite than the name.*/
    string  _title;
    
    /**Name of the stat, should be short (20 char) and R compliant (no '-', '+', ' ')*/
    string  _name;
    
    /**An argument to be passed to one of the function variable stored in the StatRecorder structure.*/
    unsigned int _arg;
    
    /**The age class for which this stat applies.*/
    age_t        _age;
    
    /**A flag specifying the way the stat values are recorded in the matrix, is initialized to FLAT by default.*/
    st_order     _ordering;
    
    
    StatRecBaseAll*  _pStatRecBaseAll;   // pointer to the ALL obejct
    
    
public:
    
    StatRecBase  ( ) : _title(), _name(), _arg(0), _age(ALL),
    _ordering(FLAT){}
    
 //   StatRecBase  (unsigned int size)
 //   { init(); }
    
    ~StatRecBase ( );
    
    /**Creates the _val matrix if its dimensions have been set and sets each elements to "NaN".*/
    void init    ( );
    
    /**Sets the recorder attributes.
     * @param T the stat title
     * @param N the stat name (headers in the output file)
     * @param Order stat table ordering flag
     * @param AGE the age class for which the stat will be recorded
     * @param ARG the argument to pass to the S function
     **/
    void set     (string T, string N, st_order Order, age_t AGE, unsigned int ARG);
    
    ///@name Accessors
    ///@{
    void setName                      (string N)            {_name = N;}
    
    /**@param i the number of rows
     @param j the number of columns */
    // void set_pStat                    (TStat_db* p)               {_pStat_db = p;}
    st_order     getOrdering           ( )                        {return _ordering;}
    string       getTitle              ( )                        {return _title;}
    string       getName               ( )                        {return _name;}
    age_t        getAge                ( )                        {return _age;}
    unsigned int getArg                ( )                        {return _arg;}
    
    double   getMean(const unsigned int& index) {return (this->*getMean_func_ptr)(index);}
    double   (StatRecBase::*getMean_func_ptr)(const unsigned int& index);
    double   getMedian(const unsigned int& index) {return (this->*getMedian_func_ptr)(index);}
    double   (StatRecBase::*getMedian_func_ptr)(const unsigned int& index);
    double   getVar(const unsigned int& index) {return (this->*getVar_func_ptr)(index);}
    double   (StatRecBase::*getVar_func_ptr)(const unsigned int& index);
    double   getMean_FLAT             (const unsigned int& gen_index);
    double   getMedian_FLAT           (const unsigned int& gen_index);
    double   getVar_FLAT              (const unsigned int& gen_index);
    double   get_GEN                  (const unsigned int& gen_index);
    double   get_RPL                  (const unsigned int& rpl_index);
    
    //    unsigned int getIndex             (unsigned int i)           {return _index[i];}
    
    
    void set_pStatRecBaseAll(StatRecBaseAll* a){_pStatRecBaseAll = a;}
    
    ///@}
    //    bool   (StatRecBase::*setVal_func_ptr)(int crnt_gen, int rpl_cntr, double value);
    //    bool   setVal_FLAT                (int crnt_gen, int rpl_cntr, double value);
    //    bool   setVal_GEN                 (int crnt_gen, int rpl_cntr, double value);
    //    bool   setVal_RPL                 (int crnt_gen, int rpl_cntr, double value);
    
};

////////////////////////////////////////////////////////////////////////////////
/**Stores the pointers to the StatHandler's stat functions.*/
template <class S>
class StatRecorder : public StatRecBase
{
private:
    /**Pointer to a 'stat getter' function of S using no argument.*/
    double (S::* _getStat)     		(void);
    double (S::* _getStatAGE)     (const age_t& AGE);
    
    /**Pointer to a 'stat getter' function of S using a bool argument.*/
    double (S::* _getStatBool) 		(bool);
    double (S::* _getStatBoolAGE) (bool, const age_t& AGE);
    
    /**Pointer to a 'stat getter' function of S using a unsigned int argument.**/
    double (S::* _getStatUI)   		(unsigned int);
    double (S::* _getStatUIAGE) 	(unsigned int, const age_t& AGE);
    
    /**The recording matrix, stores the stat values in a [_rows X _cols] matrix, i.e. first indice for row, second for column.*/
    double*     _stat;          // _stat[generation] // is not deleted here but in TStat_db
    
public:
    
    
    TStat_db* _pStat_db;
    
    StatRecBaseAll* _recAll;
    
    
    StatRecorder() : _getStat(0),    _getStatBool(0),    _getStatUI(0),
    _getStatAGE(0), _getStatBoolAGE(0), _getStatUIAGE(0), _stat(0), _recAll(0){}
    
    /**
     * @param n_rows number of rows of the stat table (usually nbr of records per replicates)
     * @param n_cols number of columns of the stat table (usually nbr of replicates)
     **/
    StatRecorder(TStat_db* pStat_db);
    
    //  ~StatRecorder() {}
    /**@brief sets the recorder attributes
     * @param Title 						the stat title
     * @param Name 							the stat name (headers in the output file)
     * @param Order 						stat table ordering flag
     * @param AGE 							age on which the stat should be processed
     * @param ARG 							the argument to pass to the S function
     * @param getSt 						function ptr to a S getter
     * @param getStBoolArg 			function ptr to a S getter with boolean argument
     * @param getStUintArg 			function ptr to a S getter with unsigned int argument
     * @param getStAGE 					function ptr to a S getter with age argument
     * @param getStBoolArgAGE 	function ptr to a S getter with boolean and age argument
     * @param getStUintArgAGE		function ptr to a S getter with unsigned int and age argument
     **/
    void set(string Title,string Name,st_order Order,age_t AGE,unsigned int ARG,
             double(S::* getSt)(void),
             double(S::* getStBoolArg)(bool)=0,
             double(S::* getStUintArg)(unsigned int)=0,
             double(S::* getStAGE)(const age_t&)=0,
             double(S::* getStBoolArgAGE)(bool, const age_t&)=0,
             double(S::* getStUintArgAGE)(unsigned int, const age_t&)=0);
    
    /**@brief launch the linked stat function and store the result in the stat table
     * @param AGE age on which the stat should be processed
     * @param crnt_gen the current running generation
     * @param rpl_cntr the current running replicate
     * @param StatHandler an instance of S
     **/
    void setStats (age_t AGE,int crnt_gen,int rpl_cntr, S* StatHandler, TMetapop* pop);  // to db
    void setStats (age_t AGE,int crnt_gen,int rpl_cntr, S* StatHandler,  ostream& FH);  // to the stream/file
   
    double*         get_stat(){return _stat;}
    void            set_stat(double* a){_stat = a;}
    void            init_stat_rec(TMetapop* pop);
    void            init_param_rec(TMetapop* pop);
    
};
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

// ----------------------------------------------------------------------------------------
// StatRecorder::set()
// ----------------------------------------------------------------------------------------
template<class S>
void StatRecorder<S>::set(string T, string N, st_order Order, age_t AGE,unsigned int ARG,
                          double (S::* getSt) (void),
                          double (S::* getStBoolArg) (bool),
                          double (S::* getStUintArg)(unsigned int),
                          double (S::* getStAGE) (const age_t&),
                          double (S::* getStBoolArgAGE) (bool, const age_t&),
                          double (S::* getStUintArgAGE)(unsigned int, const age_t&))
{
    StatRecBase::set(T, N, Order, AGE, ARG);
    
    _getStat        = getSt;
    _getStatBool    = getStBoolArg;
    _getStatUI      = getStUintArg;
    _getStatAGE     = getStAGE;
    _getStatBoolAGE = getStBoolArgAGE;
    _getStatUIAGE   = getStUintArgAGE;
}

// ----------------------------------------------------------------------------------------
// StatRecorder::setStats
// ----------------------------------------------------------------------------------------
/** this function writes the stata directly to the file (file is just opended for the output and then is closed again */
template<class S>
void StatRecorder<S>::setStats(age_t AGE, int crnt_gen, int rpl_cntr, S* StatHandler, ostream& FH)
{
    double statValue=0;
    age_t age = getAge();
    
    // test if the stat has to be computed for this age
    if(age & AGE) {
#ifdef _DEBUG
        message(" %s\n",getName().c_str());
#endif
        
        //get the value:
        if(_getStat)           		statValue = (StatHandler->*_getStat)       ();
        else if(_getStatBool)  		statValue = (StatHandler->*_getStatBool)   ((bool)getArg());
        else if(_getStatUI)    		statValue = (StatHandler->*_getStatUI)     (getArg());
        else if(_getStatAGE)      	statValue = (StatHandler->*_getStatAGE)    (age);
        else if(_getStatBoolAGE) 	statValue = (StatHandler->*_getStatBoolAGE)((bool)getArg(), age);
        else if(_getStatUIAGE)    	statValue = (StatHandler->*_getStatUIAGE)  (getArg(), age);
        else error("StatRecorder::setVal2file: no _getStat fct ptr !!\n");
        
        // write the stat directly to the file
        FH<<"\t";
        FH.width(12);
        if(statValue==my_NAN) FH << my_NANstr;
        else                  FH << statValue;
    }
}






#endif

