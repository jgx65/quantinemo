/** @file tstat_db.h
 *
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


#ifndef tstat_dbH
#define tstat_dbH


#include <string>
#include <map>
#include <list>
#include "types.h"
#include "tarray.h"
#include <vector>
using namespace std;


class TStat_db;
class StatRecBase;
class Param;

class StatRecBaseAll {
private:
public:
    string  _title;         /**The title of the stat recorder, longer and more explicite than the name.*/
    string  _name;          /**Name of the stat, should be short (20 char) and R compliant (no '-', '+', ' ')*/
    unsigned int _arg;      /**An argument to be passed to one of the function variable stored in the StatRecorder structure.*/
    age_t        _age;      /**The age class for which this stat applies.*/
    st_order     _ordering; /**A flag specifying the way the stat values are recorded in the matrix, is initialized to FLAT by default.*/
    double**     _stat;          // _val[replicate][generation]
    string**     _param;         // _val[replicate][generation]
    TStat_db* _pStat_db;
    
    
    StatRecBaseAll(TStat_db* db, StatRecBase* rec);
    StatRecBaseAll(TStat_db* db, string name, string title): _stat(0), _param(0){
        _pStat_db=db;
        _name=name;
        _title=title;
        _age=ALL;
        _ordering=PARAM;
    }

    ~StatRecBaseAll ( );
    
    
    void set_stat_pointer(double* array, unsigned int rep){
        assert(!_param); // either _param or _stat is used, but not both
        assert(_stat);
        assert(!_stat[rep]);
        _stat[rep] = array;
    }
    void set_param_pointer(string* array, unsigned int rep){
        assert(!_stat); // either _param or _stat is used, but not both
        assert(_param);
        assert(!_param[rep]);
        _param[rep] = array;
    }
    
    /**Sets the recorder attributes.
     * @param T the stat title
     * @param N the stat name (headers in the output file)
     * @param Order stat table ordering flag
     * @param AGE the age class for which the stat will be recorded
     * @param ARG the argument to pass to the S function
     **/
    //    void set     (string T, string N, st_order Order, age_t AGE, unsigned int ARG);
    
    ///@name Accessors
    ///@{
    //    void setName                      (string N)            {_name = N;}
    
    /**@param i the number of rows
       @param j the number of columns */
        st_order     getOrdering           ( )                        {return _ordering;}
        string       getTitle              ( )                        {return _title;}
        string       getName               ( )                        {return _name;}
        age_t        getAge                ( )                        {return _age;}
        unsigned int getArg                ( )                        {return _arg;}
    
        /**Returns the jth element of the ith row of the _val matrix.*/
    double getVal(unsigned int rep, unsigned int gen);
    string getParam(unsigned int rep, unsigned int gen);
    
    double   getMean(const unsigned int& index) {return (this->*getMean_func_ptr)(index);}
    double   (StatRecBaseAll::*getMean_func_ptr)(const unsigned int& index);
    double   getMedian(const unsigned int& index) {return (this->*getMedian_func_ptr)(index);}
    double   (StatRecBaseAll::*getMedian_func_ptr)(const unsigned int& index);
    double   getVar(const unsigned int& index) {return (this->*getVar_func_ptr)(index);}
    double   (StatRecBaseAll::*getVar_func_ptr)(const unsigned int& index);
    double   getMean_FLAT             (const unsigned int& gen_index);
    double   getMedian_FLAT           (const unsigned int& gen_index);
    double   getVar_FLAT              (const unsigned int& gen_index);
    double   getSum_FLAT              (const unsigned int& gen_index);
    
    //    unsigned int getIndex             (unsigned int i)           {return _index[i];}
    
    ///@}
    //    bool   (StatRecBase::*setVal_func_ptr)(int crnt_gen, int rpl_cntr, double value);
    //    bool   setVal_FLAT                (int crnt_gen, int rpl_cntr, double value);
    //    bool   setVal_GEN                 (int crnt_gen, int rpl_cntr, double value);
    //    bool   setVal_RPL                 (int crnt_gen, int rpl_cntr, double value);
    
};


////////////////////////////////////////////////////////////////////////////////
//#include "stat_rec_base.h"

class TStat_db {
private:
    unsigned int _nb_stats;//t
    unsigned int _nb_generation; // number of items to store per replicate
    unsigned int _nb_replicate;
    string _NaN;
    
    unsigned int _tot_occurrence;   // total number of values/generations to store
    unsigned int _occurrence;        // the current logtime
    
    unsigned int* _index;  // contains the stored generations _index[item]=generation; must be deleted
    
    // parameters and stats are stored in the same list: firs thte params then
    // the stats.
    list<StatRecBaseAll*> _all_stats;
    typedef list< StatRecBaseAll* >::iterator STAT_IT;
    
    int           _save_choice; // 0: all
    // 1: each rep (output immediatly)
    // 2: mean and var
    // 3: only mean
    // 4: only var
    // 5: only median
    // 6: nothing
    // 7: ABC: only mean file without header
    
    string 				_file_stats;  // file name for the stats (each generation and replicate)
    string				_file_mean;   // file name for the mean across replicates
    string				_file_median; // file name for the median across replicates
    string              _file_var;    // file name for the variance across replciates
    string              _file_legend; // file name for the legend
    
    
    
public:
    TStat_db(map<string, string> params, unsigned int rep);
    ~TStat_db();
    
    // index
    void set_index(unsigned int* i, unsigned int size) {
        _index=i;
        _tot_occurrence=size;
        _occurrence = _index[0];
    }
    void set_index(vector<unsigned int>& vec) {
        _tot_occurrence = (unsigned int)vec.size();
        assert(!_index);
        ARRAY::vector2array<unsigned int>(vec, _index, _tot_occurrence);
        _occurrence = _index[0];
    }
    
    unsigned int get_tot_occurrence(){return _tot_occurrence;}
    unsigned int get_occurrence(){return _occurrence;}
    void         set_occurrence(unsigned int i){_occurrence=i;}
    
    void            set_nb_replicate(unsigned int i){_nb_replicate = i;}
    unsigned int    get_nb_replicate(){return _nb_replicate;}
    
    StatRecBaseAll* add_statRecBase(StatRecBase* rec);
    void            add_paramRecBase(Param* pParam);
    list<StatRecBaseAll*>& get_all_stats(){return _all_stats;}
    
    void              print_headers    (ostream& FH, unsigned int order);
    void              print_value      (ostream& FH, unsigned int i, unsigned int j);
    void              print_mean       (ostream& FH, unsigned int i);
    void              print_median     (ostream& FH, unsigned int i);
    void              print_variance   (ostream& FH, unsigned int i);
    
    void printStat_each    ( );
    void printStat_mean    (bool header = true);
    void printStat_median  ( );
    void printStat_variance( );
    void printStat_legend  (unsigned int order);
    
    void printStatHeaders(ostream& FH,unsigned int order);
    void printStatLegend(ostream& FH,unsigned int order);
    void printStatValue(ostream& FH, unsigned int i, unsigned int j);
    void printStatMean(ostream& FH, unsigned int i);
    void printStatMedian(ostream& FH, unsigned int i);
    void printStatVariance(ostream& FH, unsigned int i);
    unsigned int getStatRecIndex(unsigned int i);
    
    string get_file_name_legend( )const {return _file_legend;}
    string get_file_name_stats( ) const {return _file_stats;}
    string get_file_name_mean( )  const {return _file_mean;}
    string get_file_name_median( )const {return _file_median;}
    string get_file_name_var( )   const {return _file_var;}
    int    get_save_choice( )     const {return _save_choice;}
    
    void set_file_name_legend(string s){_file_legend = s;}
    void set_file_name_stats(string s) {_file_stats = s;}
    void set_file_name_mean(string s)  {_file_mean = s;}
    void set_file_name_median(string s){_file_median = s;}
    void set_file_name_var(string s)   {_file_var = s;}
    void set_save_choice(const int& i) {_save_choice = i;}
    
    void          set_NaN(string s)  {_NaN = s;}
    string        get_NaN()          {return _NaN;}
    

    void FHwrite();
};


#endif /* defined(__quantiNemo2__stat_db__) */
