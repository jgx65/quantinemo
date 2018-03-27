/** @file tstat_db.cpp
 *
 *   Copyright (C) 2008 Samuel Neuenschwander <samuel.neuenschwander@unil.ch>
 *
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

//

#include "tstat_db.h"
#include "stat_rec_base.h"
#include "param.h"

#include <fstream>
using namespace std;


// ----------------------------------------------------------------------------------------
// StatRecBaseAll
// ----------------------------------------------------------------------------------------
StatRecBaseAll::StatRecBaseAll(TStat_db* db, StatRecBase* rec): _stat(0), _param(0)
{
    _pStat_db = db;
    _name = rec->getName();
    _title = rec->getTitle();
    _arg = rec->getArg();
    _ordering = rec->getOrdering();
    
    assert(!_stat);
    assert(!_param);
    
    ARRAY::create_1D(_stat, _pStat_db->get_nb_replicate(), (double*)NULL);
    
    switch(_ordering){
        case PARAM:  // parameter arguments
            break;
        case FLAT:  // values are stored for each generation and replicate
            getMean_func_ptr    = &StatRecBaseAll::getMean_FLAT;    // across replciates
            getMedian_func_ptr  = &StatRecBaseAll::getMedian_FLAT;  // across replciates
            getVar_func_ptr     = &StatRecBaseAll::getVar_FLAT;     // across replciates
            break;
        case GEN:   // values are stored for each generation: replicate values are added
            getMean_func_ptr    = &StatRecBaseAll::getSum_FLAT;
            getMedian_func_ptr  = &StatRecBaseAll::getSum_FLAT;
            getVar_func_ptr     = &StatRecBaseAll::getSum_FLAT;
            break;
        case RPL:   // values are stored for each replicate: generation values are added
            error("Not anymore used!\n");
//            getMean_func_ptr    = &StatRecBaseAll::get_RPL;
//            getMedian_func_ptr  = &StatRecBaseAll::get_RPL;
//            getVar_func_ptr     = &StatRecBaseAll::get_RPL;
//            break;
    }
}

// ----------------------------------------------------------------------------------------
// add_statRecBase
// ----------------------------------------------------------------------------------------
StatRecBaseAll::~StatRecBaseAll( )
{
    ARRAY::delete_2D(_stat, _pStat_db->get_nb_replicate());
    ARRAY::delete_2D(_param, _pStat_db->get_nb_replicate());
}

// ----------------------------------------------------------------------------------------
// getVal
// ----------------------------------------------------------------------------------------
/**Returns the jth element of the ith row of the _val matrix.*/
double
StatRecBaseAll::getVal(unsigned int rep, unsigned int gen)
{
    assert(gen<_pStat_db->get_tot_occurrence());
    assert(rep<_pStat_db->get_nb_replicate());
    return _stat[rep][gen];
}

// ----------------------------------------------------------------------------------------
// getVal
// ----------------------------------------------------------------------------------------
/**Returns the jth element of the ith row of the _val matrix.*/
string
StatRecBaseAll::getParam(unsigned int rep, unsigned int gen)
{
    assert(gen<_pStat_db->get_tot_occurrence());
    assert(rep<_pStat_db->get_nb_replicate());
    return _param[rep][gen];
}

// ----------------------------------------------------------------------------------------
// StatRecorder::get
// ----------------------------------------------------------------------------------------
/** values are present for each generation and replicate
 * stats across replicates, i.e. per generation
 * gen_ind: generation index and NOT the generation itself!
 */
double
StatRecBaseAll::getMean_FLAT(const unsigned int& gen_index)
{
    return ARRAY::mean2D(_stat, _pStat_db->get_nb_replicate(), gen_index);
}

double
StatRecBaseAll::getMedian_FLAT(const unsigned int& gen_index)
{
    return ARRAY::median2D(_stat, _pStat_db->get_nb_replicate(), gen_index);
}

double
StatRecBaseAll::getVar_FLAT(const unsigned int& gen_index)
{
    return ARRAY::var2D(_stat, _pStat_db->get_nb_replicate(), gen_index);
}

double
StatRecBaseAll::getSum_FLAT(const unsigned int& gen_index)
{
    return ARRAY::sum2D(_stat, _pStat_db->get_nb_replicate(), gen_index);
}

// ----------------------------------------------------------------------------------------
// TStat_db
// ----------------------------------------------------------------------------------------
TStat_db::TStat_db(map<string, string> params, unsigned int rep): _nb_stats(0), _nb_generation(0), _nb_replicate(rep), _index(0)
{
    map<string, string>::iterator iter = params.find("stat_NaN");
    if(iter==params.end()) set_NaN(my_NANstr);
    else                   set_NaN(iter->second);
}

// ----------------------------------------------------------------------------------------
// ~TStat_db
// ----------------------------------------------------------------------------------------
TStat_db::~TStat_db()
{
    if(_index) delete[] _index;
    
    for(STAT_IT cur=_all_stats.begin(), end=_all_stats.end(); cur!=end; ++cur){
        delete *cur;
    }
    _all_stats.clear();
    
}

// ----------------------------------------------------------------------------------------
// add_statRecBase
// ----------------------------------------------------------------------------------------
StatRecBaseAll*
TStat_db::add_statRecBase(StatRecBase* rec)
{
    StatRecBaseAll* pp = new StatRecBaseAll(this, rec);
    assert(pp->_stat);
    assert(!pp->_param);
    _all_stats.push_back(pp);
    return pp;
}

// ----------------------------------------------------------------------------------------
// add_statRecBase
// ----------------------------------------------------------------------------------------
void
TStat_db::add_paramRecBase(Param* pParam)
{
    StatRecBaseAll* pp = new StatRecBaseAll(this,pParam->get_name(),
        "Output of parameter agrument (only in detailed output)");
    assert(!pp->_param);
    assert(!pp->_stat);
    ARRAY::create_1D(pp->_param, _nb_replicate, (string*)NULL);
    _all_stats.push_back(pp);
}

// ----------------------------------------------------------------------------------------
// print_headers
// ----------------------------------------------------------------------------------------
void
TStat_db::print_headers(ostream& FH, unsigned int order)
{
    for(list< StatRecBaseAll* >::iterator IT=_all_stats.begin(); IT!=_all_stats.end(); ++IT) {
        if((*IT)->getOrdering() & order) {
            FH<<"\t";
            FH.width(12);
            FH<<(*IT)->getName();
        }
    }
}


// ----------------------------------------------------------------------------------------
// FHwrite
// ----------------------------------------------------------------------------------------
void
TStat_db::FHwrite()
{
    switch(get_save_choice()) {
        case 0: // all
            printStat_each();
            printStat_mean();
            printStat_variance();
            printStat_legend((GEN | FLAT | PARAM));
            break;
        case 1: // each replicate (no db)
            printStat_legend((GEN | FLAT | PARAM));
            break;
        case 2: // mean and var
            printStat_mean();
            printStat_variance();
            printStat_legend((GEN | FLAT));
            break;
        case 3: // only mean
            printStat_mean();
            printStat_legend((GEN | FLAT));
            break;
        case 4: // only var
            printStat_variance();
            printStat_legend((GEN | FLAT));
            break;
        case 5: // only median
            printStat_median();
            printStat_legend((GEN | FLAT));
            break;
        case 6:
            break; // no ouput
        case 7: // ABC: no header
            printStat_mean(false);
            break;
    }
}

// ----------------------------------------------------------------------------------------
// printStat_legend
// ----------------------------------------------------------------------------------------
/** prints the legend of the statistics (only once at generation 0) */
void
TStat_db::printStat_legend(unsigned int order)
{
#ifdef _DEBUG
    message("  TStat_db::PrintStat_legend (%s)\n", get_file_name_legend().c_str());
#endif
    
    ofstream FH(get_file_name_legend().c_str(), ios::out);
    
    if (!FH)
        error("Could not open stat legend file '%s'!\n", get_file_name_legend().c_str());
    
    // print the heading
    FH << "Legend of the used statistics\n" <<
    "*****************************\n\n";
    FH.width(20);
    FH.setf(ios::left, ios::adjustfield);
    FH << "Statistic" << "  Legend\n";
    FH << "--------------------------------------------------------------\n";
    
    // print the stat names:
    printStatLegend(FH, order);
    
    FH.close();
}

// ----------------------------------------------------------------------------------------
// printStat_each
// ----------------------------------------------------------------------------------------
/** prints the stats for each replicate separately */
void
TStat_db::printStat_each()
{
#ifdef _DEBUG
    message("  TStat_db::FHwrite (%s)\n", get_file_name_stats().c_str());
#endif
    
    ofstream FH(get_file_name_stats().c_str(), ios::out);
    if (!FH)
        error("Could not open stat output file \"%s\"!\n", get_file_name_stats().c_str());
    
    printStatHeaders(FH, (FLAT | PARAM)); // print the stat names
    
    unsigned int i, j;
    unsigned int ncols = get_nb_replicate();
    unsigned int nrows = get_tot_occurrence();
    
    for (i = 0; i < ncols; i++) {
        for (j = 0; j < nrows; j++) {
            // val is ordered rows x cols, row = generation, col = replicate
            // we write: for 1st rpl all the stat for all the gen recorded, then we move to the next rpl.
            FH.width(12);
            FH.setf(ios::left, ios::adjustfield);
            FH << setprecision(4);
            FH << (i + 1);
            FH << "\t";
            FH.width(12);
            FH << _index[j];
            
            printStatValue(FH, i, j);
            
            FH << "\n";
        }
    }
    FH.close();
}

// ----------------------------------------------------------------------------------------
// PrintStat_mean
// ----------------------------------------------------------------------------------------
/** merges the stats across replicates */
void
TStat_db::printStat_mean(bool header)
{
#ifdef _DEBUG
    message("  TStat_db::PrintStat_mean (%s)\n", get_file_name_mean().c_str());
#endif
    
    ofstream FH(get_file_name_mean().c_str(), ios::out);
    if (!FH) error("Could not open stat mean file '%s'!\n", get_file_name_mean().c_str());
    
    if (header) printStatHeaders(FH, (GEN | FLAT));
    // print the stat names
    
    unsigned int nrows = get_tot_occurrence();
    
    for (unsigned int i = 0; i < nrows; i++) {
        FH.width(12);
        FH.setf(ios::left, ios::adjustfield);
        FH << setprecision(4);
        FH << _index[i];
        
        printStatMean(FH, i);
        
        FH << "\n";
    }
    FH.close();
}

// ----------------------------------------------------------------------------------------
// PrintStat_median
// ----------------------------------------------------------------------------------------
/** merges the stats across replicates */
void
TStat_db::printStat_median()
{
#ifdef _DEBUG
    message("  TStat_db::PrintStat_median (%s)\n", get_file_name_median().c_str());
#endif
    
    ofstream FH(get_file_name_mean().c_str(), ios::out);
    if (!FH) error("Could not open stat median file '%s'!\n", get_file_name_median().c_str());
    
    // print the first row with stats headers:
    FH.width(12);
    FH.setf(ios::left, ios::adjustfield);
    FH << "generation";
    
    // print the stat names:
    printStatHeaders(FH, (GEN | FLAT));
    
    FH << "\n";
    
    // next rows:
    unsigned int nrows = get_tot_occurrence();
    
    for (unsigned int i = 0; i < nrows; i++) {
        FH.width(12);
        FH.setf(ios::left, ios::adjustfield);
        FH << setprecision(4);
        FH << _index[i];
        
        printStatMedian(FH, i);
        
        FH << "\n";
    }
    FH.close();
}

// ----------------------------------------------------------------------------------------
// PrintStat_var
// ----------------------------------------------------------------------------------------
/** merges the stats across replicates */
void
TStat_db::printStat_variance()
{
#ifdef _DEBUG
    message("  TStat_db::PrintStat_var (%s)\n", get_file_name_var().c_str());
#endif
    
    ofstream FH(get_file_name_var().c_str(), ios::out);
    if (!FH) error("Could not open stat var file '%s'!\n", get_file_name_var().c_str());
    
    printStatHeaders(FH, (GEN | FLAT));
    // print the stat names
    
    unsigned int nrows = get_tot_occurrence();
    
    // val is ordered rows x cols, row = generation, col = replicate
    // we take the variance over rows
    
    for (unsigned int i = 0; i < nrows; i++) {
        FH.width(12);
        FH.setf(ios::left, ios::adjustfield);
        FH << setprecision(4);
        FH << _index[i];
        
        printStatVariance(FH, i);
        
        FH << "\n";
    }
    FH.close();
}

//-----------------------------------------------------------------------------------------
// print_headers
// ----------------------------------------------------------------------------------------
/** prints the heading  for the stat files */
void
TStat_db::printStatHeaders(ostream& FH,unsigned int order)
{
    //print the first row with stat headers:
    FH.width(12);
    FH.setf(ios::left,ios::adjustfield);
    if(!(order & GEN)){
        FH<<"replicate"<<"\t";
        FH.width(12);
    }
    FH<<"generation";
    
    // print the name of each stat
    print_headers(FH,order);
    
    FH<<"\n";        //new row
}

//-----------------------------------------------------------------------------------------
// print_legend
// ----------------------------------------------------------------------------------------
void
TStat_db::printStatLegend(ostream& FH,unsigned int order)
{
    for(STAT_IT IT=_all_stats.begin(); IT!=_all_stats.end(); ++IT) {
        if((*IT)->getOrdering() & order) {
            FH.width(20);
            FH.setf(ios::left,ios::adjustfield);
            FH<<(*IT)->getName();
            FH<<": " << (*IT)->getTitle() << "\n";
        }
    }
}

//-----------------------------------------------------------------------------------------
// printStatValue
// ----------------------------------------------------------------------------------------
void
TStat_db::printStatValue(ostream& FH, unsigned int i, unsigned int j)
{
    double value;
    for(STAT_IT IT=_all_stats.begin(); IT!=_all_stats.end(); ++IT) {
        if((*IT)->getOrdering() == FLAT) {
            FH<<"\t";
            FH.width(12);
            value = (*IT)->getVal(i,j);
            if(value==my_NAN) FH << get_NaN();
            else              FH << value;
        }
        else if((*IT)->getOrdering() == PARAM) {  // currently used for the seed
            FH<<"\t";
            FH.width(12);
            string tt = (*IT)->getParam(i,j);
            FH << (*IT)->getParam(i,j);
        }
    }
}

//-----------------------------------------------------------------------------------------
// printStatMean
// ----------------------------------------------------------------------------------------
void
TStat_db::printStatMean(ostream& FH, unsigned int i)
{
    double value;
    for(STAT_IT IT = _all_stats.begin(); IT != _all_stats.end(); ++IT) {
        if((*IT)->getOrdering() == PARAM) continue;
        
        FH<<"\t";
        FH.width(12);
        value = (*IT)->getMean(i);
        if(value==my_NAN) FH << get_NaN();
        else              FH << value;
    }
}

//-----------------------------------------------------------------------------------------
// printStatMedian
// ----------------------------------------------------------------------------------------
void
TStat_db::printStatMedian(ostream& FH, unsigned int i)
{
    double value;
    for(STAT_IT IT = _all_stats.begin(); IT != _all_stats.end(); ++IT) {
        if((*IT)->getOrdering() == PARAM) continue;

        FH<<"\t";
        FH.width(12);
        value = (*IT)->getMedian(i);
        if(value==my_NAN) FH << get_NaN();
        else              FH << value;
    }
}

//-----------------------------------------------------------------------------------------
// printStatVariance
// ----------------------------------------------------------------------------------------
void
TStat_db::printStatVariance(ostream& FH, unsigned int i)
{
    double value;
    for(STAT_IT IT = _all_stats.begin(); IT != _all_stats.end(); ++IT) {
        if((*IT)->getOrdering() == PARAM) continue;
 
        FH<<"\t";
        FH.width(12);
        value = (*IT)->getVar(i);
        if(value==my_NAN) FH << get_NaN();
        else              FH << value;
    }
}


