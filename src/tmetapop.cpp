/** @file metapop.cpp
 *
 *   Copyright (C) 2006 Frederic Guillaume    <guillaum@zoology.ubc.ca>
 *   Copyright (C) 2008 Samuel Neuenschwander <samuel.neuenschwander@unil.ch>
 *   Copyright (C) 2018 Frederic Michaud <frederic.michaud@unil.ch>

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

#include "tmetapop.h"
#include <algorithm>

#include "ttquanti.h"
//#include "ttneutral.h"
#include "lce_breed.h"
#include "lce_disperse.h"
#include "lce_coalescence_base.h"
//#include "lce_misc.h"
#include "tstat_db.h"

#include "tsimulation.h"
#include "treplicate.h"
#include "tselection.h"

#include "stathandler.cpp"
#include <time.h>

using namespace std;


unsigned int TMetapop::getReplicates   ( ) {assert(_pSimulation); return _pSimulation->get_replicates();}
TSimulation* TMetapop::get_pSimulation(){return _pSimulation;}
TStat_db*    TMetapop::get_pStat_db(){assert(_pSimulation); return _pSimulation->stats;}
void         TMetapop::set_pStat_db(TStat_db* p){_pSimulation->stats = p;}




// ----------------------------------------------------------------------------------------
// constructor
// ----------------------------------------------------------------------------------------
TMetapop::TMetapop(TReplicate* p, unsigned int rep) :_statHandler(), _tot_sampled_patches(0), _patchNbr(0),
_generations(0), _currentGeneration(0),
_currentAge(NONE),  _sexInitRatio(0.5), _service(0), _pSelection(0),
_total_carrying_capacity(my_NAN), _pReplicate(p),
_density_threshold(0), _density_threshold_param(0),  _sampleAllOrNothingCur(0),
_current_replicate(rep), _current_index_stat_db(my_NAN), func_ptr_new_fullPatch(0),
func_ptr_new_emptyPatch(0), func_ptr_add_tempPatch(0),_last_nbSamplePatch(0),_current_nbSamplePatch(0)
{
    _pReplicate  = p;
    _pSimulation = p->get_pSimulation();
    _pSimulation->get_paramset()._popPtr = this;
    if(_pSimulation->get_paramSetKeys()) _pSimulation->get_paramSetKeys()->_popPtr = this;
    
    
    
    
    _isSampled[0] = 0;  // generation
    _isSampled[1] = 0;  // replicate
    _isSampled[2] = 0;  // age
    
#ifdef _DEBUG
    message(" Metapop::Metapop()\n");
#endif
    
    _statHandler._popPtr=this;
    if(_pSimulation->stats) _current_statRecBase  = _pSimulation->stats->get_all_stats().begin();
    
    init_paramset();
    
}

// ----------------------------------------------------------------------------------------
// init_paramset
// ----------------------------------------------------------------------------------------
void
TMetapop::init_paramset()
{
    set_paramset("metapop", "metapopulation", true, this);
    
    _paramSet->add_param("generations",INT2,true,0,my_NAN,"0",false,
                         "Number of generations to run.",0);
    
    _paramSet->add_param("coalescence",INT2,false,0,4,"0",false,
                         "Mode of simulation:\n" \
                         "  0: individual based (forward in time)\n" \
                         "  1: population based (coalescence)\n" ,3);
    
    
    // carrying capacities
    _paramSet->add_param("patch_number",INT2,false,0,my_NAN,"1",false,
                         "Number of patches of the metapopulation.",0);
    
    _paramSet->add_param("patch_capacity",INT_MAT,false,0,my_NAN,"",true,
                         "Carrying capacity of a patch.",0);
    
    _paramSet->add_param("patch_capacity_fem",INT_MAT,false,0,my_NAN,"",true,
                         "Female carrying capacity of a patch.",2);

    _paramSet->add_param("patch_capacity_mal",INT_MAT,false,0, my_NAN,"",true,
                         "Male carrying capacity of a patch.",2);
    
    
    // population initial sizes
    _paramSet->add_param("patch_ini_size",INT_MAT,false,0,my_NAN,"",false,
                         "Initial population size of a patch.",1);
    
    _paramSet->add_param("patch_ini_size_fem",INT_MAT,false,0,my_NAN,"",false,
                         "Initial population size of females of a patch.",2);
    
    _paramSet->add_param("patch_ini_size_mal",INT_MAT,false,0,my_NAN,"",false,
                         "Initial population size of males of a patch.",2);
    
    
    // environmental influence
    _paramSet->add_param("patch_ve_mean",DBL_MAT,false,my_NAN,my_NAN,"0", true,
                         "Mean environmental variance of a patch.",5);
    
    _paramSet->add_param("patch_ve_mean_fem",DBL_MAT,false,my_NAN,my_NAN,"0",true,
                         "Mean environmental variance of females of a patch.",5);
    
    _paramSet->add_param("patch_ve_mean_mal",DBL_MAT,false,my_NAN,my_NAN,"0",true,
                         "Mean environmental variance of males of a patch.",5);
    
    _paramSet->add_param("patch_ve_var",DBL_MAT,false,0,my_NAN,"0",true,
                         "Variance of the mean environmental variance of a patch.",5);
    
    _paramSet->add_param("patch_ve_var_fem",DBL_MAT,false,0,my_NAN,"0",true,
                         "Variance of the mean environmental variance of females of a patch.",5);
    
    _paramSet->add_param("patch_ve_var_mal",DBL_MAT,false,0,my_NAN,"0",true,
                         "Variance of the mean environmental variance of males of a patch.",5);
    
    // stabilizing selection
    _paramSet->add_param("patch_stab_sel_optima",DBL_MAT,false,my_NAN,my_NAN,"0",true,
                         "Stabilizing selection: the fitness optimum",0);
    
    _paramSet->add_param("patch_stab_sel_optima_fem",DBL_MAT,false,my_NAN,my_NAN,"0",true,
                         "Stabilizing selection: the fitness optimum for females",2);
    
    _paramSet->add_param("patch_stab_sel_optima_mal",DBL_MAT,false,my_NAN,my_NAN,"0",true,
                         "Stabilizing selection: the fitness optimum for males",2);
    
    _paramSet->add_param("patch_stab_sel_optima_var",DBL_MAT,false,0,my_NAN,"0",true,
                         "Stabilizing selection: the variance of the fitness optimum over time",5);
    
    _paramSet->add_param("patch_stab_sel_intensity",DBL_MAT,false,0,my_NAN,"1",true,
                         "Stabilizing selection: the selection intensity.\n" \
                         "Caution: the higher the value is, the weaker the selection pressure is.",0);
    
    _paramSet->add_param("patch_stab_sel_intensity_fem",DBL_MAT,false,0,my_NAN,"1",true,
                         "Stabilizing selection: the selection intensity for females.\n" \
                         "Caution: the higher the value is, the weaker the selection pressure is.",2);
    
    _paramSet->add_param("patch_stab_sel_intensity_mal",DBL_MAT,false,0,my_NAN,"1",true,
                         "Stabilizing selection: the selection intensity for males.\n" \
                         "Caution: the higher the value is, the weaker the selection pressure is.",2);
    
    _paramSet->add_param("patch_stab_sel_intensity_var",DBL_MAT,false,0,my_NAN,"0",true,
                         "Stabilizing selection: the variance of the fitness optimum over time",5);
    
    
    // directional selection
    _paramSet->add_param("patch_dir_sel_min",				DBL_MAT,false,0,		 1,"0",true,
                         "Directional selection: lower asymptote of selection curve.",0);
    
    _paramSet->add_param("patch_dir_sel_min_fem",			DBL_MAT,false,0,		 1,"0",true,
                         "Directional selection: lower asymptote of selection curve for females.",2);
    
    _paramSet->add_param("patch_dir_sel_min_mal",			DBL_MAT,false,0,		 1,"0",true,
                         "Directional selection: lower asymptote of selection curve for males.",2);
    
    _paramSet->add_param("patch_dir_sel_min_var",			DBL_MAT,false,0,		 my_NAN,"0",true,
                         "Directional selection: variance of the lower asymptote over time.",5);
    
    _paramSet->add_param("patch_dir_sel_max",				DBL_MAT,false,0,		 1,"1",true,
                         "Directional selection: upper asymptote of selection curve.",0);
    
    _paramSet->add_param("patch_dir_sel_max_fem",			DBL_MAT,false,0,		 1,"1",true,
                         "Directional selection: upper asymptote of selection curve for females.",2);
    
    _paramSet->add_param("patch_dir_sel_max_mal",			DBL_MAT,false,0,		 1,"1",true,
                         "Directional selection: lower asymptote of selection curve for males.",2);
    
    _paramSet->add_param("patch_dir_sel_max_var",			DBL_MAT,false,0,		 my_NAN,"0",true,
                         "Directional selection: variance of the upper asymptote over time.",5);
    
    _paramSet->add_param("patch_dir_sel_max_growth",        DBL_MAT,false,my_NAN,my_NAN,"0",true,
                         "Directional selection: phenotype with maximal growth.",0);
    
    _paramSet->add_param("patch_dir_sel_max_growth_fem",    DBL_MAT,false,my_NAN,my_NAN,"0",true,
                         "Directional selection: phenotype with maximal growth for females.",2);
    
    _paramSet->add_param("patch_dir_sel_max_growth_mal",    DBL_MAT,false,my_NAN,my_NAN,"0",true,
                         "Directional selection: phenotype with maximal growth for males.",2);
    
    _paramSet->add_param("patch_dir_sel_max_growth_var",    DBL_MAT,false,0,	 my_NAN,"0",true,
                         "Directional selection: variance of the phenotype with maximal growth over time.",5);
    
    _paramSet->add_param("patch_dir_sel_growth_rate",       DBL_MAT,false,my_NAN,my_NAN,"1",true,
                         "Directional selection: growth rate (when positive, larger values have a " \
                         "higher fitness; when negative, lower phenotypes have a higher fitness.",0);
    
    _paramSet->add_param("patch_dir_sel_growth_rate_fem",   DBL_MAT,false,my_NAN,my_NAN,"1",true,
                         "Directional selection: growth rate fir females (when positive, larger values have a " \
                         "higher fitness; when negative, lower phenotypes have a higher fitness.",2);
    
    _paramSet->add_param("patch_dir_sel_growth_rate_mal",   DBL_MAT,false,my_NAN,my_NAN,"1",true,
                         "Directional selection: growth rate for males (when positive, larger values have a " \
                         "higher fitness; when negative, lower phenotypes have a higher fitness.",2);
    
    _paramSet->add_param("patch_dir_sel_growth_rate_var",   DBL_MAT,false,0,	 my_NAN,"0",true,
                         "Directional selection: variance of the growth rate over time.",5);
    
    _paramSet->add_param("patch_dir_sel_symmetry",          DBL_MAT,false,0,	1,		"1",true,
                         "Directional selection: symmetry of the slope.",0);
    
    _paramSet->add_param("patch_dir_sel_symmetry_fem",      DBL_MAT,false,0,	1,		"1",true,
                         "Directional selection: symmetry of the slope for females.",2);
    
    _paramSet->add_param("patch_dir_sel_symmetry_mal",      DBL_MAT,false,0,	1,		"1",true,
                         "Directional selection: symmetry of the slope for males.",2);
    
    _paramSet->add_param("patch_dir_sel_symmetry_var",      DBL_MAT,false,0,	my_NAN,"0",true,
                         "Directional selection: variance of the symmetry over time.",5);
    
    
    // fitness landscape
    _paramSet->add_param("patch_fitness_landscape",         DBL_MAT,false,my_NAN,my_NAN,"1",true,
                         "Fitness landscape: Specifies the fitness values of the phenotype-fitness " \
                         "landscape (to be used with parameter 'patch_phenotype_landscape').",0);
    
    _paramSet->add_param("patch_fitness_landscape_fem",     DBL_MAT,false,my_NAN,my_NAN,"1",true,
                         "Fitness landscape: Specifies the fitness values for females of the phenotype-fitness " \
                         "landscape (to be used with parameter 'patch_phenotype_landscape').",2);
    
    _paramSet->add_param("patch_fitness_landscape_mal",     DBL_MAT,false,my_NAN,my_NAN,"1",true,
                         "Fitness landscape: Specifies the fitness values for males of the phenotype-fitness " \
                         "landscape (to be used with parameter 'patch_phenotype_landscape').",2);
    
    _paramSet->add_param("patch_phenotype_landscape",       DBL_MAT,false,my_NAN,my_NAN,"0",true,
                         "Fitness landscape: Specifies the phenotype values of the phenotype-fitness " \
                         "landscape (to be used with parameter 'patch_fitness_landscape').",0);
    
    _paramSet->add_param("patch_phenotype_landscape_fem",   DBL_MAT,false,my_NAN,my_NAN,"0",true,
                         "Fitness landscape: Specifies the phenotype values for females of the phenotype-fitness " \
                         "landscape (to be used with parameter 'patch_fitness_landscape_fem').",2);
    
    _paramSet->add_param("patch_phenotype_landscape_mal",   DBL_MAT,false,my_NAN,my_NAN,"0",true,
                         "Fitness landscape: Specifies the phenotype values for males of the phenotype-fitness " \
                         "landscape (to be used with parameter 'patch_fitness_landscape_mal').",2);
    
    
    // selection coefficient
    _paramSet->add_param("patch_coef_sel_AA",               DBL_MAT,false,my_NAN,my_NAN,"1",true,
                         "Selection coefficient: specifies the fitness for wild type AA.",0);
    
    _paramSet->add_param("patch_coef_sel_AA_fem",           DBL_MAT,false,my_NAN,my_NAN,"1",true,
                         "Selection coefficient: specifies the fitness for wild type AA of females.",2);
    
    _paramSet->add_param("patch_coef_sel_AA_mal",           DBL_MAT,false,my_NAN,my_NAN,"1",true,
                         "Selection coefficient: specifies the fitness for wild type AA of males.",2);
    
    _paramSet->add_param("patch_coef_sel_AA_var",           DBL_MAT,false,0,	 my_NAN,"0",true,
                         "Selection coefficient: Variance over time of the fitness for wild type AA.",3);
    
    _paramSet->add_param("patch_coef_sel",                  DBL_MAT,false,my_NAN,my_NAN,"0",true,
                         "Selection coefficient: The selection coefficient 's'.",0);
    
    _paramSet->add_param("patch_coef_sel_fem",              DBL_MAT,false,my_NAN,my_NAN,"0",true,
                         "Selection coefficient: The selection coefficient of females 's'.",2);
    
    _paramSet->add_param("patch_coef_sel_mal",              DBL_MAT,false,my_NAN,my_NAN,"0",true,
                         "Selection coefficient: The selection coefficient of males 's'.",2);
    
    _paramSet->add_param("patch_coef_sel_var",              DBL_MAT,false,0,		 my_NAN,"0",true,
                         "Selection coefficient: Variance over time of the selection coefficient 's'.",3);
    
    
    // how should the seleciton be defined? (0: using patch_x; 1: using quanti_x)
    _paramSet->add_param("selection_pressure_definition",   INT2,false,0,1,"1", false,
                         "Should selection be defined as:\n" \
                         "  0: patch (at the patch level 'patch_*')\n" \
                         "  1: quanti (at the quantitative trait level 'quanti_*)",0);
    
    // sample specifications
    _paramSet->add_param("sample_all_or_nothing",           INT_MAT,false,0,    1,"3",false,
                         "How to sample individuals, when not enough individuals are present:\n" \
                         "  0: when possible (output and compute stats with what is present)\n" \
                         "  1: all or nothing (if not enough individuals can be sampled, " \
                         "then don't output and compute stats)",3);
    
    _paramSet->add_param("sampled_patches",					INT_MAT,false,0,    my_NAN,"0", false,
                         "What to sample:\n" \
                         "  0: all patches are entirely sampled\n" \
                         "  integer: draw randomly 'n' patches to sample entirely\n" \
                         "  matrix: specifies explicitly the patches to sample (matrix contains the patch IDs",3);
    
    _paramSet->add_param("patch_sample_size",               DBL_MAT,false,0,    my_NAN,"NaN",true,
                         "Specifies the number of individuals to sample per patch:\n" \
                         "  integer: absolute number of individuals\n" \
                         "  decimal: relative number of individuals",2);
    
    _paramSet->add_param("patch_sample_size_fem",           DBL_MAT,false,0,    my_NAN,"NaN",true,
                         "Specifies the number of females to sample per patch:\n" \
                         "  integer: absolute number of females\n" \
                         "  decimal: relative number of females",3);
    
    _paramSet->add_param("patch_sample_size_mal",           DBL_MAT,false,3,    my_NAN,"NaN",true,
                         "Specifies the number of males to sample per patch:\n" \
                         "  integer: absolute number of males\n" \
                         "  decimal: relative number of males",5);
    
    
    _paramSet->add_param("selection_position",				INT2,   false,0,    5,		"0",false,
                         "Where does selection acts:\n" \
                         "  0: reproductive success (the fitter an individual is the more offspring it gets)\n" \
                         "  1: Not functionnal " \
                         "is the same as for 'reproductive success', but reproduction is neutral " \
                         "and the fitness of the offspring is considered)\n" \
                         "  2: juvenile survival (offspring are created and only latter selection acts)\n" \
                         "  3: adult survival\n" \
                         "  4: no selection",1);
    
    _paramSet->add_param("selection_level",                 INT2,	false,0,    3,		"0",false,
                         "How does selection acts:\n" \
                         "  0: soft (fitness is relative to the other individuals  of the patch)\n" \
                         "  1: metapop (fitness is relative to the other individuals  of the metapopulation)\n" \
                         "  2: hard (fitness is absolute)\n TODO " ,0);
    
    _paramSet->add_param("patch_mean_fitness",              INT2,	false,0,	2,		"0",false,
                         "Mean fitness of a patch id defined as:\n" \
                         "  0: the average fitness across all individuals  (females and males)\n" \
                         "  1: the average fitness across females only\n" \
                         "  2: the average fitness across males only",4);
    
    _paramSet->add_param("temporal_change_following_density",STR_MAT,false,my_NAN,my_NAN,"0",false,
                         "Still under development",5);
    
}

// ----------------------------------------------------------------------------------------
// destructor
// ----------------------------------------------------------------------------------------
TMetapop::~TMetapop()
{
    clear();
#ifdef _DEBUG
    message(" Metapop::~Metapop()\n");
#endif
}

// ----------------------------------------------------------------------------------------
// setSampleSizes
// ----------------------------------------------------------------------------------------
/** set the sample sizes
 * every time the sampling schema changes time this function is called
 */
bool TMetapop::setSampleSizes(bool set_sampleID)
{
    if(setPatchParam<double>("patch_sample_size", // parameter name (without "_fem" or "_mal")
                             &TPatch::set_N_sample,           // set a sex specific value
                             &TPatch::set_N_sample,           // set a general value
                             &TPatch::set_N_sample,           // function to set all three params at once
                             &TPatch::get_N_sample,           // get a sex specific value
                             &TPatch::get_N_sample))          // get a general value
    {
        if(set_sampleID){
            for(unsigned int i = 0; i < _patchNbr; ++i){  // no parameter set: default initialization
                if(_vPatch[i]->get_N_sample()) _vPatch[i]->set_sampleID(SAMPLED);   // patch is sampled
                else                           _vPatch[i]->set_sampleID(my_NAN);    // patch is not sampled
            }
        }
        return true;                                // it was set by the input
    }
    else{                                         // all is sampled
        if(set_sampleID){
            for(unsigned int i = 0; i < _patchNbr; ++i){  // no parameter set: default initialization
                _vPatch[i]->set_N_sample(my_NAN, my_NAN, my_NAN);
                _vPatch[i]->set_sampleID(SAMPLED);
            }                                                                  // patch is sampled
        }
        return false;                                // it was not set by the input
    }
}

// ----------------------------------------------------------------------------------------
// setCarryingCapacities
// ----------------------------------------------------------------------------------------
/** set the carrying capacities
 */
void TMetapop::setCarryingCapacities()
{
    if(!setPatchParam<unsigned int>("patch_capacity", // parameter name (without "_fem" or "_mal")
                                    &TPatch::set_K,                      // set a sex specific value
                                    &TPatch::set_K,                      // set a general value
                                    &TPatch::set_K,                      // function to set all three params at once
                                    &TPatch::get_K,                      // get a sex specific value
                                    &TPatch::get_K))                     // get a general value
    {
        error("Carrying capacity: the carrying capacity has to be set!\n");
    }
}



// ----------------------------------------------------------------------------------------
// setInitPopulationSizes
// ----------------------------------------------------------------------------------------
/** set the initial populations
 */
void TMetapop::setInitPopulationSizes()
{
    if(!setPatchParam<unsigned int>("patch_ini_size", // parameter name (without "_fem" or "_mal")
                                    &TPatch::set_N_ini,                  // set a sex specific value
                                    &TPatch::set_N_ini,                  // set a general value
                                    &TPatch::set_N_ini,                  // function to set all three params at once
                                    &TPatch::get_N_ini,                  // get a sex specific value
                                    &TPatch::get_N_ini))                 // get a general value
    {
        for(unsigned int i = 0; i < _patchNbr; ++i){    // no parameter set: default initialization
            _vPatch[i]->set_PopSizes_ini_carrying_capacity();
        }
    }
}

// ----------------------------------------------------------------------------------------
// set_sampled_patches
// ----------------------------------------------------------------------------------------
/** this function sets the sampled patches which may at any time be sampled (called once)
 * this depends on the following parameters listed in the order of priority:
 *  - patch_sample_size_fem/mal   (may vary over time)
 *  - patch_sample_size           (may vary over time)
 *  - sampled_patches             (may not vary over time)
 * patches sampled at any time have a sampleID of SAMPLED
 * never sampled patches have a sampleID of NaN
 * the parameter allGens defines if the sampled patches should be set across all times (i.e. if a parameter is temporal)
 *   by default its value is false.
 * returns true if all patches are sampled at any time
 */
bool
TMetapop::set_sampled_patches(bool allGens)
{
    // reset the patches state ///////////////////////////////////////////////////
    unsigned int i;
    vector<TPatch*>::iterator curPatch, endPatch;
    for(curPatch=_vPatch.begin(), endPatch=_vPatch.end(); curPatch!=endPatch; ++curPatch){
        (*curPatch)->set_sampleID(my_NAN);
    }
    
    // set by "patch_sample_size" ////////////////////////////////////////////////
    if(setSampleSizes(false)){ // dont set the default values if not set
        // assign the sampleID of all sampled patches
        TPatch* curPatch;
        for(i=0, _tot_sampled_patches=0; i<_patchNbr; ++i){
            curPatch = get_vPatch(i);
            if(curPatch->get_N_sample()==my_NAN || !curPatch->get_N_sample()) continue;
            curPatch->set_sampleID(SAMPLED);        // this patch is sampled
            ++_tot_sampled_patches;
        }
        if(!allGens) return _tot_sampled_patches==_patchNbr;
        
        // do the sex specific parameters change over time?
        if(_paramSet->get_param("patch_sample_size_fem")->isTemporalParam()
           || _paramSet->get_param("patch_sample_size_mal")->isTemporalParam())
        {
            map<int, string>::iterator curFem, endFem, curMal, endMal;
            Param* paramFem = _paramSet->get_param("patch_sample_size_fem");
            Param* paramMal = _paramSet->get_param("patch_sample_size_mal");
            curFem	= paramFem->get_temporal_args()->begin();
            endFem	= paramFem->get_temporal_args()->end();
            curMal	= paramMal->get_temporal_args()->begin();
            endMal	= paramMal->get_temporal_args()->end();
            
            // loop through the generations and set the new samples
            bool changed = false;
            for(i=2, ++curFem, ++curMal; i<=_generations && (curFem!=endFem || curMal!=endMal); ++i){ // for each generation (without the first one)
                if(i == (unsigned int) curFem->first){           // check for female changes
                    paramFem->update_arg(i);
                    ++curFem;
                    changed = true;
                }
                if(i == (unsigned int) curMal->first){           // check for male changes
                    paramMal->update_arg(i);
                    ++curMal;
                    changed = true;
                }
                if(!changed) continue;                           // nothing has chagned : continue
                
                // set the samples and update the sampleIDs
                if(setSampleSizes(false)){ // dont set the default values if not set
                    // assign the sampleID of all sampled patches
                    TPatch* curPatch;
                    for(i=0; i<_patchNbr; ++i){
                        curPatch = get_vPatch(i);
                        if(curPatch->get_N_sample()==my_NAN || !curPatch->get_N_sample()) continue;
                        curPatch->set_sampleID(SAMPLED);        // this patch is sampled
                    }
                }
            }
        }
        
        // does the general parameter change over time?
        else if(_paramSet->get_param("patch_sample_size")->isTemporalParam()){
            map<int, string>::iterator cur, end;
            Param* param = _paramSet->get_param("patch_sample_size");
            cur	= param->get_temporal_args()->begin();
            end	= param->get_temporal_args()->end();
            
            // loop through the generations and set the new samples
            for(i=2, ++cur; i<=_generations && cur!=end; ++i){ // for each generation (without the first one)
                if(i != (unsigned int) cur->first) continue;           // check for changes
                param->update_arg(i);
                ++cur;
                
                // set the samples and update the sampleIDs
                if(setSampleSizes(false)){ // dont set the default values if not set
                    // assign the sampleID of all sampled patches
                    TPatch* curPatch;
                    for(i=0; i<_patchNbr; ++i){
                        curPatch = get_vPatch(i);
                        if(curPatch->get_N_sample()==my_NAN || !curPatch->get_N_sample()) continue;
                        curPatch->set_sampleID(SAMPLED);        // this patch is sampled
                    }
                }
            }
        }
        
        // count the total sampled patches
        for(i=0, _tot_sampled_patches=0; i<_patchNbr; ++i){      // set the correct sampleID's
            if(get_vPatch(i)->get_sampleID()!=my_NAN) ++_tot_sampled_patches;
        }
        return _tot_sampled_patches==_patchNbr;
    }
    
    
    // set by "sampled_patches" //////////////////////////////////////////////////
    /** no sample sizes specified: check parameter "sampled_patches"
     * if matrix:       the numbers correspond to the ID of the patches to sample
     * if single value: the number correponds to the number of patches to randomly sample
     * if 0 (default):  all patches are sampled
     * Caution: input: id starts at 1 - in quantiNemo2 id starts at 0
     */
    // if the patches are explicitly defined by their ID
    if(_paramSet->get_param("sampled_patches")->is_matrix()){
        TMatrix* m = _paramSet->get_param("sampled_patches")->get_matrix();
        _tot_sampled_patches = m->get_dims(NULL);
        double* vec = m->get();
        for(i=0; i<_tot_sampled_patches; ++i){
            if(vec[i] < 1)              error("Parameter 'sampled_patches': patch index (%i) has to be positive!\n", vec[i]);
            if(vec[i] > (int)_patchNbr) error("Parameter 'sampled_patches': patch index (%i) is exceeding the number of patches (%i)!\n", vec[i], _patchNbr);
            get_vPatch((unsigned int)vec[i]-1)->set_sampleID(SAMPLED);
        }
        delete m;
        return _tot_sampled_patches==_patchNbr;
    }
    
    // it is a single value (0: default: sample all patches)
    _tot_sampled_patches = (unsigned int)_paramSet->getValue("sampled_patches");
    if(!_tot_sampled_patches) _tot_sampled_patches = _patchNbr;
    if(_tot_sampled_patches == _patchNbr){ // default: if all patches are sampled
        for(unsigned int i=0; i<_patchNbr; ++i){
            get_vPatch(i)->set_sampleID(SAMPLED);
        }
        return true;
    }
    
    // n patches are randomly selected
    if(_tot_sampled_patches > _patchNbr) error("Parameter 'sampled_patches' (%i) has to be less or equal the total number of patches (%i)!\n", _tot_sampled_patches, _patchNbr);
    unsigned int* vSample = rand().sample(_tot_sampled_patches, _patchNbr);   // draw UNIQUE numbers
    sort(vSample, vSample+_tot_sampled_patches);
    for(unsigned int i=0; i<_tot_sampled_patches; ++i){
        get_vPatch(vSample[i])->set_sampleID(SAMPLED);
    }
    delete[] vSample;
    return false;
}

// ----------------------------------------------------------------------------------------
// temporal_change
// ----------------------------------------------------------------------------------------
/** calls all "temporal_change" functions to be called */
void
TMetapop::temporal_change()
{
    // keywords
    if(_pReplicate->get_paramSetKeys()){
        map<string, Param*>* pParam = _pReplicate->get_paramSetKeys()->getTemporalParams(_currentGeneration);
        if(pParam) _pReplicate->get_paramSetKeys()->updateTemporalParams(_currentGeneration); // if it is a temporal parameter
    }
    
    // metapop
    temporal_change(_currentGeneration);
    
    // traits
    map< string,TTraitProto* >::iterator iter = getTraitPrototypes().begin();
    for(; iter != getTraitPrototypes().end(); ++iter){
        iter->second->temporal_change(_currentGeneration);
    }
    
    // LCEs
    vector<LCE*>::iterator iterLCE = _theCycle.begin();
    for(; iterLCE != _theCycle.end(); ++iterLCE){
        (*iterLCE)->temporal_change(_currentGeneration);
    }
    
    // selection
    if(_pSelection) _pSelection->temporal_change(_currentGeneration);
    
    // genetic map
    if(!isCoalescence()) _protoGenome->temporal_change(_currentGeneration);
}

// ----------------------------------------------------------------------------------------
// temporal_change
// ----------------------------------------------------------------------------------------
/** parameters which may change over time */
void
TMetapop::temporal_change(const unsigned int& gen)
{
    map<string, Param*>* pParam = _paramSet->getTemporalParams(gen);
    
    // if it is a temporal parameter
    if(pParam){
        // check if a change has to be made
        map<string, Param*>* pMap = _paramSet->updateTemporalParams(gen);
        if(pMap){
            // iterate through the map and perform the updates
            bool optima = false, intensity = false, carCap = false,
            growth_rate = false, max_growth = false, symmetry = false,
            min = false, max = false, meanVe = false, h2 = false, fit_land = false,
            stabVar = false, dirVar = false, samples = false;
            map<string, Param*>::iterator pos = pMap->begin();
            for(; pos != pMap->end(); ++pos){
                if(pos->first.find("patch_capacity") != string::npos)               	 carCap =      true;
                

                
                else if(pos->first.find("patch_sample_size") != string::npos)            samples =     true;
                
                else if(pos->first.find("patch_ve_mean") != string::npos)                meanVe =      true;
                else if(pos->first.find("patch_ve_var") != string::npos)                 h2 =          true;
                
                // stabilizing selection changes
                else if(pos->first.find("patch_stab_sel_optima_var") != string::npos
                        || pos->first.find("patch_stab_sel_intensity_var") != string::npos) stabVar=      true;
                else if(pos->first.find("patch_stab_sel_optima") != string::npos)        optima =      true;
                else if(pos->first.find("patch_stab_sel_intensity") != string::npos)     intensity =   true;
                
                // directional selection changes
                else if(pos->first.find("patch_dir_sel_growth_rate_var") != string::npos
                        || pos->first.find("patch_dir_sel_max_growth_var") != string::npos
                        || pos->first.find("patch_dir_sel_symmetry_var") != string::npos
                        || pos->first.find("patch_dir_sel_min_var") != string::npos
                        || pos->first.find("patch_dir_sel_max_var") != string::npos)         dirVar =     true;
                else if(pos->first.find("patch_dir_sel_growth_rate") != string::npos)     growth_rate= true;
                else if(pos->first.find("patch_dir_sel_max_growth") != string::npos)      max_growth = true;
                else if(pos->first.find("patch_dir_sel_symmetry") != string::npos)        symmetry =   true;
                else if(pos->first.find("patch_dir_sel_min") != string::npos)             min =        true;
                else if(pos->first.find("patch_dir_sel_max") != string::npos)             max =        true;
                
                // fitness landscape changes
                else if(pos->first.find("patch_fitness_landscape") != string::npos
                        || pos->first.find("patch_phenotype_landscape") != string::npos)     fit_land =   true;
                
                else if(pos->first == "selection_level_coef"){if(_pSelection) _pSelection->set_selection_level_coef();}
            }
            
            // make the changes only now, since several parameters may have changed (male and female...)
            if(carCap){
                setCarryingCapacities();
                if(_total_carrying_capacity != my_NAN) set_total_carrying_capacity();
            }
            
        
            
            if(samples){
                if(!_vSamplePatchUsed)	setSampleSizes(false);	// do not alter the sampleID's
                else{                                           // if the sampling is used: reset it
                    setSampleSizes(true);                         // set as well the sample IDs
                    set_sampledInds(ALL, false);                  // a newly not sampled pop has to be somehow reset (do not store the state as it will be recalled!!!)
                    _vSamplePatch.clear();
                    vector<TPatch*>::iterator curPop, endPop = get_vFullPatch().end();
                    unsigned int i=0;
                    for(curPop = get_vFullPatch().begin(); curPop!= endPop; ++curPop){
                        if((*curPop)->get_sampleID()==my_NAN) continue;
                        (*curPop)->set_sampleID(i++);
                        _vSamplePatch.push_back(*curPop);
                    }
                }
            }
            
            unsigned int linkedTraits = _vPatch[0]->get_nbLinkedTraits();
            if(meanVe) _pSelection->set_ve_mean();
            if(h2){
                set_patch_parameter(linkedTraits, "patch_ve_var", "heritability", &TPatch::set_localh2Ve);
                _pSelection->reset_Ve();
            }
            
            if(stabVar)     _pSelection->reset_selectionTypes();
            if(optima)      set_patch_parameter(linkedTraits, "patch_stab_sel_optima", "optima", &TPatch::set_localOptima);
            if(intensity)   set_patch_parameter(linkedTraits, "patch_stab_sel_intensity", "intensity", &TPatch::set_localIntensity);
            
            if(dirVar)     _pSelection->reset_selectionTypes();
            if(min)         set_patch_parameter(linkedTraits, "patch_dir_sel_min", "min", &TPatch::set_localMin);
            if(max)         set_patch_parameter(linkedTraits, "patch_dir_sel_max", "max", &TPatch::set_localMax);
            if(growth_rate) set_patch_parameter(linkedTraits, "patch_dir_sel_growth_rate", "growth rate", &TPatch::set_localGrowthRate);
            if(max_growth)  set_patch_parameter(linkedTraits, "patch_dir_sel_max_growth", "max growth", &TPatch::set_localMaxGrowth);
            if(symmetry)    set_patch_parameter(linkedTraits, "patch_dir_sel_symmetry", "symmetry", &TPatch::set_localSymmetry);
            
            if(fit_land){   // (phenotype has to be set after fitness due to sorting)
                set_patch_parameter_array(linkedTraits, "patch_fitness_landscape",   "fitness_landscape", &TPatch::set_fitnessLandscape_fitness);
                set_patch_parameter_array(linkedTraits, "patch_phenotype_landscape", "phenotype_landscape", &TPatch::set_fitnessLandscape_phenotype);
            }
        }
    }
}

// ----------------------------------------------------------------------------------------
// ini_globs
// ----------------------------------------------------------------------------------------
bool TMetapop::init_coal(map< string,TTraitProto* >& traits, map< int,LCE* >& LCEs)
{
#ifdef _DEBUG
    message("Metapop::init_coal\n");
#endif
    _generations = (unsigned int) _paramSet->getValue("generations");

    if(!_paramSet->isSet()) error("parameters in 'population' are not properly set!\n");
    set_LCE_pointers(LCEs);
    set_SexInitRatio(LCEs);
    createPopulations_coal();
    setCarryingCapacities();
    setInitPopulationSizes();
    setSampleSizes();
    set_sampled_patches(true);
    _sampleAllOrNothing = (bool)_paramSet->getValue("sample_all_or_nothing");
    if(allPatchPopulated_atInitialization()) set_func_ptr_noSample_noFull_coal();    // do not use _vFullPatch
    else                                     set_func_ptr_noSample_withFull_coal();  // use _vFullPatch
    makePrototype(traits, this);
    
    setLifeCycle_coal(LCEs);
    
    set_change_disp_rate_after_density();
    return true;
}

// ----------------------------------------------------------------------------------------
// ini_globs
// ----------------------------------------------------------------------------------------
bool TMetapop::init(map< string,TTraitProto* >& traits, map< int,LCE* >& LCEs)
{
#ifdef _DEBUG
    message("Metapop::init\n");
#endif
    
    _generations = (unsigned int) _paramSet->getValue("generations");
    
    if(!_paramSet->isSet()) error("parameters in 'population' are not properly set!\n");
    set_LCE_pointers(LCEs);
    set_SexInitRatio(LCEs);
    createPopulations();
    setCarryingCapacities();
    setInitPopulationSizes();
    setSampleSizes();
    if(set_sampled_patches(true)){	// if all patches are sampled:                // do not use _vSamplePatch
        if(allPatchPopulated_atInitialization()) set_func_ptr_noSample_noFull();    // do not use _vFullPatch
        else                                     set_func_ptr_noSample_withFull();  // use _vFullPatch
    }
    else{                                                                         // use _vSamplePatch
        if(allPatchPopulated_atInitialization()) set_func_ptr_withSample_noFull();  // do not use _vFullPatch
        else                                     set_func_ptr_withSample_withFull();// use _vFullPatch
    }
    makePrototype(traits, this);
    
    // initialize selection (if not use destroy the object)
    if(_pSelection) delete _pSelection;
    _pSelection = new TSelection();
    
    if(_paramSet->getValue("selection_pressure_definition")){
        if(!_pSelection->init2(this)){delete _pSelection; _pSelection=NULL;} // using quanti_x
    }
    else{
        if(!_pSelection->init(this)){delete _pSelection; _pSelection=NULL;}  // using patch_x (default)
    }
    
    setLifeCycle(LCEs);
    
    set_change_disp_rate_after_density();
    
    //empty and clean the RecyclingPOOL, safer...
    for(unsigned int i = 0; i < RecyclingPOOL.size(); ++i) {
        assert(RecyclingPOOL[i]);
        delete RecyclingPOOL[i];
    }
    RecyclingPOOL.clear();
    
    return true;
}

// ----------------------------------------------------------------------------------------
// set_sampledInds
// ----------------------------------------------------------------------------------------
/** for each generation with output the (stats or files) the individuals to sample have to be defined.
 * only the patches to sample are considered.
 * The parameter storeState allows to call the function without that the gen, rep and age are stored (used for temporal parameter)
 */
void
TMetapop::set_sampledInds(age_t AGE, bool storeState)
{
    _generations = (unsigned int) _paramSet->getValue("generations");
    
    // check if this generation was already sampled somehow
    if(AGE == NONE) return;                                    // no age class to update
    if(storeState){
        if((_isSampled[0]==_currentGeneration)  && (_isSampled[1]==getCurrentReplicate())){
            if(_isSampled[2]== AGE) return;                          // all the specified age classes were sampled: stop here
            
            // find out which age classes still have to be sampled
            unsigned int a, mask;
            for(a = 0, mask=1; a < NB_AGE_CLASSES; ++a, mask <<= 1){ // check each age class
                if((mask & AGE)                                        // if age class has to be sampled
                   && (mask & _isSampled[2])){                          //   and was already sampled: remove it
                    AGE ^= a;                                            // remove the age class
                }
            }
        }
        else{                                                      // nothing has been yet computed for this generation
            _isSampled[0]=_currentGeneration;                        // this generation
            _isSampled[1]=getCurrentReplicate();                         // this replicate
            _isSampled[2] = NONE;                                    // thus no age class has been sampeld yet
        }
    }
    
    assert(AGE!=NONE);
    
    
    TPatch* curPatch;
    vector<TPatch*>::iterator curPop, endPop = _popPtr->get_vSamplePatch().end();
    for(curPop=_popPtr->get_vSamplePatch().begin(); curPop!=endPop; ++curPop){
        curPatch = *curPop;
        curPatch->set_sampledInds(FEM, AGE);                    // females
        if(_sexInitRatio) curPatch->set_sampledInds(MAL, AGE);  // males
    }
    
    // update the age classes of the sampled pops
    _isSampled[2] |= AGE;                                     // add the age classes of this sampling
    
}

// ----------------------------------------------------------------------------------------
// get_current_nbPops
// ----------------------------------------------------------------------------------------
/** returns the number of current POPULATED pops of the given age class
 * Note: the vector _current_pops  contains all populated pops, thus not all
 *        pops have to be populated with all age classes
 *
 unsigned int
 TMetapop::get_current_nbPops(const age_idx& AGE)
 {
 unsigned int i, size=0;
 vector<Patch*>::iterator cur, end;
 for(cur=_current_pops.begin(), end=_current_pops.end(); cur!=end; ++cur){
 if((*cur)->sampleSize(AGE)) ++size;
 }
 return size;
 }
 */
// ----------------------------------------------------------------------------------------
// get_current_nbPops
// ----------------------------------------------------------------------------------------
/** returns the number of current POPULATED pops of the given age class
 * Note: the vector _current_pops  contains all populated pops, thus not all
 *        pops have to be populated with all age classes
 *
 unsigned int
 TMetapop::get_current_nbPops(const age_t& AGE)
 {
 unsigned int size=0;
 vector<Patch*>::iterator cur, end;
 for(cur=_current_pops.begin(), end=_current_pops.end(); cur!=end; ++cur){
 if((*cur)->sampleSize(AGE)) ++size;
 }
 return size;
 }
 */
//---------------------------------------------------------------------------
/**Removes all individual pointers of the appropriate sex and age class and flush them into the recycling pool.
 \b Note: not memory operation is performed, the total amount of memory allocated is
 left untouched.
 @param SEX the sex class of the individual
 @param AGE the index of the age class
 @param pop the pointer to the metapop for access to the recycling pool
 */
void
TMetapop::flush(sex_t SEX, age_idx AGE)
{
    vector<TPatch*>::iterator curPatch, endPatch;
    for(curPatch=get_vFullPatch().begin(), endPatch=get_vFullPatch().end(); curPatch!=endPatch; ++curPatch){
        (*curPatch)->flush(SEX, AGE);
    }
}

//---------------------------------------------------------------------------
void
TMetapop::flush(age_idx AGE)
{
    flush(FEM, AGE);
    flush(MAL, AGE);
}

//---------------------------------------------------------------------------
/**Removes all individual pointers of the appropriate sex and age class and flush them into the recycling pool.
 @param AGE an unsigned int containing the flags of the age classes to flush
 @param pop the pointer to the metapop for access to the recycling pool
 @see flush()
 */
void
TMetapop::flush(age_t AGE)
{
    for(unsigned int i = 0, mask = 1; i < NB_AGE_CLASSES; i++, mask <<= 1) {
        if(mask & AGE) {
            flush(MAL, static_cast<age_idx>(i));
            flush(FEM, static_cast<age_idx>(i));
        }
    }
}

//---------------------------------------------------------------------------
/**Removes all individual pointers of all sex and age classes and flush them into the recycling pool.*/
void
TMetapop::flush(){
    for(unsigned int i = 0; i < NB_AGE_CLASSES; i++)
    {
        flush(static_cast<age_idx>(i));
    }
}

// ----------------------------------------------------------------------------------------
// set_LCE_pointers
// ----------------------------------------------------------------------------------------
/** sets the pointer to the LCE breeding and dispersal */
void
TMetapop::set_LCE_pointers(map< int,LCE* >& LCEs)
{
    _pDisperse_LCE = NULL;
    _pBreed_LCE = NULL;
    _pCoalescence_LCE = NULL;
    map<int, LCE* >::iterator cur, end;
    for(cur=LCEs.begin(), end=LCEs.end(); cur!=end; ++cur) {
        if(!_pDisperse_LCE)   _pDisperse_LCE    = dynamic_cast<LCE_Disperse*>(cur->second);
        if(!_pBreed_LCE)      _pBreed_LCE       = dynamic_cast<LCE_Breed*>(cur->second);
        if(!_pCoalescence_LCE)_pCoalescence_LCE = dynamic_cast<LCE_Coalescence_base*>(cur->second); // is not a must!!!
    }
    assert(_pDisperse_LCE && _pBreed_LCE);       // the breed and the dispersal LCE must be present!
}

// ----------------------------------------------------------------------------------------
// buildPopulation
// ----------------------------------------------------------------------------------------
/** sets the initial sex ratio. For hermaphrodites it is 0, if both sexes are equally set it is 1
 * input parameter "sex_ratio"=nMmal/nbFem                1:=  nbFem=nbMal
 * internal parameter "_sexInitRatio"=nbMal/(nbMal+nbFem) 0.5:=nbFem=nbMal
 */
void TMetapop::set_SexInitRatio(map< int,LCE* >& LCEs)
{
    assert(_pBreed_LCE);     // the breed LCE must be present!
    if(_pBreed_LCE->get_parameter_value("mating_system") < 3
       || _pBreed_LCE->get_parameter_value("mating_system") == 6){     // for stats
        _sexInitRatio = 0;                      // hermaphrodites
    }
    else {
        _sexInitRatio = _pBreed_LCE->get_parameter_value("sex_ratio");
        if(!_sexInitRatio) error("Two sexes with a sex ratio of 0 make no sense!\n");
        _sexInitRatio /= _sexInitRatio+1;           // transform the sex ratio
    }
}

// ----------------------------------------------------------------------------------------
// buildPopulation
// ----------------------------------------------------------------------------------------
void TMetapop::createPopulations()
{
    assert(_vFullPatch.empty());
    assert(_vSamplePatch.empty());
    assert(_vPatch.empty());
    
    // get the number of populations
    _patchNbr = (unsigned int)_paramSet->getValue("patch_number");
    
    // create the patches
    _vPatch.reserve(_patchNbr);
    for(unsigned int i=0; i<_patchNbr; ++i){
        _vPatch.push_back(new TPatch(this, i));
    }
    
    set_func_pointer();
}

// ----------------------------------------------------------------------------------------
// buildPopulation
// ----------------------------------------------------------------------------------------
void TMetapop::createPopulations_coal()
{
    assert(_vFullPatch.empty());
    assert(_vSamplePatch.empty());
    assert(_vPatch.empty());
 
    // get the number of populations
    _patchNbr = (unsigned int)_paramSet->getValue("patch_number");
 
    _vPatch.reserve(_patchNbr);
    for(unsigned int i=0; i<_patchNbr; ++i){
        _vPatch.push_back(new TPatch(this, i));
        _vPatch[i]->init_coal(i);
    }
    
    set_func_pointer_coal();
}

// ----------------------------------------------------------------------------------------
// setPatchParam
// ----------------------------------------------------------------------------------------
/** this is the main fucntion to set a patch parameter. The parameter may be set
 * for each sex separately or for both sexes together.
 * The following parameters may use this functions:
 *   - patch_capacity
 *   - patch_ini_size
 *   - patch_sample_size
 * returns false if none of the parameters was explicitly set and true if it was
 */
template<typename T>
bool
TMetapop::setPatchParam(string name,                                // parameter name (without "_fem" or "_mal")
                       void (TPatch::*setSexSpecific)(T, sex_t),   // set a sex specific value
                       void (TPatch::*setGeneral)(T),              // set a general value
                       void (TPatch::*reset)(T, T, T),             // function to set all three params at once
                       T (TPatch::*getSexSpecific)(sex_t),         // get a sex specific value
                       T (TPatch::*getGeneral)())
{
    
    // hermaphrodites
    if(!_sexInitRatio){
        if(setPatchParam(name+"_fem", FEM, setSexSpecific, reset)){    // sex specific parameter
            setGeneralPatchParam_fem<T>(getSexSpecific, setGeneral);     // if only females are used
            return true;
        }
        if(setPatchParam(name, FEM, setSexSpecific, reset)){					 // general parameter
            setGeneralPatchParam_fem<T>(getSexSpecific, setGeneral);     // if only females are used
            return true;
        }
        return false;
    }
    
    // two sexes
    // sex specific settings
    if(setPatchParam(name+"_fem", FEM, setSexSpecific, reset)){                  // sex specific parameter
        if(setPatchParam(name+"_mal", MAL, setSexSpecific, reset)){
            setGeneralPatchParam_sum<T>(getSexSpecific, setGeneral);                 // if both sexes are set
            return true;
        }
        else	error("Parameter %s: only fem sex specific parameter is specified: both are required!\n", name.c_str());
    }
    if(setPatchParam(name+"_mal", MAL, setSexSpecific, reset)){
        error("Parameter %s: only mal sex specific parameter is specified: both are required!\n", name.c_str());
    }
    
    // general parameter
    if(setPatchParam(name, setGeneral, reset)){
        setSexSpecificPatchParam_sexRatio<T>(getGeneral, setSexSpecific);
        return true;
    }
    return false;
}


// ----------------------------------------------------------------------------------------
// setPatchParam
// ----------------------------------------------------------------------------------------
/** this function sets the given general parameter if set or desired
 * returns true if the parameter was set and false if the parameter was not set
 */
template<typename T>
bool
TMetapop::setPatchParam(string name, sex_t SEX,
                       void (TPatch::*pt2Func)(T, sex_t),
                       void (TPatch::*pt2reset)(T, T, T))
{
    Param* pParam = _paramSet->get_param(name);
    if(!pParam->isSet()) return false;                                                         // not set
    if(pParam->is_matrix()) setPatchParam(name, pParam->get_matrix(), SEX, pt2Func, pt2reset); // a matrix
    else                    setPatchParam(name, (T)pParam->get_value(), SEX, pt2Func);         // a single value
    return true;
}

// ----------------------------------------------------------------------------------------
// setPatchParam
// ----------------------------------------------------------------------------------------
/** the sex specific value for each patch is given by a matrix,
 * either for each patch a value (singel matrix)
 * or each patch is directly addressed (2D matrix with two columns, schema {patchID value})
 * note, that the the valuea where already reset.
 */
template<typename T>
void
TMetapop::setPatchParam(string name, TMatrix* m, sex_t SEX,
                       void (TPatch::*pt2Func)(T, sex_t),
                       void (TPatch::*pt2reset)(T, T, T))
{
    if(!m->get_oneDims() && m->getNbCols()==2){        // the patches are directly addressed ({{patchID value}{patchID value}})
        unsigned int i, id,  nbRows = m->getNbRows();
        
        // reset all patches
        if(SEX==FEM){               // do the initial setting only once, thus if males don't do it again since it has been done for females
            for(i = 0; i < _patchNbr; ++i){
                _vPatch[i]->set_localParameter<T>(0,0,0,pt2reset);
            }
        }
        
        // set the individualvalues
        for(i=0; i<nbRows; ++i){
            id = (unsigned int)m->get(i, 0) - 1;
            if(id>=_patchNbr) error("Parameter %s: the patch ID exceeds the number of patches!\n", name.c_str());
            _vPatch[id]->set_localParameter<T>((T)m->get(i,1), SEX, pt2Func);
        }
    }
    else{                         // a matrix is passed
        // get the vector and its size
        unsigned int size = m->get_dims(NULL);
        double* vec = m->get();
        
        // check the matrix sizes with the number of patches
        if(size>_patchNbr) warning("Parameter %s: The matrix length exceeds the number of patches! Not all values will be used!\n", name.c_str());
        else if(_patchNbr % size) warning("Parameter %s: The matrix length is not an entire subset of the number of patches!\n", name.c_str());
        
        for(unsigned int i = 0; i < _patchNbr; ++i){
            _vPatch[i]->set_localParameter<T>((T)vec[i % size], SEX, pt2Func);
        }
    }
    
    // clean up
    delete m;
}

// ----------------------------------------------------------------------------------------
// setPatchParam
// ----------------------------------------------------------------------------------------
// all patches have the same sex specific value
template<typename T>
void
TMetapop::setPatchParam(string name, T value, sex_t SEX,
                       void (TPatch::*pt2Func)(T, sex_t))
{
    for(unsigned int i = 0; i < _patchNbr; ++i){
        _vPatch[i]->set_localParameter(value, SEX, pt2Func);
    }
}

// ----------------------------------------------------------------------------------------
// setPatchParam
// ----------------------------------------------------------------------------------------
/** this function sets the given general parameter if set or desired
 * returns true if the parameter was set and false if the parameter was not set
 */
template<typename T>
bool
TMetapop::setPatchParam(string name,
                       void (TPatch::*pt2Func)(T),
                       void (TPatch::*pt2reset)(T, T, T))
{
    Param* pParam = _paramSet->get_param(name);
    if(!pParam->isSet()) return false;                                                    // not set
    if(pParam->is_matrix()) setPatchParam(name, pParam->get_matrix(), pt2Func, pt2reset); // a matrix
    else                    setPatchParam(name, (T)pParam->get_value(), pt2Func);         // a single value
    return true;
}

// ----------------------------------------------------------------------------------------
// setPatchParam
// ----------------------------------------------------------------------------------------
// the general patch values are defined by a matrix
template<typename T>
void
TMetapop::setPatchParam(string name, TMatrix* m,
                       void (TPatch::*pt2Func)(T),
                       void (TPatch::*pt2reset)(T, T, T))
{
    if(!m->get_oneDims() && m->getNbCols()==2){        // the patches are directly addressed ({{patchID value}{patchID value}})
        unsigned int i, id,  nbRows = m->getNbRows();
        
        // reset all patches
        for(i = 0; i < _patchNbr; ++i){
            _vPatch[i]->set_localParameter<T>(0,0,0,pt2reset);
        }
        
        // set the individualvalues
        for(i=0; i<nbRows; ++i){
            id = (unsigned int)m->get(i, 0) - 1;
            if(id>=_patchNbr) error("Parameter %s: the patch ID exceeds the number of patches!\n", name.c_str());
            _vPatch[id]->set_localParameter<T>((T)m->get(i,1), pt2Func);
        }
    }
    else{                         // a matrix is passed
        // get the vector and its size
        unsigned int size = m->get_dims(NULL);
        double* vec = m->get();
        
        // check the matrix sizes with the number of patches
        if(size>_patchNbr) warning("Parameter %s: The matrix length exceeds the number of patches! Not all values will be used!\n", name.c_str());
        else if(_patchNbr % size) warning("Parameter %s: The matrix length is not an entire subset of the number of patches!\n", name.c_str());
        
        for(unsigned int i = 0; i < _patchNbr; ++i){
            _vPatch[i]->set_localParameter<T>((T)vec[i % size], pt2Func);
        }
    }
    
    // clean up
    delete m;
}

// ----------------------------------------------------------------------------------------
// setPatchParam
// ----------------------------------------------------------------------------------------
// all patches have the same general value
template<typename T>
void
TMetapop::setPatchParam(string name,	T value,
                       void (TPatch::*pt2Func)(T))
{
    for(unsigned int i = 0; i < _patchNbr; ++i){
        _vPatch[i]->set_localParameter(value, pt2Func);
    }
}

// ----------------------------------------------------------------------------------------
// set_patch_selection_slope
// ----------------------------------------------------------------------------------------
/** sets patch specific parameters for each quantitative trait. They may be passed as:
 *  - sex specific (in this case both sexes (if there are two sexes) have to be specified)
 *  - value: same value for all patches
 *  - matrix: matrixes may be expanded to meet the number of patches and traits
 */
void TMetapop::set_patch_parameter(unsigned int nbTrait,
                                  string name,             // "patch_dir_sel_slope"
                                  string name_full,        // "slope"
                                  void (TPatch::*pt2Func)(double*, sex_t))
{
    // sex specific settings have precedence
    if(_paramSet->isSet(name+"_fem")){
        if(_paramSet->isMatrix(name+"_fem")){  // if matrix
            set_patch_value_byMatrix(nbTrait, _paramSet->getMatrix(name+"_fem"), FEM, name_full, pt2Func);
        }
        else{                                  // if single value
            set_patch_value_byValue(nbTrait, _paramSet->getValue(name+"_fem"), FEM, pt2Func);
        }
        
        if(_sexInitRatio){
            if(_paramSet->isSet(name+"_mal")){
                if(_paramSet->isMatrix(name+"_mal")){    // if matrix
                    set_patch_value_byMatrix(nbTrait, _paramSet->getMatrix(name+"_mal"), MAL, name_full, pt2Func);
                }
                else{                                    // if single value
                    set_patch_value_byValue(nbTrait, _paramSet->getValue(name+"_mal"), MAL, pt2Func);
                }
            }
            else error("Only one sex specific %s is specified: both are required!\n", name_full.c_str());
        }
    }
    else{                         // general settings
        if(_paramSet->isMatrix(name)){
            set_patch_value_byMatrix(nbTrait, _paramSet->getMatrix(name), FEM, name_full, pt2Func);
            if(_sexInitRatio) set_patch_value_byMatrix(nbTrait, _paramSet->getMatrix(name), MAL, name_full, pt2Func);
        }
        else{
            set_patch_value_byValue(nbTrait, _paramSet->getValue(name), FEM, pt2Func);
            if(_sexInitRatio) set_patch_value_byValue(nbTrait, _paramSet->getValue(name), MAL, pt2Func);
        }
    }
}

// ----------------------------------------------------------------------------------------
// set_patch_selection_slope
// ----------------------------------------------------------------------------------------
/** sets patch specific parameters for a given quantitative trait. They may be passed as:
 *  - sex specific (in this case both sexes (if there are two sexes) have to be specified)
 *  - value: same value for all patches
 *  - matrix: matrixes may be expanded to meet the number of patches and traits
 */
void TMetapop::set_patch_parameter_ofTrait(TTraitProto* pTrait,
                                          unsigned int curTrait,   // 1
                                          string trait,            // "_2"
                                          string name,             // "patch_dir_sel_slope"
                                          string name_full,        // "slope"
                                          void (TPatch::*pt2Func)(unsigned int, double, sex_t))
{
    // sex specific settings have precedence
    if(pTrait->get_parameter_isSet(name+"_fem"+trait)){
        if(pTrait->get_parameter_isMatrix(name+"_fem"+trait)){  // if matrix
            set_patch_value_byMatrix_ofTrait(curTrait, pTrait->get_parameter_matrix(name+"_fem"+trait), FEM, name_full, pt2Func);
        }
        else{                                  // if single value
            set_patch_value_byValue_ofTrait(curTrait,pTrait->get_parameter_value(name+"_fem"+trait), FEM, pt2Func);
        }
        
        if(_sexInitRatio){
            if(pTrait->get_parameter_isSet(name+"_mal"+trait)){
                if(pTrait->get_parameter_isMatrix(name+"_mal"+trait)){    // if matrix
                    set_patch_value_byMatrix_ofTrait(curTrait, pTrait->get_parameter_matrix(name+"_mal"+trait), MAL, name_full, pt2Func);
                }
                else{                                    // if single value
                    set_patch_value_byValue_ofTrait(curTrait, pTrait->get_parameter_value(name+"_mal"+trait), MAL, pt2Func);
                }
            }
            else error("Only one sex specific %s is specified: both are required!\n", name_full.c_str());
        }
    }
    else{                         // general settings
        if(pTrait->get_parameter_isMatrix(name+trait)){
            set_patch_value_byMatrix_ofTrait(curTrait, pTrait->get_parameter_matrix(name+trait), FEM, name_full, pt2Func);
            if(_sexInitRatio) set_patch_value_byMatrix_ofTrait(curTrait, pTrait->get_parameter_matrix(name+trait), MAL, name_full, pt2Func);
        }
        else{
            set_patch_value_byValue_ofTrait(curTrait, pTrait->get_parameter_value(name+trait), FEM, pt2Func);
            if(_sexInitRatio) set_patch_value_byValue_ofTrait(curTrait, pTrait->get_parameter_value(name+trait), MAL, pt2Func);
        }
    }
}

// ----------------------------------------------------------------------------------------
// set_patch_value_by_value
// ----------------------------------------------------------------------------------------
/** this function sets the same value for each patch and trait*/
void TMetapop::set_patch_value_byValue(unsigned int nbTrait,
                                      double value,
                                      sex_t SEX,
                                      void (TPatch::*pt2Func)(double*, sex_t))
{
    // get the slope for each patch/trait
    double *array = new double[nbTrait];
    for(unsigned int t = 0; t < nbTrait; ++t) {
        array[t] = value;
    }
    
    // set each patch
    for(unsigned int p = 0; p < _patchNbr; ++p) {
        _vPatch[p]->set_localParameter(array, SEX, pt2Func);
    }
    delete[] array;
}

// ----------------------------------------------------------------------------------------
// set_patch_value_by_value
// ----------------------------------------------------------------------------------------
/** this function sets the same value for each patch and trait*/
void TMetapop::set_patch_value_byValue_ofTrait(unsigned int curTrait,
                                              double value,
                                              sex_t SEX,
                                              void (TPatch::*pt2Func)(unsigned int, double, sex_t))
{
    // set each patch
    for(unsigned int p = 0; p < _patchNbr; ++p) {
        _vPatch[p]->set_localParameter_ofTrait(curTrait, value, SEX, pt2Func);
    }
}

// ----------------------------------------------------------------------------------------
// v
// ----------------------------------------------------------------------------------------
/** this function sets the value to each patch based on a matrix */
void TMetapop::set_patch_value_byMatrix(unsigned int nbTrait,
                                       TMatrix* m,
                                       sex_t SEX,
                                       string name_full,        // "slope"
                                       void (TPatch::*pt2Func)(double*, sex_t))
{
    unsigned int p, t;
    
    unsigned int dims[2];
    m->get_dims(&dims[0]);
    unsigned int count_trait = dims[0];   //nbr trait slope per patch
    unsigned int count_patch = dims[1];   //nbr of patches slope
    
    // check the number of slope per number of patches
    if(count_patch>_patchNbr) warning("There are more %s specified than patches! Only a part of the %s is considered!\n", name_full.c_str(), name_full.c_str());
    else if(_patchNbr % count_patch) warning("The number of %s is not a entire subset of the number of patches!\n", name_full.c_str());
    
    // check the number of slope per number of traits
    if(count_trait>nbTrait) warning("There are more %s specified than traits! Only a part of the %s is considered!\n", name_full.c_str(), name_full.c_str());
    else if(nbTrait % count_trait) warning("The number of %s is not a entire subset of the number of traits!\n", name_full.c_str());
    
    // create the arrays
    double* array = new double[nbTrait];
    
    // get the optima for each patch/trait
    for(p = 0; p < _patchNbr; ++p) {
        for(t = 0; t < nbTrait; ++t) {
            array[t] = m->get(t % count_trait, p % count_patch);
        }
        _vPatch[p]->set_localParameter(array, SEX, pt2Func);
    }
    
    delete[] array;
    delete m;
}

// ----------------------------------------------------------------------------------------
// set_patch_array_byMatrix_ofTrait
// ----------------------------------------------------------------------------------------
/** this function sets the array to each patch based on a matrix */
void TMetapop::set_patch_array_byMatrix_ofTrait(unsigned int curTrait,
                                               TMatrix* m,
                                               sex_t SEX,
                                               string name_full,        // "slope"
                                               void (TPatch::*pt2Func)(unsigned int, double*, unsigned int, sex_t))
{
    unsigned int p;
    
    unsigned int dims[2];
    m->get_dims(&dims[0]);
    unsigned int count_vals = dims[0];   //nbr trait slope per patch
    unsigned int count_patch = dims[1];   //nbr of patches slope
    
    // check the number of slope per number of patches
    if(count_patch>_patchNbr) warning("There are more %s specified than patches! Only a part of the %s is considered!\n", name_full.c_str(), name_full.c_str());
    else if(_patchNbr % count_patch) warning("The number of %s is not a entire subset of the number of patches!\n", name_full.c_str());
    
    // get the vals for each patch
    double* array = new double[count_vals];
    for(p = 0; p < _patchNbr; ++p) {
        m->getRowView(p % count_patch, count_vals, array);
        _vPatch[p]->set_localParameter_matrix_ofTrait(curTrait, array, count_vals, SEX, pt2Func);
    }
    
    delete[] array;
    delete m;
}

// ----------------------------------------------------------------------------------------
// set_patch_value_byMatrix_ofTrait
// ----------------------------------------------------------------------------------------
/** this function sets the value to each patch based on a matrix */
void TMetapop::set_patch_value_byMatrix_ofTrait(unsigned int curTrait,
                                               TMatrix* m,
                                               sex_t SEX,
                                               string name_full,        // "slope"
                                               void (TPatch::*pt2Func)(unsigned int t, double, sex_t))
{
    double* mat = m->get();
    unsigned int size = m->length();
    
    // check the number of slope per number of patches
    if(size>_patchNbr) warning("There are more %s specified than patches! Only a part of the %s is considered!\n", name_full.c_str(), name_full.c_str());
    else if(_patchNbr % size) warning("The number of %s is not a entire subset of the number of patches!\n", name_full.c_str());
    
    // get the optima for each patch/trait
    for(unsigned int p = 0; p < _patchNbr; ++p) {
        _vPatch[p]->set_localParameter_ofTrait(curTrait, mat[p % size], SEX, pt2Func);
    }
    //delete[] mat;
    delete m;
}

// ----------------------------------------------------------------------------------------
// set_patch_parameter_array
// ----------------------------------------------------------------------------------------
/** sets patch specific parameters for each quantitative trait. They may be passed as:
 *  - sex specific (in this case both sexes (if there are two sexes) have to be specified)
 *  - value: same value for all patches
 *  - matrix: matrixes may be expanded to meet the number of patches and traits
 */
void TMetapop::set_patch_parameter_array(unsigned int nbTrait,
                                        string name,             // "patch_dir_sel_slope"
                                        string name_full,        // "slope"
                                        void (TPatch::*pt2Func)(unsigned int, double*, unsigned int, sex_t))
{
    TTree<unsigned int, double>* tree;
    
    // sex specific settings have precedence
    if(_paramSet->isSet(name+"_fem")){
        tree = _paramSet->get_param(name+"_fem")->get_tree<double>();
        if(tree->get_depth() == 3)      set_patch_array_byMatrix(nbTrait, tree, FEM, name_full, pt2Func);
        else if(tree->get_depth() == 1) set_patch_array_byArray(nbTrait, tree, FEM, pt2Func);
        else error("The dimension of the matrix is not valid!\n", name_full.c_str());
        delete tree;
        
        if(_sexInitRatio){
            if(_paramSet->isSet(name+"_mal")){
                tree = _paramSet->get_param(name+"_mal")->get_tree<double>();
                if(tree->get_depth() == 3)      set_patch_array_byMatrix(nbTrait, tree, MAL, name_full, pt2Func);
                else if(tree->get_depth() == 1) set_patch_array_byArray(nbTrait, tree, MAL, pt2Func);
                else error("The dimension of the matrix is not valid!\n", name_full.c_str());
                delete tree;
            }
            else error("Only one sex specific %s is specified: both are required!\n", name_full.c_str());
        }
    }
    else{                         // general settings
        tree = _paramSet->get_param(name)->get_tree<double>();
        if(tree->get_depth() == 3){
            set_patch_array_byMatrix(nbTrait, tree, FEM, name_full, pt2Func);
            if(_sexInitRatio) set_patch_array_byMatrix(nbTrait, tree, MAL, name_full, pt2Func);
        }
        else if(tree->get_depth() == 1){
            set_patch_array_byArray(nbTrait, tree, FEM, pt2Func);
            if(_sexInitRatio) set_patch_array_byArray(nbTrait, tree, MAL, pt2Func);
        }
        else error("The dimension of the matrix is not valid!\n", name_full.c_str());
        delete tree;
    }
}

// ----------------------------------------------------------------------------------------
// set_patch_parameter_array
// ----------------------------------------------------------------------------------------
/** sets patch specific parameters for each quantitative trait. They may be passed as:
 *  - sex specific (in this case both sexes (if there are two sexes) have to be specified)
 *  - value: same value for all patches
 *  - matrix: matrixes may be expanded to meet the number of patches and traits
 */
void TMetapop::set_patch_parameter_array_ofTrait(TTraitProto* pTrait,
                                                unsigned int curTrait,
                                                string trait,
                                                string name,             // "patch_dir_sel_slope"
                                                string name_full,        // "slope"
                                                void (TPatch::*pt2Func)(unsigned int t, double*, unsigned int, sex_t))
{
    // sex specific settings have precedence
    TMatrix* m;
    if(pTrait->get_parameter_isSet(name+"_fem"+trait)){
        m = pTrait->get_parameter_matrix(name+"_fem"+trait);        // get the female matrix
        if(!m) error("Parameter %s is not set, but needed!\n", (name+"_fem"+trait).c_str());
        if(m->get_oneDims()) set_patch_array_byArray_ofTrait(curTrait, m, FEM, name_full, pt2Func);
        else                 set_patch_array_byMatrix_ofTrait(curTrait, m, FEM, name_full, pt2Func);
        
        if(_sexInitRatio){
            if(pTrait->get_parameter_isSet(name+"_mal"+trait)){
                m = pTrait->get_parameter_matrix(name+"_mal"+trait);    // get the male matrix
                if(!m) error("Parameter %s is not set, but needed!\n", (name+"_mal"+trait).c_str());
                if(m->get_oneDims()) set_patch_array_byArray_ofTrait(curTrait, m, MAL, name_full, pt2Func);
                else                 set_patch_array_byMatrix_ofTrait(curTrait, m, MAL, name_full, pt2Func);
            }
            else error("Only one sex specific %s is specified: both are required!\n", name_full.c_str());
        }
    }
    else{                         // general settings
        m = pTrait->get_parameter_matrix(name+trait);               // get the general matrix
        if(!m) error("Parameter %s is not set, but needed!\n", (name+trait).c_str());
        if(m->get_oneDims()) set_patch_array_byArray_ofTrait(curTrait, m, FEM, name_full, pt2Func);
        else                 set_patch_array_byMatrix_ofTrait(curTrait, m, FEM, name_full, pt2Func);
        m = pTrait->get_parameter_matrix(name+trait);               // get the general matrix
        if(_sexInitRatio){
            if(m->get_oneDims()) set_patch_array_byArray_ofTrait(curTrait, m, MAL, name_full, pt2Func);
            else                 set_patch_array_byMatrix_ofTrait(curTrait, m, MAL, name_full, pt2Func);
        }
    }
}

// ----------------------------------------------------------------------------------------
// set_patch_array_by_array
// ----------------------------------------------------------------------------------------
/** this function sets the same array for each patch and trait*/
void TMetapop::set_patch_array_byArray(unsigned int nbTrait,
                                      TTree<unsigned int, double>* m,
                                      sex_t SEX,
                                      void (TPatch::*pt2Func)(unsigned int, double*, unsigned int, sex_t))
{
    assert(m->get_depth()==1);
    
    // get the unique array
    unsigned int p, t, size= 0;
    double *array=NULL;
    m->get_elements_as_array(array, size);    // size and array are adjusted if needed
    
    // set the slope for each patch/trait
    for(t = 0; t < nbTrait; ++t) {
        for(p = 0; p < _patchNbr; ++p) {
            _vPatch[p]->set_localParameter_matrix(t, array, size, SEX, pt2Func);
        }
    }
    delete[] array;
}

// ----------------------------------------------------------------------------------------
// v
// ----------------------------------------------------------------------------------------
/** this function sets the array to each patch based on a 1D matrix */
void TMetapop::set_patch_array_byArray_ofTrait(unsigned int curTrait,
                                              TMatrix* m,
                                              sex_t SEX,
                                              string name_full,        // "slope"
                                              void (TPatch::*pt2Func)(unsigned int, double*, unsigned int, sex_t))
{
    assert(m->get_oneDims());
    
    for(unsigned int p = 0; p < _patchNbr; ++p) {
        _vPatch[p]->set_localParameter_matrix_ofTrait(curTrait, m->get(), m->length(), SEX, pt2Func);
    }
    delete m;
}

// ----------------------------------------------------------------------------------------
// set_patch_array_by_matrix
// ----------------------------------------------------------------------------------------
/** this function sets the value to each patch based on a matrix
 * tree[trait][patch][value]
 */

void TMetapop::set_patch_array_byMatrix(unsigned int nbTrait,
                                       TTree<unsigned int, double>* m,
                                       sex_t SEX,
                                       string name_full,        // "slope"
                                       void (TPatch::*pt2Func)(unsigned int, double*, unsigned int, sex_t))
{
    assert(m->get_depth()==3);
    // get and check the number of traits
    TNode<unsigned int, double>** arrayTrait = NULL;
    unsigned int count_trait = 0;              // number of traits per patch
    m->get_nodes_as_array_static(arrayTrait, count_trait);     // arrayTrait has to be deleted manually
    if(count_trait>nbTrait) warning("There are more %s specified than traits! Only a part of the %s is considered!\n", name_full.c_str(), name_full.c_str());
    else if(nbTrait % count_trait) warning("The number of %s is not a entire subset of the number of traits!\n", name_full.c_str());
    
    
    // create the arrays
    double *array = NULL;
    TNode<unsigned int, double>** arrayPatch = NULL;
    unsigned int count_patch = 0;
    unsigned int p, t, size = 0;
    TNode<unsigned int, double>* curTrait, *curPatch;
    
    // get the optima for each patch/trait
    for(t = 0; t < nbTrait; ++t) {                                 // for each trait
        curTrait = arrayTrait[t % count_trait];                       // get the current trait
        if(!curTrait){                                                // not set
            for(p = 0; p < _patchNbr; ++p){                             // for each patch
                _vPatch[p]->set_localParameter_matrix(t, NULL, 0, SEX, pt2Func);// set the values
            }
        }
        curTrait->get_nodes_as_array_static(arrayPatch, count_patch); // size and array are adjusted if needed
        if(count_patch>_patchNbr) warning("There are more %s specified than patches! Only a part of the %s is considered!\n", name_full.c_str(), name_full.c_str());
        else if(_patchNbr % count_patch) warning("The number of %s is not a entire subset of the number of patches!\n", name_full.c_str());
        
        for(p = 0; p < _patchNbr; ++p){                               // for each patch
            curPatch = arrayPatch[p % count_patch];                     // get the current patch
            if(!curPatch){                                              // not set
                _vPatch[p]->set_localParameter_matrix(t, NULL, 0, SEX, pt2Func);// set the values
                continue;
            }
            curPatch->get_elements_as_array(array, size);               // size and array are adjusted if needed (array has to be deleted manually)
            _vPatch[p]->set_localParameter_matrix(t, array, size, SEX, pt2Func);// set the values
        }
        
    }
    
    // delete the arrays
    delete[] arrayPatch;
    delete[] arrayTrait;
    delete[] array;
}

// ----------------------------------------------------------------------------------------
// setLifeCycle
// ----------------------------------------------------------------------------------------
/** */
void TMetapop::setLifeCycle_coal(map< int, LCE*>& lifeCycle)
{
    string name;
    _theCycle.clear();
    _theCycleCoal.clear();
    
    for(map< int,LCE* >::iterator LCE = lifeCycle.begin() ; LCE != lifeCycle.end(); ++LCE) {
        if(LCE->second->init(this)){
            if(LCE->second->get_rank()<100) _theCycle.push_back(LCE->second);
            else                            _theCycleCoal.push_back(LCE->second);
        }
    }
}

// ----------------------------------------------------------------------------------------
// setLifeCycle
// ----------------------------------------------------------------------------------------
/** */
void TMetapop::setLifeCycle(map< int, LCE*>& lifeCycle)
{
    string name;
    _theCycle.clear();
    
    for(map< int,LCE* >::iterator LCE = lifeCycle.begin() ; LCE != lifeCycle.end(); ++LCE) {
        if(LCE->second->init(this)){
            _theCycle.push_back(LCE->second);
        }
    }
}

// ----------------------------------------------------------------------------------------
// Loop_generation
// ----------------------------------------------------------------------------------------
/* performs the generation loop for current replicate */
void TMetapop::Loop_generation(time_t startTime)
{
    vector< LCE* >::iterator LCE;
    
#ifdef _DEBUG
    message("\n");
    message("Starting point (Patches: %i; Individuals(F/M): [%i/%i; %i/%i])\n",
            get_nbFullPatch(),
            size(FEM, OFFSx),  size(MAL, OFFSx),
            size(FEM, ADLTx),   size(MAL, ADLTx));
    assert(_popPtr->individual_container_ok());
#endif
    
    int nbOutput = _generations/10;
    if(!nbOutput) nbOutput = 1;   // if less than 10 generations are simulated
    
    
    // ------------------------------ CYCLE --------------------------------
    for(_currentGeneration = 1; _currentGeneration <= _generations; _currentGeneration++) {
        executeBeforeEachGeneration(_currentGeneration);
        
   //     if(_currentGeneration>1) temporal_change();
        temporal_change();
        // ouput to the console
        if(_currentGeneration<=10 || !(_currentGeneration % nbOutput)){
#ifdef _SHOW_MEMORY
            message("\r    replicate %i/%i [%s] %i/%i   RAM: %f MB",getCurrentReplicate()+1,getReplicates()
                    ,getElapsedTime(time(0), startTime).c_str(),_currentGeneration,_generations, process_mem_usage());
#else
            message("\r    replicate %i/%i [%s] %i/%i",getCurrentReplicate()+1,getReplicates()
                    ,getElapsedTime(time(0), startTime).c_str(),_currentGeneration,_generations);
#endif
            fflush(stdout);
        }
        
#ifdef _DEBUG
        message("\n%i. generation (%i. replicate) ***********\n", _currentGeneration, getCurrentReplicate()+1);
#endif
        
        // loop through the life cycle events
        for(LCE = _theCycle.begin(); LCE != _theCycle.end(); ++LCE) {
            (*LCE)->execute();
            _currentAge ^= (*LCE)->removeAgeClass();
            _currentAge |= (*LCE)->addAgeClass();
            
            // check if at least a pioopualtion still exists
            if(!isAlive()){
                message("\r    replicate %i/%i [%s] %i/%i -> Pop extinction !\n",getCurrentReplicate()+1,getReplicates()
                        ,getElapsedTime(time(0), startTime).c_str(),_currentGeneration,_generations);
                _currentGeneration++;
                break;
            }
        }
        
        // check if it was stopped ("isAlive()")
        if(LCE != _theCycle.end()) break;
        
        executeAfterEachGeneration(_currentGeneration);
    }
    // --------------------------- END OF CYCLE --------------------------
}

// ----------------------------------------------------------------------------------------
// Loop_generation
// ----------------------------------------------------------------------------------------
/* performs the generation loop for current replicate */
void TMetapop::Loop_generation_coal(time_t startTime)
{
    vector< LCE* >::iterator LCE;
    
#ifdef _DEBUG
    message("\n");
    message("Starting point (Patches: %i; Individuals(F/M): [%i/%i; %i/%i])\n",
            get_nbFullPatch(),
            size(FEM, OFFSx),  size(MAL, OFFSx),
            size(FEM, ADLTx),   size(MAL, ADLTx));
#endif
    
    int nbOutput = _generations/10;
    if(!nbOutput) nbOutput = 1;   // if less than 10 generations are simulated
    
    
    // ------------------------------ CYCLE --------------------------------
    for(_currentGeneration = 1; _currentGeneration <= _generations; _currentGeneration++) {
        executeBeforeEachGeneration(_currentGeneration);
        
        if(_currentGeneration) temporal_change();
        
        // ouput to the console
        if(_currentGeneration<=10 || !(_currentGeneration % nbOutput)){
#ifdef _SHOW_MEMORY
            message("\r    replicate %i/%i [%s] %i/%i   RAM: %f MB",getCurrentReplicate()+1,getReplicates()
                    ,getElapsedTime(time(0), startTime).c_str(),_currentGeneration,_generations, process_mem_usage());
#else
            message("\r    replicate %i/%i [%s] %i/%i",getCurrentReplicate()+1, getReplicates()
                    ,getElapsedTime(time(0), startTime).c_str(),_currentGeneration,_generations);
#endif
            fflush(stdout);
        }
        
#ifdef _DEBUG
        message("\n%i. generation (%i. replicate) ***********\n", _currentGeneration, getCurrentReplicate()+1);
#endif
        
        // loop through the life cycle events
        bool is_alive = isAlive();
        for(LCE = _theCycle.begin(); LCE != _theCycle.end() && is_alive; ++LCE) {
            //cout << endl << (*LCE)->get_event_name();     // DEBUGG
            (*LCE)->execute();
            _currentAge ^= (*LCE)->removeAgeClass();
            _currentAge |= (*LCE)->addAgeClass();
        }
        
        // if all pos are empty
        if(!is_alive) {
            message("\r    replicate %i/%i [%s] %i/%i -> Pop extinction !\n",getCurrentReplicate()+1,getReplicates()
                    ,getElapsedTime(time(0), startTime).c_str(),_currentGeneration,_generations);
            _currentGeneration++;
            break;
        }
        
        executeAfterEachGeneration(_currentGeneration);
    }
    // --------------------------- END OF CYCLE --------------------------
}

//------------------------------------------------------------------------------
// executeBeforeEachGeneration
//------------------------------------------------------------------------------
void
TMetapop::executeBeforeEachGeneration(const unsigned int& gen)
{
    // density dependent temporal change of the dispersal rate
    change_disp_rate_after_density(gen);
    
    // metapop
    //temporal_change(gen);
    
    // proto traits
    map< string,TTraitProto* >::iterator iter = getTraitPrototypes().begin();
    for(; iter != getTraitPrototypes().end(); ++iter){
        iter->second->executeBeforeEachGeneration(gen);
    }
    
    // LCEs
    vector<LCE*>::iterator iterLCE = _theCycle.begin();
    for(; iterLCE != _theCycle.end(); ++iterLCE){
        (*iterLCE)->executeBeforeEachGeneration(gen);
    }
    
    if(_pSelection) _pSelection->executeBeforeEachGeneration(gen);
    
    // genetic map
    _protoGenome->executeBeforeEachGeneration(gen);
    
}

//------------------------------------------------------------------------------
// executeAfterEachGeneration
//------------------------------------------------------------------------------
void
TMetapop::executeAfterEachGeneration(const unsigned int& gen)
{
    // proto traits
    map< string,TTraitProto* >::iterator iter = getTraitPrototypes().begin();
    for(; iter != getTraitPrototypes().end(); ++iter){
        iter->second->executeAfterEachGeneration(gen);
    }
    
    // LCEs
    vector<LCE*>::iterator iterLCE = _theCycle.begin();
    for(; iterLCE != _theCycle.end(); ++iterLCE){
        (*iterLCE)->executeAfterEachGeneration(gen);
    }
}

//------------------------------------------------------------------------------
// executeBeforeEachReplicate
//------------------------------------------------------------------------------
void
TMetapop::executeBeforeEachReplicate(const unsigned int& rep)
{
    // genetic map
    _protoGenome->executeBeforeEachReplicate(rep);
    
    // traits
    map< string,TTraitProto* >::iterator iterTrait = getTraitPrototypes().begin();
    for(; iterTrait != getTraitPrototypes().end(); ++iterTrait){
        iterTrait->second->executeBeforeEachReplicate(rep);
    }
    
    // selection
    if(_pSelection) _pSelection->executeBeforeEachReplicate(rep);
    if(_total_carrying_capacity != my_NAN){
        set_total_carrying_capacity();   // soft selection at the metapopulation level
    }
    
    // LCEs
    vector<LCE*>::iterator iterLCE = _theCycle.begin();
    for(; iterLCE != _theCycle.end(); ++iterLCE){
        (*iterLCE)->executeBeforeEachReplicate(rep);
    }
}

//------------------------------------------------------------------------------
// executeBeforeAfterReplicate
//------------------------------------------------------------------------------
void
TMetapop::executeAfterEachReplicate(const unsigned int& rep)
{
    // traits
    map< string,TTraitProto* >::iterator iter = getTraitPrototypes().begin();
    for(; iter != getTraitPrototypes().end(); ++iter){
        iter->second->executeAfterEachReplicate(rep);
    }
    
    // LCEs
    vector<LCE*>::iterator iterLCE = _theCycle.begin();
    for(; iterLCE != _theCycle.end(); ++iterLCE){
        (*iterLCE)->executeAfterEachReplicate(rep);
    }
}

// ----------------------------------------------------------------------------------------
// densityDependend_dispRate_change
// ----------------------------------------------------------------------------------------
/** this function allows to change the dispersal rate depending on the density
 * i.e. at the specified ratio the temoral parameter is changed
 * _density_threshold[n][3] = 0: patch, 1: density, 2: change   //n is the number of changes
 */
void
TMetapop::change_disp_rate_after_density(const int& gen)
{
    if(!_density_threshold) return;
    
    double density;
    TPatch* curPatch;
    map<int, string>::iterator iter_tempParam;
    for(unsigned int i=0; i<_density_threshold_nbChanges; ++i){
        if(_density_threshold[i][0] == my_NAN) continue;		// if already used
        curPatch = get_vPatch((unsigned int)_density_threshold[i][0]);
        if(!curPatch->get_K()) continue;										// the patch has to have a carring capacity
        density = (double)curPatch->size(ADLTx) / curPatch->get_K();
        if(density >= _density_threshold[i][1]){            // if a change has to be made
            _density_threshold[i][0] = my_NAN;              	// set to NAN to recognise used items
            
            // get the x item (it is a map so scroll through it)
            iter_tempParam = _density_threshold_param[i]->get_temporal_args()->begin();
            for(int j=0; j<_density_threshold[i][3]; ++j, ++iter_tempParam){}
            string arg  = iter_tempParam->second;
            int cur_gen = iter_tempParam->first;
            int new_gen = gen + (int)_density_threshold[i][2];
            (*_density_threshold_param[i]->get_temporal_args())[new_gen] = arg; // make the change
            _density_threshold_param[i]->get_temporal_args()->erase(iter_tempParam); // delete the old item
            
            // also the paramset has to be changed...
            multimap<int, map<string, Param*>* >* m = _density_threshold_param[i]->getParamSet()->getTemporalParams();
            multimap<int, map<string, Param*>* >::iterator pos = m->find(cur_gen);
            assert(pos != m->end());
            delete pos->second;
            m->erase(pos);               // delete the item
            _density_threshold_param[i]->getParamSet()->set_temporal_param(new_gen, _density_threshold_param[i]);  // add the new one
        }
    }
}

// ----------------------------------------------------------------------------------------
/** this function allows to change the dispersal rate depending on the density
 * i.e. at the specified ratio the temoral parameter is changed
 * _density_threshold[n][4] = 0: patch, 1: density, 2: delay, 3: change   //n is the number of changes
 */
void
TMetapop::set_change_disp_rate_after_density()
{
    Param* p = _paramSet->get_param("temporal_change_following_density");
    if(_density_threshold){delete[] _density_threshold; _density_threshold=NULL;}
    if(!p->isSet())return;
    
    // if param is specified
    if(!p->is_matrix()) error("Parameter 'temporal_change_following_density' must be a matrix!");
    TMatrixVar<string>* m = p->get_matrixVarStr();
    
    _density_threshold_nbChanges = m->getNbRows();
    Param*    curParam;
    string paramName;
    vector<string>* curVec;
    _density_threshold  = ARRAY::new_2D<double>(_density_threshold_nbChanges, 4);
    _density_threshold_param = new Param*[_density_threshold_nbChanges];
    for(unsigned int i=0; i<_density_threshold_nbChanges; ++i){
        curVec = m->get(i);
        if(curVec->size() != 5) error("Parameter 'temporal_change_following_density' must be a matrix with 5 elements per row!");
        paramName = (*curVec)[0];
        
        // check the metapop params
        curParam = _paramSet->find_param(paramName);
        
        // check the LCE params
        if(!curParam){
            vector<LCE*>::iterator  iterLCE = _theCycle.begin();
            for(;iterLCE != _theCycle.end(); ++iterLCE) {
                curParam = (*iterLCE)->get_paramset()->find_param(paramName);
                if(curParam) break;
            }
        }
        
        // check the trait params
        if(!curParam){
            map< string,TTraitProto* >::iterator iter = getTraitPrototypes().begin();
            for(; iter != getTraitPrototypes().end(); ++iter){
                curParam = iter->second->get_paramset()->find_param(paramName);
                if(curParam) break;
            }
        }
        
        // was the param found?
        if(curParam) _density_threshold_param[i] = curParam;
        else error("Parameter 'temporal_change_following_density' contains unknown parameter name '%s'!", paramName.c_str());
        
        // get the dynamic change settings
        _density_threshold[i][0] = (double)(strTo<int>((*curVec)[1])-1); // patch (here starting at 0)
        _density_threshold[i][1] = strTo<double>((*curVec)[2]);          // density
        _density_threshold[i][2] = (double)(strTo<int>((*curVec)[3]));   // delay
        _density_threshold[i][3] = (double)(strTo<int>((*curVec)[4])-1); // change (here starting at 0)
        
        // some tests
        if(_density_threshold[i][0]>=0 && _density_threshold[i][0]<_patchNbr) error("Parameter 'temporal_change_following_density' must have a valide patch id!");
        if(_density_threshold[i][2]<0) error("Parameter 'temporal_change_following_density' must have a positive delay!");
        if(!curParam->isTemporalParam() || curParam->get_temporal_args()->size()<_density_threshold[i][3])
            error("Parameter 'temporal_change_following_density' should at least contain %i temporal changes!", _density_threshold[i][3]);
    }
    
    delete m;
}

// ----------------------------------------------------------------------------------------
// get_total_carrying_capacity_bothSex
// ----------------------------------------------------------------------------------------
unsigned int
TMetapop::get_total_carrying_capacity_bothSex(){
    unsigned int size = 0;
    for(unsigned int home = 0; home < _patchNbr; ++home) {
        size += _vPatch[home]->get_K();
    }
    return size;
}

// ----------------------------------------------------------------------------------------
// get_settingsStats
// ----------------------------------------------------------------------------------------
void
TMetapop::get_settingsStats(unsigned int& nbPatchK,unsigned int& K,
                           unsigned int& NiniPatch, unsigned int& Nini,
                           unsigned int& samplePatch, unsigned int& nbSample)
{
    nbPatchK = K = NiniPatch = Nini = samplePatch = nbSample =0;
    TPatch* curPatch;
    double curSample;
    vector<TPatch*>::iterator curPop, endPop=_vPatch.end();
    for(curPop=_vPatch.begin(); curPop!=endPop; ++curPop){
        curPatch = *curPop;
        if(curPatch->get_K()){
            ++nbPatchK;                             // number of patches with K
            K += curPatch->get_K();                 // number of individuals at K
        }
        if(curPatch->get_N_ini()){
            ++NiniPatch;                            // number of initial populations
            Nini += curPatch->get_N_ini();          // number of initial individuals
        }
        if(curPatch->get_sampleID()!=my_NAN){
            ++samplePatch;                          // number of sampled populations
            curSample = curPatch->get_N_sample();
            if(curSample==my_NAN)  nbSample += curPatch->get_K();                                    // no specific sampling
            else if(curSample<1.0) nbSample += (unsigned int)my_round(curSample*curPatch->get_K());  // relative sampling
            else                   nbSample += (unsigned int)curSample;                              // absolute sampling
        }
    }
}

// ----------------------------------------------------------------------------------------
// set_total_carrying_capacity
// ----------------------------------------------------------------------------------------
unsigned int
TMetapop::get_total_carrying_capacity(const sex_t& SEX){
    unsigned int sum = 0;
    for(unsigned int home = 0; home < _patchNbr; ++home) {
        sum += _vPatch[home]->get_K(SEX);
    }
    return sum;
}

// ----------------------------------------------------------------------------------------
// get_total_iniSize_bothSex
// ----------------------------------------------------------------------------------------
unsigned int
TMetapop::get_total_iniSize_bothSex(){
    unsigned int size = 0;
    for(unsigned int home = 0; home < _patchNbr; ++home) {
        size += _vPatch[home]->get_N_ini();
    }
    return size;
}

//------------------------------------------------------------------------------
/** Metapop::reset_metapopulation */
//------------------------------------------------------------------------------
/* remove all individuals of each pop and reset the counters */
void
TMetapop::reset_metapopulation()
{
    
    for(unsigned int i = 0; i < _patchNbr; ++i) {
        _vPatch[i]->reset_counters();
        _vPatch[i]->flush();
        _vPatch[i]->set_isExtinct(true);
    }
    
    setSampleSizes(true);
    
    //flush();            // recycle the remaining individuals
    _vFullPatch.clear();
    _vSamplePatch.clear();
    
}

//------------------------------------------------------------------------------
/** Metapop::get_requiredAgeClass */
//------------------------------------------------------------------------------
/* get the required age class to start a simulation */
age_t
TMetapop::get_requiredAgeClass()
{
    age_t requiredAge = NONE;
    vector< LCE* >::iterator IT = _theCycle.begin();
    for(; IT != _theCycle.end() && requiredAge == NONE; IT++){
        requiredAge = (*IT)->requiredAgeClass();
    }
    if(requiredAge == NONE) error("The life cycle makes no sense: no age is required!\n");
    return requiredAge;
}

//------------------------------------------------------------------------------
// allPatchPopulated_atInitialization
//------------------------------------------------------------------------------
/** checks if all patches are populated at the start of the simulation
 * returns true if it is the case and false if not
 * Caution: currently it is not checked if over time carrying capacity is changing
 */
bool
TMetapop::allPatchPopulated_atInitialization()
{
    vector<TPatch*>::iterator curPop, endPop = _vPatch.end();
    for(curPop=_vPatch.begin(); curPop!=endPop; ++curPop){
        if(!(*curPop)->get_K()) return false;      // if the patch cannot be colonized
        if(!(*curPop)->get_N_ini()) return false;
    }
    return true;
}

//------------------------------------------------------------------------------
// setPopulation
//------------------------------------------------------------------------------
/** Metapop::setPopulation */
void
TMetapop::setPopulation(const age_t& requiredAge)
{
#ifdef _DEBUG
    message("Metapop::setPopulation: required age is %i ",requiredAge);
#endif
    
    for(unsigned int i = 0; i < _patchNbr; ++i){
        if(_vPatch[i]->setNewGeneration(requiredAge)) new_fullPatch(_vPatch[i]);
    }
    add_tempPatch();
    _currentAge = requiredAge;
    
#ifdef _DEBUG
    message("(Patches: %i; Individuals(F/M): [%i/%i; %i/%i])\n",
            get_nbFullPatch(),
            size(FEM, OFFSx),  size(MAL, OFFSx),
            size(FEM, ADLTx),  size(MAL, ADLTx));
    assert(_popPtr->individual_container_ok());
#endif
}

//------------------------------------------------------------------------------
// setPopulation_coal
//------------------------------------------------------------------------------
/** Metapop::setPopulation_coal */
void
TMetapop::setPopulation_coal(const age_t& requiredAge)
{
#ifdef _DEBUG
    message("Metapop::setPopulationCoal: required age is %i ",requiredAge);
#endif
    
    for(unsigned int i = 0; i < _patchNbr; ++i){
        if(_vPatch[i]->setNewGeneration_coal(requiredAge)) new_fullPatch(_vPatch[i]);
    }
    add_tempPatch();
    _currentAge = requiredAge;
    
#ifdef _DEBUG
    message("(Patches: %i; Individuals(F/M): [%i/%i; %i/%i])\n",
            get_nbFullPatch(),
            size(FEM, OFFSx),  size(MAL, OFFSx),
            size(FEM, ADLTx),   size(MAL, ADLTx));
#endif
}

//------------------------------------------------------------------------------
// setPopulation_FSTAT
//------------------------------------------------------------------------------
/** Metapop::setPopulation_FSTAT
 * returns true if the pop was set by genotypes and false if not
 * format of the fstat file as array:
 *       v[ind][0][0] = patch             // current patch (compulsory, index starting at 0)
 *       v[ind][0][1] = NaN               // currently not used
 *       v[ind][1][0] = nb_loci           // number of locus found in the file (compulsory)
 *       v[ind][1][1] = max_allele        // maximal index of the alleles (compulsory)
 *       v[ind][2][0] = age               // off: 0, adlt: 1
 *       v[ind][2][1] = sex               // mal: 0, fem: 1
 *       v[ind][3][0] = id                // id of the individual
 *       v[ind][3][1] = natal patch       // id of the individual's natal patch (index starting at 0)
 *       v[ind][4][0] = id_mom            // id of mother
 *       v[ind][4][1] = mom's patch       // id of mother's natal patch (index starting at 0)
 *       v[ind][5][0] = id_dad            // id of father
 *       v[ind][5][1] = dad's patch       // id of father's natal patch (index starting at 0)
 *       v[ind][loc+6][all] (compulsory)  // allele starting at 0
 */
bool
TMetapop::setPopulation_FSTAT(const age_t& requiredAge)
{
    // read all initial genotype files (a single file for each trait type)
    vector<unsigned int**> *vAny, *vAny2;
    string type;
    unsigned int t;
    TTraitProto* pProto;
    map<string, vector<unsigned int**>* > mFstat;         // for each Fstat file an entry
    for(t=0; t<getTraitPrototypeSize(); ++t){             // scroll through each trait
        pProto = &getTraitPrototype(t);                     // get the trait
        if(pProto->get_ini_fstat_file().empty()) continue;  // if not specified continue
        if(pProto->get_trait_index()>1) continue;           // read the file only once
        type = getTraitPrototype(t).get_type();             // get the type
        if(mFstat.find(type) != mFstat.end()) continue; 		// if alrerady read, ignore it
        mFstat[type] = read_Fstat_file(get_pSimulation()->get_fileName(pProto->get_ini_fstat_file())); // read the Fstat file
    }
    if(mFstat.empty()) return false;                      // the metapop was not initialized by Fstat files
    
    /////// the metapop will be initialized by the genotypes of the Fstat file /////////
#ifdef _DEBUG
    message("Metapop::setPopulationFSTAT: required age is %i ",requiredAge);
#endif
    
    // check if the number of loci of the Fstat file correponds to the number of loci specified in the ini-file
    // and find the correct Fstat file for each trait
    map<string, unsigned int> nbLoci;
    map<string, vector<unsigned int**>* >::iterator pos, end;
    end = mFstat.end();
    vector<unsigned int**>** trait_Fstat = new vector<unsigned int**>*[getTraitPrototypeSize()];
    for(t=0; t<getTraitPrototypeSize(); ++t){             // get the number of loci of each type of trait
        type = getTraitPrototype(t).get_type();             // get the type of trait
        nbLoci[type] += getTraitPrototype(t).get_nb_locus();// sum up the number of loci
        pos = mFstat.find(type);                            // check if a file is available for this trait
        trait_Fstat[t] = (pos==end) ? NULL : pos->second;   // set the corresponding pointer if available
    }
    for(pos = mFstat.begin(); pos != end; ++pos){ 	      // check each Fstat file
        if((*pos->second)[0][1][0] != nbLoci[pos->first])
            error("initial FSTAT file (%s): the number of loci (%i) does not correspond to the settings (%i)!\n",
                  pos->first.c_str(), (*pos->second)[0][1][0], nbLoci[pos->first]);
    }
    
    // if there are more than one Fstat file check if the files have the same basic information and get the file with most of the info (vAny)
    vAny = mFstat.begin()->second;                        // get the first file
    unsigned int nbInd = (unsigned int)vAny->size();                    // get the number of individuals
    for(pos = mFstat.begin(), ++pos; pos != end; ++pos){ 	// check if each other file has the same info
        vAny2 = pos->second;
        if(nbInd != vAny2->size()) error("initial FSTAT files: the number of individuals does not correspond between files!\n");
        if((*vAny)[0][2][0] != my_NAN){                     // check if the file has the supplement info (age checked)
            if((*vAny2)[0][2][0] != my_NAN){                  // check if the file has the supplement info (age checked)
                // check each individual for equal suplement information (both are present)
                for(unsigned int i=0, size=(unsigned int)vAny->size(); i<size; ++i){
                    if(  (*vAny)[i][0][0] != (*vAny2)[i][0][0]    // current patch: starting at 0 (compulsory)
                       || (*vAny)[i][2][0] != (*vAny2)[i][2][0]    // age          // off: 0, adlt: 1
                       || (*vAny)[i][2][1] != (*vAny2)[i][2][1]    // sex          // mal: 0, fem: 1
                       || (*vAny)[i][3][0] != (*vAny2)[i][3][0]    // id           // id of the current patch starting at 1
                       || (*vAny)[i][3][1] != (*vAny2)[i][3][1]    // natal patch  // id of the natal patch starting at 1
                       || (*vAny)[i][4][0] != (*vAny2)[i][4][0]    // id mom       // id of mother of the current patch starting at 1
                       || (*vAny)[i][4][1] != (*vAny2)[i][4][1]    // patch mom    // id of mother of the natal patch starting at 1
                       || (*vAny)[i][5][0] != (*vAny2)[i][5][0]    // id dad       // id of father of the current patch starting at 1
                       || (*vAny)[i][5][1] != (*vAny2)[i][5][1]){  // patch dad    // id of father of the natal patch starting at 1
                        error("initial FSTAT files: the charactersitics of the individuals do not match between genotype files!\n");
                    }
                }
            }
        }
        else if((*vAny2)[0][2][0] != my_NAN) vAny = vAny2;         // only 2. has info: change the 'any'
    }
    
    // for each individual
    TIndividual* cur_ind;
    age_idx age;
    sex_t sex;
    TPatch* pPatch;
    unsigned int i, l, p, curLocus, nbAllele;
    unsigned int** curGenotype;
    unsigned int** array;
    unsigned char** seq;
    unsigned char a;
    unsigned int *aMax = ARRAY::new_1D<unsigned int>(_patchNbr, (unsigned int)0);   // array for each patch to find the maximal index (only used when infos are present)
    for(i=0; i<nbInd; ++i){                               // for each individual
        array   = (*vAny)[i];                               // get the individual array which has the complete information
        
        // get the age of the individual and check if the age is desired
        if(array[2][0] != my_NAN){                          // the age is given
            age = (age_idx)array[2][0];                       // get the age
            if(  (age==OFFSx && requiredAge!=OFFSPRG)         // stop here if the individual has the wrong age
               || (age==ADLTx && requiredAge!=ADULTS)) continue;
        }
        else{                                               // the age is not given set it according to the required age
            switch(requiredAge){                              // it is assumed that all individuals have the required age
                case OFFSPRG: age = OFFSx;  break;
                case ADULTS:  age = ADLTx;  break;
            }
        }
        
        // create the individual
        cur_ind = _protoIndividual.clone();
        cur_ind->init();                                    //allocate memory for the traits sequence
        
        // set the sex
        sex = (array[2][1]!=my_NAN) ? (sex_t)array[2][1] : (rand().Uniform()<_sexInitRatio ? MAL : FEM);
        cur_ind->setSex(sex);
        
        // set the current patch and add the individeual to the patch (compulsory)
        l = array[0][0]; // get the patch index
        if(l >= _patchNbr) error("initial FSTAT files: the patch index (%i) exceeds the number of patches (%i)!\n", l+1, _patchNbr);
        pPatch = get_vPatch(l);
        cur_ind->setCurrentPatch(pPatch);
        pPatch->add(sex, age, cur_ind);     	     // add the individualto the container
        
        // set the ID of the individual
        if(array[3][0] != my_NAN && array[3][1] != my_NAN){
            if(array[3][1] >= _patchNbr) error("initial FSTAT files: the patch index (%i) exceeds the number of patches (%i)!\n", l+1, _patchNbr);
            cur_ind->setID(toStr(array[3][0]) + "_" + toStr(array[3][1]));
            cur_ind->setNatalPatch(get_vPatch(array[3][1]));    // set natal patch
            if(aMax[array[3][1]] < array[3][0]) aMax[array[3][1]] = array[3][0]; // get the maximal id of an individual
        }
        else{                                               // the ID is not given
            cur_ind->setID(pPatch->get_next_IDofIndividual());
            cur_ind->setNatalPatch(pPatch);                   // set the current instead of the natal patch
        }
        
        // set the ID of the mother
        if(array[4][0] != my_NAN && array[4][1] != my_NAN){
            if(array[4][1] >= _patchNbr) error("initial FSTAT files: the patch index (%i) exceeds the number of patches (%i)!\n", l+1, _patchNbr);
            cur_ind->setMotherID(toStr(array[4][0]) + "_" + toStr(array[4][1]+1));
        }
        else cur_ind->setMotherID("");
        
        // set ID of the father
        if(array[5][0] != my_NAN && array[5][1] != my_NAN){
            if(array[5][1] >= _patchNbr) error("initial FSTAT files: the patch index (%i) exceeds the number of patches (%i)!\n", l+1, _patchNbr);
            cur_ind->setFatherID(toStr(array[5][0]) + "_" + toStr(array[5][1]+1));
        }
        else cur_ind->setFatherID("");
        
        // set the genotype of each trait (if not given by the Fstat file allocate it randomly)
        for(t=0; t<getTraitPrototypeSize(); ++t){           // scroll through each trait
            if(!trait_Fstat[t]){                              // genotypes are not given: initialize them randomly
                cur_ind->Traits[t]->ini_sequence(pPatch);
            }
            else{                                             // gentoypes are given by the input file
                pProto = &getTraitPrototype(t);                 // get the trait
                curGenotype = (*trait_Fstat[t])[i];             // genotype is given by the Fstat file
                seq = cur_ind->Traits[t]->sequence;
                if(pProto->get_trait_index()<2) curLocus=6;     // a new type of trait starts: reset the index
                for(l = 0; l < pProto->get_nb_locus(); ++l, ++curLocus){ // for each locus
                    nbAllele = pProto->get_nb_allele(l);          // get the number of alleles
                    for(p = 0; p < ploidy; ++p){    							// for each allele
                        a = curGenotype[curLocus][p];
                        if((unsigned int)a >= nbAllele)
                            error("initial FSTAT files: the allele index (%u) exeeds the maximal number of alleles (%u, individual %u, locus %u, type %s)!",
                                  (unsigned int) a, nbAllele, i+1, curLocus-5, pProto->get_type().c_str());
                        seq[l][p] = a;
                    }
                }
            }
            cur_ind->Traits[t]->set_value();                  // set the genotype (for all traits)
        }
    }
    
    // clean up
    for(pos = mFstat.begin(); pos != end; ++pos){
        for(i=0; i<pos->second->size(); ++i){
            curGenotype = (*pos->second)[i];
            ARRAY::delete_2D(curGenotype, curGenotype[1][0]+6);
        }
        delete pos->second;
    }
    
    delete[] trait_Fstat;
    
    // adjust the individual indexer of the patch to the highest present index
    for(t=0; t<_patchNbr; ++t){                                 // for each patch
        if(aMax[t]) pPatch->set_ID_individual(++aMax[t]);
    }
    delete[] aMax;
    
    // get the populated patches
    for(i = 0; i < _patchNbr; ++i){
        if(_vPatch[i]->size(ALL)) new_fullPatch(_vPatch[i]);
    }
    add_tempPatch();
    
    
    _currentAge = requiredAge;
    
#ifdef _DEBUG
    message("(Patches: %i; Individuals(F/M): [%i/%i; %i/%i])\n",
            get_nbFullPatch(),
            size(FEM, OFFSx),  size(MAL, OFFSx),
            size(FEM, ADLTx),   size(MAL, ADLTx));
    assert(_popPtr->individual_container_ok());
#endif
    
    return true;
}

// ----------------------------------------------------------------------------------------
// reset
// ----------------------------------------------------------------------------------------
/* Metapop::reset
 * resets each Patch, all individuals are moved to the POOL
 */
void TMetapop::reset()
{
    for(unsigned int i = 0; i < _patchNbr; ++i) {
        _vPatch[i]->flush();
        _vPatch[i]->reset_ID_individual();
    }
}

// ----------------------------------------------------------------------------------------
// clear
// ----------------------------------------------------------------------------------------
void TMetapop::clear()
{
   // cout << endl << "delete patches" << endl;
    vector<TPatch*>::iterator curPop, endPop;
    for(curPop=_vPatch.begin(), endPop=_vPatch.end(); curPop!=endPop; ++curPop){
        delete *curPop;
    }
    _vPatch.clear();
    
   // cout << endl << "delete selection" << endl;
    if(_pSelection)        {delete _pSelection;          _pSelection=NULL;}
    if(_density_threshold_param) {delete[] _density_threshold_param; _density_threshold_param=NULL;}
    if(_density_threshold) ARRAY::delete_2D(_density_threshold, _density_threshold_nbChanges);
   // cout << endl << "delete patches II" << endl;
}

// ----------------------------------------------------------------------------------------
// getGenerationCounter
// ----------------------------------------------------------------------------------------
string
TMetapop::getGenerationCounter (){
    return toStr(_currentGeneration, _generations);
}
// ----------------------------------------------------------------------------------------
string
TMetapop::getGenerationCounter_ (){
    return "_" + getGenerationCounter();
}
// ----------------------------------------------------------------------------------------
string
TMetapop::getGenerationCounter_g (){
    return "_g" + getGenerationCounter();
}
// ----------------------------------------------------------------------------------------
/** soft selection (at the patch level) to K
 * the order of the Kmal and Kfem arrays are only for the populated pops
 */
void
TMetapop::regulate_selection_fitness_patch(age_idx AGE, unsigned int* Kmal, unsigned int* Kfem)
{
    assert(_pSelection);
    if(Kmal && Kfem){       // the size to regulate to is given
        vector<TPatch*>::iterator curPop, endPop;
        for(curPop=get_vFullPatch().begin(), endPop=get_vFullPatch().end(); curPop!=endPop; ++curPop, ++Kmal, ++Kfem) {
            _pSelection->set_fitness(*curPop, AGE); // compute the fitnesses, but don't sort or make yet the array cumulative
            (*curPop)->regulate_selection_fitness(*Kfem, _pSelection, FEM, AGE);
            (*curPop)->regulate_selection_fitness(*Kmal, _pSelection, MAL, AGE);
        }
    }
    else{                 	// the pops are regulated to carrying capacity
        vector<TPatch*>::iterator curPop, endPop;
        for(curPop=get_vFullPatch().begin(), endPop=get_vFullPatch().end(); curPop!=endPop; ++curPop) {
            _pSelection->set_fitness(*curPop, AGE); // compute the fitnesses, but don't sort or make yet the array cumulative
            (*curPop)->regulate_selection_fitness((*curPop)->get_KFem(), _pSelection, FEM, AGE);
            (*curPop)->regulate_selection_fitness((*curPop)->get_KMal(), _pSelection, MAL, AGE);
        }
    }
}

// ----------------------------------------------------------------------------------------
/** soft selection (at the metapopulation level) to K */
void
TMetapop::regulate_selection_fitness_metapop(age_idx AGE)
{
    // array to buffer the fitness arrays for each patch (female and male)
    double totFitnessM = 0;                              // total fitness
    double totFitnessF = 0;                              // total fitness
    double* sumFitnessM = new double[get_nbFullPatch()];         // mean fitness of each patch
    double* sumFitnessF = new double[get_nbFullPatch()];         // mean fitness of each patch
    unsigned nbInd[2];
    unsigned i, Kmal_tot=_popPtr->get_total_carrying_capacity(MAL),
    Kfem_tot=_popPtr->get_total_carrying_capacity(FEM);
    TPatchFitness** fitnessArrays[2];
    for(i=0; i<2; ++i){
        fitnessArrays[i] = new TPatchFitness*[get_nbFullPatch()];
        nbInd[i] = 0;
    }
    
    // for each patch compute the fitness of the individuals
    vector<TPatch*>::iterator curPop, endPop=get_vFullPatch().end();
    for(i=0, curPop=get_vFullPatch().begin(); curPop!=endPop; ++curPop, ++i) {
        // compute and store the fitness of the current patch (0: male, 1: female)
        _pSelection->set_fitness(*curPop, AGE);
        sumFitnessM[i]   = _pSelection->getSumFitness(MAL);
        totFitnessM      += sumFitnessM[i];
        sumFitnessF[i]   = _pSelection->getSumFitness(FEM);
        totFitnessF      += sumFitnessF[i];
        if((*curPop)->size(MAL, AGE)) fitnessArrays[MAL][i] = _pSelection->get_FitnessObject(MAL, true);
        else                          fitnessArrays[MAL][i] = NULL;
        if((*curPop)->size(FEM, AGE)) fitnessArrays[FEM][i] = _pSelection->get_FitnessObject(FEM, true);
        else                          fitnessArrays[FEM][i] = NULL;
    }
    
    // regulate each patch separately
    for(i=0, curPop=get_vFullPatch().begin(); curPop!=endPop; ++curPop, ++i) {
        // regulate the males
        if(fitnessArrays[MAL][i]){
            _pSelection->set_FitnessObject(MAL, fitnessArrays[MAL][i]);    // set the fitness arrays of the current patch
            (*curPop)->regulate_selection_fitness(my_round(sumFitnessM[i]/totFitnessM*Kmal_tot),
                                                  _pSelection, MAL, AGE);
        }
        
        // regulate the females
        if(fitnessArrays[FEM][i]){
            _pSelection->set_FitnessObject(FEM, fitnessArrays[FEM][i]);    // set the fitness arrays of the current patch
            (*curPop)->regulate_selection_fitness(my_round(sumFitnessF[i]/totFitnessF*Kfem_tot),
                                                  _pSelection, FEM, AGE);
        }
    }//end_for_nbPatch
    
    // delete the fitness arrays
    if(sumFitnessF)      delete[] sumFitnessF;
    if(sumFitnessM)      delete[] sumFitnessM;
    for(i=0; i<2; ++i){
        if(fitnessArrays[i]) delete[] fitnessArrays[i];
    }
}

// ----------------------------------------------------------------------------------------
/** hard selection to N (not to K!),
 * thus a indivdual with a fitness of 0.5 has a survival probability of 0.5
 */
void
TMetapop::regulate_selection_fitness_hard(age_idx AGE)
{
    vector<TPatch*>::iterator curPop, endPop;
    for(curPop=get_vFullPatch().begin(), endPop=get_vFullPatch().end(); curPop!=endPop; ++curPop) {
        (*curPop)->regulation_selection_hard(AGE);
    }
}

// ----------------------------------------------------------------------------------------
/** adds the _vTempPatch container to the _vFullPatch container
 * samleID is in this case the index of the patch in _vFullPatch!!!
 */
void
TMetapop::add_tempPatch_noSample_withFull()
{
    if(_vTempPatch.empty()) return;
    
    // set the sample id
    vector<TPatch*>::iterator curPop=_vTempPatch.begin(), endPop=_vTempPatch.end();
    for(unsigned int i=(unsigned int)_vFullPatch.size(); curPop!=endPop; ++curPop){
        //if((*curPop)->get_sampleID()==my_NAN) continue; // should never happen
        assert((*curPop)->get_sampleID()==SAMPLED);
        (*curPop)->set_sampleID(i++);
    }
    
    // merge the containers
    _vFullPatch.insert(_vFullPatch.end(), _vTempPatch.begin(), _vTempPatch.end());
    _vTempPatch.clear();
}

// ----------------------------------------------------------------------------------------
/** adds the _vTempPatch container to the _vFullPatch container
 */
void
TMetapop::add_tempPatch_noSample_withFull_coal()
{
    if(_vTempPatch.empty()) return;
    
    // merge the containers
    _vFullPatch.insert(_vFullPatch.end(), _vTempPatch.begin(), _vTempPatch.end());
    _vTempPatch.clear();
}

// ----------------------------------------------------------------------------------------
/** adds the _vTempPatch container to the _vFullPatch container
 */
void
TMetapop::add_tempPatch_withSample_withFull()
{
    if(_vTempPatch.empty()) return;
    
    // set the sample id
    vector<TPatch*>::iterator curPop=_vTempPatch.begin(), endPop=_vTempPatch.end();
    for(unsigned int i=(unsigned int)_vSamplePatch.size(); curPop!=endPop; ++curPop){
        if((*curPop)->get_sampleID()==my_NAN) continue;
        assert((*curPop)->get_sampleID()==SAMPLED);
        (*curPop)->set_sampleID(i++);
        _vSamplePatch.push_back(*curPop);
    }
    
    // merge the containers
    _vFullPatch.insert(_vFullPatch.end(), _vTempPatch.begin(), _vTempPatch.end());
    _vTempPatch.clear();
}

// ----------------------------------------------------------------------------------------
/** since jsut the vector _vPatch is used the sampleIid is set to the index of that patch
 */
void
TMetapop::add_tempPatch_noSample_noFull()
{
    vector<TPatch*>::iterator curPop=_vPatch.begin(), endPop=_vPatch.end();
    for(unsigned int i=0; curPop!=endPop; ++curPop){
        assert((*curPop)->get_sampleID()==SAMPLED);
        (*curPop)->set_sampleID(i++);
    }
}

// ----------------------------------------------------------------------------------------
/** since jsut the vector _vPatch is used the sampleIid is set to the index of that patch
 */
void
TMetapop::add_tempPatch_noSample_noFull_coal()
{
    return;
}

// ----------------------------------------------------------------------------------------
/** adds the _vTempPatch container to the _vFullPatch container
 * only called at initizialisation of the pops
 */
void
TMetapop::add_tempPatch_withSample_noFull()
{
    if(_vTempPatch.empty()) return;
    
    // set the sample id
    vector<TPatch*>::iterator curPop=_vTempPatch.begin(), endPop=_vTempPatch.end();
    for(unsigned int i=(unsigned int)_vSamplePatch.size(); curPop!=endPop; ++curPop){
        if((*curPop)->get_sampleID()==my_NAN) continue;
        assert((*curPop)->get_sampleID()==SAMPLED);
        (*curPop)->set_sampleID(i++);
        _vSamplePatch.push_back(*curPop);
    }
    
    // merge the containers
    _vFullPatch.insert(_vFullPatch.end(), _vTempPatch.begin(), _vTempPatch.end());
    _vTempPatch.clear();
}

// ----------------------------------------------------------------------------------------
/** newly populated patch: without _vSamplePatch and with _vFullPatch */
void
TMetapop::new_fullPatch_noSample_withFull(TPatch* curPatch)
{
    assert(curPatch->size(ALL));
    curPatch->set_isExtinct(false);
    _vTempPatch.push_back(curPatch);       // add patch to the populated patches
}

// ----------------------------------------------------------------------------------------
/** newly populated patch:  with _vSamplePatch and _vFullPatch */
void
TMetapop::new_fullPatch_withSample_withFull(TPatch* curPatch)
{
    assert(curPatch->size(ALL));
    curPatch->set_isExtinct(false);
    _vTempPatch.push_back(curPatch);       // add patch to the populated patches
}

// ----------------------------------------------------------------------------------------
/** newly populated patch: without _vSamplePatch and without _vFullPatch */
void
TMetapop::new_fullPatch_noSample_noFull(TPatch* curPatch)
{
    assert(curPatch->size(ALL));
    curPatch->set_isExtinct(false);
}

// ----------------------------------------------------------------------------------------
/** newly populated patch:  with _vSamplePatch and _vFullPatch
 * only called when pops are initialized
 */
void
TMetapop::new_fullPatch_withSample_noFull(TPatch* curPatch)
{
    assert(curPatch->size(ALL));
    curPatch->set_isExtinct(false);
    _vTempPatch.push_back(curPatch);       // add patch to the populated patches
}

// ----------------------------------------------------------------------------------------
/** newly freed patch:  without _vSamplePatch and without _vFullPatch
 * all sampled and always populated
 * curPop is addapted to point to the next element
 * endPop remains the same
 */
void
TMetapop::new_emptyPatch_noSample_noFull(vector<TPatch*>::iterator& curPop, vector<TPatch*>::iterator& endPop)
{
    ++curPop;
    assert(endPop==get_vFullPatch().end());
}

// ----------------------------------------------------------------------------------------
/** newly freed patch: with _vSamplePatch, but without _vFullPatch
 * subsampling and always populated
 * curPop is addapted to point to the next element
 * endPop remains the same
 */
void
TMetapop::new_emptyPatch_withSample_noFull(vector<TPatch*>::iterator& curPop, vector<TPatch*>::iterator& endPop)
{
    assert(!(*curPop)->size(ALL));
    (*curPop)->set_isExtinct(true);
    if((*curPop)->get_sampleID()!=my_NAN){           // is the patch sampled?
        assert((*curPop)->get_sampleID()!=SAMPLED);
        erase_vSamplePatch(*curPop);                   // the sampleID is used to delete it
    }
    ++curPop;
    assert(endPop==get_vFullPatch().end());
}

// ----------------------------------------------------------------------------------------
/** newly freed patch: without _vSamplePatch, but with _vFullPatch
 * all sampled but not always populated
 * the patch to erase has to be swapped with the last patch in the container and
 * then the last not yet used element has to be removed
 * curPop remains (deleted element is swapped with the last element)
 * endPop is decreased by one!
 */
void
TMetapop::new_emptyPatch_noSample_withFull(vector<TPatch*>::iterator& curPop, vector<TPatch*>::iterator& endPop)
{
    assert(!(*curPop)->size(ALL));
    (*curPop)->set_isExtinct(true);
    (*curPop)->nbImmigrant = 0;
    (*curPop)->nbEmigrant = 0;
    
    unsigned int sampleID = (*curPop)->get_sampleID();
    assert(sampleID!=my_NAN);
    assert(sampleID!=SAMPLED);
    assert(get_vFullPatch()[sampleID]->get_sampleID()==sampleID);
    (*curPop)->set_sampleID(SAMPLED);
    
    if(*curPop != get_vFullPatch().back()){               // if it is not the last element, swap it
        *curPop = get_vFullPatch().back();                // copy the last element to the current location
        (*curPop)->set_sampleID(sampleID);                // set its new sample id
    }
    get_vFullPatch().pop_back();                          // remove the last element
    --endPop;
    assert(endPop==get_vFullPatch().end());
}

// ----------------------------------------------------------------------------------------
/** newly freed patch: without _vSamplePatch, but with _vFullPatch
 * all sampled but not always populated
 * curPop remains (deleted element is swapped with the last element)
 * endPop is decreased by one!
 */
void
TMetapop::new_emptyPatch_noSample_withFull_coal(vector<TPatch*>::iterator& curPop, vector<TPatch*>::iterator& endPop)
{
    assert(!(*curPop)->size(ALL));
    (*curPop)->set_isExtinct(true);
    if(*curPop!= get_vFullPatch().back()){
        *curPop = get_vFullPatch().back();                 // copy the last element to the current location
    }
    get_vFullPatch().pop_back();                          // remove the last element
    --endPop;
    assert(endPop==get_vFullPatch().end());
}

// ----------------------------------------------------------------------------------------
/** newly freed patch: with _vSamplePatch and with _vFullPatch
 * subsampled with not always populated patches
 * curPop remains (deleted element is swapped with the last element)
 * endPop is decreased by one!
 */
void
TMetapop::new_emptyPatch_withSample_withFull(vector<TPatch*>::iterator& curPop, vector<TPatch*>::iterator& endPop)
{
    assert(!(*curPop)->size(ALL));
    (*curPop)->set_isExtinct(true);
    if((*curPop)->get_sampleID()!=my_NAN){          // is the patch sampled?
        assert((*curPop)->get_sampleID()!=SAMPLED);
        erase_vSamplePatch(*curPop);                 // the sampleID is needed to erase the patch!!!
    }
    (*curPop)->nbImmigrant = 0;
    (*curPop)->nbEmigrant = 0;
    
    // remove the element
    if(*curPop!= get_vFullPatch().back()){                 // if not last elemetn swap it with the last one
        *curPop = get_vFullPatch().back();                 // copy the last element to the current location
    }
    get_vFullPatch().pop_back();                          // remove the last element
    --endPop;
    assert(endPop==get_vFullPatch().end());
}

// ----------------------------------------------------------------------------------------
/** removes the current patch from the vector _vSamplePatch
 * (exchange the patch by the last one and remove the last element
 */
void
TMetapop::erase_vSamplePatch(TPatch* curPop)
{
    unsigned int sampleID = curPop->get_sampleID();
    assert(sampleID!=my_NAN);
    assert(sampleID!=SAMPLED);
    assert(_vSamplePatch[sampleID]->get_sampleID()==sampleID);
    
    _vSamplePatch[sampleID] = _vSamplePatch.back();  // copy the last element to the current location
    _vSamplePatch[sampleID]->set_sampleID(sampleID); // set its new sample id
    _vSamplePatch.pop_back();                        // remove the last element
    curPop->set_sampleID(SAMPLED);
}

// ----------------------------------------------------------------------------------------
/** add the immigrants to the db for coalescence simulations */
void
TMetapop::add_immigrants(const unsigned int& to, const unsigned int& from, const unsigned int& nb)
{
    assert(_pCoalescence_LCE);
    _pCoalescence_LCE->add_immigrants(_currentGeneration, to, from, nb);
}

// ----------------------------------------------------------------------------------------
/** add the current pop sizes to the db for coalescence simulations, but only if it is not empty */
void
TMetapop::store_popSizes()
{
    assert(_pCoalescence_LCE);
    _pCoalescence_LCE->store_popSizes(_currentGeneration);
}

// ----------------------------------------------------------------------------------------
/** get the sample sizes: the function. The function Patch::get_curN_sample()
 * cannot be used, since we want to know how many inds/pops cannot be sampled.
 */
bool
TMetapop::set_samples_coalescence()
{
    assert(_pCoalescence_LCE);
    unsigned int size, N, nbSampledPatches=0, notSampledPatches=0, notSampledInds=0, nbSampledInds=0;
    double sampleN;
    TPatch* curPatch;
    map<unsigned int, unsigned int> sampleSizes;
    vector<TPatch*>::iterator curPop, endPop = get_vPatch().end();
    for(curPop = get_vPatch().begin(); curPop!=endPop; ++curPop) {
        curPatch = *curPop;
        sampleN = curPatch->get_N_sample(FEM);  // sample size/proportion
        if(!sampleN) continue;                  // nothing is sampled: skip this patch
        N = curPatch->size(FEM, ADLTx);         // current pop size
        
        if(sampleN == my_NAN) size = N;                                  // not specified: sample all inds
        else if(sampleN < 1) size = (unsigned int) my_round(sampleN*N);
        else{                                                            // absolute number
            if(sampleN <= N) size = (unsigned int) sampleN;
            else{                                                          // if sampling exceeds the pop size
                size = N;
                notSampledInds += (unsigned int)sampleN - N;                 // not all individuals may be sampled
            }
        }
        
        if(size){
            sampleSizes[curPatch->get_ID()] = size;
            nbSampledInds += size;
            ++nbSampledPatches;
        }
        else  ++notSampledPatches;                                        // the patch should be sampled, but the pop size is zero or too small
    }
    if(sampleSizes.empty()){
        warning("All sampled patches are not populated: no stats computed!\n");
        _pCoalescence_LCE->set_samples(sampleSizes, 0);         // set the sample vector
        return false;
    }
    if(_sampleAllOrNothing && (notSampledInds || notSampledPatches)){
        warning("Some patches are not enough populated: no stats computed!\n");
        sampleSizes.clear();
        _pCoalescence_LCE->set_samples(sampleSizes, 0);         // set the sample vector
        return false;
    }
    _pCoalescence_LCE->set_samples(sampleSizes, nbSampledInds);         // set the sample vector
    if(notSampledPatches)   warning("Only %i of %i patches could be sampled: the others were not populated!\n", nbSampledPatches, nbSampledPatches + notSampledPatches);
    else if(notSampledInds) warning("All specified patches were sampled, but only %i of %i individuals could be sampled!\n", nbSampledInds, nbSampledInds + notSampledInds);
    return true;
}

// ----------------------------------------------------------------------------------------
void
TMetapop::set_func_pointer()
{
    for(unsigned int p=0; p<_patchNbr; ++p){
        _vPatch[p]->set_func_pointer();
    }
}

// ----------------------------------------------------------------------------------------
/** run the coalescence simulations scaled by the db */
void
TMetapop::run_coalescence()
{
    vector< LCE* >::iterator LCE;
    for(LCE = _theCycleCoal.begin(); LCE != _theCycleCoal.end(); ++LCE) {
        (*LCE)->execute();
        _currentAge ^= (*LCE)->removeAgeClass();
        _currentAge |= (*LCE)->addAgeClass();
    }
}

// ----------------------------------------------------------------------------------------
void
TMetapop::set_func_pointer_coal()
{
    for(unsigned int p=0; p<_patchNbr; ++p){
        _vPatch[p]->set_func_pointer_coal();
    }
}

// ----------------------------------------------------------------------------------------
void
TMetapop::set_func_ptr_noSample_noFull()
{
#ifdef _DEBUG
    message("Containers used: _vPatch\n");
#endif
    
    func_ptr_get_vSamplePatch = &TMetapop::get_vPatch;  _vSamplePatchUsed=false;
    func_ptr_get_vFullPatch   = &TMetapop::get_vPatch;  _vFullPatchUsed=false;
    
    func_ptr_add_tempPatch    = &TMetapop::add_tempPatch_noSample_noFull;
    func_ptr_new_emptyPatch   = &TMetapop::new_emptyPatch_noSample_noFull;
    func_ptr_new_fullPatch    = &TMetapop::new_fullPatch_noSample_noFull;
    
    assert(_pDisperse_LCE);
    _pDisperse_LCE->set_func_ptr_noSample_noFull();
}

// ----------------------------------------------------------------------------------------
void
TMetapop::set_func_ptr_noSample_noFull_coal()
{
#ifdef _DEBUG
    message("Containers used: _vPatch\n");
#endif
    
    func_ptr_get_vSamplePatch = &TMetapop::get_vPatch;  _vSamplePatchUsed=false;
    func_ptr_get_vFullPatch   = &TMetapop::get_vPatch;  _vFullPatchUsed=false;
    
    func_ptr_add_tempPatch    = &TMetapop::add_tempPatch_noSample_noFull_coal;
    func_ptr_new_emptyPatch   = &TMetapop::new_emptyPatch_noSample_noFull;
    func_ptr_new_fullPatch    = &TMetapop::new_fullPatch_noSample_noFull;
    
    assert(_pDisperse_LCE);
    _pDisperse_LCE->set_func_ptr_noSample_noFull();
}

// ----------------------------------------------------------------------------------------
void
TMetapop::set_func_ptr_withSample_noFull()
{
#ifdef _DEBUG
    message("Containers used: _vPatch, _vSamplePatch\n");
#endif
    
    func_ptr_get_vSamplePatch = &TMetapop::get_vSamplePatch_eff;  _vSamplePatchUsed=true;
    func_ptr_get_vFullPatch   = &TMetapop::get_vPatch;            _vFullPatchUsed=false;
    
    func_ptr_add_tempPatch    = &TMetapop::add_tempPatch_withSample_noFull;
    func_ptr_new_emptyPatch   = &TMetapop::new_emptyPatch_withSample_noFull;
    func_ptr_new_fullPatch    = &TMetapop::new_fullPatch_withSample_noFull;
    
    assert(_pDisperse_LCE);
    _pDisperse_LCE->set_func_ptr_withSample_noFull();
}

// ----------------------------------------------------------------------------------------
void
TMetapop::set_func_ptr_noSample_withFull()
{
#ifdef _DEBUG
    message("Containers used: _vPatch, _vFullPatch\n");
#endif
    
    func_ptr_get_vSamplePatch = &TMetapop::get_vFullPatch_eff;     _vSamplePatchUsed=false;
    func_ptr_get_vFullPatch   = &TMetapop::get_vFullPatch_eff;     _vFullPatchUsed=true;
    
    func_ptr_add_tempPatch    = &TMetapop::add_tempPatch_noSample_withFull;
    func_ptr_new_emptyPatch   = &TMetapop::new_emptyPatch_noSample_withFull;
    func_ptr_new_fullPatch    = &TMetapop::new_fullPatch_noSample_withFull;
    
    assert(_pDisperse_LCE);
    _pDisperse_LCE->set_func_ptr_noSample_withFull();
}

// ----------------------------------------------------------------------------------------
void
TMetapop::set_func_ptr_noSample_withFull_coal()
{
#ifdef _DEBUG
    message("Containers used: _vPatch, _vFullPatch\n");
#endif
    
    func_ptr_get_vSamplePatch = &TMetapop::get_vFullPatch_eff;    _vSamplePatchUsed=false;
    func_ptr_get_vFullPatch   = &TMetapop::get_vFullPatch_eff;    _vFullPatchUsed=true;
    
    func_ptr_add_tempPatch    = &TMetapop::add_tempPatch_noSample_withFull_coal;
    func_ptr_new_emptyPatch   = &TMetapop::new_emptyPatch_noSample_withFull_coal;
    func_ptr_new_fullPatch    = &TMetapop::new_fullPatch_noSample_withFull;
    
    assert(_pDisperse_LCE);
    _pDisperse_LCE->set_func_ptr_noSample_withFull();
}

// ----------------------------------------------------------------------------------------
void
TMetapop::set_func_ptr_withSample_withFull()
{
#ifdef _DEBUG
    message("Containers used: _vPatch, _vFullPatch, _vSamplePatch\n");
#endif
    
    func_ptr_get_vSamplePatch = &TMetapop::get_vSamplePatch_eff;    _vSamplePatchUsed=true;
    func_ptr_get_vFullPatch   = &TMetapop::get_vFullPatch_eff;      _vFullPatchUsed=true;
    
    func_ptr_add_tempPatch    = &TMetapop::add_tempPatch_withSample_withFull;
    func_ptr_new_emptyPatch   = &TMetapop::new_emptyPatch_withSample_withFull;
    func_ptr_new_fullPatch    = &TMetapop::new_fullPatch_withSample_withFull;
    
    assert(_pDisperse_LCE);
    _pDisperse_LCE->set_func_ptr_withSample_withFull();
}

// ----------------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
/** ouput of the first two columns for the stat file */
void
TMetapop::printGenRep2File(ostream& FH)
{
    FH.width(12);
    FH.setf(ios::left, ios::adjustfield);
    FH << setprecision(4);
    FH << getCurrentReplicate()+1 << "\t";
    FH.width(12);
    FH << getCurrentGeneration();
}

// ----------------------------------------------------------------------------------------
// update
// ----------------------------------------------------------------------------------------
/** updates the patche vectors of the colonized patches
	* must be called once before stats are computed
	*/
void TMetapop::update_patch_states()
{
    // get the number of patches to sample
    if(_current_nbSamplePatch) 	_last_nbSamplePatch = _current_nbSamplePatch; // the last non-zero sampled patch number
    _current_nbSamplePatch = get_nbSamplePatch();
    
    if(_current_index_stat_db==my_NAN) _current_index_stat_db=0;        // for the first round
    else ++_current_index_stat_db;
}

// ----------------------------------------------------------------------------------------
// getReplicateCounter
// ----------------------------------------------------------------------------------------
string
TMetapop::getReplicateCounter (){
    if (getReplicates()==1) return "";
    return toStr(_current_replicate+1, getReplicates());
}

// ----------------------------------------------------------------------------------------
// getReplicateCounter_
// ----------------------------------------------------------------------------------------
string
TMetapop::getReplicateCounter_ (){
    if (getReplicates()==1) return "";
    return "_" + toStr(_current_replicate+1, getReplicates());
}

// ----------------------------------------------------------------------------------------
// getReplicateCounter_r
// ----------------------------------------------------------------------------------------
string
TMetapop::getReplicateCounter_r (){
    if (getReplicates()==1) return "";
    return "_r" + toStr(_current_replicate+1, getReplicates());
}

// ----------------------------------------------------------------------------------------
// get_working_directory
// ----------------------------------------------------------------------------------------
string
TMetapop::get_working_directory(){
    assert(_pSimulation);
    return _pSimulation->get_working_directory();
}

// ----------------------------------------------------------------------------------------
// getSimfolder
// ----------------------------------------------------------------------------------------
string
TMetapop::getSimfolder(){
    assert(_pSimulation);
    return _pSimulation->get_simfolder();
}

// ----------------------------------------------------------------------------------------
// getBaseFileName
// ----------------------------------------------------------------------------------------
string
TMetapop::getBaseFileName(){
    assert(_pSimulation);
    return _pSimulation->get_basename();
}

// ----------------------------------------------------------------------------------------
// get_iniFile_directory
// ----------------------------------------------------------------------------------------
string
TMetapop::get_iniFile_directory(){
    assert(_pSimulation);
    return _pSimulation->get_iniFile_directory();
}

// ----------------------------------------------------------------------------------------
// get_exe_directory
// ----------------------------------------------------------------------------------------
string
TMetapop::get_exe_directory(){
    assert(_pSimulation);
    return _pSimulation->get_exe_directory();
}

// ----------------------------------------------------------------------------------------
// getSimfolderShort
// ----------------------------------------------------------------------------------------
string
TMetapop::getSimfolderShort (){
    assert(_pSimulation);
    return _pSimulation->get_simfolder_short();
}

// ----------------------------------------------------------------------------------------
// getSimfolderShort
// ----------------------------------------------------------------------------------------
RAND&
TMetapop::rand()
{
    assert(_pReplicate);
    return *_pReplicate->rand;
}

// ------------------------------------------------------------------------------
/** move randomly 'nbMigr' individuals  fom patch 'from_deme' to patch 'to_deme' */
void
TMetapop::move_random(sex_t SEX, age_idx from_age, unsigned int from_deme,
                     age_idx to_age, unsigned int to_deme, unsigned int nbMigr)
{
    unsigned int size = get_vPatch(from_deme)->size(SEX, from_age);
    TPatch* fromDeme = _vPatch[from_deme];
    TPatch* toDeme = _vPatch[to_deme];
    for(unsigned int i=0; i<nbMigr; ++i, --size){
        unsigned int at = rand().Uniform(size);
        toDeme->add(SEX, to_age, fromDeme->get(SEX, from_age, at));
        fromDeme->remove(SEX, from_age, at);
    }
}

// ------------------------------------------------------------------------------
/** move randomly 'nbMigr' individuals  fom patch 'from_deme' to patch 'to_deme' */
void
TMetapop::move_random(sex_t SEX, age_idx from_age, TPatch* fromDeme,
                     age_idx to_age, TPatch* toDeme, unsigned int nbMigr)
{
    unsigned int size = fromDeme->size(SEX, from_age);
    for(unsigned int i=0; i<nbMigr; ++i, --size){
        unsigned int at = rand().Uniform(size);
        toDeme->add(SEX, to_age, fromDeme->get(SEX, from_age, at));
        fromDeme->remove(SEX, from_age, at);
    }
}

// ------------------------------------------------------------------------------
/** move randomly 'nbMigr' individuals  fom patch 'from_deme' to patch 'to_deme'
 * do not delete the individual in the old pop!
 * "move with replacement"!
 */
void
TMetapop::copyMove_random_withReplacement(sex_t SEX, age_idx from_age, TPatch* fromDeme,
                                         age_idx to_age, TPatch* toDeme, unsigned int nbMigr)
{
    unsigned int size = fromDeme->size(SEX, from_age);
    assert(size);
    for(unsigned int i=0; i<nbMigr; ++i){
        unsigned int at = rand().Uniform(size);
        TIndividual* indOld = fromDeme->get(SEX, from_age, at);
        TIndividual* indCopy = copyIndividual(indOld);
        //delete indCopy;
        toDeme->add(SEX, to_age, copyIndividual(indOld));
    }
}

// ------------------------------------------------------------------------------
        
