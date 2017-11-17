/** @file tselectiontrait.h
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
//---------------------------------------------------------------------------

#ifndef tselectiontraitH
#define tselectiontraitH
//---------------------------------------------------------------------------
#include "simcomponent.h"


class TSelection;
class TTQuantiProto;
class Patch;
class Individual;


class TSelectionTrait{
protected:
    TSelection* _pSel;        // pointer to the selection object (do not delete)
    
    int _traitIndex;          // index of the quantitative trait among all traits
    int _quantiIndex;         // index of the quantitative trait among all quantitative traits
    
    // selection pressure
    unsigned int _nb_selection_params;
    double* _selection_pressure[2]; // if stabilizing selection (size = 2):
    //   0. pos: optima       (default: 0)
    //   1. pos: intensity    (default: 1)
    // if directional selection (size = 5):
    //   0. pos: min          (minimal fitness,             default: 0)
    //   1. pos: max          (maximal fitness,             default: 1)
    //   2. pos: max growth   (phenotype of maximal growth, default: 1)
    //   3. pos: growth rate  (steepness of the slope,      default: 1)
    //   4. pos: symmetry     (symmetry of the slope,       default: 1)
    // if fitness landscape  (size = 2*size_array+1)
    //   0. pos:           number of data points (size)
    //   1 -> size:        phenotype landscape
    //   size+1 -> 2*size: fitness landscape
    
    bool _selection_varies; // if false then set the pointer of _selection_pressure only to the correct patch
    // if true values vary between generations and thus have to be set in the
    // _selection_pressure array (_selection_pressure array has to be created and deleted!)
    
    
    
    
    double* _selection_sd;
    
    unsigned int _selection_model;  // 0: stabilizing selection (default)
    // 1: directional selection
    // 2: fitness landscape
    // 3: neutral
    // 4: selection coefficinet (bi-allelic)
    
    bool   _independent;
    
    
    
    // environment
    int     _environmental_model;
    double  _Ve_prop;         // Proportion of the used Ve (0: only natal patch; 1: only mating patch(default))
    double* _Ve_h2[2];        // heritability, repsectivly Ve directly  _Ve_hs[sex][patch]
    int     _Ve_model;        // 0: set VarE(0), VarE(t) constant;
    // 1: set h2(0),  VarE(t) constant;
    // 2: set h2(t), h2(t) constant, VarE(t) variable
    // 3: set H2(0),  VarE(t) constant;
    // 4: set H2(t), H2(t) constant, VarE(t) variable
    
    double _pheno;            // phenotype
    
    // current settings
    Patch*      _curPatch;
    sex_t       _curSex;
    Individual* _curInd;
    
    
    
    
    
public:
    Metapop*    _popPtr;      // pointer to the metapopulation

    TSelectionTrait();
    virtual ~TSelectionTrait(){
        if(_get_selection_pressure_func_ptr){
            delete[] _selection_pressure[FEM];
            delete[] _selection_pressure[MAL];
            delete[] _selection_sd;
            delete[] _get_selection_pressure_func_ptr;
        }
    }
    
    TTQuantiProto* _pQuantiProto;  // pointer to the corresponding quantitative trait (do not delete)
    
    void init(TSelection* s, const int& trait);
    virtual void init() = 0;
    virtual void set_ve_mean_func_ptr();
    
    double (TSelectionTrait::* _func_ptr_get_fitness)();
    double get_fitness(){return (this->*_func_ptr_get_fitness)();}
    virtual double get_fitnessFactor_trait();
    virtual double get_fitness_none() = 0;
    double get_fitnessFactor(){return get_fitness_none() * get_fitnessFactor_trait();}

    virtual int get_nb_selection_params() {return _nb_selection_params;}
    
    // functions to set the selection pressure (globally)
    void _get_selection_pressure_tot_var(Patch* patch, sex_t SEX);
    void _get_selection_pressure_tot_const(Patch* patch, sex_t SEX);
    void (TSelectionTrait::*_get_selection_pressure_tot_func_ptr)(Patch* patch, sex_t SEX);
    
    // functions to set the selection pressure (for each parameter individually)
    typedef double (TSelectionTrait::*_func_ptr)(double value, const int&);
    _func_ptr* _get_selection_pressure_func_ptr;
    double  _get_selection_pressure_var(double value, const int& i);
    double  _get_selection_pressure_const(double value, const int& i);
    
    // function pointer to get the Ve (from natal, current patch, or a mix of them)
    double (TSelectionTrait::*func_ptr_get_meanVe)();
    double (TSelectionTrait::*func_ptr_get_sdVe)();
    double  _get_null() {return 0;}
    double  _get_meanVe_current();
    double  _get_meanVe_natal();
    double  _get_meanVe_mix();
    double  _get_sdVe_current();
    double  _get_sdVe_natal();
    double  _get_sdVe_mix();
    double  get_Ve_prop()             const     {return _Ve_prop;}
    
    
    /** phenotype specific functions ********************************************/
    double   get_phenotype();
    void   (TSelectionTrait::* func_ptr_set_phenotype)(Individual* ind);
    void    set_phenotype(Individual* ind);
    void    set_phenotype_noVe(Individual* ind);
    void    set_phenotype_and_fitness_factor(Individual* ind);
    void    set_phenotype_and_fitness_factor_noVe(Individual* ind);
    void    set_fitness_factor(Individual* ind);
    double  set_get_fitness_factor(Individual* ind);
    
    void    set_Ve();
    void    set_quantiHeritability();
    
    
    // getter
    int     get_Ve_model    ( )             const  {return _Ve_model;}
    bool    get_independent ( )             const  {return _independent;}
    TTQuantiProto* get_pQuantiProto()       const  {return _pQuantiProto;}
    
    void setPhenotype_funcPointer(bool corr);
    
};
#endif
