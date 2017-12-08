/** @file ttquanti.cpp
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

#include "ttquanti.h"
#include "tmetapop.h"

#include "stathandler.cpp"
#include "tree.cpp"

#include "newmatap.h"                // need matrix applications
#include "newmatio.h"                // need matrix output routines


// ----------------------------------------------------------------------------------------
// copy constructor
// ----------------------------------------------------------------------------------------
TTraitQuanti::TTraitQuanti(const TTraitQuanti& T): _phenotype(my_NAN), _fitness(my_NAN), _fitness_factor(my_NAN), pProto(T.pProto){
    _genotype=T._genotype;
    _phenotype=T._phenotype;
    _fitness=T._fitness;
    _fitness_factor=T._fitness_factor;
    
    _copyTTraitParameters(T);  // copy the parameters of TTrait
}


// ----------------------------------------------------------------------------------------
// set_from_prototype (used during hatching)
// ----------------------------------------------------------------------------------------
void TTraitQuanti::set_from_prototype(TTraitProto* T)
{
    pProto = dynamic_cast<TTraitQuantiProto*> (T);
    pTraitProto = T;
}

// ----------------------------------------------------------------------------------------
// set_value
// ----------------------------------------------------------------------------------------
void
TTraitQuanti::set_value()
{
    _genotype = (pProto->*(pProto->get_genotype_func_ptr))(sequence);
}   // set the genotype


// ----------------------------------------------------------------------------------------
// set_fitness_factor
// ----------------------------------------------------------------------------------------
void
TTraitQuanti::set_fitness_factor()
{
    _fitness_factor = (pProto->*(pProto->get_fitnessFactor_func_ptr))(sequence);
}   // set the set_fitness_factor

// ----------------------------------------------------------------------------------------
// set_fitness_factor
// ----------------------------------------------------------------------------------------
double
TTraitQuanti::set_get_fitness_factor()
{
    _fitness_factor = (pProto->*(pProto->get_fitnessFactor_func_ptr))(sequence);
    return _fitness_factor;
}   // set the set_fitness_factor


//----------------------------------------------------------------------------------------
// operator=
// ----------------------------------------------------------------------------------------

TTraitQuanti&
TTraitQuanti::operator= (const TTrait& T)
{
    error("");
    return *this;
}

//----------------------------------------------------------------------------------------
// operator==
// ----------------------------------------------------------------------------------------
bool
TTraitQuanti::operator== (const TTrait& T)
{
    if(*pTraitProto != *T.pTraitProto) return false;
    const TTraitQuanti& TN = dynamic_cast<const TTraitQuanti&> (T);
    if(this != &TN || *pProto != *TN.pProto) return false;
    return true;
}

//----------------------------------------------------------------------------------------
// operator!=
// ----------------------------------------------------------------------------------------
bool
TTraitQuanti::operator!= (const TTrait& T)
{
    if(!((*this) == T)) return true;
    else                return false;
}

// ----------------------------------------------------------------------------------------
// Destructor
// ----------------------------------------------------------------------------------------
TTraitQuanti::~TTraitQuanti()
{
}

// ----------------------------------------------------------------------------------------
// reset
// ----------------------------------------------------------------------------------------
void
TTraitQuanti::reset()
{
    _phenotype      = my_NAN;
    _genotype       = my_NAN;
    _fitness       = my_NAN;
    _fitness_factor = my_NAN;
    
}

// ----------------------------------------------------------------------------------------
// show_up
// ----------------------------------------------------------------------------------------
void
TTraitQuanti::show_up ()
{
    message("\n  Trait's type: discretequanti\nlocus: %i\nalleles: %i\nsequence:", pProto->_nb_locus,pProto->get_nb_allele_max());
    
    for(unsigned int i = 0; (i < pProto->_nb_locus && i < 10); i++){
        message("\n              %i %i",(int)sequence[i][0],(int)sequence[i][1]);
    }
    message("\n");
}

// ----------------------------------------------------------------------------------------
// cstor
// ----------------------------------------------------------------------------------------
TTraitQuantiProto::TTraitQuantiProto( ):
_stats(0), _allelicValues(0), _allelic_file(0),
_dominanceValues(0),
_fitnessFactor_heterozygote(0), _fitnessFactor_homozygote(0), _fitnessFactor_freqDepend(0),
_fitnessFactor(0), _locusFreqs(0), _fitnessFactorTree(0), _phenoTree(0), get_genotype_func_ptr(0), _phenotyper(0), _writer(0), _genotyper(0), 
get_fitnessFactor_func_ptr(0), get_fitnessFactor2_func_ptr(0)
{
    _type = "quanti";
    ini_paramset();
}

// ----------------------------------------------------------------------------------------
// cstor
// ----------------------------------------------------------------------------------------
TTraitQuantiProto::TTraitQuantiProto(int i):
_stats(0), _allelicValues(0), _allelic_file(0),
_dominanceValues(0),
_fitnessFactor_heterozygote(0), _fitnessFactor_homozygote(0),
_fitnessFactor_freqDepend(0), _fitnessFactor(0), _locusFreqs(0),
_fitnessFactorTree(0), _phenotyper(0), _writer(0), _genotyper(0),
_phenoTree(0), get_genotype_func_ptr(0), get_fitnessFactor_func_ptr(0),
get_fitnessFactor2_func_ptr(0)
{
    _trait_index = i;
    _type = "quanti";
    ini_paramset();
}

// ----------------------------------------------------------------------------------------
// cstor
// ----------------------------------------------------------------------------------------
TTraitQuantiProto::TTraitQuantiProto(const TTraitQuantiProto& T):
_stats(0), _allelicValues(T._allelicValues), _allelic_file(T._allelic_file),
_dominanceValues(T._dominanceValues),
_fitnessFactor_heterozygote(T._fitnessFactor_heterozygote),
_fitnessFactor_homozygote(T._fitnessFactor_homozygote),
_fitnessFactor_freqDepend(T._fitnessFactor_freqDepend),
_fitnessFactor(T._fitnessFactor),
_locusFreqs(T._locusFreqs),
_fitnessFactorTree(0), _phenoTree(0),
get_genotype_func_ptr(T.get_genotype_func_ptr),
get_fitnessFactor_func_ptr(T.get_fitnessFactor_func_ptr),
get_fitnessFactor2_func_ptr(T.get_fitnessFactor2_func_ptr),
_phenotyper(0), _writer(0), _genotyper(0)
{
    _copyTraitPrototypeParameters(T);
    ini_paramset();
}

// ----------------------------------------------------------------------------------------
// dstor
// ----------------------------------------------------------------------------------------
TTraitQuantiProto::~TTraitQuantiProto ()
{
    if(_stats)          delete _stats;
    
    if(get_trait_index() <= 1){
        if(_phenotyper)     {delete _phenotyper; _phenotyper=NULL;} // _phenotyper is static!!!
        if(_genotyper)      {delete _genotyper;  _genotyper=NULL;}  // _genotyper is static!!!
        if(_writer)         {delete _writer;     _writer=NULL;}     //_writer is static!!!
    }
    
    resetTotal();
}

// -----------------------------------------------------------------------------
// ini_paramset
// -----------------------------------------------------------------------------
void
TTraitQuantiProto::ini_paramset ()
{
    TTraitProto::ini_paramset();
    string trait = get_trait_indexStr_();
    
    add_parameter("quanti_allelic_var"+trait,DBL,false,0,my_NAN,"1",false,
                  "Variance of the normal distribution from where the allelic effects are drawn randomly.",0);
    
    add_parameter("quanti_dominance_model"+trait,INT2,false,0,1,"0",false,
                  "Which dominance model to use:\n" \
                  "  0: method k: G = a1 + a2 + K*abs(a1-a2)\n" \
                  "  1: method h: G = 2[(1-h)a1 + h*a2]",2);
    
    add_parameter("quanti_dominance_mean"+trait,DBL,false,my_NAN,my_NAN,"0",false,
                  "Mean of the normal distribution from where the dominance effects are randomly drawn.",2);
    
    add_parameter("quanti_dominance_var"+trait,DBL,false,0,my_NAN,"0",false,
                  "Variance of the normal distribution from where the dominance effects are randomly drawn.",2);
    
    add_parameter("quanti_dominance_file"+trait,STR,false,my_NAN,my_NAN,"",false,
                  "Filename containing information about a locus, such as the dominance " \
                  "effects, and/or the fitness factors. " \
                  "Note, no '%' is assumed in front of the filename.",2);
    
    
    add_parameter("quanti_fitness_factor_heterozygote"+trait,DBL_MAT,false,0,my_NAN,"1",false,
                  "Fitness factor at heterozygote loci.",4);
    
    add_parameter("quanti_fitness_factor_homozygote"+trait,DBL_MAT,false,0,my_NAN,"1",false,
                  "Fitness factor at homozygote loci.",4);
    
    add_parameter("quanti_fitness_frequency_dependent"+trait, DBL_MAT,false,my_NAN,my_NAN,"0",false,
                  "TODO",5);
    
    
    add_parameter("quanti_epistatic_var"+trait,DBL,false,0,my_NAN,"0",false,
                  "Variance of the normal distribution from where the epistatic effects are randomly drawn.",4);
    
    add_parameter("quanti_epistatic_file"+trait,STR,false,my_NAN,my_NAN,"",false,
                  "Filename containing information about a genotype, such as the epistatic " \
                  "effects, and/or the fitness factors. " \
                  "Note, no '%' is assumed in front of the filename.",4);
    
    
    add_parameter("quanti_mutation_model"+trait,INT2,false,0,3,"0",false,
                  "The mutation model for QTLs:\n" \
                  "  0: RMM (Random mutation model). A new effect is drawn randomly " \
                  "in a normal distribution centered at zero.\n" \
                  "  1: IMM (Increment mutation model) A new effect is drawn randomly " \
                  "in a normal distribution centered at zero and added to the old effect.\n" \
                  "  2: KAM (K-Allele mutation model) Mutation to any other allele.\n" \
                  "  3: SSM (Single step mutation model) Mutation to next bigger or smaller effect.",0);
    
    
    add_parameter("quanti_multiplicative_model"+trait,INT2,false,0,1,"0",false,
                  "TODO",5);
    
    
    add_parameter("quanti_output"+trait,INT2,false,0,1,"0",false,
                  "Should the architecture of the trait be outputted:\n" \
                  "  0: none (no output)\n" \
                  "  1: output (allelic, dominance and epistatic files are outputted if used).",2);
    
    
    add_parameter("quanti_save_phenotype",INT2,false,0,2,"0",false,
                  "Should the phenotypes be outputted:\n" \
                  "  0: none (no output is generated)\n" \
                  "  1: standard (phenotypes are outputted in the standard FSTAT-like format)\n" \
                  "  2: extended (same as point 1, but the file contain the following six " \
                  "additional columns: the age class (1 = offspring, 2 = adult), the " \
                  "sex (0 = male; 1 = female), the ID of the individual, the ID of the " \
                  "mother, the ID of the father, and the fitness of the individual).",3);
    
    add_parameter("quanti_phenot_sex",INT2,false,0,2,"0",false,
                  "Which sex(es) should the phenotype output contain:\n" \
                  "  0: both (output includes both sexes)\n" \
                  "  1: females (output includes only female phenotypes)\n" \
                  "  2: males (output includes only male phenotypes).",3);
    
    add_parameter("quanti_phenot_age",INT2,false,0,2,"0",false,
                  "Which ag(es) should the phenotype output contain:\n" \
                  "  0: adults (output includes only adult phenotypes)\n" \
                  "  1: juveniles (output includes only juvenile phenotypes)\n" \
                  "  2: both (output includes juvenile and adult phenotypes)\n" \
                  "Note, that phenotypes are not always already available for juveniles. TODO",3);
    
    add_parameter("quanti_phenot_dir",STR,false,my_NAN,my_NAN,"",false,
                  "The directory name where phenotypes are stored.",3);
    
    add_parameter("quanti_phenot_filename",STR,false,my_NAN,my_NAN,"",false,
                  "The base name for the phenotype output files.",3);
    
    add_parameter("quanti_phenot_logtime",INT2,false,0,my_NAN,"1",true,
                  "The time interval of the phenotype output.",3);
    
    add_parameter("quanti_phenot_script",STR,false,my_NAN,my_NAN,"",false,
                  "The script which will be launched after the phenotype file is generated.",3);
    
    add_parameter("quanti_save_geno_value",INT2,false,0,2,"0",false,
                  "Should the genotpic values be outputted:\n" \
                  "  0: none (no output is generated)\n" \
                  "  1: standard (genotypic values are outputted in the standard FSTAT-like format)\n" \
                  "  2: extended (same as point 1, but the file contain the following six " \
                  "additional columns: the age class (1 = offspring, 2 = adult), the " \
                  "sex (0 = male; 1 = female), the ID of the individual, the ID of the " \
                  "mother, the ID of the father, and the fitness of the individual).",3);
    
    add_parameter("quanti_geno_value_sex",INT2,false,0,2,"0",false,
                  "Which sex(es) should the genotypic values output contain:\n" \
                  "  0: both (output includes both sexes)\n" \
                  "  1: females (output includes only female genotypic values)\n" \
                  "  2: males (output includes only male genotypic values).",3);
    
    add_parameter("quanti_geno_value_age",INT2,false,0,2,"0",false,
                  "Which ag(es) should the genotypic values output contain:\n" \
                  "  0: adults (output includes only adult genotypic value)\n" \
                  "  1: juveniles (output includes only juvenile genotypic value)\n" \
                  "  2: both (output includes juvenile and adult genotypic value).",3);
    
    add_parameter("quanti_geno_value_dir",STR,false,my_NAN,my_NAN,"",false,
                  "The directory name where genotypic values are stored.",3);
    
    add_parameter("quanti_geno_value_filename",STR,false,my_NAN,my_NAN,"",false,
                  "The base name for the genotypic values output files.",3);
    
    add_parameter("quanti_geno_value_logtime",INT2,false,0,my_NAN,"1",true,
                  "The time interval of the genotypic values output.",3);
    
    add_parameter("quanti_geno_value_script",STR,false,my_NAN,my_NAN,"",false,
                  "The script which will be launched after the genotypic values file is generated.",3);
    
    
    add_parameter("quanti_selection_model"+trait,INT2,false,0,4,"0",false,
                  "Specification of the selection pressure type:\n" \
                  "  0: neutral quantitative trait (no selection)\n" \
                  "  1: stabilizing selection\n" \
                  "  2: directional selection\n" \
                  "  3: fitness landscape\n" \
                  "  4: selection coefficient (only bi-allelic loci).",1);
    
    
    // environment
    add_parameter("quanti_heritability"+trait,DBL,false,0,my_NAN,"0",false,
                  "Specifies the environmental variance / heritability for each patch. "
                  "The applied model is specified by the parameter 'quanti_environmental_model'.",2);
    
    add_parameter("quanti_heritability_fem"+trait,DBL,false,0,my_NAN,"0",false,
                  "Specifies the environmental variance / heritability for females and for each patch. " \
                  "The applied model is specified by the parameter 'quanti_environmental_model'.",2);
    
    add_parameter("quanti_heritability_mal"+trait,DBL,false,0,my_NAN,"0",false,
                  "Specifies the environmental variance / heritability for males and for each patch. \n" \
                  " The applied model is specified by the parameter 'quanti_environmental_model'.",2);
    
    add_parameter("quanti_environmental_proportion"+trait,DBL,false,0,1,"1",false,
                  "Specifies which patch affects the environmental variance (0: only natal patch; 1: only current patch).",5);
    
    add_parameter("quanti_environmental_model"+trait,INT2,false,0,4,"0",false,
                  "How should the environmental variance / heritability be defined:\n" \
                  "  0: Ve directly set\n" \
                  "  1: Ve defined by the narrow-sense heritability (VE constant)\n" \
                  "  2: Ve defined by the narrow-sense heritability (h^2 constant)\n" \
                  "  3: Ve defined by the broad-sense heritability (Ve constant)\n" \
                  "  4: Ve defined by the broad-sense heritability (H^2 constant)\n" \
                  "The value itself is set by the parameter 'quanti_heritability'.",2);
    
    
    add_parameter("quanti_va_model"+trait,INT2,false,0,2,"0",false,
                  "How is the genetic additive variance computed:\n" \
                  "  0: for any case (computational intensive; not always computable)\n" \
                  "  1: limited to random mating\n" \
                  "  2: Va = Vg (additive genetic variance is approximated by the genetic variance).",5);
    
    
    // stabilizing selection
    add_parameter("quanti_stab_sel_optima"+trait,			  DBL_MAT,false,my_NAN,my_NAN,"0",true,
                  "Stabilizing selection: the fitness optimum.",0);
    
    add_parameter("quanti_stab_sel_optima_fem"+trait,	  DBL_MAT,false,my_NAN,my_NAN,"0",true,
                  "Stabilizing selection: the fitness optimum for females.",2);
    
    add_parameter("quanti_stab_sel_optima_mal"+trait,	  DBL_MAT,false,my_NAN,my_NAN,"0",true,
                  "Stabilizing selection: the fitness optimum for males.",2);
    
    add_parameter("quanti_stab_sel_intensity"+trait,    DBL_MAT,false,0,		 my_NAN,"1",true,
                  "Stabilizing selection: the selection intensity.\n" \
                  "Caution: the higher the value is, the weaker the selection pressure is.",0);
    
    add_parameter("quanti_stab_sel_intensity_fem"+trait,DBL_MAT,false,0,	   my_NAN,"1",true,
                  "Stabilizing selection: the selection intensity for females.\n" \
                  "Caution: the higher the value is, the weaker the selection pressure is.",2);
    
    add_parameter("quanti_stab_sel_intensity_mal"+trait,DBL_MAT,false,0,	   my_NAN,"1",true,
                  "Stabilizing selection: the selection intensity for males.\n" \
                  "Caution: the higher the value is, the weaker the selection pressure is.",2);
    
    
    // directional selection
    add_parameter("quanti_dir_sel_min"+trait,						DBL_MAT,false,0,		 1,		  "0",true,
                  "Directional selection: lower asymptote of selection curve.",0);
    
    add_parameter("quanti_dir_sel_min_fem"+trait,				DBL_MAT,false,0,		 1,			"0",true,
                  "Directional selection: lower asymptote of selection curve for females.",2);
    
    add_parameter("quanti_dir_sel_min_mal"+trait,				DBL_MAT,false,0,		 1,			"0",true,
                  "Directional selection: lower asymptote of selection curve for males.",2);
    
    add_parameter("quanti_dir_sel_min_var"+trait,				DBL_MAT,false,0,		 my_NAN,"0",true,
                  "Directional selection: variance of the lower asymptote over time.",5);
    
    add_parameter("quanti_dir_sel_max"+trait,						DBL_MAT,false,0,		 1,			"1",true,
                  "Directional selection: upper asymptote of selection curve.",0);
    
    add_parameter("quanti_dir_sel_max_fem"+trait,				DBL_MAT,false,0,		 1,			"1",true,
                  "Directional selection: upper asymptote of selection curve for females.",2);
    
    add_parameter("quanti_dir_sel_max_mal"+trait,				DBL_MAT,false,0,		 1,			"1",true,
                  "Directional selection: lower asymptote of selection curve for males.",2);
    
    add_parameter("quanti_dir_sel_max_var"+trait,				DBL_MAT,false,0,		 my_NAN,"0",true,
                  "Directional selection: variance of the upper asymptote over time.",5);
    
    add_parameter("quanti_dir_sel_max_growth"+trait,	 	DBL_MAT,false,my_NAN,my_NAN,"0",true,
                  "Directional selection: phentoype with maximal growth.",0);
   
    add_parameter("quanti_dir_sel_max_growth_fem"+trait,DBL_MAT,false,my_NAN,my_NAN,"0",true,
                  "Directional selection: phenotype with maximal growth for females.",2);
    
    add_parameter("quanti_dir_sel_max_growth_mal"+trait,DBL_MAT,false,my_NAN,my_NAN,"0",true,
                  "Directional selection: phenotype with maximal growth for males.",2);
    
    add_parameter("quanti_dir_sel_max_growth_var"+trait,DBL_MAT,false,0,	 	 my_NAN,"0",true,
                  "Directional selection: variance of the phenotype with maximal growth over time.",5);
    
    add_parameter("quanti_dir_sel_growth_rate"+trait,		DBL_MAT,false,my_NAN,my_NAN,"1",true,
                  "Directional selection: growth rate (when positive, larger values have a " \
                  "higher fitness; when negative, lower phenotypes have a higher fitness.",0);
    
    add_parameter("quanti_dir_sel_growth_rate_fem"+trait,DBL_MAT,false,my_NAN,my_NAN,"1",true,
                  "Directional selection: growth rate fir females (when positive, larger values have a " \
                  "higher fitness; when negative, lower phenotypes have a higher fitness.",2);
    
    add_parameter("quanti_dir_sel_growth_rate_mal"+trait,DBL_MAT,false,my_NAN,my_NAN,"1",true,
                  "Directional selection: growth rate for males (when positive, larger values have a " \
                  "higher fitness; when negative, lower phenotypes have a higher fitness.",2);
    
    add_parameter("quanti_dir_sel_growth_rate_var"+trait,DBL_MAT,false,0,	   my_NAN,"0",true,
                  "Directional selection: variance of the growth rate over time.",5);
    
    add_parameter("quanti_dir_sel_symmetry"+trait,			DBL_MAT,false,0,		 1,			"1",true,
                  "Directional selection: symmetry of the slope.",0);
    
    add_parameter("quanti_dir_sel_symmetry_fem"+trait,	DBL_MAT,false,0,		 1,			"1",true,
                  "Directional selection: symmetry of the slope for females.",2);
    
    add_parameter("quanti_dir_sel_symmetry_mal"+trait,	DBL_MAT,false,0,		 1,			"1",true,
                  "Directional selection: symmetry of the slope for males.",2);
    
    add_parameter("quanti_dir_sel_symmetry_var"+trait,	DBL_MAT,false,0,		 my_NAN,"0",true,
                  "Directional selection: variance of the symmetry over time.",5);
    
    
    // fitness landscape
    add_parameter("quanti_fitness_landscape"+trait,   	DBL_MAT,false,my_NAN,my_NAN,"1",true,
                  "Fitness landscape: Specifies the fitness values of the phenotype-fitness " \
                  "landscape (to be used with parameter 'patch_phenotype_landscape').",0);
    
    add_parameter("quanti_fitness_landscape_fem"+trait,	DBL_MAT,false,my_NAN,my_NAN,"1",true,
                  "Fitness landscape: Specifies the fitness values for females of the phenotype-fitness " \
                  "landscape (to be used with parameter 'patch_phenotype_landscape').",2);
    
    add_parameter("quanti_fitness_landscape_mal"+trait,	DBL_MAT,false,my_NAN,my_NAN,"1",true,
                  "Fitness landscape: Specifies the fitness values for males of the phenotype-fitness " \
                  "landscape (to be used with parameter 'patch_phenotype_landscape').",2);
    
    add_parameter("quanti_phenotype_landscape"+trait, 	DBL_MAT,false,my_NAN,my_NAN,"0",true,
                  "Fitness landscape: Specifies the phenotype values of the phenotype-fitness " \
                  "landscape (to be used with parameter 'patch_fitness_landscape').",0);
    
    add_parameter("quanti_phenotype_landscape_fem"+trait,DBL_MAT,false,my_NAN,my_NAN,"0",true,
                  "Fitness landscape: Specifies the phenotype values for females of the phenotype-fitness " \
                  "landscape (to be used with parameter 'patch_fitness_landscape_fem').",2);
    
    add_parameter("quanti_phenotype_landscape_mal"+trait,DBL_MAT,false,my_NAN,my_NAN,"0",true,
                  "Fitness landscape: Specifies the phenotype values for males of the phenotype-fitness " \
                  "landscape (to be used with parameter 'patch_fitness_landscape_mal').",2);
    
    
    // selection coefficient
    add_parameter("quanti_coef_sel_AA"+trait,						DBL_MAT,false,my_NAN,my_NAN,"1",true,
                  "Selection coefficient: specifies the fitness for wild type AA.",0);
    
    add_parameter("quanti_coef_sel_AA_fem"+trait,					DBL_MAT,false,my_NAN,my_NAN,"1",true,
                  "Selection coefficient: specifies the fitness for wild type AA of females.",2);
    
    add_parameter("quanti_coef_sel_AA_mal"+trait,					DBL_MAT,false,my_NAN,my_NAN,"1",true,
                  "Selection coefficient: specifies the fitness for wild type AA of males.",2);
    
    add_parameter("quanti_coef_sel"+trait,						    DBL_MAT,false,my_NAN,my_NAN,"0",true,
                  "Selection coefficient: The selection coefficient 's'.",0);
    
    add_parameter("quanti_coef_sel_fem"+trait,						DBL_MAT,false,my_NAN,my_NAN,"0",true,
                  "Selection coefficient: The selection coefficient of females 's'.",2);
    
    add_parameter("quanti_coef_sel_mal"+trait,						DBL_MAT,false,my_NAN,my_NAN,"0",true,
                  "Selection coefficient: The selection coefficient of males 's'.",2);
    
}

// ----------------------------------------------------------------------------------------
// reset
// ----------------------------------------------------------------------------------------
// called between replicates to reset this object
void TTraitQuantiProto::reset()
{
    TTraitProto::reset();
    
    string trait = get_trait_indexStr_();
    //----------------------------------------------------------------------------
    // read locus file
    if(_dominanceValues){
        delete_dominanceValues();
        string dominanceFile = get_parameter("quanti_dominance_file"+trait)->get_arg();
        if(!dominanceFile.empty()){      // if the parameter is set
            if(STRING::file_exists(dominanceFile, get_popPtr()->get_iniFile_directory())) read_locus_file(dominanceFile);
            else error("Dominance file '%s' (quanti) cannot be found!\n", dominanceFile.c_str());
        }
        set_dominanceValues(trait);
    }
    
    //----------------------------------------------------------------------------
    // epistatic values
    if(_phenoTree || _fitnessFactorTree){
        if(_phenoTree)        {delete _phenoTree; _phenoTree=NULL;}
        if(_fitnessFactorTree){delete _fitnessFactorTree; _fitnessFactorTree=NULL;}
        string epistaticFile = get_parameter("quanti_epistatic_file"+trait)->get_arg();
        if(!epistaticFile.empty()){       // epistatic effects set by file
            if(STRING::file_exists(epistaticFile, get_popPtr()->get_iniFile_directory())) read_genome_file(epistaticFile);
            else error("Epistatic file '%s' (quanti) cannot be found!\n", epistaticFile.c_str());
        }
        set_epistaticValues(trait);
    }
}

// ----------------------------------------------------------------------------------------
// delete_dominanceValues
// ----------------------------------------------------------------------------------------
void TTraitQuantiProto::delete_dominanceValues()
{
    if(!_dominanceValues) return;  // if not used
    for(unsigned int a, l=0; l<_nb_locus; ++l){
        for(a=0; a<_nb_allele[l]; ++a){
            delete[] _dominanceValues[l][a];
        }
        delete[] _dominanceValues[l];
    }
    delete[] _dominanceValues;
    _dominanceValues = NULL;
}

// ----------------------------------------------------------------------------------------
// delete_allelicValues
// ----------------------------------------------------------------------------------------
void TTraitQuantiProto::delete_allelicValues()
{
    if(!_allelicValues) return;
    if(_nb_locus > 1 && _allelicValues[0] != _allelicValues[1]){   // each locus has its one array
        ARRAY::delete_2D(_allelicValues, _nb_locus);
    }
    else{                                                        // all loci point to the first array
        delete[] _allelicValues[0];                                // delete the first array
        ARRAY::delete_1D(_allelicValues);
    }
}

// ----------------------------------------------------------------------------------------
// delete_fitnessFactor
// ----------------------------------------------------------------------------------------
void TTraitQuantiProto::delete_fitnessFactor()
{
    if(!_fitnessFactor) return;  // if not used
    for(unsigned int a, l=0; l<_nb_locus; ++l){
        for(a=0; a<_nb_allele[l]; ++a){
            delete[] _fitnessFactor[l][a];
        }
        delete[] _fitnessFactor[l];
    }
    delete[] _fitnessFactor;
    _fitnessFactor = NULL;
}


// ----------------------------------------------------------------------------------------
// delete_locusFreqs
// ----------------------------------------------------------------------------------------
void TTraitQuantiProto::delete_locusFreqs()
{
    if(!_locusFreqs) return;  // if not used
    delete[] _locusFreqs;
    _locusFreqs = NULL;
}

// ----------------------------------------------------------------------------------------
// reset
// ----------------------------------------------------------------------------------------
// used between simulations
void TTraitQuantiProto::resetTotal()
{
    if(_phenoTree)        {delete _phenoTree;		_phenoTree=NULL;}
    if(_fitnessFactorTree){delete _fitnessFactorTree;		_fitnessFactorTree=NULL;}
    delete_dominanceValues();
    delete_fitnessFactor();
    delete_locusFreqs();
    if(_fitnessFactor_heterozygote) {delete[] _fitnessFactor_heterozygote; _fitnessFactor_heterozygote=NULL;}
    if(_fitnessFactor_homozygote) {delete[] _fitnessFactor_homozygote; _fitnessFactor_homozygote=NULL;}
    if(_fitnessFactor_freqDepend) {delete[] _fitnessFactor_freqDepend; _fitnessFactor_freqDepend=NULL;}
    delete_allelicValues();
    
    TTraitProto::resetTotal();
}

// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
void TTraitQuantiProto::init (TMetapop* pMetapop)
{
    mut_model_t mutation_models[4] = {RMM, IMM, KAM, SSM};
    ini_base(pMetapop, mutation_models, 4);
    string trait = get_trait_indexStr_();
    
    assert(!_writer);
    assert(!_phenotyper);
    assert(!_genotyper);
    
    //----------------------------------------------------------------------------
    // get quanti specific parameters
    _output             = (int) get_parameter_value("quanti_output"+trait);
    _selection_model    = (unsigned int) get_parameter_value("quanti_selection_model"+trait);
    _Va_model           = (int) get_parameter_value("quanti_va_model"+trait);
    
    //----------------------------------------------------------------------------
    // read the allelic file if present
    string allelicFile = get_parameter("quanti_allelic_file"+trait)->get_arg();  // if empty parameter is not set
    if(!allelicFile.empty()){
        allelicFile = get_popPtr()->get_iniFile_directory()+allelicFile;
        if(STRING::file_exists(allelicFile)) read_allele_file(allelicFile);
        else error("Allelic file '%s' (quanti) cannot be found!\n", allelicFile.c_str());
        
    }
    
    //----------------------------------------------------------------------------
    // read locus file
    string dominanceFile = get_parameter("quanti_dominance_file"+trait)->get_arg();
    if(!dominanceFile.empty()){      // if the parameter is set
        dominanceFile = get_popPtr()->get_iniFile_directory()+dominanceFile;
        if(STRING::file_exists(dominanceFile)) read_locus_file(dominanceFile);
        else error("Dominance file '%s' (quanti) cannot be found!\n", dominanceFile.c_str());
    }
    set_dominanceValues(trait);
    
    //----------------------------------------------------------------------------
    // read genome file
    string epistaticFile = get_parameter("quanti_epistatic_file"+trait)->get_arg();
    if(!epistaticFile.empty()){       // epistatic effects set by file
        epistaticFile = get_popPtr()->get_iniFile_directory()+epistaticFile;
        if(STRING::file_exists(epistaticFile)) read_genome_file(epistaticFile);
        else error("Epistatic file '%s' (quanti) cannot be found!\n", epistaticFile.c_str());
    }
    set_epistaticValues(trait);
    
    //----------------------------------------------------------------------------
    // set the mutation rate probabilities
    set_allelicValues(trait);
    set_mutationFreq(trait);
    
    //----------------------------------------------------------------------------
    // set the initial allele frequencies
    set_initAlleleFreq(trait);
    _protoGenome->set_ini_sequence_model(_aLocus, _nb_locus, _ini_allele_model); // following a distribution
    
    //----------------------------------------------------------------------------
    // set the fitness factor
    set_fitnessFactor(trait, "quanti_fitness_factor_heterozygote", _fitnessFactor_heterozygote);
    set_fitnessFactor(trait, "quanti_fitness_factor_homozygote", _fitnessFactor_homozygote);
    if(_fitnessFactorTree)  get_fitnessFactor_func_ptr = &TTraitQuantiProto::get_fitnessFactor_genome;
    else if(_fitnessFactor) get_fitnessFactor_func_ptr = &TTraitQuantiProto::get_fitnessFactor_locus;
    else if(_fitnessFactor_homozygote || _fitnessFactor_heterozygote){
        get_fitnessFactor_func_ptr = &TTraitQuantiProto::get_fitnessFactor_global;
    }
    else get_fitnessFactor_func_ptr = NULL;	// fitness factor is not used
    
    
    // frequency dependent selection?
    set_fitnessFactor(trait, "quanti_fitness_frequency_dependent", _fitnessFactor_freqDepend);
    if(_fitnessFactor_freqDepend){
        get_fitnessFactor2_func_ptr = get_fitnessFactor_func_ptr;
        get_fitnessFactor_func_ptr = &TTraitQuantiProto::get_fitnessFactor_freqDepend;
    }
}

// ----------------------------------------------------------------------------------------
// temporal_change
// ----------------------------------------------------------------------------------------
void
TTraitQuantiProto::temporal_change(const unsigned int& gen)
{
    // if it is a temporal parameter
    map<string, Param*>* pParam = _paramSet->getTemporalParams(gen);
    if(pParam){
        // check if a change has to be made
        map<string, Param*>* pMap = _paramSet->updateTemporalParams(gen);
        if(pMap){
            TTraitProto::temporal_change(gen, pMap);    // base class changes?
            
            map<string, Param*>::iterator pos = pMap->begin();
            bool fit_land=false;
            string trait = get_trait_indexStr_();
            unsigned int traitIndex = get_trait_index();
            if(traitIndex) traitIndex -= 1;
            for(; pos != pMap->end(); ++pos){
                // stabilizing selection
                if(pos->first == "quanti_stab_sel_optima"+trait){
                    _popPtr->set_patch_parameter_ofTrait(this, traitIndex, trait, "quanti_stab_sel_optima",     "optima", &TPatch::set_localOptima);    // default 0
                }
                else if(pos->first == "quanti_stab_sel_intensity"+trait){
                    _popPtr->set_patch_parameter_ofTrait(this, traitIndex, trait, "quanti_stab_sel_intensity",  "intensity", &TPatch::set_localIntensity);    // default 1
                }
                
                // directional selection
                else if(pos->first == "quanti_dir_sel_growth_rate"+trait){
                    _popPtr->set_patch_parameter_ofTrait(this, traitIndex, trait, "quanti_dir_sel_growth_rate", "growth_rate", &TPatch::set_localGrowthRate);  // default 1
                }
                else if(pos->first == "quanti_dir_sel_max_growth"+trait){
                    _popPtr->set_patch_parameter_ofTrait(this, traitIndex, trait, "quanti_dir_sel_max_growth",  "max_growth",  &TPatch::set_localMaxGrowth);   // default 1
                }
                else if(pos->first == "quanti_dir_sel_symmetry"+trait){
                    _popPtr->set_patch_parameter_ofTrait(this, traitIndex, trait, "quanti_dir_sel_symmetry",    "symmetry",    &TPatch::set_localSymmetry);    // default 0.5
                }
                else if(pos->first == "quanti_dir_sel_min"+trait){
                    _popPtr->set_patch_parameter_ofTrait(this, traitIndex, trait, "quanti_dir_sel_min",         "min",         &TPatch::set_localMin);         // default 0
                }
                else if(pos->first == "quanti_dir_sel_max"+trait){
                    _popPtr->set_patch_parameter_ofTrait(this, traitIndex, trait, "quanti_dir_sel_max",         "max",         &TPatch::set_localMax);         // default 1
                }
                
                // fitness landscape (both have to be changed at the satraitIndexme time: do it later)
                else if(pos->first == "quanti_fitness_landscape"+trait
                        || pos->first == "quanti_phenotype_landscape"+trait)     fit_land =   true;
                
                if(_trait_index<=1){     // the logtimes must only be changed for the first trait!!!
                    if(pos->first == "quanti_genot_logtime" && _writer){
                        _writer->set_GenerationOccurrence((unsigned int)pos->second->get_value());
                    }
                    else if(pos->first == "quanti_phenot_logtime" && _phenotyper){
                        _phenotyper->set_GenerationOccurrence((unsigned int)pos->second->get_value());
                    }
                    else if(pos->first == "quanti_geno_value_logtime" && _genotyper){
                        _genotyper->set_GenerationOccurrence((unsigned int)pos->second->get_value());
                    }
                }
            }
            
            // temporal changes have been made: perform if necessary the last changes
            if(fit_land){   // (phenotype has to be set after fitness due to sorting)
                _popPtr->set_patch_parameter_array_ofTrait(this, traitIndex, trait, "patch_fitness_landscape",   "fitness_landscape", &TPatch::set_fitnessLandscape_fitness);
                _popPtr->set_patch_parameter_array_ofTrait(this, traitIndex, trait, "patch_phenotype_landscape", "phenotype_landscape", &TPatch::set_fitnessLandscape_phenotype);
            }
        }
    }
}

// ----------------------------------------------------------------------------------------
// loadFileServices
// ----------------------------------------------------------------------------------------
void
TTraitQuantiProto::executeAfterEachReplicate(const unsigned int& rep)
{
    if(_output){
        string dir =get_popPtr()->getSimfolder();
        print_allelic_values(dir+"allelic_values");
        print_dominance_values(dir+"dominance_values");
        print_epistatic_values(dir+"epistatic_values");
    }
}

// ----------------------------------------------------------------------------------------
// get_info
// ----------------------------------------------------------------------------------------
string
TTraitQuantiProto::get_info()
{
    string text;
    
    if(_trait_index) text = toStr(_trait_index) + ". quantitative trait: ";
    else             text = "Quantitative trait: ";
    
    text +=	get_info_locus();
    
    switch(_selection_model){
        case 0: text += " neutral trait;";            	 break;
        case 1: text += " stabilizing selection;";    	 break;
        case 2: text += " directional selection;";       break;
        case 3: text += " fitness landscape;";           break;
        case 4: text += " selection coefficient (bi-allelic);"; break;
    }
    
    return text;
}

// ----------------------------------------------------------------------------------------
// loadFileServices
// ----------------------------------------------------------------------------------------
void
TTraitQuantiProto::loadFileServices  (FileServices* loader)
{
    // genotype
    unsigned int choice = (int)get_parameter_value("quanti_save_genotype");
    assert(!_writer);
    if(choice) {
        _writer = (dynamic_cast<TTraitQuantiProto*>(_popPtr->getFirstPrototype(_type)))->_writer;
        if(!_writer){
            _writer = new TTQuantiFH();
            if(_writer->set(get_parameter("quanti_genot_logtime"),
                            get_parameter("quanti_genot_dir")->get_arg(),
                            get_parameter("quanti_genot_filename")->get_arg(),
                            get_parameter("quanti_genot_script")->get_arg(),
                            (int)get_parameter_value("quanti_genot_sex"),
                            (int)get_parameter_value("quanti_genot_age"),
                            choice,
                            this,
                            "genotype",
                            ".dat", loader, get_popPtr())){
                loader->attach(_writer);     // only the first one
            }
            else{
                delete _writer;
                _writer = NULL;
            }
        }
        else if(_writer) _writer->set(this); 	// just append the trait
    }
    
    // phenotype
    choice = (int)get_parameter_value("quanti_save_phenotype");
    assert(!_phenotyper);
    if(choice) {
        _phenotyper = (dynamic_cast<TTraitQuantiProto*>(_popPtr->getFirstPrototype(_type)))->_phenotyper;
        if(!_phenotyper){
            _phenotyper = new TTQuantiFHvalue();
            if(_phenotyper->set(get_parameter("quanti_phenot_logtime"),
                                get_parameter("quanti_phenot_dir")->get_arg(),
                                get_parameter("quanti_phenot_filename")->get_arg(),
                                get_parameter("quanti_phenot_script")->get_arg(),
                                (int)get_parameter_value("quanti_phenot_sex"),
                                (int)get_parameter_value("quanti_phenot_age"),
                                choice,
                                this,
                                "phenotype",
                                ".phe", loader, get_popPtr())){
                _phenotyper->set_getter(1); // phenotype
                loader->attach(_phenotyper);
            }
            else{
                delete _phenotyper;
                _phenotyper = NULL;
            }
        }
        else{
            _phenotyper->set(this); 	// just append the trait
        }
        
    }
    
    //genotypic value
    choice = (int)get_parameter_value("quanti_save_geno_value");
    assert(!_genotyper);
    if(choice) {
        _genotyper = (dynamic_cast<TTraitQuantiProto*>(_popPtr->getFirstPrototype(_type)))->_genotyper;
        if(!_genotyper){
            _genotyper = new TTQuantiFHvalue();
            if(_genotyper->set(get_parameter("quanti_geno_value_logtime"),
                               get_parameter("quanti_geno_value_dir")->get_arg(),
                               get_parameter("quanti_geno_value_filename")->get_arg(),
                               get_parameter("quanti_geno_value_script")->get_arg(),
                               (int)get_parameter_value("quanti_geno_value_sex"),
                               (int)get_parameter_value("quanti_geno_value_age"),
                               choice,
                               this,
                               "genotypic value",
                               ".gen", loader, get_popPtr())){
                _genotyper->set_getter(0); // genotype
                loader->attach(_genotyper);
            }
            else{
                delete _genotyper;
                _genotyper = NULL;
            }
        }
        else if(_genotyper) _genotyper->set(this); 	// just append the trait
    }
}

// ----------------------------------------------------------------------------------------
// loadStatServices
// ----------------------------------------------------------------------------------------
void
TTraitQuantiProto::loadStatServices  (StatServices* loader)
{
    if(_stats)  delete _stats;
    _stats = new TTQuantiSH(this);
    loader->attach(_stats);
}

// ----------------------------------------------------------------------------------------
// clone
// ----------------------------------------------------------------------------------------
TTraitQuanti*
TTraitQuantiProto::hatch ()
{
    TTraitQuanti* new_trait = new TTraitQuanti();
    new_trait->set_from_prototype(this);
    return new_trait;
}

// ----------------------------------------------------------------------------------------
// check_allelic_array
// ----------------------------------------------------------------------------------------
/** check the mutation probabilities for the mutation model IMM */
void
TTraitQuantiProto::check_mutationValues_IMM()
{
    for(unsigned int l = 0; l < _nb_locus; ++l){
        if(abs(_mutationFreq[l][_nb_allele[l]/2]) > 1e-4){
            warning("The mutation model IMM requiers that the mutation probability of the allele with effect=0 is zero (locus %i, allele %i): Automatically adjusted!\n",
                    l+1, 1+_nb_allele[l]/2);
        }
    }
}

// ----------------------------------------------------------------------------------------
// check_allelic_array
// ----------------------------------------------------------------------------------------
/** check the allelic effect array for the mutation model IMM */
void
TTraitQuantiProto::check_allelicValues_IMM()
{
    
    unsigned int l, a;
    for(l=0; l<_nb_locus; ++l){
        // the number of alleles has to be odd
        if(!(_nb_allele[l]%2)) error("The mutation model IMM requires an odd number of alleles!\n");
        
        // alleles have to be regularely spaced
        double step = _allelicValues[l][1] - _allelicValues[l][0];
        for(a = 2; a < _nb_allele[l]; ++a){
            if(abs(_allelicValues[l][a] - _allelicValues[l][a-1] - step) > 1e-4){
                error("Allelic values: The effects of the alleles have to be regularly distributed if mutation model IMM is used(locus %i)!\n", l);
            }
        }
        
        // the effects have to be symmetric, i.e. the middle allele has to be 0
        if(abs(_allelicValues[l][_nb_allele[l]/2]) > 1e-4){
            warning("Allelic values: For the IMM mutation model the allelic effects are assumed to be symmetrically distributed (allele %i is not zero)!\n", 1+_nb_allele[l]/2);
        }
    }
}

// ----------------------------------------------------------------------------------------
// create_regular_spaced_array
// ----------------------------------------------------------------------------------------
/** create a regular stepped array from -range to +range (inlcuding the edges)
 ** due to round-off errors the array is made symmetrically starting in the middle
 */
void
TTraitQuantiProto::create_regular_spaced_array(double* array, const int& size, double half_range)
{
    assert(array);
    double step = 2.*half_range/(size-1.);    // compute step size
    double curVal;
    double *a_up, *a_down;
    double *end = array + size;
    
    if(size%2){    // odd  number of alleles
        curVal = 0;
        a_up = a_down = array + (int)size/2;
    }
    else{          // event number of alleles
        curVal = step/2.;
        a_up = array + (int)size/2;
        a_down = a_up - 1;
    }
    
    while(1){
        *a_up   = curVal;
        *a_down = -curVal;
        
        // next element  (done so complicate due to code guard memory leak detection)
        ++a_up;
        if(a_up == end) break;
        --a_down;
        curVal += step;
    }
}

// ----------------------------------------------------------------------------------------
// set_allelic_valuesFreq
// ----------------------------------------------------------------------------------------
/** compute the frequencies for the normal distribution N(0, sd) for the values
 * of array effect_array and store them in the array
 * Caution: they do not yet sum up to 1!!!
 */
void
TTraitQuantiProto::compute_frequencies(double* array, double* effect_array, const unsigned int& size, const double& sd)
{
    assert(array);
    assert(effect_array);
    double* end = array + size;
    for(; array != end; ++array, ++effect_array){
        *array = ProbDensNorm(*effect_array, 0, sd);
    }
}

// ----------------------------------------------------------------------------------------
// set_allelic_valuesFreq
// ----------------------------------------------------------------------------------------
/** read the allelic values and their frequencies from a file (used for mutation model IMM and RMM)
 */
void
TTraitQuantiProto::read_allele_file (string filename)
{
#ifdef _DEBUG
    message("  TTraitQuantiProto::read_allele_file (%s) ...",filename.c_str());
#endif
    
    // read the file
    TMatrix mat;
    map<string, int> fileInfo;
    bool hasLocus = true;
    unsigned int i, l, a, nbValues=0;
    
    // read the file
    fileInfo = mat.read_matrix(filename);
    unsigned int nbLines = mat.getNbRows();
    unsigned int nbCols  = mat.getNbCols();
    
    const unsigned int colsSize = 5;
    unsigned int cols[colsSize] = {my_NAN,my_NAN,my_NAN,my_NAN,my_NAN};  // 0: locus; 1: allele; 2: allelic_value; 3: mut_freq; 4: ini_freq;
    string colsName[colsSize] = {"col_locus","col_allele","col_mut_freq","col_ini_freq","col_allelic_value"};
    
    // get col values form file info
    map<string, int>::iterator pos;
    for(i=0; i<colsSize; ++i){
        pos = fileInfo.find(colsName[i]);
        if(pos != fileInfo.end()){
            cols[i] = pos->second - 1;
            fileInfo.erase(pos);
        }
    }
    
    // output all not used columns
    for(pos = fileInfo.begin(); pos != fileInfo.end(); ++pos){
        warning("Allelic file for trait '%s' (%s): Column '%s' (%i) is not a valid column names!\n",
                get_type_index().c_str(), filename.c_str(), pos->first.c_str(), pos->second);
    }
    
    // check if all cols are available and also dimensions of the matrix are met and create the arrays
    if(cols[0] == my_NAN) hasLocus = false; // error("Allelic file (%s): Locus column is missing!", filename.c_str());
    if(cols[1] == my_NAN) error("Allelic file for trait '%s' (%s): Allele column is missing!", get_type_index().c_str(), filename.c_str());
    if(cols[2] != my_NAN){    // mutation frequency
        ARRAY::create_2D(_mutationFreq, _nb_locus, _nb_allele, (double)my_NAN);
    }
    if(cols[3] != my_NAN){     // initial frequencies
        set_ini_allelicArray(&mat, cols[3]);  // _initAlleleFreqCols is set
    }
    else _initAlleleFreqCols = 0;
    if(cols[4] != my_NAN) {
        ARRAY::create_2D(_allelicValues, _nb_locus, _nb_allele, (double)my_NAN);
        _allelic_file = true;
    }
    
    
    for (i=0; i<colsSize; ++i){
        if(cols[i]!=my_NAN && cols[i]>nbCols) error("Allelic file of trait '%s' (%s): File info is wrong: matrix has only %i columns and not %i!\n", get_type_index().c_str(), filename.c_str(), nbCols, cols[i]);
    }
    
    // read the table: copy the allelic effects and their frequencies
    for(i=0; i<nbLines; ++i){
        if(hasLocus) l = (int)mat.get(i,cols[0])-1;    // if locus specific values
        else         l = 0;                            // if all loci have the same settings
        a = (int)mat.get(i,cols[1])-1;    // allele
        
        // test if out of range
        if(hasLocus && l>=_nb_locus){
            error("Allelic values: Locus %i (allele %i) is out of range (only %i loci specified)!\n",
                  l+1, a+1, _nb_locus);
        }
        if(a>=_nb_allele[l]){
            error("Allelic values: Allele %i of locus %i is out of range (only %i alleles specified)!\n",
                  a+1, l+1, _nb_allele[l]);
        }
        
        // check if the combination was already read
        if(   (_allelicValues  && _allelicValues[l][a]     != my_NAN)
           || (_mutationFreq   && _mutationFreq[l][a]      != my_NAN)
           || (_initAlleleFreq && _initAlleleFreq[0][l][a] != my_NAN)){
            
            error("Allelic values: Allelic value was already specified for locus %i and allele %i!\n", l+1, a+1);
        }
        
        // set the values
        if(hasLocus){         // locus specific settings
            set_allelicValues(&mat, i, l, a, cols);
            ++nbValues;
        }
        else {                // all loci have the same settings
            for(l = 0; l < _nb_locus; ++l){
                set_allelicValues(&mat, i, l, a, cols);
                ++nbValues;
            }
        }
    }
    
    // check if we have enough values
    unsigned int totValues = 0;
    for(l=0; l<_nb_locus; ++l){
        totValues += _nb_allele[l];
    }
    if(nbValues != totValues){
        error("Allelic file: %i of %i allelic values are not set by the file '%s'!\n",
              totValues-nbValues, totValues, filename.c_str());
    }
    
    // control the arrays and make them cumulative
    if(_allelic_file && _mut_model == IMM){
        check_allelicValues_IMM();
    }
    
    if(_mutationFreq){
        if(_mut_model == IMM){
            check_mutationValues_IMM();
        }
        for(unsigned int l = 0; l < _nb_locus; ++l){
            ARRAY::make_frequency(_mutationFreq[l], _nb_allele[l]);  // adjust the sum to 1
            ARRAY::cumulative(_mutationFreq[l], _nb_allele[l]);      // make it cumulative
        }
    }
    
    if(_initAlleleFreq){
        unsigned int c, l;
        for(l = 0; l < _nb_locus; ++l){
            for(c = 0; c < _initAlleleFreqCols; ++c) {
                ARRAY::make_frequency(_initAlleleFreq[c][l], _nb_allele[l]);  // adjust the sum to 1
                ARRAY::cumulative(_initAlleleFreq[c][l], _nb_allele[l]);      // make it cumulative
            }
        }
    }
    
#ifdef _DEBUG
    message("done!\n");
#endif
}

// ----------------------------------------------------------------------------------------
// set_mutationFreq
// ----------------------------------------------------------------------------------------
/** set the probabilities to mutate to a certain allele */
void
TTraitQuantiProto::set_mutationFreq(const string& trait)
{
    if(_mutationFreq) return;      // the file was already set explicitly by the allelic file
    
    assert(_allelicValues);        // they have to be set beforehand
    
    // if the variance of the alleic effects varies between loci
    Param* p = get_parameter("quanti_allelic_var"+trait);
    if(p->is_matrix()){     // each locus has an individual mutation freq array
        // get the values
        TMatrix* m = p->get_matrix();
        unsigned int nb = m->get_dims(NULL);
        double* var = m->get();
        assert(nb<=_nb_locus); // check the number of variances with the number of loci (already done before)
        
        // create the mutation freq matrix: each locus has other settings
        ARRAY::create_2D(_mutationFreq, _nb_locus, _nb_allele);
        for(unsigned int l=0; l<_nb_locus; ++l){
            set_mutationFreq_locus(l, var[l % nb]);
        }
        delete m;
    }
    else{	// all loci have the same settings -> all loci point to the same allele array (caution, when deleting)
        // check that all loci have the same nubmer of alleles
        for(unsigned int l=1; l<_nb_locus; ++l){
            if(_nb_allele[l]!=_nb_allele[l-1]) error("Loci have the same settings, thus the number of alleles must be identical!\n");
        }
        
        ARRAY::create_1D(_mutationFreq, _nb_locus);
        _mutationFreq[0] = new double[*_nb_allele];         // create the array for the first locus
        for(unsigned int l = 1; l < _nb_locus; ++l){
            _mutationFreq[l] = _mutationFreq[0];             // all loci point to this first locus array
        }
        set_mutationFreq_locus(0, p->get_value());         // only for the first locus
    }
}

// ----------------------------------------------------------------------------------------
// set_mutationFreq_locus
// ----------------------------------------------------------------------------------------
/** set the probabilities to mutate to a certain allele for locus l*/
void
TTraitQuantiProto::set_mutationFreq_locus(const unsigned int& l, const double& var)
{
    assert(var>0);
    compute_frequencies(_mutationFreq[l], _allelicValues[l], _nb_allele[l], sqrt(var));
    
    //if IMM the "middle allele" has to have a freq of zero!
    if(_mut_model == IMM) _mutationFreq[l][_nb_allele[l]/2] = 0;
    
    // check that the frequencies sum up to 1 and make the frequencies cumulative
    ARRAY::make_frequency(_mutationFreq[l], _nb_allele[l]);
    ARRAY::cumulative(_mutationFreq[l], _nb_allele[l]);
}


// ----------------------------------------------------------------------------------------
// set_initAlleleFreq
// ----------------------------------------------------------------------------------------
/** set the initial allelic values for a maximal "polymorphism"
 * _ini_allele_model is set already either to uniform or monomorph*/
void
TTraitQuantiProto::set_initAlleleFreq(const string& trait)
{
    // explicitly given values overwrite the settings
    if(_initAlleleFreq){
        _ini_allele_model = INI_DIST;
        return;
    }
    
    // if all populations are monomorph stop
    if(_ini_allele_model == INI_MONO) return;
    
    
    _ini_allele_model = INI_DIST;   
    // if the populations are maximal polymorph: _ini_allele_model = INI_DIST;
    // the allelic effects have to be set automatically
    Param* p = get_parameter("quanti_allelic_var"+trait);
    
    // if all loci have the same polymorphism
    if(p->is_matrix()){     // each locus has another setting
        // get the values
        TMatrix* m = p->get_matrix();
        unsigned int nb = m->get_dims(NULL);
        double* var = m->get();
        assert(nb<=_nb_locus); // check the number of variances with the number of loci (already done)
        
        // create the initial allele frequency matrix: each locus has other settings
        _initAlleleFreqCols = 1;      // all patches have the same initial allele frequencies
        ARRAY::create_3D(_initAlleleFreq, _initAlleleFreqCols, _nb_locus, _nb_allele); //_initAlleleFreq[patch][locus][allele]
        for(unsigned int l=0; l<_nb_locus; ++l){
            set_initAlleleFreq_locus(0, var[l % nb]);
        }
        delete m;
    }
    else{  // all loci have the same settings -> all loci point to the same allele array (caution, when deleting)
        // check that all loci have the same number of alleles
        for(unsigned int l=1; l<_nb_locus; ++l){
            if(_nb_allele[l]!=_nb_allele[l-1]) error("Loci have the same settings, thus the number of alleles must be identical!\n");
        }
        
        // the initial frequncies are the same as the mutation probabilities for the RMM
        _initAlleleFreqCols = 1;      // all patches have the same initial allele frequencies
        ARRAY::create_2D(_initAlleleFreq, _initAlleleFreqCols, _nb_locus);  // all are the same -> use only one locus array (caution, when deleting)
        _initAlleleFreq[0][0] = new double[*_nb_allele];   // create the array for the first locus
        for(unsigned int l = 1; l < _nb_locus; ++l){
            _initAlleleFreq[0][l] = _initAlleleFreq[0][0];       // all loci point to this first locus array
        }
        set_initAlleleFreq_locus(0, p->get_value());
    }
}

// ----------------------------------------------------------------------------------------
// set_initAlleleFreq_locus
// ----------------------------------------------------------------------------------------
/** set the initial allelic values for a maximal "polymorphism" of locus l
 * _ini_allele_model is set already either to uniform or monomorph*/
void
TTraitQuantiProto::set_initAlleleFreq_locus(const unsigned int& l, const double& var)
{
    // fill the array with values
    assert(var>0);
    compute_frequencies(_initAlleleFreq[0][l], _allelicValues[l], _nb_allele[l], sqrt(var));
    
    // check that the frequencies sum up to 1 and make the frequencies cumulative
    ARRAY::make_frequency(_initAlleleFreq[0][l], _nb_allele[l]);
    ARRAY::cumulative(_initAlleleFreq[0][l], _nb_allele[l]);
}

// ----------------------------------------------------------------------------------------
// set_fitnessFactor
// ----------------------------------------------------------------------------------------
void
TTraitQuantiProto::set_fitnessFactor(const string& trait, string name, double* &array)
{
    Param* p = get_parameter(name+trait);
    if(array) {delete[] array; array=NULL;}
    if(p->is_matrix()){      // locus specific settings
        // get the values
        TMatrix* m = p->get_matrix();
        unsigned int nb = m->get_dims(NULL);
        double* values = m->get();
        
        // check the number of values with the number of loci
        if(nb>_nb_locus) warning("There are more fitness factors ('%s') than loci defined! Only a part of the fitness factors is considered!\n", name.c_str());
        else if(_nb_locus % nb) warning("The number of fitness factors ('%s') is not an entire subset of the number of loci!\n", name.c_str());
        
        // create the array
        array = new double[_nb_locus];
        
        // get the values for each locus
        for(unsigned int l=0; l<_nb_locus; ++l){
            array[l] = values[l % nb];
        }
        delete m;
    }
    else{              	// common settings for all loci
        double value = p->get_value();
        if(value != strTo<double>(p->get_default_arg())){   // create the array only if it is not the default value
            array = new double[_nb_locus];
            for(unsigned int l=0; l<_nb_locus; ++l){
                array[l] = value;
            }
        }
    }
}

// ----------------------------------------------------------------------------------------
// set_dominaceValues
// ----------------------------------------------------------------------------------------
/** set the dominance values and/or the corresponding function pointers.
 * Before calling this function the dominance file has to be read...
 */
void
TTraitQuantiProto::set_dominanceValues(const string& trait)
{
    _dominance_mean = get_parameter_value("quanti_dominance_mean"+trait);
    _dominance_sd   = sqrt(get_parameter_value("quanti_dominance_var"+trait));
    
    
    if(_dominanceValues){				// dominance effects are set explicitly by the file (however not all must be set nby the file
        get_locus_genotype_func_ptr = &TTraitQuantiProto::get_locus_genotype_dominance_array;
    }
    else if(_dominance_sd){     // dominance effects are set by a distribution (redraw at each replicate!)
        ARRAY::create_3D(_dominanceValues, _nb_locus, _nb_allele, _nb_allele, (double)my_NAN);
        get_locus_genotype_func_ptr = &TTraitQuantiProto::get_locus_genotype_dominance_array;
    }
    else if(_dominance_mean) get_locus_genotype_func_ptr = &TTraitQuantiProto::get_locus_genotype_dominance_single;
    else                     get_locus_genotype_func_ptr = &TTraitQuantiProto::get_locus_genotype_additive;
    
    if(get_parameter_value("quanti_dominance_model"+trait)==1){
        get_genotype_dominace_func_ptr=&TTraitQuantiProto::get_genotype_dominance_h;
        if(!get_parameter_isSet("quanti_dominance_mean"+trait))_dominance_mean=0.5; // if not set codominant is 0.5 and not 0
    }
    else{
        get_genotype_dominace_func_ptr=&TTraitQuantiProto::get_genotype_dominance_k;
        if(_selection_model==4) error("Parameter 'quanti_dominance_model' has to be set to 1 when using selection coefficient!");
    }
}

// ----------------------------------------------------------------------------------------
// set_epistaticValues
// ----------------------------------------------------------------------------------------
/** set the epistatic values and/or the corresponding function pointers.
 * Before calling this function the epistatic file has to be read...
 */
void
TTraitQuantiProto::set_epistaticValues(const string& trait)
{
    _epistatic_sd = sqrt(get_parameter_value("quanti_epistatic_var"+trait));
    if(_phenoTree){                 // epistatic effects are set explicitly by the file
        get_genotype_func_ptr = &TTraitQuantiProto::get_genotype_epistatic;
    }
    else if(_epistatic_sd){          // epistatic effects set by its variance (reset at each replicate!)
        get_genotype_func_ptr = &TTraitQuantiProto::get_genotype_epistatic;
        if(!_phenoTree){
            unsigned int nb_allele_max = 0;             // get the max number of alleles
            for(unsigned int l=0; l<_nb_locus; ++l){
                if(nb_allele_max<_nb_allele[l]) nb_allele_max=_nb_allele[l];
            }
            _phenoTree = new Tree<unsigned char>(_nb_locus, nb_allele_max);
        }
    }
    else get_genotype_func_ptr = &TTraitQuantiProto::get_genotype_additive;   // no epistatic effects
}

// ----------------------------------------------------------------------------------------
// set_allelicValues
// ----------------------------------------------------------------------------------------
/** set the allelic effects automatically depending on the variance of the effects
 * the variance may be set globally or for each locus separately (by a matrix)
 * the allelic effects are regularly spaced (the frequency to mutate to the
 *    different alleles varies according to the allelic effect variance.
 * if the variance is identical for all loci only a single effect array is
 *    created and all loci point to this one
 * if the variance differs among loci each locus has its own allelic effect array
 *
 * prerequisite of this function are: _mut_model, _nb_allel, _allelicValues(if set by explicitly by the file)
 */
void
TTraitQuantiProto::set_allelicValues(const string& trait)
{
    // controls
    for(unsigned int l=0; l<_nb_locus; ++l){
        if(_nb_allele[l] < 2) error("A single allele makes no sense to the simulation!\n");
        if(_mut_model == IMM){  	// controls specific for the IMM model
            if(!(_nb_allele[l]%2)) error("The IMM mutation model requires an odd number of alleles!\n");
            if(_nb_allele[l] < 51) warning("The IMM mutation model makes no sense with less than 51 alleles per locus!\n");
        }
    }
    if(_selection_model==4){	// selection coefficient: 1 QTL with two alleles and no effects
        if(_nb_locus!=1) error("Selection coefficient requires exactly one QTL per trait!\n");
        if(_nb_allele[0]!=2) error("Selection coefficient requires two alleles per QTL!\n");
        ARRAY::create_1D(_allelicValues, _nb_locus);
        _allelicValues[0] = new double[*_nb_allele];    // create the allele array for the first locus
        _allelicValues[0][0] = 0; // allele A
        _allelicValues[0][1] = 1; // allele a
        return;
    }
    
    if(_allelicValues) return;  // the allelic effects are already set explicitly by the allelic file
    
    // the allelic effects have to be fixed automatically
    Param* p = get_parameter("quanti_allelic_var"+trait);
    
    // if the variance of the allelic effects varies between loci
    if(p->is_matrix()){
        // get the values
        TMatrix* m = p->get_matrix();
        unsigned int nb = m->get_dims(NULL);
        double* var = m->get();
        
        // check the number of variances with the number of loci
        if(nb>_nb_locus) warning("There are more allelic effect variances than loci defined! Only a part of the allelic effect variances is considered!\n");
        else if(_nb_locus % nb) warning("The number of allelic effect variances is not an entire subset of the number of loci!\n");
        
        // create the allelic file: each locus has other settings
        ARRAY::create_2D(_allelicValues, _nb_locus, _nb_allele);
        
        // get the variances for each locus
        for(unsigned int l=0; l<_nb_locus; ++l){
            set_allelicValues_locus(l, var[l % nb], trait);
        }
        delete m;
    }
    else{	// all loci have the same settings -> all loci point to the same allelic effect array (caution, when deleting)
        // check that all loci have the same nubmer of alleles
        for(unsigned int l=1; l<_nb_locus; ++l){
            if(_nb_allele[l]!=_nb_allele[l-1]) error("Loci have the same settings, thus the number of alleles must be identical!\n");
        }
        
        ARRAY::create_1D(_allelicValues, _nb_locus);
        _allelicValues[0] = new double[*_nb_allele];    // create the allele array for the first locus
        for(unsigned int l = 1; l < _nb_locus; ++l){
            _allelicValues[l] = _allelicValues[0];      // all loci point to this first allele array
        }
        set_allelicValues_locus(0, p->get_value(), trait);
    }
}

// ----------------------------------------------------------------------------------------
// set_allelicValues_locus
// ----------------------------------------------------------------------------------------
/** set the allelic effects of locus l
 */
void
TTraitQuantiProto::set_allelicValues_locus(const unsigned int& l, const double& var, const string& trait)
{
    if(var<=0) error("Parameter 'quanti_allelic_var%s' must be positive!\n", trait.c_str());
    double sd = sqrt(var);
    double range;
    if(_nb_allele[l] < 6) range = sd * (_nb_allele[l]-1);
    else range = sd * ((_mut_model == IMM) ? 20 : 6);
    create_regular_spaced_array(_allelicValues[0], _nb_allele[l], range);
}

// ----------------------------------------------------------------------------------------
// set_allelicValues
// ----------------------------------------------------------------------------------------
/** set the values for a single row
 * 	cols = {"col_locus","col_allele","col_mut_freq","col_ini_freq","col_allelic_value","col_allelic_value"};
 */
void
TTraitQuantiProto::set_allelicValues(TMatrix* mat, const unsigned int& i, const unsigned int& l, const unsigned int& a, unsigned int* cols)
{
    if(_allelicValues){
        assert(cols[4]!=my_NAN);
        _allelicValues[l][a]  = mat->get(i,cols[4]); // allelic value
    }
    TTraitProto::set_allelicValues(mat, i, l, a, cols);
}

// ----------------------------------------------------------------------------------------
// read_locus_file
// ----------------------------------------------------------------------------------------
/** read dominance values and fitness factors from a file
 * not set combinations will be recessive
 */
void
TTraitQuantiProto::read_locus_file(string filename)
{
#ifdef _DEBUG
    message("  TTraitQuantiProto::read_locus_file (%s) ...",filename.c_str());
#endif
    unsigned int i, l, a1, a2;
    bool hasLocus = true;
    
    // read the file
    TMatrix mat;
    map<string, int> fileInfo = mat.read_matrix(filename);      // read the file
    unsigned int nbLines = mat.getNbRows();
    unsigned int nbCols  = mat.getNbCols();
    
    const unsigned int colsSize=5;
    unsigned int cols[5] = {my_NAN,my_NAN,my_NAN,my_NAN};  // 0: locus; 1: allele1; 2: allele2; 3: value; 4: fitness_factor
    string colsName[5]   = {"col_locus","col_allele1","col_allele2","col_dominance","col_fitness_factor"};
    
    // get col values form file info
    map<string, int>::iterator pos;
    for(i=0; i<colsSize; ++i){
        pos = fileInfo.find(colsName[i]);
        if(pos != fileInfo.end()){
            cols[i] = pos->second - 1;
            fileInfo.erase(pos);
        }
    }
    
    // ouput all not used columns
    for(pos = fileInfo.begin(); pos != fileInfo.end(); ++pos){
        warning("Dominance file for trait '%s' (%s): Column '%s' (%i) is not a valid column names!\n",
                get_type_index().c_str(), filename.c_str(), pos->first.c_str(), pos->second);
    }
    
    // check if all cols are available and also dimensions of the matrix are met adn create the arrays
    if(cols[0] == my_NAN) hasLocus = false; // error("Dominance file (%s): Locus column is missing!", filename.c_str());
    if(cols[1] == my_NAN) error("Dominace file (%s): 1. allele column is missing!", filename.c_str());
    if(cols[2] == my_NAN) error("Dominace file (%s): 2. allele column is missing!", filename.c_str());
    if(cols[3] != my_NAN) {  // dominance values
        ARRAY::create_3D(_dominanceValues, _nb_locus, _nb_allele, _nb_allele, (double) my_NAN);
    }
    if(cols[4] != my_NAN) {  // fitness factor
        ARRAY::create_3D(_fitnessFactor, _nb_locus, _nb_allele, _nb_allele, (double) my_NAN);
    }
    
    for (i=0; i<colsSize; ++i){
        if(cols[i]!=my_NAN && cols[i]>nbCols) error("Dominance file (%s): File info is wrong: matrix has only %i columns and not %i!\n", filename.c_str(), nbCols, cols[i]);
    }
    
    // read the table: copy the dominance values
    for(i=0; i<nbLines; ++i){
        if(hasLocus) l  = (int)mat.get(i,cols[0])-1;    // if locus specific values
        else         l = 0;                             // if all loci have the same settings
        a1 = (int)mat.get(i,cols[1])-1;    // allele 1
        a2 = (int)mat.get(i,cols[2])-1;    // allele 2
        
        // test if out of range
        if(l>=_nb_locus){
            error("Dominance values: Locus %i is out of range (only %i loci specified)!\n",
                  l+1, _nb_locus);
        }
        if(a1>=_nb_allele[l]){
            error("Dominance values: Allele %i of locus %i is out of range (only %i alleles specified)!\n",
                  a1+1, l+1, _nb_allele[l]);
        }
        if(a2>=_nb_allele[l]){
            error("Dominance values: Allele %i of locus %i is out of range (only %i alleles specified)!\n",
                  a2+1, l+1, _nb_allele[l]);
        }
        if(a1 > a2) swap(a1, a2);  // wrong order of the alleles: -> swap
        
        // check if the combination was already read
        if(   (_dominanceValues && _dominanceValues[l][a1][a2] != my_NAN)
           || (_fitnessFactor   && _fitnessFactor[l][a1][a2]   != my_NAN)){
            
            error("Dominance values: Value was already specified for locus %i and alleles %i and %i!\n", l+1, a1+1, a2+1);
        }
        
        // set the values
        if(hasLocus){       // if all loci have the same settings
            if(_dominanceValues) _dominanceValues[l][a1][a2] = mat.get(i,cols[3]);
            if(_fitnessFactor)   _fitnessFactor[l][a1][a2]   = mat.get(i,cols[4]);
        }
        else{               // all loci have the same settings
            for(l=0; l<_nb_locus; ++l){
                if(_dominanceValues) _dominanceValues[l][a1][a2] = mat.get(i,cols[3]);
                if(_fitnessFactor)   _fitnessFactor[l][a1][a2]   = mat.get(i,cols[4]);
            }
        }
    }
    
#ifdef _DEBUG
    message("done!\n");
#endif
}

// ----------------------------------------------------------------------------------------
// read_genome_file
// ----------------------------------------------------------------------------------------
/** read epistatic values from a file */
void
TTraitQuantiProto::read_genome_file(string filename)
{
#ifdef _DEBUG
    message("  TTraitQuantiProto::read_genome_file (%s) ...",filename.c_str());
#endif
    
    try{
        // read the file
        TMatrix mat;
        map<string, int> fileInfo;
        int i;
        unsigned int nbValues=0;
        
        // read the file
        fileInfo = mat.read_matrix(filename, 1);
        int nbLines = mat.getNbRows();
        int nbCols  = mat.getNbCols();
        map<unsigned int, map<unsigned int, string> >* pStrMatrix = mat.get_strMatrix();
        
        const int nbParam = 4;
        int colsSize=nbParam;
        int cols[nbParam]        = {my_NAN,my_NAN,my_NAN,my_NAN};  // 0: genotype; 1: epistaticVal; 2: genotypicVal; 3: fitnessFactor
        string colsName[nbParam] = {"col_genotype","col_epistatic_value","col_genotypic_value", "col_fitness_factor"};
        
        // get col values form file info
        map<string, int>::iterator pos;
        for(i=0; i<colsSize; ++i){
            pos = fileInfo.find(colsName[i]);
            if(pos != fileInfo.end()){
                cols[i] = pos->second - 1;
                fileInfo.erase(pos);        // the element is correct remove it from the list
            }
        }
        
        // ouput all not used columns
        for(pos = fileInfo.begin(); pos != fileInfo.end(); ++pos){
            warning("Epistatic file for trait '%s' (%s): Column '%s' (%i) will not be considered!\n",
                    get_type_index().c_str(), filename.c_str(), pos->first.c_str(), pos->second);
        }
        
        // check if all cols are available and also dimensions of the matrix are met and create the arrays
        if(cols[0] == my_NAN) error("Epistatic file of trait '%s' (%s): Genotype column is missing!", get_type_index().c_str(), filename.c_str());
        
        if(cols[1] != my_NAN || cols[2] != my_NAN){      // epistatic or genotypic value
            if(_phenoTree)  delete _phenoTree;
            unsigned int nb_allele_max = 0;             // get the max number of alleles
            for(unsigned int l=0; l<_nb_locus; ++l){
                if(nb_allele_max<_nb_allele[l]) nb_allele_max=_nb_allele[l];
            }
            _phenoTree = new Tree<unsigned char>(_nb_locus, nb_allele_max);
        }
        else if(_phenoTree){delete _phenoTree; _phenoTree=NULL;}
        
        if(cols[1] != my_NAN && cols[2] != my_NAN){
            warning("Epistatic file of trait '%s' (%s): Both genotypic and epistatic values are present: only the genotypic values are considered!", get_type_index().c_str(), filename.c_str());
            cols[1] = my_NAN;
        }
        
        if(cols[2] != my_NAN){
            // allelic and dominance values are not used: remove them
            _allelic_file = _dominance_file = false;
            delete_allelicValues();
            delete_dominanceValues();
        }
        
        if(cols[3] != my_NAN){        // fitness factor
            if(_fitnessFactorTree)  delete _fitnessFactorTree;
            unsigned int nb_allele_max = 0;             // get the max number of alleles
            for(unsigned int l=0; l<_nb_locus; ++l){
                if(nb_allele_max<_nb_allele[l]) nb_allele_max=_nb_allele[l];
            }
            _fitnessFactorTree = new Tree<unsigned char>(_nb_locus, nb_allele_max);
        }
        else if(_fitnessFactorTree){delete _fitnessFactorTree; _fitnessFactorTree=NULL;}
        
        for (i=0; i<colsSize; ++i){
            if(cols[i] != my_NAN && cols[i] > nbCols) error("Epistatic file of trait '%s' (%s): File info is wrong: matrix has only %i columns and not %i!\n", get_type_index().c_str(), filename.c_str(), nbCols, cols[i]);
        }
        
        
        unsigned int length, digit=my_NAN;       // number of characters representing an allele
        unsigned char** genotype = ARRAY::new_2D<unsigned char>(_nb_locus, ploidy);
        double value;
        string text, t;
        unsigned int locus;
        
        // get the values
        for(i=0; i<nbLines; ++i){
            // get the genotype
            istringstream LINE;       // make a new stream
            if(mat.get(i,cols[0]) != my_NAN) error("Could not read epistatic file of trait '%s' (%s): Genotype must be written within brackets!\n", get_type_index().c_str(), filename.c_str());
            assert((*pStrMatrix)[i].find(cols[0]) != (*pStrMatrix)[i].end());
            text = (*pStrMatrix)[i].find(cols[0])->second;
            text = text.substr(1, text.length()-2); // remove the start and ending bracket
            LINE.str(text);                         // allocate the matrix to the new stream
            
            // for each genotype
            locus = 0;
            do{
                LINE >> text;
                length = (unsigned int)text.length();
                if(length != 2*digit){
                    digit = length/2;
                    if(length != 2*digit) error("Could not read epistatic file of trait '%s' (%s): Genotype could not be read at line %i!\n", get_type_index().c_str(), filename.c_str(), i+1);
                }
                
                if(locus>=_nb_locus) error("Could not read epistatic file of trait '%s' (%s): The genotype has too many loci (%u loci expected)!\n", get_type_index().c_str(), filename.c_str(), _nb_locus);
                
                try{
                    // get first allele of locus i
                    genotype[locus][0] = (unsigned char) (strTo<unsigned int>(text.substr(0, digit))-1);
                    if((unsigned int)genotype[locus][0] >= (unsigned int)_nb_allele[locus]){
                        error("Epistatic values: Allele 1 of locus %i (line %i) is out of range (only %i alleles specified)!\n", locus+1, i+1, _nb_allele);
                    }
                    
                    // get second allele of locus i
                    genotype[locus][1] = (unsigned char) (strTo<unsigned int>(text.substr(digit, digit))-1);
                    if((unsigned int)genotype[locus][1] >= (unsigned int)_nb_allele[locus]){
                        error("Epistatic values: Allele 2 of locus %i (line %i) is out of range (only %i alleles specified)!\n", locus+1, i+1, _nb_allele);
                    }
                }
                catch(...){
                    error("File '%s' could not be read at line %i: locus %i is not a number!\n", filename.c_str(), i+1, locus+1);
                }
                ++locus;
            }while(!LINE.eof());
            
            // control if the genotype has the correct number of loci
            if(locus!=_nb_locus) error("File '%s' could not be read at line %i: wrong number of loci specified (%i instead of %i!\n", filename.c_str(), i+1, locus+1, _nb_locus);
            
            // control if the genotype was already set
            
            // set the genotype
            if(cols[1] != my_NAN){ // set the epistatic value
                if(_phenoTree->get_value(genotype) != my_NAN){
                    error("Epistatic file %s: value for genotype of line %i was already specified!\n", filename.c_str(), i+1);
                }
                value = mat.get(i,cols[1]);
                if(value == my_NAN) error("Epistatic file '%s': epistatic values are not available!\n", filename.c_str());
                _phenoTree->set_value(genotype, get_genotype_additive(genotype)+value);
            }
            if(cols[2] != my_NAN){ // set directly the genotypic value
                if(_phenoTree->get_value(genotype) != my_NAN){
                    error("Epistatic file %s: value for genotype of line %i was already specified!\n", filename.c_str(), i+1);
                }
                value = mat.get(i,cols[2]);
                if(value == my_NAN) error("Epistatic file '%s': genotypic values are not available!\n", filename.c_str());
                _phenoTree->set_value(genotype, value);
            }
            if(cols[3] != my_NAN){ // set the fitness factor
                if(_fitnessFactorTree->get_value(genotype) != my_NAN){
                    error("Epistatic file %s: value for genotype of line %i was already specified!\n", filename.c_str(), i+1);
                }
                value = mat.get(i,cols[3]);
                if(value != my_NAN) _fitnessFactorTree->set_value(genotype, value);
            }
            ++nbValues;
        }
        
        // check if we have enough values
        unsigned int totValues = 1;
        for(unsigned int l=0; l<_nb_locus; ++l){
            totValues *= _nb_allele[l]*(_nb_allele[l]+1) / 2;
        }
        if(nbValues != totValues && (cols[1] != my_NAN || cols[2] != my_NAN)){
            error("Epistatic values: %i of %i epistatic values are not set by the file '%s'!\n",
                  totValues-nbValues, totValues, filename.c_str());
        }
        
        ARRAY::delete_2D(genotype, _nb_locus);
        
    }catch(...) {error("Was not able to read file '%s'!\n", filename.c_str());}
    
#ifdef _DEBUG
    message("done!\n");
#endif
}

// ----------------------------------------------------------------------------------------
// print_allelic_values
// ----------------------------------------------------------------------------------------
/** print only the allelic values which were used */
void
TTraitQuantiProto::print_allelic_values(string name)
{
    if(!_allelicValues) return;  // if not used do not print it
    
    string rpl = _popPtr->getReplicateCounter();
    unsigned int p;
    
    // create the filename
    string filename = name
    + get_trait_indexStr_t()
    + _popPtr->getReplicateCounter_r()
    + ".txt";
#ifdef _DEBUG
    message("TTraitQuantiProto::print_allelic_values (%s)\n",filename.c_str());
#endif
    
    ofstream FILE(filename.c_str());
    
    FILE.width(12);
    FILE.setf(ios::left,ios::adjustfield);
    FILE << setprecision(4);
    
    // write title
    FILE << "#Allelic values ";
    if(get_trait_index() || !rpl.empty()){
        FILE << "(";
        if(get_trait_index())                 FILE << "Trait: " << get_trait_indexStr();
        if(get_trait_index() && !rpl.empty()) FILE << ", ";
        if(!rpl.empty())                      FILE << "Replicate: " << rpl;
        FILE << ")";
    }
    FILE << "\n################################################################\n";
    
    // write file info
    int col = 1;
    FILE << "\n[FILE_INFO]{";
    FILE << "\n    col_locus "         << col++;
    FILE << "\n    col_allele "        << col++;
    FILE << "\n    col_allelic_value " << col++;
    FILE << "\n    col_mut_freq "      << col++;
    FILE << "\n    col_ini_freq "      << col;
    FILE << "\n}";
    
    FILE << "\n\n# " << _nb_locus  << " loci"
    << "\n# max "   << get_nb_allele_max() << " alleles per locus" ;
    
    switch(_mut_model){
        case KAM: FILE << "\n# mutation model: KAM"; break;
        case SSM: FILE << "\n# mutation model: SSM"; break;
        case RMM: FILE << "\n# mutation model: RMM"; break;
        case IMM: FILE << "\n# mutation model: IMM"; break;
        case NO:  FILE << "\n# mutation model: NO"; break;
    }
    
    
    // write heading
    FILE << "\n\n#locus" << "\tallele" << "\tvalue" << "\tmut_freq" << "\tini_freq";
    
    // write the values line by line
    unsigned int l, a;
    for(l=0; l<_nb_locus; ++l){
        for(a=0; a<_nb_allele[l]; ++a){
            if(_allelicValues[l][a] != my_NAN){
                FILE << "\n" << (l+1)        // locus
                << "\t" << (a+1)        // allele
                << "\t" << _allelicValues[l][a]; // allelic value
                
                if(_mutationFreq){
                    if(a) FILE << "\t" << (_mutationFreq[l][a]-_mutationFreq[l][a-1]);  // the array is cumulative!
                    else  FILE << "\t" << _mutationFreq[l][a];
                }
                else    FILE << "\t" << "0"; 		// mutation rate is zero
                
                if(_initAlleleFreq){
                    if(_initAlleleFreqCols == 1){
                        if(a) FILE << "\t" << (_initAlleleFreq[0][l][a]-_initAlleleFreq[0][l][a-1]);  // the array is cumulative!
                        else  FILE << "\t" << _initAlleleFreq[0][l][a];
                    }
                    else{
                        FILE << "\t{";
                        for(p=0; p<_initAlleleFreqCols; ++p){
                            if(p) FILE << "\t";
                            if(a) FILE << (_initAlleleFreq[p][l][a]-_initAlleleFreq[p][l][a-1]);  // the array is cumulative!
                            else  FILE << _initAlleleFreq[p][l][a];
                        }
                        FILE << "}";
                    }
                }
                else FILE << "\t" << ((_nb_allele[l]/2)== a ? "1" : "0"); // monomorph initial populations
            }
        }
    }
    FILE.close();
}

// ----------------------------------------------------------------------------------------
// print_dominance_values
// ----------------------------------------------------------------------------------------
/** print only the dominance values which were used */
void
TTraitQuantiProto::print_dominance_values(string name)
{
    if(!(_dominanceValues || _fitnessFactor)) return;   // if not used don't make the output
    
    string rpl = _popPtr->getReplicateCounter();
    
    // create the filename
    string filename = name
    + get_trait_indexStr_t()
    + _popPtr->getReplicateCounter_r()
    + ".txt";
    
#ifdef _DEBUG
    message("TTraitQuantiProto::print_dominance_values (%s)\n",filename.c_str());
#endif
    
    ofstream FILE(filename.c_str());
    
    // write title
    FILE << "#Dominance values";
    if(get_trait_index() || !rpl.empty()){
        FILE << " (";
        if(get_trait_index())                 FILE << "Trait: " << get_trait_indexStr();
        if(get_trait_index() && !rpl.empty()) FILE << ", ";
        if(!rpl.empty())                      FILE << "Replicate: " << rpl;
        FILE << ")";
    }
    FILE << "\n################################################################\n";
    
    // write file info
    int col = 1;
    FILE << "\n[FILE_INFO]{";
    FILE << "\n    col_locus " << col++;
    FILE << "\n    col_allele1 " << col++;
    FILE << "\n    col_allele2 " << col++;
    if(_dominanceValues) FILE << "\n    col_dominance " << col++;
    if(_fitnessFactor)   FILE << "\n    col_fitness_factor " << col;
    FILE << "\n}";
    
    FILE << "\n\n# " << _nb_locus  << " loci"
    << "\n# max "   << get_nb_allele_max() << " alleles per locus"
    << "\n# vg(l) = a1 + a2 + k(a2-a1)   with a1 < a2";
    
    // write heading
    FILE << "\n\n#locus"
    << "\tall_1"
    << "\tall_2";
    if(_dominanceValues) FILE << "\tvalue";
    if(_fitnessFactor)   FILE << "\tfitness";
    
    // write the values line by line
    unsigned int l, a1, a2;
    for(l=0; l<_nb_locus; ++l){
        for(a1=0; a1<_nb_allele[l]; ++a1){
            for(a2=a1; a2<_nb_allele[l]; ++a2){
                // are there values present?
                if(   !(_dominanceValues && _dominanceValues[l][a1][a2]!=my_NAN)
                   && !(_fitnessFactor   && _fitnessFactor[l][a1][a2]  !=my_NAN)) continue;
                
                // write the allele info
                FILE << "\n" << (l+1)                             // locus
                << "\t" << (a1+1)                            // allele 1
                << "\t" << (a2+1);                           // allele 2
                
                // write dominance value if needed
                if(_dominanceValues){
                    if(_dominanceValues[l][a1][a2]!=my_NAN) FILE << "\t" << _dominanceValues[l][a1][a2];
                    else                                    FILE << "\t" << my_NANstr;
                }
                
                // write fitness factor if needed
                if(_fitnessFactor){
                    if(_fitnessFactor[l][a1][a2]!=my_NAN)   FILE << "\t" << _fitnessFactor[l][a1][a2];
                    else                                    FILE << "\t" << my_NANstr;
                }
            }
        }
    }
    
    FILE << "\n";
    FILE.close();
}

// ----------------------------------------------------------------------------------------
// print_epistatic_values
// ----------------------------------------------------------------------------------------
/** print only the epistatic values which were used */
void
TTraitQuantiProto::print_epistatic_values(string name)
{
    
    if(!(_phenoTree || _fitnessFactorTree)) return;   // if not used don't make the output
    
    string rpl = _popPtr->getReplicateCounter();
    
    // create the filename
    string filename = name
    + get_trait_indexStr_t()
    + _popPtr->getReplicateCounter_r()
    + ".txt";
    
#ifdef _DEBUG
    message("TTraitQuantiProto::print_epistatic_values (%s)\n",filename.c_str());
#endif
    
    ofstream FILE(filename.c_str());
    
    unsigned int digit = getNbDigits(get_nb_allele_max());
    
    // write title
    FILE << "#Epistatic values";
    if(get_trait_index() || !rpl.empty()){
        FILE << " (";
        if(get_trait_index())                 FILE << "Trait: " << get_trait_indexStr();
        if(get_trait_index() && !rpl.empty()) FILE << ", ";
        if(!rpl.empty())                      FILE << "Replicate: " << rpl;
        FILE << ")";
    }
    FILE << "\n################################################################\n";
    
    // write file info
    int col = 1;
    FILE << "\n[FILE_INFO]{";
    FILE << "\n    col_genotype " << col++;
    if(_phenoTree){
        if(_allelicValues){
            FILE << "\n    col_epistatic_value " << col++;
            FILE << "\n    # col_genotypic_value " << col++;
        }
        else FILE << "\n    col_genotypic_value " << col++;
    }
    if(_fitnessFactorTree)   FILE << "\n    col_fitness_factor " << col;
    FILE << "\n}";
    
    FILE << "\n\n# " << _nb_locus  << " loci"
    << "\n# max "   << get_nb_allele_max() << " alleles per locus";
    
    // write heading
    FILE << "\n\n#genotype";
    if(_phenoTree){
        if(_allelicValues) FILE << "\tepistaticVal";
        FILE << "\tgenotypicVal";
    }
    if(_fitnessFactorTree) FILE <<"\tfitnessFactor";
    
    // create the first possible genotype to explore all possible genotypes
    unsigned char** seq = ARRAY::new_2D<unsigned char>(_nb_locus, ploidy, (unsigned char) 0);
    
    double val_pheno, val_fitness;
    do{
        // get the values, i.e. check if the genotype is set
        if(_phenoTree)         val_pheno   = _phenoTree->get_value(seq);
        else                   val_pheno   = my_NAN;
        if(_fitnessFactorTree) val_fitness = _fitnessFactorTree->get_value(seq);
        else                   val_fitness = my_NAN;
        
        // if genotype is not set do not plot it
        if(val_pheno != my_NAN || val_fitness != my_NAN){
            // print the genotype
            FILE << "\n{";
            print_gentoype(FILE, seq, digit);
            FILE << "}";
            
            // print the values
            if(_phenoTree){
                if(_allelicValues) FILE << "\t" << val_pheno - get_genotype_additive(seq);
                FILE << "\t" << val_pheno;
            }
            if(_fitnessFactorTree) FILE << "\t" << val_fitness;
        }
    }while(get_next_gentoype(seq));      // get the next genotype
    
    FILE.close();
    
    ARRAY::delete_2D(seq, _nb_locus);
}

//----------------------------------------------------------------------------------------
// operator=
// ----------------------------------------------------------------------------------------
void
TTraitQuantiProto::print_gentoype(ostream& FILE, unsigned char** seq, const unsigned int& digit)
{
    unsigned int l, a;
    for(l=0; l<_nb_locus; ++l){
        if(l) FILE << " ";
        for(a=0; a<ploidy; ++a){
            FILE.fill('0');
            FILE.width(digit);
            FILE<<(unsigned int)(seq[l][a]+1);
        }
    }
}

//----------------------------------------------------------------------------------------
// operator=
// ----------------------------------------------------------------------------------------
/** generates the next possible genotype (seq is altered)
 * returns true if this was possible and false if the last possible genotype is reached
 */
bool
TTraitQuantiProto::get_next_gentoype(unsigned char** seq)
{
    int l;
    
    for(l=(int)_nb_locus-1; l>=0; --l){
        // increment the last allele and check if it is valable
        ++seq[l][1];
        if((unsigned int)seq[l][1] < (unsigned int)_nb_allele[l]) return true;
        
        // increment the second allele and reset the least allele to the same one (note 0102 = 0201!!!)
        ++seq[l][0];
        if((unsigned int)seq[l][0] < (unsigned int)_nb_allele[l]){
            seq[l][1] = seq[l][0];
            return true;
        }
        
        // the locus has to be switched: reset both alleles to zero
        seq[l][0] = seq[l][1] = 0;
    }
    return false;
}

//----------------------------------------------------------------------------------------
// getAllelicValue
// ----------------------------------------------------------------------------------------
double
TTraitQuantiProto::getAllelicValue(const unsigned int& l, const unsigned char& i)
{
    assert(_allelicValues[l][i] != my_NAN);
    return _allelicValues[l][i];
}

//----------------------------------------------------------------------------------------
// getDominanceValue
// ----------------------------------------------------------------------------------------
double
TTraitQuantiProto::getDominanceValue(const unsigned int& l, const unsigned char& a1, const unsigned char& a2)
{
    double *value;
    assert(a1<a2);
    value = &_dominanceValues[l][a1][a2];
    if(*value == my_NAN) *value = _popPtr->rand().Normal(_dominance_mean, _dominance_sd);  //get the value if not yet set
    return *value;
}

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                             ******** Prototype ********

// ----------------------------------------------------------------------------------------
// get_locus_genotype
// ----------------------------------------------------------------------------------------
double
TTraitQuantiProto::get_locus_genotype_additive(const unsigned int& l, const unsigned char& a1, const unsigned char& a2)
{
    return (getAllelicValue(l, a1) + getAllelicValue(l, a2));
}

double
TTraitQuantiProto::get_locus_genotype_dominance_array(const unsigned int& l, const unsigned char& a1, const unsigned char& a2){
    double a_1 = getAllelicValue(l, a1);    // get effect of allele 1
    double a_2 = getAllelicValue(l, a2);    // get effect of allele 2
    if(a1<a2) return get_genotype_dominance(a_1, a_2, getDominanceValue(l, a1, a2));
    if(a1>a2) return get_genotype_dominance(a_2, a_1, getDominanceValue(l, a2, a1));
    return a_1 + a_2;
}

double
TTraitQuantiProto::get_locus_genotype_dominance_single(const unsigned int& l, const unsigned char& a1, const unsigned char& a2)
{
    return get_genotype_dominance(getAllelicValue(l, a1), getAllelicValue(l, a2), _dominance_mean);
}


/* G = a1 + a2 + k(a2-a1)       // a1 < a2
 * k < -1:  underdominant
 * k = -1:  the smaller allele (a1) is dominant
 * k = 0:   purely additive
 * k = 1:   the bigger allele (a2) is dominant
 * k > 1:   overdominance
 */
double
TTraitQuantiProto::get_genotype_dominance_h(double a1, double a2, double h)
{
    if(a1<a2) return (1-h)*a1 + h*a2;
    return (1-h)*a2 + h*a1;
}

/* G = (1-h)a1 + ha2       // a1 < a2
 * h < 0:   underdominant
 * h = 0:   the smaller allele (a1) is dominant
 * h = 0.5: purely additive
 * h = 1:   the bigger allele (a2) is dominant
 * h > 1:   overdominance
 */
double
TTraitQuantiProto::get_genotype_dominance_k(double a1, double a2, double k)
{
    return a1 + a2 + k*abs(a2-a1);
}

// ----------------------------------------------------------------------------------------
// get_genotype
// ----------------------------------------------------------------------------------------
double
TTraitQuantiProto::get_genotype_additive(unsigned char** seq)
{
    double sum=0;
    for(unsigned int l=0; l<_nb_locus; ++l){
        sum += (this->*get_locus_genotype_func_ptr)(l, seq[l][0], seq[l][1]);
    }
    return sum;
}

double
TTraitQuantiProto::get_genotype_epistatic(unsigned char** seq)
{
    double val  = _phenoTree->get_value(seq);
    if(val != my_NAN) return val;
    
    // if not yet computed compute it
    val = get_genotype_additive(seq);
    val += get_popPtr()->rand().Normal(0, _epistatic_sd);
    _phenoTree->set_value(seq, val);
    return val;
}

// ----------------------------------------------------------------------------------------
// get_fitnessFactor_genome
// ----------------------------------------------------------------------------------------
/** fitness factor specified explicitly for each genome */
double
TTraitQuantiProto::get_fitnessFactor_genome(unsigned char** seq)
{
    double val  = _fitnessFactorTree->get_value(seq);
    if(val != my_NAN) return val;
    
    // if not set for the genome get the explicit fitness factor
    if(_fitnessFactor) return get_fitnessFactor_locus(seq);
    return get_fitnessFactor_global(seq);
}

// ----------------------------------------------------------------------------------------
// get_fitnessFactor_locus
// ----------------------------------------------------------------------------------------
/** fitness factor specified at the locus level */
double
TTraitQuantiProto::get_fitnessFactor_locus(unsigned char** seq)
{
    assert(_fitnessFactor);
    
    double product=1, value;
    unsigned char a1, a2;
    for(unsigned int l=0; l<_nb_locus; ++l){
        a1 = seq[l][0];     // get allele 1
        a2 = seq[l][1];     // get allele 2
        if(a1 > a2) value = _fitnessFactor[l][a2][a1];
        else        value = _fitnessFactor[l][a1][a2];
        
        if(value != my_NAN){ // value is explicitly defined
            product *= value;
        }
        else{                // generally defined
            if(a1==a2) product *= _fitnessFactor_homozygote   ? _fitnessFactor_homozygote[l]   : 1;
            else       product *= _fitnessFactor_heterozygote ? _fitnessFactor_heterozygote[l] : 1;
        }
    }
    return product;
}

// ----------------------------------------------------------------------------------------
// get_fitnessFactor_global
// ----------------------------------------------------------------------------------------
/** fitness factor specifed globally for heterozygote/homozygote loci */
double
TTraitQuantiProto::get_fitnessFactor_global(unsigned char** seq)
{
    double product=1;
    for(unsigned int l=0; l<_nb_locus; ++l){
        if(seq[l][0]==seq[l][1]) product *= _fitnessFactor_homozygote   ? _fitnessFactor_homozygote[l]   : 1;
        else                     product *= _fitnessFactor_heterozygote ? _fitnessFactor_heterozygote[l] : 1;
    }
    return product;
}

// ----------------------------------------------------------------------------------------
// get_fitnessFactor_locus_freqDepend
// ----------------------------------------------------------------------------------------
/** fitness factor for frequency dependent selection
 * Noe: _locusFreqs have to be recomputed at each generation and patch
 */
double
TTraitQuantiProto::get_fitnessFactor_freqDepend(unsigned char** seq)
{
    assert(_locusFreqs);
    
    // first get it without frequency dependent selection
    double product = get_fitnessFactor2_func_ptr ? (this->*get_fitnessFactor2_func_ptr)(seq) : 1;
    
    // and now add the frequency depend part on it (since it is multiplicate the order does not matter)
    double freqDependFactor;
    unsigned char a1, a2;
    for(unsigned int l=0; l<_nb_locus; ++l){
        freqDependFactor = _fitnessFactor_freqDepend[l];
        if(!freqDependFactor) continue; // if zero it is not set
        
        a1 = seq[l][0];     // get allele 1
        a2 = seq[l][1];     // get allele 2
        
        if(a1 < a2){
            assert(_locusFreqs[l].find(a1)!=_locusFreqs[l].end());
            assert(_locusFreqs[l][a1].find(a2)!=_locusFreqs[l][a1].end());
            product *= 1-freqDependFactor*_locusFreqs[l][a1][a2];
        }
        else{
            assert(_locusFreqs[l].find(a2)!=_locusFreqs[l].end());
            assert(_locusFreqs[l][a2].find(a1)!=_locusFreqs[l][a2].end());
            product *= 1-freqDependFactor*_locusFreqs[l][a2][a1];
        }
    }
    return product;
}

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                 ******** Phenotyper / Genotyper ********

// ----------------------------------------------------------------------------------------
// FHwrite
// ----------------------------------------------------------------------------------------
void
TTQuantiFHvalue::FHwrite ()
{
    
    // get the number of occupied patches
    if(!_popPtr->get_nbSamplePatch()) return;  // if all patches are extinct
    
    string filename = get_path() + getGenerationReplicateFileName() + get_extension();
    
#ifdef _DEBUG
    message("TTQuantiFHvalue::FHwrite (%s)\n",filename.c_str());
#endif
    
    ofstream FILE (filename.c_str(), ios::out);
    
    if(!FILE) error("could not open Phenotypes output file!!\n");
    
    // write the heading line (nbPatch must be the max patch index!!!)
    FILE << (*(--_popPtr->get_vSamplePatch().end()))->get_ID()+1 << " " << _nb_trait << "\n";
    
    // write the names of the traits
    unsigned int t, i;
    for(t=0; t<_nb_trait; ++t){
        FILE << get_name() << "_trait-" << (t+1) << "\n";
    }
    
    //FILE.width(12);
    FILE.setf(ios::left,ios::adjustfield);
    FILE.precision(4);
    
    unsigned int a, mask;
    vector<TPatch*>::iterator curPop = _popPtr->get_vSamplePatch().begin();
    vector<TPatch*>::iterator endPop = _popPtr->get_vSamplePatch().end();
    for(i=0; curPop!=endPop; ++curPop) {
        for(a = 0, mask=1; a < NB_AGE_CLASSES; ++a, mask <<= 1) {
            if(mask & _age){
                if(_sex != 2) FHwrite(static_cast<age_idx>(a), FEM, FILE, *curPop, i);
                if(_sex != 1) FHwrite(static_cast<age_idx>(a), MAL, FILE, *curPop, i);
            }
        }
    }
}

// ----------------------------------------------------------------------------------------
// FHwrite
// ----------------------------------------------------------------------------------------
void TTQuantiFHvalue::FHwrite (const age_idx& cur_age, const sex_t& cur_sex, ostream& FILE,
                               TPatch* current_patch, const int& patch_id)
{
    unsigned int nbInd = current_patch->size(cur_sex, cur_age);
    if(!nbInd) return;
    
    TIndividual* ind;
    double value;
    unsigned int j, t;
    
    for(j = 0; j < nbInd; ++j){
        ind = current_patch->get(cur_sex,cur_age, j);
        FILE << (patch_id+1);
        for(t=0; t<_nb_trait; ++t){
            value = (this->*get_value_func_ptr)(ind, t);
            FILE << "\t" << (value==my_NAN ? my_NAN : value);
        }
        if(get_format()==2){
            FILE << "\t";
            write_individual_info_to_stream(FILE, ind, cur_age, cur_sex, '\t');
        }
        FILE << "\n";
    }
}

// ----------------------------------------------------------------------------------------
// get_genotype
// ----------------------------------------------------------------------------------------
double TTQuantiFHvalue::get_genotype(TIndividual* ind, const int& t){
    return ind->getTraitGenotype(_TTidx[t]);
}

// ----------------------------------------------------------------------------------------
// get_phenotype
// ----------------------------------------------------------------------------------------
double TTQuantiFHvalue::get_phenotype(TIndividual* ind, const int& t){
    return ind->getTraitPhenotype(_TTidx[t]);
}

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                             ******** StatHandler ********

// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
bool TTQuantiSH::init()
{
    StatHandler<TTQuantiSH>::init();
    return true;
}

// ----------------------------------------------------------------------------------------
// setStatsRecorders
// ----------------------------------------------------------------------------------------
bool TTQuantiSH::setStatRecorders(const string& t)   // ex. t: q1.adlt.fst_pair
{
    if(t=="quanti"){
        add("Genetic variance between patches","q.VgB",FLAT,ADULTS,0,0,0,0,&TTQuantiSH::getVgB);
        add("Genetic variance within patches","q.VgW",FLAT,ADULTS,0,0,0,0,&TTQuantiSH::getVgW);
        add("Phenotypic variance between patches","q.VpB",FLAT,ADULTS,0,&TTQuantiSH::getVpB);
        add("Phenotypic variance within patches","q.VpW",FLAT,ADULTS,0,&TTQuantiSH::getVpW);
        return true;
    }
    
    string type = _SHLinkedTrait->get_type();
    if(t[0] != type[0]) return false;         // wrong trait type   "q"?
    
    // is their a trait specification (ex. q1)?
    string::size_type pos=1;                   // is the second position a '.'
    if(t[pos]!= '.'){                          // trait index is specified (and this has to be an integer number)
        pos = t.find('.');
        unsigned int traitID;
        try{
            traitID = strTo<unsigned int>(t.substr(1, pos-1));
        }catch(...) {error("Stat 't' cannot be read: trait specification is not a number!\n", t.c_str());}
        if(!((!_trait_index && traitID==1)               // a single trait is present
             || (traitID==_trait_index)))    return false; // wrong index when several traits are present
    }
    string token = t.substr(pos+1);                   // remove the "q.": eg. "adlt.fst_pair"
    string i = "q";
    
    // get the end
    pos = token.find('_');
    string end;                                        // "end" is empty
    if(pos != string::npos){                           // there is an end
        end = token.substr(pos+1);                       // eg. "pair"
        token.erase(pos+1);                              // eg. "adlt.fst_"
    }
    
    // group quanti
    if(token=="varA_")  return add_perPatch(end,"Additive genetic variance","q.varA",FLAT,ADULTS,0,0,0,0,0,0,&TTQuantiSH::getVarA);
    if(token=="meanG_") return add_perPatch(end,"Genetic mean","q.meanG",FLAT,ADULTS,0,0,0,0,0,0,&TTQuantiSH::getMeanG);
    if(token=="varG_")  return add_perPatch(end,"Genetic variance","q.varG",FLAT,ADULTS,0,0,0,0,0,0,&TTQuantiSH::getVarG);
    if(token=="meanP_") return add_perPatch(end,"Phenotypic mean","q.meanP",FLAT,ADULTS,0,0,0,&TTQuantiSH::getMeanP);
    if(token=="varP_")  return add_perPatch(end,"Phenotypic variance","q.varP",FLAT,ADULTS,0,0,0,&TTQuantiSH::getVarP);
    if(token=="meanGfem_") return add_perPatch(end,"Genetic mean of females","q.meanGfem",FLAT,ADULTS,0,0,0,0,0,0,&TTQuantiSH::getMeanGfem);
    if(token=="varGfem_")  return add_perPatch(end,"Genetic variance of females","q.varGfem",FLAT,ADULTS,0,0,0,0,0,0,&TTQuantiSH::getVarGfem);
    if(token=="meanPfem_") return add_perPatch(end,"Phenotypic mean of females","q.meanPfem",FLAT,ADULTS,0,0,0,&TTQuantiSH::getMeanPfem);
    if(token=="varPfem_")  return add_perPatch(end,"Phenotypic variance of females","q.varPfem",FLAT,ADULTS,0,0,0,&TTQuantiSH::getVarPfem);
    if(token=="meanGmal_") return add_perPatch(end,"Genetic mean of males","q.meanGmal",FLAT,ADULTS,0,0,0,0,0,0,&TTQuantiSH::getMeanGmal);
    if(token=="varGmal_")  return add_perPatch(end,"Genetic variance of males","q.varGmal",FLAT,ADULTS,0,0,0,0,0,0,&TTQuantiSH::getVarGmal);
    if(token=="meanPmal_") return add_perPatch(end,"Phenotypic mean of males","q.meanPmal",FLAT,ADULTS,0,0,0,&TTQuantiSH::getMeanPmal);
    if(token=="varPmal_")  return add_perPatch(end,"Phenotypic variance of males","q.varPmal",FLAT,ADULTS,0,0,0,&TTQuantiSH::getVarPmal);
    if(token=="meanW_")    return add_perPatch(end,"Fitness mean","q.meanW",FLAT,ADULTS,0,0,0,&TTQuantiSH::getMeanW);
    if(token=="varW_")     return add_perPatch(end,"Fitness variance","q.varW",FLAT,ADULTS,0,0,0,&TTQuantiSH::getVarW);
    if(token=="meanWfem_") return add_perPatch(end,"Fitness mean of males","q.meanWfem",FLAT,ADULTS,0,0,0,&TTQuantiSH::getMeanWfem);
    if(token=="varWfem_")  return add_perPatch(end,"Fitness variance of males","q.varWfem",FLAT,ADULTS,0,0,0,&TTQuantiSH::getVarWfem);
    if(token=="meanWmal_") return add_perPatch(end,"Fitness mean of males","q.meanWmal",FLAT,ADULTS,0,0,0,&TTQuantiSH::getMeanWmal);
    if(token=="varWmal_")  return add_perPatch(end,"Fitness variance of males","q.varWmal",FLAT,ADULTS,0,0,0,&TTQuantiSH::getVarWmal);
    
    if(token=="VaW")     return add("Additive genetic variance within patches","q.VaW",FLAT,ADULTS,0,0,0,0,&TTQuantiSH::getVaW);
    if(token=="VgB")     return add("Genetic variance between patches","q.VgB",FLAT,ADULTS,0,0,0,0,&TTQuantiSH::getVgB);
    if(token=="VgW")     return add("Genetic variance within patches","q.VgW",FLAT,ADULTS,0,0,0,0,&TTQuantiSH::getVgW);
    if(token=="VpB")     return add("Phenotypic variance between patches","q.VpB",FLAT,ADULTS,0,&TTQuantiSH::getVpB);
    if(token=="VpW")     return add("Phenotypic variance within patches","q.VpW",FLAT,ADULTS,0,&TTQuantiSH::getVpW);
    if(token=="qst")     return add("Qst","q.qst",FLAT,ADULTS,0,0,0,0,&TTQuantiSH::getQst);
    if(token=="qst_")    return add_pairwisePatch(end,"Qst","q.qst",FLAT,ADULTS,0,0,0,0,0,0,&TTQuantiSH::getQst_ij);
    if(token=="qst.f")   return add("Qst corrected for inbreeding","q.qst.f",FLAT,ADULTS,0,0,0,0,&TTQuantiSH::getQstF);
    if(token=="qst.f_")  return add_pairwisePatch(end,"Qst corrected for inbreeding","q.qst.f",FLAT,ADULTS,0,0,0,0,0,0,&TTQuantiSH::getQstF_ij);
    
    
    // get the age
    pos = token.find('.');                        // find the second '.'
    if(pos == string::npos) return false;
    string ageToken = token.substr(0, pos);       // eg. "adlt"
    token.erase(0, pos+1);                        // eg. "fst_"
    age_t AGE;
    string ageStr;
    if(ageToken == "adlt")     {AGE = ADULTS;  ageStr = "adult";}
    else if(ageToken == "off") {AGE = OFFSPRG; ageStr = "offsprg";}
    else return false;
    
    if(set_stat_coancestry      (token, i, type, ageToken, end, AGE, ageStr)) return true; // adults and offspring
    if(set_stat_fstat           (token, i, type, ageToken, end, AGE, ageStr)) return true; // adults and offspring
    if(set_stat_all_freq_local  (token, i, type, ageToken, end, AGE, ageStr)) return true; // adults and offspring
    if(set_stat_all_freq_global (token, i, type, ageToken, end, AGE, ageStr)) return true; // adults and offspring
    if(set_stat_locus_freq_local(token, i, type, ageToken, end, AGE, ageStr)) return true; // adults and offspring
    if(set_stat_locus_freq_global(token, i, type, ageToken, end, AGE, ageStr))return true; // adults and offspring
    if(set_stat_LD              (token, i, type, ageToken, end, AGE, ageStr)) return true; // adults and offspring
    
    return false;
}

// ----------------------------------------------------------------------------------------
// setVarAdditiveGenetic
// ----------------------------------------------------------------------------------------
/** computation of the additive genetic variance (for any case) for each patch */
void
TTQuantiSH::setVar_Va(const age_idx& AGE)
{
    // check if the table has already been computed
    if(already_computed(_computed[20], AGE)) return;
    
    if(!_varA)   _varA   = new double[get_current_nbSamplePatch()];
    else if(get_current_nbSamplePatch()!=get_last_nbSamplePatch()){ // number of sampled patches may change over time
        delete[] _varA;
        _varA   = new double[get_current_nbSamplePatch()];
    }
    set_alleleFreq(AGE);        // compute allele frequencies
    setMeanAndVar_Vg(AGE);                  // compute genetic mean and variance
    double mean; // not really used...
    
    // for each population
    vector<TPatch*>::iterator curPop = get_vSamplePatch().begin();
    vector<TPatch*>::iterator endPop = get_vSamplePatch().end();
    for(unsigned int i=0; curPop!=endPop; ++curPop, ++i){ // for each patch
        (this->*get_Va_ofPatch_func_ptr)(*curPop, AGE, mean, _varA[i], _alleleFreq_local[AGE][i]);
    }
}

// ----------------------------------------------------------------------------------------
// setVarAdditiveGenetic
// ----------------------------------------------------------------------------------------
/** computation of the additive genetic variance (only for random matings)
 * much quicker than the version for non-random matings
 * meanA is not computed. But needed for function pointer
 */
void
TTQuantiSH::get_Va_ofPatch_random_mating(TPatch* curPop, const age_idx& AGE,
                                         double& meanA, double& varA, map<unsigned char, double>* freqs)
{
    vector<TIndividual*>& curFem = curPop->get_sampled_inds(FEM, AGE);
    vector<TIndividual*>& curMal = curPop->get_sampled_inds(MAL, AGE);
    unsigned int sizeF = (unsigned int)curFem.size(),
    sizeM = (unsigned int)curMal.size();
    unsigned int size  = sizeF + sizeM;
    if(!size) {varA = my_NAN; return;}
    
    double G, meanG=0;
    double *curElem;
    TTree<unsigned char, double>       condMeanGG(3, 0);
    TTree<unsigned char, unsigned int> condMeanGGsize(3, 0);
    unsigned char aTemp[3];
    unsigned int l, a;
    unsigned char a1, a2;
    double* aGeno = new double[size];
    unsigned char** genes;
    
    // females
    vector<TIndividual*>::iterator curInd, endInd;
    for(curInd= curFem.begin(), endInd=curFem.end(); curInd!=endInd; ++curInd) {
        meanG += G = (*curInd)->getTraitGenotype(_SHLinkedTraitIndex);          // genotype
        genes = (unsigned char**)(*curInd)->getTrait(_SHLinkedTraitIndex)->get_sequence();  // sequence
        for (l = 0; l < _nb_locus; ++l){
            aTemp[0] = l;
            a1 = genes[l][0];          // 1. allele
            a2 = genes[l][1];          // 2. allele
            
            if(a1<a2){
                aTemp[1]=a1;
                aTemp[2]=a2;
            }
            else{
                aTemp[1]=a2;
                aTemp[2]=a1;
            }
            condMeanGG.get(aTemp) += G;
            condMeanGGsize.get(aTemp) += 1;
        }
    } // for each female
    
    // males
    for(curInd= curMal.begin(), endInd=curMal.end(); curInd!=endInd; ++curInd) {
        meanG += G = (*curInd)->getTraitGenotype(_SHLinkedTraitIndex);          // genotype
        genes = (unsigned char**)(*curInd)->getTrait(_SHLinkedTraitIndex)->get_sequence();  // sequence
        for (l = 0; l < _nb_locus; ++l){
            aTemp[0] = l;
            a1 = genes[l][0];          // 1. allele
            a2 = genes[l][1];          // 2. allele
            
            if(a1<a2){
                aTemp[1]=a1;
                aTemp[2]=a2;
            }
            else{
                aTemp[1]=a2;
                aTemp[2]=a1;
            }
            condMeanGG.get(aTemp) += G;
            condMeanGGsize.get(aTemp) += 1;
        }
    } // for each male
    
    // for each present genotype
    meanG /= size;
    map<unsigned char, double>* alphaStar = new map<unsigned char, double>[_nb_locus];
    for(curElem = condMeanGG.first(aTemp); curElem; curElem = condMeanGG.next(aTemp)){
        // correct the dominance effects for the sample number
        (*curElem) /= condMeanGGsize.get(aTemp);
        l  = aTemp[0];
        a1 = aTemp[1];
        a2 = aTemp[2];
        
        // compute the average excess (alphaStar) which is in this case (random mating) identical to the additive effects
        // meanG will later be substracted from alphaStar!!!
        assert(freqs[l].find(a1) != freqs[l].end() && freqs[l].find(a2) != freqs[l].end());
        alphaStar[l][a1]   += (*curElem) * freqs[l][a2];
        if(a1 != a2){                 // if the two alles are not identical ...
            alphaStar[l][a2]   += (*curElem) * freqs[l][a1];
        }
    }
    
    // compute the additive variance
    varA = 0;
    double aStar; // corrected alphaStar value (subtraction of meanG)
    map<unsigned char, double>::iterator pos, end;
    for(l=0; l<_nb_locus; ++l){
        pos = freqs[l].begin();
        end = freqs[l].end();
        for(; pos != end; ++pos){    // for each allele
            a = (unsigned int)pos->first;               // get the allele
            aStar = alphaStar[l][a] - meanG;
            varA += aStar * aStar * pos->second;        // only valable when random mating  (pos->second == freqs[l][a])
        }
    }
    varA *= 2;    // we have diploid individuals ...
    
    // delete arrays
    delete[] alphaStar;
    delete[] aGeno;
}

// ----------------------------------------------------------------------------------------
// setVarAdditiveGenetic
// ----------------------------------------------------------------------------------------
/** computation of the additive genetic variance (for any case)  */
void
TTQuantiSH::get_Va_ofPatch_regression(TPatch* curPop, const age_idx& AGE,
                                      double& meanA, double& varA,
                                      map<unsigned char, double>* freqs)
{
    vector<TIndividual*>& curFem = curPop->get_sampled_inds(FEM, AGE);
    vector<TIndividual*>& curMal = curPop->get_sampled_inds(MAL, AGE);
    unsigned int sizeF = (unsigned int)curFem.size(),
    sizeM = (unsigned int)curMal.size();
    unsigned int size  = sizeF + sizeM;
    if(!size) {varA = my_NAN; return;}
    
    unsigned int i, l;
    double* arrayG = new double[size];
    unsigned char a1, a2;
    unsigned char** genes;
    double G, meanG=0;
    unsigned int nbFailure = 0;
    
    // dynamic pairwise allele container: condMeanGG[l][a1][a2]
    map<unsigned char, map<unsigned char, double> >* condMeanGG = new map<unsigned char, map<unsigned char, double> >[_nb_locus];
    map<unsigned char, double>::iterator pos2, end2;                         // inner iterator
    map<unsigned char, map<unsigned char, double> >::iterator pos1, end1;    // outer iterator
    
    // container: condMeanGGsize[l][a1][a2]
    map<unsigned char, map<unsigned char, unsigned int> >* condMeanGGsize = new map<unsigned char, map<unsigned char, unsigned int> >[_nb_locus];
    map<unsigned char, unsigned int >::iterator pos2_size;                      // inner iterator
    map<unsigned char, map<unsigned char, unsigned int> >::iterator pos1_size;  // outer iterator
    
    
    // make matrix of predictor values (each column for an allele, each row for an individual)
    map<unsigned char, int*>*  predMatrix = new map<unsigned char, int*>[_nb_locus]; // array[locus][allele][individual]
    map<unsigned char, int*>::iterator posM;
    
    // females
    vector<TIndividual*>::iterator curInd, endInd;
    for(i=0, curInd= curFem.begin(), endInd=curFem.end(); curInd!=endInd; ++curInd, ++i) {
        // get the genotype and genotypic value
        meanG += arrayG[i] = G = (*curInd)->getTraitGenotype(_SHLinkedTraitIndex);          // genotypic value
        genes = (unsigned char**)(*curInd)->getTrait(_SHLinkedTraitIndex)->get_sequence();  // genotype
        for (l = 0; l < _nb_locus; ++l){
            a1 = genes[l][0];          // 1. allele
            a2 = genes[l][1];          // 2. allele
            
            if(a1<a2){
                condMeanGG[l][a1][a2] += G;
                ++condMeanGGsize[l][a1][a2];
            }
            else{
                condMeanGG[l][a2][a1] += G;
                ++condMeanGGsize[l][a2][a1];
            }
            
            // create the predictor matrix
            posM = predMatrix[l].find(a1); 			// allele 1 already present?
            if(posM == predMatrix[l].end()){    // if a1 apears the first time
                posM = predMatrix[l].insert(pair<unsigned char, int*>(a1, ARRAY::new_1D<int>(size, (int)0))).first;
            }
            ++posM->second[i];
            posM = predMatrix[l].find(a2); 			// allele 2 already present?
            if(posM == predMatrix[l].end()){    // if a2 apears the first time
                posM = predMatrix[l].insert(pair<unsigned char, int*>(a2, ARRAY::new_1D<int>(size, (int)0))).first;
            }
            ++posM->second[i];
        }
    } // for each indvidual
    
    // males
    for(curInd= curMal.begin(), endInd=curMal.end(); curInd!=endInd; ++curInd, ++i) {
        // get the genotype and genotypic value
        meanG += arrayG[i] = G = (*curInd)->getTraitGenotype(_SHLinkedTraitIndex);          // genotypic value
        genes = (unsigned char**)(*curInd)->getTrait(_SHLinkedTraitIndex)->get_sequence();  // genotype
        for (l = 0; l < _nb_locus; ++l){
            a1 = genes[l][0];          // 1. allele
            a2 = genes[l][1];          // 2. allele
            
            if(a1<a2){
                condMeanGG[l][a1][a2] += G;
                ++condMeanGGsize[l][a1][a2];
            }
            else{
                condMeanGG[l][a2][a1] += G;
                ++condMeanGGsize[l][a2][a1];
            }
            
            // create the predictor matrix
            posM = predMatrix[l].find(a1); 			// allele 1 already present?
            if(posM == predMatrix[l].end()){    // if a1 apears the first time
                posM = predMatrix[l].insert(pair<unsigned char, int*>(a1, ARRAY::new_1D<int>(size, (int)0))).first;
            }
            ++posM->second[i];
            posM = predMatrix[l].find(a2); 			// allele 2 already present?
            if(posM == predMatrix[l].end()){    // if a2 apears the first time
                posM = predMatrix[l].insert(pair<unsigned char, int*>(a2, ARRAY::new_1D<int>(size, (int)0))).first;
            }
            ++posM->second[i];
        }
    } // for each indvidual
    
    // correct the genotype by the mean genotype
    meanG /= size;
    for(i = 0; i < size; ++i) {
        arrayG[i] -= meanG;
    }
    
    // compute the additive effects (alpha, time consuming!)
    map<unsigned char, double>* alpha = new map<unsigned char, double>[_nb_locus];
    for(l=0; l<_nb_locus; ++l){
        if(!compute_alpha(arrayG, predMatrix[l], size, alpha[l], freqs[l])){
            if(!remove_private_alleles_compute_alpha(curPop, sizeF, sizeM, alpha[l], arrayG, AGE, freqs[l], l)){
                ++nbFailure;
            }
        }
    }
    
    // compute the average excess (alphaStar)  (loop through all present pairwise-combinations of alleles)
    map<unsigned char, double>* alphaStar = new map<unsigned char, double>[_nb_locus];
    for(l=0; l<_nb_locus; ++l){                             // for each locus
        pos1 = condMeanGG[l].begin();
        end1 = condMeanGG[l].end();
        pos1_size = condMeanGGsize[l].begin();
        for(; pos1 != end1; ++pos1, ++pos1_size){           	// for each 1. allele
            a1   = pos1->first;                                 // get the 1. allele
            pos2 = pos1->second.begin();
            end2 = pos1->second.end();
            pos2_size = pos1_size->second.begin();
            for(; pos2 != end2; ++pos2, ++pos2_size){           // for each 2. allele
                a2 = pos2->first;                                 // get the 2. allele
                
                pos2->second /= pos2_size->second;                // correct the dominance effects for the sample number
                
                // compute the average excess (alphaStar) (meanG will be later substracted!!!)
                assert(freqs[l].find(a2) != freqs[l].end());
                alphaStar[l][a1]   += pos2->second * freqs[l][a2];
                if(a1 != a2){                                     // if the two alles are not identical ...
                    assert(freqs[l].find(a1) != freqs[l].end());
                    alphaStar[l][a2]   += pos2->second * freqs[l][a1];
                }
            }
        }
    }
    
    // compute the additive variance  var = 2*sum(p a a*)
    varA = 0;
    map<unsigned char, double>::iterator curAlpha, curAlphaStar, curFreq, endFreq;
    for(l=0; l<_nb_locus; ++l){
        curAlpha     = alpha[l].begin();
        curAlphaStar = alphaStar[l].begin();
        curFreq      = freqs[l].begin();
        endFreq      = freqs[l].end();
        for(; curFreq != endFreq; ++curAlpha, ++curAlphaStar, ++curFreq){   // for each allele
            varA += curAlpha->second * (curAlphaStar->second - meanG) * curFreq->second; // any case
            // var += alphaStar[l][a]*alphaStar[l][a]*freqs[l][a1];            // limited to random mating
        }
    }
    varA *= 2;    // we have diploid individuals ...
    
    // if more than 10% of the loci were not able to be computed set Va to NAN
    if(nbFailure > 0.1*_nb_locus) varA = my_NAN;
    
    // clean up
    delete[] arrayG;
    delete[] alpha;
    delete[] alphaStar;
    delete[] condMeanGG;
    delete[] condMeanGGsize;
    
    map<unsigned char, int*>::iterator pos, end;
    for(l=0; l<_nb_locus; ++l){
        pos = predMatrix[l].begin();
        end = predMatrix[l].end();
        for(; pos != end; ++pos){
            delete[] pos->second;
        }
    }
    delete[] predMatrix;
}

// ----------------------------------------------------------------------------------------
// compute_alpha
// ----------------------------------------------------------------------------------------
/** traditional sum of squares and products method of calculation
 * compute the alpha
 * returns false if the regression was not able to compute due to
 * individuals which have two private alleles
 * singular value decomposition method
 */
bool
TTQuantiSH::compute_alpha(double* y, const map<unsigned char, int*>& x, const unsigned int& nb_ind,
                          map<unsigned char, double>& alpha, const map<unsigned char, double>& availableAllele)
{
    alpha.clear();
    map<unsigned char, double>::const_iterator pos, end = availableAllele.end();
    map<unsigned char, int*>::const_iterator posX, endX = x.end();
    unsigned int i, nbAllele = (unsigned int)availableAllele.size();  // number of alleles present in the population
    
    // test if the population is fixed for one allele
    if(nbAllele <=1){
        for(pos = availableAllele.begin(); pos != end; ++pos){
            alpha.insert(pair<unsigned char, double>(pos->first, 0));  // if singular matrix
        }
        return (nbAllele != 0); // was the computation successful or not?
    }
    
    // test if the population is fixed for two alleles
    if(nbAllele ==2){
        posX = x.begin();                   // the first allele
        for(i = 0; i<nb_ind; ++i){
            if(posX->second[i] != 1) break;   // if not all individuals have the same two alleles
        }
        if(i==nb_ind){      // if singular matrix
            for(pos = availableAllele.begin(); pos != end; ++pos){
                alpha.insert(pair<unsigned char, double>(pos->first, 0));  // if singular matrix
            }
            return true;
        }
    }
    
    // create the matrixes
    Matrix X(nb_ind, nbAllele);          // make matrix of predictor values
    ColumnVector Y(nb_ind);
    try{
        // load the predictive matrix
        for(i = 1, posX = x.begin(); posX != endX; ++posX, ++i){
            X.Column(i) << posX->second;       // the matrix index starts with 1!!!
        }
        
        // load Y values
        Y << y;
        
        // calculate estimate
        ColumnVector A = (X.t() * X).i() * (X.t() * Y);  // .i(): inverse; .t(): transpose
        
        for(i = 1, pos = availableAllele.begin(); pos != end; ++pos, ++i){
            alpha.insert(pair<unsigned char, double>(pos->first, A(i)));  // all the other values
        }
    }catch(BaseException&){
        // the problem is that "X.t() * X" can result in a singular matrix
        // thus there is no solution to the regression
        for(pos = availableAllele.begin(); pos != end; ++pos){
            alpha.insert(pair<unsigned char, double>(pos->first, 0));  // if singular matrix
        }                                                            // no additive effect at this locus
        
        // the exception has to be deleted
        delete[] BaseException::what_error;
        BaseException::what_error = NULL;
        
        return false;
    }
    return true;
}

// ----------------------------------------------------------------------------------------
// remove_private_alleles_compute_alpha
// ----------------------------------------------------------------------------------------
/** this function is called when the regression cannot be resolved
 * (i.e. when t(X)*X results in a singular matrix det(t(X)*X) = 0)
 * this function removes all individuals having two private alleles, i.e.
 * alleles which are unique to this individual.
 * Then the additive effects (alpha) are recomputed.
 * If more than 10 percent of the individuals have to be removed or the matrix is
 * again singular the computation is assumed to be unsuccessful
 * returns true if it was able to compute Va otherwise false
 */
bool
TTQuantiSH::remove_private_alleles_compute_alpha(TPatch* crnt_patch, const unsigned int& sizeF,
                                                 const unsigned int& sizeM, map<unsigned char, double>& alpha,
                                                 double* arrayG, const age_idx& age_pos,
                                                 map<unsigned char, double>& allele_freq, const unsigned int& l)
{
    // the regression was not possible to compute: remove all individuals with two private alleles
    unsigned int i;
    unsigned char* g;
    vector<unsigned char*> geno;    // vector with all "good" locus genotypes
    vector<double>   vectorG1;      // vector of the corrected genotypic values
    vector<int> rem_all;            // vector of alleles to remove
    bool stop = false;
    double private_allele_freq = 1.0/2.0/(sizeM+sizeF); 	// allele frequency of a private allele
    
    // check each individual if it has two private alleles
    for(i = 0; i < sizeF+sizeM; ++i) {
        if(i<sizeF) g = (unsigned char*)crnt_patch->get(FEM, age_pos, i)->getTrait(_SHLinkedTraitIndex)->get_sequence()[l];       // get the female
        else        g = (unsigned char*)crnt_patch->get(MAL, age_pos, i-sizeF)->getTrait(_SHLinkedTraitIndex)->get_sequence()[l]; // get the male
        
        assert(allele_freq.find(g[0]) != allele_freq.end() && allele_freq.find(g[1]) != allele_freq.end());
        if(   allele_freq[g[0]] == private_allele_freq   // check if both alleles are private
           && allele_freq[g[1]] == private_allele_freq){
            rem_all.push_back(g[0]);                        // remove these alleles
            rem_all.push_back(g[1]);
        }
        else {
            geno.push_back(g);                              // add the locus genotype
            vectorG1.push_back(arrayG[i]);
        }
    }
    
    // if more then 10% of the individuals have to be removed consider it as not computable
    if(geno.size() > 0.1*(sizeF+sizeM)) stop = true;
    else{
        // create the new available allele vector
        map<unsigned char, double> availableAllele1;
        map<unsigned char, double>::iterator pos, end;
        pos = allele_freq.begin();
        end = allele_freq.end();
        for(; pos != end; ++pos){
            if(find(rem_all.begin(), rem_all.end(), pos->first) == rem_all.end()) availableAllele1[pos->first] = pos->second;
        }
        
        // create the new predMatrix
        unsigned int size1 = (unsigned int)geno.size();
        double* arrayG1 = new double[size1];
        map<unsigned char, int*> predMatrix1;       // array[allele][individual]
        map<unsigned char, int*>::iterator posM;
        unsigned char a1, a2;
        
        for(i=0; i<size1; ++i){                     // for each individual
            a1 = geno[i][0];                        // allele 1
            posM = predMatrix1.find(a1);            // allele already present?
            if(posM == predMatrix1.end()) ARRAY::create_1D(predMatrix1[a1], size1, (int)0);
            ++posM->second[i];
            a2 = geno[i][1];                        // allele 2
            posM = predMatrix1.find(a2);            // allele already present?
            if(posM == predMatrix1.end()) ARRAY::create_1D(predMatrix1[a2], size1, (int)0);
            ++posM->second[i];
            arrayG1[i] = vectorG1[i];
        }
        
        // compute the alpha again
        if(!compute_alpha(arrayG1, predMatrix1, size1, alpha, availableAllele1)) stop = true;
    }
    
    if(stop){// the least-square regression had again no solution
        // generate a warning the first 10 times the problem occurs and then every 100 time
        static unsigned int passes;      // is initialized with 0
        static unsigned int replicate;   // needed to reset the counter when a new replicate is computed
        if(replicate != get_current_replicate()){
            passes=0;
            replicate = get_current_replicate();
        }
        ++passes;
        if(passes < 10 || !(passes%100)){
            warning("Va could not be correctly estimated (%i. time at generation %i, see manual parameter 'quanti_va_model')!\n", passes, get_current_generation());
        }
        
        return false; // unsuccessful computation of Va
    }
    return true;	  // successful computation of Va
}

// ----------------------------------------------------------------------------------------
// setMeanAndVarGenotyp
// ----------------------------------------------------------------------------------------
void
TTQuantiSH::setMeanAndVar_Vg(const age_idx& AGE)
{
    // check if the table has already been computed
    if(already_computed(_computed[21], AGE)) return;
    
    if(!_meanG) ARRAY::create_1D<double*>(_meanG, 3, NULL);
    if(!_meanG[2]) _meanG[2] = new double[get_current_nbSamplePatch()];
    else if(get_current_nbSamplePatch()!=get_last_nbSamplePatch()){ // number of sampled patches may change over time
        delete[] _meanG[2];
        _meanG[2] = new double[get_current_nbSamplePatch()];
    }
    
    if(!_varG) ARRAY::create_1D<double*>(_varG, 3, NULL);
    if(!_varG[2]) _varG[2] = new double[get_current_nbSamplePatch()];
    else if(get_current_nbSamplePatch()!=get_last_nbSamplePatch()){ // number of sampled patches may change over time
        delete[] _varG[2];
        _varG[2] = new double[get_current_nbSamplePatch()];
    }
    
    
    // for each population
    vector<TPatch*>::iterator curPop = get_vSamplePatch().begin();
    vector<TPatch*>::iterator endPop = get_vSamplePatch().end();
    for(unsigned int i=0; curPop!=endPop; ++curPop, ++i){
        setMeanAndVar_Vg_ofPatch(*curPop, AGE, _meanG[2][i], _varG[2][i]);
    }
}

// ----------------------------------------------------------------------------------------
// setMeanAndVarGenotyp
// ----------------------------------------------------------------------------------------
void
TTQuantiSH::setMeanAndVar_Vg(const age_idx& AGE, sex_t SEX)
{
    // check if the table has already been computed
    //if(already_computed(_computed[21], AGE)) return; // problem since there are two sexes
    
    if(!_meanG) ARRAY::create_1D<double*>(_meanG, 3, NULL);
    if(!_meanG[SEX]) _meanG[SEX] = new double[get_current_nbSamplePatch()];
    else if(get_current_nbSamplePatch()!=get_last_nbSamplePatch()){ // number of sampled patches may change over time
        delete[] _meanG[SEX];
        _meanG[SEX] = new double[get_current_nbSamplePatch()];
    }
    
    if(!_varG) ARRAY::create_1D<double*>(_varG, 3, NULL);
    if(!_varG[SEX]) _varG[SEX] = new double[get_current_nbSamplePatch()];
    else if(get_current_nbSamplePatch()!=get_last_nbSamplePatch()){ // number of sampled patches may change over time
        delete[] _varG[SEX];
        _varG[SEX] = new double[get_current_nbSamplePatch()];
    }
    
    
    // for each population
    vector<TPatch*>::iterator curPop = get_vSamplePatch().begin();
    vector<TPatch*>::iterator endPop = get_vSamplePatch().end();
    for(unsigned int i=0; curPop!=endPop; ++curPop, ++i){
        setMeanAndVar_Vg_ofPatch(*curPop, AGE, _meanG[SEX][i], _varG[SEX][i], SEX);
    }
}

// ----------------------------------------------------------------------------------------
// setMeanAndVarGenotyp
// ----------------------------------------------------------------------------------------
/** compute Vg across all individuals of the patch */
void
TTQuantiSH::setMeanAndVar_Vg_ofPatch_allInds(TPatch* curPop, const age_idx& AGE,
                                             double& meanG, double& varG, map<unsigned char, double>* freqs)
{
    // create a temporary array with all individuals
    vector<TIndividual*>& curFem = curPop->get_containers(FEM, AGE);
    vector<TIndividual*>& curMal = curPop->get_containers(MAL, AGE);
    unsigned int sizeF = (unsigned int)curFem.size(),
    sizeM = (unsigned int)curMal.size();
    unsigned int size  = sizeF + sizeM;
    if(!size){
        meanG = varG = my_NAN;
        return;
    }
    
    double* array = new double[size];
    vector<TIndividual*>::iterator curInd, endInd;
    unsigned int f=0;
    for(curInd=curFem.begin(), endInd=curFem.end(); curInd!=endInd; ++curInd, ++f) {
        array[f] = (*curInd)->getTraitGenotype(_SHLinkedTraitIndex);
    }
    for(curInd=curMal.begin(), endInd=curMal.end(); curInd!=endInd; ++curInd, ++f) {
        array[f] = (*curInd)->getTraitGenotype(_SHLinkedTraitIndex);
    }
    
    // compute mean and var
    meanG = ARRAY::mean(array, size);
    varG  = ARRAY::var(array, size, meanG);
    delete[] array;
}

// ----------------------------------------------------------------------------------------
// setMeanAndVarGenotyp
// ----------------------------------------------------------------------------------------
/** compute Vg for the sampled individuals */
void
TTQuantiSH::setMeanAndVar_Vg_ofPatch(TPatch* curPop, const age_idx& AGE,
                                     double& meanG, double& varG, map<unsigned char, double>* freqs)
{
    // create a temporary array with all individuals
    vector<TIndividual*>& curFem = curPop->get_sampled_inds(FEM, AGE);
    vector<TIndividual*>& curMal = curPop->get_sampled_inds(MAL, AGE);
    unsigned int sizeF = (unsigned int)curFem.size(),
    sizeM = (unsigned int)curMal.size();
    unsigned int size  = sizeF + sizeM;
    if(!size){
        meanG = varG = my_NAN;
        return;
    }
    
    double* array = new double[size];
    vector<TIndividual*>::iterator curInd, endInd;
    unsigned int f=0;
    for(curInd=curFem.begin(), endInd=curFem.end(); curInd!=endInd; ++curInd, ++f) {
        array[f] = (*curInd)->getTraitGenotype(_SHLinkedTraitIndex);
    }
    for(curInd=curMal.begin(), endInd=curMal.end(); curInd!=endInd; ++curInd, ++f) {
        array[f] = (*curInd)->getTraitGenotype(_SHLinkedTraitIndex);
    }
    
    // compute mean and var
    meanG = ARRAY::mean(array, size);
    varG  = ARRAY::var(array, size, meanG);
    delete[] array;
}

// ----------------------------------------------------------------------------------------
// setMeanAndVarGenotyp
// ----------------------------------------------------------------------------------------
void
TTQuantiSH::setMeanAndVar_Vg_ofPatch(TPatch* curPop, const age_idx& AGE,
                                     double& meanG, double& varG, sex_t SEX, map<unsigned char, double>* freqs)
{
    // create a temporary array with all individuals
    vector<TIndividual*>& cur = curPop->get_sampled_inds(SEX, AGE);
    unsigned int size = (unsigned int)cur.size();
    if(!size){
        meanG = varG = my_NAN;
        return;
    }
    
    double* array = new double[size];
    vector<TIndividual*>::iterator curInd, endInd;
    unsigned int f=0;
    for(curInd=cur.begin(), endInd=cur.end(); curInd!=endInd; ++curInd, ++f) {
        array[f] = (*curInd)->getTraitGenotype(_SHLinkedTraitIndex);
    }
    
    // compute mean and var
    meanG = ARRAY::mean(array, size);
    varG  = ARRAY::var(array, size, meanG);
    delete[] array;
}

// ----------------------------------------------------------------------------------------
// setMeanAndVarPhenotype
// ----------------------------------------------------------------------------------------
void
TTQuantiSH::setMeanAndVar_Vp()
{
    // check if the table has already been computed
    if(already_computed(_computed[22])) return;
    
    if(!_meanP) ARRAY::create_1D<double*>(_meanP, 3, NULL);
    if(!_meanP[2]) _meanP[2] = new double[get_current_nbSamplePatch()];
    else if(get_current_nbSamplePatch()!=get_last_nbSamplePatch()){ // number of sampled patches may change over time
        delete[] _meanP[2];
        _meanP[2] = new double[get_current_nbSamplePatch()];
    }
    
    if(!_varP) ARRAY::create_1D<double*>(_varP, 3, NULL);
    if(!_varP[2]) _varP[2] = new double[get_current_nbSamplePatch()];
    else if(get_current_nbSamplePatch()!=get_last_nbSamplePatch()){ // number of sampled patches may change over time
        delete[] _varP[2];
        _varP[2] = new double[get_current_nbSamplePatch()];
    }
    
    // for each population
    vector<TPatch*>::iterator curPop = get_vSamplePatch().begin();
    vector<TPatch*>::iterator endPop = get_vSamplePatch().end();
    for(unsigned int i=0; curPop!=endPop; ++curPop, ++i){ // for each patch
        setMeanAndVar_Vp_ofPatch(*curPop, _meanP[2][i], _varP[2][i]);
    }
}

// ----------------------------------------------------------------------------------------
// setMeanAndVarPhenotype
// ----------------------------------------------------------------------------------------
void
TTQuantiSH::setMeanAndVar_Vp(sex_t SEX)
{
    // check if the table has already been computed
    //if(already_computed(_computed[22])) return; problem since there are two sexes
    
    if(!_meanP) ARRAY::create_1D<double*>(_meanP, 3, NULL);
    if(!_meanP[SEX]) _meanP[SEX] = new double[get_current_nbSamplePatch()];
    else if(get_current_nbSamplePatch()!=get_last_nbSamplePatch()){ // number of sampled patches may change over time
        delete[] _meanP[SEX];
        _meanP[SEX] = new double[get_current_nbSamplePatch()];
    }
    
    if(!_varP) ARRAY::create_1D<double*>(_varP, 3, NULL);
    if(!_varP[SEX]) _varP[SEX] = new double[get_current_nbSamplePatch()];
    else if(get_current_nbSamplePatch()!=get_last_nbSamplePatch()){ // number of sampled patches may change over time
        delete[] _varP[SEX];
        _varP[SEX] = new double[get_current_nbSamplePatch()];
    }
    
    // for each population
    vector<TPatch*>::iterator curPop = get_vSamplePatch().begin();
    vector<TPatch*>::iterator endPop = get_vSamplePatch().end();
    for(unsigned int i=0; curPop!=endPop; ++curPop, ++i){ // for each patch
        setMeanAndVar_Vp_ofPatch(*curPop, _meanP[SEX][i], _varP[SEX][i], SEX);
    }
}

// ----------------------------------------------------------------------------------------
// setMeanAndVarPhenotype
// ----------------------------------------------------------------------------------------
void
TTQuantiSH::setMeanAndVar_Vp_ofPatch(TPatch* curPop, double& meanP, double& varP)
{
    vector<TIndividual*>& curFem = curPop->get_sampled_inds(FEM, ADLTx);
    vector<TIndividual*>& curMal = curPop->get_sampled_inds(MAL, ADLTx);
    unsigned int sizeF = (unsigned int)curFem.size(),
    sizeM = (unsigned int)curMal.size();
    unsigned int size  = sizeF + sizeM;
    
    // if the patch is empty or the phenotype is not yet computed -> stop
    if(!((sizeF && curFem.front()->getTraitPhenotype(_SHLinkedTraitIndex) != my_NAN)
         || (sizeM && curMal.front()->getTraitPhenotype(_SHLinkedTraitIndex) != my_NAN))){
        meanP = varP = my_NAN;
        return;
    }
    
    double* array = new double[size];
    
    vector<TIndividual*>::iterator curInd, endInd;
    unsigned int f=0;
    for(curInd=curFem.begin(), endInd=curFem.end(); curInd!=endInd; ++curInd, ++f) {
        array[f] = (*curInd)->getTraitPhenotype(_SHLinkedTraitIndex);
    }
    for(curInd=curMal.begin(), endInd=curMal.end(); curInd!=endInd; ++curInd, ++f) {
        array[f] = (*curInd)->getTraitPhenotype(_SHLinkedTraitIndex);
    }
    
    // compute mean and var
    meanP = ARRAY::mean(array, size);
    varP  = ARRAY::var(array, size, meanP);
    delete[] array;
}

// ----------------------------------------------------------------------------------------
// setMeanAndVarPhenotype
// ----------------------------------------------------------------------------------------
void
TTQuantiSH::setMeanAndVar_Vp_ofPatch(TPatch* curPop, double& meanP, double& varP, sex_t SEX)
{
    vector<TIndividual*>& cur = curPop->get_sampled_inds(SEX, ADLTx);
    unsigned int size = (unsigned int)cur.size();
    
    // if the patch is empty or the phenotype is not yet computed -> stop
    if(!(size && cur.front()->getTraitPhenotype(_SHLinkedTraitIndex) != my_NAN)){
        meanP = varP = my_NAN;
        return;
    }
    
    double* array = new double[size];
    
    vector<TIndividual*>::iterator curInd, endInd;
    unsigned int f=0;
    for(curInd=cur.begin(), endInd=cur.end(); curInd!=endInd; ++curInd, ++f) {
        array[f] = (*curInd)->getTraitPhenotype(_SHLinkedTraitIndex);
    }
    
    // compute mean and var
    meanP = ARRAY::mean(array, size);
    varP  = ARRAY::var(array, size, meanP);
    delete[] array;
}

// ----------------------------------------------------------------------------------------
// setMeanAndVarFitness
// ----------------------------------------------------------------------------------------
void
TTQuantiSH::setMeanAndVar_Wp()
{
    // check if the table has already been computed
    if(already_computed(_computed[22])) return;
    
    if(!_meanW) ARRAY::create_1D<double*>(_meanW, 3, NULL);
    if(!_meanW[2]) _meanW[2] = new double[get_current_nbSamplePatch()];
    else if(get_current_nbSamplePatch()!=get_last_nbSamplePatch()){ // number of sampled patches may change over time
        delete[] _meanW[2];
        _meanW[2] = new double[get_current_nbSamplePatch()];
    }
    
    if(!_varW) ARRAY::create_1D<double*>(_varW, 3, NULL);
    if(!_varW[2]) _varW[2] = new double[get_current_nbSamplePatch()];
    else if(get_current_nbSamplePatch()!=get_last_nbSamplePatch()){ // number of sampled patches may change over time
        delete[] _varW[2];
        _varW[2] = new double[get_current_nbSamplePatch()];
    }
    
    // for each population
    vector<TPatch*>::iterator curPop = get_vSamplePatch().begin();
    vector<TPatch*>::iterator endPop = get_vSamplePatch().end();
    for(unsigned int i=0; curPop!=endPop; ++curPop, ++i){ // for each patch
        setMeanAndVar_Wp_ofPatch(*curPop, _meanW[2][i], _varW[2][i]);
    }
}

// ----------------------------------------------------------------------------------------
// setMeanAndVarFitness
// ----------------------------------------------------------------------------------------
void
TTQuantiSH::setMeanAndVar_Wp(sex_t SEX)
{
    // check if the table has already been computed
    //if(already_computed(_computed[22])) return; problem since there are two sexes
    
    if(!_meanW) ARRAY::create_1D<double*>(_meanW, 3, NULL);
    if(!_meanW[SEX]) _meanW[SEX] = new double[get_current_nbSamplePatch()];
    else if(get_current_nbSamplePatch()!=get_last_nbSamplePatch()){ // number of sampled patches may change over time
        delete[] _meanW[SEX];
        _meanW[SEX] = new double[get_current_nbSamplePatch()];
    }
    
    if(!_varW) ARRAY::create_1D<double*>(_varW, 3, NULL);
    if(!_varW[SEX]) _varW[SEX] = new double[get_current_nbSamplePatch()];
    else if(get_current_nbSamplePatch()!=get_last_nbSamplePatch()){ // number of sampled patches may change over time
        delete[] _varW[SEX];
        _varW[SEX] = new double[get_current_nbSamplePatch()];
    }
    
    // for each population
    vector<TPatch*>::iterator curPop = get_vSamplePatch().begin();
    vector<TPatch*>::iterator endPop = get_vSamplePatch().end();
    for(unsigned int i=0; curPop!=endPop; ++curPop, ++i){ // for each patch
        setMeanAndVar_Wp_ofPatch(*curPop, _meanW[SEX][i], _varW[SEX][i], SEX);
    }
}

// ----------------------------------------------------------------------------------------
// setMeanAndVarFitness
// ----------------------------------------------------------------------------------------
void
TTQuantiSH::setMeanAndVar_Wp_ofPatch(TPatch* curPop, double& meanW, double& varW)
{
    vector<TIndividual*>& curFem = curPop->get_sampled_inds(FEM, ADLTx);
    vector<TIndividual*>& curMal = curPop->get_sampled_inds(MAL, ADLTx);
    unsigned int sizeF = (unsigned int)curFem.size(),
    sizeM = (unsigned int)curMal.size();
    unsigned int size  = sizeF + sizeM;
    
    // if the patch is empty or the phenotype is not yet computed -> stop
    if(!((sizeF && curFem.front()->getTraitFitness(_SHLinkedTraitIndex) != my_NAN)
         || (sizeM && curMal.front()->getTraitFitness(_SHLinkedTraitIndex) != my_NAN))){
        meanW = varW = my_NAN;
        return;
    }
    
    double* array = new double[size];
    
    vector<TIndividual*>::iterator curInd, endInd;
    unsigned int f=0;
    for(curInd=curFem.begin(), endInd=curFem.end(); curInd!=endInd; ++curInd, ++f) {
        array[f] = (*curInd)->getTraitFitness(_SHLinkedTraitIndex);
    }
    for(curInd=curMal.begin(), endInd=curMal.end(); curInd!=endInd; ++curInd, ++f) {
        array[f] = (*curInd)->getTraitFitness(_SHLinkedTraitIndex);
    }
    
    // compute mean and var
    meanW = ARRAY::mean(array, size);
    varW  = ARRAY::var(array, size, meanW);
    delete[] array;
}

// ----------------------------------------------------------------------------------------
// setMeanAndVarFitness
// ----------------------------------------------------------------------------------------
void
TTQuantiSH::setMeanAndVar_Wp_ofPatch(TPatch* curPop, double& meanW, double& varW, sex_t SEX)
{
    vector<TIndividual*>& cur = curPop->get_sampled_inds(SEX, ADLTx);
    unsigned int size = (unsigned int)cur.size();
    
    // if the patch is empty or the fitness is not yet computed -> stop
    if(!(size && cur.front()->getTraitFitness(_SHLinkedTraitIndex) != my_NAN)){
        meanW = varW = my_NAN;
        return;
    }
    
    double* array = new double[size];
    
    vector<TIndividual*>::iterator curInd, endInd;
    unsigned int f=0;
    for(curInd=cur.begin(), endInd=cur.end(); curInd!=endInd; ++curInd, ++f) {
        array[f] = (*curInd)->getTraitFitness(_SHLinkedTraitIndex);
    }
    
    // compute mean and var
    meanW = ARRAY::mean(array, size);
    varW  = ARRAY::var(array, size, meanW);
    delete[] array;
}

// ----------------------------------------------------------------------------------------
// getQst
// ----------------------------------------------------------------------------------------
/** calculation of QST:
 * QST = Vb/(Vb+2*h2*Vwp) = Vb/(Vb+2*Vwa)
 * h2  = narrow sense heritability (h2=Va/Vp)
 * Vb  = between pop variance of the phenotype       (variance of the phenotypic population means)
 * Vwp = within pop variance of the phenotype        (mean of the phenotypic population variances)
 * Vwa = within pop variance of the additive genetic (mean of the additive genotypic population variances)
 */
double
TTQuantiSH::getQst(const age_idx& AGE)
{
    setVar_Va(AGE);
    setMeanAndVar_Vp();
    
    double Vb = ARRAY::varUnbiased(_meanP[2], get_current_nbSamplePatch());    // round of problems
    if(Vb == my_NAN) return my_NAN;
    
    double Va = ARRAY::mean(_varA, get_current_nbSamplePatch());
    if(Va == my_NAN) return my_NAN;
    
    return  (2*Va+Vb) ? Vb/(2*Va+Vb) : my_NAN;
}

// ----------------------------------------------------------------------------------------
// setQst_perPatchPair
// ---------------------------------------------------------
/* h2  = narrow sense heritability (h2=Va/Vp)
 
 * Vb  = between pop variance of the phenotype       (variance of the phenotypic population means)
 * Vwp = within pop variance of the phenotype        (mean of the phenotypic population variances)
 * Vwa = within pop variance of the additive genetic (mean of the additive genotypic population variances)
 */
void
TTQuantiSH::setQst_perPatchPair(const age_idx& AGE)
{
    if(!_qst_matrix) ARRAY::create_2D(_qst_matrix, get_current_nbSamplePatch(), get_current_nbSamplePatch(), (double)my_NAN);
    else if(get_current_nbSamplePatch()!=get_last_nbSamplePatch()){ // number of sampled patches may change over time
        ARRAY::delete_2D(_qst_matrix, get_last_nbSamplePatch());
        ARRAY::create_2D(_qst_matrix, get_current_nbSamplePatch(), get_current_nbSamplePatch(), (double)my_NAN);
    }
    else ARRAY::reset_2D(_qst_matrix, get_current_nbSamplePatch(), get_current_nbSamplePatch(), (double)my_NAN);
    
    setVar_Va(AGE);
    setMeanAndVar_Vp();
    
    double Vb, Va;
    double array[2];		// temporary array
    
    // for each pair of pops
    unsigned int i, j;
    vector<TPatch*>::iterator curPop1, curPop2, endPop = get_vSamplePatch().end();
    for(i=0, curPop1=get_vSamplePatch().begin(); curPop1!=endPop; ++i, ++curPop1){ // for each patch
        if(_meanP[2][i] == my_NAN) continue;                   // both pops have to be populated
        array[0] = _meanP[2][i];
        
        curPop2=curPop1;
        for(j=i+1, ++curPop2; curPop2!=endPop; ++j, ++curPop2){ // for each patch
            assert(j==(*curPop2)->get_sampleID());
            if(_meanP[2][j] == my_NAN) continue;                   // both pops have to be populated
            array[1] = _meanP[2][j];
            Vb = ARRAY::var(array, 2);         		// var of means
            Va = (_varA[i] +_varA[j])/2.0;    // mean of vars
            
            _qst_matrix[i][j] = (2*Va+Vb) ? Vb/(2*Va+Vb) : my_NAN;
        }
    }
}

// ----------------------------------------------------------------------------------------
// getQst
// ----------------------------------------------------------------------------------------
/** calculation of QST with correction for inbreeding:
 * QST = (1+F)Vb/((1+F)Vb+2*h2*Vwp) = (1+F)Vb/((1+F)Vb+2*Vwa)
 * F   = inbreeding coefficient (=Fis)
 * h2  = narrow sense heritability (h2=Va/Vp)
 * Vb  = between pop variance of the phenotype       (variance of the phenotypic population means)
 * Vwp = within pop variance of the phenotype        (mean of the phenotypic population variances)
 * Vwa = within pop variance of the additive genetic (mean of the additive genotypic population variances)
 */
double
TTQuantiSH::getQstF(const age_idx& AGE)
{
    setVar_Va(AGE);
    setMeanAndVar_Vp();
    
    double Vb = ARRAY::varUnbiased(_meanP[2], get_current_nbSamplePatch());    // round off problems
    if(Vb == my_NAN) return my_NAN;
    
    double Va = ARRAY::mean(_varA, get_current_nbSamplePatch());
    if(Va == my_NAN) return my_NAN;
    
    double F = getFis(AGE);   // inbreeding coefficient
    
    return  ((1+F)*Vb+2*Va) ? (1+F)*Vb/((1+F)*Vb+2*Va) : my_NAN;
}

// ----------------------------------------------------------------------------------------
// setQst_perPatchPair
// ----------------------------------------------------------------------------------------
/** calculation of QST for all pairwise combinations with correction for inbreeding:
 * QST = (1+F)Vb/((1+F)Vb+2*h2*Vwp) = (1+F)Vb/((1+F)Vb+2*Vwa)
 * F   = inbreeding coefficient (=Fis)
 * h2  = narrow sense heritability (h2=Va/Vp)
 * Vb  = between pop variance of the phenotype       (variance of the phenotypic population means)
 * Vwp = within pop variance of the phenotype        (mean of the phenotypic population variances)
 * Vwa = within pop variance of the additive genetic (mean of the additive genotypic population variances)
 */
void
TTQuantiSH::setQstF_perPatchPair(const age_idx& AGE)
{
    if(!_qstF_matrix) ARRAY::create_2D(_qstF_matrix, get_current_nbSamplePatch(), get_current_nbSamplePatch(), (double)my_NAN);
    else if(get_current_nbSamplePatch()!=get_last_nbSamplePatch()){ // number of sampled patches may change over time
        ARRAY::delete_2D(_qstF_matrix, get_current_nbSamplePatch());
        ARRAY::create_2D(_qstF_matrix, get_current_nbSamplePatch(), get_current_nbSamplePatch(), (double)my_NAN);
    }
    else ARRAY::reset_2D(_qstF_matrix, get_current_nbSamplePatch(), get_current_nbSamplePatch(), (double)my_NAN);
    
    setVar_Va(AGE);
    setMeanAndVar_Vp();
    
    unsigned int i, j;
    double Vb, Va, F, H, hs, ho, hsnei;
    double array[2];		// temporare array
    
    // get Ho and Hs for each pop
    double * hs_pop = new double[get_current_nbSamplePatch()];
    double * ho_pop = new double[get_current_nbSamplePatch()];
    vector<TPatch*>::iterator curPop1, curPop2, endPop = get_vSamplePatch().end();
    for(i=0, curPop1=get_vSamplePatch().begin(); curPop1!=endPop; ++i, ++curPop1){
        ho_pop[i] = getHo_ofPatch(*curPop1, AGE);
        hs_pop[i] = getHs_ofPatch(*curPop1, AGE);
    }
    
    // for each pair of pops
    for(i=0, curPop1=get_vSamplePatch().begin(); curPop1!=endPop; ++i, ++curPop1){
        if(_meanP[2][i] == my_NAN) continue;                   // both pops have to be populated
        array[0] = _meanP[i][2];
        
        curPop2=curPop1;
        for(j=i+1, ++curPop2; curPop2!=endPop; ++j, ++curPop2){
            assert(j==(*curPop2)->get_sampleID());
            if(_meanP[2][j] == my_NAN) continue;                   // both pops have to be populated
            array[1] = _meanP[2][j];
            Vb = ARRAY::var(array, 2);         		// var of means
            Va = (_varA[i] +_varA[j])/2.0;        // mean of vars
            
            // inbreeding coefficient
            H = 2.0/((1.0/(*curPop1)->sampleSize(AGE))+(1.0/(*curPop2)->sampleSize(AGE))); // harmonic mean of N
            hs = (hs_pop[i] + hs_pop[j])/2.0;     // expected genetic diversity
            ho = (ho_pop[i] + ho_pop[j])/2.0;     // observed genetic diversity
            hsnei = H!=1 ? H/(H-1.0)*(hs-(ho/(2.0*H))) : hs;  //Nei's corrections
            F = hsnei ? 1.0-(ho/hsnei)    : 0;   // monomorphic = hs=0: after definition 0 and not NaN!
            
            _qstF_matrix[i][j] = ((1+F)*Vb+2*Va) ? (1+F)*Vb/((1+F)*Vb+2*Va) : my_NAN;
        }
    }
    delete[] hs_pop;
    delete[] ho_pop;
}

// ----------------------------------------------------------------------------------------
// getVgB
// ----------------------------------------------------------------------------------------
double
TTQuantiSH::getVgB(const age_idx& AGE)
{
    setMeanAndVar_Vg(AGE);
    return ARRAY::var(_meanG[2], get_current_nbSamplePatch());
}

// ----------------------------------------------------------------------------------------
// getVpB
// ----------------------------------------------------------------------------------------
double
TTQuantiSH::getVpB()
{
    setMeanAndVar_Vp();
    return ARRAY::var(_meanP[2], get_current_nbSamplePatch());
}
// ----------------------------------------------------------------------------------------
// getVaW
// ----------------------------------------------------------------------------------------
double
TTQuantiSH::getVaW(const age_idx& AGE)
{
    setVar_Va(AGE);
    return ARRAY::mean(_varA, get_current_nbSamplePatch());
}

// ----------------------------------------------------------------------------------------
// getVgW
// ----------------------------------------------------------------------------------------
double
TTQuantiSH::getVgW(const age_idx& AGE)
{
    setMeanAndVar_Vg(AGE);
    return ARRAY::mean(_varG[2], get_current_nbSamplePatch());
}

// ----------------------------------------------------------------------------------------
// getVpW
// ----------------------------------------------------------------------------------------
double
TTQuantiSH::getVpW()
{
    setMeanAndVar_Vp();
    return ARRAY::mean(_varP[2], get_current_nbSamplePatch());
}

