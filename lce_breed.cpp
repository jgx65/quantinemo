/** @file LCEbreed.cpp
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

#include "lce_breed.h"
#include "tselection.h"
#include "ttquanti.h"
#include "tselectiontrait.h"

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                             ******** LCE_Breed ********

// ----------------------------------------------------------------------------------------
// LCE_Breed_Base   NATING FUNCTIONS
// ----------------------------------------------------------------------------------------
// get a specific individual (no stochasity) /////////////////////////////////
Individual*
LCE_Breed::Index_MatingFunc   (Patch* thePatch, unsigned int& index, sex_t sex){
    return thePatch->get(sex, ADLTx, index);         // index has a meaning
}
Individual*        // for cloning
LCE_Breed::NULL_pointer   (Patch* thePatch, unsigned int& index, sex_t sex){
    return NULL;
}

// random mating /////////////////////////////////////////////////////////////
Individual*
LCE_Breed::Random_MatingFunc   (Patch* thePatch, unsigned int& index, sex_t sex){
    return thePatch->get(sex, ADLTx, get_pop_ptr()->rand().Uniform((unsigned int)_nbIndividuals[sex]));
}
Individual*
LCE_Breed::Random_Index_MatingFunc   (Patch* thePatch, unsigned int& index, sex_t sex){
    index = get_pop_ptr()->rand().Uniform((unsigned int)_nbIndividuals[sex]);
    return thePatch->get(sex, ADLTx, index);
}
Individual*
LCE_Breed::Random_S_MatingFunc (Patch* thePatch, unsigned int& index, sex_t sex){
    return _pSelection->get_RAND_mostFit(sex);
}
Individual*
LCE_Breed::Random_Index_S_MatingFunc (Patch* thePatch, unsigned int& index, sex_t sex){
    return _pSelection->get_RAND_mostFit_index(sex, index);   // index has a meaning
}

// full polygyny (one male) //////////////////////////////////////////////////
Individual*
LCE_Breed::fullPolygyny_oneMale_MatingFunc (Patch* thePatch, unsigned int& index, sex_t sex){
    return thePatch->get(sex, ADLTx, 0);    // return any individual (first one)  // to be checked !!!
}
Individual*
LCE_Breed::fullPolygyny_oneMale_S_MatingFunc (Patch* thePatch, unsigned int& index, sex_t sex){
    return _pSelection->get_mostFit(sex);         // return the fittest individual
}
Individual*
LCE_Breed::fullPolygyny_oneMale_S_MatingFunc2 (Patch* thePatch, unsigned int& index, sex_t sex){
    return _pSelection->get_RAND_mostFit(sex);    // return the randomly fittest individual
}

Individual*
LCE_Breed::fullPolygyny_manyMales_MatingFunc  (Patch* thePatch, unsigned int& index, sex_t sex){
    if(_nbIndividuals[sex]<_mating_males) return Random_MatingFunc(thePatch,index,sex);                    // random mating
    else return thePatch->get(sex, ADLTx, get_pop_ptr()->rand().Uniform(_mating_males));  // return any of the first x individuals
}
Individual*
LCE_Breed::fullPolygyny_manyMales_S_MatingFunc  (Patch* thePatch, unsigned int& index, sex_t sex){
    if(_nbIndividuals[sex]<_mating_males) return Random_S_MatingFunc(thePatch,index,sex);   // random mating
    else return _pSelection->get_RAND_mostFit_of_mostFit(sex, _mating_males);
}
Individual*
LCE_Breed::fullPolygyny_manyMales_S_MatingFunc2  (Patch* thePatch, unsigned int& index, sex_t sex){
    if(_nbIndividuals[sex]<_mating_males) return Random_S_MatingFunc(thePatch,index,sex);   // random mating
    else return _pSelection->get_RAND_mostFit_of_RAND_mostFit(sex, _mating_males);
}

// partial polygyny (one male) //////////////////////////////////////////////
Individual*
LCE_Breed::partialPolygyny_oneMale_MatingFunc (Patch* thePatch, unsigned int& index, sex_t sex){
    if(get_pop_ptr()->rand().Uniform()>_mating_proportion)return Random_MatingFunc(thePatch, index, sex);
    else return thePatch->get(sex, ADLTx, 0);   // to be checked!!!
}
Individual*
LCE_Breed::partialPolygyny_oneMale_S_MatingFunc (Patch* thePatch, unsigned int& index, sex_t sex){
    if(get_pop_ptr()->rand().Uniform()>_mating_proportion)return Random_S_MatingFunc(thePatch, index, sex);
    else return _pSelection->get_mostFit(sex);             // return the fittest individual
}
Individual*
LCE_Breed::partialPolygyny_oneMale_S_MatingFunc2 (Patch* thePatch, unsigned int& index, sex_t sex){
    if(get_pop_ptr()->rand().Uniform()>_mating_proportion)return Random_S_MatingFunc(thePatch, index, sex);
    else return _pSelection->get_RAND_mostFit(sex);        // return the ranom fittest individual
}

// partial polygyny (many males) /////////////////////////////////////////////
Individual*
LCE_Breed::partialPolygyny_manyMales_MatingFunc  (Patch* thePatch, unsigned int& index, sex_t sex){
    if(get_pop_ptr()->rand().Uniform()>_mating_proportion)return Random_MatingFunc(thePatch, index, sex);
    else return fullPolygyny_manyMales_MatingFunc(thePatch, index, sex);
}
Individual*
LCE_Breed::partialPolygyny_manyMales_S_MatingFunc  (Patch* thePatch, unsigned int& index, sex_t sex){
    if(get_pop_ptr()->rand().Uniform()>_mating_proportion)return Random_S_MatingFunc(thePatch, index, sex);
    else return fullPolygyny_manyMales_S_MatingFunc(thePatch, index, sex);
}
Individual*
LCE_Breed::partialPolygyny_manyMales_S_MatingFunc2  (Patch* thePatch, unsigned int& index, sex_t sex){
    if(get_pop_ptr()->rand().Uniform()>_mating_proportion)return Random_S_MatingFunc(thePatch, index, sex);
    else return fullPolygyny_manyMales_S_MatingFunc2(thePatch, index, sex);
}

// monogamy //////////////////////////////////////////////////////////////////
Individual*
LCE_Breed::Monogyny_MatingFunc (Patch* thePatch, unsigned int& index, sex_t sex){
    assert(_aMatingPairs[sex] && index<_aMatingPairs_size);
    if(get_pop_ptr()->rand().Uniform()>_mating_proportion || _nbIndividuals[sex]<index+1)
        return Random_MatingFunc(thePatch, index, sex);
    else
        return _aMatingPairs[sex][index];              // index has a meaning
}
Individual*
LCE_Breed::Monogyny_S_MatingFunc (Patch* thePatch, unsigned int& index, sex_t sex){
    assert(_aMatingPairs[sex] && index<_aMatingPairs_size);
    if(get_pop_ptr()->rand().Uniform()>_mating_proportion || _nbIndividuals[sex]<index+1)
        return Random_S_MatingFunc(thePatch, index, sex);
    else
        return _aMatingPairs[sex][index];             // index has a meaning
}

// one sex ///////////////////////////////////////////////////////////////////
Individual*
LCE_Breed::oneSex_notSameIndex_MatingFunc (Patch* thePatch, unsigned int& index, sex_t sex){
    if(_nbIndividuals[sex]<2) return thePatch->get(sex, ADLTx, 0);    // if there is only a single individual in the patch
    unsigned int newIndex;
    do {
        newIndex = get_pop_ptr()->rand().Uniform((unsigned int)_nbIndividuals[sex]);
    } while(index == newIndex);    // while it is the same index
    return thePatch->get(sex, ADLTx, newIndex);
}
Individual*
LCE_Breed::partialSelfing_MatingFunc (Patch* thePatch, unsigned int& index, sex_t sex){
    if(get_pop_ptr()->rand().Uniform()>_mating_proportion) return oneSex_notSameIndex_MatingFunc(thePatch, index, sex);  // random mating
    else return thePatch->get(sex, ADLTx, index);         // return the same female (index has a meaning)
}
Individual*
LCE_Breed::oneSex_notSameIndex_S_MatingFunc (Patch* thePatch, unsigned int& index, sex_t sex){
    if(_nbIndividuals[sex]<2) return thePatch->get(sex, ADLTx, 0);    // if there is only a single individual in the patch
    unsigned int newIndex, i=0;
    Individual* ind;
    do {
        ind = _pSelection->get_RAND_mostFit_index(sex, newIndex);
        if(++i > 1e5){ // for security reasons pick any individual
            ind = _pSelection->get_RAND_noFit_index(sex, newIndex);
        }
    } while(index == newIndex);  // while it is the same index
    return ind;
}
Individual*
LCE_Breed::partialSelfing_S_MatingFunc (Patch* thePatch, unsigned int& index, sex_t sex){
    if(get_pop_ptr()->rand().Uniform()>_mating_proportion) return oneSex_notSameIndex_S_MatingFunc(thePatch, index, sex); // random mating
    else return thePatch->get(sex, ADLTx, index);          // return the same female (index has a meaning)
}

// cloning ///////////////////////////////////////////////////////////////////
Individual*
LCE_Breed::partialCloning_MatingFunc (Patch* thePatch, unsigned int& index, sex_t sex){
    if(get_pop_ptr()->rand().Uniform()>_mating_proportion) return oneSex_notSameIndex_MatingFunc(thePatch, index, sex);  // random mating
    else return NULL;					 // cloning
}
Individual*
LCE_Breed::partialCloning_S_MatingFunc (Patch* thePatch, unsigned int& index, sex_t sex){
    if(get_pop_ptr()->rand().Uniform()>_mating_proportion) return oneSex_notSameIndex_S_MatingFunc(thePatch, index, sex); // random mating
    else return NULL;          // cloning
}



// ----------------------------------------------------------------------------------------
// LCE_Breed
// ----------------------------------------------------------------------------------------
LCE_Breed::LCE_Breed(int rank): LCE("breed","reproduction","",rank),
_aMatingPairs_size(0), _mating_system(0), _threshold(my_NAN), _fitness_factor_zero_isLethal(0),
_mating_proportion(1), _mean_fecundity(0), _growth_rate(0),
getMother_func_ptr(0), getFather_func_ptr(0), _maleSex(MAL), _pSelection(0){
    
    add_parameter("mating_system",INT2,false,0,6,"0", false,
                  "Available mating systems:\n" \
                  "  0: random mating (hermaphrodite)\n" \
                  "  1: selfing mating (hermaphrodite)\n" \
                  "  2: cloning mating (hermaphrodite)\n" \
                  "  3: random mating (promiscuity)\n" \
                  "  4: polygyny\n" \
                  "  5: monogamy\n" \
                  "  6: no mating/reproduction (used to compute only stats)\n" \
                  ,0);
    
    add_parameter("mating_proportion",DBL,false,0,1,"1", true,
                  "Proportion of special mating system in relation to random mating.",2);
    
    add_parameter("mating_males",INT2,false,my_NAN,my_NAN,"1", false,
                  "The number of males available for mating (only polygyny).",3);
    
    add_parameter("mean_fecundity",DBL,false,0,my_NAN,"0", true,
                  "Mean fecundity of females.",3);
    
    add_parameter("growth_rate",DBL,false,my_NAN,my_NAN,"0", true,
                  "Specifies the growth rate for logistic growth.",5);
    
    add_parameter("sex_ratio",DBL,false,0,my_NAN,"1", false,
                  "Sex ratio of males to females.",5);
    
    add_parameter("mating_nb_offspring_model",INT2,false,0,9,"0", false,
                  "How is the number of offspring computed:\n" \
                  "  0: carrying capacity (N_off = K)\n" \
                  "  1: keep number (N_off = N_adlt)\n" \
                  "  2: fecundity (N_off = Poisson(F_f*f))\n" \
                  "  3: fecundity simple (N_off = round(N_f*f))\n" \
                  "  4: fecundity binomial (N_off = floor(N_f*f) + Binomial(N_f*f - floor(N_f*f),1))\n" \
                  "  5: fecundity limited (same as 'fecundity' but " \
                  "population size may not exceed carrying capacity)\n" \
                  "  6: fecundity simple & limited (same as 'fecundtiy simple' but " \
                  "population size may not exceed carrying capacity)\n" \
                  "  7: fecundity binomial & limited (same as 'fecundity binomial' but " \
                  "population size may not exceed carrying capacity)\n" \
                  "  8: logistic regulation (N_off = N*K*(1+r)/(N*(1+r)-N+K))\n" \
                  "  9: stochastic logistic regulation (N_off = Poisson(N*K*(1+r)/(N*(1+r)-N+K))).\n" \
                  ,0);
    
    add_parameter("fitness_factor_zero_lethal",INT2,false,0,1,"0", false,
                  "How to thread a fitness of zero (0.0):\n" \
                  "  0: normal (may even reproduce)\n" \
                  "  1: lethal (a fitness factor of 0.0 is lethal).",4);
    
    
    add_parameter("sex_ratio_threshold",DBL,false,my_NAN,my_NAN,"0", false,
                  "If the parameter is set, the sex of an offspring is determined " \
                  "by the phenotype of the first quantitative trait (above " \
                  "threshold: male; below threshold: female).",4);
    
    
    add_parameter("fem_sex_allocation",STR,false,0,my_NAN,"0.5", false,
                  "Female sex allocation to sperm/pollen and ovules (0: only sperm/pollen; 1 (only ovules).",5);
    
    _aMatingPairs[MAL] = _aMatingPairs[FEM] = NULL;
}

// ----------------------------------------------------------------------------------------
// LCE_Breed::execute
// ----------------------------------------------------------------------------------------
void LCE_Breed::execute () {
#ifdef _DEBUG
    message("  LCE_Breed ");
    message("done! (Patches: %i; Individuals(F/M): [%i/%i; %i/%i])\n",
            _popPtr->get_nbFullPatch(),
            _popPtr->size(FEM, OFFSx),  _popPtr->size(MAL, OFFSx),
            _popPtr->size(FEM, ADLTx),  _popPtr->size(MAL, ADLTx));
    assert(_popPtr->individual_container_ok());
#endif
    
    (this->*breed)();
    
    if(_fitness_factor_zero_isLethal) remove_inds_zero_fitnessFactor(OFFSx);
    
    if(_threshold != my_NAN) reset_sex_after_phenotype(OFFSx);
    
#ifdef _DEBUG
    message("done! (Patches: %i; Individuals(F/M): [%i/%i; %i/%i])\n",
            _popPtr->get_nbFullPatch(),
            _popPtr->size(FEM, OFFSx),  _popPtr->size(MAL, OFFSx),
            _popPtr->size(FEM, ADLTx),  _popPtr->size(MAL, ADLTx));
    assert(_popPtr->individual_container_ok());
#endif
}

// ----------------------------------------------------------------------------------------
// temporal_change
// ----------------------------------------------------------------------------------------
/** parameters which may change over time */
void LCE_Breed::temporal_change(const unsigned int& gen)
{
    map<string, Param*>* pParam = _paramSet->getTemporalParams(gen);
    
    // if it is a temporal parameter
    if (pParam) {
        // check if a change has to be made
        map<string, Param*>* pMap = _paramSet->updateTemporalParams(gen);
        if (pMap) {
            // iterate through the map and perform the updates
            map<string, Param*>::iterator pos = pMap->begin();
            for (; pos != pMap->end(); ++pos) {
                if (pos->first == "mating_proportion")
                    _mating_proportion = pos->second->get_value();
                else if(pos->first == "mean_fecundity")
                    _mean_fecundity = pos->second->get_value();
                else if(pos->first == "growth_rate"){
                    _growth_rate = pos->second->get_value();
                    if((_nbOffspring_model==6 || _nbOffspring_model==7) &&
                       _growth_rate==my_NAN)
                        error("Mating system: The growth rate has to be specified!\n");
                }
             }
        }
    }
}

// ----------------------------------------------------------------------------------------
// LCE_Breed::init
// ----------------------------------------------------------------------------------------
bool LCE_Breed::init(Metapop* popPtr)
{
    LCE::init(popPtr);
    
    _mating_system      = (int)this->get_parameter_value("mating_system");
    _mating_proportion  = this->get_parameter_value("mating_proportion");
    _mean_fecundity     = this->get_parameter_value("mean_fecundity");
    _growth_rate        = this->get_parameter_value("growth_rate");
    _mating_males       = (int)_paramSet->getValue("mating_males");
    
    
    set_nbOffspring_model();
    
    // does selection acts at reproduction?
    _pSelection = _popPtr->get_pSelection();
    breed = &LCE_Breed::breed_neutral;                  // by default breeding is neutral
    if(_pSelection){                                    // if selection acts somewhere
        switch(_pSelection->get_selection_position()){
            case 0: // selection acts at the reproductive success
                switch(_pSelection->get_selection_level()){
                    case 0:                                                   // patch level
                    case 2: breed = &LCE_Breed::breed_soft2hard;     break;   // soft to hard selection
                    case 1: breed = &LCE_Breed::breed_soft2metapop;           // soft to metapop level
                        _popPtr->set_total_carrying_capacity();  break;
                    case 3: breed = &LCE_Breed::breed_metapop2hard;           // metapop 2 hard selection
                        _popPtr->set_total_carrying_capacity();  break;
                }
                break;
                
                
            case 1: // selection acts at the survival of the newborn
                switch(_pSelection->get_selection_level()){
                    case 0:                                                           // patch level
                    case 2: breed = &LCE_Breed::breed_offspring_soft2hard;     break; // soft to hard selection
                    case 1: breed = &LCE_Breed::breed_offspring_soft2metapop;         // soft to metapop level
                        _popPtr->set_total_carrying_capacity();            break;
                    case 3: breed = &LCE_Breed::breed_offspring_metapop2hard;         // metapop 2 hard selection
                        _popPtr->set_total_carrying_capacity();            break;
                }
                break;
        }
    }
    
    //set the mating function ptr:
    _sort[MAL] = _sort[FEM] = -1;
    if(_pSelection && _pSelection->get_selection_position()==0){ // selection acts at the reproductive success of the parents
        switch(_mating_system) {
                // with selection
            default:
            case 0:    // random mating hermaphrodite
                _maleSex = FEM;    // only one sex
                getMother_func_ptr = &LCE_Breed::Random_S_MatingFunc;
                getFather_func_ptr = &LCE_Breed::Random_S_MatingFunc;
                break;
            case 1:    // selfing hermaphrodite
                _maleSex = FEM;    // only one sex
                getMother_func_ptr = &LCE_Breed::Random_Index_S_MatingFunc;
                if(_mating_proportion == 1){
                    getFather_func_ptr = &LCE_Breed::Index_MatingFunc;        // get exactly the same individual
                }
                else {
                    getFather_func_ptr = &LCE_Breed::partialSelfing_S_MatingFunc;
                }
                break;
            case 2:    // cloning
                _maleSex = FEM;    // only one sex
                getMother_func_ptr = &LCE_Breed::Random_Index_S_MatingFunc;
                if(_mating_proportion == 1){
                    getFather_func_ptr = &LCE_Breed::NULL_pointer;        // get no father, i.e. NULL (cloning)
                }
                else {
                    getFather_func_ptr = &LCE_Breed::partialCloning_S_MatingFunc;
                }
                break;
            case 3:    // promiscuity/random mating.
                getMother_func_ptr = &LCE_Breed::Random_S_MatingFunc;
                getFather_func_ptr = &LCE_Breed::Random_S_MatingFunc;
                break;
            case 4:    // polygyny (random subset after the fitness of the males)
                _sort[MAL] = -3;
                getMother_func_ptr = &LCE_Breed::Random_S_MatingFunc;
                if(_mating_proportion == 1){
                    if(_mating_males == 1)  getFather_func_ptr = &LCE_Breed::fullPolygyny_oneMale_S_MatingFunc2;
                    else                    getFather_func_ptr = &LCE_Breed::fullPolygyny_manyMales_S_MatingFunc2;
                }
                else{
                    if(_mating_males == 1)  getFather_func_ptr = &LCE_Breed::partialPolygyny_oneMale_S_MatingFunc2;
                    else                    getFather_func_ptr = &LCE_Breed::partialPolygyny_manyMales_S_MatingFunc2;
                }
                break;
            case 5:    // monogamy
                getMother_func_ptr = &LCE_Breed::Random_Index_S_MatingFunc;
                getFather_func_ptr = &LCE_Breed::Monogyny_S_MatingFunc;
                break;
            case 6:   	// no mating and reproduction occurs (used to just compute stats from a ini genotype file)
                breed = &LCE_Breed::breed_none;
                if(_popPtr->getGenerations()>1) error("No mating/reproduction takes place thus only a single generation makes sense!\n");
                break;
        }
    }
    else{        // neutral mating (selection may act at a different stage)
        switch(_mating_system) {
            default:
            case 0:    // random mating hermaphrodite
                _maleSex = FEM;    // only one sex
                getMother_func_ptr = &LCE_Breed::Random_MatingFunc;
                getFather_func_ptr = &LCE_Breed::Random_MatingFunc;
                break;
            case 1:    // selfing hermaphrodite
                _maleSex = FEM;    // only one sex
                getMother_func_ptr = &LCE_Breed::Random_Index_MatingFunc;
                if(_mating_proportion == 1){
                    getFather_func_ptr = &LCE_Breed::Index_MatingFunc;        // get exactly the same individual
                }
                else {
                    getFather_func_ptr = &LCE_Breed::partialSelfing_MatingFunc;
                }
                break;
            case 2:    // cloning
                _maleSex = FEM;    // only one sex
                getMother_func_ptr = &LCE_Breed::Random_Index_MatingFunc;
                if(_mating_proportion == 1){
                    getFather_func_ptr = &LCE_Breed::NULL_pointer;        // get no father, i.e. NULL (cloning)
                }
                else {
                    getFather_func_ptr = &LCE_Breed::partialCloning_MatingFunc;
                }
                break;
            case 3:     // promiscuity/random mating
                getMother_func_ptr = &LCE_Breed::Random_MatingFunc;
                getFather_func_ptr = &LCE_Breed::Random_MatingFunc;
                break;
            case 4:     // polygyny
            {
                getMother_func_ptr = &LCE_Breed::Random_MatingFunc;
                if(_mating_proportion == 1){
                    if(_mating_males == 1)  getFather_func_ptr = &LCE_Breed::fullPolygyny_oneMale_MatingFunc;
                    else                    getFather_func_ptr = &LCE_Breed::fullPolygyny_manyMales_MatingFunc;
                }
                else{
                    if(_mating_males == 1)  getFather_func_ptr = &LCE_Breed::partialPolygyny_oneMale_MatingFunc;
                    else                    getFather_func_ptr = &LCE_Breed::partialPolygyny_manyMales_MatingFunc;
                }
                break;
            }
            case 5:     // monogamy
                getMother_func_ptr = &LCE_Breed::Random_Index_MatingFunc;
                getFather_func_ptr = &LCE_Breed::Monogyny_MatingFunc;
                break;
            case 6:   	// no mating and reproduction occurs (used to just compute stats from a ini genotype file)
                breed = &LCE_Breed::breed_none;
                if(_popPtr->getGenerations()>1) error("No mating/reproduction takes place thus only a single generation makes sense!\n");
                break;
        }
    }
    
    // sex ratio / number of sexes
    if(_popPtr->get_sexInitRatio()){															// what is the intial sex ratio: 0: only females
        isMatingPossible_func_ptr = &LCE_Breed::isMatingPossible_2_sex;
        // is a sex ratio specified?
        _sex_ratio = this->get_parameter_value("sex_ratio");        // input: males/females
        if(_sex_ratio<0){     // is the sex ratio within the limits?
            setSexRatio_func_ptr  = &LCE_Breed::setSexRatio_KeepSexRatio;
            getRandomSex_func_ptr = &LCE_Breed::getRandomSex_KeepSexRatio;
        }
        else{
            _sex_ratio /= _sex_ratio+1;                            // now:   males/(males+females)
            setSexRatio_func_ptr  = &LCE_Breed::setSexRatio_NoSelfing;
            getRandomSex_func_ptr = &LCE_Breed::getRandomSex_NoSelfing;
        }
    }
    else {  // hermaphrodites (a single sex is simulated)
        isMatingPossible_func_ptr = &LCE_Breed::isMatingPossible_1_sex;
        setSexRatio_func_ptr  = &LCE_Breed::setSexRatio_Selfing;
        getRandomSex_func_ptr = &LCE_Breed::getRandomSex_Selfing;
    }
    
    // if the sex is determined by the phenotype of the first quantitative trait
    if(this->get_parameter("sex_ratio_threshold")->isSet()){
        _threshold = this->get_parameter_value("sex_ratio_threshold");
        
        // we need two sexes
        if(isMatingPossible_func_ptr != &LCE_Breed::isMatingPossible_2_sex) error("Two sexes are requiered when using a sex threshold!\n");
        
        // the phenotype of the first trait has to be set just after the birth to set the sex
        if(!_pSelection) error("Sex-chromosome: At least a quantitative trait is needed to specify the sex!\n");
        if(_pSelection->get_selTrait(0)->get_Ve_prop()) error("Sex-chromosome: only the natal environment can be taken into account!\n");
    }
    else _threshold = my_NAN;       // not used, simulations are normal
    
    // is a fitness factor of zero lethal?
    _fitness_factor_zero_isLethal = (unsigned int)this->get_parameter_value("fitness_factor_zero_lethal");
    if(_fitness_factor_zero_isLethal){  // is the fitness factor indeed used?
        if(!_pSelection) _fitness_factor_zero_isLethal = false;     // selection is needed
        else{                                                       // check if fitness factors are set
            vector<int> vTraits = _popPtr->getTraitIndex("quanti");   // get all quanti traits
            vector<int>::iterator curTrait, endTrait = vTraits.end();
            for(curTrait = vTraits.begin(); curTrait!=endTrait; ++curTrait){   // for each quanti trait
                assert((dynamic_cast<TTQuantiProto*>(&_popPtr->getTraitPrototype(*curTrait))));
                if((dynamic_cast<TTQuantiProto*>(&_popPtr->getTraitPrototype(*curTrait)))->fitnessFactor_used()) break;
            }
            if(curTrait==endTrait) _fitness_factor_zero_isLethal = false; // no fitness factors set
        }
    }
    
    // evolution of female sex allocation
    _fem_sex_allocation = NULL;
    //_fem_sex_allocation_all_trait = my_NAN;
    //_fem_sex_allocation_prop = my_NAN;
    //if(_pSelection) _pSelection->set_fem_sex_allocation(1, my_NAN);
    if(_pSelection) _pSelection->set_fem_sex_allocation(NULL);
    if(this->get_parameter_isSet("fem_sex_allocation")){
        vector<int> traitID = get_pop_ptr()->getTraitIndex("quanti");	// containing the absolute index across all types of traits
        _fem_sex_allocation = new TEquation(this->get_parameter("fem_sex_allocation")->get_arg(), traitID);
        if(_pSelection) _pSelection->set_fem_sex_allocation(_fem_sex_allocation);
        isMatingPossible_func_ptr = &LCE_Breed::isMatingPossible_sex_allocation_equation;
        
        if(_pSelection && _pSelection->get_selection_position()!=0)
            error("Parameter fem_sex_allocation requires selection at the repoductive success or no selection!\n");
        
        setSexRatio_func_ptr  = &LCE_Breed::setSexRatio_Selfing; // produce just daugthers
        getRandomSex_func_ptr = &LCE_Breed::getRandomSex_Selfing;
        
        _maleSex = MAL;                 // we need hermaphrodites butneverthelss we "mis"-use the males ;-)
        
        switch(_mating_system) {        // we have to use selection at mating even for neutral simulation
            default: error("Parameter fem_sex_allocation is only valid with hermaphrodites!\n"); break;
            case 0:    // random mating hermaphrodite
                getMother_func_ptr = &LCE_Breed::Random_S_MatingFunc;
                getFather_func_ptr = &LCE_Breed::Random_S_MatingFunc;
                break;
            case 1:    // selfing hermaphrodite
                getMother_func_ptr = &LCE_Breed::Random_Index_S_MatingFunc;
                if(_mating_proportion == 1) getFather_func_ptr = &LCE_Breed::Index_MatingFunc;        // get exactly the same individual
                else  getFather_func_ptr = &LCE_Breed::partialSelfing_S_MatingFunc;
                break;
        }
    }
    
    return true;
}

// ----------------------------------------------------------------------------------------
// set_nbOffspring_model
// ----------------------------------------------------------------------------------------
// how should the number of offspring be computed?
void LCE_Breed::set_nbOffspring_model()
{
    _nbOffspring_model = (unsigned int)_paramSet->getValue("mating_nb_offspring_model");
    switch(_nbOffspring_model){
        case 0: setNbOffspring_func_ptr = &LCE_Breed::setNbOffspring_CarryCapacity;  break;
        case 1: setNbOffspring_func_ptr = &LCE_Breed::setNbOffspring_KeepNb;         break;
        case 2: setNbOffspring_func_ptr = &LCE_Breed::setNbOffspring_Fecundity;
            if(_mean_fecundity==my_NAN) error("\nMating system: The fecundity has to be specified!");
            break;
        case 3: setNbOffspring_func_ptr = &LCE_Breed::setNbOffspring_FecunditySimple;
            if(_mean_fecundity==my_NAN) error("\nMating system: The fecundity has to be specified!");
            break;
        case 4: setNbOffspring_func_ptr = &LCE_Breed::setNbOffspring_FecundityBinomial;
            if(_mean_fecundity==my_NAN) error("\nMating system: The fecundity has to be specified!");
            break;
        case 5: setNbOffspring_func_ptr = &LCE_Breed::setNbOffspring_FecundityLimited;
            if(_mean_fecundity==my_NAN) error("\nMating system: The fecundity has to be specified!");
            break;
        case 6: setNbOffspring_func_ptr = &LCE_Breed::setNbOffspring_FecunditySimpleLimited;
            if(_mean_fecundity==my_NAN) error("\nMating system: The fecundity has to be specified!");
            break;
        case 7: setNbOffspring_func_ptr = &LCE_Breed::setNbOffspring_FecundityBinomialLimited;
            if(_mean_fecundity==my_NAN) error("\nMating system: The fecundity has to be specified!");
            break;
        case 8: setNbOffspring_func_ptr = &LCE_Breed::setNbOffspring_Logistic;
            if(_growth_rate==my_NAN) error("\nMating system: The growth rate has to be specified!");
            break;
        case 9: setNbOffspring_func_ptr = &LCE_Breed::setNbOffspring_RandLogistic;
            if(_growth_rate==my_NAN) error("\nMating system: The growth rate has to be specified!");
            break;
    }
}

// ----------------------------------------------------------------------------------------
// executeBeforeEachReplicate
// ----------------------------------------------------------------------------------------
void LCE_Breed::executeBeforeEachReplicate(const unsigned int& rep)
{
    if(_threshold != my_NAN) reset_sex_after_phenotype(ADLTx); // that has to be done just before the start
    
    // if fitness factor = 0 is lethal, then fitness factor is set after breeding
    // therefore for the first generation the fitness factors of the adults have to be set in advance
    if(_fitness_factor_zero_isLethal) remove_inds_zero_fitnessFactor(ADLTx);
}

// ----------------------------------------------------------------------------------------
// breed_offspring_soft2metapop
// ----------------------------------------------------------------------------------------
/** This breeding function allows to simulate any degree of selection between soft and metapop selection.
 * First mating occurs purely neutral with nbOff = fecundity*nbFem
 * Second selection acts depending on the mating_nb_offspring_model to downregulate the population size
 * nbOff = (1-f)*Kp  + f*(Wp/Wm)*Km
 *   f=0: soft selection
 *   f=1: metapop selection
 */
void LCE_Breed::breed_offspring_soft2metapop ()
{
#ifdef _DEBUG
    message("(soft-metapop (selection on offspring)) ... ");
#endif
    
    unsigned int nbBaby, i, nbPatches = _popPtr->get_nbFullPatch(), nbMal, nbFem;
    unsigned int nbDaughters, nbSons, nbOff;
    
    // array to buffer the fitness arrays for each patch (female and male)
    double totFitness = 0;                             // total fitness
    double* sumFitness = new double[nbPatches];        // sum of the fitness of each patch
    unsigned nbInd[2];                                 // total number of adult females and males
    TPatchFitness** fitnessArrays[2];
    for(i=0; i<2; ++i){
        fitnessArrays[i] = new TPatchFitness*[nbPatches];
        nbInd[i] = 0;
    }
    
    // create the offspring with neutral mating and store their fitnesses
    vector<Patch*>::iterator curPop = _popPtr->get_vFullPatch().begin();
    vector<Patch*>::iterator endPop = _popPtr->get_vFullPatch().end();
    for(i=0; curPop!=endPop; ++curPop, ++i) {
        // check if the mating requirements are met
        if(!(this->*isMatingPossible_func_ptr)(*curPop)){
            fitnessArrays[MAL][i] = NULL;
            fitnessArrays[FEM][i] = NULL;
            continue;
        }
        
        // create the mating pairs
        preMating(*curPop);
        
        // create the offspring with neutral mating
        nbBaby = get_pop_ptr()->rand().Poisson(_mean_fecundity*_nbIndividuals[FEM]);
        (this->*setSexRatio_func_ptr)(nbBaby, nbSons, nbDaughters, _nbIndividuals[MAL], _nbIndividuals[FEM]);
        createOffspring(*curPop, nbDaughters, nbSons);
        nbInd[MAL] += _nbIndividuals[MAL];      // get total number of male adults
        nbInd[FEM] += _nbIndividuals[FEM];      // get total number of female adults
        
        // compute and store the OFFSPRING fitness of the current patch (0: male, 1: female)
        _pSelection->set_fitness(*curPop, OFFSx);  // do not sort or make the array cumulative
        sumFitness[i]   = _pSelection->getSumFitness(MAL) + _pSelection->getSumFitness(FEM);
        totFitness      += sumFitness[i];
        if(nbSons)      fitnessArrays[MAL][i] = _pSelection->get_FitnessObject(MAL, true);
        else            fitnessArrays[MAL][i] = NULL;
        if(nbDaughters) fitnessArrays[FEM][i] = _pSelection->get_FitnessObject(FEM, true);
        else            fitnessArrays[FEM][i] = NULL;
    }
    
    // compute the total number of offspring of the metapopulation (taking the adults as reference...)
    nbOff = (this->*setNbOffspring_func_ptr)(nbInd[MAL], nbInd[FEM], _popPtr->get_total_carrying_capacity());
    
    // population size regulation depending on the fitness of the offspring?
    for(i=0, curPop = _popPtr->get_vFullPatch().begin(); curPop!=endPop; ++curPop, ++i) {
        if(!fitnessArrays[FEM][i]) continue;  // female missing -> no offspring
        
        // get the number of sons/daughters to produce of this patch (in relation to the adults)
        nbFem = (*curPop)->size(FEM, ADLTx);
        nbMal = (*curPop)->size(MAL, ADLTx);
        if(totFitness) nbBaby = _pSelection->get_SoftHardSelection(nbMal+nbFem, (unsigned int)(sumFitness[i]/totFitness*nbOff));
        else           nbBaby = _pSelection->get_SoftHardSelection(nbMal+nbFem, (unsigned int)((nbMal+nbFem)/nbOff));
        (this->*setSexRatio_func_ptr)(nbBaby, nbSons, nbDaughters, nbMal, nbFem);
        
        // perform population size regulation for the males
        if(nbSons<(*curPop)->size(MAL, OFFSx)){              	// are there too many sons?
            _pSelection->set_FitnessObject(MAL, fitnessArrays[MAL][i]);  // set the fitness arrays of the current patch
            (*curPop)->regulate_selection_fitness(nbSons, _pSelection, MAL, OFFSx);
        }
        
        // perform population size regulation for the females
        if(nbDaughters<(*curPop)->size(FEM, OFFSx)){        	// are there too many daughters?
            _pSelection->set_FitnessObject(FEM, fitnessArrays[FEM][i]);  // set the fitness arrays of the current patch
            (*curPop)->regulate_selection_fitness(nbDaughters, _pSelection, FEM, OFFSx);
        }
    }//end_for_nbPatch
    
    // delete the fitness arrays
    if(sumFitness)      delete[] sumFitness;
    for(i=0; i<2; ++i){
        if(fitnessArrays[i]) delete[] fitnessArrays[i];
    }
}

// ----------------------------------------------------------------------------------------
// breed_offspring_metapop2hard
// ----------------------------------------------------------------------------------------
/** This breeding function allows to simulate any degree of selection between metapop and hard-selection.
 * First mating occurs purely neutral with nbOff = fecundity*nbFem
 * Second selection acts depending on the mating_nb_offspring_model to downregulate the population size
 * nbOff = (1-f)*(Wp/Wm)*Km  + f*Wp*Kp
 *   f=0: metapop selection
 *   f=1: hard selection
 */
void LCE_Breed::breed_offspring_metapop2hard ()
{
#ifdef _DEBUG
    message("(metapop-hard (selection on offspring)) ... ");
#endif
    
    unsigned int nbBaby, i, nbMal, nbFem, nbPatches = _popPtr->get_nbFullPatch();
    unsigned int nbDaughters, nbSons, nbOff;
    
    // array to buffer the fitness arrays for each patch (female and male)
    double totFitness = 0;                             // total fitness
    double* sumFitness = new double[nbPatches];        // sum of the fitness of each patch
    unsigned nbInd[2];                                 // total number of adult females and males
    TPatchFitness** fitnessArrays[2];
    for(i=0; i<2; ++i){
        fitnessArrays[i] = new TPatchFitness*[nbPatches];
        nbInd[i] = 0;
    }
    
    // create the offspring with neutral mating and store their fitnesses
    vector<Patch*>::iterator curPop = _popPtr->get_vFullPatch().begin();
    vector<Patch*>::iterator endPop = _popPtr->get_vFullPatch().end();
    for(i=0; curPop!=endPop; ++curPop, ++i) {
        // check if the mating requirements are met
        if(!(this->*isMatingPossible_func_ptr)(*curPop)){
            fitnessArrays[MAL][i] = NULL;
            fitnessArrays[FEM][i] = NULL;
            continue;
        }
        
        // create the mating pairs
        preMating(*curPop);
        
        // create the offspring with neutral mating
        nbBaby = get_pop_ptr()->rand().Poisson(_mean_fecundity*_nbIndividuals[FEM]);
        (this->*setSexRatio_func_ptr)(nbBaby, nbSons, nbDaughters, _nbIndividuals[MAL], _nbIndividuals[FEM]);
        createOffspring(*curPop, nbDaughters, nbSons);
        nbInd[MAL] += _nbIndividuals[MAL];      // get total number of male adults
        nbInd[FEM] += _nbIndividuals[FEM];      // get total number of female adults
        
        // compute and store the OFFSPRING fitness of the current patch (0: male, 1: female)
        _pSelection->set_fitness(*curPop, OFFSx);  // do not sort or make the array cumulative
        sumFitness[i]   = _pSelection->getSumFitness(MAL) + _pSelection->getSumFitness(FEM);
        totFitness      += sumFitness[i];
        if(nbSons)      fitnessArrays[MAL][i] = _pSelection->get_FitnessObject(MAL, true);
        if(nbDaughters) fitnessArrays[FEM][i] = _pSelection->get_FitnessObject(FEM, true);
    }
    
    // compute the total number of offspring of the metapopulation (taking the adults as reference...)
    nbOff = (this->*setNbOffspring_func_ptr)(nbInd[MAL], nbInd[FEM], _popPtr->get_total_carrying_capacity());
    
    // population size regulation depending on the fitness of the offspring?
    for(i=0, curPop = _popPtr->get_vFullPatch().begin(); curPop!=endPop; ++curPop, ++i) {
        if(!fitnessArrays[FEM][i]) continue;  // female missing -> no offspring
        
        // get the number of sons/daughters to produce of this patch (in relation to the adults)
        nbFem = (*curPop)->size(FEM, ADLTx);
        nbMal = (*curPop)->size(MAL, ADLTx);
        if(totFitness) nbBaby = _pSelection->get_SoftHardSelection((unsigned int)(sumFitness[i]/totFitness*nbOff), (unsigned int)(_pSelection->getMeanFitness()*(nbMal+nbFem)));
        else nbBaby = _pSelection->get_SoftHardSelection((unsigned int)((nbMal+nbFem)/nbOff), (unsigned int)(_pSelection->getMeanFitness()*(nbMal+nbFem)));
        (this->*setSexRatio_func_ptr)(nbBaby, nbSons, nbDaughters, nbMal, nbFem);
        
        // perform population size regulation for the males
        if(nbSons<(*curPop)->size(MAL, OFFSx)){              	// are there too many sons?
            _pSelection->set_FitnessObject(MAL, fitnessArrays[MAL][i]);  // set the fitness arrays of the current patch
            (*curPop)->regulate_selection_fitness(nbSons, _pSelection, MAL, OFFSx);
        }
        
        // perform population size regulation for the females
        if(nbDaughters<(*curPop)->size(FEM, OFFSx)){        	// are there too many daughters?
            _pSelection->set_FitnessObject(FEM, fitnessArrays[FEM][i]);  // set the fitness arrays of the current patch
            (*curPop)->regulate_selection_fitness(nbDaughters, _pSelection, FEM, OFFSx);
        }
    }//end_for_nbPatch
    
    // delete the fitness arrays
    if(sumFitness)      delete[] sumFitness;
    for(i=0; i<2; ++i){
        if(fitnessArrays[i]) delete[] fitnessArrays[i];
    }
}

// ----------------------------------------------------------------------------------------
// breed_offspring_soft2hard
// ----------------------------------------------------------------------------------------
/** This breeding function allows to simulate any degree of selection between soft and hard selection.
 * First mating occurs purely neutral with nbOff = fecundity*nbFem
 * Second selection acts depending on the mating_nb_offspring_model to downregulate the population size
 * nbOff = (1-f)*Kp  + f*Wp*Kp
 *   f=0: soft selection
 *   f=1: hard selection
 */
void LCE_Breed::breed_offspring_soft2hard ()
{
#ifdef _DEBUG
    message("(soft-hard (selection on offspring)) ... ");
#endif
    
    unsigned int nbBaby, nbSons, nbDaughters, Kp;
    
    // for each patch
    vector<Patch*>::iterator curPop = _popPtr->get_vFullPatch().begin();
    vector<Patch*>::iterator endPop = _popPtr->get_vFullPatch().end();
    for(; curPop!=endPop; ++curPop) {
        // check if the mating requirements are met
        if(!(this->*isMatingPossible_func_ptr)(*curPop)) continue;
        
        // create the mating pairs
        preMating(*curPop);
        
        // create the offspring with neutral mating depending on the fecundity
        nbBaby = get_pop_ptr()->rand().Poisson(_mean_fecundity*_nbIndividuals[FEM]);
        (this->*setSexRatio_func_ptr)(nbBaby, nbSons, nbDaughters, _nbIndividuals[MAL], _nbIndividuals[FEM]);
        createOffspring(*curPop, nbDaughters, nbSons);
        
        // compute the mean fitness of the current patch (offspring, 0: male, 1: female)
        _pSelection->set_fitness(*curPop, OFFSx); 	// do not sort or make the array cumulative
        
        // get the number of sons/daughters to produce of this patch (in relation to adults)
        Kp = (this->*setNbOffspring_func_ptr)(_nbIndividuals[MAL], _nbIndividuals[FEM], (*curPop)->get_K());
        nbBaby = _pSelection->get_SoftHardSelection(Kp, (unsigned int)(_pSelection->getMeanFitness()*Kp));
        (this->*setSexRatio_func_ptr)(nbBaby, nbSons, nbDaughters, _nbIndividuals[MAL], _nbIndividuals[FEM]);
        
        // regulate pop size
        (*curPop)->regulate_selection_fitness(nbSons,      _pSelection, MAL, OFFSx);
        (*curPop)->regulate_selection_fitness(nbDaughters, _pSelection, FEM, OFFSx);
    }//end_for_nbPatch
}

// ----------------------------------------------------------------------------------------
// breed_offspring_soft2hard
// ----------------------------------------------------------------------------------------
/** This breeding function allows to simulate any degree of selection between soft and hard selection.
 * First mating occurs purely neutral with nbOff = fecundity*nbFem
 * Second selection acts depending on the mating_nb_offspring_model to downregulate the population size
 * nbOff = (1-f)*Kp  + f*Wp*Kp
 *   f=0: soft selection
 *   f=1: hard selection
 */
void LCE_Breed::breed_offspring_soft2hard2 ()
{
#ifdef _DEBUG
    message("(soft-hard (selection on offspring)) ... ");
#endif
    
    unsigned int nbBaby, nbSons, nbDaughters, Kp;
    
    // for each patch
    vector<Patch*>::iterator curPop = _popPtr->get_vFullPatch().begin();
    vector<Patch*>::iterator endPop = _popPtr->get_vFullPatch().end();
    for(; curPop!=endPop; ++curPop) {
        // check if the mating requirements are met
        if(!(this->*isMatingPossible_func_ptr)(*curPop)) continue;
        
        // create the mating pairs
        preMating(*curPop);
        
        // create the offspring with neutral mating depending on the fecundity
        nbBaby = get_pop_ptr()->rand().Poisson(_mean_fecundity*_nbIndividuals[FEM]);
        (this->*setSexRatio_func_ptr)(nbBaby, nbSons, nbDaughters, _nbIndividuals[MAL], _nbIndividuals[FEM]);
        createOffspring(*curPop, nbDaughters, nbSons);
        
        // compute the mean fitness of the current patch (offspring, 0: male, 1: female)
        _pSelection->set_fitness(*curPop, OFFSx); 	// do not sort or make the array cumulative
        
        // get the number of sons/daughters to produce of this patch (in relation to adults)
        Kp = (this->*setNbOffspring_func_ptr)(_nbIndividuals[MAL], _nbIndividuals[FEM], (*curPop)->get_K());
        nbBaby = _pSelection->get_SoftHardSelection(Kp, (unsigned int)(_pSelection->getMeanFitness()*Kp));
        (this->*setSexRatio_func_ptr)(nbBaby, nbSons, nbDaughters, _nbIndividuals[MAL], _nbIndividuals[FEM]);
        
        // regulate pop size
        (*curPop)->regulate_selection_fitness(nbSons,      _pSelection, MAL, OFFSx);
        (*curPop)->regulate_selection_fitness(nbDaughters, _pSelection, FEM, OFFSx);
    }//end_for_nbPatch
}

// ----------------------------------------------------------------------------------------
// breed_neutral
// ----------------------------------------------------------------------------------------
/** This breeding function is purely neutral
 */
void LCE_Breed::breed_neutral ()
{
#ifdef _DEBUG
    message("(neutral) ... ");
#endif
    
    unsigned int nbBaby, nbSons, nbDaughters;
    
    // for each patch
    vector<Patch*>::iterator curPop = _popPtr->get_vFullPatch().begin();
    vector<Patch*>::iterator endPop = _popPtr->get_vFullPatch().end();
    for(; curPop!=endPop; ++curPop) {
        // check if the mating requirements are met
        if(!(this->*isMatingPossible_func_ptr)(*curPop)) continue;
        
        // create the mating pairs
        preMating(*curPop);
        
        // get the number of sons/daughters to produce of this patch
        nbBaby = (this->*setNbOffspring_func_ptr)(_nbIndividuals[MAL], _nbIndividuals[FEM], (*curPop)->get_K());
        (this->*setSexRatio_func_ptr)(nbBaby, nbSons, nbDaughters, _nbIndividuals[MAL], _nbIndividuals[FEM]);
        
        createOffspring(*curPop, nbDaughters, nbSons);
    }//end_for_nbPatch
}

// ----------------------------------------------------------------------------------------
// breed_soft2metapop
// ----------------------------------------------------------------------------------------
/** This breeding function allows to simulate any degree of selection between soft and metapop selection.
 * This fucntion is not optimal to simulate just soft selection.
 * nbOff = (1-f)*Kp  + f*(Wp/Wm)*Km
 *   f=0: soft selection
 *   f=1: metapop selection
 */
void LCE_Breed::breed_soft2metapop ()
{
#ifdef _DEBUG
    message("(soft-metapop) ... ");
#endif
    
    assert(_pSelection);
    
    unsigned int nbBaby, i, nbPatches = _popPtr->get_nbFullPatch();
    unsigned int nbDaughters, nbSons, nbOff;
    
    // array to buffer the fitness arrays for each patch (female and male)
    double totFitness = 0;                              // total fitness
    double* sumFitness = new double[nbPatches];        // mean fitness of each patch
    unsigned nbInd[2];
    TPatchFitness** fitnessArrays[2];
    for(i=0; i<2; ++i){
        fitnessArrays[i] = new TPatchFitness*[nbPatches];
        nbInd[i] = 0;
    }
    
    // for each patch compute the fitness of the individuals
    vector<Patch*>::iterator curPop = _popPtr->get_vFullPatch().begin();
    vector<Patch*>::iterator endPop = _popPtr->get_vFullPatch().end();
    for(i=0; curPop!=endPop; ++curPop, ++i) {
        
        // check if the mating requirements are met
        if(!(this->*isMatingPossible_func_ptr)(*curPop)){
            fitnessArrays[MAL][i] = NULL;
            fitnessArrays[FEM][i] = NULL;
            continue;
        }
        
        // compute and store the fitness of the current patch (0: male, 1: female)
        _pSelection->set_fitness(*curPop, ADLTx, _sort, _mating_males);
        nbInd[MAL] += _nbIndividuals[MAL];      // get total number of males
        nbInd[FEM] += _nbIndividuals[FEM];      // get total number of females
        
        // sum of the fitnesses and store the fitness arrays
        sumFitness[i]   = _pSelection->getSumFitness(MAL) + _pSelection->getSumFitness(FEM);
        totFitness      += sumFitness[i];
        if(_nbIndividuals[MAL]) fitnessArrays[MAL][i] = _pSelection->get_FitnessObject(MAL, true);
        if(_nbIndividuals[FEM]) fitnessArrays[FEM][i] = _pSelection->get_FitnessObject(FEM, true);
    }
    
    // compute the total number of offspring
    nbOff = (this->*setNbOffspring_func_ptr)(nbInd[MAL], nbInd[FEM], _popPtr->get_total_carrying_capacity());
    
    // create the offspring for each patch separately
    for(i=0, curPop = _popPtr->get_vFullPatch().begin(); curPop!=endPop; ++curPop, ++i) {
        if(!fitnessArrays[MAL][i]) continue;  // female or male missing -> no offspring
        
        // (re)set the patch parameters
        _nbIndividuals[MAL] = (*curPop)->size(MAL, ADLTx);                // reset number of males
        _nbIndividuals[FEM] = (*curPop)->size(FEM, ADLTx);                // reset number of females
        if(_nbIndividuals[MAL]) _pSelection->set_FitnessObject(MAL, fitnessArrays[MAL][i]);    // set the fitness arrays of the current patch
        if(_nbIndividuals[FEM]) _pSelection->set_FitnessObject(FEM, fitnessArrays[FEM][i]);    // set the fitness arrays of the current patch
        
        // create the mating pairs
        preMating(*curPop);
        
        // get the number of sons/daughters to produce of this patch
        if(totFitness) nbBaby = _pSelection->get_SoftHardSelection(_nbIndividuals[MAL]+_nbIndividuals[FEM], (unsigned int)(sumFitness[i]/totFitness*nbOff));
        else nbBaby = _pSelection->get_SoftHardSelection(_nbIndividuals[MAL]+_nbIndividuals[FEM], (unsigned int)((_nbIndividuals[MAL]+_nbIndividuals[FEM])/nbOff));
        (this->*setSexRatio_func_ptr)(nbBaby, nbSons, nbDaughters, _nbIndividuals[MAL], _nbIndividuals[FEM]);
        
        createOffspring(*curPop, nbDaughters, nbSons);
    }//end_for_nbPatch
    
    // delete the fitness arrays
    if(sumFitness)      delete[] sumFitness;
    for(i=0; i<2; ++i){
        if(fitnessArrays[i]) delete[] fitnessArrays[i];
    }
}

// ----------------------------------------------------------------------------------------
// LCE_Breed::breed_metapop2hard
// ----------------------------------------------------------------------------------------
/** This breeding function allows to simualte any degree of selection between metapop and
 * hard selection. However the function is not optimal for just hard selection.
 * First the fitness of all individuals is computed. The number of total offspring (of the metapop)
 * is distributed among the patches according to the mean fitness.
 * It does not make sence to use this function for a neutral case.
 * nbOff = (1-f)*(Wp/Wm)*Km  + f*Wp*Kp
 *   f=0: metapop selection
 *   f=1: hard selection
 */
void LCE_Breed::breed_metapop2hard ()
{
#ifdef _DEBUG
    message("(metapop-hard) ... ");
#endif
    
    assert(_pSelection);
    
    unsigned int nbBaby, i, nbPatches = _popPtr->get_nbFullPatch();
    unsigned int nbDaughters, nbSons, nbOff;
    
    // array to buffer the fitness arrays for each patch (female and male)
    double totFitness = 0;                              // total fitness
    double* sumFitness = new double[nbPatches];        // mean fitness of each patch
    unsigned nbInd[2];
    TPatchFitness** fitnessArrays[2];
    for(i=0; i<2; ++i){
        fitnessArrays[i] = new TPatchFitness*[nbPatches];
        nbInd[i] = 0;
    }
    
    // for each patch compute the fitness of the individuals
    vector<Patch*>::iterator curPop = _popPtr->get_vFullPatch().begin();
    vector<Patch*>::iterator endPop = _popPtr->get_vFullPatch().end();
    for(i=0; curPop!=endPop; ++curPop, ++i) {
        // check if the mating requirements are met
        if(!(this->*isMatingPossible_func_ptr)(*curPop)){
            fitnessArrays[MAL][i] = NULL;
            fitnessArrays[FEM][i] = NULL;
            continue;
        }
        
        // compute and store the fitness of the current patch (0: male, 1: female)
        _pSelection->set_fitness(*curPop, ADLTx, _sort, _mating_males);
        nbInd[MAL] += _nbIndividuals[MAL];      // get total number of males
        nbInd[FEM] += _nbIndividuals[FEM];      // get total number of females
        
        // sum of the fitnesses and store the fitness arrays
        sumFitness[i]   = _pSelection->getSumFitness(MAL) + _pSelection->getSumFitness(FEM);
        totFitness      += sumFitness[i];
        if(_nbIndividuals[MAL]) fitnessArrays[MAL][i] = _pSelection->get_FitnessObject(MAL, true);
        if(_nbIndividuals[FEM]) fitnessArrays[FEM][i] = _pSelection->get_FitnessObject(FEM, true);
    }
    
    // compute the total number of offspring
    nbOff = (this->*setNbOffspring_func_ptr)(nbInd[MAL], nbInd[FEM], _popPtr->get_total_carrying_capacity());
    
    // create the offspring for each patch separately
    for(i=0, curPop = _popPtr->get_vFullPatch().begin(); curPop!=endPop; ++curPop, ++i) {
        if(!fitnessArrays[MAL][i]) continue;  // female or male missing -> no offspring
        
        // (re)set the patch parameters
        _nbIndividuals[MAL] = (*curPop)->size(MAL, ADLTx);                // reset number of males
        _nbIndividuals[FEM] = (*curPop)->size(FEM, ADLTx);                // reset number of females
        if(_nbIndividuals[MAL]) _pSelection->set_FitnessObject(MAL, fitnessArrays[MAL][i]);    // set the fitness arrays of the current patch
        if(_nbIndividuals[FEM]) _pSelection->set_FitnessObject(FEM, fitnessArrays[FEM][i]);    // set the fitness arrays of the current patch
        
        // create the mating pairs
        preMating(*curPop);
        
        // get the number of sons/daughters to produce of this patch
        if(totFitness) nbBaby = _pSelection->get_SoftHardSelection((unsigned int)(sumFitness[i]/totFitness*nbOff), (unsigned int)(_pSelection->getMeanFitness()*(_nbIndividuals[MAL]+_nbIndividuals[FEM])));
        else nbBaby = _pSelection->get_SoftHardSelection((unsigned int)((_nbIndividuals[MAL]+_nbIndividuals[FEM])/nbOff), (unsigned int)(_pSelection->getMeanFitness()*(_nbIndividuals[MAL]+_nbIndividuals[FEM])));
        (this->*setSexRatio_func_ptr)(nbBaby, nbSons, nbDaughters, _nbIndividuals[MAL], _nbIndividuals[FEM]);
        
        createOffspring(*curPop, nbDaughters, nbSons);
    }//end_for_nbPatch
    
    // delete the fitness arrays
    if(sumFitness)      delete[] sumFitness;
    for(i=0; i<2; ++i){
        if(fitnessArrays[i]) delete[] fitnessArrays[i];
    }
}

// ----------------------------------------------------------------------------------------
// breed_soft2hard
// ----------------------------------------------------------------------------------------
/** This breeding function allows to simulate any degree of selection between soft and hard selection.
 * nbOff = (1-f)*Kp  + f*Wp*Kp
 *   f=0: soft selection
 *   f=1: hard selection
 */
void LCE_Breed::breed_soft2hard ()
{
#ifdef _DEBUG
    message("(soft-hard) ... ");
#endif
    
    unsigned int nbBaby, nbSons, nbDaughters, Kp;
    
    // for each patch
    vector<Patch*>::iterator curPop = _popPtr->get_vFullPatch().begin();
    vector<Patch*>::iterator endPop = _popPtr->get_vFullPatch().end();
    for(; curPop!=endPop; ++curPop) {
        // check if the mating requirements are met
        if(!(this->*isMatingPossible_func_ptr)(*curPop)) continue;
        
        // create the mating pairs
        preMating(*curPop);
        
        // compute the fitness of the current patch (0: male, 1: female)
        if(_pSelection) _pSelection->set_fitness(*curPop, ADLTx, _sort, _mating_males);
        
        // get the number of sons/daughters to produce of this patch
        Kp = (this->*setNbOffspring_func_ptr)(_nbIndividuals[MAL], _nbIndividuals[FEM], (*curPop)->get_K());
        nbBaby = _pSelection->get_SoftHardSelection(Kp, (unsigned int)(_pSelection->getMeanFitness()*Kp));
        (this->*setSexRatio_func_ptr)(nbBaby, nbSons, nbDaughters, _nbIndividuals[MAL], _nbIndividuals[FEM]);
        
        createOffspring(*curPop, nbDaughters, nbSons);
    }//end_for_nbPatch
}

// ----------------------------------------------------------------------------------------
// LCE_Breed
// ----------------------------------------------------------------------------------------
/** create the daugthers and sons */
void
LCE_Breed::createOffspring(Patch* cur_patch, unsigned int nbDaughters, unsigned int nbSons){
    Individual *MotherPtr, *FatherPtr;
    unsigned int index;
    
    // daughters: mate randomly a female and a male following their fitness
    while(nbDaughters) {
        MotherPtr = getMotherPtr(cur_patch, index);        // get a random male and female (according to their fitness)
        FatherPtr = getFatherPtr(cur_patch, index);
        _popPtr->makeOffsprg(MotherPtr, FatherPtr, FEM, cur_patch); // create the baby
        --nbDaughters;
    }//end_for_nbDaughters
    
    // sons: mate randomly a female and a male following their fitness
    while(nbSons) {
        MotherPtr = getMotherPtr(cur_patch, index);        // get a random male and female (according to their fitness)
        FatherPtr = getFatherPtr(cur_patch, index);
        _popPtr->makeOffsprg(MotherPtr, FatherPtr, MAL, cur_patch);   // create the baby
        --nbSons;
    }//end_for_nbSons
}

// ----------------------------------------------------------------------------------------
// preMating
// ----------------------------------------------------------------------------------------
/** this funciton is called before mating starts in a patch	*/
void
LCE_Breed::preMating(Patch* cur_patch){
    if(  getFather_func_ptr == &LCE_Breed::Monogyny_S_MatingFunc
       || getFather_func_ptr == &LCE_Breed::Monogyny_MatingFunc) create_mating_pairs(cur_patch);
}

// ----------------------------------------------------------------------------------------
// LCE_Breed
// ----------------------------------------------------------------------------------------
/** Create mating pairs for the monogamy mating system.
 * The number of paris corresponds to the number of females.
 * if(nbFemales > nbMales) males may mate with several females
 * if(nbFemales < nbMales) not all males may mate
 */
void
LCE_Breed::create_mating_pairs(Patch* cur_patch)
{
    unsigned int i, pos, nbMales;
    
    unsigned int nbPairs = cur_patch->size(FEM, ADLTx);    //number of females/ number of mating pairs
    if(nbPairs > _aMatingPairs_size){
        _aMatingPairs_size = nbPairs;
        if(_aMatingPairs[MAL]) delete[] _aMatingPairs[MAL];
        if(_aMatingPairs[FEM]) delete[] _aMatingPairs[FEM];
        _aMatingPairs[MAL] = new Individual*[_aMatingPairs_size];
        _aMatingPairs[FEM] = new Individual*[_aMatingPairs_size];
    }
    
    nbMales = cur_patch->size(MAL, ADLTx);    // number of males
    
    // while we do not have made all the pairs
    while(nbPairs>0){
        
        // put the males into a vector
        vector<Individual*> vecMales;
        for(i = 0; i<nbMales; ++i){
            vecMales.push_back(cur_patch->get(MAL, ADLTx, i));
        }
        
        // create the pairs (choose randomly a male for each female)
        for(i = 0; i<nbMales && nbPairs>0; ++i){
            --nbPairs;
            _aMatingPairs[FEM][nbPairs] = cur_patch->get(FEM, ADLTx, nbPairs); // get the female
            pos = get_pop_ptr()->rand().Uniform((unsigned int)vecMales.size());// select randomly a male
            _aMatingPairs[MAL][nbPairs] = vecMales[pos];                       // put this male to the pair
            vecMales.erase(vecMales.begin()+pos);                              // remove the male from the vector
        }
    }
}


// ----------------------------------------------------------------------------------------
// LCE_Breed
// ----------------------------------------------------------------------------------------
/** This function removes all indivduals where the fitness factor results in a fitness of zero,
 * is lethal.
 * TSelection knows all quanti traits, so we have to pass by TSelection
 */
void
LCE_Breed::remove_inds_zero_fitnessFactor(age_idx AGE)
{
    if(!_pSelection) return;
    assert(_fitness_factor_zero_isLethal);
    
    vector<Patch*>::iterator curPop = _popPtr->get_vFullPatch().begin();
    vector<Patch*>::iterator endPop = _popPtr->get_vFullPatch().end();
    Individual* ind;
    int i, nbInds;
    vector<int> vTraits = _pSelection->get_vTraits();
    vector<int>::iterator curTrait, endTrait = vTraits.end();
    for(; curPop!=endPop; ++curPop) {                                      // for each populated patch
        // check all females (starting from the back)
        nbInds = (int)(*curPop)->size(FEM, AGE);
        for (i = nbInds-1; i >= 0; --i) {                                     // for each female
            ind = (*curPop)->get(FEM, AGE, i);
            for(curTrait = vTraits.begin(); curTrait!=endTrait; ++curTrait){   // for each quanti trait
                if(!ind->set_getFitnessFactor(*curTrait)){                       // if W=0: remove individual and continue;
                    (*curPop)->recycle(FEM, AGE, i);
                    break;
                }
            }
        }
        // check all males
        nbInds = (int)(*curPop)->size(MAL, AGE);
        for (i = nbInds-1; i >= 0; --i) {                                     // for each male
            ind = (*curPop)->get(MAL, AGE, i);
            for(curTrait = vTraits.begin(); curTrait!=endTrait; ++curTrait){   // for each quanti trait
                if(!ind->set_getFitnessFactor(*curTrait)){                       // if W=0: remove individual and continue;
                    (*curPop)->recycle(MAL, AGE, i);
                    break;
                }
            }
        }
    }
}

// ----------------------------------------------------------------------------------------
// LCE_Breed
// ----------------------------------------------------------------------------------------
/** This function to set the sexes according to the phenotype of the first trait.
 * female <= threshold < male
 */
void
LCE_Breed::reset_sex_after_phenotype(age_idx AGE)
{
    Individual* ind;
    int i, sizeF, sizeM;
    assert(get_pop_ptr()->individual_container_ok());
    
    vector<Patch*>::iterator curPop = _popPtr->get_vFullPatch().begin();
    vector<Patch*>::iterator endPop = _popPtr->get_vFullPatch().end();
    for(; curPop!=endPop; ++curPop) {
        sizeF = (int)(*curPop)->size(FEM, AGE);
        sizeM = (int)(*curPop)->size(MAL, AGE);
        
        // check all females
        for (i = sizeF-1; i >= 0; --i) {
            ind = (*curPop)->get(FEM, AGE, i);
            _pSelection->get_selTrait(0)->set_phenotype(ind);  // set the phenotype
            if(_pSelection->get_selTrait(0)->get_phenotype() > _threshold){    // above threshold it should be a male
                ind->switch_sex(AGE,i);
            }
        }
        
        // check all males
        for (i = sizeM-1; i >= 0; --i) {      // note that here sizeM is not anymore cur_patch->size(MAL, ADLTx)
            ind = (*curPop)->get(MAL, AGE, i);
            _pSelection->get_selTrait(0)->set_phenotype(ind);  // set the phenotype
            if(_pSelection->get_selTrait(0)->get_phenotype() <= _threshold){    // above threshold it should be a male
                ind->switch_sex(AGE,i);
            }
        }
    }
    assert(get_pop_ptr()->individual_container_ok());
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------------------------
// LCE_Breed::execute
// ----------------------------------------------------------------------------------------
void LCE_Breed_coal::execute () {
#ifdef _DEBUG
    message("  LCE_Breed_Coal ");
#endif
    
    breed_coal();
    
#ifdef _DEBUG
    message("done! (Patches: %i; Individuals(F/M): [%i/%i; %i/%i])\n",
            _popPtr->get_nbFullPatch(),
            _popPtr->size(FEM, OFFSx),  _popPtr->size(MAL, OFFSx),
            _popPtr->size(FEM, ADLTx),  _popPtr->size(MAL, ADLTx));
#endif
}

// ----------------------------------------------------------------------------------------
// LCE_Breed::init
// ----------------------------------------------------------------------------------------
bool LCE_Breed_coal::init(Metapop* popPtr)
{
    LCE::init(popPtr);
    
    _mating_system      = (int)this->get_parameter_value("mating_system");
    _mean_fecundity     = this->get_parameter_value("mean_fecundity");
    _growth_rate        = this->get_parameter_value("growth_rate");
    
    set_nbOffspring_model();
    
    return true;
}

// ----------------------------------------------------------------------------------------
// LCE_Breed::breed_coal
// ----------------------------------------------------------------------------------------
bool LCE_Breed_coal::breed_coal()
{
#ifdef _DEBUG
    message("(breed_coal) ... ");
#endif
    
    // for each patch
    unsigned int size;
    vector<Patch*>::iterator curPop, endPop = _popPtr->get_vFullPatch().end();
    for(curPop = _popPtr->get_vFullPatch().begin(); curPop!=endPop;) {
        size = (*curPop)->size(FEM, ADLTx);
        if(size) size = (this->*setNbOffspring_func_ptr)(0, size, (*curPop)->get_K());
        (*curPop)->set_size(FEM, OFFSx, size);
        (*curPop)->set_size(FEM, ADLTx, 0);      // has to be set to 0 since dispersal is adding (and not setting) the number of inds
        if(size) ++curPop;
        else _popPtr->new_emptyPatch(curPop, endPop);   // curPop and endPop are adjusted
    }//end_for_nbPatch
    return true;
}

// ----------------------------------------------------------------------------------------
// LCE_Breed::isMatingPossible
// ----------------------------------------------------------------------------------------
/** function to set the number of females and males and checks if the conditions are met for mating */
bool LCE_Breed::isMatingPossible_1_sex(Patch* cur_patch){
    _nbIndividuals[FEM] = cur_patch->size(FEM, ADLTx);
    assert(!cur_patch->size(MAL, ADLTx));
    _nbIndividuals[MAL] = 0;            // no males are present
    if(cur_patch->size(OFFSx)) cur_patch->flush(OFFSx);
    return (_nbIndividuals[FEM]>0);
}
bool LCE_Breed::isMatingPossible_2_sex(Patch* cur_patch){
    _nbIndividuals[FEM] = cur_patch->size(FEM, ADLTx);
    _nbIndividuals[MAL] = cur_patch->size(MAL, ADLTx);
    if(cur_patch->size(OFFSx)) cur_patch->flush(OFFSx);
    return (_nbIndividuals[FEM] && _nbIndividuals[MAL]);
}

bool LCE_Breed::isMatingPossible_sex_allocation_fix(Patch* cur_patch){
    unsigned int nbInds = (unsigned int)cur_patch->get_all_inds(FEM, ADLTx).size();
    _nbIndividuals[FEM] = my_round(_fem_sex_allocation_prop*nbInds);
    _nbIndividuals[MAL] = nbInds-_nbIndividuals[FEM];
    if(cur_patch->size(OFFSx)) cur_patch->flush(OFFSx);
    return (_nbIndividuals[FEM] && _nbIndividuals[MAL]);
}

bool LCE_Breed::isMatingPossible_sex_allocation_G(Patch* cur_patch){
    double nbFem = 0;
    vector<Individual*> inds = cur_patch->get_all_inds(FEM, ADLTx);
    vector<Individual*>::iterator cur=inds.begin(), end=inds.end();
    for(;cur!=end; ++cur){
        nbFem += (*cur)->getTraitGenotype(_fem_sex_allocation_all_trait);
    }
    _nbIndividuals[FEM] = my_round(nbFem);
    _nbIndividuals[MAL] = (unsigned int)inds.size()-_nbIndividuals[FEM];
    if(cur_patch->size(OFFSx)) cur_patch->flush(OFFSx);
    return (_nbIndividuals[FEM] && _nbIndividuals[MAL]);
}

bool LCE_Breed::isMatingPossible_sex_allocation_Z(Patch* cur_patch){
    double nbFem = 0;
    vector<Individual*> inds = cur_patch->get_all_inds(FEM, ADLTx);
    vector<Individual*>::iterator cur=inds.begin(), end=inds.end();
    for(;cur!=end; ++cur){
        _pSelection->set_phenotype(*cur, _fem_sex_allocation_quanti_trait);
        nbFem += (*cur)->getTraitValue(_fem_sex_allocation_all_trait);
    }
    _nbIndividuals[FEM] = my_round(nbFem);
    _nbIndividuals[MAL] = (unsigned int)inds.size()-_nbIndividuals[FEM];
    if(cur_patch->size(OFFSx)) cur_patch->flush(OFFSx);
    return (_nbIndividuals[FEM] && _nbIndividuals[MAL]);
}

bool LCE_Breed::isMatingPossible_sex_allocation_equation(Patch* cur_patch){
    double nbFem = 0;
    vector<Individual*> inds = cur_patch->get_all_inds(FEM, ADLTx);
    vector<Individual*>::iterator cur=inds.begin(), end=inds.end();
    for(;cur!=end; ++cur){
        _pSelection->set_phenotype(*cur);
        nbFem += _fem_sex_allocation->getValue(*cur);
    }
    _nbIndividuals[FEM] = nbFem;
    _nbIndividuals[MAL] = (unsigned int)inds.size()-_nbIndividuals[FEM];
    if(cur_patch->size(OFFSx)) cur_patch->flush(OFFSx);
    return (_nbIndividuals[FEM] && _nbIndividuals[MAL]);
}

// ----------------------------------------------------------------------------------------
// LCE_Breed::setSexRatio
// ----------------------------------------------------------------------------------------
/** function to compute the number of sons/daugthers depending on the total number of childs
 * nbSons and nbDaughters are changed
 * nbBaby is the total number of new offspring
 * nbMAL and nbFEM are the number of females and males currently in the patch
 *     (nbMAL and nbFEM are used to compute the current sex ratio)
 */
void LCE_Breed::setSexRatio_Selfing(const unsigned int& nbBaby, unsigned int& nbSons, unsigned int& nbDaugthers, const unsigned int& nbMAL, const unsigned int& nbFEM){
    nbSons = 0;
    nbDaugthers = nbBaby;
}
void LCE_Breed::setSexRatio_NoSelfing(const unsigned int& nbBaby, unsigned int& nbSons, unsigned int& nbDaugthers, const unsigned int& nbMAL, const unsigned int& nbFEM){
    nbSons = my_round(get_pop_ptr()->rand().Binomial(_sex_ratio, nbBaby));
    nbDaugthers = nbBaby - nbSons;
}
void LCE_Breed::setSexRatio_KeepSexRatio(const unsigned int& nbBaby, unsigned int& nbSons, unsigned int& nbDaugthers, const unsigned int& nbMAL, const unsigned int& nbFEM){
    nbSons = my_round(nbBaby*((double)nbMAL)/(nbMAL+nbFEM));
    nbDaugthers = nbBaby - nbSons;
}

// ----------------------------------------------------------------------------------------
// LCE_Breed::getRandomSex
// ----------------------------------------------------------------------------------------
sex_t LCE_Breed::getRandomSex_Selfing(const unsigned int& nbMAL, const unsigned int& nbFEM){
    return FEM;
}
sex_t LCE_Breed::getRandomSex_NoSelfing(const unsigned int& nbMAL, const unsigned int& nbFEM){
    if(get_pop_ptr()->rand().Uniform()<_sex_ratio) return MAL;
    else                                  return FEM;
}
sex_t LCE_Breed::getRandomSex_KeepSexRatio(const unsigned int& nbMAL, const unsigned int& nbFEM){
    if(get_pop_ptr()->rand().Uniform()<(nbMAL/(double)(nbMAL+nbFEM))) return MAL;
    else                                                     return FEM;
}



/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

