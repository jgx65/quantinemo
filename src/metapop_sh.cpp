/** @file metapop_sh.cpp
 *
 *   Copyright (C) 2006 Frederic Guillaume    <guillaum@zoology.ubc.ca>
 *   Copyright (C) 2008 Samuel Neuenschwander <samuel.neuenschwander@unil.ch>
 *   Copyright (C) 2018 Frederic Michaud <frederic.michaud@unil.ch>
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


#include "metapop_sh.h"
#include "tmetapop.h"
#include "stathandler.cpp"
#include "tsimulation.h"


/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                             ******** StatHandler ********

// ----------------------------------------------------------------------------------------
// setStatsRecorders
// ----------------------------------------------------------------------------------------
bool MetapopSH::setStatRecorders(const string& t)
{
    string token = t;
    
    // get the end
    string::size_type pos = token.find('_');
    string end;                                        // "end" is empty
    if(pos != string::npos){                           // there is an end
        end = token.substr(pos+1);                       // eg. "pair"
        token.erase(pos+1);                              // eg. "adlt.fst_"
    }
    
    // patch extinction
    if(token == "ext.rate") return add("Extinction rate","ext.rate",FLAT,ALL,0,&MetapopSH::getObsrvdExtinctionRate);
    
    // Fecundity
    if(token == "fem.meanFec") return add("Mean fecundity of females","fem.meanFec",FLAT,ADULTS,1,0,&MetapopSH::getReproductiveMean);
    if(token == "fem.varFec")  return add("Variance of the fecundity of females","fem.varFec",FLAT,ADULTS,1,0,&MetapopSH::getReproductiveVar);
    if(token == "mal.meanFec") return add("Mean fecundity of males","mal.meanFec",FLAT,ADULTS,0,0,&MetapopSH::getReproductiveMean);
    if(token == "mal.varFec")  return add("Variance of the fecundity of males","mal.varFec",FLAT,ADULTS,0,0,&MetapopSH::getReproductiveVar);
    
    if(token == "fecundity"){
        add("Mean realized female fecundity","fem.real.fec",FLAT,ADULTS,1,0,&MetapopSH::getReproductiveMean);
        add("Variance of the mean realized female fecundity","fem.var.fec",FLAT,ADULTS,1,0,&MetapopSH::getReproductiveVar);
        add("Mean realized male fecundity  ","mal.real.fec",FLAT,ADULTS,0,0,&MetapopSH::getReproductiveMean);
        add("Variance of the mean realized male fecundity","mal.var.fec",FLAT,ADULTS,0,0,&MetapopSH::getReproductiveVar);
        return true;
    }
    
    // migration
    if(token == "emigrants")  return add("Mean number of emigrants/patch",token,FLAT,ADULTS,0,&MetapopSH::getMeanEmigrantPerPatch);
    if(token == "immigrants") return add("Mean number of immigrants/patch",token,FLAT,ADULTS,0,&MetapopSH::getMeanImmigrantPerPatch);
    if(token == "residents")  return add("Mean number of residents/patch",token,FLAT,ADULTS,0,&MetapopSH::getMeanResidantPerPatch);
    if(token == "immigrate")  return add("Mean effective immigration rate",token,FLAT,ADULTS,0,&MetapopSH::getMeanMigrantRatio);
    if(token == "colonisers") return add("Mean number of colonizers per extinct patch",token,FLAT,ADULTS,0,&MetapopSH::getMeanKolonisersPerPatch);
    if(token == "colon.rate") return add("Mean effective colonization rate of extinct patches",token,FLAT,ADULTS,0,&MetapopSH::getMeanKolonisersProportion);
    
    if(token == "migration"){
        add("Mean number of emigrants/patch","emigrants",FLAT,ADULTS,0,&MetapopSH::getMeanEmigrantPerPatch);
        add("Mean number of immigrants/patch","immigrants",FLAT,ADULTS,0,&MetapopSH::getMeanImmigrantPerPatch);
        add("Mean number of residents/patch","residents",FLAT,ADULTS,0,&MetapopSH::getMeanResidantPerPatch);
        add("Mean effective immigration rate","immigrate",FLAT,ADULTS,0,&MetapopSH::getMeanMigrantRatio);
        add("Mean number of colonizers per extinct patch","colonisers",FLAT,ADULTS,0,&MetapopSH::getMeanKolonisersPerPatch);
        add("Mean effective colonization rate of extinct patches","colon.rate",FLAT,ADULTS,0,&MetapopSH::getMeanKolonisersProportion);
        return true;
    }
    
    // fitness
    if(token=="VwW")     return add("Fitness variance of adults within patches ",token,FLAT,ADULTS,0,&MetapopSH::getVwW);
    if(token=="VwB")     return add("Fitness variance of adults between patches ",token,FLAT,ADULTS,0,&MetapopSH::getVwB);
    if(token=="meanW_")  return add_perPatch(end,"Fitness mean of adults","meanW",FLAT,ADULTS,0,0,0,&MetapopSH::getMeanW);
    if(token=="varW_")   return add_perPatch(end,"Fitness variance of adults","varW",FLAT,ADULTS,0,0,0,&MetapopSH::getVarW);
    
    if(token=="fitness") {
        add("Fitness variance of adults within patches ","VwW",FLAT,ADULTS,0,&MetapopSH::getVwW);
        add("Fitness variance of adults between patches ","VwB",FLAT,ADULTS,0,&MetapopSH::getVwB);
        return true;
    }
    
    // get the age
    pos = token.find('.');                        // find the first '.'
    if(pos == string::npos) return false;
    string ageToken = token.substr(0, pos);       // eg. "adlt"
    token.erase(0, pos+1);                        // eg. "fst_"
    age_t AGE;
    string ageStr;
    if(ageToken == "adlt")     {AGE = ADULTS;  ageStr = "adult";}
    else if(ageToken == "off") {AGE = OFFSPRG; ageStr = "offsprg";}
    else return false;
    
    
    if(set_stat_kinship(token, ageToken, end, AGE, ageStr)) return true;     // adults and offspring
    if(set_stat_demography(token, ageToken, end, AGE, ageStr)) return true;  // adults and offspring
    
    return false;
}

// ----------------------------------------------------------------------------------------
// set_stat_demography
// ----------------------------------------------------------------------------------------
/** fstat stat options */
bool MetapopSH::set_stat_demography(const string& t, const string& token,
                                    const string& end, const age_t& AGE, const string& ageStr)
{
    // Demography
    if(t == "nbInd")    return add("Total number of "+ageStr,token+".nbInd",FLAT,AGE,0,0,0,0,&MetapopSH::getNbInd);
    if(t == "nbFem")    return add("Total number of female "+ageStr,token+".nbFem",FLAT,AGE,0,0,0,0,&MetapopSH::getNbFem);
    if(t == "nbMal")    return add("Total number of male "+ageStr,token+".nbMal",FLAT,AGE,0,0,0,0,&MetapopSH::getNbMal);
    if(t == "meanInd")  return add("Mean number of "+ageStr,token+".meanInd",FLAT,AGE,0,0,0,0,&MetapopSH::getMeanNbInd);
    if(t == "meanFem")  return add("Mean number of female "+ageStr,token+".meanFem",FLAT,AGE,0,0,0,0,&MetapopSH::getMeanNbFem);
    if(t == "meanMal")  return add("Mean number of male "+ageStr,token+".meanMal",FLAT,AGE,0,0,0,0,&MetapopSH::getMeanNbMal);
    if(t == "sexRatio") return add("Sex ratio of "+ageStr,token+".sexRatio",FLAT,AGE,0,0,0,0,&MetapopSH::getSexRatio);
    if(t == "nbPops")   return add("Number of inhabited patches by "+ageStr,token+".nbPops",FLAT,AGE,0,0,0,0,&MetapopSH::getNbPops);
    if(t == "nbInd_")   return add_perPatch(end,"Number of individuals  "+ageStr,token+".nbInd",FLAT,AGE,0,0,0,0,0,0,&MetapopSH::getNbInd);
    if(t == "nbFem_")   return add_perPatch(end,"Number of female "+ageStr,token+".nbFem",FLAT,AGE,0,0,0,0,0,0,&MetapopSH::getNbFem);
    if(t == "nbMal_")   return add_perPatch(end,"Number of male "+ageStr,token+".nbMal",FLAT,AGE,0,0,0,0,0,0,&MetapopSH::getNbMal);
    
    
    if(t == "demo"){
        add("Total number of "+ageStr,token+".nbInd",FLAT,AGE,0,0,0,0,&MetapopSH::getNbInd);
        add("Total number of female "+ageStr,token+".nbFem",FLAT,AGE,0,0,0,0,&MetapopSH::getNbFem);
        add("Total number of male "+ageStr,token+".nbMal",FLAT,AGE,0,0,0,0,&MetapopSH::getNbMal);
        add("Mean number of "+ageStr,token+".meanInd",FLAT,AGE,0,0,0,0,&MetapopSH::getMeanNbInd);
        add("Mean number of female "+ageStr,token+".meanFem",FLAT,AGE,0,0,0,0,&MetapopSH::getMeanNbFem);
        add("Mean number of male "+ageStr,token+".meanMal",FLAT,AGE,0,0,0,0,&MetapopSH::getMeanNbMal);
        add("Sex ratio of "+ageStr,token+".sexRatio",FLAT,AGE,0,0,0,0,&MetapopSH::getSexRatio);
        add("Number of inhabited patches by "+ageStr,token+".nbPops",FLAT,AGE,0,0,0,0,&MetapopSH::getNbPops);
        return true;
    }
    
    // Demography of all pops and patches
    if(t == "nbIndTot")    return add("Total number (not just sampled) of "+ageStr,token+".nbIndTot",FLAT,AGE,0,0,0,0,&MetapopSH::getNbIndTot);
    if(t == "nbFemTot")    return add("Total number (not just sampled) of female "+ageStr,token+".nbFemTot",FLAT,AGE,0,0,0,0,&MetapopSH::getNbFemTot);
    if(t == "nbMalTot")    return add("Total number (not just sampled) of male "+ageStr,token+".nbMalTot",FLAT,AGE,0,0,0,0,&MetapopSH::getNbMalTot);
    if(t == "meanIndTot")  return add("Mean number (not just sampled) of "+ageStr,token+".meanIndTot",FLAT,AGE,0,0,0,0,&MetapopSH::getMeanNbIndTot);
    if(t == "meanFemTot")  return add("Mean number (not just sampled) of female "+ageStr,token+".meanFemTot",FLAT,AGE,0,0,0,0,&MetapopSH::getMeanNbFemTot);
    if(t == "meanMalTot")  return add("Mean number (not just sampled) of male "+ageStr,token+".meanMalTot",FLAT,AGE,0,0,0,0,&MetapopSH::getMeanNbMalTot);
    if(t == "sexRatioTot") return add("Sex ratio (not just sampled) of "+ageStr,token+".sexRatioTot",FLAT,AGE,0,0,0,0,&MetapopSH::getSexRatioTot);
    if(t == "nbPopsTot")   return add("Number (not just sampled) of inhabited patches by "+ageStr,token+".nbPopsTot",FLAT,AGE,0,0,0,0,&MetapopSH::getNbPopsTot);
    if(t == "nbIndTot_")   return add_perPatch(end,"Number (not just sampled) of individuals  "+ageStr,token+".nbIndTot",FLAT,AGE,0,0,0,0,0,0,&MetapopSH::getNbIndTot);
    if(t == "nbFemTot_")   return add_perPatch(end,"Number (not just sampled) of female "+ageStr,token+".nbFemTot",FLAT,AGE,0,0,0,0,0,0,&MetapopSH::getNbFemTot);
    if(t == "nbMalTot_")   return add_perPatch(end,"Number (not just sampled) of male "+ageStr,token+".nbMalTot",FLAT,AGE,0,0,0,0,0,0,&MetapopSH::getNbMalTot);
    
    
    if(t == "demoTot"){
        add("Total number (not just sampled) of "+ageStr,token+".nbInd",FLAT,AGE,0,0,0,0,&MetapopSH::getNbIndTot);
        add("Total number (not just sampled) of female "+ageStr,token+".nbFem",FLAT,AGE,0,0,0,0,&MetapopSH::getNbFemTot);
        add("Total number (not just sampled) of male "+ageStr,token+".nbMal",FLAT,AGE,0,0,0,0,&MetapopSH::getNbMalTot);
        add("Mean number (not just sampled) of "+ageStr,token+".meanInd",FLAT,AGE,0,0,0,0,&MetapopSH::getMeanNbIndTot);
        add("Mean number (not just sampled) of female "+ageStr,token+".meanFem",FLAT,AGE,0,0,0,0,&MetapopSH::getMeanNbFemTot);
        add("Mean number (not just sampled) of male "+ageStr,token+".meanMal",FLAT,AGE,0,0,0,0,&MetapopSH::getMeanNbMalTot);
        add("Sex ratio (not just sampled) of "+ageStr,token+".sexRatio",FLAT,AGE,0,0,0,0,&MetapopSH::getSexRatioTot);
        add("Number (not just sampled) of inhabited patches by "+ageStr,token+".nbPops",FLAT,AGE,0,0,0,0,&MetapopSH::getNbPopsTot);
        return true;
    }
    return false;
}

// ----------------------------------------------------------------------------------------
// set_stat_kinship
// ----------------------------------------------------------------------------------------
/** kinship stat options */
bool MetapopSH:: set_stat_kinship(const string& t, const string& token,
                                  const string& end, const age_t& AGE, const string& ageStr)
{
    if(t == "fsib")  return add("Mean proportion of full-sib ("+ageStr+")",token+".fsib",FLAT,AGE,3,0,0,0,0,0,&MetapopSH::getSibProportion);
    if(t == "phsib") return add("Mean proportion of paternal half-sib ("+ageStr+")",token+".phsib",FLAT,AGE,2,0,0,0,0,0,&MetapopSH::getSibProportion);
    if(t == "mhsib") return add("Mean proportion of maternal half-sib ("+ageStr+")",token+".mhsib",FLAT,AGE,1,0,0,0,0,0, &MetapopSH::getSibProportion);
    if(t == "nsib")	 return add("Mean proportion of non-sib ("+ageStr+")",token+".nsib",FLAT,AGE,0,0,0,0,0,0,&MetapopSH::getSibProportion);
    if(t == "self")  return add("Mean proportion of selfed offspring ("+ageStr+")",token+".self",FLAT,AGE,4,0,0,0,0,0,&MetapopSH::getSibProportion);
    
    if(t == "kinship"){
        add("Mean proportion of full-sib ("+ageStr+")",token+".fsib",FLAT,AGE,3,0,0,0,0,0,&MetapopSH::getSibProportion);
        add("Mean proportion of paternal half-sib ("+ageStr+")",token+".phsib",FLAT,AGE,2,0,0,0,0,0,&MetapopSH::getSibProportion);
        add("Mean proportion of maternal half-sib ("+ageStr+")",token+".mhsib",FLAT,AGE,1,0,0,0,0,0, &MetapopSH::getSibProportion);
        add("Mean proportion of non-sib ("+ageStr+")",token+".nsib",FLAT,AGE,0,0,0,0,0,0,&MetapopSH::getSibProportion);
        add("Mean proportion of selfed offspring ("+ageStr+")",token+".self",FLAT,AGE,4,0,0,0,0,0,&MetapopSH::getSibProportion);
        return true;
    }
    return false;
}
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                         ******** Migrants analysis ********
// ----------------------------------------------------------------------------------------
// getMeanEmigrantPerPatch
// ----------------------------------------------------------------------------------------
/** compute the mean number of emigrants per patch */
double MetapopSH::getMeanEmigrantPerPatch()
{
    // before any migration
    if(get_current_generation() == 1 || _popPtr->getCoalescence()){
        meanEmigrant = my_NAN; return
        meanEmigrant;
    }
    
    meanEmigrant = 0;
    unsigned int nbPatch=0;
    vector<TPatch*>::iterator curPop, endPop = get_vSamplePatch().end();
    for(curPop=get_vSamplePatch().begin(); curPop != endPop; ++curPop){
        if(!(*curPop)->size(ADLTx)) continue;     // the patch can just be colonised
        meanEmigrant += (*curPop)->nbEmigrant;
        ++nbPatch;
    }
    if(nbPatch) meanEmigrant /= nbPatch;
    else        meanEmigrant = my_NAN;
    return meanEmigrant;
}
// ----------------------------------------------------------------------------------------
// getMeanImmigrantPerPatch
// ----------------------------------------------------------------------------------------
double MetapopSH::getMeanImmigrantPerPatch()
{
    // check if the table has already been computed
    if(already_computed(_computed[0])) return meanImmigrant;
    if(get_current_generation()==1){meanImmigrant=my_NAN; return meanImmigrant;}
    
    meanImmigrant = 0;
    unsigned int nbPatch=0;
    vector<TPatch*>::iterator curPop, endPop = get_vSamplePatch().end();
    for(curPop=get_vSamplePatch().begin(); curPop != endPop; ++curPop){
        if(!(*curPop)->size(OFFSx)) continue;     // the patch can just become extinct
        meanImmigrant += (*curPop)->nbImmigrant;
        ++nbPatch;
    }
    if(nbPatch) meanImmigrant /= nbPatch;
    else        meanImmigrant = my_NAN;
    return meanImmigrant;
}

// ----------------------------------------------------------------------------------------
// getMeanResidantPerPatch
// ----------------------------------------------------------------------------------------
double MetapopSH::getMeanResidantPerPatch()
{
    // check if the table has already been computed
    if(already_computed(_computed[1])) return meanResidant;
    if(get_current_generation()==1){meanResidant=my_NAN; return meanResidant;}
    
    meanResidant = 0;
    unsigned int nbPatch=0;
    vector<TPatch*>::iterator curPop = get_vSamplePatch().begin();
    vector<TPatch*>::iterator endPop = get_vSamplePatch().end();
    for(; curPop != endPop; ++curPop){
        if(!(*curPop)->size(ADLTx)) continue;     // the patch can just be colonized
        meanResidant += (*curPop)->size(ADLTx) - (*curPop)->nbEmigrant;
        ++nbPatch;
    }
    if(nbPatch) meanResidant /= nbPatch;
    else        meanResidant = my_NAN;
    return meanResidant;
}

// ----------------------------------------------------------------------------------------
// getMeanMigrantRatio
// ----------------------------------------------------------------------------------------
double MetapopSH::getMeanMigrantRatio()
{
    // before any migration
    if(get_current_generation() == 1 || _popPtr->getCoalescence()) return my_NAN;
    
    getMeanImmigrantPerPatch(); // set the variable meanImmigrant
    getMeanResidantPerPatch();  // set the variable meanResidant
    if(meanImmigrant==my_NAN || meanResidant==my_NAN) return my_NAN;
    return ((meanImmigrant+meanResidant) ? meanImmigrant/(meanImmigrant+meanResidant) : my_NAN);
}
// ----------------------------------------------------------------------------------------
// getMeanMigrantProportion
// ----------------------------------------------------------------------------------------
double MetapopSH::getMeanKolonisersProportion()
{
    if(get_current_generation() == 1 || _popPtr->getCoalescence()) return my_NAN;
    
    double mean = 0;
    unsigned int nbFullPatch=0;
    
    TPatch* curPatch;
    vector<TPatch*>::iterator curPop = get_vSamplePatch().begin();
    vector<TPatch*>::iterator endPop = get_vSamplePatch().end();
    for(; curPop != endPop; ++curPop){
        curPatch = *curPop;
        if(curPatch->nbKolonisers != my_NAN){
            mean += (double)curPatch->nbKolonisers / curPatch->size(ADLTx);
            ++nbFullPatch;
        }
    }
    
    return  nbFullPatch ? mean/nbFullPatch : my_NAN;
}

// ----------------------------------------------------------------------------------------
// getMeanKolonisersPerPatch
// ----------------------------------------------------------------------------------------
double MetapopSH::getMeanKolonisersPerPatch()
{
    if(get_current_generation() == 1 || _popPtr->getCoalescence()) return my_NAN;
    
    double meanKolonisers = 0;
    unsigned int nbFullPatch=0;
    vector<TPatch*>::iterator curPop = get_vSamplePatch().begin();
    vector<TPatch*>::iterator endPop = get_vSamplePatch().end();
    for(; curPop != endPop; ++curPop){
        if((*curPop)->nbKolonisers != my_NAN){
            meanKolonisers += (*curPop)->nbKolonisers;
            ++nbFullPatch;
        }
    }
    return nbFullPatch ? meanKolonisers/nbFullPatch : my_NAN;
}

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                     ****** Patch extinction analysis ********

// ----------------------------------------------------------------------------------------
// setObsrvdExtinctionRate
// ----------------------------------------------------------------------------------------
double MetapopSH::getObsrvdExtinctionRate ()
{
    return get_current_nbSamplePatch() / get_nbTotSamplePatch();
}


/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                   ****** Fecundity and matings analysis ********

// ----------------------------------------------------------------------------------------
// getMeanRealizedFecundity
// ----------------------------------------------------------------------------------------
void MetapopSH::setReproductiveStats(bool sex)  // 0: MAL, 1: FEM
{
    // check if the table has already been computed
    if(already_computed(_computed[24], sex)) return;
    
    vector<double> stat;
    vector<TPatch*>::iterator curPop = get_vSamplePatch().begin();
    vector<TPatch*>::iterator endPop = get_vSamplePatch().end();
    vector<TIndividual*>::iterator curInd, endInd;
    for(; curPop != endPop; ++curPop){                                          // for each patch
        curInd = (*curPop)->get_sampled_inds((sex_t)sex, ADLTx).begin();
        endInd = (*curPop)->get_sampled_inds((sex_t)sex, ADLTx).end();
        for(; curInd!=endInd; ++curInd){                                          // for each indvidual
            stat.push_back((*curInd)->getTotRealizedFecundity());
        }
    }
    
    _mean_reprod_success = ARRAY::mean(stat);
    _var_reprod_success  = ARRAY::var(stat, _mean_reprod_success);
}

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                             ****** kinship analysis ********
// ----------------------------------------------------------------------------------------
// setKinship
// ----------------------------------------------------------------------------------------
/** compute the mean kinship */
void MetapopSH::setKinship(const age_idx& AGE)
{
    // check if the table has already been computed
    if(already_computed(_computed[2], AGE)) return;
    
    unsigned int j, tot_size=0;
    // if it is the first generation stop here
    if(get_current_generation() == 1 || _popPtr->getCoalescence()){
        for(j = 0; j < 5; ++j) {
            _sib_prop[AGE][j] = my_NAN;
        }
        return;
    }
    
    //reset counters
    for(j = 0; j < 5; ++j){
        _sib_prop[AGE][j] = 0.0;
    }
    
    vector<TPatch*>::iterator curPop = get_vSamplePatch().begin();
    vector<TPatch*>::iterator endPop = get_vSamplePatch().end();
    vector<TIndividual*>::iterator curInd1, curInd2, endInd1, endInd2;
    for(; curPop != endPop; ++curPop){                                          // for each patch
        vector<TIndividual*>& curFem = (*curPop)->get_sampled_inds(FEM, AGE);
        vector<TIndividual*>& curMal = (*curPop)->get_sampled_inds(MAL, AGE);
        tot_size += (unsigned int)(curFem.size() + curMal.size());
        
        //female-female
        for(curInd1=curFem.begin(), endInd1=curFem.end(); curInd1!=endInd1; ++curInd1){ // for each indvidual
            for(curInd2=curInd1+1; curInd2!=endInd1; ++curInd2) {
                setKinClassCounter(*curInd1, *curInd2, AGE);
            }
            if((*curInd1)->getIsSelfed()) _sib_prop[AGE][4]++;     //selfed offspring counter:
        }
        
        //male-male
        for(curInd1=curMal.begin(), endInd1=curMal.end(); curInd1!=endInd1; ++curInd1){ // for each indvidual
            for(curInd2=curInd1+1; curInd2!=endInd1; ++curInd2) {
                setKinClassCounter(*curInd1, *curInd2, AGE);
            }
            if((*curInd1)->getIsSelfed()) _sib_prop[AGE][4]++;     //selfed offspring counter:
        }
        
        //female-male
        for(curInd1=curFem.begin(), endInd1=curFem.end(); curInd1!=endInd1; ++curInd1){ // for each female
            for(curInd2=curMal.begin(), endInd2=curMal.end(); curInd2!=endInd2; ++curInd2){ // for each male
                setKinClassCounter(*curInd1, *curInd2, AGE);
            }
        }
    } //end for i < patchNbr
    
    //total number of pairwise comparisons:
    double tot = _sib_prop[AGE][0] + _sib_prop[AGE][1] + _sib_prop[AGE][2] + _sib_prop[AGE][3];
    
    if(tot){
        for(j = 0 ; j < 4; ++j){
            _sib_prop[AGE][j] /= tot;
        }
    }
    
    if(tot_size)	_sib_prop[AGE][4] /= tot_size;
    else          _sib_prop[AGE][4] = my_NAN;
}

// ----------------------------------------------------------------------------------------
// setKinClassCounter
// ----------------------------------------------------------------------------------------
/** sets the kinship */
void MetapopSH::setKinClassCounter(TIndividual *I1, TIndividual *I2, const age_idx& AGE)
{
    if(I1->getMotherID() == I2->getMotherID()){
        if(I1->getFatherID() == I2->getFatherID()) _sib_prop[AGE][3]++;   // full sibs
        else 																	     _sib_prop[AGE][1]++;   // maternal half sibs
    }
    else{
        if(I1->getFatherID() == I2->getFatherID()) _sib_prop[AGE][2]++;   // paternal half sibs
        else                                       _sib_prop[AGE][0]++;   // non sibs
    }
}
// ----------------------------------------------------------------------------------------
// stats across all inds and pops
double MetapopSH::getNbIndTot        (const age_t& AGE){return (double)get_nbInds(AGE);}
double MetapopSH::getNbFemTot        (const age_t& AGE){return (double)get_nbInds(AGE,FEM);}
double MetapopSH::getNbMalTot        (const age_t& AGE){return (double)get_nbInds(AGE,MAL);}

double MetapopSH::getMeanNbIndTot    (const age_t& AGE){return get_meanNbInds(AGE);}
double MetapopSH::getMeanNbFemTot    (const age_t& AGE){return get_meanNbInds(AGE,FEM);}
double MetapopSH::getMeanNbMalTot    (const age_t& AGE){return get_meanNbInds(AGE,MAL);}

double MetapopSH::getNbIndTot        (unsigned int i, const age_t& AGE){return (double)get_vPatch(i)->size(AGE);}
double MetapopSH::getNbFemTot        (unsigned int i, const age_t& AGE){return (double)get_vPatch(i)->size(FEM, AGE);}
double MetapopSH::getNbMalTot        (unsigned int i, const age_t& AGE){return (double)get_vPatch(i)->size(MAL, AGE);}

double MetapopSH::getNbPopsTot       (const age_t& AGE){return (double)get_nbSamplePatch(AGE);}

double MetapopSH::getSexRatioTot     (const age_t& AGE) {
    unsigned int nbFem = get_nbInds(AGE,FEM);
    return nbFem ? (double) get_nbInds(AGE,MAL)/nbFem : my_NAN;
}

// stats across sampled inds and pops
double MetapopSH::getNbInd           (const age_t& AGE){return (double)get_nbSamples(AGE);}
double MetapopSH::getNbFem           (const age_t& AGE){return (double)get_nbSamples(AGE,FEM);}
double MetapopSH::getNbMal           (const age_t& AGE){return (double)get_nbSamples(AGE,MAL);}

double MetapopSH::getMeanNbInd       (const age_t& AGE){return get_meanNbSamples(AGE);}
double MetapopSH::getMeanNbFem       (const age_t& AGE){return get_meanNbSamples(AGE,FEM);}
double MetapopSH::getMeanNbMal       (const age_t& AGE){return get_meanNbSamples(AGE,MAL);}

double MetapopSH::getNbInd           (unsigned int i, const age_t& AGE){return (double)get_vPatch(i)->sampleSize(AGE);}
double MetapopSH::getNbFem           (unsigned int i, const age_t& AGE){return (double)get_vPatch(i)->sampleSize(FEM, AGE);}
double MetapopSH::getNbMal           (unsigned int i, const age_t& AGE){return (double)get_vPatch(i)->sampleSize(MAL, AGE);}

double MetapopSH::getNbPops          (const age_t& AGE){return (double)get_nbSamplePatch(AGE);}

double MetapopSH::getSexRatio        (const age_t& AGE) {
    unsigned int nbFem =  get_nbSamples(AGE,FEM);
    return nbFem ? (double) get_nbSamples(AGE,MAL)/nbFem : my_NAN;
}
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                             ******** Fitness ********

// ----------------------------------------------------------------------------------------
// setMeanAndVarFitness
// ----------------------------------------------------------------------------------------
void
MetapopSH::setMeanAndVar_W(){
    // check if the table has already been computed
    if(already_computed(_computed[3])) return;
    
    if(!_meanW) _meanW = new double[get_current_nbSamplePatch()];
    else if(get_current_nbSamplePatch()!=get_last_nbSamplePatch()){
        delete[] _meanW;
        _meanW = new double[get_current_nbSamplePatch()];
    }
    else ARRAY::reset_1D(_meanW, get_current_nbSamplePatch(), (double)my_NAN);
    
    if(!_varW)  _varW  = new double[get_current_nbSamplePatch()];
    else if(get_current_nbSamplePatch()!=get_last_nbSamplePatch()){
        delete[] _varW;
        _varW = new double[get_current_nbSamplePatch()];
    }
    else ARRAY::reset_1D(_varW, get_current_nbSamplePatch(), (double)my_NAN);
    
    vector<TPatch*>::iterator curPop = get_vSamplePatch().begin();
    vector<TPatch*>::iterator endPop = get_vSamplePatch().end();
    for(; curPop != endPop; ++curPop){                                          // for each patch
        setMeanAndVar_W_ofPatch(*curPop, (*curPop)->get_sampleID());
    }
}

// ----------------------------------------------------------------------------------------
// setMeanAndVarFitness
// ----------------------------------------------------------------------------------------
void
MetapopSH::setMeanAndVar_W_ofPatch(TPatch* crnt_patch, const unsigned int& i){
    unsigned int sizeF = crnt_patch->size(FEM, ADLTx),
    sizeM = crnt_patch->size(MAL, ADLTx),
    size  = sizeF + sizeM;
    
    // if the patch is empty or the fitness is not yet computed -> stop
    if(!((sizeF && crnt_patch->get(FEM, ADLTx, 0)->getFitness() != my_NAN)
         ||(sizeM && crnt_patch->get(MAL, ADLTx, 0)->getFitness() != my_NAN))){
        _meanW[i] = _varW[i] = my_NAN;
        return;
    }
    
    unsigned int f, m;
    double* array = new double[size];
    
    for(f = 0; f < sizeF; ++f) {
        array[f] = crnt_patch->get(FEM, ADLTx, f)->getFitness();
    }
    for(m = 0; m < sizeM; ++m, ++f) {
        array[f] = crnt_patch->get(MAL, ADLTx, m)->getFitness();
    }
    
    // compute mean and var
    _meanW[i] = ARRAY::mean(array, size);
    _varW[i]  = ARRAY::var(array, size, _meanW[i]);
    delete[] array;
}

// ----------------------------------------------------------------------------------------
// getVwB    fitness variance between patches
// ----------------------------------------------------------------------------------------
double
MetapopSH::getVwB(){
    setMeanAndVar_W();
    return ARRAY::varUnbiased(_meanW, get_current_nbSamplePatch());
}

// ----------------------------------------------------------------------------------------
// getVwB   fitness variance within patches
// ----------------------------------------------------------------------------------------
double
MetapopSH::getVwW(){
    setMeanAndVar_W();
    return ARRAY::mean(_varW, get_current_nbSamplePatch());
}





