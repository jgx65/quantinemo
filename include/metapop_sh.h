/** @file metapopSH.h
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

#ifndef metapop_shH
#define metapop_shH

#include "stathandler.h"

class TIndividual;

/**A StatHandler for the Metapop SimComponent.*/
class MetapopSH : public StatHandler<MetapopSH> {
private:
    double meanEmigrant,meanImmigrant,meanResidant,meanDeadDisp;
    unsigned int _compEmigrant[3], _compImmigrant[3], _compKolonisers[3], _compDeadDisp[3];
    double _mean_reprod_success;
    double _var_reprod_success;
    
    // fitness
    double *_meanW, *_varW;
    
public:
    
    MetapopSH( ) : _meanW(0), _varW(0){}
    
    virtual ~MetapopSH() {
        if(_meanW)    delete[] _meanW;
        if(_varW)     delete[] _varW;
    }
    
    virtual bool setStatRecorders(const string& token);
    virtual bool set_stat_demography(const string& t,const string& token,const string& end,const age_t& AGE,const string& ageStr);
    virtual string getName() {return "MetapopSH";}
    
    ///@name Migration
    ///@{
    double getMeanEmigrantPerPatch           ();
    double getMeanImmigrantPerPatch          ();
    double getMeanMigrantRatio               ();
    double getMeanResidantPerPatch           ();
    double getMeanKolonisersProportion       ();
    double getMeanKolonisersPerPatch         ();
    ///@}
    
    ///@name Patch extinction
    ///@{
    double getObsrvdExtinctionRate           ();
    double get_isAlive                       ();
    ///@}
    
    ///@name Demography
    ///@{
    void setReproductiveStats         (bool sex);
    double getReproductiveMean        (bool sex) {setReproductiveStats(sex); return _mean_reprod_success;}
    double getReproductiveVar         (bool sex) {setReproductiveStats(sex); return _var_reprod_success;}
    
    // stats across all inds and pops
    double getNbIndTot                (const age_t& AGE);
    double getNbFemTot                (const age_t& AGE);
    double getNbMalTot                (const age_t& AGE);
    
    double getMeanNbIndTot            (const age_t& AGE);
    double getMeanNbFemTot            (const age_t& AGE);
    double getMeanNbMalTot            (const age_t& AGE);
    
    double getNbIndTot                (unsigned int i, const age_t& AGE);
    double getNbFemTot                (unsigned int i, const age_t& AGE);
    double getNbMalTot                (unsigned int i, const age_t& AGE);
    
    double getNbPopsTot               (const age_t& AGE);
    
    double getSexRatioTot             (const age_t& AGE);
    
    // stats across sampled inds and pops
    double getNbInd                   (const age_t& AGE);   // total number of sampled individuals
    double getNbFem                   (const age_t& AGE);   // total number of sampled females
    double getNbMal                   (const age_t& AGE);   // total number of sampled males
    
    double getMeanNbInd               (const age_t& AGE);   // mean number of sampled individuals across popualted pops
    double getMeanNbFem               (const age_t& AGE);   // mean number of sampled females across popualted pops
    double getMeanNbMal               (const age_t& AGE);   // mean number of sampled males across popualted pops
    
    double getNbInd                   (unsigned int i, const age_t& AGE);
    double getNbFem                   (unsigned int i, const age_t& AGE);
    double getNbMal                   (unsigned int i, const age_t& AGE);
    
    double getNbPops                  (const age_t& AGE);
    
    double getSexRatio                (const age_t& AGE);
    ///@}
    
    ///@name Kinship
    ///@{
    void   setKinship                 (const age_idx& AGE);
    void   setKinClassCounter         (TIndividual *I1, TIndividual *I2, const age_idx& AGE);
    double getSibProportion           (unsigned int i, const age_t& AGE) {
        age_idx curAGE = age_t2idx(AGE);
        setKinship(curAGE);
        return _sib_prop[curAGE][i];
    }
    bool   set_stat_kinship(const string& t,const string& token,const string& end,const age_t& AGE,const string& ageStr);
    ///@}
    
    // fitness
    // fitness
    double getMeanW     (unsigned int i)  {
        setMeanAndVar_W();
        return _meanW[i];
    }
    double getVarW      (unsigned int i)  {
        setMeanAndVar_W();
        return _varW[i];
    }
    double getVwB();
    double getVwW();
    void setMeanAndVar_W();
    void setMeanAndVar_W_ofPatch(Patch* crnt_patch, const unsigned int& i);
};

#endif

