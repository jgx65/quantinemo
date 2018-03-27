/** @file tlocus.h
 *
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


//---------------------------------------------------------------------------

#ifndef tlocusH
#define tlocusH
//---------------------------------------------------------------------------
#include "types.h"
class TTraitProto;
class TPatch;
class TMetapop;
//---------------------------------------------------------------------------
/** this class contains the information of a locus position charachterized by the chromosme and
	* the genetic distance on in cM (separately for both sexes)
	* unlinked loci: chrom has a value of NaN and pos the index of the unlinked locus of the trait (0. umlinked locus, 1. unlinked locus,...)
	* Pleiotropic loci point to the same object
 */
class TLocusPosition{
public:
    unsigned int chrom;    // index of the chromosome (starts with 0)
    double pos[2];         // position on the chromosome in cM from the beggining of the chromosome
    
    TLocusPosition(){chrom = pos[FEM] = pos[MAL] = my_NAN;}
    TLocusPosition(const unsigned int& c, const double& f){
        chrom = c; pos[FEM] = pos[MAL] = f;
    }
    TLocusPosition(const unsigned int& c, const double& d, const sex_t& s){
        chrom = c; pos[FEM] = pos[MAL] = my_NAN; pos[s] = d;
    }
    TLocusPosition(const unsigned int& c, const double& f, const double& m){
        chrom = c; pos[FEM] = f; pos[MAL] = m;
    }
    
    bool isUnlinked(){return chrom == my_NAN;}
    
    bool operator<(const TLocusPosition& l) const {
        if(chrom < l.chrom) return true;
        if(chrom > l.chrom) return false;
        if(pos[FEM] < l.pos[FEM]){
            if(pos[MAL] <= l.pos[MAL]) return true;
            throw("Locus positions: loci positions have to be in increasing distance for both sexes!\n");
        }
        if(pos[MAL] < l.pos[MAL]){
            if(pos[FEM] <= l.pos[FEM]) return true;
            throw("Locus positions: loci positions have to be in increasing distance for both sexes!\n");
        }
        return false;
    }
    bool operator>(const TLocusPosition& l) const {return !(*this<=l);}
    
    bool operator<=(const TLocusPosition& l) const {
        if(chrom < l.chrom) return true;
        if(chrom > l.chrom) return false;
        if(pos[FEM] <= l.pos[FEM]){
            if(pos[MAL] <= l.pos[MAL]) return true;
            throw("Locus positions: loci positions have to be in increasing distance for both sexes!\n");
        }
        return false;
    }
    bool operator>=(const TLocusPosition& l) const {return !(*this<l);}
    
    bool operator==(const TLocusPosition& l) const {
        if(chrom != l.chrom) return false;
        if(pos[FEM] == l.pos[FEM] && pos[MAL] == l.pos[MAL]) return true;
        return false;
    }
    bool operator!=(const TLocusPosition& l) const {return !(*this==l);}
    
    
};

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
/** object for each locus containing the locus specific information */
class TLocus{
private:
    TLocusPosition*   _locusPos;             // locus coordinates on the genome (NULL if unlinked)
    unsigned int _nb_allele;            // maximal number of alleles
    unsigned int _locus_id_trait;       // the locus id of the trait
    unsigned int _locus_id_tot;         // the locus id of the genome
    
    double* _allelic_effect;            // if quantitative trait this pointer points to the alleic effects
    TTraitProto* _pTrait;               // pointer to the corresponding trait
    TMetapop* _popPtr;
    
    // mutations
    double _mut_rate;
    double* _mutationFreq;       		// _mutationFreq[allele]
    mut_model_t _mut_model;             // mutation model
    void (TLocus::*_mut_model_func_ptr)(unsigned char* seq);
    void  _mutate_none                 (unsigned char* seq){ }
    void  _mutate_KAM                  (unsigned char* seq);
    void  _mutate_SSM                  (unsigned char* seq);
    void  _mutate_RMM                  (unsigned char* seq);
    void  _mutate_IMM                  (unsigned char* seq);
    bool (TLocus::*_mut_model2_func_ptr)(unsigned char* seq, const double& ran);
    bool  _mutate_RMM                   (unsigned char* seq, const double& ran);
    bool  _mutate_IMM                   (unsigned char* seq, const double& ran);
    
    // sequence initialization
    double** _ini_frequencies; 	// frequencies to initialize a locus _ini_frequencies[patch][allele](only used for _ini_sequence_dist())
    void (TLocus::*_ini_seq_model_func_ptr)(unsigned char* seq, TPatch* patch, unsigned int size);
    void  _ini_sequence_monomorph      (unsigned char* seq, TPatch* patch, unsigned int size=ploidy);
    void  _ini_sequence_uniform        (unsigned char* seq, TPatch* patch, unsigned int size=ploidy);
    void  _ini_sequence_dist           (unsigned char* seq, TPatch* patch, unsigned int size=ploidy);
    
    // mutational correlation
private:
    double* _MNR_param;                 // NULL if not used
    
public:
    double* get_MNR_param(){return _MNR_param;}
    
    TLocus();
    ~TLocus();
    
    void set_locusPosition(TLocusPosition* l){_locusPos = l;}
    void set_mutationRate(const double& m){_mut_rate = m;}
    void set_mutationModel(const mut_model_t& m);
    void set_nb_allele(const unsigned int& n){_nb_allele = n;}
    void set_locus_id_trait(const unsigned int& n){_locus_id_trait = n;}
    void set_locus_id_tot(const unsigned int& n);
    void set_traitPointer(TTraitProto* p);
    void set_ini_sequence_model(const ini_model_t& model);
    void set_ini_sequence_model(const unsigned int& model);
    void set_ini_frequencies(double** pointer){_ini_frequencies=pointer;}
    
    unsigned int 	get_chromosomePosition()  {return _locusPos ? _locusPos->chrom : my_NAN;}
    double      	get_locusPosition(sex_t s){return _locusPos ? _locusPos->pos[s]: my_NAN;}
    TLocusPosition*  	get_locusPosition()       {return _locusPos;}
    unsigned int 	get_nb_allele()			  {return _nb_allele;}
    unsigned int 	get_locus_id_trait()	  {return _locus_id_trait;}
    unsigned int 	get_locus_id_tot()		  {return _locus_id_tot;}
    double          get_mut_rate()            {return _mut_rate;}
    mut_model_t     get_mutationModel()       {return _mut_model;}
    TTraitProto*    get_pTrait()              {return _pTrait;}
    
    void mutate(unsigned char* seq);
    void mutate_now(unsigned char* seq);
    bool mutate_now(unsigned char* seq, const double& rand);
    void ini_sequence(unsigned char* seq, TPatch* p, unsigned int s=ploidy){(this->*_ini_seq_model_func_ptr)(seq, p, s);}
    
    void reset();
    
};

#endif
