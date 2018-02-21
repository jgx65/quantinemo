/** @file stathandler.h
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

#ifndef stathandlerH
#define stathandlerH

#include "stathandlerbase.h"
#include "stat_rec_base.h"

class TTraitProto;
class TIndividual;
class TPatch;
/**A class to compute and store the summary statistics associated with a SimComponent.
 * The template type must be the type of the class that declares the methods linked into the
 * StatRecorder elements. */
template <class SH>
class StatHandler: public StatHandlerBase {
    
private:
    double get_coancestry(TPatch*, const sex_t&, TPatch*, const sex_t&, const age_idx&);
    
protected:
    /**The list of stat recorders.*/
	list<StatRecorder<SH>*> _recorders;
    
    TTraitProto* _SHLinkedTrait;
    unsigned int _SHLinkedTraitIndex;
	unsigned int _trait_index; // for multiple instansiations of the corresponding trait
    
    typedef typename list< StatRecorder<SH>* >::iterator REC_IT;
    
	// allele frequencies / counts
	unsigned int _nb_patch;
	unsigned int _nb_locus;
	unsigned int* _nb_allele;         // number of allele per locus (do not delete this array, it is just pointing to another array)
	unsigned int _nb_allele_max;      // the highest number of alleles across loci
    
	map<unsigned char, double>*** _alleleFreq_local;  // [age][patch][locus][allele]
	map<unsigned char, double>**  _alleleFreq_global; // [age][locus][allele]
	map<unsigned char, map<unsigned char, double> >*** _locusFreq_local;  // [age][patch][locus][allele1][allele2]
	map<unsigned char, map<unsigned char, double> >**  _locusFreq_global; // [age][locus][allele1][allele2]
    
	map<unsigned char, map<unsigned char, double> >* get_genotypeFreq(const age_idx& AGE, TPatch* p,
                                                                      const unsigned int& l1, const unsigned int& l2);
    
    
	bool already_computed(unsigned int* array, age_t AGE=NONE);
	bool already_computed(unsigned int* array, const age_idx& AGE);
    
	bool set_stat_coancestry     (string t, string i, string trait, string token, string end, age_t, string);
	bool set_stat_fstat          (string t, string i, string trait, string token, string end, age_t, string);
	bool set_stat_all_freq_local (string t, string i, string trait, string token, string end, age_t, string);
	bool set_stat_all_freq_global(string t, string i, string trait, string token, string end, age_t, string);
	bool set_stat_locus_freq_local (string t, string i, string trait, string token, string end, age_t, string);
	bool set_stat_locus_freq_global(string t, string i, string trait, string token, string end, age_t, string);
	bool set_stat_LD             (string t, string i, string trait, string token, string end, age_t, string);
    
    void set_stat_all_freq_local_text(string t, string i, string trait,
                                      string ageStr, unsigned int p,
                                      unsigned int l, unsigned int a,
                                      string& n, string & _i);
    void set_stat_all_freq_global_text(string t, string i, string trait,
                                       string ageStr,
                                       unsigned int l, unsigned int a,
                                       string& stat, string & text);
    void set_stat_locus_freq_local_text(string t, string i, string trait,
                                   string ageStr, unsigned int p,
                                   unsigned int l, unsigned int a1,
                                   unsigned char a2,
                                        string& stat, string & text);
    void set_stat_locus_freq_global_text(string t, string i, string trait,
                                   string ageStr,
                                   unsigned int l, unsigned int a1,
                                   unsigned char a2,
                                         string& stat, string & text);
    
	///@}
    ///@name Coancestries
	///@{
	double   Theta_FF[NB_AGE_CLASSES], Theta_MM[NB_AGE_CLASSES], Theta_FM[NB_AGE_CLASSES];
	double   _mean_theta[NB_AGE_CLASSES], _mean_alpha[NB_AGE_CLASSES];
	double*** _coa_matrix;            				// _coa_matrix[age][patch1][patch2]
    
    /**Kinship classes proportions*/
	double _sib_prop[NB_AGE_CLASSES][5];    // 0: non sib, 1: maternal half sibs; 2: paternal half sibs; 3: full sib; 4: selfed
	double _sib_coa[NB_AGE_CLASSES][5];     // 0: non sib, 1: maternal half sibs; 2: paternal half sibs; 3: full sib; 4: selfed
    
    
    /**Gives the coancestry (probability of identity by state) of two gene sequences.
     The probability returned is the averageproportion of identical alleles per locus between the two sequences.
     @param ind1 first sequence, treated as of type (unsigned char**)
     @param ind2 second sequence, treated as of type (unsigned char**)
     */
    double Coancestry               (unsigned char** seq1, unsigned char** seq2);
    
    /**Computes the within and between patches coancestry coefficients.
     @param age_pos the age class index
     @param dim the dimension of the matrix to fill:
     - 1 = the diagonal (i.e. the within patch coancestries or theta's)
     - 2 = the upper half (i.e. the between patch coancestries or alpha's)
     - 3 = both
     */
	void   setCoaMatrixTheta             (const age_idx& AGE);
	void   setCoaMatrixAlpha             (const age_idx& AGE);
	void   setSexspecific_Theta          (const age_idx& AGE);
    
    /**Gets the given coancestry coefficient from the coancestry matrix.
     @param i combination of the row and column indexes (see setCoaMatrixRecorders()).
     \note the upper half and the diagonal of the matrix are filled, other positions are set to 0.
     */
	double getCoaTheta (unsigned int i, const age_t& AGE){
		age_idx curAge=age_t2idx(AGE);
        unsigned int sampleID1 = get_vPatch(i)->get_sampleID();
       //if this statistic cannot be computed, return my_Nan
        if(sampleID1==my_NAN||sampleID1==SAMPLED){
            return my_NAN;
        }
        setCoaMatrixTheta(curAge);
        return _coa_matrix[curAge][sampleID1][sampleID1];           // diagonal

    }
	double getCoaAlpha (unsigned int i, const age_t& AGE){
		age_idx curAge=age_t2idx(AGE);
		setCoaMatrixAlpha(curAge);
		unsigned int id1, id2;
		if(!id2sampleIDs(i, id1, id2))return my_NAN;	// if at least one patch is not sampled
		return _coa_matrix[curAge][id1][id2];           // pairwise
	}
	double getMeanTheta(const age_t& AGE){
		age_idx curAge = age_t2idx(AGE);
		setCoaMatrixTheta(curAge);
		return _mean_theta[curAge];
	}
	double getMeanAlpha(const age_t& AGE){
		age_idx curAge = age_t2idx(AGE);
		setCoaMatrixAlpha(curAge);
		return _mean_alpha[curAge];
	}
    
	/**Gives the mean within females coancestry coefficient.*/
	double getTheta_FF              (const age_t& AGE){
		age_idx curAge = age_t2idx(AGE);
		setSexspecific_Theta(curAge);
		return Theta_FF[curAge];
	}
    
	/**Gives the mean within males coancestry coefficient.*/
	double getTheta_MM              (const age_t& AGE){
		age_idx curAge = age_t2idx(AGE);
		setSexspecific_Theta(curAge);
		return Theta_MM[curAge];
	}
    
	/**Gives the mean between males and females coancestry coefficient.*/
	double getTheta_FM              (const age_t& AGE){
		age_idx curAge = age_t2idx(AGE);
		setSexspecific_Theta(curAge);
		return Theta_FM[curAge];
	}
    
	void   setSibStats              (const age_idx& AGE);
	void   setSibCoa                (TIndividual *I1, TIndividual *I2, const age_idx& AGE);
	double getSibProportions        (unsigned int i, const age_t& AGE) {
		age_idx curAGE = age_t2idx(AGE);
		setSibStats(curAGE);
		return _sib_prop[curAGE][i];
	}
	double getSibCoaMeans           (unsigned int i, const age_t& AGE) {
		age_idx curAGE = age_t2idx(AGE);
		setSibStats(curAGE);
		return _sib_coa[curAGE][i];
	}
	///@}
    
    
    //fstat:
	double _ho[NB_AGE_CLASSES], _hs[NB_AGE_CLASSES], _ht[NB_AGE_CLASSES],
    _hsnei[NB_AGE_CLASSES], _htnei[NB_AGE_CLASSES],
    _fst[NB_AGE_CLASSES], _fis[NB_AGE_CLASSES], _fit[NB_AGE_CLASSES];
	double _fst_wc[NB_AGE_CLASSES], _fis_wc[NB_AGE_CLASSES], _fit_wc[NB_AGE_CLASSES];
	double*** _fst_matrix_wc, ***_fst_matrix;
    
	double **_hsnei_locus, **_htnei_locus, **_fst_locus, **_fis_locus, **_fit_locus,  // for Nei and Chesser
    **_ho_locus, **_hs_locus, **_ht_locus;
    
	double **_fst_WC_locus, **_fis_WC_locus, **_fit_WC_locus;                       // Weir and Cockerham
    
	double **_het0_locus, **_het1_locus, **_fst_bn_locus;                           // fdist: Beaumont & Nichols 1996
    
    
	unsigned int** _computed;    // computed[stat][option]: 0: gen; 1: rep; 2: age
	unsigned int _computed_size;
    
public:
    
	StatHandler( ):_SHLinkedTrait(0),  _trait_index(0),
    _nb_allele(0), _nb_allele_max(0),
    _alleleFreq_local(0), _alleleFreq_global(0),
    _locusFreq_local(0), _locusFreq_global(0),
    _coa_matrix(0), _fst_matrix_wc(0),_fst_matrix(0),
    _hsnei_locus(0), _htnei_locus(0), _fst_locus(0), _fis_locus(0), _fit_locus(0),
    _ho_locus(0), _hs_locus(0), _ht_locus(0), _fst_WC_locus(0), _fis_WC_locus(0), _fit_WC_locus(0),
    _het0_locus(0), _het1_locus(0), _fst_bn_locus(0){
        
		_computed_size = 35;
		_computed = new unsigned int*[_computed_size];
		for(unsigned int i=0; i<_computed_size; ++i){
			_computed[i] = new unsigned int[3];
			_computed[i][0] = 0;
			_computed[i][1] = 0;
			_computed[i][2] = 0;
		}
	}
    
	void set(TTraitProto* TT);
    
	virtual ~StatHandler ( );
    
	virtual bool init();
    
	virtual string getName() {return "StatHandler";}
    
	/**Empties the _recorders list, they are destroyed in StatHandlerBase::reset().*/
	virtual void clear   ( )  {_recorders.clear();}
    
	virtual unsigned int get_nb_stats ( ){return (unsigned int)_recorders.size();}
    
	/**Computes the stats by executing the function variables stored in the StatRecorder's.*/
    void execute ( );              // if stats are stored in db
    void execute (ostream& FH);    // if stats are directly written to file
    void execute_param ( );              // if stats are stored in db
    void execute_param (ostream& FH);    // if stats are directly written to file
    
	/** Adds a StatRecorder to the list (for each pairwise patch combination: postfix: _pi-j) */
	bool add_pairwisePatch (const string& end, const string& Title, const string& Name,
                            const st_order& Order, const age_t& AGE, const unsigned int& ARG,
                            double(SH::* getStat)(void),
                            double(SH::* getStatBoolArg)(bool)=0,
                            double(SH::* getStatUintArg)(unsigned int)=0,
                            double(SH::* getStatAGE)(const age_t&)=0,
                            double(SH::* getStatBoolArgAGE)(bool, const age_t&)=0,
                            double(SH::* getStatUintArgAGE)(unsigned int, const age_t&)=0);
    
    bool add_pairwisePatch_ij (const string& end, const string& Title, const string& Name,
                                       const st_order& Order, const age_t& AGE, const unsigned int& ARG,
                                       double(SH::* getStat)(void),
                                       double(SH::* getStatBoolArg)(bool)=0,
                                       double(SH::* getStatUintArg)(unsigned int)=0,
                                       double(SH::* getStatAGE)(const age_t&)=0,
                                       double(SH::* getStatBoolArgAGE)(bool, const age_t&)=0,
                                       double(SH::* getStatUintArgAGE)(unsigned int, const age_t&)=0);
    
	/** Adds a StatRecorder to the list (for each patch: postfix: _pi) */
	bool add_perPatch (const string& end, const string& Title, const string& Name,
                       const st_order& Order, const age_t& AGE, const unsigned int& ARG,
                       double(SH::* getStat)(void),
                       double(SH::* getStatBoolArg)(bool)=0,
                       double(SH::* getStatUintArg)(unsigned int)=0,
                       double(SH::* getStatAGE)(const age_t&)=0,
                       double(SH::* getStatBoolArgAGE)(bool, const age_t&)=0,
                       double(SH::* getStatUintArgAGE)(unsigned int, const age_t&)=0);
    
	bool add_perPatch_i (const string& end, const string& Title, const string& Name,
                         const st_order& Order, const age_t& AGE, const unsigned int& ARG,
                         double(SH::* getStat)(void),
                         double(SH::* getStatBoolArg)(bool)=0,
                         double(SH::* getStatUintArg)(unsigned int)=0,
                         double(SH::* getStatAGE)(const age_t&)=0,
                         double(SH::* getStatBoolArgAGE)(bool, const age_t&)=0,
                         double(SH::* getStatUintArgAGE)(unsigned int, const age_t&)=0);
    
    bool add_perPatch_andLocus (const string& end, const string& Title, const string& Name,
                                        const st_order& Order, const age_t& AGE, const unsigned int& ARG,
                                        double(SH::* getStat)(void),
                                        double(SH::* getStatBoolArg)(bool)=0,
                                        double(SH::* getStatUintArg)(unsigned int)=0,
                                        double(SH::* getStatAGE)(const age_t&)=0,
                                        double(SH::* getStatBoolArgAGE)(bool, const age_t&)=0,
                                        double(SH::* getStatUintArgAGE)(unsigned int, const age_t&)=0);
    
	/** Adds a StatRecorder to the list (for each locus: postfix: _li) */
	bool add_perLocus (const string& end, const string& Title, const string& Name,
                       const st_order& Order, const age_t& AGE, const unsigned int& ARG,
                       double(SH::* getStat)(void),
                       double(SH::* getStatBoolArg)(bool)=0,
                       double(SH::* getStatUintArg)(unsigned int)=0,
                       double(SH::* getStatAGE)(const age_t&)=0,
                       double(SH::* getStatBoolArgAGE)(bool, const age_t&)=0,
                       double(SH::* getStatUintArgAGE)(unsigned int, const age_t&)=0);
    
	bool add_perLocus_i (const string& end, const string& Title, const string& Name,
                         const st_order& Order, const age_t& AGE, const unsigned int& ARG,
                         double(SH::* getStat)(void),
                         double(SH::* getStatBoolArg)(bool)=0,
                         double(SH::* getStatUintArg)(unsigned int)=0,
                         double(SH::* getStatAGE)(const age_t&)=0,
                         double(SH::* getStatBoolArgAGE)(bool, const age_t&)=0,
                         double(SH::* getStatUintArgAGE)(unsigned int, const age_t&)=0);
    
	/** Adds a StatRecorder to the list (for each pairwise locus combination: postfix: _li-j) */
	bool add_pairwiseLocus (const string& end, const string& Title, const string& Name,
                            const st_order& Order, const age_t& AGE, const unsigned int& ARG,
                            double(SH::* getStat)(void),
                            double(SH::* getStatBoolArg)(bool)=0,
                            double(SH::* getStatUintArg)(unsigned int)=0,
                            double(SH::* getStatAGE)(const age_t&)=0,
                            double(SH::* getStatBoolArgAGE)(bool, const age_t&)=0,
                            double(SH::* getStatUintArgAGE)(unsigned int, const age_t&)=0);
    
	bool add_pairwiseLocus_ij (const string& end, const string& Title, const string& Name,
                               const st_order& Order, const age_t& AGE, const unsigned int& ARG,
                               double(SH::* getStat)(void),
                               double(SH::* getStatBoolArg)(bool)=0,
                               double(SH::* getStatUintArg)(unsigned int)=0,
                               double(SH::* getStatAGE)(const age_t&)=0,
                               double(SH::* getStatBoolArgAGE)(bool, const age_t&)=0,
                               double(SH::* getStatUintArgAGE)(unsigned int, const age_t&)=0);
    
	bool add_pairwiseLocus_perPatch (const string& end, const string& Title, const string& Name,
                                     const st_order& Order, const age_t& AGE, const unsigned int& ARG,
                                     double(SH::* getStat)(void),
                                     double(SH::* getStatBoolArg)(bool)=0,
                                     double(SH::* getStatUintArg)(unsigned int)=0,
                                     double(SH::* getStatAGE)(const age_t&)=0,
                                     double(SH::* getStatBoolArgAGE)(bool, const age_t&)=0,
                                     double(SH::* getStatUintArgAGE)(unsigned int, const age_t&)=0);
    
	bool add_pairwiseLocus_perPatch_p (const string& end, const string& Title, const string& Name,
                                       const st_order& Order, const age_t& AGE, const unsigned int& ARG,
                                       double(SH::* getStat)(void),
                                       double(SH::* getStatBoolArg)(bool)=0,
                                       double(SH::* getStatUintArg)(unsigned int)=0,
                                       double(SH::* getStatAGE)(const age_t&)=0,
                                       double(SH::* getStatBoolArgAGE)(bool, const age_t&)=0,
                                       double(SH::* getStatUintArgAGE)(unsigned int, const age_t&)=0);
    
	bool add_pairwisePatch_perLocus (const string& end, const string& Title, const string& Name,
                                     const st_order& Order, const age_t& AGE, const unsigned int& ARG,
                                     double(SH::* getStat)(void),
                                     double(SH::* getStatBoolArg)(bool)=0,
                                     double(SH::* getStatUintArg)(unsigned int)=0,
                                     double(SH::* getStatAGE)(const age_t&)=0,
                                     double(SH::* getStatBoolArgAGE)(bool, const age_t&)=0,
                                     double(SH::* getStatUintArgAGE)(unsigned int, const age_t&)=0);
    
	bool add_pairwisePatch_perLocus_l (const string& end, const string& Title, const string& Name,
                                       const st_order& Order, const age_t& AGE, const unsigned int& ARG,
                                       double(SH::* getStat)(void),
                                       double(SH::* getStatBoolArg)(bool)=0,
                                       double(SH::* getStatUintArg)(unsigned int)=0,
                                       double(SH::* getStatAGE)(const age_t&)=0,
                                       double(SH::* getStatBoolArgAGE)(bool, const age_t&)=0,
                                       double(SH::* getStatUintArgAGE)(unsigned int, const age_t&)=0);
    
    /** Adds a StatRecorder to the list (a single stat) */
    bool add (const string& Title, string Name,
              const st_order& Order, const age_t& AGE, const unsigned int& ARG,
              double(SH::* getStat)(void),
              double(SH::* getStatBoolArg)(bool)=0,
              double(SH::* getStatUintArg)(unsigned int)=0,
              double(SH::* getStatAGE)(const age_t&)=0,
              double(SH::* getStatBoolArgAGE)(bool, const age_t&)=0,
              double(SH::* getStatUintArgAGE)(unsigned int, const age_t&)=0);
    
    
	// allele frequencies
	double  get_allele_freq_local(unsigned int i, const age_t& AGE);
	double  get_allele_freq_global(unsigned int i, const age_t& AGE);
	bool    set_alleleFreq(const age_idx& AGE);
	unsigned int set_alleleFreq_ofPatch(TPatch* crnt_patch, const age_idx& AGE, map<unsigned char, double>*& freqs, map<unsigned char, double>*& global_freqs);
    unsigned int set_alleleFreq_ofPatch_allInds(TPatch* crnt_patch, const age_idx& AGE,
                                                map<unsigned char, double>*& freqs,
                                                map<unsigned char, double>*& global_freqs);
	void    set_alleleFreq_global_ofLocus(const unsigned int& l,
                                          map<unsigned char, double>** localFreqs, map<unsigned char, double>* globalFreqs,
                                          unsigned int* popSizes, const unsigned int& tot_popSize, const unsigned int& nbPatch);
	void    set_alleleFreq_global(map<unsigned char, double>** localFreqs, map<unsigned char, double>* globalFreqs,
                                  unsigned int* popSizes, const unsigned int& tot_popSize, const unsigned int& nbPatch);
    
	double  getHarmonicMean_ofPopSize(const age_idx& AGE, const vector<TPatch*>& aPatch, unsigned int& nbPopFull); // for any number of patches
    
    // locus genotype frequencies
	double  get_locus_freq_local(unsigned int i, const age_t& AGE);
	double  get_locus_freq_global(unsigned int i, const age_t& AGE);
	bool    set_locusFreq(const age_idx& AGE);
	unsigned int set_locusFreq_ofPatch(TPatch* crnt_patch, const age_idx& AGE,
                                       map<unsigned char, map<unsigned char, double> >*& freqs,
                                       map<unsigned char, map<unsigned char, double> >*& global_freqs);
    unsigned int get_locusGenotypeFreqs_ofPatch(TPatch * curPatch, const age_idx & AGE,
                                                map<unsigned char, map<unsigned char, double> >*& freqs);
    unsigned int get_locusGenotypeFreqs_ofPatch(TPatch * curPatch, const age_idx & AGE,
                                                unsigned int traitID,
                                                map<unsigned char, map<unsigned char, double> >*& freqs);
    unsigned int get_locusGenotypeFreqs_ofPatch_allInds(TPatch * curPatch, const age_idx & AGE,
                                                unsigned int traitID,
                                                map<unsigned char, map<unsigned char, double> >*& freqs);
    unsigned int get_locusGenotypeCounts_ofPatch(TPatch * curPatch, const age_idx & AGE,
                                                 unsigned int traitID,
                                                map<unsigned char, map<unsigned char, double> >*& freqs);
    unsigned int get_locusGenotypeCounts_ofPatch_andSex(TPatch * curPatch, const age_idx & AGE, sex_t SEX,
                                                        unsigned int traitID,
                                                        map<unsigned char, map<unsigned char, double> >*& freqs);
    
	///@name F-stats:
	///@{
    
	/** Computes the weighted within and between patch Fst's as well as the overall Fst (Theta) */
	void   setFstat_Weir_Cockerham          (const age_idx& AGE);      // across loci
	void   setFstat_Weir_Cockerham(const age_idx& AGE, vector<TPatch*>& aPatch,
                                   map<unsigned char, double>** alleleFreqs,
                                   map<unsigned char, double>* alleleFreqsGlobal,
                                   double& fst, double& fis, double &fit);
	void   setFstat_Weir_Cockerham_perLocus (const age_idx& AGE);      // for each locus separately
	void   setFstat_Weir_Cockerham_perPatchPair(const age_idx& AGE);
	double getFstat_Weir_Cockerham_perPatchPair(const age_idx& AGE, TPatch* pop1, TPatch* p2);
	double getFstat_Weir_Cockerham_perPatchPair_andLocus(const age_idx& AGE, TPatch* pop1, TPatch* p2, unsigned int l);
	void   get_sigma_of_locus_Weir_Cockerham(double& sigma_a2,
                                             double& sigma_b2, double& sigma_w2, const unsigned int& nbPatch, const unsigned int& nbPop,
                                             map<unsigned char, double>** alleleFreqs, map<unsigned char, double>* alleleFreqsGlobal,
                                             const unsigned int& l, const unsigned int& tot_size, unsigned int* pop_sizes,
                                             map<unsigned char, double>** ho_patch_allele, const double& nc);
    
    /** F-statistics following Weir & Cockeram (1984) */
	double getFst_WC (const age_t& AGE) {age_idx curAGE=age_t2idx(AGE); setFstat_Weir_Cockerham(curAGE); return _fst_wc[curAGE];}
	double getFit_WC (const age_t& AGE) {age_idx curAGE=age_t2idx(AGE); setFstat_Weir_Cockerham(curAGE); return _fit_wc[curAGE];}
	double getFis_WC (const age_t& AGE) {age_idx curAGE=age_t2idx(AGE); setFstat_Weir_Cockerham(curAGE); return _fis_wc[curAGE];}
    
	//////////////////////////////////////////////////////////////////////////////
	/** To find the corresponding locus, patch, or whatever a single integer is used.
     * These functions here allow to extract fromv this number the corresponding ids if
     * several quantities are needed, eg pairwise patches.
     */
    
	/* 2 args to id */
	inline unsigned int toID(unsigned int id1, unsigned int id2,
                             unsigned int tot1, unsigned int tot2)
	{
		return id1 + tot1*id2;
	}
	/* id to 2 args */
	inline void fromID(unsigned int in, unsigned int& id1, unsigned int& id2,
                       unsigned int tot1, unsigned int tot2)
	{
		assert(in < tot1*tot2);
		id2 = (unsigned int) in/tot1;
		id1 = in - id2*tot1;               //in%tot1;
	}
    
	/* 3 args to id */
	inline unsigned int toID(unsigned int id1, unsigned int id2, unsigned int id3,
                             unsigned int tot1, unsigned int tot2, unsigned int tot3)
	{
        return toID(id1, id2, tot1, tot2) + tot1*tot2*id3;
	}
    
	/* id to 3 args */
	inline void fromID(unsigned int in, unsigned int& id1, unsigned int& id2, unsigned int& id3,
                       unsigned int tot1, unsigned int tot2, unsigned int tot3)
	{
		assert(in < tot1*tot2*tot3);
		unsigned int sum = tot1*tot2;
		id3 = (unsigned int) in/sum;
        fromID(in-id3*sum, id1, id2, tot1, tot2);
	}
    
    
	/* 4 args to id */
	inline unsigned int toID(unsigned int id1, unsigned int id2, unsigned int id3, unsigned int id4,
                             unsigned int tot1, unsigned int tot2, unsigned int tot3, unsigned int tot4)
	{
        return toID(id1, id2, id3, tot1, tot2, tot3) + tot1*tot2*tot3*id4;
	}
    
	/* id to 4 args */
	inline void fromID(unsigned int in, unsigned int& id1, unsigned int& id2, unsigned int& id3, unsigned int& id4,
                       unsigned int tot1, unsigned int tot2, unsigned int tot3, unsigned int tot4)
	{
		assert(in < tot1*tot2*tot3*tot4);
		unsigned int sum = tot1*tot2*tot3;
		id4 = (unsigned int) in/sum;
        fromID(in-id4*sum, id1, id2, id3, tot1, tot2, tot3);
	}
    
    unsigned int extract_id(string& text, char first, unsigned max, bool charOblig=true);
   
    /*	computes from "in" the two corresponding sampleID's or reverse.
     * returns false if one or the other patch are not populated
     */
	inline bool id2sampleIDs(const unsigned int& in, unsigned int& sampleID1, unsigned int& sampleID2)
	{
		fromID(in, sampleID1, sampleID2, _nb_patch, _nb_patch);
		sampleID1 = get_vPatch(sampleID1)->get_sampleID();
		sampleID2 = get_vPatch(sampleID2)->get_sampleID();
		if(sampleID1==SAMPLED || sampleID2==SAMPLED) return false; // currently sampled patch, but currently not populated
		if(sampleID1==my_NAN  || sampleID2==my_NAN)  return false; // it is a temporal sample parameter and patch is currently not sampled
		if(sampleID1>sampleID2) swap(sampleID1, sampleID2);
		return true;  // return true if both patches are populated
	}
    
    
    
	inline bool id2sampleIDs(const unsigned int& in, unsigned int& sampleID1, unsigned int& sampleID2, unsigned int& locus)
	{
		fromID(in, sampleID1, sampleID2, locus, _nb_patch, _nb_patch, _nb_locus);
		sampleID1 = get_vPatch(sampleID1)->get_sampleID();
		sampleID2 = get_vPatch(sampleID2)->get_sampleID();
		if(sampleID1==SAMPLED || sampleID2==SAMPLED) return false; // currently sampled patch, but currently not populated
		if(sampleID1==my_NAN  || sampleID2==my_NAN)  return false; // it is a temporal sample parameter and patch is currently not sampled
		if(sampleID1>sampleID2) swap(sampleID1, sampleID2);
		return true;  // return true if both patches are populated
	}
    
    
	double getFst_WC_ij(unsigned int i, const age_t& AGE){  // all pairs are computed
		age_idx curAge = age_t2idx(AGE);
		setFstat_Weir_Cockerham_perPatchPair(curAge);
		unsigned int id1, id2;
		if(!id2sampleIDs(i, id1, id2))return my_NAN;	// if at least one patch is not sampled
		return _fst_matrix_wc[curAge][id1][id2];
	}
	double getFst_WC_ij_single(unsigned int i, const age_t& AGE){  // just this pair is computed
		unsigned int sample1, sample2;
		fromID(i, sample1, sample2, _nb_patch, _nb_patch);
		return getFstat_Weir_Cockerham_perPatchPair(age_t2idx(AGE), get_vPatch(sample1), get_vPatch(sample2));
	}
	double getFst_WC_ij_l_single(unsigned int i, const age_t& AGE){  // just this pair is computed
		unsigned int id1, id2, l;
		if(!id2sampleIDs(i, id1, id2, l))return my_NAN;	// if at least one patch is not sampled
		return getFstat_Weir_Cockerham_perPatchPair_andLocus(age_t2idx(AGE), get_vPatch(id1), get_vPatch(id2), l);
	}
    
	double getFst_WC_ofLocus  (unsigned int i, const age_t& AGE) {
		age_idx curAge = age_t2idx(AGE);
		setFstat_Weir_Cockerham_perLocus(curAge);
		return _fst_WC_locus[curAge][i];
	}
	double getFis_WC_ofLocus  (unsigned int i, const age_t& AGE) {
		age_idx curAge = age_t2idx(AGE);
		setFstat_Weir_Cockerham_perLocus(curAge);
		return _fis_WC_locus[curAge][i];
	}
	double getFit_WC_ofLocus  (unsigned int i, const age_t& AGE) {
		age_idx curAge = age_t2idx(AGE);
		setFstat_Weir_Cockerham_perLocus(curAge);
		return _fit_WC_locus[curAge][i];
	}
    
    
    /** F-statistics following Nei and Chesser 1983 */
	void   setFstat_Nei_Chesser              (const age_idx& AGE);    // F-statistics globally across patch and loci
	void   setFstat_Nei_Chesser_perPatchPair (const age_idx& AGE);    // F-statistics per pair of patches and across loci
	double getFstat_Nei_Chesser_perPatchPair (const age_idx& AGE, TPatch* p1, TPatch* p2);    // F-statistics per pair of patches and across loci
	double getFstat_Nei_Chesser_perPatchPair_andLocus (const age_idx& AGE, TPatch* p1, TPatch* p2, unsigned int l);    // F-statistics per pair of patches and across loci
	void   setFstat_Nei_Chesser_perLocus     (const age_idx& AGE);    // F-statistics for each locus separately across patches
	double get_fixed_loci_global (const age_t& AGE);
	void   set_fixed_loci_perPatch(const age_t& AGE);
    
	double getHo                    (const age_idx& AGE);
	double getHo                    (const age_t& AGE) {return getHo(age_t2idx(AGE));}
	double*getHo_perLocus           (const age_idx& AGE);
	double getHo_ofLocus            (unsigned int i, const age_t& AGE){return getHo_perLocus(age_t2idx(AGE))[i];}
	double getHo_ofPatch            (TPatch* curPatch, const age_idx& AGE);
	double getHo_ofPatch            (unsigned int i, const age_t& AGE){return getHo_ofPatch(get_vPatch(i), age_t2idx(AGE));}
	double*getHo_ofPatch_andLocus    (const age_idx& AGE, TPatch* cur_patch, double* array=NULL);
	double getHo_ofPatch_andLocus    (const age_idx& AGE, TPatch* cur_patch, const unsigned int& l);
	double getHo_ofPatch_andLocus(unsigned int i, const age_t& AGE){
		unsigned int id, loc;
		fromID(i, id, loc, _nb_patch, _nb_locus);
		return getHo_ofPatch_andLocus(age_t2idx(AGE), get_vPatch(id), loc);
	}
	map<unsigned char, double>* getHo_ofPatchperAllele  (const age_idx& AGE, TPatch* cur_patch, map<unsigned char, double>* array=NULL);
    
	double getHs                    (const age_idx& AGE);
	double getHs                    (const age_t& AGE){return getHs(age_t2idx(AGE));}
	double getHs_ofPatch            (TPatch* p, const age_idx& AGE);
	double getHs_ofPatch            (unsigned int p, const age_t& AGE){return getHs_ofPatch(get_vPatch(p), age_t2idx(AGE));}
	double getHsUnbiased_ofPatch    (unsigned int p, const age_t& AGE);
	double getHsUnbiased_ofPatch_andLocus(const age_idx& AGE, TPatch* p, const unsigned int& l);
	double getHsUnbiased_ofPatch_andLocus(unsigned int i, const age_t& AGE){
		unsigned int id, loc;
		fromID(i, id, loc, _nb_patch, _nb_locus);
		return getHsUnbiased_ofPatch_andLocus(age_t2idx(AGE), get_vPatch(id), loc);
	}
	double getHs_ofLocus            (const unsigned int& l, const age_idx& AGE);
	double getHs_ofLocus            (unsigned int l, const age_t& AGE){return getHs_ofLocus(l, age_t2idx(AGE));}
	double getHs_ofPatch_andLocus   (const age_idx& AGE, TPatch* p, const unsigned int& l);
	double getHs_ofPatch_andLocus   (unsigned int i, const age_t& AGE){
		unsigned int id, loc;
		fromID(i, id, loc, _nb_patch, _nb_locus);
		return getHs_ofPatch_andLocus(age_t2idx(AGE), get_vPatch(id), loc);
	}
    
	double getHt                    (const age_idx& AGE);
	double getHt                    (const age_t& AGE){return getHt(age_t2idx(AGE));}
	double getHt                    (const age_idx& AGE, const unsigned int& p1, const unsigned int& p2);   // total Ht of the 2 patches
	double getHt                    (const age_idx& AGE, TPatch* p1, TPatch* p2);                             // total Ht of the 2 patches
	double getHt_ofLocus            (const age_idx& AGE, const unsigned int& l);
	double getHt_ofLocus            (unsigned int l, const age_t& AGE){return getHt_ofLocus(age_t2idx(AGE), l);}
	double getHt_ofLocus            (const age_idx& AGE, const unsigned int& l, const unsigned int& p1, const unsigned int& p2); // total Ht of the 2 patches
	double getHt_ofLocus            (const age_idx& AGE, const unsigned int& l, TPatch* p1, TPatch* p2); // total Ht of the 2 patches
    
	// allelic richness
	unsigned int get_minSampleSize  (const age_idx& AGE);
    
	double getRs                    (const age_idx& AGE);
	double getRs                    (const age_t& AGE){return getRs(age_t2idx(AGE));}
	double getRs_ofPatch            (TPatch* p, const age_idx& AGE, unsigned int minN=my_NAN);
	double getRs_ofPatch            (unsigned int p, const age_t& AGE){return getRs_ofPatch(get_vPatch(p), age_t2idx(AGE));}
	double getRs_ofLocus            (const unsigned int& l, const age_idx& AGE, unsigned int minN=my_NAN);
	double getRs_ofLocus            (unsigned int l, const age_t& AGE){return getRs_ofLocus(l, age_t2idx(AGE));}
	double getRs_ofPatch_andLocus   (const age_idx& AGE, TPatch* p, const unsigned int& l, unsigned int minN);
	double getRs_ofPatch_andLocus   (unsigned int i, const age_t& AGE){
		age_idx AGEidx = age_t2idx(AGE);
		unsigned int sample1, locus, minN = get_minSampleSize(AGEidx);
		fromID(i, sample1, locus, _nb_patch, _nb_locus);
		return getRs_ofPatch_andLocus(AGEidx, get_vPatch(sample1), locus, minN);
	}
    
	double getRt                    (const age_idx& AGE);
	double getRt                    (const age_t& AGE){return getRt(age_t2idx(AGE));}
	double getRt_ofLocus            (unsigned int i, const age_t& AGE) {return getRt_ofLocus(i, age_t2idx(AGE));}
	double getRt_ofLocus            (const unsigned int& l, const age_idx& AGE, unsigned int minN=my_NAN);
    
	// allelic range
	double getRange                    (const age_idx& AGE);
	double getRange                    (const age_t& AGE){return getRange(age_t2idx(AGE));}
	double getRange_ofPatch            (TPatch* p, const age_idx& AGE);
	double getRange_ofPatch            (unsigned int p, const age_t& AGE){return getRange_ofPatch(get_vPatch(p), age_t2idx(AGE));}
	double getRange_ofLocus            (const unsigned int& l, const age_idx& AGE);
	double getRange_ofLocus            (unsigned int l, const age_t& AGE){return getRange_ofLocus(l, age_t2idx(AGE));}
	unsigned int getRange_ofPatch_andLocus(const age_idx& AGE, TPatch* p, const unsigned int& l);
	double getRange_ofPatch_andLocus   (unsigned int i, const age_t& AGE){
		unsigned int sample, locus;
		fromID(i, sample, locus, _nb_patch, _nb_locus);
		return (double)getRange_ofPatch_andLocus(age_t2idx(AGE), get_vPatch(sample), locus);
	}
    
	double getRangeTot                 (const age_idx& AGE);
	double getRangeTot                 (const age_t& AGE){return getRangeTot(age_t2idx(AGE));}
	double getRangeTot_ofLocus         (unsigned int i, const age_t& AGE) {return getRangeTot_ofLocus(i, age_t2idx(AGE));}
	unsigned int getRangeTot_ofLocus   (const unsigned int& l, const age_idx& AGE);
    
	// garza williamson stat
	double getGW                    (const age_idx& AGE);
	double getGW                    (const age_t& AGE){return getGW(age_t2idx(AGE));}
	double getGW_ofPatch            (TPatch* p, const age_idx& AGE);
	double getGW_ofPatch            (unsigned int p, const age_t& AGE){return getGW_ofPatch(get_vPatch(p), age_t2idx(AGE));}
	double getGW_ofLocus            (const unsigned int& l, const age_idx& AGE);
	double getGW_ofLocus            (unsigned int l, const age_t& AGE){return getGW_ofLocus(l, age_t2idx(AGE));}
	double getGW_ofPatch_andLocus(const age_idx& AGE, TPatch* p, const unsigned int& l);
	double getGW_ofPatch_andLocus   (unsigned int i, const age_t& AGE){
		unsigned int sample, locus;
		fromID(i, sample, locus, _nb_patch, _nb_locus);
		return getGW_ofPatch_andLocus(age_t2idx(AGE), get_vPatch(sample), locus);
	}
    
	double getGWTot                 (const age_idx& AGE){return (double)getNbAlleleTot(AGE)/(1+getRangeTot(AGE));}
	double getGWTot                 (const age_t& AGE){return getGWTot(age_t2idx(AGE));}
	double getGWTot_ofLocus         (unsigned int i, const age_t& AGE) {return getGWTot_ofLocus(i, age_t2idx(AGE));}
	double getGWTot_ofLocus         (const unsigned int& l, const age_idx& AGE);
    
	// mean number of alleles per locus and patch/metapopulation
	double getNbAllele              (const age_t& AGE);
	double getNbAllele_ofPatch      (TPatch* p, const age_idx& AGE);
	double getNbAllele_ofPatch      (unsigned int i, const age_t& AGE);
	double getNbAllele_ofLocus      (unsigned int i, const age_t& AGE);
	unsigned int getNbAllele_ofPatch_andLocus(const age_idx& AGE, TPatch* p, const unsigned int& l);
	unsigned int getNbAllele_ofPatch_andLocus(const age_t& AGE, unsigned int p, const unsigned int& l){
		return getNbAllele_ofPatch_andLocus(age_t2idx(AGE), p, l);
	}
	double getNbAllele_ofPatch_andLocus   (unsigned int i, const age_t& AGE){
		unsigned int sample, locus;
		fromID(i, sample, locus, _nb_patch, _nb_locus);
		return (double)getNbAllele_ofPatch_andLocus(age_t2idx(AGE), get_vPatch(sample), locus);
	}
	double getNbAlleleTot           (const age_t& AGE){return getNbAlleleTot(age_t2idx(AGE));}
	double getNbAlleleTot           (const age_idx& AGE);
	unsigned int getNbAlleleTot_ofLocus(const age_idx & AGE, const unsigned int& l);
	double getNbAlleleTot_ofLocus   (unsigned int i, const age_t& AGE){
		return (double)getNbAlleleTot_ofLocus(age_t2idx(AGE), i);
	}
    
	// mean number of fixed loci
	double getNbFixedLocus              (const age_t& AGE);
	unsigned int getNbFixedLocus_ofPatch(TPatch* p, const age_idx& AGE);
	double getNbFixedLocus_ofPatch      (unsigned int i, const age_t& AGE);
	double getNbFixedLocus_ofLocus      (unsigned int i, const age_t& AGE);
	unsigned int getNbFixedLocus_ofPatch_andLocus(const age_idx& AGE, TPatch* p, const unsigned int& l);
	unsigned int getNbFixedLocus_ofPatch_andLocus(const age_t& AGE, TPatch* p, const unsigned int& l){
		return getNbFixedLocus_ofPatch_andLocus(age_t2idx(AGE), p, l);
	}
	double getNbFixedLocus_ofPatch_andLocus   (unsigned int i, const age_t& AGE){
		unsigned int sample, locus;
		fromID(i, sample, locus, _nb_patch, _nb_locus);
		return (double)getNbFixedLocus_ofPatch_andLocus(age_t2idx(AGE), get_vPatch(sample), locus);
	}
	double getNbFixedLocusTot           (const age_t& AGE){return getNbFixedLocusTot(age_t2idx(AGE));}
	unsigned int getNbFixedLocusTot     (const age_idx& AGE);
	unsigned int getNbFixedLocusTot_ofLocus(const age_idx & AGE, const unsigned int& l);
	double getNbFixedLocusTot_ofLocus   (unsigned int i, const age_t& AGE){
		return (double)getNbFixedLocusTot_ofLocus(age_t2idx(AGE), i);
	}
    
    
	double getHsnei(const age_t& AGE) {
		age_idx curAGE=age_t2idx(AGE);
		setFstat_Nei_Chesser(curAGE);
		return _hsnei[curAGE];
	}
	double getHtnei(const age_t& AGE) {
		age_idx curAGE=age_t2idx(AGE);
		setFstat_Nei_Chesser(curAGE);
		return _htnei[curAGE];
	}
    
	double getFst(const age_t& AGE) {
		age_idx curAGE=age_t2idx(AGE);
		setFstat_Nei_Chesser(curAGE);
		return _fst[curAGE];
	}
	double getFis(const age_t& AGE) {
		age_idx curAGE=age_t2idx(AGE);
		setFstat_Nei_Chesser(curAGE);
		return _fis[curAGE];
	}
	double getFit(const age_t& AGE) {
		age_idx curAGE=age_t2idx(AGE);
		setFstat_Nei_Chesser(curAGE);
		return _fit[curAGE];
	}
    
	double getHsnei_ofLocus(unsigned int i, const age_t& AGE) {
		age_idx curAge = age_t2idx(AGE);
		setFstat_Nei_Chesser_perLocus(curAge);
		return _hsnei_locus[curAge][i];
	}
	double getHtnei_ofLocus(unsigned int i, const age_t& AGE) {
		age_idx curAge = age_t2idx(AGE);
		setFstat_Nei_Chesser_perLocus(curAge);
		return _htnei_locus[curAge][i];
	}
	double getFst_ofLocus(unsigned int i, const age_t& AGE) {
		age_idx curAge = age_t2idx(AGE);
		setFstat_Nei_Chesser_perLocus(curAge);
		return _fst_locus[curAge][i];
	}
	double getFis_ofLocus(unsigned int i, const age_t& AGE) {
		age_idx curAge = age_t2idx(AGE);
		setFstat_Nei_Chesser_perLocus(curAge);
		return _fis_locus[curAge][i];
	}
	double getFit_ofLocus(unsigned int i, const age_t& AGE) {
		age_idx curAge = age_t2idx(AGE);
		setFstat_Nei_Chesser_perLocus(curAge);
		return _fit_locus[curAge][i];
	}
    
	/**Accessor to the Fst matrix as set by setFstMatrix().*/
	double getFst_ij(unsigned int i, const age_t& AGE){  // all pairs are computed
		age_idx curAge = age_t2idx(AGE);
		setFstat_Nei_Chesser_perPatchPair(curAge);
		unsigned int id1, id2;
		if(!id2sampleIDs(i, id1, id2))return my_NAN;	// if at least one patch is not sampled
		return _fst_matrix[curAge][id1][id2];
	}
	double getFst_ij_single(unsigned int i, const age_t& AGE){ // just this pair is computed
		unsigned int sample1, sample2;
		fromID(i, sample1, sample2, _nb_patch, _nb_patch);
		return getFstat_Nei_Chesser_perPatchPair(age_t2idx(AGE), get_vPatch(sample1), get_vPatch(sample2));
	}
	double getFst_ij_l_single(unsigned int i, const age_t& AGE){  // just this pair is computed
		unsigned int id1, id2, l;
		if(!id2sampleIDs(i, id1, id2, l))return my_NAN;	// if at least one patch is not sampled
		return getFstat_Nei_Chesser_perPatchPair_andLocus(age_t2idx(AGE), get_vPatch(id1), get_vPatch(id2), l);
	}
    
    // fdist (Beaumont & Nichols, 1996)
	double getHet0_ofLocus(unsigned int i, const age_t& AGE) {
		age_idx curAge = age_t2idx(AGE);
		setFstat_Beaumont_Nichols_perLocus(curAge);
		return _het0_locus[curAge][i];
	}
	double getHet1_ofLocus(unsigned int i, const age_t& AGE) {
		age_idx curAge = age_t2idx(AGE);
		setFstat_Beaumont_Nichols_perLocus(curAge);
		return _het1_locus[curAge][i];
	}
	double getFst_fdist_ofLocus(unsigned int i, const age_t& AGE) {
		age_idx curAge = age_t2idx(AGE);
		setFstat_Beaumont_Nichols_perLocus(curAge);
		return _fst_bn_locus[curAge][i];
	}
	void setFstat_Beaumont_Nichols_perLocus(const age_idx& AGE);
    
	// linkage disequilibrium
	double getDprime                (const age_idx& AGE, const unsigned int& l1, const unsigned int& l2);
	double getDprime_ofPatch        (const age_idx& AGE, TPatch* p, const unsigned int& l1, const unsigned int& l2);
	double getDstar                 (const age_idx& AGE, const unsigned int& l1, const unsigned int& l2);
	double getDstar_ofPatch         (const age_idx& AGE, TPatch* p, const unsigned int& l1, const unsigned int& l2);
	double getR2                    (const age_idx& AGE, const unsigned int& l1, const unsigned int& l2);
	double getR2_ofPatch            (const age_idx& AGE, TPatch* p, const unsigned int& l1, const unsigned int& l2);
	double getChi2                  (const age_idx& AGE, const unsigned int& l1, const unsigned int& l2);
	double getChi2_ofPatch          (const age_idx& AGE, TPatch* p, const unsigned int& l1, const unsigned int& l2);
	double getDprime_ij             (unsigned int i, const age_t& AGE){
		unsigned int loc1, loc2;
		fromID(i, loc1, loc2, _nb_locus, _nb_locus);
		return getDprime(age_t2idx(AGE), loc1, loc2);
	}
	double getDstar_ij(unsigned int i, const age_t& AGE){
		unsigned int loc1, loc2;
		fromID(i, loc1, loc2, _nb_locus, _nb_locus);
		return getDstar (age_t2idx(AGE), loc1, loc2);
	}
	double getR2_ij(unsigned int i, const age_t& AGE){
		unsigned int loc1, loc2;
		fromID(i, loc1, loc2, _nb_locus, _nb_locus);
		return getR2    (age_t2idx(AGE), loc1, loc2);
	}
	double getChi2_ij(unsigned int i, const age_t& AGE){
		unsigned int loc1, loc2;
		fromID(i, loc1, loc2, _nb_locus, _nb_locus);
		return getChi2  (age_t2idx(AGE), loc1, loc2);
	}
	double getDprime_ij_ofPatch(unsigned int i, const age_t& AGE){
		unsigned int p,l1, l2;
		fromID(i, p, l1, l2, _nb_patch, _nb_locus, _nb_locus);
		return getDprime_ofPatch(age_t2idx(AGE), get_vPatch(p), l1, l2);
	}
	double getDstar_ij_ofPatch(unsigned int i, const age_t& AGE){
		unsigned int p,l1, l2;
		fromID(i, p, l1, l2, _nb_patch, _nb_locus, _nb_locus);
		return getDstar_ofPatch(age_t2idx(AGE), get_vPatch(p), l1, l2);
	}
	double getR2_ij_ofPatch(unsigned int i, const age_t& AGE){
		unsigned int p,l1, l2;
		fromID(i, p, l1, l2, _nb_patch, _nb_locus, _nb_locus);
		return getR2_ofPatch(age_t2idx(AGE), get_vPatch(p), l1, l2);
	}
	double getChi2_ij_ofPatch(unsigned int i, const age_t& AGE){
		unsigned int p,l1, l2;
		fromID(i, p, l1, l2, _nb_patch, _nb_locus, _nb_locus);
		return getChi2_ofPatch(age_t2idx(AGE), get_vPatch(p), l1, l2);
	}
    
	// get the total population size of the sampled patches
	unsigned int get_nbInds(const age_idx& AGE, const sex_t& SEX);
	unsigned int get_nbInds(const age_idx& AGE){return get_nbInds(AGE,FEM)+get_nbInds(AGE,MAL);}
	unsigned int get_nbInds(const age_t& AGE, const sex_t& SEX);
	unsigned int get_nbInds(const age_t& AGE){return get_nbInds(AGE,FEM)+get_nbInds(AGE,MAL);}
    
	// get the mean population size of the sampled patches
	double get_meanNbInds(const age_idx& AGE, const sex_t& SEX);
	double get_meanNbInds(const age_idx& AGE);
	double get_meanNbInds(const age_t& AGE, const sex_t& SEX);
	double get_meanNbInds(const age_t& AGE);
    
	// get the total sampled size of the sampled patches
	unsigned int get_nbSamples(const age_idx& AGE, const sex_t& SEX);
	unsigned int get_nbSamples(const age_idx& AGE){return get_nbSamples(AGE,FEM)+get_nbSamples(AGE,MAL);}
	unsigned int get_nbSamples(const age_t& AGE, const sex_t& SEX);
	unsigned int get_nbSamples(const age_t& AGE){return get_nbSamples(AGE,FEM)+get_nbSamples(AGE,MAL);}
    
	// get the mean population size of the sampled patches
	double get_meanNbSamples(const age_idx& AGE, const sex_t& SEX);
	double get_meanNbSamples(const age_idx& AGE);
	double get_meanNbSamples(const age_t& AGE, const sex_t& SEX);
	double get_meanNbSamples(const age_t& AGE);
    
	// chord distance following Cavalli_Sforza & Bodmer (1971)
	double getChordDist2_perPatchPair(const age_idx& AGE, TPatch* pop1, TPatch* pop2);
	double getChordDist_perPatchPair(const age_idx& AGE, TPatch* pop1, TPatch* pop2);
	double getChordDist(const age_t& AGE);   // average chord distance
	double getChordDist_ij(unsigned int i, const age_t& AGE){  // just this pair is computed
		unsigned int sample1, sample2;
		fromID(i, sample1, sample2, _nb_patch, _nb_patch);
		return getChordDist_perPatchPair(age_t2idx(AGE), get_vPatch(sample1), get_vPatch(sample2));
	}
};

#endif //STATHANDLER_H

