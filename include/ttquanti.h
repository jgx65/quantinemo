/** @file ttquanti.h
 *
 *   Copyright (C) 2006 Frederic Guillaume    <guillaum@zoology.ubc.ca>
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

#ifndef ttquantiH
#define ttquantiH

#include "ttrait.h"
#include "filehandler.h"
#include "stathandler.h"
#include "tree.h"

class TTQuantiFHvalue;
class TTraitQuantiProto;

/******************************************************************************/
/******************************************************************************/
/**Quantitative trait coded on discrete allelic values.
 * Alleles coded on \<char\>, max number of alleles per locus is 256.*/
class TTraitQuanti : public TTrait{
	friend class TTraitQuantiProto; // we allow to access these parameters from TTraitQuantiProto directly
private:
	double _genotype;             // this is the genotype (Vg)
	double _phenotype;            // this is the phenotype (Vp = Vg + Ve)
	double _fitness;              // this is the fitness component
	double _fitness_factor;       // by default 1

	TTraitQuantiProto* pProto;

public:

	TTraitQuanti (): _phenotype(my_NAN), _fitness(my_NAN), _fitness_factor(my_NAN), pProto(0){}

    TTraitQuanti(const TTraitQuanti& T);
    
    virtual ~TTraitQuanti ();

	virtual void set_from_prototype(TTraitProto* T);

	///@}
	///@name Implementations
	///@{
	virtual TTraitQuanti& operator= (const TTrait& T);
	virtual bool    operator== (const TTrait& T);
	virtual bool    operator!= (const TTrait& T);
	virtual void    reset                ();
	virtual void*   set_trait            (void* value)           {return NULL;}
	virtual void**  get_sequence         ()  const              {return (void**)sequence;}
	virtual void    set_sequence         (void** seq){
		reset();
		sequence = (unsigned char**)seq;
	}
	virtual void    set_value            (); // set the genotype
	virtual void    set_value            (double val)          {_phenotype = val+_genotype;}   // set the phenotype
	virtual double  get_value            ()					   {return _phenotype;}
	virtual double  get_genotype         ()					   {return _genotype;}
	virtual double  get_phenotype        ()                    {return _phenotype;}
	virtual double  get_fitness          ()                    {return _fitness;}
	virtual double  get_fitnessFactor    ()                    {return _fitness_factor;}
	virtual void    set_fitness          (double value)        {_fitness = value;}
	virtual void    set_fitness_factor   () ;
	virtual double  set_get_fitness_factor() ;


	virtual void    show_up              ();
	virtual TTraitQuanti*  clone     ()                      {return new TTraitQuanti(*this);}

	///@}
};



class TTQuantiSH;
class TTQuantiFH;
class TTQuantiFHvalue;
/******************************************************************************/
/******************************************************************************/
/**Prototype class for the TTQuanti trait class.
 **/
class TTraitQuantiProto : public TTraitProto {
	friend class TTraitQuanti; // we allow to access these parameters from TTQuanti directly
private:
	int     _output;              //0: nothing, 1: allelic values
//	int     _Ve_model;            // environmental model: 0: Ve directly, 1: h2; 2: h2 at every gen; 3: H2; 4: H" at every gen
//	double  _Ve_prop;             // proportion of the environment from the current versus the natal patch
	unsigned int     _selection_model;     // 0: neutral; 1: stabilizing, 2: directional 3: fitness landscape, 4: selection coeffcient
	int     _Va_model;            // 0: always valid, but slow, 1: limited to random mating, but fast

public:
    double getAllelicValue(const unsigned int& l, const unsigned char& i);
    double getDominanceValue(const unsigned int& l, const unsigned char& a1, const unsigned char& a2);

	TTQuantiSH*             _stats;
	TTQuantiFH*      _writer;        // only one instance per type of trait: the first trait (index 1) owns it
	TTQuantiFHvalue* _phenotyper;    // only one instance per type of trait: the first trait (index 1) owns it
    TTQuantiFHvalue* _genotyper;     // only one instance per type of trait: the first trait (index 1) owns it

	double**             _allelicValues;      // _allelicValues[locus][allele]
	bool                 _allelic_file;       // is a alleic file passed?

	double***            _dominanceValues;    // only for dominance effect    _dominanceValues[locus][allele1][allele2]
	double               _dominance_mean;
	double               _dominance_sd;       // 0 if not used!
	bool                 _dominance_file;     // is a dominance file passed?

	double*              _fitnessFactor_heterozygote;   // _fitnessFactor_heterozygote[locus]
	double*              _fitnessFactor_homozygote;     // _fitnessFactor_homozygote[locus]
	double*              _fitnessFactor_freqDepend;     // _fitnessFactor_freqDepend[locus]
	double***            _fitnessFactor;                // _fitnessFactor[locus][allele1][allele2]
	map<unsigned char, map< unsigned char, double> >*  _locusFreqs; // _locusFreqs[locus][allele1][allele2] used for _fitnessFactor_freqDepend
	Tree<unsigned char>* _fitnessFactorTree;            //  fitness factor defined for the entire genome

	Tree<unsigned char>* _phenoTree;          // only for epistatic effect
	double               _epistatic_sd;

	// determination of the genotype
	double  (TTraitQuantiProto::* get_locus_genotype_func_ptr)(const unsigned int& l, const unsigned char& a1, const unsigned char& a2);
	double  get_locus_genotype_additive (const unsigned int& l, const unsigned char& a1, const unsigned char& a2);
	double  get_locus_genotype_dominance_single (const unsigned int& l, const unsigned char& a1, const unsigned char& a2);
	double  get_locus_genotype_dominance_array (const unsigned int& l, const unsigned char& a1, const unsigned char& a2);

	double (TTraitQuantiProto::* get_genotype_dominace_func_ptr)(double a1, double a2, double dom);
	double get_genotype_dominance(double a1, double a2, double h){
		return (this->*get_genotype_dominace_func_ptr)(a1, a2, h);
	}
	double get_genotype_dominance_h(double a1, double a2, double h);
	double get_genotype_dominance_k(double a1, double a2, double k);

	double get_dominance_mean(){return _dominance_mean;}

	double  (TTraitQuantiProto::* get_genotype_func_ptr)(unsigned char** seq);
	double  get_genotype_additive (unsigned char** seq);
	double  get_genotype_epistatic (unsigned char** seq);

	double* get_fitnessFactor_heterozygote(){return _fitnessFactor_heterozygote;}
	double* get_fitnessFactor_homozygote()  {return _fitnessFactor_homozygote;}
	double*** get_fitnessFactor_array()  {return _fitnessFactor;}
	bool    fitnessFactor_used() {return (get_fitnessFactor_func_ptr != NULL || get_fitnessFactor2_func_ptr != NULL);}
	bool    fitnessFactor_freqDep_used() {return (get_fitnessFactor2_func_ptr != NULL);}
	double  (TTraitQuantiProto::* get_fitnessFactor_func_ptr)(unsigned char** seq);
	double  (TTraitQuantiProto::* get_fitnessFactor2_func_ptr)(unsigned char** seq);
	double  get_fitnessFactor_genome(unsigned char** seq);
	double  get_fitnessFactor_locus(unsigned char** seq);
	double  get_fitnessFactor_global(unsigned char** seq);
    double  get_fitnessFactor_freqDepend(unsigned char** seq);
    double* get_fitnessFactor_freqDepend(){return _fitnessFactor_freqDepend;}
    
    map<unsigned char, map< unsigned char, double> >*& get_locusFreqs() {return _locusFreqs;}


	void    check_allelicValues_IMM();
	void    check_mutationValues_IMM();

	void    ini_paramset             ( );   // before the normal ini()

public:

	TTraitQuantiProto ( );
	TTraitQuantiProto (int i);
	TTraitQuantiProto(const TTraitQuantiProto& T);

	~TTraitQuantiProto ( );

	//implementation of TTraitProto:
	virtual void                     init (TMetapop* pMetapop);

	virtual void                      reset ();
	virtual void                      resetTotal ();
	void                              delete_dominanceValues();
	void                              delete_allelicValues();
	void                              delete_fitnessFactor();
    void                              delete_locusFreqs();

	virtual TTraitQuanti*          hatch ();

	virtual TTraitQuantiProto*      clone () {return new TTraitQuantiProto(*this);}

	//implementation of SimComponent:
	virtual void loadFileServices ( FileServices* loader );

	virtual void loadStatServices ( StatServices* loader );

	void read_allele_file(string name);
	void set_fitnessFactor(const string& trait, string name, double* &array);
	void set_allelicValues(const string& trait);
	void set_allelicValues_locus(const unsigned int& l, const double& var, const string& trait);
	void set_allelicValues(TMatrix* m, const unsigned int& i, const unsigned int& l, const unsigned int& a, unsigned int* cols);
	void set_dominanceValues(const string& trait);
	void set_epistaticValues(const string& traitS);
	void set_mutationFreq(const string& trait);
	void set_mutationFreq_locus(const unsigned int& l, const double& var);
	void set_initAlleleFreq(const string& trait);
	void set_initAlleleFreq_locus(const unsigned int& l, const double& var);
	void read_locus_file(string name);
	void read_genome_file(string name);

	void print_allelic_values (string name);
	void print_dominance_values (string name);
	void print_epistatic_values (string name);
	void print_gentoype(ostream& FILE, unsigned char** seq, const unsigned int& digit);
	bool get_next_gentoype(unsigned char** seq);

	virtual void temporal_change(const unsigned int& gen);
	virtual void executeAfterEachReplicate(const unsigned int& rep);

	unsigned int get_selection_model()     const     {return _selection_model;}
	int     get_Va_model()            const     {return _Va_model;}

	string  get_info             ();

	void    create_regular_spaced_array(double* array, const int& size, double half_range);
	void    compute_frequencies(double* array, double* effect_array, const unsigned int& size, const double& sd);

};




/******************************************************************************/
/******************************************************************************/
/** A file handler to save the discrete quantitative trait genotypes in a FSTAT-like text file.
 *  The file extension is ".dat".
 */
class TTQuantiFH: public FileHandler {
public:

	TTQuantiFH (){}

	virtual ~TTQuantiFH ( ) { }

};


/******************************************************************************/
/******************************************************************************/
/** File handler used to save the phenotypes 0por genotypes generated from the
 * discrete quantitative trait. The file extension is
 * ".phe". !!!
 */
class TTQuantiFHvalue: public FileHandler {

private:
    virtual void FHwrite (const age_idx& cur_age, const sex_t& cur_sex, ostream& FILE,
                          TPatch* current_patch, const int& patch_id);
    
    double (TTQuantiFHvalue::*get_value_func_ptr)(TIndividual* ind, const int& t);
    double get_genotype(TIndividual* ind, const int& t);
    double get_phenotype(TIndividual* ind, const int& t);
    

public:

TTQuantiFHvalue(){}

    virtual  void set_getter(int i){
        switch(i){
            case 0: get_value_func_ptr= &TTQuantiFHvalue::get_genotype;
                set_name("genotypic value");
                set_extension(".gen");
                break;
            case 1: get_value_func_ptr= &TTQuantiFHvalue::get_phenotype;
                set_name("phenotype");
                set_extension(".phe");
                break;
            default: error("TTQuantiFHvalue::set_getter: Not valid option for getter!\n");
        }
    }
    
virtual ~TTQuantiFHvalue ( ) { }

virtual void FHwrite  ();
};

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                               ****** TTQuantiSH ******

/*_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\*/
/**The stat handler for neutral markers. */
class TTQuantiSH: public StatHandler<TTQuantiSH> {
private:


	double *_varA, *_varA2; //, _h2;

	double 	**_meanG, **_varG, **_meanP, **_varP, **_meanW, **_varW;		// 0: MAL; 1: FEM; 2: both


	double** _qst_matrix;        // Qst matrix for pairwise values
	double** _qstF_matrix;       // Qst matrix for pairwise values corrected for inbreeding

//	double Theta_FF, Theta_MM, Theta_FM;
//	double _mean_theta, _mean_alpha;

public:

	TTQuantiSH (TTraitQuantiProto* TT) :  _varA(0), _varA2(0),/* _h2(0),*/ _meanG(0), _varG(0),
	_meanP(0), _varP(0), _meanW(0), _varW(0), _qst_matrix(0), _qstF_matrix(0)
{
		set(TT);

		// how the additive genetice variance is computed
		switch(TT->get_Va_model()){
		case 0: get_Va_ofPatch_func_ptr = &TTQuantiSH::get_Va_ofPatch_regression;    break; // any case
		case 1: get_Va_ofPatch_func_ptr = &TTQuantiSH::get_Va_ofPatch_random_mating; break; // random mating
		case 2: get_Va_ofPatch_func_ptr = &TTQuantiSH::setMeanAndVar_Vg_ofPatch;     break; // Va = Vg
		}
}

	virtual ~TTQuantiSH ( )
	{
		if(_varA)     delete[] _varA;
		if(_varA2)    delete[] _varA2;

		ARRAY::delete_2D(_meanG, 3);
		ARRAY::delete_2D(_varG, 3);
		ARRAY::delete_2D(_meanP, 3);
		ARRAY::delete_2D(_varP, 3);
		ARRAY::delete_2D(_meanW, 3);
		ARRAY::delete_2D(_varW, 3);
        assert(_popPtr);
		ARRAY::delete_2D(_qst_matrix, get_current_nbSamplePatch());
		ARRAY::delete_2D(_qstF_matrix, get_current_nbSamplePatch());

	}

	virtual bool init ( ) ;

	virtual bool setStatRecorders (const string& token);
	virtual string getName() {return "QuantiSH";}


	// variance components
	bool compute_alpha(double* y, const map<unsigned char, int*>& x, const unsigned int& nb_ind,
			map<unsigned char, double>& alpha, const map<unsigned char, double>& availableAllele);
	bool remove_private_alleles_compute_alpha(TPatch* crnt_patch, const unsigned int& sizeF,
			const unsigned int& sizeM, map<unsigned char, double>& alpha,
			double* arrayG, const age_idx& age_pos,
			map<unsigned char, double>& allele_freq, const unsigned int& l);
	void (TTQuantiSH::*get_Va_ofPatch_func_ptr)(TPatch*, const age_idx&, double&, double&, map<unsigned char, double>*);
	void get_Va_ofPatch(TPatch* p, const age_idx& a, double& m, double& v, map<unsigned char, double>* f){(this->*get_Va_ofPatch_func_ptr)(p,a,m,v,f);}
	void get_Va_ofPatch_regression(TPatch*, const age_idx&, double&, double&, map<unsigned char, double>*);    // any case, but slower
	void get_Va_ofPatch_random_mating(TPatch*, const age_idx&, double&, double&, map<unsigned char, double>*); // only random mating

	double getVarA      (unsigned int i, const age_t& AGE){return getVarA(i, age_t2idx(AGE));}
	double getVarA      (unsigned int i, const age_idx& AGE)  {
		setVar_Va(AGE);
		return _varA[i];
	}

	// genotypic value
	double getMeanG     (unsigned int i, const age_t& AGE){return getMeanG(i, age_t2idx(AGE));}
	double getMeanG     (unsigned int i, const age_idx& AGE)  {
		setMeanAndVar_Vg(AGE);
		return _meanG[2][i];
	}
	double getMeanGfem  (unsigned int i, const age_t& AGE){return getMeanGfem(i, age_t2idx(AGE));}
	double getMeanGfem  (unsigned int i, const age_idx& AGE)  {
		setMeanAndVar_Vg(AGE, FEM);
		return _meanG[FEM][i];
	}
	double getMeanGmal  (unsigned int i, const age_t& AGE){return getMeanGmal(i, age_t2idx(AGE));}
	double getMeanGmal  (unsigned int i, const age_idx& AGE)  {
		setMeanAndVar_Vg(AGE, MAL);
		return _meanG[MAL][i];
	}
	double getVarG      (unsigned int i, const age_t& AGE){return getVarG(i, age_t2idx(AGE));}
	double getVarG      (unsigned int i, const age_idx& AGE)  {
		setMeanAndVar_Vg(AGE);
		return _varG[2][i];
	}
	double getVarGfem   (unsigned int i, const age_t& AGE){return getVarGfem(i, age_t2idx(AGE));}
	double getVarGfem   (unsigned int i, const age_idx& AGE)  {
		setMeanAndVar_Vg(AGE, FEM);
		return _varG[FEM][i];
	}
	double getVarGmal   (unsigned int i, const age_t& AGE){return getVarGmal(i, age_t2idx(AGE));}
	double getVarGmal   (unsigned int i, const age_idx& AGE)  {
		setMeanAndVar_Vg(AGE, MAL);
		return _varG[MAL][i];
	}

	// phenotypic value
	double getMeanP     (unsigned int i)  {
		setMeanAndVar_Vp();
		return _meanP[2][i];
	}
	double getMeanPfem  (unsigned int i)  {
		setMeanAndVar_Vp(FEM);
		return _meanP[FEM][i];
	}
	double getMeanPmal  (unsigned int i)  {
		setMeanAndVar_Vp(MAL);
		return _meanP[MAL][i];
	}
	double getVarP      (unsigned int i)  {
		setMeanAndVar_Vp();
		return _varP[2][i];
	}
	double getVarPfem   (unsigned int i)  {
		setMeanAndVar_Vp(FEM);
		return _varP[FEM][i];
	}
	double getVarPmal   (unsigned int i)  {
		setMeanAndVar_Vp(MAL);
		return _varP[MAL][i];
	}

	// fitness value
	double getMeanW     (unsigned int i)  {
		setMeanAndVar_Wp();
		return _meanW[2][i];
	}
	double getMeanWfem  (unsigned int i)  {
		setMeanAndVar_Wp(FEM);
		return _meanW[FEM][i];
	}
	double getMeanWmal  (unsigned int i)  {
		setMeanAndVar_Wp(MAL);
		return _meanW[MAL][i];
	}
	double getVarW      (unsigned int i)  {
		setMeanAndVar_Wp();
		return _varW[2][i];
	}
	double getVarWfem   (unsigned int i)  {
		setMeanAndVar_Wp(FEM);
		return _varW[FEM][i];
	}
	double getVarWmal   (unsigned int i)  {
		setMeanAndVar_Wp(MAL);
		return _varW[MAL][i];
	}

	// QST
	double getQst       (const age_t& AGE){return getQst(age_t2idx(AGE));};
	double getQst       (const age_idx& AGE);
	void setQst_perPatchPair(const age_idx& AGE);
	double getQst_ij    (unsigned int i, const age_t& AGE){return getQst_ij(i, age_t2idx(AGE));};
	double getQst_ij    (unsigned int i, const age_idx& AGE){
		setQst_perPatchPair(AGE);
		unsigned int id1, id2;
		if(!id2sampleIDs(i, id1, id2))return my_NAN;	// if at least one patch is not sampled
		return _qst_matrix[id1][id2];
	}

	// QST corrected for inbreeding
	double getQstF      (const age_t& AGE){return getQstF(age_t2idx(AGE));};
	double getQstF      (const age_idx& AGE);
	void setQstF_perPatchPair(const age_idx& AGE);
	double getQstF_ij   (unsigned int i, const age_t& AGE){return getQstF_ij(i, age_t2idx(AGE));};
	double getQstF_ij   (unsigned int i, const age_idx& AGE){
		setQstF_perPatchPair(AGE);
		unsigned int id1, id2;
		if(!id2sampleIDs(i, id1, id2))return my_NAN;	// if at least one patch is not sampled
		return _qstF_matrix[id1][id2];
	}


	// between variance
	double getVgB       (const age_t& AGE){return getVgB(age_t2idx(AGE));};
	double getVgB       (const age_idx& AGE);
	double getVpB       ();

	// within variance
	double getVaW       (const age_t& AGE){return getVaW(age_t2idx(AGE));};
	double getVaW       (const age_idx& AGE);
	double getVgW       (const age_t& AGE){return getVgW(age_t2idx(AGE));};
	double getVgW       (const age_idx& AGE);
	double getVpW       ();

	void   setVar_Va(const age_idx& AGE);
	void   setMeanAndVar_Vg(const age_idx& AGE);
	void   setMeanAndVar_Vg(const age_idx& AGE, sex_t SEX);
	void   setMeanAndVar_Vg_ofPatch_allInds(TPatch*, const age_idx&, double&, double&, map<unsigned char, double>* freqs=NULL);
	void   setMeanAndVar_Vg_ofPatch(TPatch*, const age_idx&, double&, double&, map<unsigned char, double>* freqs=NULL);
	void   setMeanAndVar_Vg_ofPatch(TPatch*, const age_idx&, double&, double&, sex_t SEX, map<unsigned char, double>* freqs=NULL);
	void   setMeanAndVar_Vp();
	void   setMeanAndVar_Vp(sex_t SEX);
	void   setMeanAndVar_Vp_ofPatch(TPatch* crnt_patch, double& meanP, double& varP);
	void   setMeanAndVar_Vp_ofPatch(TPatch* crnt_patch, double& meanP, double& varP, sex_t SEX);
	void   setMeanAndVar_Wp();
	void   setMeanAndVar_Wp(sex_t SEX);
	void   setMeanAndVar_Wp_ofPatch(TPatch* crnt_patch, double& meanP, double& varP);
	void   setMeanAndVar_Wp_ofPatch(TPatch* crnt_patch, double& meanP, double& varP, sex_t SEX);
};

#endif //TTDISCRETEQUANTI_H

