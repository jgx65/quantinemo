/** @file geneticmap.h
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

#ifndef tgenomeH
#define tgenomeH

#include "tlocus.h"
#include "simcomponent.h"

class TIndividual;
class TGenomeProto;
class TPatch;
class FileServices;
class StatServices;


//---------------------------------------------------------------------------
/** TGenome: object initialized for each individual
	* contains the sequence of the individual
	*/

class TGenome{
private:
	TGenomeProto* _protoGenome; // pointer to the prototype (don't delete it)
	unsigned char** sequence;  // sequence[locus][allele]; allele 0: from mother; allele 1: from father

public:
	TGenome();
	TGenome(const TGenome& T): _protoGenome(T._protoGenome), sequence(0){}
	~TGenome();
	void clear();
	void inherit(TIndividual* mother, TIndividual* father);
	void mutate();
	void ini_sequence(TPatch* p);
	void ini_sequence(unsigned char** seq);
	void ini_sequence_dadFirst(unsigned char** seq);
	void ini_sequence(unsigned char* seq_mum, unsigned char* seq_dad);
	void set_protoGenome(TGenomeProto* p){_protoGenome = p;}
	void create_sequence();
	void clone(const TGenome& gen);

	TGenomeProto* get_protoGenome() const {return _protoGenome;}
	unsigned char** get_sequence() const {return sequence;}

	TGenome& operator=(const TGenome& g);
};


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
/** TGenomeProto: a single object is initialized (contained by the metapop)
	* contains the super chromosome with the positions of all loci
	*/

class TMetapop;

class TGenomeProto: public SimComponent {
friend class TGenome;

protected:
	// super chromosome parameters (do not change over time, but may differ among female and male)
	// first the genetic map is set up for females and then if necessary the position for males is checked
	// note that the order of the loci on the genome must be identical for males and females
	// only the genetic distance between loci may differ
	double            _chromosomeLength[2];  // _chromosomeLength[SEX]: total length of the super chromosome
	double*           _locus_position_tot[2];// _locus_position_tot[SEX]: super chromosome with all loci positions
	unsigned int      _nb_chromosome;        // number of chromosomes (_chromosomeSizeVector.size(), unlinked loci are also on a chromosome!)

	TLocus**          _locus_tot;            // _locus_tot[locus]: locus array containing a pointer to each locus (first linked loci, followed by unlinked loci)
	unsigned int      _nb_locus_tot;		 // total number of loci
	unsigned int      _nb_locus_linked;      // number of linked loci
	unsigned int      _nb_locus_unlinked;	 // number of unlinked loci
	unsigned int*     _nb_locus_per_chromosome;// cumulative number of loci per chromosome _nb_locus_per_chromosome[chrom]

	// the actual chromosome to use (allows to use a recombination factor)
	double            _chromosomeLength_temp[2];  	// total length of the super chromosome
	double*           _locus_position_tot_temp[2];	// super chromosome with all loci positions (only delete it when recombination factor is used)
	double* 		  _recombination_factor_chrom[2]; // used if quantitative traits define the recombination factors
    
    // parameters just used to set up the super chromosome
    double*           _chromosomeSize[2];    // _chromosomeSize[SEX]: array of the cumulative length of the chromosomes (only used to set up the super chromosome);
	map<string, vector<TLocusPos*> >                   _locus_pos;     // _locus_pos[type][TLocusPos*]
	map<string, map<unsigned int, vector<TLocus*> > >  _locus_vector;  // _locus_vector[type][locus_id][TLocus*]

	// help arrays to simulate pleiotrophy
	// Idea: Each locus of each trait are listed in the array "_locus_tot" sorted by their genome position.
	//       The same is true for pleiotropic loci. If pleiotropic loci are simulated some additional parameters
	//       are needed to perform the mutation as efficient as possible. Note, that the
	//       parameter "_locus_tot" is not modified. The parameter "_locus_pleiotropic" somehow replaces the
	//       parameter "_locus_tot" and points to the first locus within the parameter "_locus_tot" of a pleiotropic
	//       set of loci. This new parameter is used for mutation functions if pleiotropic loci
	//       are present.
    //       Note, that pleiotrophy has no effect on inheritance as they behave the same as to loci at the same location!
	unsigned int*   _locus_pleiotropic;              // an array of size _locus_pleiotropic+1
                                                     // contains the starting position of each pleiotropic loci set
                                                     // the last element contains the first position beyond the last locus
	unsigned int     _nb_locus_pleiotropic;          // size of the _locus_pleiotropic array -1
    
    ////////////////////
	bool 			  _sex_specific_genome;
	double 			  _mut_rate_mean;				 // a mean mutation rate (used if the mutation rate is the same for all loci)
    
	bool _change_mutationRates;                      // flag to know when ini_mutate() has to be called
    
    
	virtual void loadFileServices ( FileServices* loader ) {}
	virtual void loadStatServices ( StatServices* loader ) {}

	// inheritance functions
	void (TGenomeProto::*_inherit_func_ptr)(TIndividual* mother, TIndividual* father, unsigned char** child);
	void _inherit_linked  (TIndividual* mother, TIndividual* father, unsigned char** child);
	void _inherit_unlinked(TIndividual* mother, TIndividual* father, unsigned char** child);
	void _inherit_mixed   (TIndividual* mother, TIndividual* father, unsigned char** child);
    
	// recombination factor
	void ini_recombination_factor();
	bool ini_recombination_factor(string param_name, sex_t SEX);
	void ini_recombination_qtrait(string param_name, sex_t SEX, double* vec, unsigned int size, TMatrix* m);
	void (TGenomeProto::*_recombine_func_ptr[2])(TIndividual* parent, unsigned char** child, int index);
	void _recombine_normal   (TIndividual* parent, unsigned char** child, int index = 0);
	void _recombine_qtrait   (TIndividual* parent, unsigned char** child, int index = 0);
    
	typedef double (TGenomeProto::*_func_ptr)(TIndividual* parent, sex_t SEX, unsigned int chrom);
	_func_ptr* _recombination_chrom_func_ptr[2];
	double _recombination_chrom_factor(TIndividual* parent, sex_t SEX, unsigned int chrom);
	double _recombination_chrom_qtraitZ(TIndividual* parent, sex_t SEX, unsigned int chrom);
	double _recombination_chrom_qtraitG(TIndividual* parent, sex_t SEX, unsigned int chrom);
    
    
	// mutation functions
	void  (TGenomeProto::*_mutate_func_ptr)(unsigned char** seq);
	void  _mutate_zero_mutation_rate   (unsigned char** seq) { }
	void  _mutate_equal_mutation_rate  (unsigned char** seq);
	void  _mutate_equal_mutation_rate_pleiotrophy(unsigned char** seq);
	void  _mutate_equal_mutation_rate_pleiotrophy_correl(unsigned char** seq);
	void  _mutate_unequal_mutation_rate(unsigned char** seq);
	void  _mutate_unequal_mutation_rate_pleiotrophy(unsigned char** seq);
	void  _mutate_unequal_mutation_rate_pleiotrophy_correl(unsigned char** seq);
    
	void ini();
	void ini_paramset();
    

public:
    
    static unsigned int MUTATION_TRIAL;        // the number of trials for the RMM model before giving up mutation
	TMetapop* _popPtr;

	void set_genome_positions(string type, TTraitProto* pTrait);
	bool set_genome_positions_fix();                                // global (genome)
	bool set_genome_positions_fix(string type, TTraitProto* pTrait);// per trait (ntrl_gneome)
	bool set_genome_positions_none(string type);
	void set_genome_positions_matrix(TMatrixVar<double>* matrix, sex_t SEX, string type, bool delMatrix=true);

	TGenomeProto();
	~TGenomeProto();

	void ini_all();

	void set_metapop_ptr(TMetapop* p){_popPtr = p;}

	TMatrixVar<double>* drawGeneticMapRandom(TMatrix* matrix, unsigned int& nbLocus);

	void inherit(TIndividual* mother, TIndividual* father, unsigned char** child);

	void mutate(unsigned char** seq);
	void ini_mutate();
	void set_mutation_of_locus(const unsigned int& l, const double& rate, const mut_model_t& model);

	void ini_sequence(unsigned char** seq, TPatch* patch);
	void set_ini_sequence_model(TLocus* aLocus, const unsigned int& size, const unsigned int& model);
	void set_ini_sequence_model(TLocus* aLocus, const unsigned int& size, const ini_model_t& model);

	/** creates the super chromosome based on the single chromosome */
	void ini_genetic_map();
	void printGeneticMapInfo();
	unsigned int getNbChromosome(){return _nb_chromosome;}

	void set_protoGenomePtr(TGenome* p){p->set_protoGenome(this);}

	void print_genome();
	void print_genome(ofstream& FILE, const string& type);
	string print_genome_sex(const string& type, sex_t SEX=FEM);
	void print_locus_index_tot(ofstream& FILE, TTraitProto* pTrait);
	void print_locus_index_trait(ofstream& FILE, TTraitProto* pTrait);
    
    void print_genome_tot(ofstream& FILE);
    string print_genome_tot_sex(sex_t SEX);
    void print_genome_trait(ofstream& FILE, const string& type);
    string print_genome_trait_sex(sex_t SEX, const string& type);
    

	void executeBeforeEachReplicate(const unsigned int& rep){}
	void executeAfterEachReplicate(const unsigned int& rep){}

	void executeBeforeEachGeneration(const unsigned int& gen){}
	void executeAfterEachGeneration(const unsigned int& gen){}

	void temporal_change(const unsigned int& gen);

	void clear();
    
    bool is_sexSpecificMap(){return _locus_position_tot[FEM] != _locus_position_tot[MAL];}

	void add_locus_position(const string&,const sex_t&,const unsigned int&,const double&, const unsigned int&);
	void add_locus2genome(const string& t, const unsigned int& i, TLocus* l);
	void reset_locus_vector(){_nb_locus_tot=0; _locus_vector.clear(); _sex_specific_genome=false;}
	void create_locus_tot();
	void create_locus_tot_global();
	void create_locus_tot_perTrait();
	void create_locus_position_tot();
	void create_pleiotrophy();
	void set_change_mutationRates(bool b){_change_mutationRates=b;}
	void check_pleiotrophic_loci();

	unsigned int  get_nb_locus_tot() {return _nb_locus_tot;}
    unsigned int* get_nb_locus_per_chromosome() {return _nb_locus_per_chromosome;}
    unsigned int  get_nb_chromosome() {return _nb_chromosome;}
	double*       get_locus_position_tot(sex_t SEX) {return _locus_position_tot[SEX];}
	TLocus&       get_locus_tot(const unsigned int& l) {return *_locus_tot[l];}
	TLocus**      get_locus_tot() {return _locus_tot;}
	double        get_chromosomeLength(sex_t SEX) {return _chromosomeLength[SEX];}
    unsigned int  get_nb_locus_unlinked(){return _nb_locus_unlinked;}
    unsigned int  get_nb_locus_linked(){return _nb_locus_linked;}
};
#endif
