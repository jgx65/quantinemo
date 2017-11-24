/** @file ttrait.h
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

#ifndef ttraitH
#define ttraitH

#include "types.h"
#include "simcomponent.h"
#include "tgenome.h"

class Patch;
class TIndividual;

class TTraitProto;

/**Interface for all trait types, declares all basic trait operations.
 * This pure abstract class declares the traits interface. The precise genetic architecture
 * of the trait is not defined and is up to the designer of the trait. It uses void pointers to
 * allow users to define their own structure. Existing traits use various types for their genes
 * sequence, like \c char, \c double or \c bitset. It is up to the trait's user to know what kind
 * of structure they expect.
 * The trait objects are contained in the Individual class. Their parameters are set by their
 * TTraitProto.
 **/

class TTrait{
    
    /** Genetic map for the recombination:
     * Idea: if we know how long the chromosomes (in cM) are it is possible to compute
     * the mean number of recombination events. Then it is possible to compute
     * randomly the positions of the recombination events. The connection points between the
     * chromosomes are treated as recombination events with the probability of 50%.
     * Each trait knows the positions of its loci. This allows to make the inheritance within
     * each trait separately.
     
     * Needed elements for all the trait together (static):
     *   - vector of the cumulative lengths of the chromosomes (all chromosomes are added together)
     *   - vector with the positions of the recombinations (f and m separately)
     *   - index for the starting chromosome (f and m separately)
     
     * Needed parameters for each trait separately:
     *   - vector of the positions ( from the start) of each loci
     
     * Important:
     *   - The recombination function has to be called only once for each generation of a child.
     *   - After the initial population is initialized the genetic map has to be initialized.
     *   - Before the traits are deleted the genetic map has to be deleted.
     */
    
public:
    
    
    unsigned char** sequence; 			// pointer to the allele of the genetic map sequence[locus][allele]
    
    virtual void* get_allele(const unsigned int& loc, const unsigned int& all)  const;
    
public:
    TTraitProto* pTraitProto;
    
    // the type has to be set by the inherited class
    TTrait() : sequence(0){
    }
    
    TTrait(const TTrait& T) : sequence(0){
        _copyTTraitParameters(T);  // copy the parameters of TTrait
    }
    
    // end genetic map
    
protected:
    void _copyTTraitParameters(const TTrait& T); // for the copy constructor
    
public:
    /**Called to allocate the trait's genotypic sequences. Called each time a new Individual is created (\c TIndividual::init())**/
    virtual   void            ini(TIndividual* ind);
    
    /**Called at the end of each simulation/replicate, deallocates sequence memory. **/
    virtual   void            reset () = 0;
    
    
    /** Mutation procedure, perform mutations on the genes sequence. **/
    //  virtual   void            mutate ();
    
    /** Called to set the phenotypic to a particular value or to give context-dependent value(s) to the trait.
     * @param value the value passed to the trait
     **/
    virtual   void*           set_trait (void* value) = 0;
    
    virtual   void            set_fitness_factor(){}
    virtual   void            set_fitness(double value){}
    virtual   double          set_get_fitness_factor(){return my_NAN;}
    
    virtual   void            set_from_prototype(TTraitProto* T) = 0;
    
    /**Called at the start of each replicate, sets the initial genotypes. Called by \c TIndividual::create(). **/
    virtual   void            ini_sequence (Patch* patch);
    
    /** Called to set the sequence pointer to an existing trait
     * @param seq the existing sequence pointer
     **/
    virtual   void            set_sequence (void** seq) = 0;
    
    /**Tells the trait to set its phenotype from genotype, should be used instead of get_value().**/
    virtual   void            set_value () = 0;
    virtual   void            set_value (double value) = 0;
    
    /** Genotype to phenotype mapper.
     * @return the phenotype computed from the genotype
     **/
    virtual   double          get_value () = 0;
    virtual   double          get_phenotype (){return my_NAN;}
    virtual   double          get_fitness (){return my_NAN;}
    virtual   double          get_fitnessFactor(){return 1;}
    virtual   double          get_genotype (){return my_NAN;}
    
    /** sequence accessor.
     * @return the sequence pointer
     **/
    virtual   void**          get_sequence () const = 0;
    
    /** Called to read one allele value at a particular locus.
     * @return the allelic value at position 'all' at locus 'loc'
     * @param loc locus position in the sequence
     * @param all which allele we want to read the value from
     **/
    //  virtual   void*           get_allele   (const unsigned int& loc, const unsigned int& all) const;
    
    /** Writes some info to stdout. **/
    virtual   void            show_up  () = 0;
    
    /** Returns a copy of itself.
     \b Note: call the copy constructor of the trait which should only copy the parameters values
     not the complete state of the trait (i.e. shallow copy). The copy of the sequence data is
     made through the assignment operator!**/
    virtual   TTrait*         clone () = 0;
    
    ///@name Operators
    ///@{
    /**Copies the complete state of the trait from right to left side of the operator, sequence
     data included.*/
    virtual   TTrait& operator= (const TTrait&) = 0;
    
    /**Checks for parameters equivalence, not genetic equivalence.*/
    virtual   bool operator== (const TTrait&) = 0;
    virtual   bool operator!= (const TTrait&) = 0;
    
    ///@}
    virtual~TTrait ( );
};




//------------------------------------------------------------------------------
/**TTrait setter.
 * Encapsulates the methods to set the traits parameters and to generate traits. Also stores
 * the position of the trait in the individuals trait table.
 * This class manages the file and stat handlers through its inheritance of the SimComponent
 * interface.
 */
class TTraitProto : public SimComponent {
    friend class TTrait;
    
protected:
    
    TMetapop* _popPtr;       /**The ptr to the current TMetapop.*/
    
    double***            _initAlleleFreq;     // _initAlleleFreq[patch][locus][allele] is cumulative
    unsigned int         _initAlleleFreqCols; // the number of columns used per allele, i.e. per patch
    double**             _mutationFreq;       // _mutationFreq[locus][allele] is cumulative
    mut_model_t          _mut_model;          // RMM, IMM
    ini_model_t          _ini_allele_model;   // INI_UNIF, INI_MONO, INI_DIST
    void set_ini_allelicArray (TMatrix* mat, const int& i);
    void set_allelicValues(TMatrix* mat, const unsigned int& i, const unsigned int& l, const unsigned int& a, unsigned int* cols);
    
    vector<unsigned int**>* _ini_genotypes;   // vector of initial genotypes read from FSTAT file:
    // vector[ind], 1 value: patch, then [locus][allele]
    // locus and allele start already by 0
    
    TGenomeProto*        _protoGenome;		  // pointer to the genetic map (dont' delete it)
    unsigned int*        _locus_index;		  // array containing the locus index of the genetic map
    
public:
    //parameters:
    unsigned int*    _nb_allele;              // number of alleles per locus
    unsigned int     _nb_locus;
    double           _mut_rate_mean;
    string           _type;
    unsigned int     _absolute_index;      	  // absolute trait number across all the traits
    string           _ini_fstat_file;
    
    
    TLocus* _aLocus;
    
    TLocus* get_aLocus(){return _aLocus;}
    
protected:
    
    // parameters concerning multiple instances of the same class
    unsigned int _trait_index;                // number of the instantiated class (separately for ntrl and quanti)
    // 0: if a single instantiation
    // 1 or higher:  if several traits
    
    void _copyTraitPrototypeParameters(const TTraitProto& T); // for the copy constructor
    
public:
    TTraitProto(): _popPtr(0),
    _initAlleleFreq(0), _mutationFreq(0),
    _protoGenome(0), _locus_index(0),
    _nb_allele(0), _nb_locus(0),
    _absolute_index(0), _aLocus(0), _trait_index(0){}
    
    
    ~TTraitProto();
    
    virtual   TTraitProto&  operator=   (const TTraitProto& T){return setTTraitProto(T);}
    virtual   bool          operator==  (const TTraitProto& T){return isEqualTTraitProto(T);}
    virtual   bool          operator!=  (const TTraitProto& T){return !isEqualTTraitProto(T);}
    
    /** for the derived classes to compare directly the paramters in this base class */
    virtual TTraitProto& setTTraitProto(const TTraitProto&);     // for operator=
    virtual bool isEqualTTraitProto(const TTraitProto&);         // for operator==
    virtual bool isUnequalTTraitProto(const TTraitProto&);       // for operator!=
    
    double (TTraitProto::*get_h2)(void);
    
    
    ///@name Setters
    ///@{
    virtual void    set_type              (const string& x)   {_type=x;}
    virtual void    set_nb_allele         (const unsigned int& x, const unsigned int& l){_nb_allele[l]=x;}
    virtual void    set_nb_locus          (const int& x)      {_nb_locus=x;}
    virtual void    set_absolute_index    (const unsigned int& x) {_absolute_index=x;}
    virtual void    set_trait_index       (const unsigned int& x) {_trait_index=x;}
    ///@}
    
    ///@name Getters
    ///@{
    virtual TMetapop*get_popPtr()                     const     {return _popPtr;}
    virtual string  get_type ()                      const     {return _type;}
    virtual string  get_type_index()                 const     {return _type+get_trait_indexStr_();}
    virtual unsigned int get_nb_allele(const unsigned int& l)const     {return _nb_allele[l];}
    virtual unsigned int*get_nb_allele()             const     {return _nb_allele;}
    virtual unsigned int get_nb_locus()              const     {return _nb_locus;}
    virtual unsigned int get_absolute_index()        const     {return _absolute_index;}
    virtual double* get_mutationFreq(const unsigned int& l) {return _mutationFreq[l];}
    virtual double** get_mutationFreq()                     {return _mutationFreq;}
    virtual string get_ini_fstat_file()              const     {return _ini_fstat_file;}
    unsigned int get_nb_allele_max();
    ///@}
    
    virtual string get_info(){return "not set";}
    virtual string get_info_locus();
    
    virtual unsigned int get_trait_index()      const;
    virtual string  get_trait_indexStr()        const;
    virtual string  get_trait_indexStr_()       const;
    virtual string  get_trait_indexStr_t()      const;
    
    virtual void  set_mutation_model(const string& trait,  mut_model_t* model, const unsigned int& size);
    virtual void  set_mutation_rates(const string& trait);
    virtual void  add_loci2genome(const string& trait);
    unsigned int  get_locus_max_idx_of_trait_type();
    
    virtual void  add_locus_index(const unsigned int& trait_index, const unsigned int& geome_index);
    virtual double* get_ini_alleleFreq(Patch* patch, const unsigned& locus);
    
    /**Inits the parameters, called by \c IndFactory::makePrototype().*/
    virtual   void            init (TMetapop* pMetapop) = 0;
    
    /**Re-inits the parameters, called by \c IndFactory::makePrototype().*/
    virtual   void            reset (){};            // between replicates
    virtual   void            resetTotal ();
    
    /**Creates the trait of which it is the prototype, called by \c IndFactory::makePrototype().*/
    virtual   TTrait*         hatch () = 0;
    
    /**Returns a itself. \b Note: call the copy constructor and only copy the parameters state.*/
    virtual   TTraitProto* clone () = 0;
    void    ini_base(TMetapop* pMetapop, mut_model_t* mutation_models, unsigned int nbModels);
    void    ini_paramset();
    
    /** These tasks are performed before/after each replicate */
    virtual void executeBeforeEachReplicate(const unsigned int& rep){}
    virtual void executeAfterEachReplicate(const unsigned int& rep) {}
    
    virtual void executeBeforeEachReplicate(const unsigned int& rep, map<string, Param*>* pMap){} // called by the derived classes
    virtual void executeAfterEachReplicate(const unsigned int& rep, map<string, Param*>* pMap){}  // called by the derived classes
    
    /** These tasks are performed just before/after each generation */
    virtual void executeBeforeEachGeneration(const unsigned int& gen) {}
    virtual void executeAfterEachGeneration(const unsigned int& gen) {}
    
    virtual void executeBeforeEachGeneration(const unsigned int& gen, map<string, Param*>* pMap){} // called by the derived classes
    virtual void executeAfterEachGeneration(const unsigned int& gen, map<string, Param*>* pMap){}  // called by the derived classes
    
    virtual void temporal_change(const unsigned int& gen){}  // default function
    virtual void temporal_change(const unsigned int& gen, map<string, Param*>* pMap);  // called by the derived class
    
    virtual void read_allele_file (string filename);
    
    /** get any to be specified model value */
    virtual double get_model_value(){return my_NAN;}
};



#endif //TTRAIT_H

