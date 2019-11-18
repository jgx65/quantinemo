/** @file lce_coalescence_base.h
 *
 *   Copyright (C) 2010 Samuel Neuenschwander <samuel.neuenschwander@unil.ch>
 *   Copyright (C) 2018 Frederic Michaud <frederic.a.michaud@gmail.com>

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
// ---------------------------------------------------------------------------

#ifndef lce_coalescence_baseH
#define lce_coalescence_baseH

#include <map>
//#include "genealogy.h"
//#include "functions.h"
//#include "simcomponent.h"
#include "stathandler.h"
#include "filehandler.h"
#include "coal_deme.h"
#include "lifecycleevent.h"

using namespace std;


class LCE_CoalescenceSH;
class LCE_CoalescenceFH;
class LCE_CoalescenceFH_mrca;
class LCE_CoalescenceFH_lineages;
class LCE_CoalescenceFH_popSizes;

class LCE_Coalescence_base : public LCE {
private:
    unsigned int             _nbSamples;         // number of samples (haploid)
    map<unsigned int, unsigned int> _iniSamples; // container with the initial samples
    ALLELE**          _seq;               // _seq[sample][locus]
    double                   _model_threshold;   // 0: just multiple coalescence events
    unsigned int             _nbGen;             // number of generations to store (time of the demographic simulation)
    unsigned int             _nbLocus;           // total number of loci
    dbCoalescence*           _dbCoalescence;     // db to store the pop sizes and migration numbers
    LCE_CoalescenceSH*       _stats;
    
    // stats
    LCE_CoalescenceFH*       _writer;
    LCE_CoalescenceFH_mrca*  _writerMRCA;
    LCE_CoalescenceFH_lineages*  _writerLineages;
    LCE_CoalescenceFH_popSizes*  _writerPopSizes;
    
    unsigned long*           _vMRCA;             // vector of the time to the MRCA of each locus




public:
    LCE_Coalescence_base(int rank = my_NAN);
    virtual ~LCE_Coalescence_base();
    
    
    virtual LCE_Coalescence_base*  clone() {return new LCE_Coalescence_base(11);}
    virtual age_t removeAgeClass () {return 0;}
    virtual age_t addAgeClass () {return 0;}
    virtual age_t requiredAgeClass () {return ADULTS;}
    
    virtual void execute();
   
    virtual void loadStatServices(StatServices* loader);
    virtual void loadFileServices ( FileServices* loader );
       
    virtual void temporal_change(const unsigned int& gen);

    
    void set_up_metapop();
    void set_up_metapop_empty();
    void init(TMetapop* popPtr, const unsigned int& i);
    void ini_paramset();
    void set_samples(const map<unsigned int, unsigned int>& m, const unsigned int& nb);
    
    inline unsigned int get_nbSamples(){return _nbSamples;}
    inline map<unsigned int, unsigned int>& get_iniSamples(){return _iniSamples;} // container with the initial samples
    inline double get_model_threshold(){return _model_threshold;}
    inline dbCoalescence& get_dbCoalescence(unsigned int gen){return _dbCoalescence[gen];}
    inline dbCoalescence* get_dbCoalescence(){return _dbCoalescence;}
    inline unsigned int get_nbGen(){return _nbGen;}
    inline unsigned int get_nbLocus(){return _nbLocus;}
    inline ALLELE** get_seq(){return _seq;}
    inline unsigned long* get_vMRCA() {return _vMRCA;}

    LCE_CoalescenceFH*           get_writer(){return _writer;}
    LCE_CoalescenceFH_mrca*      get_writerMRCA(){return _writerMRCA;}
    LCE_CoalescenceFH_lineages*  get_writerLineages(){return _writerLineages;}
    LCE_CoalescenceFH_popSizes*  get_writerPopSizes(){return _writerPopSizes;}

    
    inline void set_vMRCA(unsigned int l, unsigned long t){assert(l<_nbLocus); _vMRCA[l]=t;}

    
    // functions dealing with the database
    void add_immigrants(const unsigned int& gen, const unsigned int& to, const unsigned int& from, const unsigned int& nb);
    void set_popSize(const unsigned int& gen, const unsigned int& i, const unsigned int& size);
    void store_popSizes(const unsigned int& curGen);
};

void run_coal_locus(LCE_Coalescence_base* pCoal, unsigned int from,
                            unsigned int to, unsigned int curThread);


////////////////////////////////////////////////////////////////////////////////
class LCE_CoalescenceSH  : public StatHandler<LCE_CoalescenceSH> {
private:
    LCE_Coalescence_base*  pCoal;
    
public:
    
    LCE_CoalescenceSH (LCE_Coalescence_base* TT): pCoal(TT) {}
    LCE_CoalescenceSH(){}
    virtual ~LCE_CoalescenceSH ( ){}
    
    double get_MRCA_mean() {return ARRAY::mean(pCoal->get_vMRCA(), pCoal->get_nbLocus());}
    double get_MRCA_var()  {return ARRAY::var(pCoal->get_vMRCA(), pCoal->get_nbLocus());}
    
    bool init ( );
    
    bool setStatRecorders (const string& t);
    
    string getName() {return "CoalescenceSH";}
};

////////////////////////////////////////////////////////////////////////////////
/**A file handler to save the neutral markers genotypes in the FSTAT format*/
class LCE_CoalescenceFH: public FileHandler {
private:
    LCE_Coalescence_base*  _pCoal;
    string  _NEXUS_filename;     // file name is stored in order not to recreate it all the time
    string* _locusTree;          // for each locus the tree as text
    
public:
    void    write_NEXUS_tree(unsigned int l, LCE_Coalescence* pCoal);
    void    write_NEXUS();
    void    set_pCoal(LCE_Coalescence_base* p){_pCoal = p;}
    void    create_locusTree(){_locusTree = new string[_pCoal->get_nbLocus()];}
    void    update(){}; // needed to supresss the Fstat writing
    
    LCE_CoalescenceFH () : _pCoal(0), _locusTree(0) {}
    virtual ~LCE_CoalescenceFH ( ) {if(_locusTree) delete[] _locusTree;}
};

////////////////////////////////////////////////////////////////////////////////
/**A file handler to save the times to the MRCAs*/
class LCE_CoalescenceFH_mrca: public FileHandler {
private:
    LCE_Coalescence_base*  _pCoal;
public:
    void    write_mrca();
    void    set_pCoal(LCE_Coalescence_base* p){_pCoal = p;}
    void    update(){}; // needed to supresss the Fstat writing
    
    LCE_CoalescenceFH_mrca () : _pCoal(0) {}
    virtual ~LCE_CoalescenceFH_mrca ( ) { }
};


////////////////////////////////////////////////////////////////////////////////
/**A file handler to save the lineages and population sizes over time*/
class LCE_CoalescenceFH_lineages: public FileHandler {
private:
    LCE_Coalescence_base*  _pCoal;
    string* _lineage_filename; // file name for each lcous is stored in order not to recreate it all the time
public:
    void    write_lineages_start(unsigned int locus);
    void    write_current_lineages(unsigned int nbGen, unsigned long curGen,
                                     dbCoalescence* curDBs, LCE_Coalescence* pCoal);
    void    write_current_lineages(unsigned int nbGen, unsigned long curGen,
                                     LCE_Coalescence* pCoal);
    void    write_lineages_end();
    void    set_pCoal(LCE_Coalescence_base* p){_pCoal = p;}
    void    create_lineage_filename(){_lineage_filename = new string[_pCoal->get_nbLocus()];}
    void    update(){}; // needed to supresss the Fstat writing
    
    LCE_CoalescenceFH_lineages () : _pCoal(0), _lineage_filename(0) {}
    virtual ~LCE_CoalescenceFH_lineages ( ) {delete[] _lineage_filename;}
};

////////////////////////////////////////////////////////////////////////////////
/**A file handler to store over time the population sizes for the given patches */
class LCE_CoalescenceFH_popSizes: public FileHandler {
private:
    LCE_Coalescence_base*  _pCoal;
    // string            _filename;
public:
    void    write_popSizes();
    void    set_pCoal(LCE_Coalescence_base* p){_pCoal = p;}
    void    update(){}; // needed to supresss the Fstat writing
    
    LCE_CoalescenceFH_popSizes () : _pCoal(0) {}
    virtual ~LCE_CoalescenceFH_popSizes ( ) { }
};

////////////////////////////////////////////////////////////////////////////////



#endif

