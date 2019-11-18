/** @file lce_coalescence.h
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

#ifndef lce_coalescenceH
#define lce_coalescenceH

//#include <map>
#include "genealogy.h"
//#include "functions.h"
//#include "simcomponent.h"
//#include "stathandler.h"
#include "filehandler.h"
#include "coal_deme.h"
#include "lifecycleevent.h"
using namespace std;

// ---------------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////////
class TLocus;
class TMetapop;
class LCE_Coalescence_base;

////////////////////////////////////////////////////////////////////////////////
class LCE_Coalescence {
protected:
    LCE_Coalescence_base*    _pCoalBase;
    
    ALLELE**         _seq;           // belongs to LCE_Coalescnece_base (don't delete it)
    
    unsigned int            _curLocus;      // the locus id (staerting at 0) which is simulated
    unsigned int            _endLocus;      // the after last locus to simulate
    
    
    dbCoalescence*           _dbCoalescence;     // db to store the pop sizes and migration numbers (don't delete it here)
    
    map<PATCH_ID, TDeme*>    _demes;             // container to store all active demes, i.e. patches: Demes[patch_id]
    vector<TDeme*>           _demesLineages;     // container with pointer to all demes with lineages
    unsigned int             _num_lineages;      // total number of lineages
    unsigned long            _curGen;            // current generation going backward in time (0 is the last generation after demographic simulation)
    
    unsigned long            _timeDivergence;    // time of divergence (parameter "divergence_time")
    
    map<PATCH_ID, vector<pair<PATCH_ID, double> > > _dispRates; // dispersal rates for before onset of sim and divergence time
    double 					_coal_migr_rateIn;
    double 					_coal_migr_rateOut;
    double 					_coal_migr_rateCorner;
    double 					_coal_migr_rate;
    double 					_coal_migr_rate_propagule;
    
    unsigned long 			_tmrca; 			// time to the MRCA
    TTNode* 				_mrca;				// the MRCA, normally the last node of geneTree, however due to recombination and reuse of nodes this may not anymore be the case.
    
    unsigned long			_treeSize; 			// the total length of the tree (at each coalescence event the childern branches are added to the length)
    
    unsigned int			_nbNodes;			// number of nodes created in total (used to have an increment of the NodeID)
    
    double                   _model_threshold;   // 0: just multiple coalescence events
    // 1: just single coalescence event
    
    bool                    allPopulated();                          // test
    bool                    curDemesAsDB(dbCoalescence* curDBs);      // test
    void                    plotDemes();                             // test
    void                    plotMigrDB(map<unsigned int, unsigned int>& curMigr); // test
    
    
    
public:
    LCE_Coalescence(LCE_Coalescence_base* p);
    
    ~LCE_Coalescence();
    
    LCE_Coalescence_base* get_pCoalBase(){return _pCoalBase;}
    
    virtual void    init(LCE_Coalescence_base* pBase);
    
    TMetapop*        get_pop_ptr();
    
    //    virtual void    ini_paramset();
    //	virtual void    init_base(Metapop* popPtr, const unsigned int& i);
    //	virtual void    init(Metapop* popPtr, const unsigned int& i);
    virtual void    run_coalescence(unsigned int from, unsigned int to);
    //	virtual void    set_up_metapop();
    
    virtual void    set_dispRates();
    
    void 			set_mrca(TTNode* n){_mrca=n; _tmrca=_mrca->time;}
    TTNode*			get_mrca(){return _mrca;}
    unsigned long   get_tmrca(){return _tmrca;}
    
    void log2console(double diffTime, unsigned int l, unsigned int nbLocus, unsigned int period);
    
    virtual void    execute(){}
    
    LCE_Coalescence*  clone () {
        return new LCE_Coalescence(_pCoalBase);
    }
    
    age_t removeAgeClass () {return 0;}
    age_t addAgeClass () {return 0;}
    age_t requiredAgeClass () {return ADULTS;}
    
    TTNode**	 _sampleNodes;			// array of all sample nodes
    
    virtual bool evacuate_lineages(TDeme* curDeme, bool all){return (this->*func_ptr_evacuate_lineages)(curDeme, all);}
    bool (LCE_Coalescence::*func_ptr_evacuate_lineages)(TDeme* curDeme, bool all);
    virtual bool evacuate_lineages_stepping_stone(TDeme* curDeme, bool all=false){return curDeme->evacuate_lineages_stepping_stone(this, all);}
    virtual bool evacuate_lineages_island(TDeme* curDeme, bool all=false)        {return curDeme->evacuate_lineages_island(this, all);}
    virtual bool evacuate_lineages_disp_matrix(TDeme* curDeme, bool all=false)   {return curDeme->evacuate_lineages_disp_matrix(this, all);}
    virtual bool evacuate_stepping_stone(TDeme* origDeme, unsigned int id, unsigned int& lineages, unsigned int steps);
    
    virtual void set_samples();
    virtual void set_ini_demeSizes();
    
    virtual void perform_popSize_regulation(dbCoalescence* curDB);
    virtual void perform_migration(dbCoalescence* curDB);
    
    virtual void migration(TDeme* curDeme, vector<pair<PATCH_ID, POP_SIZE> >& curMigr,
                           const unsigned int& demeSize, vector<TDeme*>& newDemes);
    virtual void perform_migration_before();
    void (LCE_Coalescence::*func_ptr_migration_before)(TDeme* curDeme, vector<pair<PATCH_ID, double> >& curMigr,
                                                       const unsigned int& demeSize, map<PATCH_ID, TDeme*>& newDemes);
    virtual void migration_before_matrix(TDeme* curDeme, vector<pair<PATCH_ID, double> >& curMigr,
                                         const unsigned int& demeSize, map<PATCH_ID, TDeme*>& newDemes);
    virtual void migration_before_island(TDeme* curDeme, vector<pair<PATCH_ID, double> >& curMigr,
                                         const unsigned int& demeSize, map<PATCH_ID, TDeme*>& newDemes);
    virtual void migration_before_propagule(TDeme* curDeme, vector<pair<PATCH_ID, double> >& curMigr,
                                            const unsigned int& demeSize, map<PATCH_ID, TDeme*>& newDemes);
    virtual void migration_before_1D_SS(TDeme* curDeme, vector<pair<PATCH_ID, double> >& curMigr,
                                        const unsigned int& demeSize, map<PATCH_ID, TDeme*>& newDemes);
    virtual void migration_before_2D_SS_4Neighbour(TDeme* curDeme, vector<pair<PATCH_ID, double> >& curMigr,
                                                   const unsigned int& demeSize, map<PATCH_ID, TDeme*>& newDemes);
    virtual void migration_before_2D_SS_8Neighbour(TDeme* curDeme, vector<pair<PATCH_ID, double> >& curMigr,
                                                   const unsigned int& demeSize, map<PATCH_ID, TDeme*>& newDemes);
    
    void (LCE_Coalescence::*func_ptr_perform_coalescence)();
    virtual void perform_coalescence();
    virtual void perform_coalescence_singleEvent();       // model 1
    virtual void perform_coalescence_multipleEvent();     // model 2
    virtual void perform_coalescence_mixedEvent();        // model 3
    
    
    virtual unsigned long build_tree();
    virtual void merge_demes();
    
    virtual bool check_db(const unsigned int& gen);
    
    virtual void generateGeneticData(ALLELE** seq, TLocus** genome, const unsigned int& curLoc);
    virtual void clear_tree();
    virtual void clear_tree(TTNode*);
    
    
    //getter
    inline unsigned int get_num_lineages()   {return _num_lineages;}
    inline map<PATCH_ID, TDeme*>& get_demes(){return _demes;}
    
    inline TDeme* get_demes(PATCH_ID id){		// if not present create the deme
        map<PATCH_ID, TDeme*>::iterator pos = _demes.find(id);
        if(pos != _demes.end()) return pos->second;
        TDeme * curDeme = new TDeme(get_pop_ptr(), id);
        _demes.insert(pair<PATCH_ID, TDeme*>(id, curDeme));
        return curDeme;
    }
    
    inline TDeme* get_demes(PATCH_ID id, map<PATCH_ID, TDeme*>& newDemes){		// if not present create the deme
        // is the deme already present
        map<PATCH_ID, TDeme*>::iterator pos = _demes.find(id);
        if(pos != _demes.end()) return pos->second;
        
        // or has the deme already been added
        pos = newDemes.find(id);
        if(pos != newDemes.end()) return pos->second;
        
        // no, then create a new one
        TDeme * curDeme = new TDeme(get_pop_ptr(), id);
        newDemes.insert(pair<PATCH_ID, TDeme*>(id, curDeme));
        return curDeme;
    }
    
    double get_model_threshold();
    
    
    inline vector<TDeme*>& get_demesLineages(){return _demesLineages;}
    inline unsigned long   get_curGen()         {return _curGen;}
    inline dbCoalescence&  get_dbCoalescence(const unsigned int& gen){return _dbCoalescence[gen];}
    
    
    // total tree size functions
    void				test_tree(bool alsoSize=false);
    void				print_tree(TTNode* pNode);
    unsigned long		compute_treeSize();
    void				set_treeSize(unsigned long i){_treeSize=i;}
    unsigned long		get_treeSize(){return _treeSize;}
    void				add2treeSize(long i){_treeSize+=i;}
    
    unsigned int&		get_set_nbNodes(){return _nbNodes;}
    
    void                remove_mrca_from_deme();
    
    unsigned long build_tree_simple();
    
    void migration(TDeme* curDeme, vector<pair<PATCH_ID, POP_SIZE> >& curMigr,
                   const unsigned int& demeSize, map<PATCH_ID, TDeme*>& newDemes);
    
    unsigned int get_curLocus(){return _curLocus;}
    
    
};


//////////////////////////////////////////////////////////////////////////////////
//// coalescence version based on LCE_CoaelscenceII supporting recombination
//class LCE_CoalescenceRecomb : public LCE_Coalescence {
//private:
//    //	map<PATCH_ID, TDeme*>     _demes;             // used, but stores just demes with lineages
//    //	vector<TDeme*>           _demesLineages;     // not used
//    
//public:
//    multimap<unsigned long, TTNode*> _vecNodes;// store the nodes below the recombination vecNodes<genOfRecomb, node>
//    multimap<unsigned long, TTNode*>::iterator _vecNodes_pos; // iterators to loop through vecNodes
//    multimap<unsigned long, TTNode*>::iterator _vecNodes_end;
//    
//    LCE_CoalescenceRecomb(int rank):LCE_Coalescence(rank){}
//    
//    ~LCE_CoalescenceRecomb(){}
//    
//    LCE_CoalescenceRecomb*  clone () {
//        return new LCE_CoalescenceRecomb(11);
//    }
//    
//    //	void perform_popSize_regulation(dbCoalescence* curDBs);
//    unsigned long build_tree();
//    unsigned long build_tree_simple();
//    //	void set_samples();
//    virtual void init(LCE_Coalescence_base* pBase);
//    void run_coalescence();
//    //	void set_up_metapop();
//    
//    virtual void recombine_tree(double locus_dist);
//    
//    unsigned int setup_tree_atGen(unsigned long curGen);
//    
//    // for the first tree with recombination
//    void (LCE_CoalescenceRecomb::*func_ptr_perform_coalescenceII)();
//    void perform_coalescenceII();
//    void perform_coalescence_singleEventII();       // model 1
//    void perform_coalescence_multipleEventII();     // model 2
//    void perform_coalescence_mixedEventII();        // model 3
//    
//    // for the follow up trees with recombination (for a vector of recombNodes: map<PATCH_ID, vector<TTNode*> >)
//    void (LCE_CoalescenceRecomb::*func_ptr_perform_coalescenceIII)(map<PATCH_ID, vector<TTNode*> >& recombNodes);
//    void perform_coalescenceIII(map<PATCH_ID, vector<TTNode*> >& recombNodes);
//    void perform_coalescence_singleEventIII(map<PATCH_ID, vector<TTNode*> >& recombNodes);       // model 1
//    void perform_coalescence_multipleEventIII(map<PATCH_ID, vector<TTNode*> >& recombNodes);     // model 2
//    void perform_coalescence_mixedEventIII(map<PATCH_ID, vector<TTNode*> >& recombNodes);        // model 3
//    
//    // for the follow up trees with recombination (for a single recombNode: TDeme* ; TTNode*)
//    void (LCE_CoalescenceRecomb::*func_ptr_perform_coalescenceIV)(TDeme* recombDeme, TTNode*& recombNode);
//    void perform_coalescenceIV(TDeme* recombDeme, TTNode*& recombNode);
//    void perform_coalescence_singleEventIV(TDeme* recombDeme, TTNode*& recombNode);       // model 1
//    void perform_coalescence_multipleEventIV(TDeme* recombDeme, TTNode*& recombNode);     // model 2
//    void perform_coalescence_mixedEventIV(TDeme* recombDeme, TTNode*& recombNode);        // model 3
//    
//    void perform_popSize_regulationIII(dbCoalescence* curDBs);
//    
//    void perform_migration_recomb(dbCoalescence* DB);
//    void perform_migration_before_recomb();
//    void perform_migration_recomb(dbCoalescence* DB, map<PATCH_ID, vector<TTNode*> >& recombDemes);
//    void perform_migration_recomb(dbCoalescence* DB, PATCH_ID& recombDemeID, TTNode*& recombNode);
//    void migrationRecomb(TDeme* curDeme, vector<pair<unsigned int, POP_SIZE> >& curMigr,
//                         const unsigned int& demeSize, map<PATCH_ID, TDeme*>& newDemes);
//    void migrationRecomb(TDeme* curDeme,
//                         vector<pair<PATCH_ID, POP_SIZE> >& curMigr,
//                         const unsigned int& demeSize,
//                         map<PATCH_ID, TDeme*>& newDemes,
//                         map<PATCH_ID, vector<TTNode*> >& recombDemes,
//                         map<PATCH_ID, vector<TTNode*> >& newRecombDemes);
//    void migrationRecomb(TDeme* curDeme,
//                         vector<pair<PATCH_ID, POP_SIZE> >& curMigr,
//                         const unsigned int& demeSize,
//                         map<PATCH_ID, TDeme*>& newDemes,
//                         PATCH_ID& recombDemeID,  TTNode*& recombNode);
//    
//    void update_recombDemes(multimap<unsigned long, TTNode*>::iterator& vecNodes_pos,
//                            multimap<unsigned long, TTNode*>::iterator& vecNodes_end,
//                            map<PATCH_ID, vector<TTNode*> >& recombDemes);
//    void init_recomb(PATCH_ID recombDemeID, TTNode* recombNode, unsigned long time);
//    
//    unsigned int update_lineages(unsigned long curGen);
//    void update_oldTree_coalescence(unsigned long curGen, bool onlyOnce = false);
//    void update_oldTree_migration(unsigned long curGen, PATCH_ID recombDemeID);
//    
//    void merge_demes_recomb();
//    void merge_demes(PATCH_ID& recombDemeID, TTNode* recombNode);
//    void merge_demes();
//    void remove_empty_demes(map<PATCH_ID, vector<TTNode*> >& recombDemes);
//    void remove_empty_demes(PATCH_ID recombDemeID);
//    
//    void recombine_tree_log(map<PATCH_ID, vector<TTNode*> >& recombDemes);
//    void recombine_tree_log(PATCH_ID recombDemeID, TTNode* recombNode);
//    unsigned int update_for_recomb_after_MRCA(map<PATCH_ID, vector<TTNode*> >& recombDemes);
//    unsigned int update_for_recomb_after_MRCA(PATCH_ID recombDemeID, TTNode*& recombNode, unsigned long time);
//};
//
//////////////////////////////////////////////////////////////////////////////////
//// coalescence version based on LCE_CoaelscenceII supporting recombination
//// coalescence are performed sequentially
//class LCE_CoalescenceRecombII : public LCE_CoalescenceRecomb {
//private:
//public:
//    LCE_CoalescenceRecombII(int rank):LCE_CoalescenceRecomb(rank){}
//    
//    ~LCE_CoalescenceRecombII(){}
//    
//    LCE_CoalescenceRecombII*  clone () {
//        return new LCE_CoalescenceRecombII(11);
//    }
//    
//    void init(LCE_Coalescence_base* pBase);
//    void recombine_tree(double locus_dist);
//    void recombine_tree2(double locus_dist);
//    
//    bool loopGen_duringDemo_beforeMRCA(PATCH_ID& recombDemeID, TTNode*& recombNode);
//    bool loopGen_duringDemo_afterMRCA(PATCH_ID& recombDemeID, TTNode*& recombNode);
//    bool loopGen_afterDivergence_beforeMRCA(PATCH_ID& recombDemeID, TTNode*& recombNode);
//    bool loopGen_afterDivergence_afterMRCA(PATCH_ID& recombDemeID, TTNode*& recombNode);
//    bool loopGen_beforeDivergence_beforeMRCA(PATCH_ID& recombDemeID, TTNode*& recombNode);
//    bool loopGen_beforeDivergence_afterMRCA(PATCH_ID& recombDemeID, TTNode*& recombNode);
//};




#endif
