/** @file tmetapop.h
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
 *   along with quantiNemo.  If not, see <hrttp://www.gnu.org/licenses/>.
 */

#ifndef metapopH
#define metapopH


#include "patch.h"

#ifdef _RANDOM_C11
#include "randomC11.h"
#else
#include "random.h"
#endif

class LCE;
class LCE_Breed;
class LCE_Disperse;
class LCE_Coalescence_base;
class TSimulation;
class TReplicate;

//CLASS METAPOP

/**Top class of the metapopulation structure, contains patches, implements the replicate and generation loops.
 * This class implements the two main loops of a simulation, the replicate and the generation loops. The replicate
 * loop iterates the generation loop which itself iterates the life cycle loop composed of the LCEs selected by the user.
 
 * The basic design for the metapopulation structure is a top-down chain of responsibility where the Metapop class
 * takes care of the patches it contains which are themselves concerned by the management of their individual containers.
 * The Individual class is only concerned by the management of its traits. Thereby, a metapopulation can be viewed
 * as an interleaving of containers where one container class takes care of its directly contained class only,
 * without knowledge of the upward container state.
 
 * The Metapop class thus implements methods used to manage and get information from the patches, to manage the parameters
 * necessary to build a population and to get the population state information.
 
 * <b>Population states:</b> given by the number and the position of the individuals (both spatially in demes and temporally
 * in age class containers). The Metapop::_currentAge flag is set according to the age state of the metapopulation.
 * The life cycle events modify that state by moving individuals between individuals containers within the metatpopulation.
 * The Metapop::_currentAge flag can contain the following age class bits as defined in types.h.
 * The OFFSPRNG (=1) age class bit is set whenever the offspring containers are not empty. The ADULTS (=2) age class bit is set
 * whenever the adult containers are not empty. The ALL (=3) age class is the addition of the previous tags and NONE (=0) is the negation of them.
 * The individual containers are stored in the patches and handled through the Patch class interface. Each age class is represented
 * by two containers, one for the males (index 0) and the other for the females (index 1). These containers are store in a table
 * and are accessed through their age class index as defined by the age_idx enum (see types.h). These indexes are as follows:
 * The OFFSPRNG age class has index OFFSx = 0, and the ADULTS age class has index ADLTx = 2.
 *
 */
class TMetapop: public TSimComponent, public IndFactory
{
private:
    
    /**The life cycle events*/
    vector < LCE* > _theCycle;
    vector < LCE* > _theCycleCoal;  // used for the stats and file output after the coalescence
    
    /**The stat handler for the population stats*/
    MetapopSH _statHandler;
    
    vector<TPatch*>            _vPatch;          // all patches (contains the patches and does not change over time!)
    vector<TPatch*>            _vFullPatch;      // only the populated patches (adapted during sim)
    vector<TPatch*>            _vSamplePatch;    // only sampled and populated patches [sampleID] (adapted during sim)
    
    vector<TPatch*>            _vTempPatch;      // temp vector for newly populated patches
    
    unsigned int 			  _tot_sampled_patches; // the total/max number of patches to sample (set at the start of the sim)
    
    vector<TPatch*>&        (TMetapop::*func_ptr_get_vFullPatch)();
    inline vector<TPatch*>& get_vFullPatch_eff()     {return _vFullPatch;}
    
    vector<TPatch*>&        (TMetapop::*func_ptr_get_vSamplePatch)();
    inline vector<TPatch*>& get_vSamplePatch_eff()     {return _vSamplePatch;}
    
public:
    // functions dealing with the patch containers ///////////////////////////////
    /** the following containers with the same structure are available:
     * - _vPatch:
     *       - contains ALL patches (the patches have to be deleted)
     *       - does not change over time (the patchID corresponds to the index in the vector)
     * - _vFullPatch:
     *       - contains just POPULATED patches
     *       - the order of the patches is meaningless
     *       - performance issue:
     *          - replace the vector by _vPatch if:
     *             - all patches with carrying capacity are populated at the start
     *             - patches have never a carrying capacity of 0
     *       - is dynamically adjusted:
     *          - uses the vector _vTempPatch to add newly populated patches to this vector
     *            since otherwise the iterator would be invalidated during looping
     *          - newly populated patches are added at disepersal() to _vTempPatch and at postDispersal() to _vFullPatch
     *          - newly freed patches are removed at postDiserpsal(), extinction(), aging()
     * - _vSamplePatch:
     *       - contains just POPULATED AND SAMPLED patches
     *       - the sampleID corresponds to the index in the vector (sampleID: NaN: not sampled; SAMPLED: sampled, but not populated: other: index in vector
     *       - is dynamically adjusted at the same stages as the vector _vFullPatch
     *       - performance issue:
     *          - replace the vector by _vFullPatch (and thus perhaps by _vPatch) if:
     *             - for coalescence simulations
     *             - if all patches are sampled
     */
    
    // all patches (does not change over time)
    inline vector<TPatch*>& get_vPatch()         {return _vPatch;}
    inline TPatch*          get_vPatch(const unsigned int& i){return _vPatch[i];}
    inline unsigned int    get_nbPatch()        {return (unsigned int)_vPatch.size();}
    
    // populated patches (changes over time)
    inline vector<TPatch*>& get_vFullPatch()     {return (this->*func_ptr_get_vFullPatch)();}
    inline unsigned int    get_nbFullPatch()    {return (unsigned int)get_vFullPatch().size();}
    
    // temporarily populated patches (changes over time)
    inline vector<TPatch*>& get_vTempPatch()     {return _vTempPatch;}
    inline unsigned int    get_nbTempPatch()    {return (unsigned int)_vTempPatch.size();}
    
    // sampled AND populated patches (changes over time)
    inline vector<TPatch*>& get_vSamplePatch()   {return (this->*func_ptr_get_vSamplePatch)();}
    inline unsigned int    get_nbSamplePatch()  {return (unsigned int)get_vSamplePatch().size();}  // current number of
    inline unsigned int    get_nbTotSamplePatch(){return _tot_sampled_patches;} // maximum number of patches to sample
    void                   erase_vSamplePatch(TPatch*);
    
    double                 get_fullPatch_ratio(){return (double)get_nbFullPatch()/get_nbPatch();}
    
    // functions to set the function pointers of the containers
    void set_func_ptr_noSample_noFull();
    void set_func_ptr_noSample_withFull();
    void set_func_ptr_noSample_noFull_coal();
    void set_func_ptr_noSample_withFull_coal();
    void set_func_ptr_withSample_noFull();
    void set_func_ptr_withSample_withFull();
    
    vector<LCE*>& get_theCycle() {return _theCycle;}
    vector<LCE*>& get_theCycleCoal() {return _theCycleCoal;}
    
private:
    // pointers to LCEs
    LCE_Breed* _pBreed_LCE;             /** pointer to the breed LCE */
    LCE_Disperse* _pDisperse_LCE;       /** pointer to the disperse LCE */
    LCE_Coalescence_base* _pCoalescence_LCE;	/** pointer to the coalescence LCE */
    TSimulation* _pSimulation;
    TReplicate* _pReplicate;
    
    // output (stats and files)
    age_t        _sampled_age_class;    // all the age classes used for any output
    unsigned int _isSampled[3];         // to prevent the multiple calling of the function set_sampledInds() //_sampledInds[0:gen; 1:rep; 2: age]
    
    //parameters:
    unsigned int _patchNbr;             /**Number of patches in the population.*/
    unsigned int _generations;          /**Number of generations to iterate.*/
    
    //counters:
    unsigned int _currentGeneration;    /**The current generation in the generation loop, starts at 1.*/
    age_t        _currentAge;          	/**The current age class, might be changed by the LCEs.*/
    
    /** initial sex ratio */  // param input: sex_ratio = nbMal/nbFem
    double  _sexInitRatio;   // internaly (variable _seInitRatio: nbMal/(nbMal+nbFem): 0: only females => hermaphrodite; 0.5: f = m
    
    FileServices* _service;
    
    TSelection* 	_pSelection;   // pointer to the selection stuff if needed
    
    unsigned int  _total_carrying_capacity;
    
    /** parameters for denstiy dependent dispersal rate */
    double** _density_threshold;					// allows to change the disp rate depending on the density (0: patch, 1: density, 2: change)
    Param**  _density_threshold_param;    // pointer to the param
    unsigned int _density_threshold_nbChanges;
    bool     _vSamplePatchUsed;
    bool     _vFullPatchUsed;
    bool     _sampleAllOrNothing;       // 0: compute stats always if possible; 1: compute only stats if entire sampling can be performed
    bool     _sampleAllOrNothingCur;    // same as above but for the current settings (sampling schema may change)
        
public:
    // stat parameters unique for a replicate (was static, but this does not work for mulithreading)
    unsigned int _current_replicate;        // the current replicate (starting at 0) (does not change during replicate sim!!!)
    //unsigned int _current_generation;       // the current generation (starting at 0)
    unsigned int _current_index_stat_db;    // the current index for the array _val of the stat_db where to store the stat
    unsigned int _last_nbSamplePatch;    /** the number of sampled AND populated patches of the last NON-ZERO stats computation */
    unsigned int _current_nbSamplePatch; /** the number of sampled AND populated patches of this generation */
    void update_patch_states();
    
    list<StatRecBaseAll*>::iterator _current_statRecBase;
    
    TReplicate* get_pReplicate(){return _pReplicate;}

    
    ////////////////////////////////////////////////////////////////////////////
    
    bool  get_sampleAllOrNothing(){return _sampleAllOrNothing;}
    
    TMetapop(TReplicate* p, unsigned int rep);
    virtual ~TMetapop();
    void init_paramset();
    
    /**Inits the population parameters from the ParamSet and builds the pop (adds patches), the prototypes and the life cycle.
     Called at the start of each simulation, resets the individual garbage collector.
     @param traits the map of the selected trait prototypes from the SimBuilder
     @param LCEs the map of the selected LCEs from the SimBuilder
     @callgraph
     */
    bool init( map< string, TTraitProto* >& traits, map< int, LCE* >& LCEs );
    
    void set_sampledInds(age_t AGE, bool storeState=true); // sets the current sampling numbers (called at each generation with output)
    
    void          flush      (sex_t SEX, age_idx AGE);
    void          flush      (age_idx AGE);
    void          flush      (age_t AGE);
    void          flush      ();
    
    /** function which changes parameters over time */
    void temporal_change();
    void temporal_change(const unsigned int& gen);
    
    /**Called to empty the patches, individuals are move to the garbage collector.*/
    void reset();
    
    /**Called at the end of each simulation, empties the pop and the garbage collector, Individuals are destroyed.*/
    void clear();
    
    /**Creates the list of LCEs from the selected ones in the SimBuilder instance.
     \b Note: the LCEs inserted in the list are the LCE templates stored in the SimBuilder. They aren't copied or cloned.
     **/
    void setLifeCycle(map< int, LCE* >& lifeCycle);
    
    ///@name Main loops
    ///@{
    /**Replicate loop, iterates the life cycle \a _replicates times.
     @callgraph
     */
    
    /** function which is executed before/after each replicate */
    void executeBeforeEachReplicate(const unsigned int& rep);
    void executeAfterEachReplicate(const unsigned int& rep);
    
    /** function which is executed before/after each generation */
    void executeBeforeEachGeneration(const unsigned int& gen);
    void executeAfterEachGeneration(const unsigned int& gen);
    
    /**Life cycle loop, executes the list of LCEs \a _generations times.
     @param startTime the starting time of the current replicate.
     @callgraph
     */
    void Loop_generation(time_t startTime);
    
    ///@}
    
    void createPopulations();
    
    ///@name Population builders
    ///@{
    
    /** function to set the initial population sizes */
    void setInitPopulationSizes();
    void setCarryingCapacities();
    bool setSampleSizes(bool set_sampleID=true);
    bool set_sampled_patches(bool allGens = false);
    
    template<typename T>
    bool setPatchParam(string name,                                // parameter name (without "_fem" or "_mal")
                       void (TPatch::*setSexSpecific)(T, sex_t), // set a sex specific value
                       void (TPatch::*setGeneral)(T),            // set a general value
                       void (TPatch::*reset)(T, T, T),           // function to set all three params at once
                       T (TPatch::*getSexSpecific)(sex_t),       // get a sex specific value
                       T (TPatch::*getGeneral)());               // get a general value
    
    // sex specific parameter
    template<typename T>
    bool setPatchParam(string name, sex_t SEX,
                       void (TPatch::*pt2Func)(T, sex_t),
                       void (TPatch::*pt2reset)(T, T, T));
    template<typename T>
    void setPatchParam(string name, TMatrix* m, sex_t SEX,
                       void (TPatch::*pt2Func)(T, sex_t),
                       void (TPatch::*pt2reset)(T, T, T));
    template<typename T>
    void setPatchParam(string name, T val, sex_t SEX,
                       void (TPatch::*pt2Func)(T, sex_t));
    
    // general parameter
    template<typename T>
    bool setPatchParam(string name, void (TPatch::*pt2Func)(T) ,
                       void (TPatch::*pt2reset)(T, T, T));
    template<typename T>
    void setPatchParam(string name, TMatrix* m, void (TPatch::*pt2Func)(T),
                       void (TPatch::*pt2reset)(T, T, T));
    template<typename T>
    void setPatchParam(string name, T val, void (TPatch::*pt2Func)(T));
    
    template <typename T>
    void setSexSpecificPatchParam_sexRatio(T (TPatch::*getter)(), void (TPatch::*setter)(T, sex_t)){
        for(unsigned int i = 0; i < _patchNbr; ++i){
            _vPatch[i]->setSexSpecificPatchParam_sexRatio(getter, setter, _sexInitRatio);
        }
    }
    
    template <typename T>
    void setGeneralPatchParam_fem(T (TPatch::*getter)(sex_t), void (TPatch::*setter)(T)){
        for(unsigned int i = 0; i < _patchNbr; ++i){
            _vPatch[i]->setGeneralPatchParam_fem(getter, setter);
        }
    }
    
    template <typename T>
    void setGeneralPatchParam_sum(T (TPatch::*getter)(sex_t), void (TPatch::*setter)(T)){
        for(unsigned int i = 0; i < _patchNbr; ++i){
            _vPatch[i]->setGeneralPatchParam_sum(getter, setter);
        }
    }
    
    
    /**Sets the first generation of each replicates.
     @callgraph
     */
    void reset_metapopulation();
    age_t get_requiredAgeClass();
    void setPopulation (const age_t& AGE);        // random initialization
    bool setPopulation_FSTAT (const age_t& AGE);  // if genotypes are given by FSTAT  file
    bool allPatchPopulated_atInitialization();
    ///@}
    
    /**Sets the Patch trait optima.
     @param optima a matrix containing the patches optima and intensity (see manual). */
    void set_patch_parameter(unsigned int nbTrait, string name, string name_full,
                             void (TPatch::*pt2Func)(double*, sex_t));
    void set_patch_value_byValue(unsigned int nbTrait, double value, sex_t SEX,
                                 void (TPatch::*pt2Func)(double*, sex_t));
    void set_patch_value_byMatrix(unsigned int nbTrait, TMatrix* m, sex_t SEX, string name,
                                  void (TPatch::*pt2Func)(double*, sex_t));
    void set_patch_parameter_array(unsigned int nbTrait, string name, string name_full,
                                   void (TPatch::*pt2Func)(unsigned int, double*, unsigned int, sex_t));
    void set_patch_array_byArray(unsigned int nbTrait, TTree<unsigned int, double>* m, sex_t SEX,
                                 void (TPatch::*pt2Func)(unsigned int, double*, unsigned int, sex_t));
    void set_patch_array_byMatrix(unsigned int nbTrait, TTree<unsigned int, double>* m, sex_t SEX, string name,
                                  void (TPatch::*pt2Func)(unsigned int, double*, unsigned int, sex_t));
    
    void set_patch_parameter_ofTrait(TTraitProto* pTrait, unsigned int curTrait, string trait, string name,
                                     string name_full, void (TPatch::*pt2Func)(unsigned int, double, sex_t));
    void set_patch_value_byValue_ofTrait(unsigned int curTrait, double value, sex_t SEX,
                                         void (TPatch::*pt2Func)(unsigned int, double, sex_t));
    void set_patch_value_byMatrix_ofTrait(unsigned int curTrait, TMatrix* m, sex_t SEX, string name,
                                          void (TPatch::*pt2Func)(unsigned int, double, sex_t));
    void set_patch_parameter_array_ofTrait(TTraitProto* pTrait, unsigned int curTrait, string trait, string name, string name_full,
                                           void (TPatch::*pt2Func)(unsigned int, double*, unsigned int, sex_t));
    void set_patch_array_byArray_ofTrait(unsigned int curTrait, TMatrix* m, sex_t SEX, string name,
                                         void (TPatch::*pt2Func)(unsigned int, double*, unsigned int, sex_t));
    void set_patch_array_byMatrix_ofTrait(unsigned int curTrait, TMatrix* m, sex_t SEX, string name,
                                          void (TPatch::*pt2Func)(unsigned int, double*, unsigned int, sex_t));
    
    ///@name Implementations
    ///@{
    //SimComponent implementation:
    virtual void loadFileServices ( FileServices* loader ) {}
    
    virtual void loadStatServices ( StatServices* loader ) {loader->attach(&_statHandler);}
    
  //  vector<TIndividual*>& getIndividuals(const sex_t& SEX, const age_idx& AG);
    unsigned int getGenerations        ( ) {return _generations;}
    unsigned int getPatchNbr           ( ) {return _patchNbr;}
    
    /**Returns the mean number of generations performed per replicate.*/
    LCE_Disperse* get_pDisperse_LCE    ( ) {return _pDisperse_LCE;}
    LCE_Breed*    get_pBreed_LCE       ( ) {return _pBreed_LCE;}
    double        get_sexInitRatio     ( ) {return _sexInitRatio;}
    TGenomeProto* get_protoGenome      ( ) {return _protoGenome;}
    unsigned int* get_isSampled        ( ) {return _isSampled;}
    TSimulation*    get_pSimulation();
    TStat_db*       get_pStat_db();
    void            set_pStat_db(TStat_db* p);
        
    
    void set_LCE_pointers(map< int,LCE* >& LCEs);
    
    ///@}
    
    ///@name Population state interface
    ///@{
    unsigned int getCurrentReplicate(){return _current_replicate;}
    unsigned int getReplicates   ( );
    void         setCurrentGeneration  (const unsigned int& i) {_currentGeneration = i;}
    unsigned int getCurrentGeneration  ( ) {return _currentGeneration;}
    age_t        getCurrentAge         ( ) {return _currentAge;}
    
    /**
     * @return the current replicate counter string
     */
    string  getReplicateCounter ();
    string  getReplicateCounter_ ();
    string  getReplicateCounter_r ();
    
    // gat parameters of TSimulation
    string get_working_directory();
    string get_iniFile_directory();
    string get_exe_directory();
    string getBaseFileName();
    string getSimfolder();
    string getSimfolderShort();


    bool individual_container_ok();
    
    /**Sets the age flag.
     @param age the current age. */
    void setCurrentAge(age_t age) {_currentAge = age;}
    
    /**Checks if the population still contains at least one individual in any sex or age class.*/
    bool isAlive( ) {return size() != 0;}
    
    /**Get the total number of individuals present in the population, all sex and age classes together.*/
    unsigned int size( ) {return size(ALL);}
    
    /**Interface to get the size of a praticular age and sex class(es).
     @param AGE age class flags
     @param SEX sex class
     */
    unsigned int size ( sex_t SEX, age_t AGE );
    unsigned int size ( sex_t SEX, age_idx AGE );
    
    /**Interface to get the size of a praticular age and sex class within a patch.
     @param AGE age class flags
     @param SEX sex class
     @param deme the focal patch
     */
    unsigned int size (sex_t SEX, age_t AGE, unsigned int deme);
    unsigned int size (sex_t SEX, age_idx AGE, unsigned int deme);
    
    /**Simplified interface to get the size of both sexes of the appropriate age class(es) in the whole population.
     @param AGE age class flags
     */
    unsigned int size ( age_t AGE ){
        return size( FEM, AGE ) + size( MAL, AGE );
    }
    unsigned int size ( age_idx AGE ){
        return size( FEM, AGE ) + size( MAL, AGE );
    }
    
    /**Simplified interface to get the size of both sexes of the appropriate age class(es) in one patch.
     @param AGE age class flags
     @param deme the focal deme
     */
    unsigned int size ( age_t AGE, unsigned int deme ){
        return size( FEM, AGE, deme ) + size( MAL, AGE, deme );
    }
    unsigned int size ( age_idx AGE, unsigned int deme ){
        return size( FEM, AGE, deme ) + size( MAL, AGE, deme );
    }
    
    /**Get the total number of individuals present in the population, all sex and age classes together.*/
    unsigned int sampleSize( ) {return sampleSize(ALL);}
    
    /**Interface to get the size of a praticular age and sex class(es).
     @param AGE age class flags
     @param SEX sex class
     */
    unsigned int sampleSize ( sex_t SEX, age_t AGE );
    unsigned int sampleSize ( sex_t SEX, age_idx AGE );
    
    /**Interface to get the size of a praticular age and sex class within a patch.
     @param AGE age class flags
     @param SEX sex class
     @param deme the focal patch
     */
    unsigned int sampleSize (sex_t SEX, age_t AGE, unsigned int deme);
    unsigned int sampleSize (sex_t SEX, age_idx AGE, unsigned int deme);
    
    /**Simplified interface to get the size of both sexes of the appropriate age class(es) in the whole population.
     @param AGE age class flags
     */
    unsigned int sampleSize ( age_t AGE ){
        return sampleSize( FEM, AGE ) + sampleSize( MAL, AGE );
    }
    unsigned int sampleSize ( age_idx AGE ){
        return sampleSize( FEM, AGE ) + sampleSize( MAL, AGE );
    }
    
    /**Simplified interface to get the size of both sexes of the appropriate age class(es) in one patch.
     @param AGE age class flags
     @param deme the focal deme
     */
    unsigned int sampleSize ( age_t AGE, unsigned int deme ){
        return sampleSize( FEM, AGE, deme ) + sampleSize( MAL, AGE, deme );
    }
    unsigned int sampleSize ( age_idx AGE, unsigned int deme ){
        return sampleSize( FEM, AGE, deme ) + sampleSize( MAL, AGE, deme );
    }
    
    /**Returns a pointer to the appropriate individual.
     @param SEX sex class container index
     @param AGE age class container index
     @param at the index of the individual in its container
     @param deme the patch where to grab the individual*/
    TIndividual* get (const sex_t& SEX, const age_idx& AGE, unsigned int at, unsigned int deme);
    
    /**Moves an individual from a deme to an other one, both demes sizes are modified.
     @param SEX sex class container index
     @param from_age age class container index in the deme of origin
     @param from_deme index of the deme of origin
     @param to_age age class container index in the destination deme
     @param to_deme index of the destination deme
     @param at index of the focal individual in the 'from' deme
     */
    void move (sex_t SEX, age_idx from_age, unsigned int from_deme,
               age_idx to_age, unsigned int to_deme, unsigned int at);
    void move_random (sex_t SEX, age_idx from_age, unsigned int from_deme,
                      age_idx to_age, unsigned int to_deme, unsigned int nbMigr);
    void move_random (sex_t SEX, age_idx from_age, TPatch* fromDeme,
                      age_idx to_age, TPatch* toDeme, unsigned int nbMigr);
    void copyMove_random_withReplacement (sex_t SEX, age_idx from_age, TPatch* fromDeme,
                                          age_idx to_age, TPatch* toDeme, unsigned int nbMigr);
    
    /**
     * @return the current generation counter string
     */
    string  getGenerationCounter ();
    string  getGenerationCounter_ ();
    string  getGenerationCounter_g ();
    
    /** returns tre if only one sex is used for the simulations */
    void    set_SexInitRatio(map< int,LCE* >& LCEs);
    
    void            set_service(FileServices* ptr){_service=ptr;}
    FileServices*   get_service(){return _service;}
    
    void    set_generations(unsigned int i) {_generations = i;}
    
    void    regulate_selection_fitness_patch(age_idx AGE, unsigned int* Kmal=NULL, unsigned int* Kfem=NULL);
    void    regulate_selection_fitness_metapop(age_idx AGE);
    void    regulate_selection_fitness_hard(age_idx AGE);
    
    void set_total_carrying_capacity(){_total_carrying_capacity = get_total_carrying_capacity_bothSex();}
    unsigned int get_total_carrying_capacity_bothSex();
    unsigned int get_total_carrying_capacity(const sex_t& SEX);
    unsigned int get_total_carrying_capacity(){return _total_carrying_capacity;}
    unsigned int get_total_iniSize_bothSex();
    void get_settingsStats(unsigned int& nbPatch,unsigned int& K, unsigned int& NiniPatch, unsigned int& Nini, unsigned int& samplePatch, unsigned int& nbSample);
    
    TSelection* get_pSelection()              {return _pSelection;}
    
    /** functions for changing dispersal rate following pop density of a certain patch */
    void change_disp_rate_after_density(const int& gen);
    void set_change_disp_rate_after_density();
    
    /** newly populated patch */
    void (TMetapop::*func_ptr_new_fullPatch)(TPatch* curPatch);
    inline void new_fullPatch(TPatch* curPatch){(this->*TMetapop::func_ptr_new_fullPatch)(curPatch);}
    inline void new_fullPatch_noSample_noFull(TPatch* curPatch);
    inline void new_fullPatch_withSample_noFull(TPatch* curPatch);
    inline void new_fullPatch_noSample_withFull(TPatch* curPatch);
    inline void new_fullPatch_withSample_withFull(TPatch* curPatch);
    
    /** newly freed patch */
    void (TMetapop::*func_ptr_new_emptyPatch)(vector<TPatch*>::iterator& curPop, vector<TPatch*>::iterator& endPop);
    inline void new_emptyPatch(vector<TPatch*>::iterator& curPop, vector<TPatch*>::iterator& endPop){ return (this->*TMetapop::func_ptr_new_emptyPatch)(curPop, endPop);}
    inline void new_emptyPatch_withSample_withFull(vector<TPatch*>::iterator& curPop, vector<TPatch*>::iterator& endPop);
    inline void new_emptyPatch_withSample_noFull(vector<TPatch*>::iterator& curPop, vector<TPatch*>::iterator& endPop);
    inline void new_emptyPatch_noSample_noFull(vector<TPatch*>::iterator& curPop, vector<TPatch*>::iterator& endPop);
    inline void new_emptyPatch_noSample_withFull(vector<TPatch*>::iterator& curPop, vector<TPatch*>::iterator& endPop);
    inline void new_emptyPatch_noSample_noFull_coal(vector<TPatch*>::iterator& curPop, vector<TPatch*>::iterator& endPop);
    inline void new_emptyPatch_noSample_withFull_coal(vector<TPatch*>::iterator& curPop, vector<TPatch*>::iterator& endPop);
    
    // add the _vTempPatch container to the _vFullPatch container
    void (TMetapop::*func_ptr_add_tempPatch)();
    void add_tempPatch(){(this->*func_ptr_add_tempPatch)();}
    void add_tempPatch_withSample_withFull();
    void add_tempPatch_withSample_noFull();
    void add_tempPatch_noSample_withFull();
    void add_tempPatch_noSample_noFull();
    void add_tempPatch_noSample_withFull_coal();
    void add_tempPatch_noSample_noFull_coal();
    
    /** all functions used for the coalescence simulation */
    bool init_coal( map< string, TTraitProto* >& traits, map< int, LCE* >& LCEs );
    void setLifeCycle_coal(map< int, LCE* >& lifeCycle);
    void Loop_generation_coal(time_t startTime);
    void createPopulations_coal();
    void setPopulation_coal (const age_t& AGE);
    bool isCoalescence(){return (bool)_pCoalescence_LCE;}
    LCE_Coalescence_base* getCoalescence(){return _pCoalescence_LCE;}
    void add_immigrants(const unsigned int& to, const unsigned int& from, const unsigned int& nb);
    void store_popSizes();
    bool set_samples_coalescence();
    void run_coalescence();
    void set_func_pointer_coal();
    void set_func_pointer();
    
    void printGenRep2File(ostream& FH);
    
    RAND& rand();
    
    ///@}
};


inline unsigned int TMetapop::size ( sex_t SEX, age_t AGE )
{
    unsigned int size=0;
    vector<TPatch*>::iterator curPop, endPop;
    for(curPop=get_vFullPatch().begin(), endPop=get_vFullPatch().end(); curPop!=endPop; ++curPop) {
        size += (*curPop)->size(SEX, AGE);
    }
    return size;
}

inline unsigned int TMetapop::size (sex_t SEX, age_idx AGE, unsigned int deme)
{
    return get_vPatch(deme)->size(SEX, AGE);
}

inline unsigned int TMetapop::size ( sex_t SEX, age_idx AGE )
{
    unsigned int size=0;
    vector<TPatch*>::iterator curPop, endPop;
    for(curPop=get_vFullPatch().begin(), endPop=get_vFullPatch().end(); curPop!=endPop; ++curPop) {
        size += (*curPop)->size(SEX, AGE);
    }
    return size;
}

inline unsigned int TMetapop::size (sex_t SEX, age_t AGE, unsigned int deme){
    return get_vPatch(deme)->size(SEX, AGE);
}


inline unsigned int TMetapop::sampleSize ( sex_t SEX, age_t AGE )
{
    unsigned int size=0;
    vector<TPatch*>::iterator curPop, endPop;
    for(curPop=get_vFullPatch().begin(), endPop=get_vFullPatch().end(); curPop!=endPop; ++curPop) {
        size += (*curPop)->sampleSize(SEX, AGE);
    }
    return size;
}

inline unsigned int TMetapop::sampleSize (sex_t SEX, age_idx AGE, unsigned int deme)
{
    return get_vPatch(deme)->sampleSize(SEX, AGE);
}

inline unsigned int TMetapop::sampleSize ( sex_t SEX, age_idx AGE )
{
    unsigned int size=0;
    vector<TPatch*>::iterator curPop, endPop;
    for(curPop=get_vFullPatch().begin(), endPop=get_vFullPatch().end(); curPop!=endPop; ++curPop) {
        size += (*curPop)->sampleSize(SEX, AGE);
    }
    return size;
}

inline unsigned int TMetapop::sampleSize (sex_t SEX, age_t AGE, unsigned int deme){
    return get_vPatch(deme)->sampleSize(SEX, AGE);
}


inline TIndividual* TMetapop::get (const sex_t& SEX, const age_idx& AGE, unsigned int at, unsigned int deme){
    return get_vPatch(deme)->get(SEX, AGE, at);
}

inline void TMetapop::move (sex_t SEX, age_idx from_age, unsigned int from_deme,
                           age_idx to_age, unsigned int to_deme, unsigned int at){
    _vPatch[to_deme]->add(SEX, to_age, get(SEX, from_age, at, from_deme));
    _vPatch[from_deme]->remove(SEX, from_age, at);
    
}

inline bool TMetapop::individual_container_ok ()
{
    vector<TPatch*>::iterator curPop, endPop;
    for(curPop=get_vFullPatch().begin(), endPop=get_vFullPatch().end(); curPop!=endPop; ++curPop) {
         if(!(*curPop)->individual_container_ok()) return false;
    }
    return true;
}

#endif
