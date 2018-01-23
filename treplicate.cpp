/** @file treplicate.cpp
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

#include "treplicate.h"

#include "tsimulation.h"
#include "ttneutral.h"
#include "ttquanti.h"
#include "lce_disperse.h"
#include "lce_coalescence_base.h"
#include "lce_breed.h"
#include "lce_misc.h"
#include "lce_regulation.h"
#include "lce_extinction.h"

using namespace std;


//----------------------------------------------------------------------------------------
// TReplicate
// ----------------------------------------------------------------------------------------
TReplicate::TReplicate(TSimulation* ptr, map<string, string> params,
                       map<string, string> keys, unsigned int rep,
                       unsigned int curThread, unsigned int nbThreads):
_current_replicate(rep), _generations(0), _FileServices(0), _StatServices(0),
_pSimulation(0), _pStat_db(0), _thePop(0), _current_thread(curThread),
_nbThreads(nbThreads), rand(0){
#ifdef _DEBUG
    message(" TReplicate::TReplicate()\n");
#endif
    _pSimulation = ptr;
    assert(_pSimulation->randEngines);
    rand = &(_pSimulation->randEngines[_current_thread]);     // attribut the
    
    init_paramset();
    
    Param* pParam = _paramSet.get_param("seed");
    pParam->set_arg(_pSimulation->_vSeeds[_current_thread]);

    run_replicate(params, keys, _current_replicate);          // run the replicate

}

//----------------------------------------------------------------------------------------
// TReplicate
// ----------------------------------------------------------------------------------------
/** constructor used to initalize that sim:
 - test the parameters and write log file
 - initialize the stats and set up stats_db
 */
TReplicate::TReplicate(TSimulation* ptr): _current_replicate(my_NAN), _generations(0),
_FileServices(0), _StatServices(0), _pSimulation(0), _pStat_db(0), _thePop(0), rand(0){
#ifdef _DEBUG
    message(" TReplicate::TReplicate()\n");
#endif
    _pSimulation = ptr;
    init_paramset();
}

//----------------------------------------------------------------------------------------
// ~TReplicate
// ----------------------------------------------------------------------------------------
TReplicate::~TReplicate()
{
    if(_FileServices && _current_replicate) delete _FileServices; // don't delete the service of the first replicate here
    if(_StatServices) delete _StatServices;
    if(_thePop) delete _thePop;
}

//----------------------------------------------------------------------------------------
// init_paramset
// ----------------------------------------------------------------------------------------
void
TReplicate::init_paramset()
{
    // copy the paramSet of TSimulation
    _paramSet = _pSimulation->get_paramset();
    
}

// ----------------------------------------------------------------------------------------
// loadDefaultTemplates
// ----------------------------------------------------------------------------------------
/** load all traits and LCEs of a certain simulation */
void
TReplicate::loadDefaultTemplates(map<string, string>& inputParams)
{
#ifdef _DEBUG
    message("\nSimRunner::loadDefaultTemplates: ...\n");
#endif
    // generate all traits objects
    generate_traits(inputParams);
    
    // create the life cycle events (LCE)
    if (_isCoalescence) generate_LCE_coalescence();
    else generate_LCE_individual();
    
    // get all parameters
    build_allParams();
    
#ifdef _DEBUG
    message("\nSimRunner::loadDefaultTemplates: ... done\n");
#endif
}

// ----------------------------------------------------------------------------------------
// generate_traits
// ----------------------------------------------------------------------------------------
/** generate all trait objects (for the current number of trait!) */
void
TReplicate::generate_traits(map<string, string>& inputParams)
{
    assert(_TTrait_Templates.empty());
    
    //add TTraitProto
    unsigned int nb = get_nbTraits("quanti", inputParams);
    if (_isCoalescence && nb)
        error("Coalescence simulations do not allow to simulate quantitative traits!\n");
    switch (nb) {
        case 0:
            break;
        case 1:
            add_trait_template(new TTQuantiProto());
            break;
        default:
            for (unsigned int i = 0; i < nb; ++i) {
                add_trait_template(new TTQuantiProto(i + 1));
            }
            break;
    }
    
    // add neutral traits
    nb = get_nbTraits("ntrl", inputParams);
    switch (nb) {
        case 0:
            break;
        case 1:
            add_trait_template(new TTNeutralProto());
            break;
        default:
            for (unsigned int i = 0; i < nb; ++i) {
                add_trait_template(new TTNeutralProto(i + 1));
            }
            break;
    }
}

// ----------------------------------------------------------------------------------------
// generate_LCE_individual
// ----------------------------------------------------------------------------------------
/** specification of the life cycle */
void
TReplicate::generate_LCE_individual()
{
    assert(_LCE_Templates.empty());
    add_LCE_template(new LCE_Breed(1));
    add_LCE_template(new LCE_StatServiceNotifier(2));
    add_LCE_template(new LCE_FileServicesNotifier(3));
    add_LCE_template(new LCE_Aging(4));
    add_LCE_template(new LCE_Regulation(OFFSx, 5));
    add_LCE_template(new LCE_Disperse(6));
    add_LCE_template(new LCE_Regulation(ADLTx, 7));
    add_LCE_template(new LCE_Extinction(8));
}

// ----------------------------------------------------------------------------------------
// generate_LCE_coalescence
// ----------------------------------------------------------------------------------------
/** specification of the life cycle for the coalescence simulation */
void
TReplicate::generate_LCE_coalescence()
{
    assert(_LCE_Templates.empty());
    add_LCE_template(new LCE_Breed_coal(1));
    add_LCE_template(new LCE_Regulation(OFFSx, 2));
    add_LCE_template(new LCE_DisperseCoalescence(3));
    add_LCE_template(new LCE_store_popSizes(4));
    add_LCE_template(new LCE_Regulation(ADLTx, 5));
    add_LCE_template(new LCE_Extinction(6));
    
    // all these LCEs (rank>100) are executed after the "life"
    assert(_isCoalescence);
    add_LCE_template(new LCE_Coalescence_base(101));
    add_LCE_template(new LCE_StatServiceNotifier(102));
    add_LCE_template(new LCE_FileServicesNotifier(103));
}


// ----------------------------------------------------------------------------------------
/** get the number of traits of the given trait type */
unsigned int
TReplicate::get_nbTraits(string name, map<string, string>& inputParams)
{
    // if the number of loci is not specified no trait will be inizialized
    map<string, string>::iterator pos = inputParams.find(name + "_loci"); // general param
    if (pos == inputParams.end()) pos = inputParams.find(name + "_loci_1"); // or at least for the frist trait
    if (pos == inputParams.end()) return 0; // the number of loci is not specified
    if (pos == inputParams.end()) return 0; // the number of loci is not specified
    try {
        if (!strTo<unsigned int>(pos->second)) return 0;
    }
    catch(...) {
        error("The parameter '%s_loci' (%s) must be an integer!\n",
              name.c_str(), pos->second.c_str());
    }
    
    // get the number of traits
    pos = inputParams.find(name + "_nb_trait");
    if (pos == inputParams.end()) return 1; // the number of traits is not specified
    unsigned int nb;
    try {
        nb = strTo<unsigned int>(pos->second);
    }
    catch(...) {
        error(
              "The parameter '%s_nb_trait' must have a positive number as argument!\n",
              name.c_str());
    }
    return nb;                          // return the specified number of traits
}


// ----------------------------------------------------------------------------------------
/** is it a normal or a coalescent simulation */
unsigned int
TReplicate::isCoalescence(map<string, string>& inputParams)
{
    map<string, string>::iterator pos = inputParams.find("coalescence"); // general param
    if (pos == inputParams.end()) return false; // the parameter is not specified
    try {
        return strTo<bool>(pos->second);
    }
    catch(...) {
        error("The parameter 'coalescence' requires an integer as argument (%s)!\n",
              pos->second.c_str());
    }
    return false;
}

// ----------------------------------------------------------------------------------------
// run_replicate_thread
// ----------------------------------------------------------------------------------------
clock_t
TReplicate::print_start_replicate(unsigned int currentReplicate)
{
#ifdef _DEBUG
#ifdef _SHOW_MEMORY
    message("\n\r**** replicate %i/%i [%s] ****   RAM: %f MB\n",currentReplicate+1, _pSimulation->get_replicates()
            ,getElapsedTime(0).c_str(), process_mem_usage());
#else
    message("\n\r**** replicate %i/%i [%s] ****\n",currentReplicate+1,_pSimulation->get_replicates()
            ,getElapsedTime(0).c_str());
#endif
#else
#ifdef _SHOW_MEMORY
    message("\n\r    replicate %i/%i [%s]   RAM: %f MB",currentReplicate+1,_pSimulation->get_replicates()
            ,getElapsedTime(0).c_str(), process_mem_usage());
#else
    message("\n\r    replicate %i/%i [%s]",currentReplicate+1,_pSimulation->get_replicates()
            ,getElapsedTime(0).c_str());
#endif
#endif
    fflush(stdout);
    return clock();
}

// ----------------------------------------------------------------------------------------
// run_replicate
// ----------------------------------------------------------------------------------------
/** run a single replicate */
bool
TReplicate::run_replicate(map<string, string>& params, map<string, string>& keys,
                          unsigned int currentReplicate)
{
    // coalescence?
    _isCoalescence = isCoalescence(params);
    if(_isCoalescence) return run_replicate_coal(params, keys, currentReplicate);
    return run_replicate_ind(params, keys, currentReplicate);
    
}

// ----------------------------------------------------------------------------------------
// run_replicate
// ----------------------------------------------------------------------------------------
/** run a single replicate */
bool
TReplicate::run_replicate_ind(map<string, string>& params, map<string, string>& keys,
                          unsigned int currentReplicate)
{
    time_t start = time(0);
    
    // create metapop
    assert(!_thePop);
    _thePop = new Metapop(this, currentReplicate);
    _components.push_back(_thePop);
    _components.push_back(_thePop->get_protoGenome());
    
    
    // create all templates and set the parameters
    loadDefaultTemplates(params);
    
    //clear and build the lists of traits, LCEs, stats and files handlers, etc.
    if (!setup(params, keys, _thePop)) return false;
    
    _FileServices->set_params(_currentParams);
    
    //run the simulation
    //build initial populations for the current replicate (new adult generation)
    //reset_metapopulation();
    age_t age_required = _thePop->get_requiredAgeClass();
    if(!_thePop->setPopulation_FSTAT(age_required)) _thePop->setPopulation(age_required);
    if(!_thePop->size(age_required)) error("Metapopulation is empty: not able to start any simulation!\n");
    
    _thePop->executeBeforeEachReplicate(currentReplicate);
    //--------------------------- GENERATION LOOP ---------------------------
    
    _thePop->Loop_generation(start);
    
    //-----------------------------------------------------------------------
    
    _thePop->executeAfterEachReplicate(currentReplicate);
    
    _pSimulation->set_genLength(_thePop->getCurrentGeneration() - 1, currentReplicate);
    _pSimulation->set_elapsedTimeRepl(time(0) - start, currentReplicate);
    
    if(!_thePop->isAlive() && _thePop->getCurrentGeneration() != _thePop->getGenerations()) {
        _thePop->setCurrentGeneration(_thePop->getGenerations());
        if(_thePop->get_service()) _thePop->get_service()->notify();
    }
    
    // clean up
    reset_sim_components();
    _currentParams.clear();
    _currentTraits.clear();
    _currentLifeCycle.clear();
    delete _StatServices; _StatServices = NULL;
    delete _FileServices; _FileServices = NULL;
    
    return true;
}

// ----------------------------------------------------------------------------------------
// run_replicate
// ----------------------------------------------------------------------------------------
/** run a single replicate */
bool
TReplicate::run_replicate_coal(map<string, string>& params, map<string, string>& keys,
                          unsigned int currentReplicate)
{
    time_t start = time(0);
    
    // create metapop
    assert(!_thePop);
    _thePop = new Metapop(this, currentReplicate);
    _components.push_back(_thePop);
    _components.push_back(_thePop->get_protoGenome());
    
    // coalescence?
    _isCoalescence = isCoalescence(params);
    
    // create all templates and set the parameters
    loadDefaultTemplates(params);
    
    //clear and build the lists of traits, LCEs, stats and files handlers, etc.
    if (!setup(params, keys, _thePop)) return false;
    
    _FileServices->set_params(_currentParams);
    
    //run the simulation
    //build initial populations for the current replicate (new adult generation)
    //reset_metapopulation();
    _thePop->set_func_pointer_coal();
    _thePop->getCoalescence()->init(_thePop, _thePop->getGenerations());

    age_t age_required = _thePop->get_requiredAgeClass();
    _thePop->setPopulation_coal(age_required);
    if(!_thePop->size(age_required)) error("Metapopulation is empty: not able to start any simulation!\n");
    
    _thePop->executeBeforeEachReplicate(currentReplicate);
    //--------------------------- GENERATION LOOP ---------------------------
    
    _thePop->Loop_generation_coal(start);
    _thePop->set_samples_coalescence();
    _thePop->run_coalescence();
    
    //-----------------------------------------------------------------------
    
    _thePop->executeAfterEachReplicate(currentReplicate);
    
    _pSimulation->set_genLength(_thePop->getCurrentGeneration() - 1, currentReplicate);
    _pSimulation->set_elapsedTimeRepl(time(0) - start, currentReplicate);
    
    if(!_thePop->isAlive() && _thePop->getCurrentGeneration() != _thePop->getGenerations()) {
        _thePop->setCurrentGeneration(_thePop->getGenerations());
        if(_thePop->get_service()) _thePop->get_service()->notify();
    }
    
    // clean up
    reset_sim_components();
    _currentParams.clear();
    _currentTraits.clear();
    _currentLifeCycle.clear();
    delete _StatServices; _StatServices = NULL;
    delete _FileServices; _FileServices = NULL;
    
    return true;
}

// ----------------------------------------------------------------------------------------
// test_replicate_and_setUpStats
// ----------------------------------------------------------------------------------------
/** test all parameters and set up the stats */
bool
TReplicate::test_replicate_and_setUpStats(map<string, string>& params,
                                          map<string, string>& keys)
{
    // create metapop
    assert(!_thePop);
    _thePop = new Metapop(this, my_NAN);
    _components.push_back(_thePop);
    _components.push_back(_thePop->get_protoGenome());
    
    // coalescence?
    _isCoalescence = isCoalescence(params);
    
    // create all templates and set the parameters
    loadDefaultTemplates(params);
    
    //clear and build the lists of traits, LCEs, stats and files handlers, etc.
    if (!setup(params, keys, _thePop)) return false;
    _FileServices->set_params(_currentParams);
    
    return true;
}

// ----------------------------------------------------------------------------------------
// help
// ----------------------------------------------------------------------------------------
/** test all parameters and set up the stats */
void
TReplicate::print_help(ostream& os, unsigned int wide1, unsigned int wide2,
                       char fill, unsigned int importance, string arg)
{
    if(arg==""){
        os << "\n\nAvailable parameter sets:";
    }
    else{
        os << "\n\nParameters:\n";
        os << left << setw(wide1-2) << setfill(fill) << "\"--parameter";
        os << setfill(' ') << "  [data type] [range]  (default value)  {importance of parameter}\"\n";
        os << setw(wide1-2) << "\"" << "  description\"";
    }

    // simulation
    get_paramset().print_help(os, wide1, wide2, fill, importance, arg);
    
    Metapop pop(this,1);
    pop.get_paramset()->print_help(os, wide1, wide2, fill, importance, arg);
    
    TTQuantiProto  quanti;
    quanti.get_paramset()->print_help(os, wide1, wide2, fill, importance, arg);
    
    TTNeutralProto  ntrl;
    ntrl.get_paramset()->print_help(os, wide1, wide2, fill, importance, arg);
    
    LCE_Breed breed(1);
    breed.get_paramset()->print_help(os, wide1, wide2, fill, importance, arg);
    
    LCE_Aging age(1);
    age.get_paramset()->print_help(os, wide1, wide2, fill, importance, arg);
    
    TGenomeProto genome;
    genome.get_paramset()->print_help(os, wide1, wide2, fill, importance, arg);
    
    LCE_StatServiceNotifier stat(1);
    stat.get_paramset()->print_help(os, wide1, wide2, fill, importance, arg);
    
    LCE_FileServicesNotifier ser(1);
    ser.get_paramset()->print_help(os, wide1, wide2, fill, importance, arg);
    
    LCE_Regulation reg(OFFSx, 5);
    reg.get_paramset()->print_help(os, wide1, wide2, fill, importance, arg);
  
    LCE_Disperse disp(6);
    disp.get_paramset()->print_help(os, wide1, wide2, fill, importance, arg);
    
    LCE_Regulation reg2(ADLTx, 7);
    reg2.get_paramset()->print_help(os, wide1, wide2, fill, importance, arg);

    LCE_Extinction ext(8);
    ext.get_paramset()->print_help(os, wide1, wide2, fill, importance, arg);
    
    LCE_store_popSizes store(4);
    store.get_paramset()->print_help(os, wide1, wide2, fill, importance, arg);
    
    LCE_Coalescence_base coal(101);
    coal.get_paramset()->print_help(os, wide1, wide2, fill, importance, arg);
//    
//    LCE_CoalescenceRecomb coalR(101);
//    coalR.get_paramset()->print_help(os, wide1, wide2, fill, importance, arg);
//    
//    LCE_CoalescenceRecombII coalR2(101);
//    coalR2.get_paramset()->print_help(os, wide1, wide2, fill, importance, arg);
    
if(arg==""){
        os << "\n\nTo list the parameters of a parameter set please enter the " \
               "name of the parameter set, e.g";
        os << "\n         quantiNemo2 --help quanti";
        os << "\nor type";
        os << "\n         quantiNemo2 --help all";
        os << "\nto get the parameters of all parameter sets.";
    }
    os << "\n";
}

//----------------------------------------------------------------------------------------
// setup
// ----------------------------------------------------------------------------------------
bool
TReplicate::setup(map<string, string>& simparams, map<string, string>& simkeywords,
                  Metapop* thePop)
{
#ifdef _DEBUG
    message("SimRunner::run:building current params\n");
#endif
    
    //build the list of parameters from the record:
    if (!this->build_currentParams(simparams, simkeywords, this)) {
        error("SimRunner::run:couldn't build current params\n");
        return false;
    }
    
    //init the sim and pop params
    assert(!_FileServices);
    _FileServices = new FileServices;
    
    assert(!_StatServices);
    _StatServices = new StatServices;
    if (!build_pop(thePop)) return false;
    
    // check if the sequence of life cycle events makes sense
    checkLCEconsistency();
    
    //load the stats and files handlers
    register_all(thePop);
    
    // build the lists of stat recorders
    build_stat_recorders(thePop);
    
    print_info(thePop);
    
   return true;
}

//----------------------------------------------------------------------------------------
// register_all
// ----------------------------------------------------------------------------------------
//StatServices: build the lists of stat recorders: (cannot be done before the traits and LCEs are registered)
void
TReplicate::build_stat_recorders(Metapop* thePop)
{
    if(_pSimulation->stats && _StatServices->init()) return;
    
    // if no stats are computed remove this LCE
    delete _pSimulation->stats;
    _pSimulation->stats = NULL;
    
    
    if (thePop->isCoalescence()) {
        vector<LCE*>::iterator curLCE = thePop->get_theCycleCoal().begin(),
        endLCE = thePop->get_theCycleCoal().end();
        for (; curLCE != endLCE; ++curLCE) {
            if ((*curLCE)->get_event_name() != "save_stats") continue; // remove the stat LCE
            thePop->get_theCycleCoal().erase(curLCE);
            break;
        }
    }
    else {
        vector<LCE*>::iterator curLCE = thePop->get_theCycle().begin(),
        endLCE = thePop->get_theCycle().end();
        for (; curLCE != endLCE; ++curLCE) {
            if ((*curLCE)->get_event_name() != "save_stats") continue; // remove the stat LCE
            thePop->get_theCycle().erase(curLCE);
            break;
        }
    }
}

//----------------------------------------------------------------------------------------
// register_all
// ----------------------------------------------------------------------------------------
void
TReplicate::register_all(Metapop* thePop)
{
    // regsiter the metapop
    register_services(thePop);
    _FileServices->set_pop_ptr(thePop);
    _StatServices->set_pop_ptr(thePop);
    
    // register the traits
    map<string, TTraitProto*>::iterator trait = _currentTraits.begin();
    for (; trait != _currentTraits.end(); ++trait) {
        register_services(trait->second);
    }
    
    // register the life cycle events
    map<int, LCE*>::iterator iterLCE = _currentLifeCycle.begin();
    for (; iterLCE != _currentLifeCycle.end(); ++iterLCE) {
        register_services(iterLCE->second);
    }
}

//----------------------------------------------------------------------------------------
// print_info
// ----------------------------------------------------------------------------------------
void
TReplicate::print_info(Metapop* thePop)
{
    if(thePop->_current_replicate!=my_NAN) return;
    
    
    message("\nSETTINGS (%s)", STRING::get_file_name(_paramSet.getArg("_settings_file")).c_str());
    
    // simulation settings
    message("\n  Simulation:");
    switch (_isCoalescence) {
        case 0:
            message("\n    individual-based");
            break;
        case 1:
            message("\n    population-based (coalescence, speed optimization)");
            break;
        case 2:
            message("\n    population-based (coalescence, memory optimization)");
            break;
    }
    message("\n    %i generations", thePop->getGenerations());
    message("\n    %i replicates", _pSimulation->get_replicates());
    message("\n    %i. of %i batch simulations", 1+_pSimulation->_current_sim, _pSimulation->_nbSims);
    message("\n    %i threads used", _pSimulation->_threads);
    if (_pSimulation->get_fixedSeed()) message("\n    random generator initialized by the user");
    else message("\n    random generator initialized by the time");
    
    // traits
    message("\n\n  Loaded traits: ");
    unsigned int nb = 0;
    list<TTraitProto*>::iterator cur = _TTrait_Templates.begin(), end = _TTrait_Templates.end();
    for (; cur != end; ++cur) {
        if ((*cur)->get_paramset()->isSet()) { // _currentTraits cannot be used since there the traits are sorted alphabetically
            message("\n    %s", (*cur)->get_info().c_str());
            ++nb;
        }
    }
    if (!nb) message("\n    CAUTION: no traits (genetic data) are simulated!");
    
    // life cycle events
    message("\n\n  Life cycle sequence:");
    vector<LCE*>::iterator curLCE, endLCE;
    unsigned int i = 1;
    for (curLCE = thePop->get_theCycle().begin(), endLCE =
         thePop->get_theCycle().end(); curLCE != endLCE; ++curLCE) {
        message("\n    %i. %s", i++, (*curLCE)->get_event_name().c_str());
    }
    char c = 'A';
    for (curLCE = thePop->get_theCycleCoal().begin(), endLCE =
         thePop->get_theCycleCoal().end(); curLCE != endLCE; ++curLCE) {
        message("\n    %c. %s", c++, (*curLCE)->get_event_name().c_str());
    }
    
    // metapopulation
    unsigned int nbPatch = thePop->getPatchNbr(), nbPatchK, nbCarCap,
    nbSamplePatch, nbSample, nbIniPatch, nbIni;
    thePop->get_settingsStats(nbPatchK, nbCarCap, nbIniPatch, nbIni,
                              nbSamplePatch, nbSample);
    message("\n\n  Metapopulation:");
    if (nbPatchK == nbPatch) message("\n    %i patches (K_tot=%i)", nbPatch,
                                     nbCarCap);
    else message("\n    %i patches (K_tot=%i, nbPatch(K>0):%i)", nbPatch,
                 nbCarCap, nbPatchK);
    if (nbPatchK != nbIniPatch || nbCarCap != nbIni)
        message("\n    %i initial populations (N_tot=%i)", nbIniPatch, nbIni);
    if (nbPatchK != nbSamplePatch || nbCarCap != nbSample)
        message("\n    %i sampled patches (N_sample=%i)", nbSamplePatch,
                nbSample);
    message("\n    Migration model: %s",
            thePop->get_pDisperse_LCE()->get_disp_model_str().c_str());
    message("\n    Mating system: %s",
            thePop->get_pBreed_LCE()->getMatingSystem_str().c_str());
    
    // genetic map
    message("\n\n  Genetic map:");
    thePop->get_protoGenome()->printGeneticMapInfo();
    
    // output
    message("\n\n  Output:");
    list< FileHandler* >::iterator curFile = _FileServices->get_children().begin();
    list< FileHandler* >::iterator endFile = _FileServices->get_children().end();
    if (curFile == endFile) message("\n    CAUTION: no output is generated!");
    else if (_pSimulation->get_simfolder_short().empty()) message("\n    (to folder './')");
    else message("\n    (to folder '%s')", _pSimulation->get_simfolder_short().c_str());
    for (; curFile != endFile; ++curFile) {
        if ((*curFile)->get_trait().empty()) {            // that is not a trait
            if ((*curFile)->get_name() == "statistics") { // that is the stat file
                if(_StatServices->get_nb_params()){
                    message("\n    %s (%i parameters & %i statistics) computed at %i generations",
                            (*curFile)->get_name().c_str(),
                            _StatServices->get_nb_params(),
                            _StatServices->get_nb_stats(),
                            _pSimulation->stats->get_tot_occurrence());
                }
                else{
                    message("\n    %s (%i statistics) computed at %i generations",
                            (*curFile)->get_name().c_str(),
                            _StatServices->get_nb_stats(),
                            _pSimulation->stats->get_tot_occurrence());
                }
            }
            else {                                           // any other file
                message("\n    %s computed at %i generations",
                        (*curFile)->get_name().c_str(),
                        (*curFile)->get_tot_occurrence());
            }
        }
        else {                                   // all other files
            message("\n    %s (%s) computed at %i generations",
                    (*curFile)->get_name().c_str(),
                    (*curFile)->get_trait().front()->get_type().c_str(),
                    (*curFile)->get_tot_occurrence());
        }
    }
    
    message("\n");
}

//----------------------------------------------------------------------------------------
// register_services
// ----------------------------------------------------------------------------------------
void
TReplicate::register_services(SimComponent* cmpt)
{
    _FileServices->load(cmpt);
    _StatServices->load(cmpt);
}

//----------------------------------------------------------------------------------------
// build_pop
// ----------------------------------------------------------------------------------------
bool
TReplicate::build_pop(Metapop* thePop)
{
    if (_isCoalescence)
        return thePop->init_coal(build_currentTraits(), build_currentLifeCycle());
    return thePop->init(build_currentTraits(), build_currentLifeCycle());
}






