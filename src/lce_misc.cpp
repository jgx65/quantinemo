/** @file LCEmisc.cpp
 *
 *   Copyright (C) 2006 Frederic Guillaume    <guillaum@zoology.ubc.ca>
 *   Copyright (C) 2008 Samuel Neuenschwander <samuel.neuenschwander@unil.ch>
 *   Copyright (C) 2018 Frederic Michaud <frederic.michaud@unil.ch>
 *
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

#include "lce_misc.h"
#include "stathandler.cpp"

/* _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/ */

// ******** LCE_FileServicesNotifier ********

// ----------------------------------------------------------------------------------------
// LCE_FileServicesNotifier::execute
// ----------------------------------------------------------------------------------------
void LCE_FileServicesNotifier::execute() {
#ifdef _DEBUG
    message("  LCE_FileServicesNotifier ... ");
#endif
    _service->notify();
    
#ifdef _DEBUG
    message("done! (gen: %i rpl: %i)\n",
            this->get_pop_ptr()->getCurrentGeneration()+1,
            this->get_pop_ptr()->getCurrentReplicate()+1);
#endif
}

/* _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/ */

// ******** LCE_StatServiceNotifier ********

// -----------------------------------------------------------------------------
// LCE_StatServiceNotifier
// -----------------------------------------------------------------------------
LCE_StatServiceNotifier::LCE_StatServiceNotifier(int rank) :
LCE("save_stats", "statistics", "", rank), _service(0) {
    
    add_parameter("stat", STR_MAT, true, my_NAN, my_NAN, "",false,
                        "List of statistics to compute.",0);
    
    add_parameter("param",STR_MAT, false,my_NAN, my_NAN, "", false,
                        "List of parameter arguments to output in the statistic file. " \
                        "Note, the output is listed only in the detailed statistic file.",0);
    
    add_parameter("stat_save", INT2, false, 0, 7, "0",false,
                        "How should the statistics be stored:\n" \
                        "  0: all (detailed, mean and variance across replicates)\n" \
                        "  1: detailed (detailed output per replicate)\n" \
                        "  2: summed up (mean and variance across replicates)\n" \
                        "  3: mean (mean across replicates)\n" \
                        "  4: variance (variance across replicates)\n" \
                        "  5: median (median across replicates)\n" \
                        "  6: none (no statistics are reported)\n",0);
    
    add_parameter("stat_log_time", INT2, false, 0, my_NAN, "1", true,
                        "The time interval at which statistics are recorded.",0);
    
    add_parameter("stat_dir", STR, false, my_NAN, my_NAN, "",false,
                        "Directory to store the statistic files.", 1);
    
    add_parameter("stat_filename", STR, false, my_NAN, my_NAN, "", false,
                        "Base file name of the statistic files.", 1);
}

// -----------------------------------------------------------------------------
// loadStatServices
// -----------------------------------------------------------------------------
void LCE_StatServiceNotifier::loadStatServices(StatServices* loader)
{
    _service = loader;
    if(_popPtr->get_pStat_db())_fileHandler.set_statService(loader);
    _service->set_statArg(_stat_arg);
    _service->set_paramArg(_param_arg);
    if(_popPtr->get_pStat_db()) _fileHandler.set_save_choice(_save_choiceTemp);
    loader->attach(&_statHandler);
}

// -----------------------------------------------------------------------------
// LCE_StatServiceNotifier::execute
// -----------------------------------------------------------------------------
void LCE_StatServiceNotifier::execute()
{
#ifdef _DEBUG
    message("  LCE_StatServiceNotifier ... ");
#endif
    
    assert(_service);
    if (!(_popPtr->getCurrentGeneration() % _pStat_db->get_occurrence())){
        _service->update_states();
        _service->notify_params();
        _service->notify();
    }
    
#ifdef _DEBUG
    message("done! (occurrence %i)\n", _pStat_db->get_occurrence());
#endif
}

// -----------------------------------------------------------------------------
// init
// -----------------------------------------------------------------------------
bool LCE_StatServiceNotifier::init(TMetapop* popPtr) {
    LCE::init(popPtr);
    
    _pStat_db = get_pop_ptr()->get_pStat_db();
    
    _stat_arg = _paramSet->getArg("stat");
    _param_arg = _paramSet->getArg("param");
    _save_choiceTemp = (unsigned int)_paramSet->getValue("stat_save");
    
    if (_save_choiceTemp == 6) { // no stats are written, LCE is not used
        _paramSet->set_isSet(false);
        delete get_pop_ptr()->get_pStat_db();
        _popPtr->set_pStat_db(NULL);
        return false;
    }
    
    // get the number of logtime occurrences
    unsigned int gen = _popPtr->getGenerations();
    vector<unsigned int> vIndex = get_logtime_occurences(_paramSet->get_param("stat_log_time"), gen);
    
    // if stats are never computed
    if (vIndex.empty()) {
        warning("The log time of the summary statistics is larger than the number of generations (%i): no summary statistics will be computed!\n", gen);
        _paramSet->set_isSet(false);
        return false;
    }

    // copy vector to array
    if(popPtr->_current_replicate==my_NAN){
        _pStat_db->set_index(vIndex);
    }
    
    _fileHandler.set(NULL, gen, // is called only at the last generation
                     _paramSet->getArg("stat_dir"),
                     _paramSet->getArg("stat_filename"),
                     "", 0, 0, 0, NULL, "statistics", ".txt", NULL, get_pop_ptr());
    return true;
}


// ----------------------------------------------------------------------------------------
// set_logtime_occurences
// ----------------------------------------------------------------------------------------
/** return a vector of the loged time for the given logtime */
vector<unsigned int>
LCE_StatServiceNotifier::get_logtime_occurences(Param* pParam, unsigned int totGen)
{
    vector<unsigned int> vIndex;
    if(pParam->isTemporalParam()){
        map<int, string>* pArgs = pParam->get_temporal_args();
        map<int, string>::iterator pos = pArgs->begin();  // set the iterator
        
        // loop through the generations and sum up the occurences
        for(unsigned int occ, i=1; i<=totGen; ++i){
            assert(i!=1 || (i==1 && i==(unsigned int) pos->first)); // first one has to be set
            if(i == (unsigned int) pos->first){                      // new temporal parameter
                occ = strTo<unsigned int>(pos->second);     // change occurence
                ++pos;                                  // go the the next position
            }
            if(!((i) % occ)) vIndex.push_back(i); //++ _tot_occurrence;
        }
    }
    else {     // no temporal change
        unsigned int occ = (unsigned int) pParam->get_value();
        for(unsigned int i=1; i<=totGen; ++i){
            if(!((i) % occ)) vIndex.push_back(i); //++ _tot_occurrence;
        }
    }
    
    // coalescence simulations
    if (_popPtr->isCoalescence() && (vIndex.size() != 1 || vIndex[0] != totGen)) {
        warning("Coaelscence: statistics may only be computed for the last generation: logtime adjusted!\n");
        vIndex.clear(); // simplest way to achieve it, although not the most effective...
        vIndex.push_back(totGen);
    }
    
    return vIndex;
}


// -----------------------------------------------------------------------------
// temporal_change
// -----------------------------------------------------------------------------
void LCE_StatServiceNotifier::temporal_change(const unsigned int& gen)
{
    // temporal parameters
    map<string, Param*> *pParam = _paramSet->getTemporalParams(gen);

    if (pParam) {
        // check if a change has to be made
        map<string, Param*> *pMap = _paramSet->updateTemporalParams(gen);
        if (pMap) {
            // iterate through the map and performe the updates
            map<string, Param*>::iterator pos = pMap->begin();
            for (; pos != pMap->end(); ++pos) {
                if (pos->first == "stat_log_time" && _service) {
                    _pStat_db->set_occurrence((unsigned int)pos->second->get_value());
                }
            }
        }
    }
}

// -----------------------------------------------------------------------------
// get_isAlive
// -----------------------------------------------------------------------------
double LCE_StatSH::get_isAlive() {
    return(double)_popPtr->isAlive();
}

/* _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/ */

// ******** LCE_Aging ********

// ----------------------------------------------------------------------------------------
// LCE_Aging::execute
// ----------------------------------------------------------------------------------------
/** everything moves by a generation. The oldest age class will be flushed.
 */
void LCE_Aging::execute() {
#ifdef _DEBUG
    message("  LCE_Aging ... ");
#endif
    
    // remove all adults
    vector<TPatch*>::iterator curPop, endPop = _popPtr->get_vFullPatch().end();
    for (curPop = _popPtr->get_vFullPatch().begin(); curPop != endPop; ) {
        (*curPop)->flush(ADLTx); // remove all adults
        
        // check if the patch is still populated:
        if ((*curPop)->size(OFFSx)) ++curPop;
        else _popPtr->new_emptyPatch(curPop, endPop);   // remove if now empty (curPop and endPop are adjusted)
    }
    
#ifdef _DEBUG
    message("done! (Patches: %i; Individuals(F/M): [%i/%i; %i/%i])\n",
            _popPtr->get_nbFullPatch(), _popPtr->size(FEM, OFFSx),
            _popPtr->size(MAL, OFFSx), _popPtr->size(FEM, ADLTx),
            _popPtr->size(MAL, ADLTx));
    assert(_popPtr->individual_container_ok());
#endif
}

/* _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/ */

// ******** LCE_store_popSizes ********

// -----------------------------------------------------------------------------
// LCE_store_popSizes::execute
// -----------------------------------------------------------------------------
void LCE_store_popSizes::execute() {
#ifdef _DEBUG
    message("  LCE_store_popSizes ... ");
#endif
    
    _popPtr->store_popSizes();
    
#ifdef _DEBUG
    message("done!\n");
#endif
}

// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// LCE_store_popSizes::init
// -----------------------------------------------------------------------------
bool LCE_store_popSizes::init(TMetapop* popPtr) {
    LCE::init(popPtr);
    return true;
}



