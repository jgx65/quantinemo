/** @file lce_coalescence.cpp
 *
 *   Copyright (C) 2010 Samuel Neuenschwander <samuel.neuenschwander@unil.ch>
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
// ---------------------------------------------------------------------------
#include "lce_coalescence.h"
#include "lce_coalescence_base.h"
//#include <algorithm>
//#include "stathandler.cpp"
#include "lce_disperse.h"
//#include <limits>

using namespace std;

// ----------------------------------------------------------------------------------------
double LCE_Coalescence::get_model_threshold(){assert(_pCoalBase); return _pCoalBase->get_model_threshold();}
TMetapop* LCE_Coalescence::get_pop_ptr(){assert(_pCoalBase); return _pCoalBase->get_pop_ptr();}


// ----------------------------------------------------------------------------------------
/** initialize the class */
void LCE_Coalescence::init(LCE_Coalescence_base* pBase)
{
    _pCoalBase = pBase;
    _model_threshold = _pCoalBase->get_model_threshold();
    _dbCoalescence = _pCoalBase->get_dbCoalescence();
    _seq = _pCoalBase->get_seq();
    
    switch (get_pop_ptr()->get_pDisperse_LCE()->get_disp_model()) {
        case 0: // island model
        case 1: // island propagule model
            func_ptr_evacuate_lineages = &LCE_Coalescence::evacuate_lineages_island;
            break;
        case 2: // 1D stepping stone model
        case 3: // 2D stepping stone model
            func_ptr_evacuate_lineages = &LCE_Coalescence::evacuate_lineages_stepping_stone;
            break;
    }
    
    set_dispRates();
    if (_model_threshold == 0) func_ptr_perform_coalescence =
        &LCE_Coalescence::perform_coalescence_multipleEvent;
    else if (_model_threshold == 1e6) func_ptr_perform_coalescence =
        &LCE_Coalescence::perform_coalescence_singleEvent;
    else func_ptr_perform_coalescence = &LCE_Coalescence::perform_coalescence_mixedEvent;

}

// ----------------------------------------------------------------------------------------
/** add the current pop sizes to the db for coalescence simulations, but only if it is not empty */
bool LCE_Coalescence::check_db(const unsigned int& gen)
{
	return true; //_dbTemp.check_db();
}

// ----------------------------------------------------------------------------------------
/** set the initial deme sizes depending on the db (last entry in the db)
 */
void LCE_Coalescence::set_ini_demeSizes()
{
	_demes.clear();

	map<PATCH_ID, dbPop>::iterator curPos, endPos;
	curPos = _dbCoalescence[_pCoalBase->get_nbGen() - 1].get_patches().begin();
	endPos = _dbCoalescence[_pCoalBase->get_nbGen() - 1].get_patches().end();
	for (; curPos != endPos; ++curPos) {
		assert(curPos->second.get_popSize());
		_demes[(unsigned int) curPos->first]->set_deme_size(
				(unsigned int) curPos->second.get_popSize(), (unsigned int) curPos->first);
	}
}

// ----------------------------------------------------------------------------------------
/** set the sample sizes (lineages) AND the pop sizes of these demes.
 * Caution: sample size is diploid, lineages is haploid!!!
 * before this function is called the functions
 *  - set_samples(const map<unsigned int, unsigned int>& m, const unsigned int)
 *  - set_ini_demeSizes()
 * have to be called!!!
 * sample sizes are already adjusted to the current pop sizes
 */
void LCE_Coalescence::set_samples()
{
	assert(!_pCoalBase->get_iniSamples().empty());

	_demes.clear();
	clear_tree();
	ARRAY::create_1D<TTNode*>(_sampleNodes, _pCoalBase->get_nbSamples());

	map<PATCH_ID, dbPop>& db = _dbCoalescence[_pCoalBase->get_nbGen() - 1].get_patches(); // the database with the pop sizes
	unsigned int patchID;
	TDeme* curDeme;
	_nbNodes=0;
	map<unsigned int, unsigned int>::iterator curSample, endSample = _pCoalBase->get_iniSamples().end(); // the sample sizes
	for (curSample = _pCoalBase->get_iniSamples().begin(); curSample != endSample; ++curSample) {
		patchID = (unsigned int) curSample->first;
		curDeme = get_demes(patchID);
		curDeme->set_sample_size(curSample->second, _sampleNodes, _nbNodes, patchID);
		assert(db.find(patchID)!=db.end());
		assert((unsigned int)db[patchID].get_popSize()>=curSample->second);
		curDeme->set_deme_size((unsigned int) db[patchID].get_popSize(), patchID);
	}
}

// ----------------------------------------------------------------------------------------
bool LCE_Coalescence::allPopulated()
{
    map<PATCH_ID, TDeme*>::iterator curDeme, endDeme;
    for (curDeme = _demes.begin(), endDeme = _demes.end(); curDeme != endDeme; ++curDeme) {
        if (!curDeme->second->get_deme_size()) return false;
    }
    return true;
}

// ----------------------------------------------------------------------------------------
/** builds the three efficiently by just adjusting the deme sizes of the demes containing
 * lineages. Thus it is not possible to store the lineages_file!!!
 */
unsigned long LCE_Coalescence::build_tree_simple()
{
    // Initialize the demes
    set_samples();
    
    _num_lineages = _pCoalBase->get_nbSamples();
    _timeDivergence = (unsigned long) _pCoalBase->get_parameter("divergence_time")->get_value(); // get the divergence time
    if(_timeDivergence<_pCoalBase->get_nbGen()) _timeDivergence = _pCoalBase->get_nbGen(); // divergence is like a migration event after coalescence event
    
#ifdef _DEBUG
    message("\n  Coalescence simulations: %i lineages:", _num_lineages);
#endif
    
    // build the tree over the time period of the demographic simulation
    dbCoalescence* curDB;
    for (_curGen = 1; _curGen <= (unsigned long) _pCoalBase->get_nbGen(); ++_curGen) { // for each generation (starting with 1) in the db going backward in time
        curDB = &_dbCoalescence[_pCoalBase->get_nbGen() - _curGen]; // starting from the back
        perform_popSize_regulation(curDB);     // Step 1: adapt population sizes
        perform_coalescence();                      // Step 3: Coalescence round
        if(_num_lineages==1) break;
        perform_migration(curDB);                  	// Step 2: migration
    }
    
    // if the MRCA has not yet been found continue searching it
    if (_num_lineages > 1) {
#ifdef _DEBUG
        message("\n    At starting time (%u) there are %i lineages remaining.", _curGen, _num_lineages);
#endif
        
        // build the tree over the time period before the demographic simulation
        for (; _curGen < _timeDivergence; ++_curGen) {	// for each generation (starting with 1) in the db going backward in time
            perform_coalescence();                  // Step 3: Coalescence round
            if(_num_lineages==1) break;
            perform_migration_before();                 // Step 2: migration
        }
        
        // if the MRCA has not yet been found continue searching it
        if (_num_lineages > 1) {
#ifdef _DEBUG
            message("\n    At divergence time (%u) there are %i lineages remaining.", _curGen, _num_lineages);
#endif
            
            // divergence is reached: merge all demes, if the new pop size is not specified the pop
            // size will be the total number of individuals of the metapop
            merge_demes();
            // continue until the last gene
            for (; _curGen < ULONG_MAX && _num_lineages > 1; ++_curGen) {// for each generation (starting with 1) in the db going backward in time
                perform_coalescence();              // Step 3: Coalescence round
            }
        }
    }
    
    if (_num_lineages > 1) error("Coalescence:: Coalescence process did not converge!\n");
    
#ifdef _DEBUG
    message("\n    MRCA reached at time (%u).", _curGen);
#endif
    
    _tmrca = _mrca->time;		// get the time to the MRCA
    _pCoalBase->set_vMRCA(_curLocus, _tmrca);
    
    // remove the last lineage (the node is already in the three)
    remove_mrca_from_deme();
    
    return _curGen;
}

// ----------------------------------------------------------------------------------------
bool LCE_Coalescence::curDemesAsDB(dbCoalescence* curDBs)
{
    /*
     // are the demes of db present in _demes?
     vector<pair<PATCH_ID, POP_SIZE> >::iterator curDB, endDB;
     for(curDB=curDBs.get_popSize().begin(), endDB=curDBs.get_popSize().end(); curDB!=endDB; ++curDB){
     if(_demes.find((unsigned int)curDB->first) == _demes.end()) return false;
     }
     
     
     // are the demes of _demes present in the db?
     map<PATCH_ID, TDeme*>::iterator curDeme, endDeme;
     for(curDeme=_demes.begin(), endDeme=_demes.end(); curDeme!=endDeme; ++curDeme){
     if(curDBs.get_popSize().find((PATCH_ID)curDeme->first) == curDBs.get_popSize().end()) return false;
     }
     */
    return true;
}

// ----------------------------------------------------------------------------------------
void LCE_Coalescence::plotDemes()
{
    map<PATCH_ID, TDeme*>::iterator curDeme, endDeme;
    cout << "\n*** ";
    for (curDeme = _demes.begin(), endDeme = _demes.end(); curDeme != endDeme; ++curDeme) {
        if (curDeme->second->get_lineages())
            cout << curDeme->first << " (" << curDeme->second->get_deme_size() << "; "
            << curDeme->second->get_lineages() << "); ";
    }
}

// ------------------------------------------------------------------------------
/** merge all remaining lineages in a single (first) deme
 * the new pop size may be specified by the parameter "divergence_pop_size"
 * if not specified the new pop size is the sum of all individuals of the metapop at time 0
 * if the new pop size is smaller than the number of lineages, the pop size is readjusted
 */
void LCE_Coalescence::merge_demes()
{
    assert(_curGen==_timeDivergence);
    
    // merge the lineages (nodes) of all demes
    map<PATCH_ID, TDeme*>::iterator curDeme, endDeme;
    curDeme = _demes.begin();
    endDeme = _demes.end();
    TDeme* firstDeme = curDeme->second;
    for (++curDeme; curDeme != endDeme; ++curDeme) { // we start at the second deme...
        assert(curDeme->second->get_lineages());
        firstDeme->get_chainNodeList().insert(
                                              firstDeme->get_chainNodeList().end(), // add the lineages to the first deme
                                              curDeme->second->get_chainNodeList().begin(),
                                              curDeme->second->get_chainNodeList().end());
        curDeme->second->clear_chainNodeList(); // remove the lineages for safety reasons
    }
    firstDeme->update_lineages(); // set the number of lineages
    
    // remove all demes, but the first one
    curDeme = _demes.begin();
    _demes.erase(++curDeme, _demes.end());
    
    // adapt the pop size of the first deme (number of lineages, deme size)
    unsigned int popSize = (unsigned int) _pCoalBase->get_parameter("divergence_pop_size")->get_value();
    if (popSize) {                                      // pop size is specified
        if (2 * popSize < firstDeme->get_lineages()) firstDeme->set_deme_size(
                                                                              (unsigned int) ((double) firstDeme->get_lineages() / 2)); // pop size is too small
        else firstDeme->set_deme_size(popSize);
    }
    else {                // pop size is not specified => adjust to metapop size
        // get the pop size of the entire metapopulation at generation 1
        unsigned int sumPopSize = 0;
        map<PATCH_ID, dbPop>::iterator curDeme, endDeme = _dbCoalescence[0].get_patches().end(); // the database with the pop sizes
        for (curDeme = _dbCoalescence[0].get_patches().begin(); curDeme != endDeme; ++curDeme) {
            sumPopSize += curDeme->second.get_popSize();
        }
        firstDeme->set_deme_size(sumPopSize); // merged pop size is the sum of all pop sizes at t=0
    }
    assert(2*firstDeme->get_deme_size() >= firstDeme->get_lineages());
}

// ----------------------------------------------------------------------------------------
void LCE_Coalescence::plotMigrDB(map<unsigned int, unsigned int>& curMigr)
{
    map<unsigned int, unsigned int>::iterator curDeme, endDeme;
    cout << "\n*** ";
    for (curDeme = curMigr.begin(), endDeme = curMigr.end(); curDeme != endDeme; ++curDeme) {
        cout << curDeme->first << " (" << curDeme->second << "); ";
    }
}

// ----------------------------------------------------------------------------------------
/** perform migration
 * Caution "to" and "from" are in forward time sense!!!
 */
void LCE_Coalescence::perform_migration(dbCoalescence* DB)
{
    // for each deme with lineages: migrate lineages
    // if a lineage arrives in a new deme this deme is added to temp vector and
    // only added after all migrations to _demes, since otherwise the iterator is invalidated
    map<PATCH_ID, TDeme*> newDemes; //temp container to take the newly "colonized by lineages" demes
    map<PATCH_ID, TDeme*>::iterator curDeme, endDeme;
    //	map<PATCH_ID, dbPop>::iterator curDB, endDB = DB->get_patches().end();
    for (curDeme = _demes.begin(), endDeme = _demes.end(); curDeme != endDeme; ++curDeme) {
        assert(DB->get_patches().find(curDeme->first)!=DB->get_patches().end()); // the element has to exist
        dbPop& dbDeme = DB->get_patches()[curDeme->first]; // get the current deme
        if (dbDeme.get_immigrants().empty()) continue; // if no migration happened
        migration(curDeme->second, dbDeme.get_immigrants(), curDeme->second->get_deme_size(),
                  newDemes);
    }
    
    if (!newDemes.empty()){
        _demes.insert(newDemes.begin(), newDemes.end()); // merge the containers
    }
    
    // for each deme with lineages: update the number of lineages
    // if a deme does not anymore contain lineages it is removed from _demes
    // endDeme is not anymore the same as above, and may vary during looping!!!
    for (curDeme = _demes.begin(); curDeme != _demes.end();) {
        if (curDeme->second->post_migration_lineages()) ++curDeme; // update the number of lineages
        else {
            delete curDeme->second;
            _demes.erase(curDeme++);
        }
    }
}

// ----------------------------------------------------------------------------------------
/** perform migration
 * Caution "to" and "from" are in forward time sense!!!
 * it is possible that lineages migrate (backward in time) to a currently empty deme
 */
void LCE_Coalescence::migration(TDeme* curDeme, vector<pair<PATCH_ID, POP_SIZE> >& curMigr,
                                const unsigned int& demeSize, vector<TDeme*>& newDemes)
{
    assert(curDeme->get_lineages()); // get the number of lineages which may be send
    
    // get the total number of immigrants
    unsigned int totMigr = 0;
    vector<pair<PATCH_ID, POP_SIZE> >::iterator curFrom, endFrom;
    for (curFrom = curMigr.begin(), endFrom = curMigr.end(); curFrom != endFrom; ++curFrom) {
        totMigr += (unsigned int) curFrom->second; // the from deme may be not existent at the moment!
    }
    
    // get the number of lineages to migrate
    assert(demeSize);
    assert(((double)totMigr)/demeSize>=0);
    assert(((double)totMigr)/demeSize<=1);
    unsigned int migrLineages = get_pop_ptr()->rand().Binomial(((double) totMigr) / demeSize,
                                                               curDeme->get_lineages()); // draw the number of lineages to migrate
    if (!migrLineages) return; // although migration happened no lineages will migrate
    
    // send lineages to each neighbor (multinomial...)
    double sum_m = 1.0, cur_m;
    unsigned int nbEmigrants;
    for (curFrom = curMigr.begin(), endFrom = curMigr.end(); migrLineages; ++curFrom) {
        cur_m = (double) curFrom->second / totMigr; // get migration rate of the current neighbor
        nbEmigrants = get_pop_ptr()->rand().Binomial(cur_m / sum_m, migrLineages); // get number of lineages to emigrate to this neighbour
        sum_m -= cur_m;                     // adjust the sum of migration rates
        if (nbEmigrants) {
            curDeme->migrate(curFrom->first, nbEmigrants, *this, newDemes); // perform migration (maybe the deme has to be created)!!!
            migrLineages -= nbEmigrants;
        }
    }
}

// ----------------------------------------------------------------------------------------
/** perform migration
 * Caution "to" and "from" are in forward time sense!!!
 * it is possible that lineages migrate (backward in time) to a currently empty deme
 */
void LCE_Coalescence::migration(TDeme* curDeme, vector<pair<PATCH_ID, POP_SIZE> >& curMigr,
                                const unsigned int& demeSize, map<PATCH_ID, TDeme*>& newDemes)
{
    assert(curDeme->get_lineages()); // get the number of lineages which may be send
    
    // get the total number of immigrants
    unsigned int totMigr = 0;
    vector<pair<PATCH_ID, POP_SIZE> >::iterator curFrom, endFrom;
    for (curFrom = curMigr.begin(), endFrom = curMigr.end(); curFrom != endFrom; ++curFrom) {
        totMigr += (unsigned int) curFrom->second; // the from deme may be not existent at the moment!
    }
    
    // get the number of lineages to migrate
    assert(demeSize);
    assert(((double)totMigr)/demeSize>=0);
    assert(((double)totMigr)/demeSize<=1);
    unsigned int migrLineages = get_pop_ptr()->rand().Binomial(((double) totMigr) / demeSize,
                                                               curDeme->get_lineages()); // draw the number of lineages to migrate
    if (!migrLineages) return; // although migration happened no lineages will migrate
    
    // send lineages to each neighbor (multinomial...)
    // store them in a separate container and merge them only after a round
    double sum_m = 1.0, cur_m;
    unsigned int nbEmigrants;
    for (curFrom = curMigr.begin(), endFrom = curMigr.end(); migrLineages; ++curFrom) {
        cur_m = (double) curFrom->second / totMigr; // get migration rate of the current neighbor
        nbEmigrants = get_pop_ptr()->rand().Binomial(cur_m / sum_m, migrLineages); // get number of lineages to emigrate to this neighbour
        sum_m -= cur_m;                     // adjust the sum of migration rates
        if (nbEmigrants) {
            curDeme->migrate(curFrom->first, nbEmigrants, *this, _demes, newDemes); // perform migration (maybe the deme has to be created)!!!
            migrLineages -= nbEmigrants;
        }
    }
}

// ----------------------------------------------------------------------------------------
/** perform migration before the onset of demographic simulation until divergence
 * Caution "to" and "from" are in forward time sense!!!
 */
void LCE_Coalescence::perform_migration_before()
{
    // for each deme with lineages: migrate lineages
    // if a lineage arrives in a new deme this deme is added to temp vector and
    // only added after all migrations to _demes, since otherwise the iterator is invalidated
    map<PATCH_ID, TDeme*> newDemes; //temp container to take the newly "colonized by lineages" demes
    map<PATCH_ID, TDeme*>::iterator curDeme, endDeme;
    map<PATCH_ID, vector<pair<PATCH_ID, double> > >::iterator pos;
    for (curDeme = _demes.begin(), endDeme = _demes.end(); curDeme != endDeme; ++curDeme) {
        pos = _dispRates.find(curDeme->first);           // get the current deme
        if (pos == _dispRates.end()) continue;       // if no migration happened
        (this->*func_ptr_migration_before)(curDeme->second, pos->second,
                                           curDeme->second->get_deme_size(), newDemes);
    }
    
    if (!newDemes.empty()) _demes.insert(newDemes.begin(), newDemes.end()); // merge the containers
    
    // for each deme with lineages: update the number of lineages
    // if a deme does not anymore contain lineages it is removed from _demes
    // endDeme is not anymore the same as above, and may vary during looping!!!
    for (curDeme = _demes.begin(); curDeme != _demes.end();) {
        if (curDeme->second->post_migration_lineages()) ++curDeme; // update the number of lineages
        else {
            delete curDeme->second;
            _demes.erase(curDeme++);
        }
    }
}

// ----------------------------------------------------------------------------------------
/** perform migration
 * Caution "to" and "from" are in forward time sense!!!
 * it is possible that lineages migrate (backward in time) to a currently empty deme
 */
void LCE_Coalescence::migration_before_island(TDeme* curDeme,
                                              vector<pair<PATCH_ID, double> >& curMigr, const unsigned int& demeSize,
                                              map<PATCH_ID, TDeme*>& newDemes)
{
    assert(curDeme->get_lineages()); // get the number of lineages which may be send
    
    // get the total emigration rate
    double totMigr = 0;
    vector<pair<PATCH_ID, double> >::iterator curFrom, endFrom;
    for (curFrom = curMigr.begin(), endFrom = curMigr.end(); curFrom != endFrom; ++curFrom) {
        totMigr += curFrom->second; // the from deme may be not existent at the moment!
    }
    
    // get the number of lineages to migrate
    assert(demeSize);
    assert(totMigr>=0);
    assert(totMigr<=1);
    unsigned int migrLineages = get_pop_ptr()->rand().Binomial(totMigr, curDeme->get_lineages()); // draw the number of lineages to migrate
    if (!migrLineages) return; // although migration happened no lineages will migrate
    
    // send lineages to each neighbor (multinomial...)
    double sum_m = 1.0, cur_m;
    unsigned int nbEmigrants;
    for (curFrom = curMigr.begin(), endFrom = curMigr.end(); migrLineages; ++curFrom) {
        cur_m = curFrom->second;   // get migration rate of the current neighbor
        nbEmigrants = get_pop_ptr()->rand().Binomial(cur_m / sum_m, migrLineages); // get number of lineages to emigrate to this neighbour
        sum_m -= cur_m;                     // adjust the sum of migration rates
        if (nbEmigrants) {
            curDeme->migrate(curFrom->first, nbEmigrants, *this, _demes, newDemes); // perform migration (maybe the deme has to be created)!!!
            migrLineages -= nbEmigrants;
        }
    }
}

// ----------------------------------------------------------------------------------------
/** perform migration
 * Caution "to" and "from" are in forward time sense!!!
 * it is possible that lineages migrate (backward in time) to a currently empty deme
 */
void LCE_Coalescence::migration_before_propagule(TDeme* curDeme,
                                                 vector<pair<PATCH_ID, double> >& curMigr, const unsigned int& demeSize,
                                                 map<PATCH_ID, TDeme*>& newDemes)
{
    assert(curDeme->get_lineages()); // get the number of lineages which may be send
    
    // get the total emigration rate
    double totMigr = 0;
    vector<pair<PATCH_ID, double> >::iterator curFrom, endFrom;
    for (curFrom = curMigr.begin(), endFrom = curMigr.end(); curFrom != endFrom; ++curFrom) {
        totMigr += curFrom->second; // the from deme may be not existent at the moment!
    }
    
    // get the number of lineages to migrate
    assert(demeSize);
    assert(totMigr>=0);
    assert(totMigr<=1);
    unsigned int migrLineages = get_pop_ptr()->rand().Binomial(totMigr, curDeme->get_lineages()); // draw the number of lineages to migrate
    if (!migrLineages) return; // although migration happened no lineages will migrate
    
    // send lineages to each neighbor (multinomial...)
    double sum_m = 1.0, cur_m;
    unsigned int nbEmigrants;
    for (curFrom = curMigr.begin(), endFrom = curMigr.end(); migrLineages; ++curFrom) {
        cur_m = curFrom->second;   // get migration rate of the current neighbor
        nbEmigrants = get_pop_ptr()->rand().Binomial(cur_m / sum_m, migrLineages); // get number of lineages to emigrate to this neighbour
        sum_m -= cur_m;                     // adjust the sum of migration rates
        if (nbEmigrants) {
            curDeme->migrate(curFrom->first, nbEmigrants, *this, _demes, newDemes); // perform migration (maybe the deme has to be created)!!!
            migrLineages -= nbEmigrants;
        }
    }
}

// ----------------------------------------------------------------------------------------
/** perform migration
 * Caution "to" and "from" are in forward time sense!!!
 * it is possible that lineages migrate (backward in time) to a currently empty deme
 */
void LCE_Coalescence::migration_before_matrix(TDeme* curDeme,
                                              vector<pair<PATCH_ID, double> >& curMigr, const unsigned int& demeSize,
                                              map<PATCH_ID, TDeme*>& newDemes)
{
    assert(curDeme->get_lineages()); // get the number of lineages which may be send
    
    // get the total emigration rate
    double totMigr = 0;
    vector<pair<PATCH_ID, double> >::iterator curFrom, endFrom;
    for (curFrom = curMigr.begin(), endFrom = curMigr.end(); curFrom != endFrom; ++curFrom) {
        totMigr += curFrom->second; // the from deme may be not existent at the moment!
    }
    
    // get the number of lineages to migrate
    assert(demeSize);
    assert(totMigr>=0);
    assert(totMigr<=1);
    unsigned int migrLineages = get_pop_ptr()->rand().Binomial(totMigr, curDeme->get_lineages()); // draw the number of lineages to migrate
    if (!migrLineages) return; // although migration happened no lineages will migrate
    
    // send lineages to each neighbor (multinomial...)
    double sum_m = 1.0, cur_m;
    unsigned int nbEmigrants;
    for (curFrom = curMigr.begin(), endFrom = curMigr.end(); migrLineages; ++curFrom) {
        cur_m = curFrom->second;   // get migration rate of the current neighbor
        nbEmigrants = get_pop_ptr()->rand().Binomial(cur_m / sum_m, migrLineages); // get number of lineages to emigrate to this neighbour
        sum_m -= cur_m;                     // adjust the sum of migration rates
        if (nbEmigrants) {
            curDeme->migrate(curFrom->first, nbEmigrants, *this, _demes, newDemes); // perform migration (maybe the deme has to be created)!!!
            migrLineages -= nbEmigrants;
        }
    }
}

// ----------------------------------------------------------------------------------------
/** perform migration
 * Caution "to" and "from" are in forward time sense!!!
 * it is possible that lineages migrate (backward in time) to a currently empty deme
 */
void LCE_Coalescence::migration_before_1D_SS(TDeme* curDeme,
                                             vector<pair<PATCH_ID, double> >& curMigr, const unsigned int& demeSize,
                                             map<PATCH_ID, TDeme*>& newDemes)
{
    assert(curDeme->get_lineages()); // get the number of lineages which may be send
    
    // get the total emigration rate
    double totMigr = 0;
    vector<pair<PATCH_ID, double> >::iterator curFrom, endFrom;
    for (curFrom = curMigr.begin(), endFrom = curMigr.end(); curFrom != endFrom; ++curFrom) {
        totMigr += curFrom->second; // the from deme may be not existent at the moment!
    }
    
    // get the number of lineages to migrate
    assert(demeSize);
    assert(totMigr>=0);
    assert(totMigr<=1);
    unsigned int migrLineages = get_pop_ptr()->rand().Binomial(totMigr, curDeme->get_lineages()); // draw the number of lineages to migrate
    if (!migrLineages) return; // although migration happened no lineages will migrate
    
    // send lineages to each neighbor (multinomial...)
    double sum_m = 1.0, cur_m;
    unsigned int nbEmigrants;
    for (curFrom = curMigr.begin(), endFrom = curMigr.end(); migrLineages; ++curFrom) {
        cur_m = curFrom->second;   // get migration rate of the current neighbor
        nbEmigrants = get_pop_ptr()->rand().Binomial(cur_m / sum_m, migrLineages); // get number of lineages to emigrate to this neighbour
        sum_m -= cur_m;                     // adjust the sum of migration rates
        if (nbEmigrants) {
            curDeme->migrate(curFrom->first, nbEmigrants, *this, _demes, newDemes); // perform migration (maybe the deme has to be created)!!!
            migrLineages -= nbEmigrants;
        }
    }
}

// ----------------------------------------------------------------------------------------
/** perform migration
 * Caution "to" and "from" are in forward time sense!!!
 * it is possible that lineages migrate (backward in time) to a currently empty deme
 */
void LCE_Coalescence::migration_before_2D_SS_4Neighbour(TDeme* curDeme,
                                                        vector<pair<PATCH_ID, double> >& curMigr, const unsigned int& demeSize,
                                                        map<PATCH_ID, TDeme*>& newDemes)
{
}

// ----------------------------------------------------------------------------------------
/** perform migration
 * Caution "to" and "from" are in forward time sense!!!
 * it is possible that lineages migrate (backward in time) to a currently empty deme
 */
void LCE_Coalescence::migration_before_2D_SS_8Neighbour(TDeme* curDeme,
                                                        vector<pair<PATCH_ID, double> >& curMigr, const unsigned int& demeSize,
                                                        map<PATCH_ID, TDeme*>& newDemes)
{
}

// ----------------------------------------------------------------------------------------
/** sets the dispersal rate container if one is specified, otherwise it is empty
 * caution this are immigration rates and not emigration rates!!!
 * dispRates[to]<from, migr_rate>
 */
void LCE_Coalescence::set_dispRates()
{
    // get parameter
    Param* p = _pCoalBase->get_parameter("coalescence_dispersal_rate");
    if (!p->isSet()) return;                 // no dispersal is set => stop here
    
    // rest params
    _dispRates.clear();
    _coal_migr_rate = 0;
    _coal_migr_rateIn = 0;
    _coal_migr_rateOut = 0;
    _coal_migr_rateCorner = 0;
    _coal_migr_rate_propagule = 0;
    
    if (p->is_matrix()) {                       // a dispersal matrix is passed
        func_ptr_migration_before = &LCE_Coalescence::migration_before_matrix;
        TMatrix* m = p->get_matrix();
        unsigned int i, j, nrows, ncols;
        double m_rate;
        nrows = m->getNbRows();
        ncols = m->getNbCols();
        if (nrows == ncols)
            error("Parameter 'coalescence_dispersal_rate' has a matrix with unique sizes!");
        if (nrows != get_pop_ptr()->get_nbPatch())
            error("Parameter 'coalescence_dispersal_rate' has a matrix with wrong dimensions!");
        vector<pair<PATCH_ID, double> >* curPatch;
        for (i = 0; i < ncols; ++i) {
            curPatch = &_dispRates[i + 1];
            for (j = 0; j < nrows; ++j) {
                if (i == j) continue;                  // remaining individuals
                m_rate = m->get(j, i);              // get migration rate
                if (!m_rate) continue;               // no migration
                curPatch->push_back(pair<PATCH_ID, double>(j + 1, m_rate));
            }
            if (curPatch->empty()) _dispRates.erase(i + 1); // no migration to this patch => remove it
        }
        delete m;
    }
    else {           // a single value is given
        _coal_migr_rate = p->get_value();
        switch (get_pop_ptr()->get_pDisperse_LCE()->get_disp_model()) {
            case 0: // island model
                _coal_migr_rate /= get_pop_ptr()->get_nbPatch() - 1;
                func_ptr_migration_before = &LCE_Coalescence::migration_before_island;
                break;
            case 1: { // island propagule model
                double propagule_prob = _pCoalBase->get_parameter_value("_coalescence_dispersal_propagule_prob");
                _coal_migr_rate_propagule = _coal_migr_rate * propagule_prob;
                _coal_migr_rate *= (1.0 - propagule_prob) / (get_pop_ptr()->get_nbPatch() - 2);
                func_ptr_migration_before = &LCE_Coalescence::migration_before_propagule;
                break;
            }
            case 2: { // 1D stepping stone model
                int factorIn, factorOut;
                switch (get_pop_ptr()->get_pDisperse_LCE()->get_border_model()) {	// edge effect
                    default:
                    case 0:
                        factorIn = 2;
                        factorOut = 2;
                        break;        	// circle
                    case 1:
                        factorIn = 1;
                        factorOut = 0;
                        break;					// reflecting
                    case 2:
                        factorIn = 2;
                        factorOut = -2;
                        break;			// absorbing (negative number means removing individuals)
                }
                
                _coal_migr_rateIn = _coal_migr_rate / factorIn;
                _coal_migr_rateOut = factorOut ? _coal_migr_rate / factorOut : 0;
                _coal_migr_rate /= 2;
                func_ptr_migration_before = &LCE_Coalescence::migration_before_1D_SS;
                break;
            }
            case 3: // 2D stepping stone model
                int factorIn4, factorOut4, factorInCorner4;
                int factorIn8, factorOut8, factorInCorner8;
                switch (get_pop_ptr()->get_pDisperse_LCE()->get_border_model()) {	// edge effect
                    default:
                    case 0: // torus
                        factorIn4 = 4;
                        factorIn8 = 8;
                        factorOut4 = 4;
                        factorOut8 = 8;
                        factorInCorner4 = 4;
                        factorInCorner8 = 8;
                        break;
                    case 1: // reflecting
                        factorIn4 = 3;
                        factorIn8 = 5;
                        factorOut4 = 0;
                        factorOut8 = 0;
                        factorInCorner4 = 2;
                        factorInCorner8 = 3;
                        break;
                    case 2: // absorbing
                        factorIn4 = 4;
                        factorIn8 = 8;
                        factorOut4 = -4;
                        factorOut8 = -8; // negative number means removing individuals
                        factorInCorner4 = 4;
                        factorInCorner8 = 8;
                        break;
                }
                
                // number of neighbors
                if (get_pop_ptr()->get_pDisperse_LCE()->get_lattice_range() == 0) { // 4 neighbors
                    func_ptr_migration_before = &LCE_Coalescence::migration_before_2D_SS_4Neighbour;
                    _coal_migr_rateIn = _coal_migr_rate / factorIn4;
                    _coal_migr_rateOut = factorOut4 ? _coal_migr_rate / factorOut4 : 0;
                    _coal_migr_rateCorner = _coal_migr_rate / factorInCorner4;
                    _coal_migr_rate /= 4;
                }
                else {                          // 8 neighbors
                    func_ptr_migration_before = &LCE_Coalescence::migration_before_2D_SS_8Neighbour;
                    _coal_migr_rateIn = _coal_migr_rate / factorIn8;
                    _coal_migr_rateOut = factorOut8 ? _coal_migr_rate / factorOut8 : 0;
                    _coal_migr_rateCorner = _coal_migr_rate / factorInCorner8;
                    _coal_migr_rate /= 8;
                }
                break;
        }
    }
}

// ----------------------------------------------------------------------------------------
void LCE_Coalescence::perform_coalescence()
{
    (this->*func_ptr_perform_coalescence)();
}


// ----------------------------------------------------------------------------------------
void LCE_Coalescence::perform_coalescence_singleEvent()
{
    map<PATCH_ID, TDeme*>::iterator curDeme, endDeme;
    for (curDeme = _demes.begin(), endDeme = _demes.end(); curDeme != endDeme; ++curDeme) {
#ifdef _DEBUG
        unsigned int nbLineages=_num_lineages;
#endif
        
        _num_lineages -= curDeme->second->singleCoalescentEvent(_curGen, *this);
        
#ifdef _DEBUG
        if(nbLineages!=_num_lineages) message("\n  Coalescence at generation: %u (%i lineages; %i demes)", _curGen, _num_lineages, _demes.size());
#endif
    }
}

// ----------------------------------------------------------------------------------------
void LCE_Coalescence::perform_coalescence_multipleEvent()
{
    map<PATCH_ID, TDeme*>::iterator curDeme, endDeme;
    for (curDeme = _demes.begin(), endDeme = _demes.end(); curDeme != endDeme; ++curDeme) {
#ifdef _DEBUG
        unsigned int nbLineages=_num_lineages;
#endif
        
        _num_lineages -= curDeme->second->multipleCoalescentEvents(_curGen, *this);
        
#ifdef _DEBUG
        if(nbLineages!=_num_lineages) message("\n  Coalescence at generation: %u (%i lineages; %i demes)", _curGen, _num_lineages, _demes.size());
#endif
    }
}

// ----------------------------------------------------------------------------------------
void LCE_Coalescence::perform_coalescence_mixedEvent()
{
    map<PATCH_ID, TDeme*>::iterator curDeme, endDeme;
    for (curDeme = _demes.begin(), endDeme = _demes.end(); curDeme != endDeme; ++curDeme) {
#ifdef _DEBUG
        unsigned int nbLineages=_num_lineages;
#endif
        
        _num_lineages -= curDeme->second->mixedCoalescentEvents(_curGen, *this);
        
#ifdef _DEBUG
        if(nbLineages!=_num_lineages) message("\n  Coalescence at generation: %u (%i lineages; %i demes)", _curGen, _num_lineages, _demes.size());
#endif
    }
}

// ----------------------------------------------------------------------------------------
/** adjust the deme size for the current generation
 * NEW: only the demes with lineages are updated!!! (more efficient)
 */
void LCE_Coalescence::perform_popSize_regulation(dbCoalescence* curDBs)
{
    map<PATCH_ID, TDeme*>::iterator curDeme, endDeme;
    map<PATCH_ID, dbPop>& db = curDBs->get_patches();
    map<PATCH_ID, dbPop>::iterator curDB, endDB = db.end();
    for (curDeme = _demes.begin(), endDeme = _demes.end(); curDeme != endDeme; ++curDeme) { // for each patch with lineages
        curDB = db.find(curDeme->first);
        if (curDB == endDB) {   // the patch is now empty
            if (!(this->*func_ptr_evacuate_lineages)(curDeme->second, true)) { // check if there are still lineages: evacuate them
                warning("Coalescence: Problem: not able to evacuate all lineages!!!\n");
                curDeme->second->get_lineages();
            }
        }
        curDeme->second->set_deme_size(curDB->second.get_popSize()); // set/change the new pop size
    }
}

// ------------------------------------------------------------------------------
// constructor
// ------------------------------------------------------------------------------
LCE_Coalescence::LCE_Coalescence(LCE_Coalescence_base* p)
{
    _dbCoalescence = NULL;
    _coal_migr_rate = 0;
    _coal_migr_rateIn = 0;
    _coal_migr_rateOut = 0;
    _coal_migr_rateCorner = 0;
    _coal_migr_rate_propagule = 0;
    _treeSize = my_NAN;
    _sampleNodes = NULL;
    _mrca = NULL;
    
    init(p);

}

// ------------------------------------------------------------------------------
/** destructor */
LCE_Coalescence::~LCE_Coalescence()
{
    clear_tree();
    
    map<PATCH_ID, TDeme*>::iterator cur, end=_demes.end();
    for(cur=_demes.begin(); cur!=end; ++cur){
        delete cur->second;
    }
}

// ------------------------------------------------------------------------------
/** recursively a deme is searched which can take the supernumerous lineages
 * false is returned if this is not possible, since all demes are full
 */
bool LCE_Coalescence::evacuate_stepping_stone(TDeme* origDeme, unsigned int id,
                                              unsigned int& lineages, unsigned int steps)
{
    if (!lineages) return true;  // stop if there are no more lineages to spread
    vector<unsigned int> vNeighbours = get_pop_ptr()->get_pDisperse_LCE()->get_neighbours(id); // get the neighboring demes of this deme
    unsigned int i, size = (unsigned int)vNeighbours.size();            // number of neighbors
    unsigned int* aNeighbours = new unsigned int[size];
    ARRAY::vector2array(vNeighbours, aNeighbours, size);
    get_pop_ptr()->rand().randomize(aNeighbours, size); // randomize the order of the neighbors
    
    // spread the lineages on the neighbors if they have free places
    map<PATCH_ID, TDeme*>::iterator cur, end = get_demes().end();
    TDeme* curNeighbour;
    for (i = 0; i < size && lineages; ++i) {
        cur = get_demes().find(aNeighbours[i]);
        if (cur == end) continue; // stop here if this neighbor is not present, i.e. not populated
        curNeighbour = cur->second;
        while (lineages && curNeighbour->get_lineages() < 2 * curNeighbour->get_deme_size()) {
            origDeme->migrate(curNeighbour->get_id(), 1, *this); // migrate a single lineage
            --lineages;
        }
    }
    
    // continue the spreading recursively from the neighboring demes (for each lineage separately)
    if (steps) {
        for (i = 0; lineages && i < size; ++i) {  // for each neighbor
            if (evacuate_stepping_stone(origDeme, aNeighbours[i], lineages, steps - 1)) {
                delete[] aNeighbours;
                return true;
            }
        }
    }
    delete[] aNeighbours;
    return (lineages == 0); // return true, if all lineages were spread and false if this was not possible
}

// ------------------------------------------------------------------------------
/** we have the tree: now generate the genetic data
 * seq[ind]
 */
void LCE_Coalescence::generateGeneticData(ALLELE** seq, TLocus** genome,
                                          const unsigned int& curLoc)
{
    // initialize the MRCA sequence
    (*genome)->ini_sequence(&_mrca->sequence, NULL, 1);
    
    // mutate recursively
    if (_mrca->desc1) _mrca->desc1->mutate(get_pop_ptr(), *_mrca, genome + curLoc);
    if (_mrca->desc2) _mrca->desc2->mutate(get_pop_ptr(), *_mrca, genome + curLoc);
    
    // get the genetic data
    for (unsigned int i = 0; i < _pCoalBase->get_nbSamples(); ++i) {
        seq[i][curLoc] = _sampleNodes[i]->sequence;
    }
}

// ----------------------------------------------------------------------------------------
// run_coalescence
// ----------------------------------------------------------------------------------------
void LCE_Coalescence::run_coalescence(unsigned int from, unsigned int to)
{
#ifdef _DEBUG
    message("\nLCE_Coalescence:");
#endif
    
    _curLocus = from;
    _endLocus = to;
    
    // coalescence simulations
    time_t startTime = time(0);
    int nbOutput = _pCoalBase->get_nbLocus() / 10;
    if (!nbOutput) nbOutput = 1;   // if less than 10 generations are simulated
    log2console(0, 0, _pCoalBase->get_nbLocus(), nbOutput);	// output to the console
    
    if (_pCoalBase->get_writer()) _pCoalBase->get_writer()->create_locusTree();
    for (; _curLocus < _endLocus; ++_curLocus) {
        // generate the coalescence tree
        if (_pCoalBase->get_writerLineages()) {
            _pCoalBase->get_writerLineages()->write_lineages_start(_curLocus);
            build_tree();                                         //////////////
            _pCoalBase->get_writerLineages()->write_lineages_end();
        }
        else build_tree_simple(); // faster since it does not check for temporal parameters or any other output
        if (_pCoalBase->get_writer()) _pCoalBase->get_writer()->write_NEXUS_tree(_curLocus, this);
        
        test_tree();
        
        // sprinkle mutations on the tree
        generateGeneticData(_seq, _pCoalBase->get_pop_ptr()->get_protoGenome()->get_locus_tot(), _curLocus); //////////////
        
        // clean up
        clear_tree();
        log2console(time(0) - startTime, _curLocus+1, _pCoalBase->get_nbLocus(), nbOutput);	// output to the console
    }
}

// ----------------------------------------------------------------------------------------
// set_up_metapop
// ----------------------------------------------------------------------------------------
unsigned long LCE_Coalescence::compute_treeSize()
{
    if(!_mrca) return 0;
    unsigned long size = 0;
    if (_mrca->desc1) _mrca->desc1->compute_treeSize(size);
    if (_mrca->desc2) _mrca->desc2->compute_treeSize(size);
    return size;
}

//----------------------------------------------------------------------------------------
// remove_mrca_from_deme
// ----------------------------------------------------------------------------------------
/** remove the MRCA node from the demes ChainList
 */
void LCE_Coalescence::remove_mrca_from_deme()
{
    assert(_mrca);
    if(_mrca->disp_vec.empty()){
        // remove the last lineage (the node is already in the tree)
        map<PATCH_ID, TDeme*>::iterator curPos, endPos;
        for (curPos = _demes.begin(), endPos = _demes.end(); curPos != endPos; ++curPos) {
            if (curPos->second->get_lineages()) break;
        }
        curPos->second->clear_chainNodeList();
        curPos->second->update_lineages();
    }
    else {
        assert(_mrca->disp_vec.size()==1);
        assert(get_demes(_mrca->disp_vec.front().second)->get_chainNodeList().size()==1);
        assert(get_demes(_mrca->disp_vec.front().second)->get_chainNodeList().front()==_mrca);
        TDeme* curDeme = get_demes(_mrca->disp_vec.front().second);
        curDeme->get_chainNodeList().clear();
        curDeme->update_lineages();
    }
}

// ----------------------------------------------------------------------------------------
// log2console
// ----------------------------------------------------------------------------------------
/** output to the console
 * l is the current locus (index starting at 1)
 * */
void LCE_Coalescence::log2console(double diffTime, unsigned int l, unsigned int nbLocus,
                                  unsigned int period)
{
    if (l <= 10 || !(l % period)) {
#ifdef _DEBUG
#ifdef _SHOW_MEMORY
        message("\n\r****   coalescence [%s] %i/%i   RAM: %f MB  ****",getElapsedTime(diffTime).c_str(),l,nbLocus, process_mem_usage());
#else
        message("\n\r****   coalescence [%s] %i/%i  ****",getElapsedTime(diffTime).c_str(),l,nbLocus);
#endif
#else
#ifdef _SHOW_MEMORY
        message("\r      coalescence [%s] %i/%i   RAM: %f MB",getElapsedTime(diffTime).c_str(),l,nbLocus, process_mem_usage());
#else
        message("\r      coalescence [%s] %i/%i", getElapsedTime(diffTime).c_str(), l, nbLocus);
#endif
#endif
        fflush(stdout);
    }
}

// ----------------------------------------------------------------------------------------
void LCE_Coalescence::clear_tree()
{
//	test_tree();
	if(!_mrca) return;
    
//    assert(_mrca->disp_vec.size()==1);
//    assert(get_demes(_mrca->disp_vec.front().second)->get_chainNodeList().size()==1);
//    assert(get_demes(_mrca->disp_vec.front().first)->get_chainNodeList().front()==_mrca);
//    //get_demes(_mrca->disp_vec.front().first)->get_chainNodeList.clear();
	clear_tree(_mrca);
	delete [] _sampleNodes;
    _sampleNodes = NULL;
    _mrca = NULL;
}

// ----------------------------------------------------------------------------------------
void LCE_Coalescence::clear_tree(TTNode* p)
{
	assert(p);
	if(p->desc1) clear_tree(p->desc1);
	if(p->desc2) clear_tree(p->desc2);
	delete p;

}

// ----------------------------------------------------------------------------------------
/** builds the three efficiently by just adjusting the deme sizes of the demes containing
 * lineages. Thus it is not possible to store the lineages_file!!!
 */
unsigned long LCE_Coalescence::build_tree()
{
    _pCoalBase->temporal_change(_pCoalBase->get_nbGen());
    if (_pCoalBase->get_writerLineages()) {
        Param* curParam = _pCoalBase->get_parameter("coalescence_lineages_logtime");
        curParam->set_arg_for_gen(_pCoalBase->get_nbGen());
        _pCoalBase->get_writerLineages()->set_GenerationOccurrence((unsigned int) curParam->get_value());
    }
    
    // Initialize the demes
    //set_ini_demeSizes();
    set_samples();
    
    _num_lineages = _pCoalBase->get_nbSamples();
    _timeDivergence = (unsigned long) _pCoalBase->get_parameter("divergence_time")->get_value(); // get the divergence time
    if(_timeDivergence<_pCoalBase->get_nbGen()) _timeDivergence = _pCoalBase->get_nbGen(); // divergence is like a migration event after coalescence event
    
    
#ifdef _DEBUG
    message("\n  Coalescence simulations: %i lineages:", _num_lineages);
#endif
    
    // build the tree over the time period of the demographic simulation
    dbCoalescence* curDB;
    for (_curGen = 0; _curGen < (unsigned long) _pCoalBase->get_nbGen(); ++_curGen) {// for each generation (starting with 1) in the db going backward in time
        curDB = &_dbCoalescence[_pCoalBase->get_nbGen() - _curGen - 1]; // starting from the back
        perform_popSize_regulation(curDB);     // Step 1: adapt population sizes
        if (_pCoalBase->get_writerLineages()) _pCoalBase->get_writerLineages()->write_current_lineages(_pCoalBase->get_nbGen(), _curGen, curDB, this);
        perform_coalescence();                      // Step 3: Coalescence round
        if(_num_lineages==1) break;
        perform_migration(curDB);                  	// Step 2: migration
        _pCoalBase->temporal_change(_pCoalBase->get_nbGen() - (unsigned int)_curGen);
    }
    
    // if the MRCA has not yet be found continue searching it
    if (_num_lineages > 1) {
#ifdef _DEBUG
        message("\n    At starting time (%u) there are %i lineages remaining.", _curGen, _num_lineages);
#endif
        
        // build the tree over the time period before the demographic simulation
        for (; _curGen < _timeDivergence; ++_curGen) {	// for each generation (starting with 1) in the db going backward in time
            if (_pCoalBase->get_writerLineages()) _pCoalBase->get_writerLineages()->write_current_lineages(_pCoalBase->get_nbGen(), _curGen, this);
            perform_coalescence();                  // Step 3: Coalescence round
            if(_num_lineages==1) break;
            perform_migration_before();                 // Step 2: migration
            _pCoalBase->temporal_change(_pCoalBase->get_nbGen() - (unsigned int)_curGen);
        }
        
        // if the MRCA has not yet be found continue searching it
        if (_num_lineages > 1) {
#ifdef _DEBUG
            message("\n    At divergence time (%u) there are %i lineages remaining.", _curGen, _num_lineages);
#endif
            
            // divergence is reached: merge all demes, if the new pop size is not specified the pop
            // size will be the total number of individual sof the metapop
            merge_demes();
            // continue until the last gene
            for (; _curGen < ULONG_MAX; ++_curGen) {// for each generation (starting with 1) in the db going backward in time
                if (_pCoalBase->get_writerLineages()) _pCoalBase->get_writerLineages()->write_current_lineages(_pCoalBase->get_nbGen(), _curGen, this);
                perform_coalescence();              // Step 3: Coalescence round
                if(_num_lineages==1) break;
                _pCoalBase->temporal_change(_pCoalBase->get_nbGen() - (unsigned int)_curGen);
            }
        }
    }
    
    if (_num_lineages > 1) error("Coalescence:: Coalescence process did not converge!\n");
    
#ifdef _DEBUG
    message("\n    MRCA reached at time (%u).", _curGen);
#endif
    
    _tmrca = _mrca->time;		// get the time to the MRCA
    _pCoalBase->set_vMRCA(_curLocus,_tmrca);
    
    
    // remove the last lineage (the node is already in the three)
    remove_mrca_from_deme();
    
    return _curGen;
}


// ----------------------------------------------------------------------------------------
// test_tree
// ----------------------------------------------------------------------------------------
/** test integrity of tree (using assert)
 */
void
LCE_Coalescence::test_tree(bool alsoSize)
{
    if(!_sampleNodes) return;
    
    TTNode *curNode, *ancestor;
    vector<pair<unsigned int, unsigned int> >::iterator cur, end;
    unsigned int lastTime, lastID;
    
#ifdef _DEBUG
    message("\n#################################\n");
#endif
    
    // for each sampled node
    for(unsigned int i=0; i<_pCoalBase->get_nbSamples(); ++i){
        curNode=_sampleNodes[i];
        
#ifdef _DEBUG
        message("\ncurrentNode: %i (%i; ", curNode->ID_Node, curNode->time);
        if(curNode->desc1) message("%i & ", curNode->desc1->ID_Node);
        else message("NULL & ");
        if(curNode->desc2) message("%i; ancestor: ", curNode->desc2->ID_Node);
        else message("NULL; ancestor: ");
        if(curNode->ancestor) message("%i)", curNode->ancestor->ID_Node);
        else message("NULL)");
#endif
        
        // test that it is a sample node
        assert(curNode->time==0);
        assert(curNode->desc1==NULL);
        assert(curNode->desc2==NULL);
        
        // loop through all nodes until the MRCA
        while(curNode->ancestor){
            
            ancestor=curNode->ancestor;
            assert(curNode->time <= ancestor->time);
            assert(ancestor->desc1 == curNode || ancestor->desc2 == curNode);
            assert(ancestor->desc1 != ancestor->desc2);
            //			assert(curNode->time == curNode->disp_vec.front().first);
            /*            assert((curNode->disp_vec.back().first < curNode->ancestor->time)     migration happenes after coaelscence
             || ((curNode->disp_vec.back().first == curNode->ancestor->time)
             && (curNode->disp_vec.size()==1                                except if the entry is the colescence event
             || curNode->ancestor->time==_timeDivergence)));           or the divergence time
             */
//            assert((curNode->disp_vec.back().first < curNode->ancestor->time)    // migration happenes after coaelscence
//                   || ((curNode->disp_vec.back().first == curNode->ancestor->time)
//                       && curNode->disp_vec.size()==1));                               // except if the entry is the colescence event
            
            // test the vector disp_vec
            lastTime = 0, lastID=my_NAN;
            for(cur=curNode->disp_vec.begin(), end=curNode->disp_vec.end(); cur!=end; ++cur){
                assert(cur->first>lastTime || lastTime==0 ||
                       (cur->first==lastTime && cur==++curNode->disp_vec.begin())); // migration is after coalescence
                assert(cur->second!=lastID);
                lastTime = cur->first;
                lastID = cur->second;
            }
            
            // the last migration has to lead to the coalescence deme
            //assert(curNode->disp_vec.back().second == curNode->ancestor->disp_vec.front().second);
            
            // move to the next one
            curNode = ancestor;
        }
        
        // it is the ancestor
        assert(_tmrca<=_curGen || _mrca == curNode);
    }
    
    unsigned long treeSize = compute_treeSize();
    if(alsoSize) assert(_treeSize == treeSize);
    else assert(_tmrca<=_curGen || _treeSize == treeSize);
    
}

// ----------------------------------------------------------------------------------------
// print_tree
// ----------------------------------------------------------------------------------------
/** print the tree to the console
 */
void
LCE_Coalescence::print_tree(TTNode* curNode)
{
    if(!curNode) return;
    
    print_tree(curNode->desc1);
    print_tree(curNode->desc2);
    
    cout << "\n\n********** Node: " << curNode->ID_Node + 1 <<
    " (Time: " << curNode->time <<
    "; Ancestor: ";
    if(curNode->ancestor) cout << curNode->ancestor->ID_Node + 1; else cout << "NULL";
    cout << "; Desc1: ";
    if(curNode->desc1) cout << curNode->desc1->ID_Node + 1; else cout << "NULL";
    cout << "; Desc2: ";
    if(curNode->desc2) cout << curNode->desc2->ID_Node + 1; else cout << "NULL";
    cout << ")";
    
    /*
     *
     for (unsigned int f = 0; f < curNode->disp_vec.size(); ++f) {
     cout << "\n  time: " << curNode->disp_vec[f].first << "\tdeme: "
					<< curNode->disp_vec[f].second;
     }
     */
    cout << flush;
    
}

