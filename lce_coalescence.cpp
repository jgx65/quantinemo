/** @file lce_coalescence.cpp
 *
 *   Copyright (C) 2010 Samuel Neuenschwander <samuel.neuenschwander@unil.ch>
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
Metapop* LCE_Coalescence::get_pop_ptr(){assert(_pCoalBase); return _pCoalBase->get_pop_ptr();}


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
    for (_curGen = 0; _curGen < (unsigned long) _pCoalBase->get_nbGen(); ++_curGen) { // for each generation (starting with 1) in the db going backward in time
        curDB = &_dbCoalescence[_pCoalBase->get_nbGen() - _curGen - 1]; // starting from the back
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
void LCE_Coalescence::generateGeneticData(unsigned char** seq, TLocus** genome,
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

// ----------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------
//// ----------------------------------------------------------------------------------------
///** adjust the deme size for the current generation
// * NEW: only the demes with lineages are updated!!! (more efficient)
// * this is for recombination for subsequent loci: also the recombining nodes are taken into account
// */
//void LCE_CoalescenceRecomb::perform_popSize_regulationIII(dbCoalescence* curDBs)
//{
//	map<PATCH_ID, TDeme*>::iterator curDeme, endDeme;
//	map<PATCH_ID, dbPop>& db = curDBs->get_patches();
//	map<PATCH_ID, dbPop>::iterator curDB, endDB = db.end();
//	for (curDeme = _demes.begin(), endDeme = _demes.end(); curDeme != endDeme; ++curDeme) { // for each patch with lineages
//		curDB = db.find(curDeme->first);
//		if (curDB == endDB) {   // the patch is now empty
//			if (!(this->*func_ptr_evacuate_lineages)(curDeme->second, true)) { // check if there are still lineages: evacuate them
//				warning("Coalescence: Problem: not able to evacuate all lineages!!!\n");
//				curDeme->second->get_lineages();
//			}
//		}
//		curDeme->second->set_deme_size(curDB->second.get_popSize()); // set/change the new pop size
//	}
//}
//
//// ----------------------------------------------------------------------------------------
///** perform migration
// * Caution "to" and "from" are in forward time sense!!!
// * Note, when emigrating the number of lineages is not adjusted, but the node added,
// * only afterwards the number of lineages is adapted to the immigrants. This to avoid
// * newly immigrated genes emigrate in the same generation
// */
//void LCE_CoalescenceRecomb::perform_migration_recomb(dbCoalescence* DB)
//{
//	// for each deme with lineages: migrate lineages
//	// if a lineage arrives in a new deme this deme is added to temp vector and
//	// only added after all migrations to _demes, since otherwise the iterator is invalidated
//	map<PATCH_ID, TDeme*> newDemes; //temp container to take the newly "colonized by lineages" demes
//	map<PATCH_ID, TDeme*>::iterator curDeme, endDeme;
////	map<PATCH_ID, dbPop>::iterator curDB, endDB = DB->get_patches().end();
//	for (curDeme = _demes.begin(), endDeme = _demes.end(); curDeme != endDeme; ++curDeme) {
//		assert(DB->get_patches().find(curDeme->first)!=DB->get_patches().end()); // the element has to exist
//		dbPop& dbDeme = DB->get_patches()[curDeme->first]; // get the current deme
//		if (dbDeme.get_immigrants().empty()) continue; // if no migration happened
//		migration(curDeme->second, dbDeme.get_immigrants(), curDeme->second->get_deme_size(),
//				newDemes);
//	}
//
//	if (!newDemes.empty()) _demes.insert(newDemes.begin(), newDemes.end()); // merge the containers
//
//	// for each deme with lineages: update the number of lineages
//	// if a deme does not anymore contain lineages it is removed from _demes
//	// endDeme is not anymore the same as above, and may vary during looping!!!
//	for (curDeme = _demes.begin(); curDeme != _demes.end();) {
//		if (curDeme->second->post_migration_lineages()) ++curDeme; // update the number of lineages
//		else{
//			delete curDeme->second;
//			_demes.erase(curDeme++);
//		}
//	}
//}
//
//// ----------------------------------------------------------------------------------------
///** perform migration
// * Caution "to" and "from" are in forward time sense!!!
// * Note, when emigrating the number of lineages is not adjusted, but the node added,
// * only afterwards the number of lineages is adapted to the immigrants. This to avoid
// * that newly immigrated genes emigrate in the same generation
// * only recombining lineages are considered!!!
// */
//void LCE_CoalescenceRecomb::perform_migration_recomb(dbCoalescence* DB, map<PATCH_ID, vector<TTNode*> >& recombDemes)
//{
//	// for each deme with lineages: migrate lineages
//	// if a lineage arrives in a new deme this deme is added to temp vector and
//	// only added after all migrations to _demes, since otherwise the iterator is invalidated
//	map<PATCH_ID, TDeme*> newDemes; //temp container to take the newly "colonized by lineages" demes
//	map<PATCH_ID, vector<TTNode*> > newRecombDemes; //temp container to take the newly "colonized by lineages" demes
//	map<PATCH_ID, vector<TTNode*> >::iterator curDemeID, endDemeID;
//	for (curDemeID = recombDemes.begin(), endDemeID = recombDemes.end(); curDemeID != endDemeID; ++curDemeID) {
//		assert(DB->get_patches().find(curDemeID->first)!=DB->get_patches().end()); // the element has to exist
//		dbPop& dbDeme = DB->get_patches()[curDemeID->first]; // get the current deme
//		if (dbDeme.get_immigrants().empty()) continue; // if no migration happened
//		assert(_demes.find(curDemeID->first)!=_demes.end());
//		migrationRecomb(_demes[curDemeID->first], dbDeme.get_immigrants(),
//				_demes[curDemeID->first]->get_deme_size(),
//				newDemes, recombDemes, newRecombDemes);
//	}
//
//	if (!newDemes.empty()){
//		_demes.insert(newDemes.begin(), newDemes.end()); // merge the containers
//	}
//    
//    if (!newRecombDemes.empty()){
//        recombDemes.insert(newRecombDemes.begin(), newRecombDemes.end()); // merge the containers
//    }
//
//	// for each deme with lineages: update the number of lineages
//	// empty demes are NOT deleted as they may contain new recombining lineages
//	map<PATCH_ID, TDeme*>::iterator curDeme, endDeme;
//	for (curDeme = _demes.begin(); curDeme != _demes.end(); ++curDeme) {
//		curDeme->second->post_migration_lineages();
//	}
//}
//
//// ----------------------------------------------------------------------------------------
///** perform migration
// * Caution "to" and "from" are in forward time sense!!!
// * Note, when emigrating the number of lineages is not adjusted, but the node added,
// * only afterwards the number of lineages is adapted to the immigrants. This to avoid
// * that newly immigrated genes emigrate in the same generation
// * only recombining lineages are considered!!!
// */
//void LCE_CoalescenceRecomb::perform_migration_recomb(dbCoalescence* DB, PATCH_ID& recombDemeID, TTNode*& recombNode)
//{
//	// for each deme with lineages: migrate lineages
//	// if a lineage arrives in a new deme this deme is added to temp vector and
//	// only added after all migrations to _demes, since otherwise the iterator is invalidated
//    assert(DB->get_patches().find(recombDemeID)!=DB->get_patches().end()); // the element has to exist
//    dbPop& dbDeme = DB->get_patches()[recombDemeID]; // get the current deme
//    if (dbDeme.get_immigrants().empty()) return; // if no migration happened
//    
//	map<PATCH_ID, TDeme*> newDemes; //temp container to take the newly "colonized by lineages" demes
//    assert(_demes.find(recombDemeID)!=_demes.end());
//    TDeme* pDeme = _demes[recombDemeID];
//    migrationRecomb(pDeme, dbDeme.get_immigrants(), pDeme->get_deme_size(),
//                    newDemes, recombDemeID, recombNode);
//    
//	if (!newDemes.empty()){
//		_demes.insert(newDemes.begin(), newDemes.end()); // merge the containers
//	}
//    
//	// for each deme with lineages: update the number of lineages
//	// empty demes are NOT deleted as they may contain new recombining lineages
//	map<PATCH_ID, TDeme*>::iterator curDeme, endDeme;
//	for (curDeme = _demes.begin(); curDeme != _demes.end(); ++curDeme) {
//		curDeme->second->post_migration_lineages();
//	}
//}
//
//// ----------------------------------------------------------------------------------------
///** perform migration before the onset of demographic simulation until divergence
// * Caution "to" and "from" are in forward time sense!!!
// */
//void LCE_CoalescenceRecomb::perform_migration_before_recomb()
//{
//	// for each deme with lineages: migrate lineages
//	// if a lineage arrives in a new deme this deme is added to temp vector and
//	// only added after all migrations to _demes, since otherwise the iterator is invalidated
//	map<PATCH_ID, TDeme*> newDemes; //temp container to take the newly "colonized by lineages" demes
//	map<PATCH_ID, TDeme*>::iterator curDeme, endDeme;
//	map<PATCH_ID, vector<pair<PATCH_ID, double> > >::iterator pos;
//	for (curDeme = _demes.begin(), endDeme = _demes.end(); curDeme != endDeme; ++curDeme) {
//		pos = _dispRates.find(curDeme->first);           // get the current deme
//		if (pos == _dispRates.end()) continue;       // if no migration happened
//		(this->*func_ptr_migration_before)(curDeme->second, pos->second,
//				curDeme->second->get_deme_size(), newDemes);
//	}
//
//	if (!newDemes.empty()) _demes.insert(newDemes.begin(), newDemes.end()); // merge the containers
//
//	// for each deme with lineages: update the number of lineages
//	// if a deme does not anymore contain lineages it is removed from _demes
//	// endDeme is not anymore the same as above, and may vary during looping!!!
//	for (curDeme = _demes.begin(); curDeme != _demes.end();) {
//		if (curDeme->second->post_migration_lineages()) ++curDeme; // update the number of lineages
//		else _demes.erase(curDeme++);
//	}
//}
//
//// ----------------------------------------------------------------------------------------
///** perform migration
// * Caution "to" and "from" are in forward time sense!!!
// * it is possible that lineages migrate (backward in time) to a currently empty deme
// * only recombining genes are migrated!!!
// */
//void LCE_CoalescenceRecomb::migrationRecomb(TDeme* curDeme,
//		vector<pair<PATCH_ID, POP_SIZE> >& curMigr,
//		const unsigned int& demeSize,
//		map<PATCH_ID, TDeme*>& newDemes,
//		map<PATCH_ID, vector<TTNode*> >& recombDemes,
//		map<PATCH_ID, vector<TTNode*> >& newRecombDemes)
//{
//	assert(recombDemes.find(curDeme->get_id())!=recombDemes.end());
//	vector<TTNode*>& curVecNodes = recombDemes[curDeme->get_id()];
//	assert(!curVecNodes.empty());
//
//	// get the total number of immigrants
//	unsigned int totMigr = 0;
//	vector<pair<PATCH_ID, POP_SIZE> >::iterator curFrom, endFrom;
//	for (curFrom = curMigr.begin(), endFrom = curMigr.end(); curFrom != endFrom; ++curFrom) {
//		totMigr += (unsigned int) curFrom->second; // the from deme may be not existent at the moment!
//	}
//
//	// get the number of lineages to migrate
//	assert(demeSize);
//	assert(((double)totMigr)/demeSize>=0);
//	assert(((double)totMigr)/demeSize<=1);
//	unsigned int migrLineages = get_pop_ptr()->rand().Binomial(((double) totMigr) / demeSize,
//			curVecNodes.size()); // draw the number of lineages to migrate
//	if (!migrLineages) return; // although migration happened no lineages will migrate
//
//	// send lineages to each neighbor (multinomial...)
//	// store them in a separate container and merge them only after a round
//	double sum_m = 1.0, cur_m;
//	unsigned int randNode, nbEmigrants;
//	for (curFrom = curMigr.begin(), endFrom = curMigr.end(); migrLineages; ++curFrom) {
//		cur_m = (double) curFrom->second / totMigr; // get migration rate of the current neighbor
//		nbEmigrants = get_pop_ptr()->rand().Binomial(cur_m / sum_m, migrLineages); // get number of lineages to emigrate to this neighbour
//		sum_m -= cur_m;                     // adjust the sum of migration rates
//		if (nbEmigrants) {
//			migrLineages -= nbEmigrants;
//
//			// check if the deme already exists
//			vector<TTNode*>*  sinkVecNodes;
//			if(_demes.find(curFrom->first)==_demes.end()){
//				newDemes[curFrom->first] = new TDeme(get_pop_ptr(), curFrom->first);
//				sinkVecNodes = &newRecombDemes[curFrom->first];
//			}
//			else sinkVecNodes = &recombDemes[curFrom->first];
//
//			for(;nbEmigrants; --nbEmigrants){
//				randNode = get_pop_ptr()->rand().Uniform((unsigned int)curVecNodes.size());
//				TTNode* curNode = curVecNodes[randNode];				// draw randomly a node
//				sinkVecNodes->push_back(curNode);        				// add the lineage to the new deme
//				curVecNodes.erase(curVecNodes.begin()+randNode);       	// remove the lineage from the current deme
//			}
//		}
//	}
//}
//
//// ----------------------------------------------------------------------------------------
///** perform migration
// * Caution "to" and "from" are in forward time sense!!!
// * it is possible that lineages migrate (backward in time) to a currently empty deme
// * only recombining genes are migrated!!!
// * update disp_vec_temp
// */
//void LCE_CoalescenceRecomb::migrationRecomb(TDeme* curDeme,
//                                            vector<pair<PATCH_ID, POP_SIZE> >& curMigr,
//                                            const unsigned int& demeSize,
//                                            map<PATCH_ID, TDeme*>& newDemes,
//                                            PATCH_ID& recombDemeID,
//                                            TTNode*& recombNode)
//{
//	assert(recombNode);
//    
//	// get the total number of immigrants
//	unsigned int totMigr = 0;
//	vector<pair<PATCH_ID, POP_SIZE> >::iterator curFrom, endFrom;
//	for (curFrom = curMigr.begin(), endFrom = curMigr.end(); curFrom != endFrom; ++curFrom) {
//		totMigr += (unsigned int) curFrom->second; // the from deme may be not existent at the moment!
//	}
//    
//	// get the number of lineages to migrate
//	assert(demeSize);
//	assert(((double)totMigr)/demeSize>=0);
//	assert(((double)totMigr)/demeSize<=1);
//	unsigned int migrLineages = get_pop_ptr()->rand().Binomial(((double) totMigr) / demeSize, 1); // draw the number of lineages to migrate
//	if (!migrLineages) return; // although migration happened no lineages will migrate
//    
//    assert(migrLineages==1);
//    
//	// send lineages to each neighbor (multinomial...)
//	// store them in a separate container and merge them only after a round
//	double sum_m = 1.0, cur_m;
//	unsigned int nbEmigrants;
//	for (curFrom = curMigr.begin(), endFrom = curMigr.end(); migrLineages; ++curFrom) {
//		cur_m = (double) curFrom->second / totMigr; // get migration rate of the current neighbor
//		nbEmigrants = get_pop_ptr()->rand().Binomial(cur_m / sum_m, migrLineages); // get number of lineages to emigrate to this neighbour
//		sum_m -= cur_m;                     // adjust the sum of migration rates
//		if (nbEmigrants) {
//			// check if the deme already exists
//			if(_demes.find(curFrom->first)==_demes.end()){
//				newDemes[curFrom->first] = new TDeme(get_pop_ptr(), curFrom->first);   // add the newley colonized deme
//			}
//#ifdef _DEBUG
//            cout << "\n  **** recomb migration node " << recombNode->ID_Node+1 << "' (deme "
//            << curDeme->get_id()+1 << " => " << curFrom->first+1 << ")" << flush;
//#endif
//			recombDemeID = curFrom->first;  // update the recombDemeID to the new location
//            recombNode->disp_vec_temp.push_back(pair<unsigned int, unsigned int>(_curGen, curFrom->first));
//            
//            break; // we have just a single node to migrate
//		}
//	}
//    assert(recombDemeID==recombNode->get_currentDeme_temp(_curGen+1));
//}
//
//// ----------------------------------------------------------------------------------------
//void LCE_CoalescenceRecomb::perform_coalescenceII()
//{
//	(this->*func_ptr_perform_coalescenceII)();
//}
//
//// ----------------------------------------------------------------------------------------
//void LCE_CoalescenceRecomb::perform_coalescence_singleEventII()
//{
//	map<PATCH_ID, TDeme*>::iterator curDeme, endDeme;
//	for (curDeme = _demes.begin(), endDeme = _demes.end(); curDeme != endDeme; ++curDeme) {
//#ifdef _DEBUG
//		unsigned int nbLineages=_num_lineages;
//#endif
//
//		_num_lineages -= curDeme->second->singleCoalescentEventII(_curGen, *this);
//
//#ifdef _DEBUG
//		if(nbLineages!=_num_lineages) message("\n  Coalescence at generation: %u (%i lineages; %i demes)", _curGen, _num_lineages, _demes.size());
//#endif
//	}
//}
//
//// ----------------------------------------------------------------------------------------
//void LCE_CoalescenceRecomb::perform_coalescence_multipleEventII()
//{
//	map<PATCH_ID, TDeme*>::iterator curDeme, endDeme;
//	for (curDeme = _demes.begin(), endDeme = _demes.end(); curDeme != endDeme; ++curDeme) {
//#ifdef _DEBUG
//		unsigned int nbLineages=_num_lineages;
//#endif
//
//		_num_lineages -= curDeme->second->multipleCoalescentEventsII(_curGen, *this);
//
//#ifdef _DEBUG
//		if(nbLineages!=_num_lineages) message("\n  Coalescence at generation: %u (%i lineages; %i demes)", _curGen, _num_lineages, _demes.size());
//#endif
//	}
//}
//
//// ----------------------------------------------------------------------------------------
//void LCE_CoalescenceRecomb::perform_coalescence_mixedEventII()
//{
//	map<PATCH_ID, TDeme*>::iterator curDeme, endDeme;
//	for (curDeme = _demes.begin(), endDeme = _demes.end(); curDeme != endDeme; ++curDeme) {
//#ifdef _DEBUG
//		unsigned int nbLineages=_num_lineages;
//#endif
//
//		_num_lineages -= curDeme->second->mixedCoalescentEventsII(_curGen, *this);
//
//#ifdef _DEBUG
//		if(nbLineages!=_num_lineages) message("\n  Coalescence at generation: %u (%i lineages; %i demes)", _curGen, _num_lineages, _demes.size());
//#endif
//	}
//}
//
//// ----------------------------------------------------------------------------------------
//// ----------------------------------------------------------------------------------------
//void LCE_CoalescenceRecomb::perform_coalescenceIII(map<PATCH_ID, vector<TTNode*> >& recombNodes)
//{
//	(this->*func_ptr_perform_coalescenceIII)(recombNodes);
//}
//
//// ----------------------------------------------------------------------------------------
///** iterate only over the active demes (with active nodes)
// */
//void LCE_CoalescenceRecomb::perform_coalescence_singleEventIII(map<PATCH_ID, vector<TTNode*> >& activeNodes)
//{
//	map<PATCH_ID, vector<TTNode*> >::iterator curDeme, endDeme;
//	for (curDeme = activeNodes.begin(), endDeme = activeNodes.end(); curDeme != endDeme;
//			++curDeme) {
//#ifdef _DEBUG
//		unsigned int nbLineages=_num_lineages;
//#endif
//
//		_num_lineages -= _demes[curDeme->first]->singleCoalescentEventIII(_curGen, *this, curDeme->second);
//
//#ifdef _DEBUG
//		if(nbLineages!=_num_lineages) message("\n  Coalescence at generation: %u (%i lineages; %i demes)", _curGen, _num_lineages, _demes.size());
//#endif
//	}
//}
//
//// ----------------------------------------------------------------------------------------
///** iterate only over the active demes (with active nodes)
// */
//void LCE_CoalescenceRecomb::perform_coalescence_multipleEventIII(map<PATCH_ID, vector<TTNode*> >& activeNodes)
//{
//	map<PATCH_ID, vector<TTNode*> >::iterator curDeme, endDeme;
//	for (curDeme = activeNodes.begin(), endDeme = activeNodes.end(); curDeme != endDeme;
//         ++curDeme) {
//#ifdef _DEBUG
//		unsigned int nbLineages=_num_lineages;
//#endif
//        
//		_num_lineages -= _demes[curDeme->first]->multipleCoalescentEventsIII(_curGen, *this, curDeme->second);
//        
//#ifdef _DEBUG
//		if(nbLineages!=_num_lineages) message("\n  Coalescence at generation: %u (%i lineages; %i demes)", _curGen, _num_lineages, _demes.size());
//#endif
//	}
//}
//
//// ----------------------------------------------------------------------------------------
///** iterate only over the active demes (with active nodes)
// */
//void LCE_CoalescenceRecomb::perform_coalescence_mixedEventIII(map<PATCH_ID, vector<TTNode*> >& recombNodes)
//{
//	map<PATCH_ID, vector<TTNode*> >::iterator curDeme, endDeme;
//	for (curDeme = recombNodes.begin(), endDeme = recombNodes.end(); curDeme != endDeme;
//         ++curDeme) {
//#ifdef _DEBUG
//		unsigned int nbLineages=_num_lineages;
//#endif
//		assert(_demes.find(curDeme->first) != _demes.end());
//		_num_lineages -= _demes[curDeme->first]->mixedCoalescentEventsIII(_curGen, *this, curDeme->second);
//        
//#ifdef _DEBUG
//		if(nbLineages!=_num_lineages) message("\n  Coalescence at generation: %u (%i lineages; %i demes)", _curGen, _num_lineages, _demes.size());
//#endif
//	}
//}
//
//// ----------------------------------------------------------------------------------------
//// ----------------------------------------------------------------------------------------
//void LCE_CoalescenceRecomb::perform_coalescenceIV(TDeme* recombDeme, TTNode*& recombNode)
//{
//	(this->*func_ptr_perform_coalescenceIV)(recombDeme, recombNode);
//}
//
//// ----------------------------------------------------------------------------------------
///** iterate only over the active demes (with active sinlge node)
// */
//void LCE_CoalescenceRecomb::perform_coalescence_singleEventIV(TDeme* recombDeme, TTNode*& recombNode)
//{
//#ifdef _DEBUG
//    unsigned int nbLineages=_num_lineages;
//#endif
//    
//    _num_lineages -= recombDeme->singleCoalescentEventIV(_curGen, *this, recombNode);
//    
//#ifdef _DEBUG
//    if(nbLineages!=_num_lineages) message("\n  Coalescence at generation: %u (%i lineages; %i demes)", _curGen, _num_lineages, _demes.size());
//#endif
//}
//
//// ----------------------------------------------------------------------------------------
///** iterate only over the active demes (with active single node)
// */
//void LCE_CoalescenceRecomb::perform_coalescence_multipleEventIV(TDeme* recombDeme, TTNode*& recombNode)
//{
//#ifdef _DEBUG
//    unsigned int nbLineages=_num_lineages;
//#endif
//    
//    _num_lineages -= recombDeme->multipleCoalescentEventsIV(_curGen, *this, recombNode);
//    
//#ifdef _DEBUG
//    if(nbLineages!=_num_lineages) message("\n  Coalescence at generation: %u (%i lineages; %i demes)", _curGen, _num_lineages, _demes.size());
//#endif
//}
//
//// ----------------------------------------------------------------------------------------
///** iterate only over the active demes (with active single node)
// */
//void LCE_CoalescenceRecomb::perform_coalescence_mixedEventIV(TDeme* recombDeme, TTNode*& recombNode)
//{
//#ifdef _DEBUG
//    unsigned int nbLineages=_num_lineages;
//#endif
//    _num_lineages -= recombDeme->mixedCoalescentEventsIV(_curGen, *this, recombNode);
//    
//#ifdef _DEBUG
//    if(nbLineages!=_num_lineages) message("\n  Coalescence at generation: %u (%i lineages; %i demes)", _curGen, _num_lineages, _demes.size());
//#endif
//}
//
//// ----------------------------------------------------------------------------------------
//void LCE_CoalescenceRecomb::init(LCE_Coalescence_base* pBase)
//{
//    _pCoalBase = pBase;
//    _model_threshold = _pCoalBase->get_model_threshold();
//    _dbCoalescence = _pCoalBase->get_dbCoalescence();
//    
//    switch (_popPtr->get_pDisperse_LCE()->get_disp_model()) {
//        case 0: // island model
//        case 1: // island propagule model
//            func_ptr_evacuate_lineages = &LCE_Coalescence::evacuate_lineages_island;
//            break;
//        case 2: // 1D stepping stone model
//        case 3: // 2D stepping stone model
//            func_ptr_evacuate_lineages = &LCE_Coalescence::evacuate_lineages_stepping_stone;
//            break;
//    }
//    
//    set_dispRates();
//    
//	if (_model_threshold == 0) {
//		func_ptr_perform_coalescence = &LCE_Coalescence::perform_coalescence_multipleEvent;
//		func_ptr_perform_coalescenceII = &LCE_CoalescenceRecomb::perform_coalescence_multipleEventII;
//		func_ptr_perform_coalescenceIII = &LCE_CoalescenceRecomb::perform_coalescence_multipleEventIII;
//	}
//	else if (_model_threshold == 1e6) {
//		func_ptr_perform_coalescence = &LCE_Coalescence::perform_coalescence_singleEvent;
//		func_ptr_perform_coalescenceII = &LCE_CoalescenceRecomb::perform_coalescence_singleEventII;
//		func_ptr_perform_coalescenceIII = &LCE_CoalescenceRecomb::perform_coalescence_singleEventIII;
//	}
//	else {
//		func_ptr_perform_coalescence = &LCE_Coalescence::perform_coalescence_mixedEvent;
//		func_ptr_perform_coalescenceII = &LCE_CoalescenceRecomb::perform_coalescence_mixedEventII;
//		func_ptr_perform_coalescenceIII = &LCE_CoalescenceRecomb::perform_coalescence_mixedEventIII;
//	}
//}
//
//// ----------------------------------------------------------------------------------------
///** builds the three efficiently by just adjusting the deme sizes of the demes containing
// * lineages. Thus it is not possible to store the lineages_file!!!
// */
//unsigned long LCE_CoalescenceRecomb::build_tree()
//{
//	temporal_change(_nbGen);
//	if (_writerLineages) {
//		Param* curParam = _paramSet->get_param("coalescence_lineages_logtime");
//		curParam->set_arg_for_gen(_nbGen);
//		_writerLineages->set_GenerationOccurrence((unsigned int) curParam->get_value());
//	}
//
//	// Initialize the demes
//	//set_ini_demeSizes();
//	set_samples();
//
//	_num_lineages = _pCoalBase->get_nbSamples();
//    _timeDivergence = (unsigned long) get_parameter("divergence_time")->get_value(); // get the divergence time
//    if(_timeDivergence<_nbGen) _timeDivergence = _nbGen; // divergence is like a migration event after coalescence event
//
//
//#ifdef _DEBUG
//	message("\n  Coalescence simulations: %i lineages:", _num_lineages);
//#endif
//
//	// build the tree over the time period of the demographic simulation
//	dbCoalescence* curDB;
//	for (_curGen = 0; _curGen < (unsigned long) _nbGen; ++_curGen) {// for each generation (starting with 1) in the db going backward in time
//		curDB = &_dbCoalescence[_nbGen - _curGen - 1]; // starting from the back
//		perform_popSize_regulation(curDB);     // Step 1: adapt population sizes
//		if (_writerLineages) _writerLineages->write_current_lineagesII(_nbGen, _curGen, curDB);
//		perform_coalescenceII();                      // Step 3: Coalescence round
//        if(_num_lineages==1) break;
//		perform_migration(curDB);                  	// Step 2: migration
//		temporal_change(_nbGen - (unsigned int)_curGen);
//	}
//
//	// if the MRCA has not yet be found continue searching it
//	if (_num_lineages > 1) {
//#ifdef _DEBUG
//		message("\n    At starting time (%u) there are %i lineages remaining.", _curGen, _num_lineages);
//#endif
//
//		// build the tree over the time period before the demographic simulation
//		for (; _curGen < _timeDivergence; ++_curGen) {	// for each generation (starting with 1) in the db going backward in time
//			if (_writerLineages) _writerLineages->write_current_lineagesII(_nbGen, _curGen);
//			perform_coalescenceII();                  // Step 3: Coalescence round
//            if(_num_lineages==1) break;
//			perform_migration_before();                 // Step 2: migration
//			temporal_change(_nbGen - (unsigned int)_curGen);
//		}
//
//		// if the MRCA has not yet be found continue searching it
//		if (_num_lineages > 1) {
//#ifdef _DEBUG
//			message("\n    At divergence time (%u) there are %i lineages remaining.", _curGen, _num_lineages);
//#endif
//
//			// divergence is reached: merge all demes, if the new pop size is not specified the pop
//			// size will be the total number of individual sof the metapop
//			merge_demes();
//			// continue until the last gene
//			for (; _curGen < ULONG_MAX; ++_curGen) {// for each generation (starting with 1) in the db going backward in time
//				if (_writerLineages) _writerLineages->write_current_lineagesII(_nbGen, _curGen);
//				perform_coalescenceII();              // Step 3: Coalescence round
//                if(_num_lineages==1) break;
//				temporal_change(_nbGen - (unsigned int)_curGen);
//			}
//		}
//	}
//
//	if (_num_lineages > 1) error("Coalescence:: Coalescence process did not converge!\n");
//
//#ifdef _DEBUG
//	message("\n    MRCA reached at time (%u).", _curGen);
//#endif
//
//	_tmrca = _mrca->time;		// get the time to the MRCA
//	_vMRCA.push_back(_tmrca);
//
//	// remove the last lineage (the node is already in the three)
//	remove_mrca_from_deme();
//
//	return _curGen;
//}
//
//// ----------------------------------------------------------------------------------------
///** builds the three efficiently by just adjusting the deme sizes of the demes containing
// * lineages. Thus it is not possible to store the lineages_file!!!
// * adds the patchID to each node in order to retake the position
// */
//unsigned long LCE_CoalescenceRecomb::build_tree_simple()
//{
//	// Initialize the demes
//	set_samples();
//
//	_num_lineages = _pCoalBase->get_nbSamples();
//	_curGen = 0;
//    _timeDivergence = (unsigned long) get_parameter("divergence_time")->get_value(); // get the divergence time
//    if(_timeDivergence<_nbGen) _timeDivergence = _nbGen; // divergence is like a migration event after coalescence event
//
//
//#ifdef _DEBUG
//	message("\n  Coalescence simulations:");
//	message("\n    Generation: %u (%i lineages; %i demes)", _curGen, _num_lineages, _demes.size());
//#endif
//
//	// build the tree over the time period of the demographic simulation
//	dbCoalescence* curDB;
//	for (; _curGen < (unsigned long) _nbGen; ++_curGen) { // for each generation (starting with 1) in the db going backward in time
//		cout << endl << _curGen << ": " << _demes.size() << " demes; " << _num_lineages << " lineages" << endl;
//		curDB = &_dbCoalescence[_nbGen - _curGen - 1]; // starting from the back
//		perform_popSize_regulation(curDB);     // Step 1: adapt population sizes
//		perform_coalescenceII();                    // Step 3: Coalescence round
//        if(_num_lineages==1) break;
//		perform_migration(curDB);                   // Step 2: migration
//	}
//
//	// if the MRCA has not yet been found continue searching it
//	if (_num_lineages > 1) {
//#ifdef _DEBUG
//		message("\n    At starting time (%u) there are %i lineages remaining.", _curGen, _num_lineages);
//#endif
//
//		// build the tree over the time period before the demographic simulation
//		for (; _curGen < _timeDivergence; ++_curGen) {	// for each generation (starting with 1) in the db going backward in time
//			perform_coalescenceII();                // Step 3: Coalescence round
//            if(_num_lineages==1) break;
//			perform_migration_before();                 // Step 2: migration
//			cout << endl << _curGen << ": " <<  _demes.size() << " demes; " << _num_lineages << " lineages" << endl;
//		}
//
//		// if the MRCA has not yet been found continue searching it
//		if (_num_lineages > 1) {
//#ifdef _DEBUG
//			message("\n    At divergence time (%u) there are %i lineages remaining.", _curGen, _num_lineages);
//#endif
//			// divergence is reached: merge all demes, if the new pop size is not specified the pop
//			// size will be the total number of individual sof the metapop
//            merge_demes(); // merge demes is a migration event and is at the end of the generation
//			// continue until the last gene
//			for (; _curGen < ULONG_MAX && _num_lineages > 1; ++_curGen) {// for each generation (starting with 1) in the db going backward in time
//				perform_coalescenceII();            // Step 3: Coalescence round
//			}
//		}
//	}
//
//	if (_num_lineages > 1) error("Coalescence:: Coalescence process did not converge!\n");
//
//#ifdef _DEBUG
//	message("\n    MRCA reached at time (%u).", _curGen);
//#endif
//
//	_tmrca = _mrca->time;		// get the time to the MRCA
//	_vMRCA.push_back(_tmrca);
//
//
//	// remove the last lineage (the node is already in the three)
//	remove_mrca_from_deme();
//
//	return _curGen;
//}
//
//// ------------------------------------------------------------------------------
///** merge all remaining lineages in a single (first) deme
// * the new pop size may be specified by the parameter "divergence_pop_size"
// * if not specified the new pop size is the sum of all individuals of the metapop at time 0
// * if the new pop size is smaller than the number of lineages, the pop size is readjusted
// */
//void LCE_CoalescenceRecomb::merge_demes()
//{
//    assert(_curGen==_timeDivergence);
//
//	// merge the lineages (nodes) of all demes
//	assert(!_demes.empty());
//	map<PATCH_ID, TDeme*>::iterator curDeme, endDeme;
//	curDeme = _demes.begin();
//	endDeme = _demes.end();
//	TDeme* firstDeme = curDeme->second;
//	unsigned int firstDemeID = firstDeme->get_id();
//	vector<TTNode*>::iterator curNode, endNode;
//	for (++curDeme; curDeme != endDeme; ++curDeme) { // we start at the second deme...
//		//assert(curDeme->second->get_lineages() || curDeme->first==recombDemeID);
//        
//		// enter the migration
//        curNode = curDeme->second->get_chainNodeList().begin();
//		endNode = curDeme->second->get_chainNodeList().end();
//		for (; curNode != endNode; ++curNode) {
//#ifdef _DEBUG
//            cout << "\n  **** merge demes: node " << (*curNode)->ID_Node+1
//            << " (deme " << curDeme->second->get_id()+1 << " => " << firstDemeID+1 << ")" << flush;
//#endif
//
//            // check that the node not just emigrated to it (exception if it is the coalescence entry)
//            if((*curNode)->disp_vec.back().first==_timeDivergence-1 && (*curNode)->disp_vec.size()!=1){
//                (*curNode)->disp_vec.pop_back();
//                if((*curNode)->disp_vec.back().second!=firstDemeID){
//                    (*curNode)->disp_vec.push_back(pair<unsigned int, unsigned int>(_timeDivergence-1, firstDemeID));
//                }
//            }
//            else (*curNode)->disp_vec.push_back(pair<unsigned int, unsigned int>(_timeDivergence-1, firstDemeID));
//		}
//        
//		firstDeme->get_chainNodeList().insert(firstDeme->get_chainNodeList().end(), // add the lineages to the first deme
//                                              curDeme->second->get_chainNodeList().begin(),
//                                              curDeme->second->get_chainNodeList().end());
//		curDeme->second->clear_chainNodeList(); // remove the lineages for safety reasons
//	}
//	firstDeme->update_lineages(); // set the number of lineages
//    
//	// remove all demes, but the first one
//	curDeme = _demes.begin();
//	_demes.erase(++curDeme, _demes.end());
//    
//	// adapt the pop size of the first deme (number of lineages, deme size)
//	unsigned int popSize = (unsigned int) get_parameter("divergence_pop_size")->get_value();
//	if (popSize) {                                      // pop size is specified
//		if (2 * popSize < firstDeme->get_lineages()) firstDeme->set_deme_size((unsigned int) ((double) firstDeme->get_lineages() / 2)); // pop size is too small
//		else firstDeme->set_deme_size(popSize);
//	}
//	else {                // pop size is not specified => adjust to metapop size
//		// get the pop size of the entire metapopulation at generation 1
//		unsigned int sumPopSize = 0;
//		map<PATCH_ID, dbPop>::iterator curDeme, endDeme = _dbCoalescence[0].get_patches().end(); // the database with the pop sizes
//		for (curDeme = _dbCoalescence[0].get_patches().begin(); curDeme != endDeme; ++curDeme) {
//			sumPopSize += curDeme->second.get_popSize();
//		}
//		firstDeme->set_deme_size(sumPopSize); // merged pop size is the sum of all pop sizes at t=0
//	}
//	assert(2*firstDeme->get_deme_size() >= firstDeme->get_lineages());
//}
//
//// ------------------------------------------------------------------------------
///** merge all remaining lineages in a single (first) deme
// * the new pop size may be specified by the parameter "divergence_pop_size"
// * if not specified the new pop size is the sum of all individuals of the metapop at time 0
// * if the new pop size is smaller than the number of lineages, the pop size is readjusted
// */
//void LCE_CoalescenceRecomb::merge_demes(PATCH_ID& recombDemeID, TTNode* recombNode)
//{
//#ifdef _DEBUG
//    cout << "\n  **** divergence time reached (gen " << _timeDivergence << "): merge demes:"  << flush;
//#endif
//    //merge_demes();
//    assert(_curGen>=_timeDivergence);
//    
//    //perform_migration_before();                                 // Step 2: migration
//    //update_oldTree_migration(_timeDivergence, recombDemeID);
//    assert(_demes.size()==1 ||
//           (_demes.size()==2 && _demes.begin()->first==recombDemeID && !_demes.begin()->second->get_lineages()) ||
//           (_demes.size()==2 && _demes.rbegin()->first==recombDemeID && !_demes.rbegin()->second->get_lineages()));
//    
//    // get the deme with the old lineages
//    TDeme* firstDeme = _demes.begin()->second->get_lineages() ? _demes.begin()->second : _demes.rbegin()->second;
//
//	// adapt the pop size of the first deme (number of lineages, deme size)
//	unsigned int popSize = (unsigned int) get_parameter("divergence_pop_size")->get_value();
//	if (popSize) {                                      // pop size is specified
//		if (2 * popSize < firstDeme->get_lineages()) firstDeme->set_deme_size((unsigned int) ((double) firstDeme->get_lineages() / 2)); // pop size is too small
//		else firstDeme->set_deme_size(popSize);
//	}
//	else {                // pop size is not specified => adjust to metapop size
//		// get the pop size of the entire metapopulation at generation 1
//		unsigned int sumPopSize = 0;
//		map<PATCH_ID, dbPop>::iterator curDeme, endDeme = _dbCoalescence[0].get_patches().end(); // the database with the pop sizes
//		for (curDeme = _dbCoalescence[0].get_patches().begin(); curDeme != endDeme; ++curDeme) {
//			sumPopSize += curDeme->second.get_popSize();
//		}
//		firstDeme->set_deme_size(sumPopSize); // merged pop size is the sum of all pop sizes at t=0
//	}
//	assert(2*firstDeme->get_deme_size() >= firstDeme->get_lineages());
//
//    // does the recombining node has to be migrated?
//    if(recombNode && recombDemeID != firstDeme->get_id()){
//#ifdef _DEBUG
//        cout << "\n  **** merge demes: recomb migration node " << recombNode->ID_Node+1
//        << "* (deme " << recombDemeID+1 << " => " << firstDeme->get_id()+1 << ")" << flush;
//#endif
//        recombDemeID = firstDeme->get_id(); // update the recombDemeID
//        
//        // check that the node not just emigrated to it (exception if it is the coalescence entry)
//        if(recombNode->disp_vec_temp.back().first==_timeDivergence-1 && recombNode->disp_vec_temp.size()!=1){
//            recombNode->disp_vec_temp.pop_back();
//            if(recombNode->disp_vec_temp.back().second!=recombDemeID){
//                recombNode->disp_vec_temp.push_back(pair<unsigned int, unsigned int>(_timeDivergence-1, recombDemeID));
//            }
//        }
//        else recombNode->disp_vec_temp.push_back(pair<unsigned int, unsigned int>(_timeDivergence-1, recombDemeID));
//    }
//    
//    // delete empty demes
//    map<PATCH_ID, TDeme*>::iterator curDeme, endDeme;
//    for(curDeme = _demes.begin(), endDeme = _demes.end(); curDeme!=endDeme;){
//        if(curDeme->second->get_lineages()) ++curDeme;
//        else{
//            delete curDeme->second;
//			_demes.erase(curDeme++);
//        }
//    }
//    
//    
//    assert(_demes.size()==1);
//}
//
//// ------------------------------------------------------------------------------
///** in contrast to the function merge_demes() this function just updates the deme size
// */
//void LCE_CoalescenceRecomb::merge_demes_recomb()
//{
//    assert(_curGen==_timeDivergence);
//
//	// adapt the pop size of the first deme (number of lineages, deme size)
//	TDeme* firstDeme = _demes.begin()->second;
//	unsigned int popSize = (unsigned int) get_parameter("divergence_pop_size")->get_value();
//	if (popSize) {                                      // pop size is specified
//		if (2 * popSize < firstDeme->get_lineages()) firstDeme->set_deme_size(
//				(unsigned int) ((double) firstDeme->get_lineages() / 2)); // pop size is too small
//		else firstDeme->set_deme_size(popSize);
//	}
//	else {                // pop size is not specified => adjust to metapop size
//		// get the pop size of the entire metapopulation at generation 1
//		unsigned int sumPopSize = 0;
//		map<PATCH_ID, dbPop>::iterator curDeme, endDeme = _dbCoalescence[0].get_patches().end(); // the database with the pop sizes
//		for (curDeme = _dbCoalescence[0].get_patches().begin(); curDeme != endDeme; ++curDeme) {
//			sumPopSize += curDeme->second.get_popSize();
//		}
//		firstDeme->set_deme_size(sumPopSize); // merged pop size is the sum of all pop sizes at t=0
//	}
//	assert(2*firstDeme->get_deme_size() >= firstDeme->get_lineages());
//}
//
//// ----------------------------------------------------------------------------------------
//// run_coalescence
//// ----------------------------------------------------------------------------------------
//void LCE_CoalescenceRecomb::run_coalescence()
//{
//#ifdef _DEBUG
//	message("\nLCE_CoalescenceRecomb:");
//#endif
//
//	// coalescence simulations
//	_vMRCA.clear();
//	time_t startTime = time(0);
//	unsigned int c, l, nbChrom = _popPtr->get_protoGenome()->getNbChromosome();
//	unsigned int nbLocusTot = _popPtr->get_protoGenome()->get_nb_locus_tot();
//	unsigned int* nbLocus_perChrom = _popPtr->get_protoGenome()->get_nb_locus_per_chromosome(); // cumulative number of loci per chromosome
//	double* locus_position_tot = _popPtr->get_protoGenome()->get_locus_position_tot(FEM);
//	int nbOutput = nbLocusTot / 10;
//	if (!nbOutput) nbOutput = 1;   		// if less than 10 loci are simulated
//	log2console(0, 0, nbLocusTot, nbOutput);	// output to the console
//
//	// generate the data
//	if (_writerPopSizes) _writerPopSizes->write_popSizes();
//	if (_writer) _writer->write_NEXUS_heading();
////	ARRAY::create_2D(_seq, _pCoalBase->get_nbSamples(), nbLocusTot);
//	bool first = true; // informs if the tree has been used (in case not all chromosomes contain loci
//
//	// for all linked loci
//	for (c = 0, l = 0; c < nbChrom; ++c) {   			// for each chromosome
//		for (first = true; l < nbLocus_perChrom[c]; ++l) { // for each locus of this chromosome
//			if (first) { // generate the tree for the first locus of this chromosome
//				if (_writerLineages) {
//					_writerLineages->write_lineages_start(l);
//					build_tree();                                 //////////////
//					_writerLineages->write_lineages_end();
//				}
//				else build_tree_simple(); // faster since it does not check for temporal parameters or any other output
//				first = false;
//			}
//			else {			// recombine tree for all other loci (except first)
//				test_tree();
//                recombine_tree(locus_position_tot[l]-locus_position_tot[l-1]);						//////////////
//			}
//			if (_writer) _writer->write_NEXUS_body(l);			// write coalescence tree if desired
//
//			// generate genetic data (sprinkle mutations)
//			generateGeneticData(_seq, _popPtr->get_protoGenome()->get_locus_tot(), l); //////////////
//			log2console(time(0) - startTime, l + 1, nbLocusTot, nbOutput); // output to the console
//			print_tree(_mrca);
//
//		}
//		test_tree();
//		clear_tree();
//	}
//
//	// for all unlinked loci (placed after the linked loci)
//	for (; l < nbLocusTot; ++l) {
//		// generate the coalescence tree
//		if (_writerLineages) {
//			_writerLineages->write_lineages_start(l);
//			LCE_Coalescence::build_tree();                                         //////////////
//			_writerLineages->write_lineages_end();
//		}
//		else LCE_Coalescence::build_tree_simple(); // faster since it does not check for temporal parameters or any other output
//		if (_writer) _writer->write_NEXUS_body(l); // write coalescence tree if desired
//
//		// sprinkle mutations on the tree
//		generateGeneticData(_seq, _popPtr->get_protoGenome()->get_locus_tot(), l); //////////////
//		clear_tree();
//		log2console(time(0) - startTime, l + 1, nbLocusTot, nbOutput); // output to the console
//	}
//
//}
//
//// ----------------------------------------------------------------------------------------
//// recombine_tree
//// ----------------------------------------------------------------------------------------
///** alter tree following the number of recombinations */
//void LCE_CoalescenceRecomb::recombine_tree(double locus_dist)
//{
//	if (!locus_dist) return;  // no recombinations => tree does not change
//    
//	// draw the number of recombination events and place them on the tree
//	vector<unsigned long> vecRecombs;
//	_treeSize = compute_treeSize();	// total size of the tree
//    // draw the recombination positions if linked loci are present
//    for (int i = get_pop_ptr()->rand().Poisson(_treeSize * locus_dist / 100); i > 0; --i) {
//		vecRecombs.push_back(get_pop_ptr()->rand().Uniform(_treeSize));
//	}
//    if(vecRecombs.empty()) return; // no recombination
//    
//	sort(vecRecombs.begin(), vecRecombs.end()); // sort the vector
//    
//#ifdef _DEBUG
//	cout << "\n\n****** Recombining positions:";
//    cout << "\n    in total " << vecRecombs.size() << " events,";
//#endif
//    
//	// find recursively the recombination positions, i.e. nodes below (towards the present)
//	//    the recombination positions in the container vecNodes
//	_vecNodes.clear();
//	unsigned long curSize = 0;
//	if (_mrca->desc1) _mrca->desc1->get_recombinationPos(curSize, _vecNodes, vecRecombs);
//	if (_mrca->desc2) _mrca->desc2->get_recombinationPos(curSize, _vecNodes, vecRecombs);
//	_vecNodes_pos = _vecNodes.begin(); // iterators to loop through vecNodes
//	_vecNodes_end = _vecNodes.end();
//    
//#ifdef _DEBUG
//    cout << "\n    resulting in the following " << _vecNodes.size() << " effectiv events:";
//	for(;_vecNodes_pos!=_vecNodes_end; ++ _vecNodes_pos){
//		cout << "\n    Generation: " << _vecNodes_pos->first << "; Node: " << _vecNodes_pos->second->ID_Node+1;
//	}
//	cout << flush;
//#endif
//    
//	dbCoalescence* curDB;					// pointer to the actual db entry
//	map<PATCH_ID, vector<TTNode*> > recombDemes;	// recombDemes<PATCH_ID, vector of demes with a recombining node in it>
//	map<PATCH_ID, TDeme*>::iterator cur, end;
//    
//	// set deme sizes (it will done automatically if the tree has to be reconstructed during the demographic simulation)
//	unsigned long timeDiv = (unsigned long) get_parameter("divergence_time")->get_value(); // get the divergence time
//    
//	// perform recombination
//	// restart the building of the tree at the first recombination event (going backwards in time)
//	// then go generation by generation backwards in time and allow the gene of the recombination
//	// to migrate and coalesce to any other gene (note, that the branch of the recombining gene remains
//	// in the game!). When the gene coalesces, generate a new node, reconnect the lower and upper node
//	// of the coalescing gene as well as the one of the recombining gene. All nodes active with a recombining
//	// gene are kept in the container activeNodes. Continue until all recombining genes have re-coalesced.
//	_vecNodes_pos = _vecNodes.begin(); // iterators to loop through vecNodes
//	_vecNodes_end = _vecNodes.end();
//	while(_vecNodes_pos!=_vecNodes_end){
//		// set up simulation
//		_curGen = _vecNodes_pos->first;// tree has to be rebuild after this generation (backward in time)
//		_num_lineages = setup_tree_atGen(_curGen);	// setup the _demes for this generation time
//        
//		if (_curGen > _nbGen) {						// before demographic simulation
//			if (_curGen > timeDiv) merge_demes_recomb();			// before divergence time
//			else perform_popSize_regulation(&_dbCoalescence[0]); 	// before demographic simulation
//		}
//        
//		update_recombDemes(_vecNodes_pos, _vecNodes_end, recombDemes);
//        
//		////////////////////////////////////////////////////////////////////////////
//		// during demographic simulation
//		// before MRCA
//		for (; _curGen < (unsigned long) _nbGen && _curGen <= _tmrca; ++_curGen) { // for each generation (starting with 1) in the db going backward in time
//			curDB = &_dbCoalescence[_nbGen - _curGen - 1]; 	// starting from the back
//			_num_lineages = update_lineages(_curGen);						// update the present tree lineages
//			update_recombDemes(_vecNodes_pos, _vecNodes_end, recombDemes);	// update the recombining demes (recombDemes)
//			remove_empty_demes(recombDemes);
//			perform_popSize_regulationIII(curDB);     			// Step 1: adapt population sizes
//#ifdef _DEBUG
//			recombine_tree_log(recombDemes);				// temp console output for development
//#endif
//			perform_coalescenceIII(recombDemes);			// Step 3: Coalescence round
//			update_recombDemes(_vecNodes_pos, _vecNodes_end, recombDemes);	// has to be removed : TODO
//			if(recombDemes.empty()) break;
//			perform_migration_recomb(curDB, recombDemes);                   	// Step 2: migration
//		}
//		if(recombDemes.empty()) continue;
//        
//		// after MRCA
//		if(_curGen < (unsigned long) _nbGen){
//			_num_lineages = update_for_recomb_after_MRCA(recombDemes);
//			for (; _curGen < (unsigned long) _nbGen; ++_curGen) {// for each generation (starting with 1) in the db going backward in time
//				curDB = &_dbCoalescence[_nbGen - _curGen - 1]; // starting from the back
//				perform_popSize_regulation(curDB);     // Step 1: adapt population sizes
//#ifdef _DEBUG
//				recombine_tree_log(recombDemes);				// temp console output for development
//#endif
//				perform_coalescenceII();                      // Step 3: Coalescence round
//				if(_num_lineages==1) break;
//				perform_migration(curDB);                  	// Step 2: migration
//				temporal_change(_nbGen - (unsigned int)_curGen);
//			}
//			if(_num_lineages==1) continue;
//		}
//        
//        
//		////////////////////////////////////////////////////////////////////////////
//		// before demographic simulation
//		// before MRCA
//		if(!recombDemes.empty()){
//			for (; _curGen < timeDiv && _curGen <= _tmrca && !recombDemes.empty(); ++_curGen) { // for each generation (starting with 1) in the db going backward in time
//				_num_lineages = update_lineages(_curGen);						// update the present tree lineages
//				update_recombDemes(_vecNodes_pos, _vecNodes_end, recombDemes);
//				remove_empty_demes(recombDemes);
//#ifdef _DEBUG
//				recombine_tree_log(recombDemes);				// temp console output for development
//#endif
//				perform_coalescenceIII(recombDemes);		// Step 3: Coalescence round
//				perform_migration_before_recomb();                 // Step 2: migration
//			}
//			if(recombDemes.empty()) continue;
//			if(_curGen < timeDiv) _num_lineages = update_for_recomb_after_MRCA(recombDemes);
//		}
//        
//		// after MRCA
//		if(_curGen < timeDiv && !recombDemes.empty()){
//			for (; _curGen < timeDiv && _num_lineages > 1; ++_curGen) {	// for each generation (starting with 1) in the db going backward in time
//#ifdef _DEBUG
//				recombine_tree_log(recombDemes);				// temp console output for development
//#endif
//				perform_coalescenceII();                // Step 3: Coalescence round
//				perform_migration_before();                 // Step 2: migration
//				cout << endl << _curGen << ": " <<  _demes.size() << " demes; " << _num_lineages << " lineages" << endl;
//			}
//			if(_num_lineages==1) continue;
//		}
//        
//		////////////////////////////////////////////////////////////////////////////
//		// before divergence time
//        
//		// divergence is reached: merge all demes, if the new pop size is not specified the pop
//		// size will be the total number of individual sof the metapop
//		merge_demes();
//        
//		// before MRCA
//		if(!recombDemes.empty()){
//			for (; _curGen < ULONG_MAX && _curGen <= _tmrca && !recombDemes.empty(); ++_curGen) { // for each generation (starting with 1) in the db going backward in time
//				_num_lineages = update_lineages(_curGen);						// update the present tree lineages
//				update_recombDemes(_vecNodes_pos, _vecNodes_end, recombDemes);
//				remove_empty_demes(recombDemes);
//#ifdef _DEBUG
//				recombine_tree_log(recombDemes);				// temp console output for development
//#endif
//				perform_coalescenceIII(recombDemes);		// Step 3: Coalescence round
//			}
//			if(recombDemes.empty()) continue;
//			if(_curGen < ULONG_MAX) _num_lineages = update_for_recomb_after_MRCA(recombDemes);
//		}
//        
//		// after MRCA
//		if(_curGen < ULONG_MAX && _num_lineages > 1){
//			for (; _curGen < ULONG_MAX && _num_lineages > 1; ++_curGen) {// for each generation (starting with 1) in the db going backward in time
//#ifdef _DEBUG
//				recombine_tree_log(recombDemes);				// temp console output for development
//#endif
//				perform_coalescenceII();            // Step 3: Coalescence round
//			}
//		}
//	}
//    // remove the last lineage (the node is already in the three)
//    remove_mrca_from_deme();
//}
//
//// ----------------------------------------------------------------------------------------
//// update_for_recomb_after_MRCA
//// ----------------------------------------------------------------------------------------
///** we have passed the MRCA: set simulation up as for a normal tree:
// * add the recombining nodes to the normal nodes
// *
// */
//unsigned int
//LCE_CoalescenceRecomb::update_for_recomb_after_MRCA(map<PATCH_ID, vector<TTNode*> >& recombDemes)
//{
//#ifdef _DEBUG
//	cout << "\n**** The previous MRCA is reached: set up for normal coalescence";
//
//	// update the number of lineages per deme
//	map<PATCH_ID, TDeme*>::iterator curDeme, endDeme;
//	for (curDeme = _demes.begin(), endDeme = _demes.end(); curDeme != endDeme; ++curDeme) {
//		assert(curDeme->second->get_id()==curDeme->first);
//		cout << "\nDeme " << curDeme->first+1 << ": " << curDeme->second->get_chainNodeList().size() << flush;
//	}
//#endif
//
//	// remove obsolete nodes
//	TTNode *obsolete, *other, *node;
//	TDeme* pDeme;
//
//
//	// add the recombining lineages to _demes (_lineages is updated below)
//	vector<TTNode*>::iterator curNode, endNode;
//	map<PATCH_ID, vector<TTNode*> >::iterator pos, end;
//	unsigned int nbLineages = 0;
//	for(pos=recombDemes.begin(), end=recombDemes.end(); pos!=end; ++pos){
//		// for each recombining node
//		for(curNode=pos->second.begin(), endNode=pos->second.end(); curNode!=endNode; ++curNode){
//			node=*curNode;
//
//			assert(node->ancestor); // it cannot be the MRCA
//
//			// 1. remove obsolete old ancestor
//			obsolete = node->ancestor;
//			other =  (obsolete->desc1==node) ? obsolete->desc2 : obsolete->desc1; // which child is the other one
//			other->ancestor = obsolete->ancestor;
//            if(obsolete->ancestor->desc1==obsolete) obsolete->ancestor->desc1 = other;
//            else obsolete->ancestor->desc2 = other;
//			other->disp_vec.insert(other->disp_vec.end(), obsolete->disp_vec.begin(), obsolete->disp_vec.end());
//			obsolete->desc1 = obsolete->desc2 = obsolete->ancestor = NULL;	// reset all for security reasons
//
//			// 2. move disp_vec_temp to disp_vec
//			node->disp_vec = node->disp_vec_temp;
//			node->disp_vec_temp.clear();
//			node->ancestor = NULL;
//
//			// if the obsolete is the MRCA, replace it by other
//			if(!obsolete->ancestor){
//				// since we are at the tMRCA the other node will be in the same deme,
//				// just a single node, namely the MRCA will be present
//				assert(_demes.find(obsolete->disp_vec.front().second)!=_demes.end());	// it has to be present
//				pDeme = _demes[obsolete->disp_vec.front().second];		// it is the first entry
//				assert(pDeme->get_chainNodeList().size()==1);
//				assert(pDeme->get_chainNodeList().front()->ID_Node==obsolete->ID_Node);
//				pDeme->get_chainNodeList().clear();
//				pDeme->get_chainNodeList().push_back(other);
//			}
//		}
//
//		nbLineages += pos->second.size();
//
//		// add the recombining nodes to the old lineages
//		pDeme = get_demes(pos->first);
//		pDeme->get_chainNodeList().insert(get_demes(pos->first)->get_chainNodeList().end(),
//				pos->second.begin(), pos->second.end());
//		pDeme->update_lineages();
//	}
//
//	recombDemes.clear();
//
//	/*
//	 // update the number of lineages per deme
//	 map<PATCH_ID, TDeme*>::iterator curDeme, endDeme;
//	 unsigned int nbLineages = 0;
//	 for (curDeme = _demes.begin(), endDeme = _demes.end(); curDeme != endDeme; ++curDeme) {
//		 nbLineages += curDeme->second->update_get_lineages();
//		 assert(curDeme->second->get_id()==curDeme->first);
//	 }
//	 */
//
//	return nbLineages+1;		// +1: of the old tree only the MRCA is remaining
//}
//
//// ----------------------------------------------------------------------------------------
//// update_for_recomb_after_MRCA
//// ----------------------------------------------------------------------------------------
///** we have passed the MRCA: set simulation up as for a normal tree:
// * add the recombining nodes to the normal nodes
// * remove the current MRCA
// *
// * problme if it goes behind the MRCA:
// *  - recombining node is a child of MRCA: the MRCA has to be removed, thus one has to stop at _curGen_tmrca-1
// *  - recombining node is not a child of MRCA: the MRCA remains a node, thu onbe has to stop at _curGen=_tmrca
// */
//unsigned int
//LCE_CoalescenceRecomb::update_for_recomb_after_MRCA(PATCH_ID recombDemeID,TTNode*& recombNode, unsigned long time)
//{
//	cout << "\n**** The previous MRCA is reached (gen " << time << "): set up for normal coalescence" << flush;
//    
//    // perhaps there are more than one coalescence events at  tMRCA
//    unsigned int nbLins = recombNode->ancestor->ancestor ? 1 : 2;
//    while(get_num_lineages()>nbLins){
//        update_oldTree_coalescence(_tmrca, true);   // it may be that we are only at _tmrca-1 and not at _tmrca
//    }
//    
//#ifdef _DEBUG
//    
//	// update the number of lineages per deme
//	unsigned int nbLineages = 0;
//	map<PATCH_ID, TDeme*>::iterator curDeme, endDeme;
//	for (curDeme = _demes.begin(), endDeme = _demes.end(); curDeme != endDeme; ++curDeme) {
//		assert(curDeme->second->get_id()==curDeme->first);
//        nbLineages += curDeme->second->get_chainNodeList().size();
//		cout << "\nDeme " << curDeme->first+1 << ": " << curDeme->second->get_chainNodeList().size() << flush;
//	}
//    assert(get_num_lineages()==nbLineages);
//    assert((recombNode->ancestor==_mrca && nbLineages==2) || nbLineages==1);
//    assert(time>=_nbGen || recombDemeID == recombNode->get_currentDeme_temp(time));
//#endif
//    
//	// remove obsolete nodes
//	TTNode *obsolete, *other;
//	TDeme* pDeme;
//    
//    // 1. remove obsolete old ancestor:
//    obsolete = recombNode->ancestor;
//    other =  (obsolete->desc1==recombNode) ? obsolete->desc2 : obsolete->desc1; // which child is the other one
//    other->ancestor = obsolete->ancestor;
//
//    vector<pair<unsigned int, unsigned int> >::iterator start = obsolete->disp_vec.begin();
//    if(other->disp_vec.back().second== obsolete->disp_vec.front().second) ++start; // ignore the first entry
//    other->disp_vec.insert(other->disp_vec.end(), start, obsolete->disp_vec.end());
//
//    
//    
//    // if the obsolete is the MRCA, replace it by other
//    if(!obsolete->ancestor){
//        // correct tree size (remove active branche sizes)
//        add2treeSize(-2-(time - other->time)-(time - recombNode->time));
//
//        // since we are at the tMRCA the other node will be in the same deme,
//        // just a single node, namely the MRCA will be present
//        assert(obsolete->disp_vec.size()==1);                       // only a single entry must be present
//        assert(_demes.find(obsolete->disp_vec.front().second)!=_demes.end());	// it has to be present
//        pDeme = _demes[obsolete->disp_vec.front().second];		// it is the first entry
//        if(recombNode->ancestor==_mrca){        // the recombining node is a child of the ancestor: remove MRCA
//            assert(pDeme->get_chainNodeList().size()==2);
//            assert(pDeme->get_chainNodeList().front()->ID_Node==obsolete->desc1->ID_Node ||
//                   pDeme->get_chainNodeList().front()->ID_Node==obsolete->desc2->ID_Node);
//            assert(pDeme->get_chainNodeList().back()->ID_Node==obsolete->desc1->ID_Node ||
//                   pDeme->get_chainNodeList().back()->ID_Node==obsolete->desc2->ID_Node);
//        }
//        else{
//            assert(pDeme->get_chainNodeList().size()==1);
//            assert(pDeme->get_chainNodeList().front()->ID_Node==obsolete->ID_Node);
//        }
//        pDeme->get_chainNodeList().clear();
//        pDeme->get_chainNodeList().push_back(other);
//        pDeme->update_lineages();
//    }
//    else{
//        // correct tree size (remove active branche sizes)
//        add2treeSize(-(obsolete->time - recombNode->time));
//        
//        if(obsolete->ancestor->desc1==obsolete) obsolete->ancestor->desc1 = other;
//        else obsolete->ancestor->desc2 = other;
//    }
//    delete obsolete;        // delete the obsolete node
//    
//    
//    // 2. move disp_vec_temp to disp_vec
//    recombNode->disp_vec = recombNode->disp_vec_temp;
//    recombNode->disp_vec_temp.clear();
//    recombNode->ancestor = NULL;
//    
//    // add the recombining nodes to the old lineages
//    pDeme = get_demes(recombDemeID);
//    pDeme->get_chainNodeList().push_back(recombNode);
//    pDeme->update_lineages();
//    
//    
//    recombNode = NULL;
//    return 2;
//}
//
//// ----------------------------------------------------------------------------------------
//// recombine_tree_log
//// ----------------------------------------------------------------------------------------
//void
//LCE_CoalescenceRecomb::recombine_tree_log(map<PATCH_ID, vector<TTNode*> >& recombDemes)
//{
//	test_tree();
//
//	map<PATCH_ID, vector<TTNode*> >::iterator curRecomb, endRecomb;
//	unsigned int size=0;
//	for(curRecomb=recombDemes.begin(), endRecomb=recombDemes.end(); curRecomb!=endRecomb; ++curRecomb){
//		size += curRecomb->second.size();
//	}
//	bool first;
//
//	cout << "\n\n" << _curGen <<
//			": # recombining demes: " << recombDemes.size() <<
//			"; # lineages: " << _num_lineages <<
//			"; # recombLineages: " << size << endl;
//	map<PATCH_ID, TDeme*>::iterator cur, end;
//
//	for(cur=_demes.begin(); cur!=_demes.end(); ++cur){
//		cout << "  patch: " << cur->first+1 <<
//				"; deme_size: " << cur->second->get_deme_size() <<
//				"; lineages: " << cur->second->get_lineages() << " (";
//
//		vector<TTNode*>::iterator curVec, endVec;
//		for(first = true, curVec=cur->second->get_chainNodeList().begin(),
//				endVec=cur->second->get_chainNodeList().end(); curVec!=endVec; ++curVec){
//			if(!first) cout << "; ";
//			else first = false;
//			cout << (*curVec)->ID_Node+1;
//		}
//
//		cout << "); recombNodes: ";
//
//		map<PATCH_ID, vector<TTNode*> >::iterator curRecombDeme = recombDemes.find(cur->first);
//		if(curRecombDeme == recombDemes.end() || curRecombDeme->second.empty() ) cout << " none";
//		else{
//			cout << curRecombDeme->second.size() << " (";
//			vector<TTNode*>::iterator rec_pos = curRecombDeme->second.begin();
//			vector<TTNode*>::iterator rec_end = curRecombDeme->second.end();
//			for(first=true; rec_pos!=rec_end; ++rec_pos){
//				if(!first) cout << "; ";
//				else first = false;
//				cout << (*rec_pos)->ID_Node+1 << "*";
//			}
//			cout << ")";
//		}
//		cout << endl;
//	}
//}
//
//// ----------------------------------------------------------------------------------------
//// recombine_tree_log
//// ----------------------------------------------------------------------------------------
///** recombDeme is NULL if no recombining Node is present
// */
//void
//LCE_CoalescenceRecomb::recombine_tree_log(PATCH_ID recombDemeID, TTNode* recombNode)
//{
//	//test_tree();
//    
//	cout << "\n\n" << _curGen <<
//    ": # recombining demes: " << (recombNode ? "1" : "0") <<
//    "; # lineages: " << _num_lineages << " + " << (recombNode ? "1'" : "0'") << endl;
//	
//    map<PATCH_ID, TDeme*>::iterator cur, end;
//    bool first;
//	for(cur=_demes.begin(); cur!=_demes.end(); ++cur){
//		cout << "  patch: " << cur->first+1 <<
//        "; deme_size: " << cur->second->get_deme_size() <<
//        "; lineages: " << cur->second->get_lineages() << " + " <<
//        ((cur->first == recombDemeID && recombNode) ? "1'" : "0'") << " (";
//        
//		vector<TTNode*>::iterator curVec, endVec;
//		for(first = true, curVec=cur->second->get_chainNodeList().begin(),
//            endVec=cur->second->get_chainNodeList().end(); curVec!=endVec; ++curVec){
//			if(!first) cout << "; ";
//			else first = false;
//			cout << (*curVec)->ID_Node+1;
//		}
//        
//        if(cur->first == recombDemeID && recombNode){
//            if(cur->second->get_lineages()) cout << "; ";
//            cout << recombNode->ID_Node+1 << "'";
//        }
//        cout << ")" << endl;
//	}
//}
//
//// ----------------------------------------------------------------------------------------
//// update_recombDemes
//// ----------------------------------------------------------------------------------------
///** for each generation update the "recombining demes" (the demes with lineages (nodes) under recombination)
// * recombNodes<demeID, vector of recombining nodes>
// * modified vectors:
// *  - recombDemes
// *  - _demes (are added if a recombining lineage emigrates to a deme not yet having lineages)
// *  Demes without any recombining node are NOT removed
// */
//void LCE_CoalescenceRecomb::update_recombDemes(
//		multimap<unsigned long, TTNode*>::iterator& vecNodes_pos,
//		multimap<unsigned long, TTNode*>::iterator& vecNodes_end,
//		map<PATCH_ID, vector<TTNode*> >& recombDemes)
//{
//	// REMOVE all recombining demes without any recombining node
//	map<PATCH_ID, vector<TTNode*> >::iterator curDeme, endDeme;
//	for (curDeme = recombDemes.begin(), endDeme = recombDemes.end(); curDeme != endDeme;) {
//		if (!curDeme->second.empty()) ++curDeme;
//		else recombDemes.erase(curDeme++);
//	}
//
//	// ADD new become recombining demes to recombDemes and potentially also to _demes
//	unsigned int demeID;
//	TTNode *curNode;
//	assert(vecNodes_pos==vecNodes_end || vecNodes_pos->first >= _curGen);
//	for (; vecNodes_pos->first == _curGen && vecNodes_pos != vecNodes_end; ++vecNodes_pos) {
//		curNode = vecNodes_pos->second;
//		demeID = curNode->get_currentDeme(_curGen);
//		recombDemes[demeID].push_back(vecNodes_pos->second);
//		get_demes(demeID);		// if not present insert the deme to _demes
//		curNode->recombTime = RECOMB;
//
//		// copy the dispersal vector before the recombination event
//		// find first up to where the vector has to be copied (i.e. the lower part)
//		assert(curNode->disp_vec_temp.empty());
//        vector<pair<unsigned int, unsigned int> >::iterator cur = curNode->disp_vec.begin();
//        vector<pair<unsigned int, unsigned int> >::iterator end = curNode->disp_vec.end();;
//        while(cur!=end && cur->first<=_curGen) {
//            ++cur;
//        }
//        curNode->disp_vec_temp.insert(curNode->disp_vec_temp.begin(), curNode->disp_vec.begin(), cur);	// copy the part of the vector (the last element (cur) is not added!!!)
//        assert(curNode->disp_vec_temp.back().first<_curGen);
//	}
//}
//
//// ----------------------------------------------------------------------------------------
//// update_recombDemes
//// ----------------------------------------------------------------------------------------
///** initalize the recomination of the given node
// * check if the deme holding the linage is inthe _demes list
// * copy the part of the disp_vector before recombination to the temp vector
// */
//void LCE_CoalescenceRecomb::init_recomb(PATCH_ID recombDemeID, TTNode* recombNode, unsigned long time)
//{
//    get_demes(recombDemeID);		// if not present insert the deme to _demes
//    recombNode->recombTime = RECOMB;
//    
//    // copy the dispersal vector before the recombination event
//    // find first up to where the vector has to be copied (i.e. the lower part)
//    assert(recombNode->disp_vec_temp.empty());
//    assert(recombNode->disp_vec.front().first<=time);
//    assert(recombNode->ancestor->time>time);
//    vector<pair<unsigned int, unsigned int> >::iterator cur = recombNode->disp_vec.begin();
//    vector<pair<unsigned int, unsigned int> >::iterator end = recombNode->disp_vec.end();;
//	while(cur!=end && cur->first<=time) {
//		++cur;
//	}
//    recombNode->disp_vec_temp.insert(recombNode->disp_vec_temp.begin(), recombNode->disp_vec.begin(), cur);
//    assert(recombNode->disp_vec_temp.back().first<=time);
//}
//
//// ----------------------------------------------------------------------------------------
//// set_up_metapop
//// ----------------------------------------------------------------------------------------
///**set up _demes which contains all demes with lineages  at time curGen (time starts at the present)
// * thereby, find first all nodes having the ancestor passing the time curGen
// * returns the total number of lineages
// */
//unsigned int LCE_CoalescenceRecomb::setup_tree_atGen(unsigned long curGen)
//{
//	// search all nodes for the given time add them to the corresponding demes
//	_demes.clear();             // map<PATCH_ID, TDeme*>
//	TTNode* curNode;
//	map<unsigned int, TTNode*> tempList; // temporal vector with all nodes for the given time
//	for (unsigned int i = 0; i < _pCoalBase->get_nbSamples(); ++i) {		// for each sampled node
//
//		// find the nodes below (towards the present) the given recombination time
//		curNode = _sampleNodes[i];
//		while ((unsigned int) curNode->ancestor->time <= curGen) {
//			curNode = curNode->ancestor;
//		}
//
//		// add the node only if not yet present in the list
//		if (tempList.find(curNode->ID_Node) == tempList.end()) {
//			tempList[curNode->ID_Node] = curNode;
//			assert(curNode->get_currentDeme(curGen)!=my_NAN);
//
//			// insert the deme if not present and add the node
//			get_demes(curNode->get_currentDeme_afterMigration(curGen))->get_chainNodeList().push_back(curNode);
//		}
//	}
//	tempList.clear();
//
//	// update the number of lineages per deme
//	map<PATCH_ID, TDeme*>::iterator curDeme, endDeme;
//	unsigned int nbLineages = 0;
//	for (curDeme = _demes.begin(), endDeme = _demes.end(); curDeme != endDeme; ++curDeme) {
//		nbLineages += curDeme->second->update_get_lineages();;
//		assert(curDeme->second->get_id()==curDeme->first);
//	}
//
//	return nbLineages;
//}
//
//// ----------------------------------------------------------------------------------------
//// update_demes
//// ----------------------------------------------------------------------------------------
///**the "real" lineages of the "old" tree are updated for migration and coalescence
// * the container _demes and its chainNodeList is updated
// * returned is the total number of lineages (of the "old" tree)
// * add, but DON'T remove empty demes from _demes, since they may contain newly
// *     recombining lineages (is done later)
// */
//unsigned int
//LCE_CoalescenceRecomb::update_lineages(unsigned long time)
//{
//	map<PATCH_ID, TDeme*>::iterator curDeme, endDeme;
//	map<PATCH_ID, TDeme*> newDemes; // container to add newly colonized demes
//	vector<TTNode*>::iterator curNode, endNode;
//	TDeme* pDeme;
//	for (curDeme = _demes.begin(), endDeme = _demes.end(); curDeme != endDeme; ++curDeme) {
//		pDeme = curDeme->second;
//		// for each current lineage in this deme
//		for(curNode=pDeme->get_chainNodeList().begin(),
//				endNode=pDeme->get_chainNodeList().end();
//				curNode!=endNode;){
// 			if((*curNode)->ancestor->time==time){							// check if the node coalesced
//				// add the ancestor node only if it is desc1 (to avoid adding it twice)
//				if((*curNode)->ancestor->desc1==*curNode){
//					pDeme->get_chainNodeList_temp().push_back((*curNode)->ancestor);	// deme does not change
//					(*curNode)->ancestor->get_currentDeme(time);			// set iterators
//				}
//				curNode = pDeme->get_chainNodeList().erase(curNode);		// remove the node from the old deme
//				endNode = pDeme->get_chainNodeList().end();					// end has to be reset!!!!
//			}
//			else{															// check if the node migrated
//				// we have to check the next position to see if a migration happened
//				if((*curNode)->disp_vec_next!=(*curNode)->disp_vec.end()
//						&& (*curNode)->disp_vec_next->first==time){
//					TDeme* newDeme = get_demes((*curNode)->get_currentDeme(time, (*curNode)->disp_vec_next), newDemes);
//					newDeme->get_chainNodeList_temp().push_back(*curNode); 	// add the node to the new deme
//					curNode = pDeme->get_chainNodeList().erase(curNode);	// remove the node from the old deme
//					endNode = pDeme->get_chainNodeList().end();				// end has to be reset!!!!
//				}
//				else {
//					++curNode;						// no change
//				}
//			}
//		}
//	}
//
//	// add the newly colonized demes to _demes
//	if (!newDemes.empty()) _demes.insert(newDemes.begin(), newDemes.end()); // merge the containers
//
//	// update the number of lineages per deme
//	unsigned int nbLineages = 0;
//	for (curDeme = _demes.begin(), endDeme = _demes.end(); curDeme != endDeme; ++curDeme) {
//		nbLineages += curDeme->second->merge_chainNodeLists();
//	}
//
//	return nbLineages;
//}
//
//
//
//// ----------------------------------------------------------------------------------------
//// update_demes
//// ----------------------------------------------------------------------------------------
///**the old tree is updated for coalescence events
// * _demes does not change, but chainNodeList and _num_lineages
// * returned is the total number of lineages (of the "old" tree)
// * no demes are added or removed!
// */
//void
//LCE_CoalescenceRecomb::update_oldTree_coalescence(unsigned long time, bool onlyOnce)
//{
//	map<PATCH_ID, TDeme*>::iterator curDeme, endDeme;
//	vector<TTNode*>::iterator curNode, endNode;
//    bool multiCoal=false;       // if mutliple coalescence events occur at the same generation
//	TDeme* pDeme;
//	for (curDeme = _demes.begin(), endDeme = _demes.end(); curDeme != endDeme; ++curDeme) {
//		pDeme = curDeme->second;
//		// for each current lineage in this deme
//		for(curNode=pDeme->get_chainNodeList().begin(),
//            endNode=pDeme->get_chainNodeList().end();
//            curNode!=endNode;){
// 			if((*curNode)->ancestor && (*curNode)->ancestor->time==time){   // check if the node coalesced
//                vector<TTNode*>::iterator desc1_pos = pDeme->getNodePosInChainNodeList((*curNode)->ancestor->desc1);
//                vector<TTNode*>::iterator desc2_pos = pDeme->getNodePosInChainNodeList((*curNode)->ancestor->desc2);
//                
//                if(desc1_pos != endNode && desc2_pos != endNode){
//					// add the new node to the deme
//                    pDeme->get_chainNodeList_temp().push_back((*curNode)->ancestor);	// deme does not change
//					(*curNode)->ancestor->get_currentDeme(time);			// set iterators
//                    --_num_lineages;
//#ifdef _DEBUG
//                    cout << "\n  **** old tree coalescence of " << (*desc1_pos)->ID_Node+1
//                         << " + " << (*desc2_pos)->ID_Node+1
//                         << " => " << (*curNode)->ancestor->ID_Node+1 << flush;
//#endif
//                    
//                    if((*curNode)->ancestor->ancestor && (*curNode)->ancestor->ancestor->time==time) multiCoal=true;
//
//                    // remove both nodes from the old deme (remove first the node towards the end of the vector)
//                    if(desc1_pos>desc2_pos) swap(desc1_pos, desc2_pos);
//                    curNode = pDeme->get_chainNodeList().erase(desc2_pos);
//                    curNode = pDeme->get_chainNodeList().erase(desc1_pos);
//                    endNode = pDeme->get_chainNodeList().end();					// end has to be reset!!!!
//				}
//                else ++curNode;
//			}
//			else ++curNode;						// no change
//		}
//	}
//    
//	// update the number of lineages per deme
//	for (curDeme = _demes.begin(), endDeme = _demes.end(); curDeme != endDeme; ++curDeme) {
//		curDeme->second->merge_chainNodeLists();
//	}
//    
//    if(multiCoal && !onlyOnce) update_oldTree_coalescence(time); // if multiple coalescences repeat this step
//}
//
//
//
//// ----------------------------------------------------------------------------------------
//// update_demes
//// ----------------------------------------------------------------------------------------
///**the old tree is updated for migration events
// * _dmes and also chainNodeList may change, but not _num_lineages
// * add, but DON'T remove empty demes from _demes, since they may contain newly
// *     recombining lineages (is done later)
// */
//void
//LCE_CoalescenceRecomb::update_oldTree_migration(unsigned long time, PATCH_ID recombDemeID)
//{
//	map<PATCH_ID, TDeme*>::iterator curDeme, endDeme;
//	map<PATCH_ID, TDeme*> newDemes; // container to add newly colonized demes
//	vector<TTNode*>::iterator curNode, endNode;
//	TDeme* pDeme;
//	for (curDeme = _demes.begin(), endDeme = _demes.end(); curDeme != endDeme; ++curDeme) {
//		pDeme = curDeme->second;
//		// for each current lineage in this deme
//		for(curNode=pDeme->get_chainNodeList().begin(),
//            endNode=pDeme->get_chainNodeList().end();
//            curNode!=endNode;){
//            // check if the node migrated
//            // we have to check the next position to see if a migration happened
//            if((*curNode)->disp_vec_next!=(*curNode)->disp_vec.end()
//               && (*curNode)->disp_vec_next->first==time){
//                TDeme* newDeme = get_demes((*curNode)->get_currentDeme(time, (*curNode)->disp_vec_next), newDemes);
//                newDeme->get_chainNodeList_temp().push_back(*curNode); 	// add the node to the new deme
//#ifdef _DEBUG
//                cout << "\n  **** old tree migration node " << (*curNode)->ID_Node+1
//                << " (deme " << pDeme->get_id()+1
//                << " => " << newDeme->get_id()+1 << ")" << flush;
//#endif
//                curNode = pDeme->get_chainNodeList().erase(curNode);	// remove the node from the old deme
//                endNode = pDeme->get_chainNodeList().end();				// end has to be reset!!!!
//            }
//            else ++curNode;						// no change
//		}
//	}
//    
//	// add the newly colonized demes to _demes
//	if (!newDemes.empty()) _demes.insert(newDemes.begin(), newDemes.end()); // merge the containers
//    
//	// for each deme (needs a second round!!!, ednDeme iterator remains always valid)
//	for (curDeme = _demes.begin(), endDeme = _demes.end(); curDeme != endDeme;){
//        curDeme->second->merge_chainNodeLists();                        // update the number of lineages per deme
//		
//        if(curDeme->second->get_lineages()){                            // check if it has lineages
//            ++curDeme;
//            continue;
//        }
//        
//        if(curDeme->first==recombDemeID){                               // check if it has recombining nodes
//			++curDeme;
//			continue;
//		}
//        
//		_demes.erase(curDeme++);                                        // it is not anymore "populated" => remove it
//	}
//}
//
//
//// ----------------------------------------------------------------------------------------
//// remove_empty_demes
//// ----------------------------------------------------------------------------------------
///** function checks _demes and removes all demes which have no old AND no new lineage.
// */
//void
//LCE_CoalescenceRecomb::remove_empty_demes(map<PATCH_ID, vector<TTNode*> >& recombDemes)
//{
//	map<PATCH_ID, TDeme*>::iterator curDeme, endDeme;
//	map<PATCH_ID, vector<TTNode*> >::iterator curRecomb = recombDemes.begin();
//	map<PATCH_ID, vector<TTNode*> >::iterator endRecomb = recombDemes.end();
//	for (curDeme = _demes.begin(), endDeme = _demes.end(); curDeme != endDeme;) {
//		if(curDeme->second->get_lineages()){							// it has for sure old lineages
//			if(curRecomb!=endRecomb && curDeme->first==curRecomb->first) ++curRecomb;	// it has also new lineages
//			++curDeme;
//			continue;
//		}
//		if(curRecomb!=endRecomb && curDeme->first==curRecomb->first){					// it has new lineages
//			++curRecomb;
//			++curDeme;
//			continue;
//		}
//
//		// it is not anymore "populated" => remove it
//		_demes.erase(curDeme++);
//	}
//
//}
//
//
//// ----------------------------------------------------------------------------------------
//// remove_empty_demes
//// ----------------------------------------------------------------------------------------
///** function checks _demes and removes all demes which have no old AND no new lineage.
// */
//void
//LCE_CoalescenceRecomb::remove_empty_demes(PATCH_ID recombDemeID)
//{
//	map<PATCH_ID, TDeme*>::iterator curDeme, endDeme;
//	for (curDeme = _demes.begin(), endDeme = _demes.end(); curDeme != endDeme;) {
//		// check if it has lineages
//        if(curDeme->second->get_lineages()){
//            ++curDeme;
//            continue;
//        }
//
//        // check if it has recombining nodes
//        if(curDeme->first==recombDemeID){
//			++curDeme;
//			continue;
//		}
//        
//		// it is not anymore "populated" => remove it
//		_demes.erase(curDeme++);
//	}
//    
//}
//
//
//// ----------------------------------------------------------------------------------------
//// ----------------------------------------------------------------------------------------
//// recombine_tree
//// ----------------------------------------------------------------------------------------
///** alter tree following the number of recombinations
// * The drawn recombination events are performed sequentially
// * during demographic sims the deme migrations are stored in disp_vec
// * divergence has no entry in disp_vec
// * in a generation coalescence occurs before migration!!!
// */
//void LCE_CoalescenceRecombII::recombine_tree(double locus_dist)
//{
//    if (!locus_dist) return;  // no recombinations => tree does not change
//    
//#ifdef _DEBUG
//    cout << "\n\n****** Recombining positions:";
//#endif
//    
//    double recomb_rate = cM2r(locus_dist); // convert cM to r
//    double cur_pos=0;
//    vector<unsigned long> vecRecombs;
//    _treeSize = compute_treeSize();	// total size of the tree
//    unsigned long recombSize, curSize;   // the
//    TTNode* recombNode;
//    unsigned long lastGen;
//    map<PATCH_ID, TDeme*>::iterator cur, end;
//    _timeDivergence = (unsigned long) get_parameter("divergence_time")->get_value(); // get the divergence time
//    if(_timeDivergence<_nbGen) _timeDivergence = _nbGen; // divergence is like a migration event after coalescence event
//    
//    // perform recombination
//    // restart the building of the tree at the first recombination event (going backwards in time)
//    // then go generation by generation backwards in time and allow the gene of the recombination
//    // to migrate and coalesce to any other gene (note, that the branch of the recombining gene remains
//    // in the game!). When the gene coalesces, generate a new node, reconnect the lower and upper node
//    // of the coalescing gene as well as the one of the recombining gene. All nodes active with a recombining
//    // gene are kept in the container activeNodes. Continue until all recombining genes have re-coalesced.
//    for(unsigned int i=0; 1; ++i){
//        test_tree(true);
//        print_tree(_mrca);
//        
//        // draw the nex recombination position
//        cur_pos += get_pop_ptr()->rand().Exponential(_treeSize*recomb_rate);
//        if(cur_pos>recomb_rate) break;   // stop if the length is reached
//        
//        // draw the recombination position on the tree
//        recombSize = get_pop_ptr()->rand().Uniform(_treeSize);
//        curSize = 0;       // cumulative size
//        if (_mrca->desc1) recombNode = _mrca->desc1->get_recombinationPosII(curSize, recombSize, _curGen);
//        if (!recombNode && _mrca->desc2) recombNode = _mrca->desc2->get_recombinationPosII(curSize, recombSize, _curGen);
//        assert(recombNode);
//        
//        // if recombination happenes on time 0 after a coalescence event, this is problematic
//        if(recombNode->time==_curGen) ++_curGen;
//        
//#ifdef _DEBUG
//        cout << "\n\n    " << i+1 << ". Recombination at generation: " << _curGen << "; Node: "
//        << recombNode->ID_Node+1 << " (tMRCA: " << _tmrca
//        << "; treeLength: "<< get_treeSize() << ")" << flush;
//#endif
//        
//        // set up simulation (for _curGen-1 as the timestamp is at the end of the "generation events")
//        _num_lineages = setup_tree_atGen(_curGen-1);	// setup the _demes for this generation time
//        PATCH_ID recombDemeID = recombNode->get_currentDeme_afterMigration(_curGen-1);
//        init_recomb(recombDemeID, recombNode, _curGen-1); // init recombination
//        
//        // the time where the sims have to change from before to after MRCA depends on the location of the recombination event
//        // if a child of the MRCA recombines, then change just after Tmrca otherwise just after Tmrca-1 (since after, there is an offset of 1)
//        unsigned long time_removeMRCA = recombNode->ancestor==_mrca ? _tmrca : _tmrca+1;
//        
//        ////////////////////////////////////////////////////////////////////////////
//        // during demographic simulation
//        if(_curGen>=_nbGen){    // set the last population sizes
//            dbCoalescence* curDB = &_dbCoalescence[0];              // starting from the back
//            perform_popSize_regulationIII(curDB);                       // Step 1: adapt population sizes
//        }
//        else{
//            // before MRCA
//            lastGen = _nbGen > time_removeMRCA? time_removeMRCA : _nbGen;
//            for (; _curGen < lastGen; ++_curGen) { // for each generation (starting with 1) in the db going backward in time
//                if(!loopGen_duringDemo_beforeMRCA(recombDemeID, recombNode)) goto endOfLoop;
//            }
//            if(_curGen==time_removeMRCA) _num_lineages = update_for_recomb_after_MRCA(recombDemeID, recombNode, _curGen-1); // is at the end of the generation event
//            
//            // after MRCA
//            for (; _curGen < (unsigned long) _nbGen; ++_curGen) {       // for each generation (starting with 1) in the db going backward in time
//                if(!loopGen_duringDemo_afterMRCA(recombDemeID, recombNode)) goto endOfLoop;
//            }
//            
//            ////////////////////////////////////////////////////////////////////////////
//            // before demographic simulation, but after divergence time
//            // before MRCA
//            if(recombNode){
//                lastGen = _timeDivergence > time_removeMRCA? time_removeMRCA : _timeDivergence;
//                for (; _curGen < lastGen; ++_curGen) { // for each generation (starting with 1) in the db going backward in time
//                    if(!loopGen_afterDivergence_beforeMRCA(recombDemeID, recombNode)) goto endOfLoop;
//                }
//                if(_curGen==time_removeMRCA) _num_lineages = update_for_recomb_after_MRCA(recombDemeID, recombNode, _curGen-1);
//            }
//            
//            // after MRCA
//            for (; _curGen < _timeDivergence; ++_curGen) {     // for each generation (starting with 1) in the db going backward in time
//                if(!loopGen_afterDivergence_afterMRCA(recombDemeID, recombNode)) goto endOfLoop;
//            }
//        }
//        
//        ////////////////////////////////////////////////////////////////////////////
//        // before divergence time
//        
//        // divergence is reached: merge all demes, if the new pop size is not specified the pop
//        // size will be the total number of individuals of the metapop
//        //if(_curGen==_nbGen) perform_migration_before();           // Step 2: migration
//        
//        
//        //assert(_curGen==_timeDivergence);
//        if(recombNode) merge_demes(recombDemeID, recombNode);
//        else merge_demes();
//        
//        // before MRCA
//        if(recombNode){
//            lastGen = ULONG_MAX > time_removeMRCA? time_removeMRCA : ULONG_MAX;
//            for (; _curGen < lastGen; ++_curGen) { // for each generation (starting with 1) in the db going backward in time
//                if(!loopGen_beforeDivergence_beforeMRCA(recombDemeID, recombNode)) goto endOfLoop;
//            }
//            if(_curGen==time_removeMRCA) _num_lineages = update_for_recomb_after_MRCA(recombDemeID, recombNode, _curGen-1);
//        }
//        
//        // after MRCA
//        if(_curGen < ULONG_MAX){
//            for (; _curGen < ULONG_MAX; ++_curGen) {// for each generation (starting with 1) in the db going backward in time
//                if(!loopGen_beforeDivergence_afterMRCA(recombDemeID, recombNode)) goto endOfLoop;
//            }
//        }
//        
//        // remove the last lineage (the node is already in the tree)
//        remove_mrca_from_deme();
//        
//    endOfLoop:
//        {
//            cout << "\nRecombination finished" << endl;
//        }
//    }
//}
//
//// ----------------------------------------------------------------------------------------
//// loopGen_duringDemo_beforeMRCA
//// ----------------------------------------------------------------------------------------
///** during demographic simulation and before MRCA
// */
//bool
//LCE_CoalescenceRecombII::loopGen_duringDemo_beforeMRCA(PATCH_ID& recombDemeID, TTNode*& recombNode)
//{
//    dbCoalescence* curDB = &_dbCoalescence[_nbGen - _curGen - 1];              // starting from the back
//    perform_popSize_regulationIII(curDB);                       // Step 1: adapt population sizes
//#ifdef _DEBUG
//    recombine_tree_log(recombDemeID, recombNode);				// temp console output for development
//#endif
//    perform_coalescenceIV(_demes[recombDemeID], recombNode);	// Step 3: Coalescence round
//    if(!recombNode) return false;
//    update_oldTree_coalescence(_curGen);						// update the present tree lineages
//    
//    perform_migration_recomb(curDB, recombDemeID, recombNode);  // Step 2: migration
//    update_oldTree_migration(_curGen, recombDemeID);
//    return true;
//}
//
//// ----------------------------------------------------------------------------------------
//// loopGen_duringDemo_afterMRCA
//// ----------------------------------------------------------------------------------------
///** during demogprhic simulation and after MRCA
// */
//bool
//LCE_CoalescenceRecombII::loopGen_duringDemo_afterMRCA(PATCH_ID& recombDemeID, TTNode*& recombNode)
//{
//    dbCoalescence* curDB = &_dbCoalescence[_nbGen - _curGen - 1];          // starting from the back
//    perform_popSize_regulation(curDB);                      // Step 1: adapt population sizes
//#ifdef _DEBUG
//    recombine_tree_log(recombDemeID, recombNode);			// temp console output for development
//#endif
//    perform_coalescenceII();                                // Step 3: Coalescence round
//    if(_num_lineages==1) return false;
//    perform_migration(curDB);                               // Step 2: migration
//    temporal_change(_nbGen - (unsigned int)_curGen);
//    return true;
//}
//
//// ----------------------------------------------------------------------------------------
//// loopGen_afterDivergence_beforeMRCA
//// ----------------------------------------------------------------------------------------
///** before demographic simulation, but after divergence time and before MRCA
// */
//bool
//LCE_CoalescenceRecombII::loopGen_afterDivergence_beforeMRCA(PATCH_ID& recombDemeID, TTNode*& recombNode)
//{
//#ifdef _DEBUG
//    recombine_tree_log(recombDemeID, recombNode);			// temp console output for development
//#endif
//    perform_coalescenceIV(_demes[recombDemeID], recombNode);// Step 3: Coalescence round
//    if(!recombNode) return false;
//    update_oldTree_coalescence(_curGen);					// update the present tree lineages
//    
//    perform_migration_before_recomb();                      // Step 2: migration
//    update_oldTree_migration(_curGen, recombDemeID);
//    return true;
//}
//
//// ----------------------------------------------------------------------------------------
//// loopGen_afterDivergence_afterMRCA
//// ----------------------------------------------------------------------------------------
///** before demographic simulation, but after divergence time and after MRCA
// */
//bool
//LCE_CoalescenceRecombII::loopGen_afterDivergence_afterMRCA(PATCH_ID& recombDemeID, TTNode*& recombNode)
//{
//#ifdef _DEBUG
//    recombine_tree_log(recombDemeID, recombNode);				// temp console output for development
//#endif
//    perform_coalescenceII();                                    // Step 3: Coalescence round
//    if(_num_lineages==1) return false;
//    perform_migration_before();                                 // Step 2: migration
//    return true;
//}
//
//// ----------------------------------------------------------------------------------------
//// loopGen_beforeDivergence_beforeMRCA
//// ----------------------------------------------------------------------------------------
///** before demographic simulation, but before divergence time and before MRCA
// */
//bool
//LCE_CoalescenceRecombII::loopGen_beforeDivergence_beforeMRCA(PATCH_ID& recombDemeID, TTNode*& recombNode)
//{
//#ifdef _DEBUG
//    recombine_tree_log(recombDemeID, recombNode);			// temp console output for development
//#endif
//    perform_coalescenceIV(_demes[recombDemeID], recombNode);// Step 3: Coalescence round
//    if(!recombNode) return false;
//    update_oldTree_coalescence(_curGen);					// update the present tree lineages
//    return true;
//}
//
//// ----------------------------------------------------------------------------------------
//// loopGen_beforeDivergence_afterMRCA
//// ----------------------------------------------------------------------------------------
///** before demographic simulation, but before divergence time and after MRCA
// */
//bool
//LCE_CoalescenceRecombII::loopGen_beforeDivergence_afterMRCA(PATCH_ID& recombDemeID, TTNode*& recombNode)
//{
//#ifdef _DEBUG
//    recombine_tree_log(recombDemeID, recombNode);			// temp console output for development
//#endif
//    perform_coalescenceII();                                // Step 3: Coalescence round
//    return (_num_lineages!=1);
//}
//
//// ----------------------------------------------------------------------------------------
//void LCE_CoalescenceRecombII::init(LCE_Coalescence_base* pBase)
//{
//    _pCoalBase = pBase;
//    _model_threshold = _pCoalBase->get_model_threshold();
//    _dbCoalescence = _pCoalBase->get_dbCoalescence();
//
//    switch (_popPtr->get_pDisperse_LCE()->get_disp_model()) {
//        case 0: // island model
//        case 1: // island propagule model
//            func_ptr_evacuate_lineages = &LCE_Coalescence::evacuate_lineages_island;
//            break;
//        case 2: // 1D stepping stone model
//        case 3: // 2D stepping stone model
//            func_ptr_evacuate_lineages = &LCE_Coalescence::evacuate_lineages_stepping_stone;
//            break;
//    }
//    
//    set_dispRates();
//
//    if (_dbCoalescence) {
//        delete[] _dbCoalescence;
//        _dbCoalescence = NULL;
//    }
//    _dbCoalescence = new dbCoalescence[_nbGen];
//    
//    if (_model_threshold == 0) {
//        func_ptr_perform_coalescence = &LCE_Coalescence::perform_coalescence_multipleEvent;
//        func_ptr_perform_coalescenceII = &LCE_CoalescenceRecomb::perform_coalescence_multipleEventII;
//        func_ptr_perform_coalescenceIV = &LCE_CoalescenceRecomb::perform_coalescence_multipleEventIV;
//    }
//    else if (_model_threshold == 1e6) {
//        func_ptr_perform_coalescence = &LCE_Coalescence::perform_coalescence_singleEvent;
//        func_ptr_perform_coalescenceII = &LCE_CoalescenceRecomb::perform_coalescence_singleEventII;
//        func_ptr_perform_coalescenceIV = &LCE_CoalescenceRecomb::perform_coalescence_singleEventIV;
//    }
//    else {
//        func_ptr_perform_coalescence = &LCE_Coalescence::perform_coalescence_mixedEvent;
//        func_ptr_perform_coalescenceII = &LCE_CoalescenceRecomb::perform_coalescence_mixedEventII;
//        func_ptr_perform_coalescenceIV = &LCE_CoalescenceRecomb::perform_coalescence_mixedEventIV;
//    }
//}
//
//// ----------------------------------------------------------------------------------------
//// recombine_tree
//// ----------------------------------------------------------------------------------------
///** alter tree following the number of recombinations
// * The drawn recombination events are performed sequentially
// */
//void LCE_CoalescenceRecombII::recombine_tree2(double locus_dist)
//{
//    if (!locus_dist) return;  // no recombinations => tree does not change
//    
//#ifdef _DEBUG
//    cout << "\n\n****** Recombining positions:";
//#endif
//    
//    double recomb_rate = cM2r(locus_dist); // convert cM to r
//    double cur_pos=0;
//    vector<unsigned long> vecRecombs;
//    _treeSize = compute_treeSize();	// total size of the tree
//    unsigned long recombSize, curSize;   // the
//    TTNode* recombNode;
//    dbCoalescence* curDB;					// pointer to the actual db entry
//    map<PATCH_ID, TDeme*>::iterator cur, end;
//    _timeDivergence = (unsigned long) get_parameter("divergence_time")->get_value(); // get the divergence time
//    if(_nbGen>_timeDivergence) _timeDivergence = _nbGen;
//    
//    // perform recombination
//    // restart the building of the tree at the first recombination event (going backwards in time)
//    // then go generation by generation backwards in time and allow the gene of the recombination
//    // to migrate and coalesce to any other gene (note, that the branch of the recombining gene remains
//    // in the game!). When the gene coalesces, generate a new node, reconnect the lower and upper node
//    // of the coalescing gene as well as the one of the recombining gene. All nodes active with a recombining
//    // gene are kept in the container activeNodes. Continue until all recombining genes have re-coalesced.
//    for(unsigned int i=0; 1; ++i){
//        test_tree(true);
//        print_tree(_mrca);
//        
//        // draw the nex recombination position
//        cur_pos += get_pop_ptr()->rand().Exponential(_treeSize*recomb_rate);
//        if(cur_pos>recomb_rate) break;   // stop if the length is reached
//        
//        // draw the recombination position on the tree
//        recombSize = get_pop_ptr()->rand().Uniform(_treeSize);
//        curSize = 0;       // cumulative size
//        if (_mrca->desc1) recombNode = _mrca->desc1->get_recombinationPosII(curSize, recombSize, _curGen);
//        if (!recombNode && _mrca->desc2) recombNode = _mrca->desc2->get_recombinationPosII(curSize, recombSize, _curGen);
//        assert(recombNode);
//        
//        // if recombination happenes on time 0 after a coalescence event, this is problematic
//        if(recombNode->time==_curGen) ++_curGen;
//        
//#ifdef _DEBUG
//        cout << "\n\n    " << i+1 << ". Recombination at generation: " << _curGen << "; Node: "
//        << recombNode->ID_Node+1 << " (tMRCA: " << _tmrca
//        << "; treeLength: "<< get_treeSize() << ")" << flush;
//#endif
//        
//        // set up simulation (for _curGen-1 as the timestamp is at the end of the "generation events")
//        _num_lineages = setup_tree_atGen(_curGen-1);	// setup the _demes for this generation time
//        PATCH_ID recombDemeID = recombNode->get_currentDeme_afterMigration(_curGen-1);
//        init_recomb(recombDemeID, recombNode, _curGen-1); // init recombination
//        
//        //		if (_curGen > _nbGen) {						// before demographic simulation
//        //	if (_curGen > timeDiv) merge_demes_recomb();			// before divergence time
//        //	else perform_popSize_regulation(&_dbCoalescence[0]); 	// before demographic simulation
//        //}
//        
//        // problme it if goes behind the MRCA:
//        //  - recombining node is a child of MRCA: the MRCA has to be removed, thus one has to stop at _curGen<_tmrca
//        //  - recombining node is not a child of MRCA: the MRCA remains a node, thu onbe has to stop at _curGen==_tmrca
//        ////////////////////////////////////////////////////////////////////////////
//        // during demographic simulation
//        // before MRCA
//        for (; _curGen < (unsigned long) _nbGen &&
//             (_curGen<_tmrca || (_curGen==_tmrca && recombNode->ancestor!=_mrca)); // stop before tMRCA if the child of the MRCA is recombining, otherwise
//             ++_curGen) { // for each generation (starting with 1) in the db going backward in time
//            curDB = &_dbCoalescence[_nbGen - _curGen - 1];              // starting from the back
//            perform_popSize_regulationIII(curDB);                       // Step 1: adapt population sizes
//#ifdef _DEBUG
//            recombine_tree_log(recombDemeID, recombNode);				// temp console output for development
//#endif
//            perform_coalescenceIV(_demes[recombDemeID], recombNode);	// Step 3: Coalescence round
//            if(!recombNode) goto endOfLoop;
//            update_oldTree_coalescence(_curGen);						// update the present tree lineages
//            
//            perform_migration_recomb(curDB, recombDemeID, recombNode);  // Step 2: migration
//            update_oldTree_migration(_curGen, recombDemeID);
//        }
//        
//        // after MRCA
//        if(_curGen < (unsigned long) _nbGen){
//            _num_lineages = update_for_recomb_after_MRCA(recombDemeID, recombNode, _curGen);
//            for (; _curGen < (unsigned long) _nbGen; ++_curGen) {       // for each generation (starting with 1) in the db going backward in time
//                curDB = &_dbCoalescence[_nbGen - _curGen - 1];          // starting from the back
//                perform_popSize_regulation(curDB);                      // Step 1: adapt population sizes
//#ifdef _DEBUG
//                recombine_tree_log(recombDemeID, recombNode);			// temp console output for development
//#endif
//                perform_coalescenceII();                                // Step 3: Coalescence round
//                if(_num_lineages==1) goto endOfLoop;
//                perform_migration(curDB);                               // Step 2: migration
//                temporal_change(_nbGen - (unsigned int)_curGen);
//            }
//        }
//        
//        
//        ////////////////////////////////////////////////////////////////////////////
//        // before demographic simulation, but after divergence time
//        // before MRCA
//        if(recombNode){
//            for (; _curGen < _timeDivergence &&
//                 (_curGen<_tmrca || (_curGen==_tmrca && recombNode->ancestor!=_mrca)); // stop before tMRCA if the child of the MRCA is recombining, otherwise
//                 ++_curGen) { // for each generation (starting with 1) in the db going backward in time
//#ifdef _DEBUG
//                recombine_tree_log(recombDemeID, recombNode);			// temp console output for development
//#endif
//                perform_coalescenceIV(_demes[recombDemeID], recombNode);// Step 3: Coalescence round
//                if(!recombNode) goto endOfLoop;
//                update_oldTree_coalescence(_curGen);					// update the present tree lineages
//                
//                perform_migration_before_recomb();                      // Step 2: migration
//                update_oldTree_migration(_curGen, recombDemeID);
//            }
//            if(_curGen < _timeDivergence) _num_lineages = update_for_recomb_after_MRCA(recombDemeID, recombNode, _curGen);
//        }
//        
//        // after MRCA
//        for (; _curGen < _timeDivergence; ++_curGen) {     // for each generation (starting with 1) in the db going backward in time
//#ifdef _DEBUG
//            recombine_tree_log(recombDemeID, recombNode);				// temp console output for development
//#endif
//            perform_coalescenceII();                                    // Step 3: Coalescence round
//            if(_num_lineages==1) goto endOfLoop;
//            perform_migration_before();                                 // Step 2: migration
//        }
//        
//        ////////////////////////////////////////////////////////////////////////////
//        // before divergence time
//        
//        // divergence is reached: merge all demes, if the new pop size is not specified the pop
//        // size will be the total number of individuals of the metapop
//        if(_curGen==_nbGen) perform_migration_before();                                 // Step 2: migration
//        if(recombNode) merge_demes(recombDemeID, recombNode);
//        else merge_demes();
//        
//        // before MRCA
//        if(recombNode){
//            for (; _curGen < ULONG_MAX &&
//                 (_curGen<_tmrca || (_curGen==_tmrca && recombNode->ancestor!=_mrca)); // stop before tMRCA if the child of the MRCA is recombining, otherwise
//                 ++_curGen) { // for each generation (starting with 1) in the db going backward in time
//#ifdef _DEBUG
//                recombine_tree_log(recombDemeID, recombNode);			// temp console output for development
//#endif
//                perform_coalescenceIV(_demes[recombDemeID], recombNode);// Step 3: Coalescence round
//                if(!recombNode) goto endOfLoop;
//                update_oldTree_coalescence(_curGen);					// update the present tree lineages
//            }
//            if(_curGen < ULONG_MAX) _num_lineages = update_for_recomb_after_MRCA(recombDemeID, recombNode, _curGen);
//        }
//        
//        // after MRCA
//        if(_curGen < ULONG_MAX){
//            for (; _curGen < ULONG_MAX; ++_curGen) {// for each generation (starting with 1) in the db going backward in time
//#ifdef _DEBUG
//                recombine_tree_log(recombDemeID, recombNode);			// temp console output for development
//#endif
//                perform_coalescenceII();                                // Step 3: Coalescence round
//                if(_num_lineages==1) goto endOfLoop;
//            }
//        }
//        
//        // remove the last lineage (the node is already in the tree)
//        remove_mrca_from_deme();
//        
//    endOfLoop:
//        {
//            cout << "\nRecombination finished" << endl;
//        }        
//    }
//}
//
//
