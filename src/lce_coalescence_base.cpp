/** @file lce_coalescence_base.cpp
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
#include "lce_coalescence_base.h"
#include "treplicate.h"
#include "lce_coalescence.h"
#include "stathandler.cpp"
#include "lce_disperse.h"
#include "version.h"


//#include <limits>
//#include <algorithm>

using namespace std;


// ----------------------------------------------------------------------------------------
// constructor
// ----------------------------------------------------------------------------------------
LCE_Coalescence_base::LCE_Coalescence_base(int rank): LCE("coalesce", "coalescence", "", rank)
{
    _seq = NULL;
    _stats = NULL;
    _writer = NULL;
    _writerMRCA = NULL;
    _writerLineages = NULL;
    _writerPopSizes = NULL;
    _dbCoalescence = NULL;
    _vMRCA = NULL;



    ini_paramset();
}

// ----------------------------------------------------------------------------------------
// destructor
// ----------------------------------------------------------------------------------------
LCE_Coalescence_base::~LCE_Coalescence_base()
{
    if (_stats) delete _stats;
    if (_writer) delete _writer;
    if (_writerMRCA) delete _writerMRCA;
    if (_writerLineages) delete _writerLineages;
    if (_writerPopSizes) delete _writerPopSizes;
    if (_vMRCA) delete[] _vMRCA;
    if (_dbCoalescence) delete[] _dbCoalescence;


}

// ----------------------------------------------------------------------------------------
// loadFileServices
// ----------------------------------------------------------------------------------------
void LCE_Coalescence_base::loadFileServices(FileServices* loader)
{
    //writer tree
    unsigned int choice = (unsigned int) get_parameter_value("coalescence_save_tree");
    if (choice) {
        if (!_writer) _writer = new LCE_CoalescenceFH();
        _writer->set((unsigned int) 0,          // occurrence
                     get_parameter("coalescence_tree_dir")->get_arg(),
                     get_parameter("coalescence_tree_filename")->get_arg(),
                     get_parameter("coalescence_tree_script")->get_arg(),
                     (int) 1, // female
                     (int) 1, // adults
                     choice,
                     NULL, "coalescence_tree", ".tree", loader, get_pop_ptr());
        
        _writer->set_pCoal(this);
        loader->attach(_writer); // should not be exectued at each generation
    }
    else if (_writer) {
        delete _writer;
        _writer = NULL;
    }
    
    //writer mrca
    choice = (unsigned int) get_parameter_value("coalescence_save_mrca");
    if (choice) {
        if (!_writerMRCA) _writerMRCA = new LCE_CoalescenceFH_mrca();
        _writerMRCA->set((unsigned int) 0,          // occurrence
                         get_parameter("coalescence_mrca_dir")->get_arg(),
                         get_parameter("coalescence_mrca_filename")->get_arg(),
                         get_parameter("coalescence_mrca_script")->get_arg(),
                         (int) 1, // female
                         (int) 1, // adults
                         choice,
                         NULL, "coalescence_mrca", ".mrca", loader, get_pop_ptr());
        
        _writerMRCA->set_pCoal(this);
        loader->attach(_writerMRCA); // should not be exectued at each generation
    }
    else if (_writerMRCA) {
        delete _writerMRCA;
        _writerMRCA = NULL;
    }
    
    //writer lineages
    choice = (unsigned int) get_parameter_value("coalescence_save_lineages");
    if (choice) {
        if (!_writerLineages) _writerLineages = new LCE_CoalescenceFH_lineages();
        _writerLineages->set(get_parameter("coalescence_lineages_logtime"), // occurrence
                             get_parameter("coalescence_lineages_dir")->get_arg(),
                             get_parameter("coalescence_lineages_filename")->get_arg(),
                             get_parameter("coalescence_lineages_script")->get_arg(),
                             (int) 1,     // female
                             (int) 1,     // adults
                             (unsigned int) get_parameter_value("coalescence_lineages_logtime2"),
                             NULL, "coalescence_lineages", ".lin", loader, get_pop_ptr());
        
        _writerLineages->set_pCoal(this);
        loader->attach(_writerLineages); // should not be exectued at each generation
    }
    else if (_writerLineages) {
        delete _writerLineages;
        _writerLineages = NULL;
    }
    
    //writer pop sizes
    choice = (unsigned int) get_parameter_value("coalescence_save_pop_sizes");
    if (choice) {
        if (!_writerPopSizes) _writerPopSizes = new LCE_CoalescenceFH_popSizes();
        _writerPopSizes->set(get_parameter("coalescence_pop_sizes_logtime"), // occurrence
                             get_parameter("coalescence_pop_sizes_dir")->get_arg(),
                             get_parameter("coalescence_pop_sizes_filename")->get_arg(),
                             get_parameter("coalescence_pop_sizes_script")->get_arg(),
                             (int) 1,     // female
                             (int) 1,     // adults
                             choice,     // format: 2: no output
                             NULL, "coalescence_pop_sizes", ".txt", loader, get_pop_ptr());
        
        _writerPopSizes->set_pCoal(this);
        loader->attach(_writerPopSizes); // should not be exectued at each generation
    }
    else if (_writerPopSizes) {
        delete _writerPopSizes;
        _writerPopSizes = NULL;
    }
}


// ----------------------------------------------------------------------------------------
// loadStatServices
// ----------------------------------------------------------------------------------------
void LCE_Coalescence_base::loadStatServices(StatServices* loader)
{
    LCE_CoalescenceSH ee;
    // allocate the stat handler
    assert(!_stats);
    _stats = new LCE_CoalescenceSH(this);
    loader->attach(_stats);
}

// ----------------------------------------------------------------------------------------
// execute
// ----------------------------------------------------------------------------------------
void LCE_Coalescence_base::execute()
{
    if (_writerPopSizes) _writerPopSizes->write_popSizes();
    if (_writerLineages) _writerLineages->create_lineage_filename();

    if (_nbSamples) {
        ARRAY::create_2D(_seq, get_nbSamples(), _nbLocus);
        
        assert(_popPtr->get_pReplicate()->get_isCoalescence());
        assert(get_pop_ptr());
        assert(get_pop_ptr()->get_protoGenome());
        TGenomeProto* pGenome = get_pop_ptr()->get_protoGenome();
        unsigned int  nbChromosome = pGenome->get_nb_chromosome();
        unsigned int  nbUnlinked = pGenome->get_nb_locus_unlinked();
        
        unsigned int nbLocusThreads = nbChromosome + nbUnlinked;
        
        
        run_coal_locus(this, 0, nbLocusThreads, 0);
        // post coalescence
        delete[] _dbCoalescence; _dbCoalescence = NULL;
        if (_writer) _writer->write_NEXUS();
        if (_writerMRCA) _writerMRCA->write_mrca();

        set_up_metapop();
    }
    else set_up_metapop_empty();
}

// ----------------------------------------------------------------------------------------
void run_coal_locus(LCE_Coalescence_base* pCoalBase, unsigned int from,
                    unsigned int to, unsigned int curThread)
{
    switch (pCoalBase->get_pop_ptr()->get_pReplicate()->get_isCoalescence()) {
        case 1:{
            LCE_Coalescence* pCoal = new LCE_Coalescence(pCoalBase);
            pCoal->run_coalescence(from, to);
            break; // new version here the deme sizes of the demes with lineages are just updated
        }
    }
}

// ----------------------------------------------------------------------------------------
// set_up_metapop
// ----------------------------------------------------------------------------------------
void LCE_Coalescence_base::set_up_metapop()
{
#ifdef _DEBUG
    message("\n\n Set up metapopulation (%i individuals):\n", _nbSamples/2);
#endif
    
    // reset metapopulation after coalescence
    _popPtr->reset_metapopulation();
    _popPtr->set_func_pointer();
    _popPtr->set_func_ptr_noSample_withFull();
    
    // decrement the current generation as it was incremented once too much
    _popPtr->setCurrentGeneration(_popPtr->getCurrentGeneration() - 1);
    
    // set up the metapopulation with the simulated genotypes
    map<unsigned int, unsigned int>::iterator curPos, endPos;
    unsigned char** curSeq = _seq;  // iterator
    unsigned int i;
    Patch* curPatch;
    for (curPos = _iniSamples.begin(), endPos = _iniSamples.end(); curPos != endPos; ++curPos) {
        curPatch = _popPtr->get_vPatch(curPos->first);
        for (i = curPos->second; i > 0; --i, ++curSeq, ++curSeq) { // generate each individual (the seq has to be twice iterated)
            TIndividual* new_ind = _popPtr->makeNewIndividual(NULL, NULL, FEM, curPatch);
            new_ind->create_dadFirst(curSeq);
            curPatch->add(FEM, ADLTx, new_ind);
        }
        _popPtr->new_fullPatch(curPatch);
    }
    _popPtr->add_tempPatch();
    _popPtr->set_sampledInds(ADULTS);
    
    assert(curSeq == _seq+get_nbSamples());
    ARRAY::delete_2D(_seq, get_nbSamples());
}

// ----------------------------------------------------------------------------------------
// set_up_metapop_empty
// ----------------------------------------------------------------------------------------
/** if metapop is empty */
void LCE_Coalescence_base::set_up_metapop_empty()
{
    assert(!_nbSamples);
#ifdef _DEBUG
    message("\n\n Set up metapopulation (%i individuals):\n", _nbSamples/2);
#endif
    
    // reset metapopulation after coalescence
    _popPtr->reset_metapopulation();
    _popPtr->set_func_pointer();
    _popPtr->set_func_ptr_noSample_withFull();
    
    // decrement the current generation as it was incremented once too much
    _popPtr->setCurrentGeneration(_popPtr->getCurrentGeneration() - 1);
    
    // set up the metapopulation with the simulated genotypes
    _popPtr->add_tempPatch();
    _popPtr->set_sampledInds(ADULTS);
}

// ----------------------------------------------------------------------------------------
/** set the sample sizes map (_iniSamples[patchId],size).
 * Caution: sample size is diploid, lineage is haploid!!!
 */
void LCE_Coalescence_base::set_samples(const map<unsigned int, unsigned int>& m, const unsigned int& nb)
{
    _iniSamples = m;
    _nbSamples = 2 * nb; // diploid -> haploid
    assert(_iniSamples.size()==_nbSamples/2);
}

// ----------------------------------------------------------------------------------------
// temporal_change
// ----------------------------------------------------------------------------------------
/** parameters which may change over time */
void LCE_Coalescence_base::temporal_change(const unsigned int& gen)
{
    map<string, Param*>* pParam = _paramSet->getTemporalParams(gen);
    
    // if it is a temporal parameter
    if (pParam) {
        // check if a change has to be made
        map<string, Param*>* pMap = _paramSet->updateTemporalParams(gen);
        if (pMap) {
            // iterate through the map and perform the updates
            map<string, Param*>::iterator pos = pMap->begin();
            for (; pos != pMap->end(); ++pos) {
                if (pos->first == "coalescence_lineages_logtime" && _writerLineages) {
                    if (gen > 1) pos->second->set_arg_for_gen(gen - 1); // as it is a backward in time simulation
                    _writerLineages->set_GenerationOccurrence(
                                                              (unsigned int) pos->second->get_value());
                }
                else if (pos->first == "coalescence_lineages_logtime2" && _writerLineages) {
                    _writerLineages->set_format((unsigned int) pos->second->get_value());
                }
                else if (pos->first == "coalescence_pop_sizes_logtime" && _writerPopSizes) {
                    _writerPopSizes->set_GenerationOccurrence((unsigned int) pos->second->get_value());
                }
            }
        }
    }
}

// ----------------------------------------------------------------------------------------
void LCE_Coalescence_base::init(TMetapop* popPtr, const unsigned int& i)
{
    _popPtr = popPtr;
    _nbGen = i;
    _nbLocus = _popPtr->get_protoGenome()->get_nb_locus_tot();
    
  //  clear_tree();
    if (_writer) _writer->set_GenerationOccurrence(_nbGen); // previously _nbGen is not known
    if (_writerMRCA) _writerMRCA->set_GenerationOccurrence(_nbGen); // previously _nbGen is not known
    
    _model_threshold = get_parameter_value("coalescence_model_threshold");
    
    
    assert(!_dbCoalescence);
    _dbCoalescence = new dbCoalescence[_nbGen];
    assert(!_vMRCA);
    _vMRCA = new unsigned long[_nbLocus];
    
}

// ----------------------------------------------------------------------------------------
// ini_paramset
// ----------------------------------------------------------------------------------------
void LCE_Coalescence_base::ini_paramset()
{
    add_parameter("coalescence_model_threshold", DBL, false, 0, my_NAN, "0.1", false,
                        "The probability threshold above which the approximation p=n(n-1)/(4N) is " \
                        "not applicable as more than 1 coalescence event may occur per generation.",0);
    
    add_parameter("divergence_time", INT2, false, 0, my_NAN, "0", false,
                        "The divergence time of all remaining populations (at this " \
                        "time all remaining popualtions are merged).",0);
    
    add_parameter("divergence_pop_size", INT2, false, 0, my_NAN, "0", false,
                        "The population size before divergence. If not set it is the " \
                        "the new population size is the sum of all remaining populations.",0);
    
    
    add_parameter("coalescence_save_tree", INT2, false, 0, 2, "0", false,
                        "Should the coalescence tree be saved (Nexus format):\n" \
                        "  0: none\n" \
                        "  1: scaled by coalescence time\n" \
                        "  2: scaled by number of mutations",3);
    
    add_parameter("coalescence_tree_dir", STR, false, my_NAN, my_NAN, "", false,
                        "The directory name where the coalescence tree files are stored.",3);
    
    add_parameter("coalescence_tree_filename", STR, false, my_NAN, my_NAN, "", false,
                        "The base name for the coalescence tree files.",3);
    
    add_parameter("coalescence_tree_script", STR, false, my_NAN, my_NAN, "", false,
                        "The script which will be launched after the coalescence tree file is generated.",3);
    
    
    add_parameter("coalescence_save_mrca", INT2, false, 0, 1, "0", false,
                        "Should the MRCA tree be saved:\n" \
                        "  0: none\n" \
                        "  1: output",3);
    
    add_parameter("coalescence_mrca_dir", STR, false, my_NAN, my_NAN, "", false,
                        "The directory name where the MRCA files are stored.",3);
    
    add_parameter("coalescence_mrca_filename", STR, false, my_NAN, my_NAN, "", false,
                        "The base name for the MRCA files.",3);
    
    add_parameter("coalescence_mrca_script", STR, false, my_NAN, my_NAN, "", false,
                        "The script which will be launched after the MRCA file is generated.",3);
    
    add_parameter("coalescence_save_lineages", INT2, false, 0, 1, "0", false,
                        "Should the lineages be saved:\n" \
                        "  0: none\n" \
                        "  1: output",3);
    
    add_parameter("coalescence_lineages_logtime", INT2, false, 0, my_NAN, "1", true,
                        "The time interval of the lineages output during the demographic simulation.",3); // during demographic sims
    
    add_parameter("coalescence_lineages_logtime2", INT2, false, 0, my_NAN, "1", true,
                        "The time interval of the lineages output before (forward-in-time) " \
                        "the demographic simulation.",3); // before demographic sims
    
    add_parameter("coalescence_lineages_dir", STR, false, my_NAN, my_NAN, "", false,
                        "The directory name where the lineages files are stored.",3);
    
    add_parameter("coalescence_lineages_filename", STR, false, my_NAN, my_NAN, "", false,
                        "The base name for the lineages files.",3);
    
    add_parameter("coalescence_lineages_script", STR, false, my_NAN, my_NAN, "", false,
                        "The script which will be launched after the lineages file is generated.",3);
    
    
    add_parameter("coalescence_save_pop_sizes", INT2, false, 0, 1, "0", false,
                        "Should the population sizes per patch be saved:\n" \
                        "  0: none\n" \
                        "  1: output",3);
    
    add_parameter("coalescence_pop_sizes_logtime", INT2, false, 0, my_NAN, "1", true,
                        "The time interval of the population sizes output.",3);
    
    add_parameter("coalescence_pop_sizes_dir", STR, false, my_NAN, my_NAN, "", false,
                        "The directory name where the population size files are stored.",3);
    
    add_parameter("coalescence_pop_sizes_filename", STR, false, my_NAN, my_NAN, "", false,
                        "The base name for the population size files.",3);
    
    add_parameter("coalescence_pop_sizes_script", STR, false, my_NAN, my_NAN, "", false,
                        "The script which will be launched after the population size file is generated.",3);
    
    add_parameter("coalescence_pop_sizes_of_patch", INT_MAT, false, my_NAN, my_NAN, "", false,
                        "Allows defining the patches for with the population sizes are ouputted " \
                        "(by default all patches).",3);// list of patch ids to list (default: all)
    
    
    add_parameter("coalescence_dispersal_rate", INT_MAT, false, my_NAN, my_NAN, "", false,
                        "TODO",5);
    
    add_parameter("dispersal_propagule_prob", DBL, false, 0, 1, "1", true,
                        "TODO",5);
}

// ----------------------------------------------------------------------------------------
/** add the immigrants to the db for coalescence simulations */
void LCE_Coalescence_base::add_immigrants(const unsigned int& gen, const unsigned int& to,
                                     const unsigned int& from, const unsigned int& nb)
{
    _dbCoalescence[gen - 1].add_immigrants(to, from, nb);
}

// ----------------------------------------------------------------------------------------
/** add the current pop sizes to the db for coalescence simulations, but only if it is not empty */
void LCE_Coalescence_base::set_popSize(const unsigned int& gen, const unsigned int& i,
                                  const unsigned int& size)
{
    assert(size);
    _dbCoalescence[gen - 1].set_popSize(i, size);
}
// ----------------------------------------------------------------------------------------
/** add the current pop sizes to the db for coalescence simulations, but only if it is not empty */
void LCE_Coalescence_base::store_popSizes(const unsigned int& curGen)
{
    unsigned int size;
    vector<Patch*>::iterator curPop, endPop = get_pop_ptr()->get_vFullPatch().end();
    for (curPop = get_pop_ptr()->get_vFullPatch().begin(); curPop != endPop; ++curPop) {
        size = (*curPop)->size(FEM, ADLTx);
        if (size) set_popSize(curGen, (*curPop)->get_ID(), size);
    }
}

// ----------------------------------------------------------------------------------------
// write tree recursively
void LCE_CoalescenceFH::write_NEXUS_tree(unsigned int l, LCE_Coalescence* pCoal)
{
    //    // reset the node id's to be smooth
    //    multimap<unsigned int, TTNode*> nodeTimes;
    //	_pCoal->get_mrca()->get_nodeTimes(nodeTimes); // _format: 0: no output, 1: coal time; 2: nbMutations
    //    multimap<unsigned int, TTNode*>::iterator curNode, endNode;
    //    unsigned int i=0;
    //
    //    // jump the sample nodes (id does not change)
    //    for(i=0, curNode=nodeTimes.begin(); i<_pCoal->get_nbSamples(); ++i) ++curNode;
    //
    //    // adjust the id's
    //    for(endNode=nodeTimes.end(); curNode!=endNode; ++curNode){
    //        if(curNode->first==0) continue;     //sample node id's never change
    //        curNode->second->ID_Node = i;
    //        ++i;
    //    }
    
    
    ostringstream FILE;
    pCoal->get_mrca()->write_tree(FILE, _format); // _format: 0: no output, 1: coal time; 2: nbMutations
    _locusTree[l] = FILE.str();
}

// ----------------------------------------------------------------------------------------
/** write the trees to file (trees are stored in _locusTree) */
void LCE_CoalescenceFH::write_NEXUS()
{
    assert(_locusTree);
    
    // open the file
    _NEXUS_filename = get_path() + getReplicateFileName() + get_extension();
#ifdef _DEBUG
    message("LCE_CoalescenceFH::write_NEXUS (%s)\n",_NEXUS_filename.c_str());
#endif
    ofstream FILE(_NEXUS_filename.c_str(), ios::out);
    if (!FILE)
        error("Could not open coalescence tree output file '%s'!\n", _NEXUS_filename.c_str());
    
    char curTime[20];
    time_t tt = time(NULL);
    strftime(curTime, 20, "%d-%m-%Y %H:%M:%S", localtime(&tt));
    
    // write heading
    FILE << "#NEXUS" << "\n\n[Treefile generated by quantiNemo2"
    << "\n    version " << RELEASE
    << "." << REVISION << "." << MINOR_VERSION << TEMP_VERSION << " [" << VERSION_DATE
    << "; " << VERSION_TIME << "]" << "\n    File created the " << curTime;
    
    if (_format == 1) FILE << "\n    Branches scaled by the coalescence time";
    else FILE << "\n    Branches scaled by the number of mutations";
    
    FILE << "\n]" << "\n\nBegin trees;";
    
    // write translation
    FILE << "\n\nTranslate";
    unsigned int patchID, i = 1, j;
    map<unsigned int, unsigned int>::iterator curSample, endSample = _pCoal->get_iniSamples().end(); // the sample sizes
    for (curSample = _pCoal->get_iniSamples().begin(); curSample != endSample; ++curSample) {
        patchID = (unsigned int) curSample->first + 1;
        for (j = 1; j <= curSample->second; ++j, i += 2) { // for each individual (starting at 1)
            FILE << "\n    " << i << "  patch" << patchID << "_ind" << j << "a,";
            FILE << "\n    " << i + 1 << "  patch" << patchID << "_ind" << j << "b";
            FILE << (i == _pCoal->get_nbSamples() - 1 ? ";" : ","); // the last one has to be a ";"
        }
    }
    FILE << "\n";

    
    // write the tree of each locus
    for(unsigned int l=0; l<_pCoal->get_nbLocus(); ++l){
        FILE << "\ntree LOCUS_" << l + 1 << " = " << _locusTree[l] << ";";
    }
    
    // write ending
    FILE << "\n\nEnd;";
    FILE.close();
    
    // call an external program and pass the name of this file
    if (!_script.empty()) execute_script(_script, _NEXUS_filename);
}

// ----------------------------------------------------------------------------------------
void LCE_CoalescenceFH_mrca::write_mrca()
{
    assert(_pCoal->get_vMRCA());
    
    // open the file
    string filename = get_path() + getReplicateFileName() + get_extension();
#ifdef _DEBUG
    message("LCE_CoalescenceFH_mrca::write_mrca (%s)\n", filename.c_str());
#endif
    ofstream FILE(filename.c_str(), ios::out);
    if (!FILE) error("Could not open coalescence mrca output file '%s'!\n", filename.c_str());
    
    vector<unsigned long>::iterator cur, end;
    for (unsigned int l=0; l<_pCoal->get_nbLocus(); ++l){
        FILE << _pCoal->get_vMRCA()[l] << "\n";
    }
    
    FILE.close();
    
    // call an external program and pass the name of this file
    if (!_script.empty()) execute_script(_script, filename);
}

// ----------------------------------------------------------------------------------------
void LCE_CoalescenceFH_lineages::write_lineages_start(unsigned int locus)
{
    // open the file
    _lineage_filename[locus] = get_path() + getReplicateFileName() + "_l" + toStr(locus + 1)
    + get_extension();
#ifdef _DEBUG
    message("LCE_CoalescenceFH_lineages::write_lineages (%s)\n", _lineage_filename[locus].c_str());
#endif
    ofstream FILE(_lineage_filename[locus].c_str(), ios::out);
    if (!FILE)
        error("Could not open coalescence lineages output file '%s'!\n", _lineage_filename[locus].c_str());
    FILE.close();
}

// ----------------------------------------------------------------------------------------
/** write lineages DURING the demographic simulation  */
void LCE_CoalescenceFH_lineages::write_current_lineages(unsigned int nbGen,
                                                        unsigned long curGen,
                                                        dbCoalescence* curDBs,
                                                        LCE_Coalescence* pCoal)
{
    if (!get_GenerationOccurrence()) return;
    unsigned int diff = nbGen - (unsigned int)curGen;       // nbGen > curGen
    if (diff % get_GenerationOccurrence()) return;
    
    // open the file
    ofstream FILE(_lineage_filename[pCoal->get_curLocus()].c_str(), ios::app);
    if (!FILE)
        error("Could not open coalescence lineages output file '%s'!\n",
              _lineage_filename[pCoal->get_curLocus()].c_str());
    
    if (nbGen < curGen) FILE << "-" << diff; // the generation
    else FILE << diff;
    FILE << "   " << pCoal->get_num_lineages();
    
    map<PATCH_ID, dbPop>::iterator curDB, endDB = curDBs->get_patches().end();
    map<PATCH_ID, TDeme*>::iterator curDeme, endDeme = pCoal->get_demes().end();
    for (curDB = curDBs->get_patches().begin(); curDB != endDB; ++curDB) {
        curDeme = pCoal->get_demes().find(curDB->first); // check if this patch has lineages
        FILE << "      " << curDB->first + 1 << "   " << curDB->second.get_popSize() << "   "
        << (curDeme == endDeme ? 0 : curDeme->second->get_lineages());
    }
    FILE << "\n";
    FILE.close();
}

// ----------------------------------------------------------------------------------------
/** write lineages BEFORE the demographic simulation  */
void LCE_CoalescenceFH_lineages::write_current_lineages(unsigned int nbGen,
                                                          unsigned long curGen,
                                                          LCE_Coalescence* pCoal)
{
    if (!get_format()) return;
    unsigned int diff = (unsigned int)curGen - nbGen;       // curGen > nbGen
    if (diff % get_format()) return; // _format is used since we have here two logtimes
    
    // open the file
    ofstream FILE(_lineage_filename[pCoal->get_curLocus()].c_str(), ios::app);
    if (!FILE)
        error("Could not open coalescence lineages output file '%s'!\n",
              _lineage_filename[pCoal->get_curLocus()].c_str());
    
    if (nbGen < curGen) FILE << "-" << diff; // the generation
    else FILE << diff;
    FILE << "   " << pCoal->get_num_lineages();
    
    map<PATCH_ID, TDeme*>::iterator curDeme, endDeme = pCoal->get_demes().end();
    for (curDeme = pCoal->get_demes().begin(); curDeme != endDeme; ++curDeme) {
        FILE << "      " << curDeme->first + 1 << "   " << curDeme->second->get_deme_size() << "   "
        << curDeme->second->get_lineages();
    }
    FILE << "\n";
    FILE.close();
}

// ----------------------------------------------------------------------------------------
// write_lineages_end
// ----------------------------------------------------------------------------------------
void LCE_CoalescenceFH_lineages::write_lineages_end()
{
    // call an external program and pass the name of this file
    if (!_script.empty()) execute_script(_script, _filename);
}

// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
void LCE_CoalescenceFH_popSizes::write_popSizes()
{
    // open the file
    _filename = get_path() + getReplicateFileName() + get_extension();
#ifdef _DEBUG
    message("LCE_CoalescenceFH_popSizes::write_popSizes (%s)\n", _filename.c_str());
#endif
    ofstream FILE(_filename.c_str(), ios::out);
    if (!FILE)
        error("Could not open coalescence population size output file '%s'!\n",
              _filename.c_str());
    
    // are the patches to list specified (default no => list all patches!)
    Param* pParam = _pCoal->get_parameter("coalescence_pop_sizes_of_patch");
    dbCoalescence* curDB;
    map<PATCH_ID, dbPop>::iterator curPop, endPop;
    if (pParam->isSet()) {					// list only the specified patches
        TMatrix* m = pParam->get_as_matrix();		// get the specified patches
        double* array = m->get();                   // index starts at 1
        unsigned int i, curGen, size = m->length();
        
        // write heading
        FILE << "generation";
        for (i = 0; i < size; ++i) {
            if (array[i] > _popPtr->get_nbPatch())
                error("Parameter coalescence_pop_sizes_pops: patch id '%i' is out of range!\n",
                      (int) array[i]);
            FILE << "\t" << "nbIndTot_p" + toStr(array[i]);
        }
        
        // write pop sizes for the specified generation times
        for (curGen = 1; curGen <= _pCoal->get_nbGen(); ++curGen) {	// for each generation (starting with 1) in the db going forward in time
            _pCoal->temporal_change(curGen);
            if (!get_GenerationOccurrence() || curGen % get_GenerationOccurrence()) return;
            FILE << "\n" << curGen;
            curDB = &_pCoal->get_dbCoalescence(curGen - 1);
            endPop = curDB->get_patches().end();
            for (i = 0; i < size; ++i) {
                curPop = curDB->get_patches().find(array[i] - 1);
                if (curPop == endPop) FILE << "\t" << 0;       // patch is empty
                else FILE << "\t" << curPop->second.get_popSize(); // patch is populated
            }
        }
        delete m;
    }
    else {                         // list all patches
        // write heading
        unsigned int i, curGen;
        FILE << "generation";
        for (i = 0; i < _popPtr->get_nbPatch(); ++i) {
            FILE << "\t" << "nbIndTot_p" + toStr(i + 1);
        }
        
        // write pop sizes for the specified generation times
        for (curGen = 1; curGen <= _pCoal->get_nbGen(); ++curGen) {	// for each generation (starting with 1) in the db going forward in time
            _pCoal->temporal_change(curGen);
            if (!get_GenerationOccurrence() || curGen % get_GenerationOccurrence()) return;
            FILE << "\n" << curGen;
            curDB = &_pCoal->get_dbCoalescence(curGen - 1);
            endPop = curDB->get_patches().end();
            for (i = 0; i < _popPtr->get_nbPatch(); ++i) {
                curPop = curDB->get_patches().find(i);
                if (curPop == endPop) FILE << "\t" << 0;       // patch is empty
                else FILE << "\t" << curPop->second.get_popSize(); // patch is populated
            }
        }
    }
    FILE.close();
}

// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
bool LCE_CoalescenceSH::init ( )
{
    StatHandler<LCE_CoalescenceSH>::init();
    return true;
}

// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
bool LCE_CoalescenceSH::setStatRecorders (const string& t)
{
    if (t == "mrca") {
        add("Mean time to the MRCA", "mrca.mean", FLAT, ADULTS, 0,
            &LCE_CoalescenceSH::get_MRCA_mean);
        add("Variance of the time to the MRCA", "mrca.var", FLAT, ADULTS, 0,
            &LCE_CoalescenceSH::get_MRCA_var);
        return true;
    }
    if (t == "mrca.mean")
        return add("Mean time to the MRCA", "mrca.mean", FLAT, ADULTS, 0,
                   &LCE_CoalescenceSH::get_MRCA_mean);
    if (t == "mrca.var")
        return add("Variance of the time to the MRCA", "mrca.var", FLAT, ADULTS, 0,
                   &LCE_CoalescenceSH::get_MRCA_var);
    return false;
}


