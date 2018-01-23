/** @file filehandler.cpp
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


#include "filehandler.h"
#include "tmetapop.h"
#include "version.h"
#include "tsimulation.h"
#include <errno.h>
using namespace std;

// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
bool FileHandler::init()
{
    //check if the occurrences exceed the pop parameters:
    if (_GenerationOccurrence > _popPtr->getGenerations()) {
        _GenerationOccurrence = _popPtr->getGenerations();
    }
    
    //check if the basefilename is already used on disk:
    string filename = get_path() + getGenerationReplicateFileName()
    + get_extension();
    ifstream ifExist;
    ifExist.setstate(ios::failbit);
    
    ifExist.open(filename.c_str(), ios::in);
    
    if (ifExist.is_open() && !_popPtr->get_pSimulation()->getOverwriteFiles()) {
        warning("a file named \"%s\" exists already in folder \"%s\" !!\n",
                (getGenerationReplicateFileName() + get_extension()).c_str(),
                get_path().c_str());
        ifExist.close();
        return false;
    }
    return true;
}

// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
bool FileHandler::set (Param* gen_occ, unsigned int gen_occ_int, string path,
                       string filename, string script, int sex, int age, int format,
                       TTraitProto* trait, string name, string ext,
                       FileServices* loader, TMetapop* ptr)
{
    assert(ptr);
    _popPtr = ptr;
    set_service(loader);
    set_GenerationParam(gen_occ);
    set_GenerationOccurrence(gen_occ ? (unsigned int)gen_occ->get_value() : gen_occ_int);
    set_path(path);
    set_filename(filename);
    set_script(script);
    set_sex(sex);
    set_age(age);
    set_format(format);
    set(trait);
    set_name(name);
    set_extension(ext);
    
    if(_name=="genotype"){
        switch (_format) {
            default: break;
            case 1:
                writeGenotype_func_ptr=&FileHandler::write_Fstat;
                _format = 0;
                break;
            case 2:
                writeGenotype_func_ptr=&FileHandler::write_Fstat;
                _format = 1; // extended
                break;
            case 3:
                writeGenotype_func_ptr=&FileHandler::write_Arlequin;
                _format = 0;
                break;
            case 4:
                writeGenotype_func_ptr=&FileHandler::write_Arlequin;
                _format = 1; // extended
                break;
            case 5:
                writeGenotype_func_ptr=&FileHandler::write_Plink_ped;
                _format = 0;
                write_Plink_map(_format);  // only once
                _quantiTraitIndexes = _popPtr->getTraitIndex("quanti");   // get all quanti traits
                break;
            case 6:
                writeGenotype_func_ptr=&FileHandler::write_Plink_ped;
                _format = 1;    // extended
                write_Plink_map(_format);  // only once
                _quantiTraitIndexes = _popPtr->getTraitIndex("quanti");   // get all quanti traits
                break;
        }
    }
    
    if(_GenerationParam && !_trait.empty() && !get_tot_occurrence()) return false;
    return true;
}

// ----------------------------------------------------------------------------------------
// get_filename
// ----------------------------------------------------------------------------------------
string& FileHandler::get_filename()
{
    _current_filename = _path + _filename + _popPtr->getReplicateCounter_()
    + _extension;
    return _current_filename;
}

// ----------------------------------------------------------------------------------------
// update
// ----------------------------------------------------------------------------------------
void FileHandler::update()
{
    if (!(_popPtr->getCurrentGeneration() % _GenerationOccurrence)
        || (_popPtr->getCurrentGeneration() == _GenerationOccurrence)) {
        FHwrite();
    }
}

// ----------------------------------------------------------------------------------------
// FHwrite
// ----------------------------------------------------------------------------------------
/** Writes genotype in a FSTAT-like format, the age class and the sex are added after the pop id.*/
void FileHandler::FHwrite()
{
    _popPtr->set_sampledInds(_age);      // update the sampling
    
    (this->*writeGenotype_func_ptr)(_format);
}

// ----------------------------------------------------------------------------------------
// write_Fstat
// ----------------------------------------------------------------------------------------
/** Writes genotype in a FSTAT-like format, the age class and the sex are added after the pop id.*/
void FileHandler::write_Fstat(bool extened)
{
    unsigned int t, l;
    
    // get the number of occupied patches
    if (!_popPtr->get_nbSamplePatch()) return; // if no patches are sampled
    
    // get the number of digits to visualize
    unsigned int max_allele = 0; // find the max number of alleles across all traits
    for (t = 0; t < _nb_trait; ++t) {
        for (l = 0; l < _trait[t]->get_nb_locus(); ++l) {
            if (max_allele < _trait[t]->get_nb_allele(l))
                max_allele = _trait[t]->get_nb_allele(l);
        }
    }
    
    // open the file
    string filename = get_path() + getGenerationReplicateFileName()
    + get_extension();
#ifdef _DEBUG
    message("FileHandler::FHwrite (%s)\n",filename.c_str());
#endif
    ofstream FILE(filename.c_str(), ios::out);
    if (!FILE)
        error("Could not open FSTAT output file '%s'!\n", filename.c_str());
    
    // write the heading line
    unsigned int position = getNbDigits(max_allele);
    unsigned int tot_nb_loci = 0;
    for (t = 0; t < _nb_trait; ++t) {
        tot_nb_loci += _trait[t]->get_nb_locus();
    }
    unsigned int maxID = (*_popPtr->get_vSamplePatch().rbegin())->get_ID()
    + 1;
    FILE << maxID << " " << tot_nb_loci << " " << max_allele << " " << position
    << "\n";
    
    // write the names of the loci
    string type = _trait[0]->_type;
    for (t = 0; t < _nb_trait; ++t) {
        for (l = 0; l < _trait[t]->get_nb_locus(); ++l) {
            FILE << type[0] << (t + 1) << "_l" << (l + 1) << "\n";
        }
    }
    
    // write all individuals of all patches
    unsigned int nbPatchDigits = getNbDigits(maxID);
    unsigned int a, mask;
    vector<TPatch*>::iterator curPop = _popPtr->get_vSamplePatch().begin();
    vector<TPatch*>::iterator endPop = _popPtr->get_vSamplePatch().end();
    for (; curPop != endPop; ++curPop) {
        for (a = 0, mask = 1; a < NB_AGE_CLASSES; ++a, mask <<= 1) {
            if (mask & _age) {
                if (_sex != 2)
                    write_Fstat(static_cast<age_idx>(a), FEM, FILE, *curPop,
                                nbPatchDigits, position, extened);
                if (_sex != 1)
                    write_Fstat(static_cast<age_idx>(a), MAL, FILE, *curPop,
                                nbPatchDigits, position, extened);
            }
        }
    }
    
    FILE.close();
    
    // call an external program and pass the name of this file
    if (!_script.empty()) execute_script(_script, filename);
}

// ----------------------------------------------------------------------------------------
// write_Fstat
// ----------------------------------------------------------------------------------------
void FileHandler::write_Fstat(const age_idx& AGE, const sex_t& SEX,
                              ostream& FILE, TPatch* curPop, const int& nbPatchDigit,
                              const int& position, bool extened)
{
    
    unsigned char** seq;
    unsigned int nb_locus, k, l, t;
    
    vector<TIndividual*>::iterator curInd =
    curPop->get_sampled_inds(SEX, AGE).begin();
    vector<TIndividual*>::iterator endInd =
    curPop->get_sampled_inds(SEX, AGE).end();
    for (; curInd != endInd; ++curInd) {
        FILE << setfill('0') << setw(nbPatchDigit) << curPop->get_ID() + 1
        << setfill(' ') << " ";
        for (t = 0; t < _nb_trait; ++t) { // for multiple instanciations of a trait
            nb_locus = _trait[t]->get_nb_locus();
            seq =
            (unsigned char**) (*curInd)->getTrait(_TTidx[t])->get_sequence();
            
            for (k = 0; k < nb_locus; ++k) {
                for (l = 0; l < ploidy; ++l) {
                    FILE.fill('0');
                    FILE.width(position);
                    FILE << (unsigned int) (seq[k][l] + 1);
                }
                FILE << " ";
            }
        }
        if (extened) write_individual_info_to_stream(FILE, *curInd, AGE, SEX, ' ');
        FILE << "\n";
    }
}

// ----------------------------------------------------------------------------------------
// FHwrite
// ----------------------------------------------------------------------------------------
/** writes the information of the individual to the stream */
void FileHandler::write_individual_info_to_stream(ostream& FILE,
                                                  TIndividual* ind, const age_idx& cur_age, const sex_t& cur_sex, char sep)
{
    FILE << cur_age + 1 << sep << cur_sex << sep << ind->getID() << sep
    << (ind->getMotherID() != "" ? ind->getMotherID() : my_NANstr )
    << sep
    << (ind->getFatherID() != "" ? ind->getFatherID() : my_NANstr );
    if (ind->getFitness() != my_NAN) FILE << sep << ind->getFitness();
    else FILE << sep << my_NANstr;
}

// ----------------------------------------------------------------------------------------
// write_Arlequin
// ----------------------------------------------------------------------------------------
/** Writes genotype in a Areqluin format, the age class and the sex are added after the pop id.*/
void FileHandler::write_Arlequin(bool extened)
{
    unsigned int t, l;
    
    // get the number of occupied patches
    if (!_popPtr->get_nbSamplePatch()) return; // if no patches are sampled
    
    // get the number of digits to visualize
    unsigned int max_allele = 0; // find the max number of alleles across all traits
    for (t = 0; t < _nb_trait; ++t) {
        for (l = 0; l < _trait[t]->get_nb_locus(); ++l) {
            if (max_allele < _trait[t]->get_nb_allele(l))
                max_allele = _trait[t]->get_nb_allele(l);
        }
    }
    
    // open the file
    string filename = get_path() + getGenerationReplicateFileName() + ".arp";
#ifdef _DEBUG
    message("FileHandler::ARPwrite (%s)\n",filename.c_str());
#endif
    ofstream FILE(filename.c_str(), ios::out);
    if (!FILE)
        error("Could not open Arlequin output file '%s'!\n", filename.c_str());
    
    char curTime[20];
    time_t tt = time(NULL);
    strftime(curTime, 20, "%d-%m-%Y %H:%M:%S", localtime(&tt));
    
    // write the heading
    FILE << "[Profile]"
    << "\n  Title=\"Genetic data simulated with quantiNemo2\""
    << "\n  # quantiNemo 2 v" << RELEASE << "." << REVISION << "."
    << MINOR_VERSION << TEMP_VERSION << "[" << VERSION_DATE << "; "
    << VERSION_TIME << "]" << "\n  # File created the " << curTime
    << "\n  NbSamples=" << _popPtr->get_nbSamplePatch()
    << "\n  GenotypicData=1" << "\n  GameticPhase=1"
    << "\n  RecessiveData=0" << "\n  DataType=STANDARD"
    << "\n  LocusSeparator=WHITESPACE" << "\n\n[Data]"
    << "\n  [[Samples]]";
    
    unsigned int a, mask, nbInd, nbIntDigit;
    unsigned int position = getNbDigits(max_allele);
    
    vector<TPatch*>::iterator curPop, endPop;
    for (curPop = _popPtr->get_vSamplePatch().begin(), endPop =
         _popPtr->get_vSamplePatch().end(); curPop != endPop;
         ++curPop) {
        // get the number of individuals
        for (nbInd = 0, a = 0, mask = 1; a < NB_AGE_CLASSES; ++a, mask <<= 1) {
            if (mask & _age) {
                if (_sex != 2)
                    nbInd += (*curPop)->sampleSize(FEM,
                                                   static_cast<age_idx>(a));
                if (_sex != 1)
                    nbInd += (*curPop)->sampleSize(MAL,
                                                   static_cast<age_idx>(a));
            }
        }
        if (!nbInd) continue;                 // stop here if the patch is empty
        
        FILE << "\n    SampleName=\"pop_" << (*curPop)->get_ID() + 1 << "\""
        << "\n    SampleSize=" << nbInd << "\n    SampleData={";
        
        // write each individual
        nbIntDigit = getNbDigits(nbInd);
        for (a = 0, mask = 1; a < NB_AGE_CLASSES; ++a, mask <<= 1) {
            if (mask & _age) {
                if (_sex != 2)
                    write_Arlequin(static_cast<age_idx>(a), FEM, FILE, *curPop,
                                   nbIntDigit, position, extened);
                if (_sex != 1)
                    write_Arlequin(static_cast<age_idx>(a), MAL, FILE, *curPop,
                                   nbIntDigit, position, extened);
            }
        }
        
        FILE << "\n    }\n";
    }
    FILE.close();
    
    // call an external program and pass the name of this file
    if (!_script.empty()) execute_script(_script, filename);
}

// ----------------------------------------------------------------------------------------
// write_Arlequin
// ----------------------------------------------------------------------------------------
void FileHandler::write_Arlequin(const age_idx& AGE, const sex_t& SEX,
                                 ostream& FILE, TPatch* curPop, const unsigned int& nbIndDigit,
                                 const int& position, bool extened)
{
    unsigned char** seq;
    unsigned int nb_locus, k, l, t;
    
    vector<TIndividual*>::iterator curInd =
    curPop->get_sampled_inds(SEX, AGE).begin();
    vector<TIndividual*>::iterator endInd =
    curPop->get_sampled_inds(SEX, AGE).end();
    for (unsigned int i = 1; curInd != endInd; ++curInd, ++i) {
        for (l = 0; l < ploidy; ++l) {
            if (!l) FILE << "\n      Ind_" << setfill('0') << setw(nbIndDigit)
                << i << " 1";
            else FILE << "\n          " << setfill(' ') << setw(nbIndDigit)
                << "" << "  ";
            for (t = 0; t < _nb_trait; ++t) { // for multiple instanciations of a trait
                nb_locus = _trait[t]->get_nb_locus();
                seq =
                (unsigned char**) (*curInd)->getTrait(_TTidx[t])->get_sequence();
                
                for (k = 0; k < nb_locus; ++k) {
                    FILE << " ";
                    FILE.fill('0');
                    FILE.width(position);
                    FILE << (unsigned int) (seq[k][l] + 1);
                }
            }
            if (!l && extened) {
                FILE << "  # ";
                write_individual_info_to_stream(FILE, *curInd, AGE, SEX, ' ');
            }
        }
    }
}

// ----------------------------------------------------------------------------------------
// write_Plink_ped
// ----------------------------------------------------------------------------------------
/** Writes genotype in a Plink format. Format
 - Family ID
 - Sample ID
 - Paternal ID
 - Maternal ID
 - Sex (1=male; 2=female; other=unknown)
 - Affection (0=unknown; 1=unaffected; 2=affected): individual fitness
 - Genotypes (space or tab separated, 2 for each marker. 0=missing)
*/
void FileHandler::write_Plink_ped(bool extended)
{
    char sep=' ';
    
    // get the number of occupied patches
    if (!_popPtr->get_nbSamplePatch()) return; // if no patches are sampled
    
    // open the file
    ofstream FILE;
    bool appended;
    unsigned int curGen = _popPtr->getCurrentGeneration();
    if(last_generation == curGen-1){    // concatenate output
        FILE.open(_filename_ped.c_str(), std::ios_base::app);
        if (!FILE) error("Could not open Plink (*.ped) output file '%s'!\n", _filename_ped.c_str());
        appended = true;
    }
    else{
        _filename_ped = get_path() + getGenerationReplicateFileName() + ".ped";
        FILE.open(_filename_ped.c_str(), ios::out);
        if (!FILE) error("Could not open Plink (*.ped) output file '%s'!\n", _filename_ped.c_str());
        appended = false;
        
        // write a comment heading line
        FILE << "# PLINK .ped file created by quantiNemo 2 v" << RELEASE << "."
        << REVISION << "." << MINOR_VERSION << TEMP_VERSION << "[" << VERSION_DATE << "; "
        << VERSION_TIME << "]" << "\n";

        // print time stamp
        char curTime[20];
        time_t tt = time(NULL);
        strftime(curTime, 20, "%d-%m-%Y %H:%M:%S", localtime(&tt));
        FILE << "# file created at " << curTime;
        if(extended) FILE << " (listing also multi-allele loci)";
        FILE  << "\n";
    }
    last_generation = curGen;
#ifdef _DEBUG
    message("FileHandler::FHwrite (%s)\n",_filename_ped.c_str());
#endif
    
    // write all individuals of all patches
    unsigned int maxID = (*_popPtr->get_vSamplePatch().rbegin())->get_ID() + 1;
    unsigned int nbPatchDigits = getNbDigits(maxID);
    unsigned int a, mask;
    vector<TPatch*>::iterator curPop = _popPtr->get_vSamplePatch().begin();
    vector<TPatch*>::iterator endPop = _popPtr->get_vSamplePatch().end();
    for (; curPop != endPop; ++curPop) {
        for (a = 0, mask = 1; a < NB_AGE_CLASSES; ++a, mask <<= 1) {
            if (mask & _age) {
                if (_sex != 2)
                    write_Plink_ped(static_cast<age_idx>(a), FEM, FILE, *curPop,
                                nbPatchDigits, 1, sep, extended, appended);
                if (_sex != 1)
                    write_Plink_ped(static_cast<age_idx>(a), MAL, FILE, *curPop,
                                nbPatchDigits, 1, sep, extended, appended);
            }
        }
    }
    
    FILE.close();
    
    // call an external program and pass the name of this file
    if (!_script.empty()) execute_script(_script, _filename_ped);
    
    // write all phenotypes to a speparate file
    write_Plink_pheno(extended, appended);
}

// ----------------------------------------------------------------------------------------
// write_Plink_ped
// ----------------------------------------------------------------------------------------
/** format:
 - Family ID
 - Sample ID
 - Paternal ID
 - Maternal ID
 - Sex (1=male; 2=female; other=unknown)
 - Affection (0=unknown; 1=unaffected; 2=affected): individual fitness
 - Genotypes (space or tab separated, 2 for each marker. 0=missing)
 */
void FileHandler::write_Plink_ped(const age_idx& AGE, const sex_t& SEX,
                              ostream& FILE, TPatch* curPop, const int& nbPatchDigit,
                              const int& position, char sep, bool extended, bool append)
{
    
    unsigned char** seq;
    unsigned int nb_locus, k, l, t;
    
    vector<TIndividual*>::iterator curInd = curPop->get_sampled_inds(SEX, AGE).begin();
    vector<TIndividual*>::iterator endInd = curPop->get_sampled_inds(SEX, AGE).end();
    for (; curInd != endInd; ++curInd) {
        FILE << setfill('0') << setw(nbPatchDigit) << curPop->get_ID() + 1 << setfill(' ') << sep;
        FILE << (*curInd)->getID() << sep;
		FILE << ((append && (*curInd)->getFather()) ? (*curInd)->getFatherID() : (string)"0" ) << sep;
		FILE << ((append && (*curInd)->getMother()) ? (*curInd)->getMotherID() : (string)"0" ) << sep;
		FILE << (*curInd)->getSex() + 1 << sep;
        FILE << ((*curInd)->getFitness() != my_NAN ? (*curInd)->getFitness() : -9 );

        for (t = 0; t < _nb_trait; ++t) { // for multiple instanciations of a trait
            nb_locus = _trait[t]->get_nb_locus();
            seq = (unsigned char**) (*curInd)->getTrait(_TTidx[t])->get_sequence();
            
            for (k = 0; k < nb_locus; ++k) {
                if(_trait[t]->get_nb_allele(k)!=2 && !extended) continue; // not a SNP
                FILE << sep;
                for (l = 0; l < ploidy; ++l) {
                    FILE << sep;
                    FILE.fill('0');
                    FILE.width(position);
                    FILE << (unsigned int) (seq[k][l] + 1);
                }
            }
        }
        FILE << "\n";
    }
}

// ----------------------------------------------------------------------------------------
// write_Plink_pheno
// ----------------------------------------------------------------------------------------
/** Writes alternative phenotypeifle (Plink). Fomrat:
 - Family ID
 - Sample ID
 - Phenotype 1
 - Phenotype 2
 ...
 */
void FileHandler::write_Plink_pheno(bool extended, bool appended)
{
    char sep=' ';
    
    // get the number of occupied patches
    if (!_popPtr->get_nbSamplePatch()) return; // if no patches are sampled
    if (_quantiTraitIndexes.empty()) return;   // no quantitative traits are simulated
    
    // open the file
    ofstream FILE;
    if(appended){    // concatenate output
        FILE.open(_filename_pheno.c_str(), std::ios_base::app);
        if (!FILE) error("Could not open Plink (*.pheno) output file '%s'!\n", _filename_pheno.c_str());
    }
    else{
        _filename_pheno = get_path() + getGenerationReplicateFileName() + ".pheno";
        FILE.open(_filename_pheno.c_str(), ios::out);
        if (!FILE) error("Could not open Plink (*.pheno) output file '%s'!\n", _filename_pheno.c_str());
        
        // write a comment heading line
        FILE << "# PLINK .pheno file created by quantiNemo 2 v" << RELEASE << "."
        << REVISION << "." << MINOR_VERSION << TEMP_VERSION << "[" << VERSION_DATE << "; "
        << VERSION_TIME << "]" << "\n";
        
        // print time stamp
        char curTime[20];
        time_t tt = time(NULL);
        strftime(curTime, 20, "%d-%m-%Y %H:%M:%S", localtime(&tt));
        FILE << "# file created at " << curTime << "\n";
    }
#ifdef _DEBUG
    message("FileHandler::FHwrite (%s)\n",_filename_pheno.c_str());
#endif
    
    // write all individuals of all patches
    unsigned int maxID = (*_popPtr->get_vSamplePatch().rbegin())->get_ID() + 1;
    unsigned int nbPatchDigits = getNbDigits(maxID);
    unsigned int a, mask;
    vector<TPatch*>::iterator curPop = _popPtr->get_vSamplePatch().begin();
    vector<TPatch*>::iterator endPop = _popPtr->get_vSamplePatch().end();
    for (; curPop != endPop; ++curPop) {
        for (a = 0, mask = 1; a < NB_AGE_CLASSES; ++a, mask <<= 1) {
            if (mask & _age) {
                if (_sex != 2)
                    write_Plink_pheno(static_cast<age_idx>(a), FEM, FILE, *curPop,
                                    nbPatchDigits, 1, sep);
                if (_sex != 1)
                    write_Plink_pheno(static_cast<age_idx>(a), MAL, FILE, *curPop,
                                    nbPatchDigits, 1, sep);
            }
        }
    }
    
    FILE.close();
    
    // call an external program and pass the name of this file
    if (!_script.empty()) execute_script(_script, _filename_pheno);
}

// ----------------------------------------------------------------------------------------
// write_Plink_pheno
// ----------------------------------------------------------------------------------------
/** format:
 - Family ID
 - Sample ID
 - Phenotype 1
 - Phenotype 2
 ...
 */
void FileHandler::write_Plink_pheno(const age_idx& AGE, const sex_t& SEX,
                                  ostream& FILE, TPatch* curPop, const int& nbPatchDigit,
                                  const int& position, char sep)
{
    double pheno;
    vector<TIndividual*>::iterator curInd = curPop->get_sampled_inds(SEX, AGE).begin();
    vector<TIndividual*>::iterator endInd = curPop->get_sampled_inds(SEX, AGE).end();
    vector<int>::iterator curTrait, endTrait=_quantiTraitIndexes.end();
    for (; curInd != endInd; ++curInd) {
        FILE << setfill('0') << setw(nbPatchDigit) << curPop->get_ID() + 1 << setfill(' ') << sep;
        FILE << (*curInd)->getID() << sep;
        
        
        for(curTrait=_quantiTraitIndexes.begin(); curTrait!=endTrait; ++curTrait){
            pheno = (*curInd)->getTraitPhenotype(*curTrait);  // first qTrait
            FILE << sep;
            FILE << (pheno == my_NAN ? -9 : pheno);
        }
        FILE << "\n";
    }
}

// ----------------------------------------------------------------------------------------
// write_Plink
// ----------------------------------------------------------------------------------------
/** Writes genotype in a Plink format. Format
 - chromosome (1-22, X, Y or 0 if unplaced)
 - rs# or snp identifier
 - Genetic distance (morgans)
 - Base-pair position (bp units)
 */
void FileHandler::write_Plink_map(bool extended)
{
    // get the number of occupied patches
    if (!_popPtr->get_nbSamplePatch()) return; // if no patches are sampled
    
    // open the file
    string filename = get_path() + getReplicateFileName();
    
    
    if(_popPtr->get_protoGenome()->is_sexSpecificMap()){
        write_Plink_map(filename + "_fem.map", FEM, extended);
        write_Plink_map(filename + "_mal.map", MAL, extended);
    }
    else{
        write_Plink_map(filename + ".map", FEM, extended);
    }
}

// ----------------------------------------------------------------------------------------
// write_Plink
// ----------------------------------------------------------------------------------------
/** Writes genotype in a Plink format. Format
 - chromosome (1-22, X, Y or 0 if unplaced)
 - rs# or snp identifier
 - Genetic distance (morgans)
 - Base-pair position (bp units)
 */
void FileHandler::write_Plink_map(string filename, sex_t SEX, bool extended)
{
    // get the number of occupied patches
    assert(_popPtr->get_nbSamplePatch());
    
#ifdef _DEBUG
    message("FileHandler::FHwrite (%s)\n",filename.c_str());
#endif
    ofstream FILE(filename.c_str(), ios::out);
    if (!FILE)
        error("Could not open Plink (*.map) output file '%s'!\n", filename.c_str());
    
    // write a comment heading line
    FILE<< "# PLINK .map file created by quantiNemo 2 v" << RELEASE << "."
    << REVISION << "." << MINOR_VERSION << TEMP_VERSION << "[" << VERSION_DATE << "; "
    << VERSION_TIME << "]\n";
    
    // print time stamp
    char curTime[20];
    time_t tt = time(NULL);
    strftime(curTime, 20, "%d-%m-%Y %H:%M:%S", localtime(&tt));
    FILE << "# file created at " << curTime;
    if(extended) FILE << " (listing also multi-allele loci)";
    
    unsigned int t, l, nb_locus;
    char sep= ' ';
    TLocus* curLocus;
    
    // get the max conromosome index (to be able to correctly label the unlinked loci
    TGenomeProto* pGenome = _popPtr->get_protoGenome();
    unsigned int unlinkedChrom = 0;
    if(pGenome->get_nb_locus_linked()) unlinkedChrom = pGenome->get_locus_tot(pGenome->get_nb_locus_linked()-1).get_chromosomePosition()+1;
    
    for (t = 0; t < _nb_trait; ++t) { // for multiple instanciations of a trait
        nb_locus = _trait[t]->get_nb_locus();
        
        for (l = 0; l < nb_locus; ++l) {
            curLocus = &_trait[t]->get_aLocus()[l];
            if(curLocus->get_nb_allele() != 2 && !extended) continue; // not a SNP
            
            if(curLocus->get_chromosomePosition()==my_NAN){ // unlinked locus
                FILE << "\n" << ++unlinkedChrom
                << sep << curLocus->get_locus_id_tot()+1
                << sep << 0
                << sep << 0;
            }
            else{                                           // linked locus
                FILE << "\n" << curLocus->get_chromosomePosition()+1
                << sep << curLocus->get_locus_id_tot()+1
                << sep << curLocus->get_locusPosition(SEX)/100
                << sep << 0;
            }
        }
    }
    FILE << "\n";
    
    FILE.close();
    
    // call an external program and pass the name of this file
    if (!_script.empty()) execute_script(_script, filename);
}

// ----------------------------------------------------------------------------------------
/** launch the script at the end of all simulations */
void FileHandler::execute_script(string scriptname, string filename)
{
    runScript(_popPtr->get_pSimulation()->get_fileName(scriptname),
              _popPtr->get_pSimulation()->get_fileName(filename));

}

// ----------------------------------------------------------------------------------------
// get_tot_occurrence
// ----------------------------------------------------------------------------------------
/** returns the number of times the stats are really computed (taking into account temporal parameters) */
unsigned int FileHandler::get_tot_occurrence()
{
    return (unsigned int)get_logtime_occurences(_GenerationParam, _popPtr->getGenerations(),
                                  _popPtr->isCoalescence()).size();
}

// ----------------------------------------------------------------------------------------
// getReplicateFileName
// ----------------------------------------------------------------------------------------
string&
FileHandler::getReplicateFileName()
{
    _rep_filename = _filename + _popPtr->getReplicateCounter_r();
    
    return _rep_filename;
}

// ----------------------------------------------------------------------------------------
// getGenerationReplicateFileName
// ----------------------------------------------------------------------------------------
string FileHandler::getGenerationReplicateFileName()
{
    string name = _filename + _popPtr->getGenerationCounter_g()
    + _popPtr->getReplicateCounter_r();
    return name;
}

// ----------------------------------------------------------------------------------------
// set_filename
// ----------------------------------------------------------------------------------------
void FileHandler::set_filename(string s)
{
    if(s.empty() && _service) _filename = _popPtr->getBaseFileName(); // not specified: use the basename
    else _filename = s;
}

// ----------------------------------------------------------------------------------------
// set_filename
// ----------------------------------------------------------------------------------------
void FileHandler::set_path (string& path)
{
    _path = path;
    if(!_path.empty() && _path[_path.length()-1] != SEP) _path += SEP;
    _short_path = _path;
    if(!_popPtr->getSimfolder().empty()){
        _path = _popPtr->getSimfolder()+_path;
    }
    check_path(_path);
}
