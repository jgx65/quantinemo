/** @file geneticmap.cpp
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
//---------------------------------------------------------------------------
#include "tgenome.h"
#include "tmetapop.h"


unsigned int TGenomeProto::MUTATION_TRIAL;
//---------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------
// TGenome
// ----------------------------------------------------------------------------------------
TGenome::TGenome()
{
	sequence = NULL;
	_protoGenome = NULL;
}

// ----------------------------------------------------------------------------------------
// ~TGenome
// ----------------------------------------------------------------------------------------
TGenome::~TGenome()
{
	clear();
}

// ----------------------------------------------------------------------------------------
// clear
// ----------------------------------------------------------------------------------------
void TGenome::clear()
{
	if (sequence) {
		assert(_protoGenome);
		ARRAY::delete_2D(sequence, _protoGenome->_nb_locus_tot);
	}
}

// ----------------------------------------------------------------------------------------
// inherit
// ----------------------------------------------------------------------------------------
void TGenome::inherit(TIndividual* mother, TIndividual* father)
{
	_protoGenome->inherit(mother, father, sequence);
}

// ----------------------------------------------------------------------------------------
// mutate
// ----------------------------------------------------------------------------------------
void TGenome::mutate()
{
	_protoGenome->mutate(sequence);
}

// ----------------------------------------------------------------------------------------
// ini_sequence
// ----------------------------------------------------------------------------------------
void TGenome::ini_sequence(TPatch* patch)
{
	_protoGenome->ini_sequence(sequence, patch);
}

// ----------------------------------------------------------------------------------------
// ini_sequence
// ----------------------------------------------------------------------------------------
/** initialize the sequence with a given genotype */
void TGenome::ini_sequence(unsigned char* seq_mum, unsigned char* seq_dad)
{
	for (unsigned int l = 0; l < _protoGenome->_nb_locus_tot; ++l) {
		sequence[l][0] = seq_mum[l];
		sequence[l][1] = seq_dad[l];
	}
}

// ----------------------------------------------------------------------------------------
// ini_sequence
// ----------------------------------------------------------------------------------------
/** initialize the sequence with a given genotype
 * the input is a double sequence: first mum, than dad
 */
void TGenome::ini_sequence(unsigned char** seq)
{
	for (unsigned int l = 0; l < _protoGenome->_nb_locus_tot; ++l) {
		sequence[l][0] = seq[0][l];
		sequence[l][1] = seq[1][l];
	}
}

// ----------------------------------------------------------------------------------------
// ini_sequence_dadFirst
// ----------------------------------------------------------------------------------------
/** initialize the sequence with a given genotype
 * the input is a double sequence: first dad, than mum
 */
void TGenome::ini_sequence_dadFirst(unsigned char** seq)
{
	for (unsigned int l = 0; l < _protoGenome->_nb_locus_tot; ++l) {
		sequence[l][0] = seq[1][l];
		sequence[l][1] = seq[0][l];
	}
}

// ----------------------------------------------------------------------------------------
// create_sequence
// ----------------------------------------------------------------------------------------
void TGenome::create_sequence()
{
	assert(!sequence);
	ARRAY::create_2D(sequence, _protoGenome->_nb_locus_tot, ploidy);
}

// ----------------------------------------------------------------------------------------
// clone
// ----------------------------------------------------------------------------------------
void TGenome::clone(const TGenome& gen)
{
	_protoGenome = gen._protoGenome;
	assert(!sequence);
	ARRAY::create_2D(sequence, _protoGenome->_nb_locus_tot, ploidy);
}

// ----------------------------------------------------------------------------------------
// copy
// ----------------------------------------------------------------------------------------
TGenome&
TGenome::operator=(const TGenome& g)
{
	_protoGenome = g._protoGenome;
	unsigned int l, a, nbLocus = _protoGenome->_nb_locus_tot;
	if (!sequence) ARRAY::create_2D(sequence, nbLocus, ploidy);
	for (l = 0; l < nbLocus; ++l) {
		for (a = 0; a < ploidy; ++a) {
			sequence[l][a] = g.sequence[l][a];
		}
	}
	return *this;
}

// ----------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------
// constructor
// ----------------------------------------------------------------------------------------
TGenomeProto::TGenomeProto()
{
#ifdef _DEBUG
	message(" TGenomeProto::TGenomeProto\n");
#endif
	reset_locus_vector();
	ini();
	ini_paramset();
}

// ----------------------------------------------------------------------------------------
// ini
// ----------------------------------------------------------------------------------------
void TGenomeProto::ini()
{
    
	_locus_tot = NULL;
	_nb_locus_linked = 0;
	_nb_locus_unlinked = 0;
	_mut_rate_mean = 0;
	_change_mutationRates = true;
	_nb_chromosome = 0;
	_nb_locus_per_chromosome = NULL;
    
	for (unsigned int s = 0; s < 2; ++s) {		// for sex specific settings
		_recombination_factor_chrom[s] = NULL;
		_recombination_chrom_func_ptr[s] = NULL;
		_chromosomeSize[s] = NULL;
		_locus_position_tot[s] = NULL;
		_locus_position_tot_temp[s] = NULL;
		_chromosomeLength[s] = 0;
        
	}
    
	_locus_pleiotropic = NULL;
	_nb_locus_pleiotropic = 0;
	_mutate_func_ptr = &TGenomeProto::_mutate_zero_mutation_rate;
}

// ----------------------------------------------------------------------------------------
// ini
// ----------------------------------------------------------------------------------------
void TGenomeProto::ini_paramset()
{
	set_paramset("genome", "genetic map", false, _popPtr);
    
	// genetic map output (used for debugging)
	add_parameter("_genetic_map_output", INT2, false, 0, 1, "0", false,
                  "Output the genetic map for all traits:\n" \
                  "  0: no\n" \
                  "  1: yes",5);
    
	// factor to change the recombination rate
    add_parameter("recombination_factor", MAT, false, 0, my_NAN, "1", false,
                  "Allows modifying the recombination rate of the genetic map, " \
                  "by multiplying the genetic distance with the factor.",0);
    
    add_parameter("recombination_factor_fem", MAT, false, 0, my_NAN, "1", false,
                  "Allows modifying the recombination rate of the female genetic map, " \
                  "by multiplying the genetic distance with the factor.",0);
    
    add_parameter("recombination_factor_mal", MAT, false, 0, my_NAN, "1", false,
                  "Allows modifying the recombination rate of the male genetic map, " \
                  "by multiplying the genetic distance with the factor.",0);
    
    
	// genome positions
    add_parameter("genome", MAT, false, my_NAN, my_NAN, "", false,
                  "The genetic map in cM.",5);
    
    add_parameter("genome_fem", MAT, false, my_NAN, my_NAN, "", false,
                  "The genetic map for females in cM.",5);
    
    add_parameter("genome_mal", MAT, false, my_NAN, my_NAN, "", false,
                  "The genetic map for males in cM.",5);
    

    add_parameter("mutation_trial", INT2, false, 1, my_NAN, "1e4", false,
                  "The number of trials to perform to draw a new allele at a mutation event " \
                  " before resigning.",5);
    
}

// ----------------------------------------------------------------------------------------
// destructor
// ----------------------------------------------------------------------------------------
TGenomeProto::~TGenomeProto()
{
#ifdef _DEBUG
	message(" TGenomeProto::~TGenomeProto\n");
#endif
    
	// delete each locus position
	map<string, vector<TLocusPos*> >::iterator curType, endType;
	vector<TLocusPos*>::iterator curPos, endPos;
	
    if(_locus_pos.find("global")==_locus_pos.end()){    // per trait defined map
        for (curType = _locus_pos.begin(), endType = _locus_pos.end();
             curType != endType; ++curType) {
            for (curPos = curType->second.begin(), endPos = curType->second.end();
                 curPos != endPos; ++curPos) {
                delete *curPos;
            }
        }
    }
    else{                                               // globally defined map
        for (curPos = _locus_pos["global"].begin(), endPos = _locus_pos["global"].end();
             curPos != endPos; ++curPos) {
            delete *curPos;
        }
    }
    
	clear();
}

//------------------------------------------------------------------------------
// print_genome
//------------------------------------------------------------------------------
/** outputs the genetic map for all the traits */
void TGenomeProto::print_genome()
{
	if (!get_parameter_value("_genetic_map_output") || !_nb_chromosome) return;
    
	// create the file name
	string file = _popPtr->getSimfolder() + "genetic_map"
    + _popPtr->getReplicateCounter_r() + ".txt";
    
#ifdef _DEBUG
	message("TGenomeProto::print_genome (%s)\n",file.c_str());
#endif
    
	// open the file
	ofstream FILE(file.c_str());
	if (!FILE) error("opening file '%s' to write the genetic map!\n", file.c_str());
    
	// print the heading
	FILE << "#Genetic map";
	if (_popPtr->getReplicates() > 1) FILE << " (Replicate: " << _popPtr->getCurrentReplicate()+1 << ")";
	FILE << "\n#" << setfill('=') << setw(30) << "" << setfill(' ');

    // first way to specify it
    FILE << "\n#Globally defined genetic map:";
    print_genome_tot(FILE);

    map<string, TTraitProto*>::iterator cur=_popPtr->getTraitPrototypes().begin();
    map<string, TTraitProto*>::iterator end=_popPtr->getTraitPrototypes().end();
	for (; cur != end; ++cur) {
        print_locus_index_tot(FILE, cur->second);
    }

    FILE << "\n";
    
    //second way to specify it
    FILE << "\n#Per trait defined genetic map:";
	// print the genetic map of each trait on a separate line
    map<string, vector<TLocusPos*> >::iterator curTrait = _locus_pos.begin();
    map<string, vector<TLocusPos*> >::iterator endTrait = _locus_pos.end();
	for (curTrait = _locus_pos.begin(), endTrait = _locus_pos.end(); curTrait != endTrait; ++curTrait) {
		print_genome_trait(FILE, curTrait->first);
	}

    cur=_popPtr->getTraitPrototypes().begin();
    end=_popPtr->getTraitPrototypes().end();
	for (; cur != end; ++cur) {
        print_locus_index_trait(FILE, cur->second);
    }
    

    FILE.close();
}

// ----------------------------------------------------------------------------------------
// print_genome
// ----------------------------------------------------------------------------------------
/** print the genome of the trait to the stream */
void TGenomeProto::print_genome_tot(ofstream& FILE)
{
	if (!is_sexSpecificMap()) { // a single genetic map
		string map = print_genome_tot_sex(FEM);
		if (map.size() != 4) FILE << "\n" << setw(25) << left << "genome " << map;
	}
	else {                                           // sex specific genetic map
		string map = print_genome_tot_sex(FEM);
		if (map.size() != 4) FILE << "\n" << setw(25) << left << "genome_fem " << map; // female genetic map
		map = print_genome_tot_sex(MAL);
		if (map.size() != 4) FILE << "\n" << setw(25) << left << "genome_mal " << map; // male genetic map
	}
}

// ----------------------------------------------------------------------------------------
// print_genome
// ----------------------------------------------------------------------------------------
/** print the genome of the trait to the stream */
string TGenomeProto::print_genome_tot_sex(sex_t SEX)
{
	TLocus *curLocus;
	TLocusPos* curLocPos = NULL;
	unsigned int curChrom, lastChrom = my_NAN;
	ostringstream out;
	out << "{{";
	for (unsigned int l = 0; l < _nb_locus_linked; ++l) { // for each linked locus
		curLocus = _locus_tot[l];                       // get the current locus
		if (curLocPos == curLocus->get_locusPosition()) continue; // pointer comparison: this position was already printed
		curLocPos = curLocus->get_locusPosition();		// get the new pointer
		curChrom = curLocus->get_chromosomePosition();
		assert(curChrom != my_NAN);
		if (curChrom != lastChrom) {                // start of a new chromosome
			if (lastChrom == my_NAN) { // if it is the first chromosome in the genome
				if (curChrom) out << curChrom + 1 << ": "; // if one or several chromosomes are jumped
			}
			else {                                        // not fist chromosome
				out << "}{";
				if (curChrom > lastChrom + 1) out << curChrom + 1 << ": "; // if one or several chromosomes are jumped
			}
			lastChrom = curChrom;                 // update the chromosome index
		}
		else out << " ";                   // same chromosome => we need a space
		out << curLocus->get_locusPosition(SEX);           // print the position
	}
	out << "}}";
	return out.str();
}


// ----------------------------------------------------------------------------------------
// print_genome
// ----------------------------------------------------------------------------------------
/** print the genome of the trait to the stream */
void TGenomeProto::print_genome_trait(ofstream& FILE, const string& type)
{
	if (_locus_position_tot[FEM] == _locus_position_tot[MAL]) { // a single genetic map
		string map = print_genome_trait_sex(FEM, type);
		if (map.size() != 4) FILE << "\n" << setw(25) << left << type + "_genome" << map;
	}
	else {                                           // sex specific genetic map
		string map = print_genome_trait_sex(FEM, type);
		if (map.size() != 4) FILE << "\n" << setw(25) << left << type + "_genome_fem" << map; // female genetic map
		map = print_genome_trait_sex(MAL, type);
		if (map.size() != 4) FILE << "\n" << setw(25) << left << type + "_genome_mal" << map; // male genetic map
	}
}

// ----------------------------------------------------------------------------------------
/** returns a string of the genetic map of the trait. */
string TGenomeProto::print_genome_trait_sex(sex_t SEX, const string& type)
{
	TLocus *curLocus;
	TLocusPos* curLocPos = NULL;
	unsigned int curChrom, lastChrom = my_NAN;
	ostringstream out;
	out << "{{";
	for (unsigned int l = 0; l < _nb_locus_linked; ++l) { // for each linked locus
		if (_locus_tot[l]->get_pTrait()->get_type() != type) continue; // not the correct type
		curLocus = _locus_tot[l];                       // get the current locus
		if (curLocPos == curLocus->get_locusPosition()) continue; // pointer comparison: this position was already printed
		curLocPos = curLocus->get_locusPosition();		// get the new pointer
		curChrom = curLocus->get_chromosomePosition();
		assert(curChrom != my_NAN);
		if (curChrom != lastChrom) {                // start of a new chromosome
			if (lastChrom == my_NAN) { // if it is the first chromosome in the genome
				if (curChrom) out << curChrom + 1 << ": "; // if one or several chromosomes are jumped
			}
			else {                                        // not fist chromosome
				out << "}{";
				if (curChrom > lastChrom + 1) out << curChrom + 1 << ": "; // if one or several chromosomes are jumped
			}
			lastChrom = curChrom;                 // update the chromosome index
		}
		else out << " ";                   // same chromosome => we need a space
		out << curLocus->get_locusPosition(SEX);           // print the position
	}
	out << "}}";
	return out.str();
}

// ----------------------------------------------------------------------------------------
/** returns a string of the genetic map of the trait. */
void TGenomeProto::print_locus_index_trait(ofstream& out, TTraitProto* pTrait)
{
    // generate a temp map for the conversion of totID to traitID
    string type = pTrait->get_type();
    TLocusPos* curPos = _locus_tot[0]->get_locusPosition();
    bool firstI = true;
    map<unsigned int, unsigned int> tot2trait;      //tot2trait[totID]=traitID
    for(unsigned int traitID=0, l=0; l<_nb_locus_linked; ++l){
        if(_locus_tot[l]->get_pTrait()->get_type() != type) continue;
        if(curPos==_locus_tot[l]->get_locusPosition() || firstI){
            firstI = false;
            tot2trait[l] = traitID;
        }
        else{
            curPos=_locus_tot[l]->get_locusPosition();
            tot2trait[l] = ++traitID;
        }
    }

    TLocus* aLocus = pTrait->get_aLocus();
    unsigned int nbLocus = pTrait->get_nb_locus();
    bool first = true;
	out << "\n" << setw(25) << left << pTrait->get_type_index() + "_locus_index" << "{";
    for(unsigned int i=0; i<nbLocus; ++i){
		if (aLocus[i].get_chromosomePosition()==my_NAN) continue; // not linked
		if(!first) out << " ";
        
        unsigned int pos = aLocus[i].get_locus_id_tot();
        assert(tot2trait.find(pos)!=tot2trait.end());
        out << tot2trait[pos]+1;
        first = false;
    }
	out << "}";
}

// ----------------------------------------------------------------------------------------
/** returns a string of the genetic map of the trait. */
void TGenomeProto::print_locus_index_tot(ofstream& out, TTraitProto* pTrait)
{
    TLocus* aLocus = pTrait->get_aLocus();
    unsigned int nbLocus = pTrait->get_nb_locus();
    bool first = true;
	out << "\n" << setw(25) << left << pTrait->get_type_index() + "_locus_index" << "{";
    for(unsigned int i=0; i<nbLocus; ++i){
		if (aLocus[i].get_chromosomePosition()==my_NAN) continue; // not linked
		if(!first) out << " ";
        
        unsigned int pos = aLocus[i].get_locus_id_tot();
        if(_locus_pleiotropic) pos = ARRAY::search_smallerEqual(_locus_pleiotropic, _nb_locus_pleiotropic, pos);
        
        out << pos+1;
        first = false;
    }
	out << "}";
}

// ----------------------------------------------------------------------------------------
// set_genome_positions
// ----------------------------------------------------------------------------------------
/** this function sets the entire genome of each type of markers.
 */
void TGenomeProto::set_genome_positions(string type, TTraitProto* pTrait)
{
    set_genome_positions_fix();                         // general
    if(!set_genome_positions_fix(type, pTrait)){        // type specific
        set_genome_positions_none(type);                // no map
    }
}

// ----------------------------------------------------------------------------------------
// set_genome_positions_fix
// ----------------------------------------------------------------------------------------
/** set the loci positions explicitly.
 * returns true if done and false if not
 */
bool TGenomeProto::set_genome_positions_fix()
{
    string type ="global";
    if(_locus_pos.find(type)!= _locus_pos.end()) return true; // chekc if already done
    
    // check if there are sex specific settings
    Param *pParamF = get_parameter("genome_fem");
    Param *pParamM = get_parameter("genome_mal");
    Param *pParam  = get_parameter("genome");

	if (pParamF->isSet()){
        assert(_locus_pos.find(type)== _locus_pos.end());
        set_genome_positions_matrix(pParamF->get_matrixVar(), FEM, type);     // female
        
		if (_popPtr->get_sexInitRatio()) {
			if (!pParamM->isSet()) error("Only one sex specific genome is specified: both are required!\n");
			set_genome_positions_matrix(pParamM->get_matrixVar(), MAL, type);   // male
		}
        else set_genome_positions_matrix(pParamF->get_matrixVar(), FEM, type);     // female
		return true;
	}
    if (pParamM->isSet()) error("Only one sex specific genome is specified: both are required!\n");
    
	// check if there is a general setting
	if (pParam->isSet()) {
        assert(_locus_pos.find(type)== _locus_pos.end());
        set_genome_positions_matrix(pParam->get_matrixVar(), FEM, type);
		return true;
	}
    
	return false; // not set
}

// ----------------------------------------------------------------------------------------
// set_genome_positions_fix
// ----------------------------------------------------------------------------------------
/** set the loci positions explicitly.
 * returns true if done and false if not
 */
bool TGenomeProto::set_genome_positions_fix(string type, TTraitProto* pTrait)
{
    // check if there are sex specific settings
    Param *pParamF = pTrait->get_parameter(type + "_genome_fem");
    Param *pParamM = pTrait->get_parameter(type + "_genome_mal");
    Param *pParam  = pTrait->get_parameter(type + "_genome");
    
	if (pParamF->isSet()){
        assert(_locus_pos.find(type)== _locus_pos.end());
        
        // check if it is globally defined
        if(_locus_pos.find("global")!=_locus_pos.end())
            error("Genetic map may be either defined globally or per trait, but not in both ways!\n");
        
        set_genome_positions_matrix(pParamF->get_matrixVar(), FEM, type);     // female
        
		if (_popPtr->get_sexInitRatio()) {
			if (!pParamM->isSet()) error("Only one sex specific genome is specified: both are required!\n");
			set_genome_positions_matrix(pParamM->get_matrixVar(), MAL, type);   // male
		}
        else set_genome_positions_matrix(pParamF->get_matrixVar(), FEM, type);     // female
		return true;
	}
    if (pParamM->isSet()) error("Only one sex specific genome is specified: both are required!\n");
    
	// check if there is a general setting
	if (pParam->isSet()) {
        assert(_locus_pos.find(type)== _locus_pos.end());
        
        // check if it is globally defined
        if(_locus_pos.find("global")!=_locus_pos.end())
            error("Genetic map may be either defined globally or per trait, but not in both ways!\n");
        
        set_genome_positions_matrix(pParam->get_matrixVar(), FEM, type);
		return true;
	}
    
    
    // overtake the global definition if present
    if(_locus_pos.find("global")!=_locus_pos.end()){
        _locus_pos[type]=_locus_pos["global"];
        return true;
    }
    
	return false; // not set
}

// ----------------------------------------------------------------------------------------
// set_genome_positions_none
// ----------------------------------------------------------------------------------------
/** if no genetic map is present a virtual locus position has to be set with position NaN
 */
bool TGenomeProto::set_genome_positions_none(string type)      // no genetic map
{
	add_locus_position(type, FEM, my_NAN, my_NAN, 0); // and a Nan position
	return true;
}

// ----------------------------------------------------------------------------------------
// set_genome_positions
// ----------------------------------------------------------------------------------------
/** this function sets the genetic map depending on the parameters.
 * if delMatrix is true (default) the matrix is deleted at the end
 */
void TGenomeProto::set_genome_positions_matrix(TMatrixVar<double>* matrix, sex_t SEX,
                                               string type, bool delMatrix)
{
	unsigned int c, l, p, size;
	unsigned int nbChromosome = matrix->getNbRows();
    
	// for each chromosome
	for (c = 0, l = 0; c < nbChromosome; ++c) {
		size = matrix->getNbCols(c);
		for (p = 0; p < size; ++p, ++l) {
			add_locus_position(type, SEX, c, matrix->get(c, p), l);
		}
	}
    
	// add an unlinked position
	add_locus_position(type, SEX, my_NAN, my_NAN, l);
    
	if (delMatrix) delete matrix;
}

// ----------------------------------------------------------------------------------------
// ini_genetic_map
//------------------------------------------------------------------------------
/** The super chromosome is created based on the positions of the loci in the vector _locus_vector
 */
void TGenomeProto::ini_genetic_map()
{
    
	assert(!_locus_tot);
    
#ifdef _DEBUG
	message("TGenomeProto::ini_genetic_map\n");
#endif
    
	// reset the genetic map
	clear();
	ini();
    
	create_locus_tot();      // _locus_vector => _locus_tot
	if (_nb_locus_linked) create_locus_position_tot(); // _locus_tot    => _locus_position_tot
	create_pleiotrophy(); // _locus_tot    => _locus_pleiotropic
	
    // set the inheritance function
    if(_nb_locus_linked){
        if (_nb_locus_unlinked) _inherit_func_ptr = &TGenomeProto::_inherit_mixed;
        else _inherit_func_ptr = &TGenomeProto::_inherit_linked;
    }
    else _inherit_func_ptr = &TGenomeProto::_inherit_unlinked;

	// dump the genome to file if requested
	print_genome();
}

//------------------------------------------------------------------------------
/** Check if pleiotropic loci are present and create _loucs_pleiotropic if yes
 * Perevioulsy the super chromosome _locus_tot was created
 *  unsigned int*   _locus_pleiotropic:
 *    - size: _locus_pleiotropic+1
 *    - contains the starting position of each real locus, i.e. each set of pleiotropic locus.
 *      This results that there is an entry for each locus, expect if they are pleiotroph
 *      Note, that only quantitative traits may show pleiotrophy
 *    - the last element contains the first position after the last locus
 */
void TGenomeProto::create_pleiotrophy()
{
	string type;
    
	// check if there are pleiotropic loci
	vector<unsigned int> locus_pleiotropic;	// temp vector since we do not know the size yet
	TLocusPos* curPos = NULL;
	unsigned int l;
	for (l = 0; l < _nb_locus_linked; ++l) {
		if (curPos == _locus_tot[l]->get_locusPosition()) continue; // it is a pleiotropic locus: continue
        curPos = _locus_tot[l]->get_locusPosition(); // overtake the new position
        locus_pleiotropic.push_back(l);
	}
    
    // if no pleiotropic loci are present stop here
    if(_nb_locus_linked == (unsigned int) locus_pleiotropic.size()) return;
    
    // we have pleiotropic loci:
    
	// add the index after the end (to be able to use "<")
	locus_pleiotropic.push_back(l);

    
	// copy the vector to an array
	_nb_locus_pleiotropic = (unsigned int)locus_pleiotropic.size();
	assert(!_locus_pleiotropic);
	_locus_pleiotropic = new unsigned int[_nb_locus_pleiotropic];
	ARRAY::vector2array<unsigned int>(locus_pleiotropic, _locus_pleiotropic, _nb_locus_pleiotropic);
	_nb_locus_pleiotropic -= 1; 	// there is one pleiotropic locus fewer
    
    // test if pleiotropic loci are valid
    check_pleiotrophic_loci();
}

//------------------------------------------------------------------------------
/** prints to the console the genetic map */
void TGenomeProto::printGeneticMapInfo()
{
	// find the number of linked/unlinked and pleiotropic/non-pleiotropic loci
	unsigned int l, nb_linked_pleiotropic = 0, nb_linked = 0;
	if (_locus_pleiotropic) {
		for (l = 0; l < _nb_locus_pleiotropic; ++l) {
			if (_locus_pleiotropic[l] + 1 == _locus_pleiotropic[l + 1]) ++nb_linked;
            else ++nb_linked_pleiotropic;
		}
	}
	else nb_linked = _nb_locus_linked;
    
	// linked loci
	if (nb_linked || nb_linked_pleiotropic) {       // linked loci
		message("\n    %u chromosomes of size", _nb_chromosome);
		if (_locus_position_tot[FEM] == _locus_position_tot[MAL]) {	// common genetic map
			double start = 0;
			for (unsigned int c = 0; c < _nb_chromosome; ++c) {
				if (c) message(",");
				message(" %g", _chromosomeSize[FEM][c] - start); // size is cumulative!
				start = _chromosomeSize[FEM][c];
			}
			message("cM (in total %gcM)", _chromosomeLength[FEM]);
		}
		else {                                       // sex specific genetic map
			double tot_size = 0, sizeF, sizeM, size;
			for (unsigned int c = 0; c < _nb_chromosome; ++c) {
				sizeM = c ? _chromosomeSize[MAL][c] - _chromosomeSize[MAL][c - 1] :
                _chromosomeSize[MAL][c]; // size is cumulative!
				sizeF = c ? _chromosomeSize[FEM][c] - _chromosomeSize[FEM][c - 1] :
                _chromosomeSize[FEM][c]; // size is cumulative!
				size = sizeM > sizeF ? sizeM : sizeF; // which chromosome is bigger?
				if (c) message(",");
				message(" %g", size);
				tot_size += size;
			}
			message("cM (in total %gcM, sex specific)", tot_size);
		}
		message("\n    %u linked loci", nb_linked + nb_linked_pleiotropic);
		if (nb_linked_pleiotropic) message(" (thereof %u pleiotropic loci)", nb_linked_pleiotropic);
	}
	
    // unlinked loci
    if (_nb_locus_unlinked) {
		message("\n    %u independent loci", _nb_locus_unlinked);
	}
    
    if(!_nb_locus_unlinked && !nb_linked) message("\n    CAUTION: no loci are simulated!");
    cout << flush;
}

// ----------------------------------------------------------------------------------------
// inherit
// ----------------------------------------------------------------------------------------
/** public function for inheritance */
void TGenomeProto::inherit(TIndividual* mother, TIndividual* father, unsigned char** child)
{
	(this->*_inherit_func_ptr)(mother, father, child);
}

// ----------------------------------------------------------------------------------------
// _inherit_unlinked
// ----------------------------------------------------------------------------------------
/** inheritance with just unlinked loci
 * the genome containes first the linked loci, followed by the unlinked loci
 */
void TGenomeProto::_inherit_unlinked(TIndividual* mother, TIndividual* father, unsigned char** child)
{
	unsigned char** mother_seq = mother->genome.get_sequence();
	unsigned char** father_seq = father->genome.get_sequence();
	for (unsigned int l = _nb_locus_linked; l < _nb_locus_tot; ++l) {
		child[l][0] = mother_seq[l][_popPtr->rand().Bool()];
		child[l][1] = father_seq[l][_popPtr->rand().Bool()];
	}
}

// ----------------------------------------------------------------------------------------
// _inherit_linked
// ----------------------------------------------------------------------------------------
/** inheritance with just linked loci */
void TGenomeProto::_inherit_linked(TIndividual* mother, TIndividual* father, unsigned char** child)
{
	(this->*_recombine_func_ptr[FEM])(mother, child, 0);
	(this->*_recombine_func_ptr[MAL])(father, child, 1);
}

// ----------------------------------------------------------------------------------------
// _inherit_mixed
// ----------------------------------------------------------------------------------------
/** inheritance where linked an unlinked loci are present
 * the genome contains first the linked loci, followed by the unlinked loci
 */
void TGenomeProto::_inherit_mixed(TIndividual* mother, TIndividual* father,
                                  unsigned char** child)
{
	_inherit_linked(mother, father, child);
	_inherit_unlinked(mother, father, child);
}

// ----------------------------------------------------------------------------------------
// _recombine
// ----------------------------------------------------------------------------------------
/** perform a recombination: the array gamete will contain the index of the chromosome to be used
 * The function is the same if pleiotropic loci are present or not.
 * - _chromosomeLength_temp: used to specify the total number of recombinations
 * - _locus_position_tot_temp: vector with the locus positions
 * */
void TGenomeProto::_recombine_normal(TIndividual* parent, unsigned char** child, int index)
{
	unsigned char** parent_seq = parent->genome.get_sequence();
	sex_t SEX = parent->getSex();
    
	// draw the recombination positions
	vector<double> vecRecombs;
	for (int i = _popPtr->rand().Poisson(_chromosomeLength_temp[SEX] / 100); i > 0; --i) {
		vecRecombs.push_back(_chromosomeLength_temp[SEX] * _popPtr->rand().Uniform());
	}
	sort(vecRecombs.begin(), vecRecombs.end());  	 	// sort the vector
	vecRecombs.push_back(_chromosomeLength_temp[SEX]); // add the last position of the chromosome as recombination event
    
	// inherit and recombine
	unsigned int c, l, curChrom;
	vector<double>::iterator nextRecomb = vecRecombs.begin(); // get the position of the first recombination event
	for (l = 0, c = 0; c < _nb_chromosome; ++c) {         // for each chromosome
		curChrom = (unsigned int) _popPtr->rand().Bool(); // draw randomly the starting gamete for the chromosome
		for (; l < _nb_locus_per_chromosome[c]; ++l) { // for each locus on this chromosome
			while (_locus_position_tot_temp[SEX][l] > *nextRecomb) { // the last position of recomb is always bigger or equal
				curChrom = !curChrom;
				++nextRecomb;
			}
			child[l][index] = parent_seq[l][curChrom]; // set the gamete for the locus
		}
	}
}

// ----------------------------------------------------------------------------------------
// _recombine
// ----------------------------------------------------------------------------------------
/** perform a recombination: the array gamete will contain the index of the chromosome to be used
 * The function is the same if pleiotropic loci are present or not.
 * - _chromosomeLength: used to specify the total number of recombinations
 * - _locus_position_tot: vector with the locus positions
 * Note that _chromosomeLength_temp and _locus_position_tot_temp are not used.
 * Each chromosome is treated separately
 * The recombination factor solely affects the number of recombination, thus the
 * positions afterwards can be drawn on the original chromosomes
 * "A chromosome exactly 100 centiMorgan long would be expected to experience a
 * single crossover per meiosis, on average."
 *
 * */
void TGenomeProto::_recombine_qtrait(TIndividual* parent, unsigned char** child, int index)
{
	unsigned char** parent_seq = parent->genome.get_sequence();
	sex_t SEX = parent->getSex();
    
	// inherit and recombine for each chromosome SEPARATELY
	unsigned int c, l, curChrom;
	double curFactor;
	vector<double> vecRecombs;
	vector<double>::iterator nextRecomb; // position of the next recombination even
	for (l = 0, c = 0; c < _nb_chromosome; ++c) {         // for each chromosome
		// draw the recombination positions
		vecRecombs.clear();
		curFactor = (this->*_recombination_chrom_func_ptr[SEX][c])(parent, SEX,
                                                                   c);
		if (curFactor < 0) error("Recombination_factor cannot be negative!\n");
		for (int i = _popPtr->rand().Poisson(
                                          curFactor * _chromosomeSize[SEX][c] / 100); i > 0; --i) {
			vecRecombs.push_back(_chromosomeSize[SEX][c] * _popPtr->rand().Uniform());// doesen't matter that the positions are drawn on a different sized chromosome
		}
		sort(vecRecombs.begin(), vecRecombs.end());  	 	// sort the vector
		vecRecombs.push_back(_chromosomeSize[SEX][c]); // add the last position of the chromosome as recombination event
        
		// perform the recombination
		nextRecomb = vecRecombs.begin();
		curChrom = (unsigned int) _popPtr->rand().Bool(); // draw randomly the starting gamete for the chromosome
		for (; l < _nb_locus_per_chromosome[c]; ++l) { // for each locus on this chromosome
			while (_locus_position_tot[SEX][l] > *nextRecomb) { // the last position of recomb is always bigger or equal
				curChrom = !curChrom;
				++nextRecomb;
			}
			child[l][index] = parent_seq[l][curChrom]; // set the gamete for the locus
		}
	}
}

// ----------------------------------------------------------------------------------------
// ini_all
// ----------------------------------------------------------------------------------------
/** initialization of the different genome aspects
 */
void TGenomeProto::ini_all()
{
    MUTATION_TRIAL = (unsigned int)this->get_parameter_value("mutation_trial");
	ini_genetic_map();
	ini_recombination_factor();
	ini_mutate(); // this function has to be called after ini_genetic_map() has been executed
}

// ----------------------------------------------------------------------------------------
// ini_recombination_factor
// ----------------------------------------------------------------------------------------
/** initialization of the recombination factor
 */
void TGenomeProto::ini_recombination_factor()
{
	delete_sex_specific_array(_locus_position_tot_temp[FEM],
                              _locus_position_tot_temp[MAL]);
    
	// sex specific recombination factors used?
	if (ini_recombination_factor("recombination_factor_fem", FEM)) {
		if (_popPtr->get_sexInitRatio()) {
			// get the male values
			if (!ini_recombination_factor("recombination_factor_mal", MAL)) {
				error("Only one sex specific recombination_factor is specified: both are required!\n");
			}
		}
		else {
			_recombine_func_ptr[MAL] = _recombine_func_ptr[FEM];
			_chromosomeLength_temp[MAL] = _chromosomeLength_temp[FEM];
			_locus_position_tot_temp[MAL] = _locus_position_tot_temp[FEM];
            
		}
	}
	else {		// general recombination factor used?
		if (ini_recombination_factor("recombination_factor", FEM)) {
			ini_recombination_factor("recombination_factor", MAL);// male and female have the same recombination factor
		}
		else {	// no recombination factor used
			_recombine_func_ptr[FEM] = &TGenomeProto::_recombine_normal;
			_chromosomeLength_temp[FEM] = _chromosomeLength[FEM];
			_locus_position_tot_temp[FEM] = _locus_position_tot[FEM];
            
			_recombine_func_ptr[MAL] = &TGenomeProto::_recombine_normal;
			_chromosomeLength_temp[MAL] = _chromosomeLength[MAL];
			_locus_position_tot_temp[MAL] = _locus_position_tot[MAL];
		}
	}
	if (!_recombination_factor_chrom[MAL]) _recombination_factor_chrom[MAL] =
        _recombination_factor_chrom[FEM];
	if (!_recombination_chrom_func_ptr[MAL])
		_recombination_chrom_func_ptr[MAL] = _recombination_chrom_func_ptr[FEM];
}

// ----------------------------------------------------------------------------------------
// ini_recombination_factor
// ----------------------------------------------------------------------------------------
/** initialization of the recombination factor
 * returns true if the parameter is set and false if not
 * sets the following parameters
 *  - _locus_position_tot_temp[SEX][l]: loci positions on the super chromosome
 *  - _chromosomeLength_temp[SEX]: total size of the super chromosome
 * */
bool TGenomeProto::ini_recombination_factor(string param_name, sex_t SEX)
{
	Param* param = get_parameter(param_name);
	if (!param->isSet()) return false;
    
	if (param->is_matrix()) {	// recombination factor varies per chromosome
		TMatrix* m = param->get_as_matrix(true);
		double* vec = m->get();
		unsigned int size = m->length();

		// check the number of values per number of chromosomes
		if (size > _nb_chromosome) warning("There are more %s specified than chromosomes! Only a part of the %s is considered!\n",
                                           param_name.c_str(), param_name.c_str());
		else if (_nb_chromosome % size)
			warning("The number of %s is not an entire subset of the number of chromosomes!\n",
					param_name.c_str());
        
		// if strings were found within the matrix
        if (m->get_strMatrix()){
			ini_recombination_qtrait(param_name, SEX, vec, size, m);
            delete m;
            return true;
        }
        
		// get the values
        double recombFactor, startChrom_temp = 0, startChrom = 0;
        _locus_position_tot_temp[SEX] = new double[_nb_locus_linked];
        for (unsigned int l = 0, c = 0; c < _nb_chromosome; ++c) {// for each chromosome
            recombFactor = vec[c % size];
            for (; l < _nb_locus_per_chromosome[c]; ++l) {	// for each locus
                _locus_position_tot_temp[SEX][l] = startChrom_temp
                + recombFactor
                * (_locus_position_tot[SEX][l] - startChrom);
            }
            startChrom_temp = _locus_position_tot_temp[SEX][l - 1];
            startChrom = _locus_position_tot[SEX][l - 1];
        }
        _chromosomeLength_temp[SEX] = startChrom_temp;
        _recombine_func_ptr[SEX] = &TGenomeProto::_recombine_normal;
        delete m;
        return true;
	}
	
    // single value
    double recombFactor = param->get_value();
    unsigned int l;
    _locus_position_tot_temp[SEX] = new double[_nb_locus_linked];
    for (l = 0; l < _nb_locus_linked; ++l) {		// for each locus
        _locus_position_tot_temp[SEX][l] = recombFactor
        * _locus_position_tot[SEX][l];
    }
    _chromosomeLength_temp[SEX] = _locus_position_tot_temp[SEX][l - 1];
    _recombine_func_ptr[SEX] = &TGenomeProto::_recombine_normal;
	return true;
}

// ----------------------------------------------------------------------------------------
// ini_recombination_qtrait
// ----------------------------------------------------------------------------------------
/** initialization of the recombination factor if at least one recombination factor
 * is defined by a recombination factor
 * sets the following parameters
 *  - _locus_position_tot_temp[SEX][l]: loci positions on the super chromosome
 *  - _chromosomeLength_temp[SEX]: total size of the super chromosome
 * */
void TGenomeProto::ini_recombination_qtrait(string param_name, sex_t SEX,
                                            double* vec, unsigned int size, TMatrix* m)
{
	map<unsigned int, map<unsigned int, string> >* m_str = m->get_strMatrix();
	unsigned int pos, nbCol = m->getNbCols();
	double recombFactor;
	string text;
	vector<int> traitID = _popPtr->getTraitIndex("quanti");	// containing the absolute index across all types of traits
	_recombination_factor_chrom[SEX] = new double[_nb_chromosome];
	_recombination_chrom_func_ptr[SEX] = new _func_ptr[_nb_chromosome];
	for (unsigned int c = 0; c < _nb_chromosome; ++c) {	// for each chromosome
		pos = c % size;
		recombFactor = vec[pos];
		if (recombFactor == my_STR) {// a trait specifies the recombination factor
			text = (*m_str)[(unsigned int) pos / nbCol][pos % nbCol];
			switch (text[0]) {
				case 'Z': {	// the phenotype specifies the recombination factor
					_recombination_chrom_func_ptr[SEX][c] =
                    &TGenomeProto::_recombination_chrom_qtraitZ;
					break;
				}
				case 'G': {	// the genotypic value specifies the recombination factor
					_recombination_chrom_func_ptr[SEX][c] =
                    &TGenomeProto::_recombination_chrom_qtraitG;
					break;
				}
				default:
					error("The parameter %s has not readable characters!\n",
                          param_name.c_str());
					break;
			}
            
			// get the absolute trait index (caution may have different types of traits)
			pos = strTo<unsigned int>(text.substr(1)) - 1;	// starting at 1
			if (pos >= traitID.size())
				error("The parameter %s has a quantitative trait index which is out of range!\n",
                      param_name.c_str());
			_recombination_factor_chrom[SEX][c] = traitID[pos];
		}
		else {				// a fixed number specifies the recombination factor
			_recombination_chrom_func_ptr[SEX][c] =
            &TGenomeProto::_recombination_chrom_factor;
			_recombination_factor_chrom[SEX][c] = recombFactor;
		}
	}
	_recombine_func_ptr[SEX] = &TGenomeProto::_recombine_qtrait;
}

// ----------------------------------------------------------------------------------------
// _recombination_chrom_factor
// ----------------------------------------------------------------------------------------
/** a single recombination factor per chromosome is present
 */
double TGenomeProto::_recombination_chrom_factor(TIndividual* parent, sex_t SEX,
                                                 unsigned int chrom)
{
	return _recombination_factor_chrom[SEX][chrom];
}

// ----------------------------------------------------------------------------------------
// _recombination_chrom_qtraitZ
// ----------------------------------------------------------------------------------------
/** recombination factor is specified by the phenotype of trait X
 */
double TGenomeProto::_recombination_chrom_qtraitZ(TIndividual* parent, sex_t SEX,
                                                  unsigned int chrom)
{
	return parent->getTraitPhenotype(_recombination_factor_chrom[SEX][chrom]);
}

// ----------------------------------------------------------------------------------------
// _recombination_chrom_qtraitG
// ----------------------------------------------------------------------------------------
/** recombination factor is specified by the genotypic value of trait X
 */
double TGenomeProto::_recombination_chrom_qtraitG(TIndividual* parent, sex_t SEX,
                                                  unsigned int chrom)
{
	return parent->getTraitGenotype(_recombination_factor_chrom[SEX][chrom]);
}

// ----------------------------------------------------------------------------------------
// ini_mutate
// ----------------------------------------------------------------------------------------
/** initialization of the mutation function
 * this function has to be called after ini_genetic_map() has been executed
 */
void TGenomeProto::ini_mutate()
{
	if (!_change_mutationRates || !_nb_locus_tot) return;
	_change_mutationRates = false;
    
	// check if the mutation rate changes among loci
	double mut_rate = _locus_tot[0]->get_mut_rate(); // mutation rate of the first locus
	unsigned int l;
	for (l = 1; l < _nb_locus_tot; ++l) {
		if (mut_rate != _locus_tot[l]->get_mut_rate()) { // check if the rate changes among loci
			if (_locus_pleiotropic) _mutate_func_ptr = &TGenomeProto::_mutate_unequal_mutation_rate_pleiotrophy;
			else _mutate_func_ptr = &TGenomeProto::_mutate_unequal_mutation_rate;
			return;
		}
	}
    
	// all loci have the same mutation rate
	_mut_rate_mean = mut_rate;
	if (!_mut_rate_mean) _mutate_func_ptr =
        &TGenomeProto::_mutate_zero_mutation_rate; 	// no mutation
	else {
		if (_locus_pleiotropic) _mutate_func_ptr =
            &TGenomeProto::_mutate_equal_mutation_rate_pleiotrophy; // same mutation rate
		else _mutate_func_ptr = &TGenomeProto::_mutate_equal_mutation_rate; // same mutation rate
	}
}

// ----------------------------------------------------------------------------------------
//
// ----------------------------------------------------------------------------------------
/** check if the loci of a pleiotropic set meet the requirements:
 * - equal trait type
 * - equal mutation rate
 * - equal mutation model
 */
void TGenomeProto::check_pleiotrophic_loci()
{
	if (!_locus_pleiotropic) return;

	unsigned int cur, end;
	double mut_rate;
	mut_model_t mut_model;
	for (unsigned int l = 0; l < _nb_locus_pleiotropic; ++l) { // for each locus
		end = _locus_pleiotropic[l + 1];
		if (_locus_pleiotropic[l] == end - 1) continue;       // not a pleiotropic locus: continue
        
		// we have a pleiotropic set => test it
        cur = _locus_pleiotropic[l];
		mut_rate = _locus_tot[cur]->get_mut_rate();        // get the first mutation rate
		mut_model = _locus_tot[cur]->get_mutationModel();  // get the first mutation model
		for (; cur < end; ++cur) {          // for each successiv pleiotropic locus
            if (_locus_tot[cur]->get_pTrait()->get_type()!="quanti")
				error("Genome: pleiotropic loci must belong to a quantitatvie trait type!\n");
			if (mut_rate != _locus_tot[cur]->get_mut_rate())
				error("Genome: pleiotropic loci must have the same mutation rate!\n");
			if (mut_model != _locus_tot[cur]->get_mutationModel())
				error("Genome: pleiotropic loci must have the same mutation model!\n");
		}
	}
}

// ----------------------------------------------------------------------------------------
// mutate the entire sequence
// ----------------------------------------------------------------------------------------
/** public function for mutate of the entire sequence
 * the function redirects the call to either of the following functions
 *  - _mutate_zero_mutation_rate
 *  - _mutate_equal_mutation_rate
 *  - _mutate_unequal_mutation_rate
 */
void TGenomeProto::mutate(unsigned char** seq)
{
	(this->*_mutate_func_ptr)(seq);
}

// ----------------------------------------------------------------------------------------
/** if the mutation rate is identical for all loci perform a "global" mutation
 * draw randomly the locus and allele to mutate and then call the corresponding mutation function
 */
void TGenomeProto::_mutate_equal_mutation_rate(unsigned char** seq)
{
	unsigned int NbMut, l;
	for (NbMut = _popPtr->rand().Poisson(ploidy * _nb_locus_tot * _mut_rate_mean); NbMut != 0; --NbMut) {
		l = _popPtr->rand().Uniform(_nb_locus_tot);
		_locus_tot[l]->mutate_now(&seq[l][_popPtr->rand().Uniform((unsigned int)ploidy)]); // a mutation has to occur
	}
}

// ----------------------------------------------------------------------------------------
/** if the mutation rate is identical for all loci perform a "global" mutation
 * draw randomly the locus and allele to mutate and then call the corresponding mutation function
 */
void TGenomeProto::_mutate_equal_mutation_rate_pleiotrophy(unsigned char** seq)
{
	unsigned int NbMut, l, curLocus, end;
	for (NbMut = _popPtr->rand().Poisson(ploidy * _nb_locus_pleiotropic * _mut_rate_mean); NbMut != 0; --NbMut) {
		l = _popPtr->rand().Uniform(_nb_locus_pleiotropic);
		curLocus = _locus_pleiotropic[l];           // get the starting position
		end = _locus_pleiotropic[l + 1];            // get the after last position
		for (; curLocus < end; ++curLocus) {
			_locus_tot[curLocus]->mutate_now(&seq[curLocus][_popPtr->rand().Uniform((unsigned int)ploidy)]); // a mutation has to occur
		}
	}
}

// ----------------------------------------------------------------------------------------
/** if the mutation rate is identical for all loci perform a "global" mutation
 * draw randomly the locus and allele to mutate and then call the corresponding mutation function
 */
void TGenomeProto::_mutate_equal_mutation_rate_pleiotrophy_correl(unsigned char** seq)
{
    assert(_locus_pleiotropic);
	unsigned int NbMut, l, curPos, endPos;
	double* MNR_params, *deviates;
	for (NbMut = _popPtr->rand().Poisson(ploidy * _nb_locus_pleiotropic * _mut_rate_mean); NbMut != 0;--NbMut) {
		l = _popPtr->rand().Uniform(_nb_locus_pleiotropic);
		curPos = _locus_pleiotropic[l];              // get the starting position
		endPos = _locus_pleiotropic[l + 1];         // get the after last position
		MNR_params = _locus_tot[curPos]->get_MNR_param();
        
		if (MNR_params) {    // if the mutations are correlated
			deviates = _popPtr->rand().get_MNR_deviates_fromParam(MNR_params, endPos - curPos);   // draw the random numbers
			for (; curPos < endPos; ++curPos) {
				_locus_tot[curPos]->mutate_now(&seq[curPos][_popPtr->rand().Uniform((unsigned int)ploidy)],*deviates++); // a mutation has to occur
			}
		}
		else {	// if they are uncorrelated the mutations are independent
			for (; curPos < endPos; ++curPos) {
				_locus_tot[curPos]->mutate_now(&seq[curPos][_popPtr->rand().Uniform((unsigned int)ploidy)]); // a mutation has to occur
			}
		}
	}
}

// ----------------------------------------------------------------------------------------
/** if the mutation rate varies between loci the probability to mutate has to be
 * checked for each single locus
 */
void TGenomeProto::_mutate_unequal_mutation_rate(unsigned char** seq)
{
	for (unsigned int l = 0; l < _nb_locus_tot; ++l) {
		_locus_tot[l]->mutate(seq[l]);       // check if a mutation appears here
	}
}

// ----------------------------------------------------------------------------------------
/** if the mutation rate varies between loci the probability to mutate has to be
 * checked for each single locus
 */
void TGenomeProto::_mutate_unequal_mutation_rate_pleiotrophy(
                                                             unsigned char** seq)
{
    assert(_locus_pleiotropic);
	unsigned int NbMut, l, curLocus, end;
	for (l = 0; l < _nb_locus_pleiotropic; ++l) {
		for (NbMut = _popPtr->rand().Poisson(ploidy * _locus_tot[l]->get_mut_rate()); NbMut != 0; --NbMut) {
			curLocus = _locus_pleiotropic[l];
			end = _locus_pleiotropic[l + 1];
			for (; curLocus < end; ++curLocus) {
				_locus_tot[curLocus]->mutate_now(&seq[curLocus][_popPtr->rand().Uniform((unsigned int)ploidy)]); // a mutation has to occur
			}
		}
	}
}

// ----------------------------------------------------------------------------------------
/** if the mutation rate varies between loci the probability to mutate has to be
 * checked for each single locus
 * The mutation between the different loci of a pleiotropic set are correlated,
 * thus the random deviates are drawn simultaneously from the multivariate normal distribution
 */
void TGenomeProto::_mutate_unequal_mutation_rate_pleiotrophy_correl(unsigned char** seq)
{
    assert(_locus_pleiotropic);
	unsigned int NbMut, l, curPos, endPos;
	double* MNR_params, *deviates;
	for (l = 0; l < _nb_locus_pleiotropic; ++l) {
		for (NbMut = _popPtr->rand().Poisson(ploidy * _locus_tot[l]->get_mut_rate()); NbMut != 0; --NbMut) {
			curPos = _locus_pleiotropic[l];
			endPos = _locus_pleiotropic[l + 1];
			MNR_params = _locus_tot[curPos]->get_MNR_param();
            
			if (MNR_params) {    // if the mutations are correlated
				deviates = _popPtr->rand().get_MNR_deviates_fromParam(MNR_params,
                                                                   endPos - curPos);   // draw the random numbers
				for (; curPos < endPos; ++curPos) {
					_locus_tot[curPos]->mutate_now(&seq[curPos][_popPtr->rand().Uniform((unsigned int)ploidy)], *deviates++); // a mutation has to occur
				}
			}
			else {	// if they are uncorrelated the mutations are independent
				for (; curPos < endPos; ++curPos) {
					_locus_tot[curPos]->mutate_now(&seq[curPos][_popPtr->rand().Uniform((unsigned int)ploidy)]); // a mutation has to occur
				}
			}
		}
	}
}

// ----------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------
// set_ini_sequence_of_locus
// ----------------------------------------------------------------------------------------
/** set the sequence initialization model for all loci of the trait called by TTraitProto*/
void TGenomeProto::set_ini_sequence_model(TLocus* aLocus,
                                          const unsigned int& size, const unsigned int& model)
{
	for (unsigned int l = 0; l < size; ++l) {
		aLocus[l].set_ini_sequence_model(model);
	}
}

// ----------------------------------------------------------------------------------------
/** set the sequence initialization model for all loci of the trait called by TTraitProto*/
void TGenomeProto::set_ini_sequence_model(TLocus* aLocus,
                                          const unsigned int& size, const ini_model_t& model)
{
	for (unsigned int l = 0; l < size; ++l) {
		aLocus[l].set_ini_sequence_model(model);
	}
}

// ----------------------------------------------------------------------------------------
// ini_sequence
// ----------------------------------------------------------------------------------------
/** initialize the entire sequence */
void TGenomeProto::ini_sequence(unsigned char** seq, TPatch* patch)
{
	assert(seq);
    
	for (unsigned int l = 0; l < _nb_locus_tot; ++l) {
		_locus_tot[l]->ini_sequence(seq[l], patch);
	}
}

// ----------------------------------------------------------------------------------------
// executeBefortemporal_changeEachReplicate
// ----------------------------------------------------------------------------------------
void TGenomeProto::temporal_change(const unsigned int& gen)
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
				if (pos->first.find("recombination_factor") != string::npos
                    || pos->first.find("recombination_factor_fem") != string::npos
                    || pos->first.find("recombination_factor_mal") != string::npos) {
					ini_recombination_factor();
				}
			}
		}
	}
    
	ini_mutate();
}

// ----------------------------------------------------------------------------------------
// clear
// ----------------------------------------------------------------------------------------
void TGenomeProto::clear()
{
	delete_sex_specific_array(_recombination_chrom_func_ptr[FEM],
                              _recombination_chrom_func_ptr[MAL]);
	delete_sex_specific_array(_recombination_factor_chrom[FEM],
                              _recombination_factor_chrom[MAL]);
	delete_sex_specific_array(_locus_position_tot[FEM],
                              _locus_position_tot[MAL]);
	delete_sex_specific_array(_chromosomeSize[FEM], _chromosomeSize[MAL]);
	if (_locus_tot) {
		delete[] _locus_tot;
		_locus_tot = NULL;
	}
	if (_locus_pleiotropic) {
		delete[] _locus_pleiotropic;
		_locus_pleiotropic = NULL;
	}
	if (_nb_locus_per_chromosome) {
		delete[] _nb_locus_per_chromosome;
		_nb_locus_per_chromosome = NULL;
	}
    
	// only delete array if recombination factor was used
	if (_recombine_func_ptr[FEM] != &TGenomeProto::_recombine_normal) {
		if (_locus_position_tot_temp[FEM]) {
			if (_locus_position_tot_temp[FEM]
                == _locus_position_tot_temp[MAL]) {
				delete[] _locus_position_tot_temp[FEM];
				_locus_position_tot_temp[FEM] = _locus_position_tot_temp[MAL] =
                NULL;
			}
			delete[] _locus_position_tot_temp[FEM];
			_locus_position_tot_temp[FEM] = NULL;
		}
	}
    
	if (_recombine_func_ptr[MAL] != &TGenomeProto::_recombine_normal) {
		if (_locus_position_tot_temp[MAL]) {
			delete[] _locus_position_tot_temp[MAL];
			_locus_position_tot_temp[MAL] = NULL;
		}
	}
}


// ----------------------------------------------------------------------------------------
// create_locus_tot
// ----------------------------------------------------------------------------------------
void TGenomeProto::create_locus_tot()
{
    if(_locus_pos.find("global")==_locus_pos.end()) create_locus_tot_perTrait();
    else create_locus_tot_global();
}

// ----------------------------------------------------------------------------------------
// create_locus_tot
// ----------------------------------------------------------------------------------------
/** aligns all loci from the vector _locus_vector to the array _locus_tot
 * and adds the genetic position to the locus.
 * The vector _locus_vector is cleared at the end.
 * _locus_pos[type][TLocusPos*]
 * _locus_vector[type][locus_id][TLocus*]
 * returns true if pleiotropic loci are present and false if this is not the case
 * since the locus position index are separately for each type of trait first a
 * temp vector has to be created before the final _locus_tot vector may be created
 * The following objects are created/filled:
 * - _locus_tot[TLocus*]
 * - _nb_locus_linked: unsigned int
 * - _nb_locus_unlinked: unsigned int
 * - _nb_chromosome: unsigned int
 * - _nb_locus_per_chromosome[chromosome]: cumulative number of loci per chromosome
 */
void TGenomeProto::create_locus_tot_perTrait()
{
    assert(_locus_pos.find("global")==_locus_pos.end());
    
	// reset class variables
	if (_locus_tot) delete[] _locus_tot;
	_locus_tot = new TLocus*[_nb_locus_tot];
	_nb_locus_linked = 0;            // reset the number of linked loci
    
	// create some temp containers
	map<TLocusPos, vector<TLocus*> > locus_tot_trait; // used to align the different trait types
	vector<unsigned int> locus_pos_used; // container to check if all locus positions are used
	map<string, map<unsigned int, vector<TLocus*> > >::iterator curType, endType; // _locus_vector[type][locusID]
	map<unsigned int, vector<TLocus*> >::iterator curLocIdx, endLocIdx;
	vector<TLocus*>::iterator curLocus, endLocus;
	map<string, vector<TLocusPos*> >::iterator curGenomeMap;
	vector<TLocusPos*>* curGenome = &(_locus_pos.find("global")->second);  // get the current genome;
	TLocusPos* curLocusPos;
	unsigned int l, posIdx;
	string type;                     // ntrl/quanti
    
	curType = _locus_vector.begin();
	endType = _locus_vector.end();
	for (; curType != endType; ++curType) {            // for each type of trait
		type = curType->first;                             // "quanti" or "ntrl"
		assert(_locus_pos.find(type) != _locus_pos.end()); // at least the NaN position must be available for this trait type
        curGenome = &(_locus_pos.find(type)->second);  // get the current genome;
		locus_pos_used.clear();
		locus_pos_used.assign(curGenome->size(), 0);
		curLocIdx = curType->second.begin();
		endLocIdx = curType->second.end();
		for (; curLocIdx != endLocIdx; ++curLocIdx) {    // for each locus index
			curLocus = curLocIdx->second.begin();
			endLocus = curLocIdx->second.end();
            
			// get the locus index and find the correct position
			posIdx = curLocIdx->first; // locus position index of the current locus
			if (posIdx == my_NAN) {                // locus position is unlinked
				curLocusPos = *(*curGenome).rbegin(); // get the last position , it is the my_NaN
				++locus_pos_used[curGenome->size() - 1]; // count the number of times this locus is used
				_nb_locus_unlinked += curLocIdx->second.size(); // get number of unlinked loci (unlinked loci have as index NaN!)
			}
			else {                                   // locus position is linked
				if (posIdx >= curGenome->size() - 1) { // there is a problem with the index (was no genome specified?)
					if (curGenome->size() == 1) error("Genome (%s): no genome is specified for this trait!\n", type.c_str());
					else error("Genome (%s): locus position index '%u' of trait type '%s' exceeds the genome range!\n", type.c_str(), posIdx + 1, curType->first.c_str());
                
				}
				curLocusPos = (*curGenome)[posIdx]; // get the current locus position object
				++locus_pos_used[posIdx]; // count the number of times this locus is used
			}
            
			// set the genetic position of the locus depending on its locus index
			for (; curLocus != endLocus; ++curLocus) { // for each pleiotropic locus (locus of this index)
				(*curLocus)->set_locusPosition(curLocusPos); // set the genetic position to the locus
				locus_tot_trait[*curLocusPos].push_back(*curLocus); // add the locus to the trait vector
			}
		}
        
		// check if all locus positions are used and erase unused positions (start from the end!!!)
		vector<unsigned int> posIndex;
		for (l = (unsigned int)locus_pos_used.size(); l; --l) {
			if (locus_pos_used[l - 1]) continue; // locus position was used at least once
			if (l != locus_pos_used.size()) posIndex.push_back(l); // the last position is the unlinked position which may not always be used!
			delete (*curGenome)[l - 1]; // remove this locus position as it is not used
			(*curGenome).erase(curGenome->begin() + l - 1); // remove this locus position as it is not used
		}
		if (!posIndex.empty()) { // if not all positions are used plot a warning
			string text;
			for (l = (unsigned int)posIndex.size(); l; --l) {
				text += " " + toStr(posIndex[l - 1]);
			}
			warning("No loci were attributed to the positions%s of the %s genome!\n", text.c_str(), type.c_str());
		}
	}
	_nb_locus_linked = _nb_locus_tot - _nb_locus_unlinked; // set the number of linked loci
    
	// merge the different types of traits and create the final array _locus_tot
	map<TLocusPos, vector<TLocus*> >::iterator curPos = locus_tot_trait.begin(),
    endPos = locus_tot_trait.end();
	vector<TLocus*>::iterator curLoc, endLoc;
	unsigned int curChrom = curPos->first.chrom; // temp vector to count the number of cumulative loci per chromosome
	vector<unsigned int> nb_locus_per_chrom;
	for (l = 0; curPos != endPos; ++curPos) {         // for each locus position
                                                      // if a new chromosome starts
		if (curChrom != curPos->first.chrom) {
			curChrom = curPos->first.chrom;                // overtake the index
			nb_locus_per_chrom.push_back(l); // add the cumulative number of loci of this chromosome
		}
		// for each locus of this position
		for (curLoc = curPos->second.begin(), endLoc = curPos->second.end();
             curLoc != endLoc; ++curLoc, ++l) {
			_locus_tot[l] = *curLoc;
			(*curLoc)->set_locus_id_tot(l);    // set the index within the locus
		}
	}
	if (curChrom != my_NAN) nb_locus_per_chrom.push_back(l); // the last chromosome (caution about unlinked loci at the end!!!)
    
	if (l != _nb_locus_tot) error("Genome: the number of loci is incorrect!\n");
	if (_nb_locus_tot == _nb_locus_unlinked) _nb_chromosome = 0; // no chromosomes
	else ARRAY::vector2array<unsigned int>(nb_locus_per_chrom, _nb_locus_per_chromosome, _nb_chromosome);
}

// ----------------------------------------------------------------------------------------
// create_locus_tot
// ----------------------------------------------------------------------------------------
/** aligns all loci from the vector _locus_vector to the array _locus_tot
 * and adds the genetic position to the locus.
 * The vector _locus_vector is cleared at the end.
 * _locus_pos[type][TLocusPos*]
 * _locus_vector[type][locus_id][TLocus*]
 * since the locus position index are separately for each type of trait first a
 * temp vector has to be created before the final _locus_tot vector may be created
 * The following objects are created/filled:
 * - _locus_tot[TLocus*]
 * - _nb_locus_linked: unsigned int
 * - _nb_locus_unlinked: unsigned int
 * - _nb_chromosome: unsigned int
 * - _nb_locus_per_chromosome[chromosome]: cumulative number of loci per chromosome
 */
void TGenomeProto::create_locus_tot_global()
{
    assert(_locus_pos.find("global")!=_locus_pos.end());

	// reset class variables
	assert(!_locus_tot);
	_locus_tot = new TLocus*[_nb_locus_tot];
	_nb_locus_linked = 0;            // reset the number of linked loci
    
	// create some temp containers
	map<TLocusPos, vector<TLocus*> > locus_tot_trait; // used to align the different trait types
	vector<unsigned int> locus_pos_used; // container to check if all locus positions are used and to evaluate if pleiotropic loci are present
    
	map<string, map<unsigned int, vector<TLocus*> > >::iterator curType, endType;
	map<unsigned int, vector<TLocus*> >::iterator curLocIdx, endLocIdx;
	vector<TLocus*>::iterator curLocus, endLocus;
	map<string, vector<TLocusPos*> >::iterator curGenomeMap;
	vector<TLocusPos*>* curGenome = &(_locus_pos.find("global")->second);  // get the current genome;
    locus_pos_used.assign(curGenome->size(), 0);
	TLocusPos* curLocusPos;
	unsigned int l, posIdx;
    
	curType = _locus_vector.begin();
	endType = _locus_vector.end();
	for (; curType != endType; ++curType) {            // for each type of trait
		curLocIdx = curType->second.begin();
		endLocIdx = curType->second.end();
		for (; curLocIdx != endLocIdx; ++curLocIdx) {    // for each locus index
			curLocus = curLocIdx->second.begin();
			endLocus = curLocIdx->second.end();
            
			// get the locus index and find the correct position
			posIdx = curLocIdx->first; // locus position index of the current locus
			if (posIdx == my_NAN) {                // locus position is unlinked
				posIdx = (unsigned int) curGenome->size() - 1;
                assert( (*curGenome)[posIdx]->pos[FEM]=my_NAN);
				_nb_locus_unlinked += curLocIdx->second.size(); // get number of unlinked loci (unlinked loci have as index NaN!)
			}
			else {                                   // locus position is linked
				if (posIdx >= curGenome->size() - 1) { // there is a problem with the index (was no genome specified?)
					if (curGenome->size() == 1) error("Genome (global): no genome is specified for this trait!\n");
					else error("Genome (global): locus position index '%u' of trait type '%s' exceeds the genome range!\n", posIdx + 1, curType->first.c_str());
				}
            }
				
            // add the locus to the map
            curLocusPos = (*curGenome)[posIdx]; // get the current locus position object
            ++locus_pos_used[posIdx]; // count the number of times this locus is used
			         
			// set the genetic position of the locus depending on its locus index
			for (; curLocus != endLocus; ++curLocus) { // for each pleiotropic locus (locus of this index)
				(*curLocus)->set_locusPosition(curLocusPos); // set the genetic position to the locus
				locus_tot_trait[*curLocusPos].push_back(*curLocus); // add the locus to the trait vector
                
                
//                cout << endl << curType->first << "\t" << curLocIdx->first << "\t" << (*curLocus)->get_locus_id_trait() << "\t" << *curLocus << "\t" << (*curLocus)->get_locusPosition() << "\t" << (*curLocus)->get_chromosomePosition() << "\t" << (*curLocus)->get_locusPosition(FEM);
			}
//            cout << endl;
		}
    }
    
    
    // check if all locus positions are used and erase unused positions (start from the end!!!)
    vector<unsigned int> posIndex;
    for (l = (unsigned int)locus_pos_used.size(); l; --l) {
        if (locus_pos_used[l - 1]) continue; // locus position was used at least once
        if (l != locus_pos_used.size()) posIndex.push_back(l); // the last position is the unlinked position which may not always be used!
        delete (*curGenome)[l - 1]; // remove this locus position as it is not used
        (*curGenome).erase(curGenome->begin() + l - 1); // remove this locus position as it is not used
    }
    if (!posIndex.empty()) { // if not all positions are used plot a warning
        string text;
        for (l = (unsigned int)posIndex.size(); l; --l) {
            text += " " + toStr(posIndex[l - 1]);
        }
        warning("No loci were attributed to the positions%s of the global genome!\n", text.c_str());
    }
    
    _nb_locus_linked = _nb_locus_tot - _nb_locus_unlinked; // set the number of linked loci
    
	// merge the different types of traits and create the final array _locus_tot
	map<TLocusPos, vector<TLocus*> >::iterator curPos = locus_tot_trait.begin();
    map<TLocusPos, vector<TLocus*> >::iterator endPos = locus_tot_trait.end();
	vector<TLocus*>::iterator curLoc, endLoc;
	unsigned int curChrom = curPos->first.chrom; // temp vector to count the number of cumulative loci per chromosome
	vector<unsigned int> nb_locus_per_chrom;
	for (l = 0; curPos != endPos; ++curPos) {         // for each locus position
		// if a new chromosome starts
		if (curChrom != curPos->first.chrom) {
			curChrom = curPos->first.chrom;                // overtake the index
			nb_locus_per_chrom.push_back(l); // add the cumulative number of loci of this chromosome
		}
        
		// add each locus of this position to _locus_tot
		for (curLoc = curPos->second.begin(), endLoc = curPos->second.end(); curLoc != endLoc; ++curLoc, ++l) {
			_locus_tot[l] = *curLoc;
			(*curLoc)->set_locus_id_tot(l);    // set the index within the locus
		}
	}
	if (curChrom != my_NAN) nb_locus_per_chrom.push_back(l); // the last chromosome (caution about unlinked loci at the end!!!)
    
	if (l != _nb_locus_tot) error("Genome: the number of loci is incorrect!\n");
	if (_nb_locus_tot == _nb_locus_unlinked) _nb_chromosome = 0; // no chromosomes
	else ARRAY::vector2array<unsigned int>(nb_locus_per_chrom, _nb_locus_per_chromosome, _nb_chromosome);
}

// ----------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------
// create_locus_position_tot
// ----------------------------------------------------------------------------------------
/** create the following objects:
 *  - _locus_position_tot[SEX][l]: loci positions on the super chromosome
 *  - _chromosomeSize[SEX][c]: cumulative size of the chromosomes
 *  - _chromosomeLength[SEX]: total length of the super chromosome
 * the objects are created based on the array _locus_tot
 */
void TGenomeProto::create_locus_position_tot()
{
	if (!_nb_chromosome) return;
    
	// reset the arrays
	delete_sex_specific_array(_chromosomeSize[FEM], _chromosomeSize[MAL]);
	delete_sex_specific_array(_locus_position_tot[FEM],
                              _locus_position_tot[MAL]);
	delete_sex_specific_array(_locus_position_tot_temp[FEM],
                              _locus_position_tot_temp[MAL]);
    
	// create the arrays
	unsigned int nbSex;
	if (_sex_specific_genome) {
		for (unsigned s = 0; s < 2; ++s) {
			_locus_position_tot[s] = new double[_nb_locus_linked];
			_chromosomeSize[s] = new double[_nb_chromosome];
		}
		nbSex = 2;
	}
	else {
		_locus_position_tot[FEM] = new double[_nb_locus_linked];
		_chromosomeSize[FEM] = new double[_nb_chromosome];
		_locus_position_tot[MAL] = _locus_position_tot[FEM];
		_chromosomeSize[MAL] = _chromosomeSize[FEM];
		nbSex = 1;
	}
    
	// get the cumulative position for each locus
	double curSize[2];
	curSize[FEM] = curSize[MAL] = 0;
	double startChrom[2];
	startChrom[FEM] = startChrom[MAL] = 0;
	unsigned int s, l, c = 0, curChrom = _locus_tot[0]->get_chromosomePosition(); // first chromosome
	for (l = 0; l < _nb_locus_linked; ++l) {            // for each linked locus
		assert(_locus_tot[l]->get_chromosomePosition() != my_NAN); // the locus must be linked
		if (curChrom != _locus_tot[l]->get_chromosomePosition() && l) { // the chromosome has changed
			curChrom = _locus_tot[l]->get_chromosomePosition(); // get the current chromosome
			for (s = 0; s < nbSex; ++s) {
				_chromosomeSize[s][c] = l ? _locus_position_tot[s][l - 1] : 0; // set the chromosome length for each sex (cumulative)
				startChrom[s] = _chromosomeSize[s][c];
			}
			++c;							// increment the chromosome counter
		}
		for (s = 0; s < nbSex; ++s) {     // set the locus position for each sex
			_locus_position_tot[s][l] = startChrom[s] + _locus_tot[l]->get_locusPosition((sex_t) s);
		}
	}
	assert(c==_nb_chromosome-1 && l==_nb_locus_linked);
    
	// set the last chromosome size and the total length of the super chromosome
	for (s = 0; s < 2; ++s) { // has to be done always for both sexes (because of _chromosomeLength)!!
		_chromosomeLength[s] = _chromosomeSize[s][c] = _locus_position_tot[s][l-1];
	}
}

// -----------------------------------------------------------------------------------
// set_locus_positions
// ----------------------------------------------------------------------------------------
/** the locus position is added to the vector. It is evaluated if the order of the
 * loci corresponds to their genetic map.
 * this function is called by the protoTrait
 * _locus_pos is created
 */
void TGenomeProto::add_locus_position(const string& t, // type of trait (quanti/ntrl/global)
                                      const sex_t& s,          // sex (FEM/MAL)
                                      const unsigned int& c,   // chromosome
                                      const double& d,         // position on the chromosome in cM
                                      const unsigned int& l)   // locus index
{
	if (l < _locus_pos[t].size()) { // if this locus position is already present, then this means that the second sex is set
		TLocusPos* curPos = _locus_pos[t][l];
		if (curPos->chrom != c) error("Genome (%s): the sex specific genome does not match each other!\n", t.c_str()); // wrong chromosome
		if (l && _locus_pos[t][l - 1]->chrom == c // locus pos must increase within a chromosome
            && _locus_pos[t][l - 1]->pos[s] > d)
			error("Genome: the sex specific genome does not match each other!\n"); // wrong order of the loci
		curPos->pos[s] = d;	// set the new position
		_sex_specific_genome |= (curPos->pos[FEM] != curPos->pos[MAL]); // is genetic map different for females and males!
		assert((curPos->pos[FEM] != my_NAN && curPos->pos[MAL] != my_NAN) // both unlinked
               || (curPos->pos[FEM] == my_NAN && curPos->pos[MAL] == my_NAN));// or both linked
	}
	else {                              // add this locus position to the vector
		TLocusPos* curPos = new TLocusPos(c, d, d); // set the same postion for both sexes (it may be changed later...)
		if (l) {		// perform some checks (if it is not the first locus)
			TLocusPos* lastPos = *_locus_pos[t].rbegin();
			if (*lastPos > *curPos) error("Genome (%s): the loci positions have to be in increasing distance!\n", t.c_str()); // wrong order of the loci
		}
		_locus_pos[t].push_back(curPos);
	}
}

// -----------------------------------------------------------------------------------
// add_locus2genome
// ----------------------------------------------------------------------------------------
/** adds a locus to the genome. Function is called by TTraitProto 
 * t: trait type
 * i: locusID
 * l: TLocus
 */
void TGenomeProto::add_locus2genome(const string& t, const unsigned int& i, TLocus* l)
{
	++_nb_locus_tot;
	_locus_vector[t][i].push_back(l);
}


