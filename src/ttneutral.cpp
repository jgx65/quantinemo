/** @file ttneutral.cpp
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
 *   along with quantiNemo.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "ttneutral.h"
#include "tmetapop.h"
#include "stathandler.cpp"

// ----------------------------------------------------------------------------------------
// cstor
// ----------------------------------------------------------------------------------------
TTraitNeutralProto::TTraitNeutralProto() :
		_stats(0), _writer(0)
{
	_type = "ntrl";
	ini_paramset();
}

// ----------------------------------------------------------------------------------------
// cstor
// ----------------------------------------------------------------------------------------
TTraitNeutralProto::TTraitNeutralProto(int i) :
		_stats(0), _writer(0)
{
	_trait_index = i;
	_type = "ntrl";
	ini_paramset();
}

// ----------------------------------------------------------------------------------------
// cstor
// ----------------------------------------------------------------------------------------
TTraitNeutralProto::TTraitNeutralProto(const TTraitNeutralProto& T) :
		_stats(0), _writer(0)
{
	_copyTraitPrototypeParameters(T);
	ini_paramset();
}

// ----------------------------------------------------------------------------------------
// ini_paramset
// ----------------------------------------------------------------------------------------
void TTraitNeutralProto::ini_paramset()
{
	TTraitProto::ini_paramset();
	string trait = get_trait_indexStr_();
    add_parameter("ntrl_mutation_model" + trait, INT2, false, 0, 1, "0", false,
                  "The mutation model for neutral markers:\n" \
                  "  0: KAM (K-Allele mutation model) Mutation to any other allele.\n" \
                  "  1: SSM (Single step mutation model) Mutation to next bigger or smaller effect.",0);
}

// ----------------------------------------------------------------------------------------
// dstor
// ----------------------------------------------------------------------------------------
TTraitNeutralProto::~TTraitNeutralProto()
{
	if (_stats) delete _stats;
	if (get_trait_index() <= 1) {
		delete _writer; // writer will be deleted in testRepl
		_writer = NULL;
	}         //_writer is "static"!!!

	resetTotal();
}

// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
void TTraitNeutralProto::init(TMetapop* pMetapop)
{
	mut_model_t mutation_models[2] = { KAM, SSM };
	ini_base(pMetapop, mutation_models, 2);

	string trait = get_trait_indexStr_();

    assert(!_writer);

	// read the allelic file
	string allelicFile = get_parameter("ntrl_allelic_file" + trait)->get_arg(); // if empty parameter is not set
	if (!allelicFile.empty()) {
        allelicFile = get_popPtr()->get_iniFile_directory() + allelicFile;
		if (STRING::file_exists(allelicFile)) read_allele_file(allelicFile);
		else error("Allelic file '%s' (ntrl) cannot be found!\n",
				allelicFile.c_str());
	}

	//----------------------------------------------------------------------------
	// how to initialize alleles
	if (_initAlleleFreq) _ini_allele_model = INI_DIST;
	_protoGenome->set_ini_sequence_model(_aLocus, _nb_locus, _ini_allele_model); // following a distribution
}

// ----------------------------------------------------------------------------------------
// temporal_change
// ----------------------------------------------------------------------------------------
/** parameters which may change over time */
void TTraitNeutralProto::temporal_change(const unsigned int& gen)
{
	map<string, Param*>* pParam = _paramSet->getTemporalParams(gen);

	// if it is a temporal parameter
	if (pParam) {
		// check if a change has to be made
		map<string, Param*>* pMap = _paramSet->updateTemporalParams(gen);
		if (pMap) {
			TTraitProto::temporal_change(gen, pMap);      // base class changes?

			// iterate through the map and perform the updates
			map<string, Param*>::iterator pos = pMap->begin();
			for (; pos != pMap->end(); ++pos) {
				if (_trait_index <= 1 && pos->first == "ntrl_genot_logtime"
						&& _writer) { // the logtimes must only be changed for the first trait!!!
					_writer->set_GenerationOccurrence(
							(unsigned int) pos->second->get_value());
				}
			}
		}
	}
}

// ----------------------------------------------------------------------------------------
// get_info
// ----------------------------------------------------------------------------------------
string TTraitNeutralProto::get_info()
{
	string text;

	if (_trait_index) text = toStr(_trait_index) + ". neutral marker type: ";
	else text = "Neutral marker type: ";

	text += get_info_locus();

	return text;
}
// ----------------------------------------------------------------------------------------
// loadFileServices
// ----------------------------------------------------------------------------------------
void TTraitNeutralProto::loadFileServices(FileServices* loader)
{
	//writer
	unsigned int choice = (unsigned int) get_parameter_value("ntrl_save_genotype");
    assert(!_writer);
    if(choice) {
        _writer = (dynamic_cast<TTraitNeutralProto*>(_popPtr->getFirstPrototype(_type)))->_writer;
        if(!_writer){
            _writer = new TTNeutralFH();
			if (_writer->set(get_parameter("ntrl_genot_logtime"),
					get_parameter("ntrl_genot_dir")->get_arg(),
					get_parameter("ntrl_genot_filename")->get_arg(),
					get_parameter("ntrl_genot_script")->get_arg(),
					(int) get_parameter_value("ntrl_genot_sex"),
					(int) get_parameter_value("ntrl_genot_age"),
					choice, this, "genotype", ".dat", loader, get_popPtr())) {
				loader->attach(_writer);
			}
			else {
				delete _writer;
				_writer = NULL;
			}
		}
		else if (_writer) _writer->set(this); 		// append this trait
	}
}

// ----------------------------------------------------------------------------------------
// loadStatServices
// ----------------------------------------------------------------------------------------
void TTraitNeutralProto::loadStatServices(StatServices* loader)
{
	//allocate the stat handler
	if (_stats) delete _stats;
	_stats = new TTNeutralSH(this);
	loader->attach(_stats);
}

// ----------------------------------------------------------------------------------------
// clone
// ----------------------------------------------------------------------------------------
TTraitNeutral* TTraitNeutralProto::hatch()
{
	TTraitNeutral* new_trait = new TTraitNeutral();
	new_trait->set_from_prototype(this);
	return new_trait;
}
//----------------------------------------------------------------------------------------
// operator=
// ----------------------------------------------------------------------------------------
TTraitNeutral&
TTraitNeutral::operator=(const TTrait& T)
{
	error("");
	return *this;
}

//----------------------------------------------------------------------------------------
// operator==
// ----------------------------------------------------------------------------------------
void*
TTraitNeutral::get_allele(const unsigned int& loc, const unsigned int& all) const
{
	return (!(loc < pProto->_nb_locus) || !(all < ploidy) ?
			0 : (void*) &sequence[loc][all]);
}

//----------------------------------------------------------------------------------------
// operator==
// ----------------------------------------------------------------------------------------
bool TTraitNeutral::operator==(const TTrait& T)
{
	if (*pTraitProto != *T.pTraitProto) return false;

	const TTraitNeutral& TN = dynamic_cast<const TTraitNeutral&>(T);

	if (this != &TN || *pProto != *TN.pProto) return false;

	return true;
}

//----------------------------------------------------------------------------------------
// operator!=
// ----------------------------------------------------------------------------------------
bool TTraitNeutral::operator!=(const TTrait& T)
{
	if (!((*this) == T)) return true;
	else return false;
}

// ----------------------------------------------------------------------------------------
// Destructor
// ----------------------------------------------------------------------------------------
TTraitNeutral::~TTraitNeutral()
{
}

// ----------------------------------------------------------------------------------------
// set_from_prototype (used during hatching)
// ----------------------------------------------------------------------------------------
void TTraitNeutral::set_from_prototype(TTraitProto* T)
{
	pProto = dynamic_cast<TTraitNeutralProto*>(T);
	pTraitProto = T;
}

// ----------------------------------------------------------------------------------------
// show_up
// ----------------------------------------------------------------------------------------
void TTraitNeutral::show_up()
{
	message(
			"\n  Trait's type: ntrl\n\
       locus: %i\n\
     alleles: %i\n\
		sequence:",
			pProto->_nb_locus, pProto->get_nb_allele_max());

	for (unsigned int i = 0; (i < pProto->_nb_locus && i < 10); i++) {
		message("\n              %i %i", (int) sequence[i][0],
				(int) sequence[i][1]);
	}
	message("\n");
}

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                             ******** StatHandler ********/
// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
bool TTNeutralSH::init()
{
	StatHandler<TTNeutralSH>::init();
	return true;
}

// ----------------------------------------------------------------------------------------
// setStatsRecorders
// ----------------------------------------------------------------------------------------
/** this function invstigates if the given stat option "t" is present in this trait
 * example: n.adlt.nbLoc
 * - i:      "n"
 * - token:  "adlt"/"off"
 * - AGE:    ADULTS/OFFSPRG
 * - ageStr: "adult"/"offsprg"
 */
bool TTNeutralSH::setStatRecorders(const string& t)
{
	string type = _SHLinkedTrait->get_type();
	if (t[0] != type[0]) return false;         // wrong trait type   "q"?

	// is their a trait specification (ex. q1)?
	string::size_type pos = 1;                   // is the second position a '.'
	if (t[pos] != '.') { // trait index is specified (and this has to be an integer number)
		pos = t.find('.');
		unsigned int traitID;
		try {
			traitID = strTo<unsigned int>(t.substr(1, pos - 1));
		}
		catch(...) {
			error(
					"Stat 't' cannot be read: trait specification is not a number!\n",
					t.c_str());
		}
		if (!((!_trait_index && traitID == 1)       // a single trait is present
		|| (traitID == _trait_index))) return false; // wrong index when several traits are present
	}
	string token = t.substr(pos + 1);    // remove the "q.": eg. "adlt.fst_pair"
	string i = "n";

	// get the end
	pos = token.find('_');
	string end;                                        // "end" is empty
	if (pos != string::npos) {                           // there is an end
		end = token.substr(pos + 1);                       // eg. "pair"
		token.erase(pos + 1);                              // eg. "adlt.fst_"
	}

	// get the age
	pos = token.find('.');                        // find the second '.'
	if (pos == string::npos) return false;
	string ageToken = token.substr(0, pos);       // eg. "adlt"
	token.erase(0, pos + 1);                        // eg. "fst_"
	age_t AGE;
	string ageStr;
	if (ageToken == "adlt") {
		AGE = ADULTS;
		ageStr = "adult";
	}
	else if (ageToken == "off") {
		AGE = OFFSPRG;
		ageStr = "offsprg";
	}
	else return false;

	if (set_stat_coancestry(token, i, type, ageToken, end, AGE, ageStr)) return true; // adults and offspring
	if (set_stat_fstat(token, i, type, ageToken, end, AGE, ageStr)) return true; // adults and offspring
	if (set_stat_all_freq_local(token, i, type, ageToken, end, AGE, ageStr)) return true; // adults and offspring
	if (set_stat_all_freq_global(token, i, type, ageToken, end, AGE, ageStr))return true; // adults and offspring
	if (set_stat_locus_freq_local(token, i, type, ageToken, end, AGE, ageStr)) return true; // adults and offspring
	if (set_stat_locus_freq_global(token, i, type, ageToken, end, AGE, ageStr))return true; // adults and offspring
	if (set_stat_LD(token, i, type, ageToken, end, AGE, ageStr)) return true; // adults and offspring

	return false;
}
