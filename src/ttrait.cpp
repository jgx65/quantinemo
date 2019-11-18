/** @file ttrait.cpp
 *
 *   Copyright (C) 2006 Frederic Guillaume    <guillaum@zoology.ubc.ca>
 *   Copyright (C) 2008 Samuel Neuenschwander <samuel.neuenschwander@unil.ch>
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

#include "ttrait.h"
#include "patch.h"
#include "tmetapop.h"


//#include <algorithm>
using namespace std;

//------------------------------------------------------------------------------
TTraitProto::~TTraitProto(){
	resetTotal();
}

//------------------------------------------------------------------------------
/** same as the copy constructor */
// ----------------------------------------------------------------------------------------
void
TTraitProto::_copyTraitPrototypeParameters(const TTraitProto& T){
	error("Not implemented!\n");
}

//------------------------------------------------------------------------------
/** same as the copy constructor */
// ----------------------------------------------------------------------------------------
void
TTraitProto::resetTotal(){
	if(_initAlleleFreq){
		if(_nb_locus == 1 || _initAlleleFreq[0][0] != _initAlleleFreq[0][1]){   // each locus has its one array
			ARRAY::delete_3D(_initAlleleFreq, _initAlleleFreqCols, _nb_locus);
		}
		else{                                                        // all loci point to the first array
			delete[] _initAlleleFreq[0][0];
			ARRAY::delete_2D(_initAlleleFreq, _initAlleleFreqCols);
		}
	}

	if(_mutationFreq){
		if(_nb_locus == 1 || _mutationFreq[0] != _mutationFreq[1]){   // each locus has its one array
			ARRAY::delete_2D(_mutationFreq, _nb_locus);
		}
		else{                                                        // all loci point to the first array
			delete[] _mutationFreq[0];
			delete[] _mutationFreq;
			_mutationFreq = NULL;
		}
	}

	if(_locus_index) {delete[] _locus_index; _locus_index=NULL;}
	if(_nb_allele) {delete[] _nb_allele; _nb_allele=NULL;}
	if(_aLocus){delete[] _aLocus; _aLocus=NULL;}
}

// ----------------------------------------------------------------------------------------
// get_trait_index
// ----------------------------------------------------------------------------------------
double*
TTraitProto::get_ini_alleleFreq(TPatch* patch, const unsigned& locus){
	return _initAlleleFreq[patch->get_ID()%_initAlleleFreqCols][locus];
}

// ----------------------------------------------------------------------------------------
/** get_trait_index */
// ----------------------------------------------------------------------------------------
unsigned int
TTraitProto::get_trait_index() const {
	return _trait_index;
}

string
TTraitProto::get_trait_indexStr() const {
	if(!_trait_index) return "";
	return toStr(_trait_index);
}

string
TTraitProto::get_trait_indexStr_() const {
	if(!_trait_index) return "";
	return ("_" + toStr(_trait_index));
}

string
TTraitProto::get_trait_indexStr_t() const {
	if(!_trait_index) return "";
	return ("_t" + toStr(_trait_index));
}

// ----------------------------------------------------------------------------------------
// get_nb_allele_max
// ----------------------------------------------------------------------------------------
/** get the maximal number of possible alleles of the trait*/
unsigned int
TTraitProto::get_nb_allele_max()
{
	unsigned int size = _nb_allele[0];
	for(unsigned int i=1; i<_nb_locus; ++i){
		if(size<_nb_allele[i]) size = _nb_allele[i];
	}
	return size;
}

// ----------------------------------------------------------------------------------------
// set_ini_allelicArray
// ----------------------------------------------------------------------------------------
/** set up the initial frequency array/ get the dimensions (nothing else!)
 * set the parameter _initAlleleFreqCols and the coresponding array _initAlleleFreq
 * to find the dimension scroll through all values and look for arrays and get the largest array
 */
void
TTraitProto::set_ini_allelicArray (TMatrix* mat, const int& i)
{
	// check if the frequencies are identical for all patches
	map<unsigned int, map<unsigned int, string> >* pStrMatrix;
    map<unsigned int, map<unsigned int, string> >::iterator pos1;
    map<unsigned int, string>::iterator pos2;
    unsigned int curSize, nbLines=mat->getNbRows();
    unsigned int nbPatch = _popPtr->getPatchNbr();
    _initAlleleFreqCols = 1;
    for(unsigned int l=0; l<nbLines; ++l){
        if(mat->get(l, i) != my_NAN) continue;  // no matrix: all pops have the same init freqs

        // get the matrix and check for misreadings
		pStrMatrix = mat->get_strMatrix();
        if(!pStrMatrix) error("Could not read initial frequencies!\n");
        pos1=pStrMatrix->find(l);
        if(pos1==pStrMatrix->end()) error("Could not read initial frequencies!\n");
        pos2 = pos1->second.find(i);
        if(pos2==pos1->second.end()) error("Could not read initial frequencies!\n");
		
        //find the size of the matrix
		curSize = (unsigned int)STRING::strMatrix2vector<double>(pos2->second).size();
        if(nbPatch % curSize) warning("The number of initial allele frequencies (line %i) is not an entire subset of the number of patches!\n", l+1);
        
        if(curSize<_initAlleleFreqCols) continue;
        _initAlleleFreqCols=curSize;
	}
    
    // check the number of carrying capacities with the number of patches
    if(_initAlleleFreqCols>nbPatch) warning("There are more initial allele frequencies defined than patches! Only a part of the initial allele frequencies is considered!\n");

	ARRAY::create_3D(_initAlleleFreq, _initAlleleFreqCols, _nb_locus, _nb_allele, (double)my_NAN);
}

// ----------------------------------------------------------------------------------------
// set_allelic_valuesFreq
// ----------------------------------------------------------------------------------------
/** read the allelic values and their frequencies from a file (used for mutation model IMM and RMM)
 */
void
TTraitProto::read_allele_file (string filename)
{
#ifdef _DEBUG
	message(" TTraitProto::read_allele_file(%s)...\n", filename.c_str());
#endif
	// read the file
	TMatrix mat;
	map<string, int> fileInfo;
	bool hasLocus = true;
	unsigned int i, l, a, nbValues=0;

	// read the file
	fileInfo = mat.read_matrix(filename);
	unsigned int nbLines = mat.getNbRows();
	unsigned int nbCols  = mat.getNbCols();

	const unsigned int colsSize = 4;
	unsigned int cols[colsSize]        = {my_NAN,my_NAN,my_NAN,my_NAN};  // 0: locus; 1: allele; 2: mut_freq; 3: ini_freq;
	string colsName[colsSize] = {"col_locus","col_allele","col_mut_freq", "col_ini_freq"};
    

	// get col values form file info
	map<string, int>::iterator pos;
	for(i=0; i<colsSize; ++i){
		pos = fileInfo.find(colsName[i]);
		if(pos != fileInfo.end()){
			cols[i] = pos->second - 1;
			fileInfo.erase(pos);
		}
	}

	// output all not used columns
	for(pos = fileInfo.begin(); pos != fileInfo.end(); ++pos){
		warning("Allelic file for trait '%s' (%s): Column '%s' (%i) will not be considered!\n",
				get_type_index().c_str(), filename.c_str(), pos->first.c_str(), pos->second);
	}

	// check if all cols are available and also dimensions of the matrix are met and create the arrays
	if(cols[0] == my_NAN) hasLocus = false; // error("Allelic file (%s): Locus column is missing!", filename.c_str());
	if(cols[1] == my_NAN) error("Allelic file (%s): Allele column is missing!", filename.c_str());
	if(cols[2] != my_NAN){    // mutation frequency
		ARRAY::create_1D(_mutationFreq, _nb_locus);
		for(l=0; l<_nb_locus; ++l){
			_mutationFreq[l]= ARRAY::new_1D<double>(_nb_allele[l], (double)my_NAN);
		}
	}

	if(cols[3] != my_NAN){     // initial frequencies
		set_ini_allelicArray (&mat, cols[3]);
	}
    else _initAlleleFreqCols=0;

	for (i=0; i<colsSize; ++i){
		if(cols[i]!=my_NAN && cols[i] > nbCols) error("Allelic file (%s): File info is wrong: matrix has only %i columns and not %i!\n", filename.c_str(), nbCols, cols[i]);
	}

	// copy the allelic effects and their frequencies
	for(i=0; i<nbLines; ++i){
		if(hasLocus) l = (int)mat.get(i,cols[0])-1;    // if locus specific values
		else         l = 0;                            // if all loci have the same settings
		a = (int)mat.get(i,cols[1])-1;    // allele

		// test if out of range
		if(hasLocus && l>_nb_locus){
			error("Allelic values: Locus %i (allele %i) is out of range (only %i loci specified)!\n",
					l+1, a+1, _nb_locus);
		}
		if(a>_nb_allele[l]){
			error("Allelic values: Allele %i of locus %i is out of range (only %i alleles specified)!\n",
					a+1, l+1, _nb_allele[l]);
		}

		// check if the combination was already read
		if(   (_mutationFreq  && _mutationFreq[l][a] != my_NAN)
           || (_initAlleleFreq && _initAlleleFreq[0][l][a] != my_NAN)){
			error("Allelic values: Allelic value was already specified for locus %i and allele %i!\n", l+1, a+1);
		}

		// set the values
		if(hasLocus){         // locus specific settings
            set_allelicValues(&mat, i, l, a, cols);
			++nbValues;
		}
		else {                // all loci have the same settings
			for(l = 0; l < _nb_locus; ++l){
                set_allelicValues(&mat, i, l, a, cols);
				++nbValues;
			}
		}
	}

	// check if we have enough values
	unsigned int totValues = 0;
	for(l=0; l<_nb_locus; ++l){
		totValues += _nb_allele[l];
	}
	if(nbValues != totValues){
		error("Allelic file: %i of %i allelic values are not set by the file '%s'!\n",
				totValues-nbValues, totValues, filename.c_str());
	}

	// control the arrays and make them cumulative
	if(_mutationFreq){
		for(unsigned int l = 0; l < _nb_locus; ++l){
			ARRAY::make_frequency(_mutationFreq[l], _nb_allele[l]);  // adjust the sum to 1
			ARRAY::cumulative(_mutationFreq[l], _nb_allele[l]);      // make it cumulative
		}
	}

	if(_initAlleleFreqCols){
		unsigned int c, l;
        for(l = 0; l < _nb_locus; ++l){
            for(c = 0; c < _initAlleleFreqCols; ++c) {
				ARRAY::make_frequency(_initAlleleFreq[c][l], _nb_allele[l]);  // adjust the sum to 1
				ARRAY::cumulative(_initAlleleFreq[c][l], _nb_allele[l]);      // make it cumulative
			}
		}
	}
}

// ----------------------------------------------------------------------------------------
// set_mutation_model
// ----------------------------------------------------------------------------------------
void
TTraitProto::set_mutation_model(const string& trait, mut_model_t* model, const unsigned int& size){
	unsigned int index = (unsigned int)get_parameter_value(_type+"_mutation_model"+trait);
	if(index>=size) error("The mutation model is not available!\n");
	_mut_model = model[index];
	for(unsigned int l=0; l<_nb_locus; ++l){
		_aLocus[l].set_mutationModel(_mut_model);
	}
}

// ----------------------------------------------------------------------------------------
// set_mutation_rates
// ----------------------------------------------------------------------------------------
/** this function sets the mutation rates dependig on the parameters.
 */
void
TTraitProto::set_mutation_rates(const string& trait)
{
	_protoGenome->set_change_mutationRates(true);

	// if the single mutation rates are given as a matrix
	if(get_parameter(_type+"_mutation_rate"+trait)->is_matrix()){
		TMatrix* m = get_parameter(_type+"_mutation_rate"+trait)->get_matrix();
		unsigned int nb = m->get_dims(NULL);
		double* v = m->get();
		// check the number of mutation rates with the number of patches
		if(nb>_nb_locus) warning("There are more mutation rates than loci defined! Only a part of the mutation rates is considered!\n");
		else if(_nb_locus % nb) warning("The number of defined mutation rates is not an entire subset of the number of loci!\n");
		for(unsigned int l=0; l<_nb_locus; ++l){
			_aLocus[l].set_mutationRate(v[l % nb]);
		}
		delete m;
		return;
	}

	// the mutation rate is the same for all loci
	double mut_rate = get_parameter_value(_type+"_mutation_rate"+trait);
	for(unsigned int l=0; l<_nb_locus; ++l){
		_aLocus[l].set_mutationRate(mut_rate);
	}
}


// ----------------------------------------------------------------------------------------
// set_allelicValues
// ----------------------------------------------------------------------------------------
/** set the values for a single row 
 * 	cols = {"col_locus","col_allele","col_mut_freq","col_ini_freq","col_allelic_value"};
*/
void
TTraitProto::set_allelicValues(TMatrix* mat, const unsigned int& i, const unsigned int& l, const unsigned int& a, unsigned int* cols)
{
	if(_mutationFreq) {
        assert(cols[2]!=my_NAN);
        _mutationFreq[l][a]  = mat->get(i,cols[2]);   // mut  frequency
    }
	if(_initAlleleFreq){
        unsigned int curCol=cols[3];
        assert(curCol!=my_NAN);
        double value=mat->get(i,curCol);
        if(value!=my_NAN){
            for(unsigned int p=0; p<_initAlleleFreqCols; ++p){
                _initAlleleFreq[p][l][a] = value; // init frequencies per population
            }
        }
        else{
            // get the matrix and (check for misreadings have already been done)
            map<unsigned int, map<unsigned int, string> >* pStrMatrix = mat->get_strMatrix();
            assert(pStrMatrix);
            map<unsigned int, map<unsigned int, string> >::iterator pos1=pStrMatrix->find(i);
            assert(pos1!=pStrMatrix->end());
            map<unsigned int, string>::iterator pos2 = pos1->second.find(curCol);
            assert(pos2!=pos1->second.end());

            vector<double> vec=STRING::strMatrix2vector<double>(pos2->second);
            for(unsigned int p=0; p<_initAlleleFreqCols; ++p){
                _initAlleleFreq[p][l][a] = vec[p % vec.size()]; // init frequencies per population
            }
        }
    }
}

// ----------------------------------------------------------------------------------------
// set_locus_index
// ----------------------------------------------------------------------------------------
/** this function is called by TGenomeProto when the super chromsome is made
 * it sets the index of the locus corresponding to this trait
 */
void
TTraitProto::add_locus_index(const unsigned int& trait_index, const unsigned int& genome_index)
{
	_locus_index[trait_index] = genome_index;
}

// ----------------------------------------------------------------------------------------
// add_loci2genome
// ----------------------------------------------------------------------------------------
/** adds the loci to the genome depending on their locus index
 * in the settings file the index starts with 1 and within quantiNemo with 0
 */
void
TTraitProto::add_loci2genome(const string& trait)
{
	Param* p = get_parameter(_type+"_locus_index"+trait);
    unsigned int l = 0;
	if(p->isSet()){                                  // a matrix is passed
		TMatrix* m = p->get_matrix();
		double* v = m->get();
		unsigned int curIndex, size = m->length();
		if(size > _nb_locus){
			warning("Genome (%s): there are too many locus indexes specified. Only the first ones will be used!\n", _type.c_str());
			size = _nb_locus;
        }
		for(; l<size; ++l){
			curIndex = (unsigned int)v[l];
			_protoGenome->add_locus2genome(_type, (curIndex==my_NAN ? my_NAN : curIndex-1), &_aLocus[l]);
		}
		delete m;
	}

	// pass the not specified loci as unlinked
	for(;l<_nb_locus; ++l){
		_protoGenome->add_locus2genome(_type, my_NAN, &_aLocus[l]);
	}
}
// ----------------------------------------------------------------------------------------
//
// ----------------------------------------------------------------------------------------
/** to set the genetic map randomly one has to know the maximal locus index. This
 * function parses all proto traits of this type and investigates the locus indexes if present
 */
unsigned int
TTraitProto::get_locus_max_idx_of_trait_type()
{
	unsigned int max=0;
	map<string,TTraitProto*>::iterator curTrait = _popPtr->getTraitPrototypes().begin();
	map<string,TTraitProto*>::iterator endTrait = _popPtr->getTraitPrototypes().end();
	for(; curTrait != endTrait; ++curTrait){
		if(curTrait->second->get_type() != _type) continue;
		Param* curParam = curTrait->second->get_parameter(_type+"_locus_index"+curTrait->second->get_trait_indexStr_());
		if(!curParam->isSet()) continue;
		TMatrix* m = curParam->get_matrix();
		double* v = m->get();
		unsigned int size = m->get_dims(NULL);
		for(unsigned int l=0; l<size; ++l){
			if(max < (unsigned int)v[l]) max = (unsigned int)v[l];
		}
		delete m;
	}
	return max;
}

// ----------------------------------------------------------------------------------------
// get_info_locus
// ----------------------------------------------------------------------------------------
/** returns the info for locus and alleles */
string
TTraitProto::get_info_locus()
{
	string text =	toStr(_nb_locus) + " loci; ";
	unsigned int i, nbAllele=_nb_allele[0];
	for(i=1; i<_nb_locus; ++i){
		if(nbAllele!=_nb_allele[i]) break;
	}
	if(i==_nb_locus) text += toStr(nbAllele) + " alleles;";
	else             text += "variable number of alleles;";

	return text;
}

// ----------------------------------------------------------------------------------------
// temporal_change
// ----------------------------------------------------------------------------------------
/** parameters which may change over time */
void
TTraitProto::temporal_change(const unsigned int& gen, map<string, Param*>* pMap)
{
	assert(pMap);

	string trait = get_trait_indexStr_();

	// iterate through the map and performe the updates
	map<string, Param*>::iterator pos = pMap->begin();
	for(; pos != pMap->end(); ++pos){
		if(  pos->first == _type+"_mutation_rate"+trait) set_mutation_rates(trait);
	}
}

// ----------------------------------------------------------------------------------------
// ini_paramset
// ----------------------------------------------------------------------------------------
void
TTraitProto::ini_paramset()
{
	string trait = get_trait_indexStr_();
    string trait_temp = (_type=="quanti" ? "quantitative traits" : "neutral traits");

    set_paramset(_type+trait, trait_temp, false, _popPtr);

	add_parameter(_type+"_nb_trait",INT2,false,0,my_NAN,"1", false,
                  "Number of " + trait_temp + " to simulate.", 5);

    add_parameter(_type+"_loci"+trait,INT2,true,0,my_NAN,"0", false,
                  "Number of loci to simulate per " + trait_temp + ".", 5);
    
    add_parameter(_type+"_all"+trait,INT_MAT,false,1,255,"255", false,
                  "Maximum number of alleles per locus.", 5);
    

    add_parameter(_type+"_mutation_rate"+trait,DBL_MAT,false,0,1,"0", true,
                  "Mutation rate per locus and generation.", 5);
    
    add_parameter(_type+"_save_genotype",INT2,false,0,6,"0", false,
                  "Should the genotypes be outputted:\n" \
                  "  0: none (no output is generated)\n" \
                  "  1: FSTAT format (Goudet, 1995)\n" \
                  "  2: FSTAT extended format (same as point 1, but with " \
                  "additional individual information)\n" \
                  "  3: Arlequin format (Excoffier, 2010)\n" \
                  "  4: Arlequin format extended (same as point 3, but with " \
                  "additional individual information as comment)\n" \
                  "  5: PLINK format (Purcell, 2007)\n" \
                  "  6: PLINK format extended (same as point 5, but also " \
                  "multi-allele loci are listed).",3);
    
    add_parameter(_type+"_genot_sex",INT2,false,0,2,"0", false,
                  "Which sex(es) should the genotype output contain:\n" \
                  "  0: both (output includes both sexes)\n" \
                  "  1: females (output includes only female genotypes)\n" \
                  "  2: males (output includes only male genotypes).",3);
    
    add_parameter(_type+"_genot_age",INT2,false,0,2,"0", false,
                  "Which ag(es) should the genotype output contain:\n" \
                  "  0: adults (output includes only adult genotypes)\n" \
                  "  1: juveniles (output includes only juvenile genotypes)\n" \
                  "  2: both (output includes juvenile and adult genotypes)",3);
    
    add_parameter(_type+"_genot_dir",STR,false,0,0,"", false,
                  "The directory name where genotypes are stored.",3);
    
    add_parameter(_type+"_genot_filename",STR,false,0,0,"", false,
                  "The base name for the genotype output files.",3);
    
    add_parameter(_type+"_genot_logtime",INT2,false,0,my_NAN ,"1",true,
                  "The time interval of the genotype output.",3);
    
    add_parameter(_type+"_genot_script",STR,false,my_NAN,my_NAN,"", false,
                  "The script which will be launched after the genotype file is generated.",3);
    

	// genetic map (identical for traits of the given type)
    add_parameter(_type+"_genome"           ,MAT,false,my_NAN,my_NAN,"", false,
                  "The genetic map in cM for " + trait_temp + ".", 3);
    
    add_parameter(_type+"_genome_fem"       ,MAT,false,my_NAN,my_NAN,"", false,
                  "The genetic map for females in cM for " + trait_temp + ". " \
                  "Note, the order must be identical between male and female map.", 3);
    
    add_parameter(_type+"_genome_mal"       ,MAT,false,my_NAN,my_NAN,"", false,
                  "The genetic map for males in cM for " + trait_temp + ". " \
                  "Note, the order must be identical between male and female map.", 3);
    
    add_parameter(_type+"_locus_index"+trait,MAT,false,my_NAN,my_NAN,"", false,
                  "Assignment of the loci to a genetic map position.", 3); // for each trait

	// genotype file as input?
    add_parameter(_type+"_ini_genotypes",STR,false,my_NAN,my_NAN,"", false,
                  "Filename containing the initial genotypes in FSTAT format.", 3);
    
    
    add_parameter(_type+"_allelic_file" + trait, STR, false, my_NAN, my_NAN, "", false,
                  "Filename containing allelic information, such as " +
                  (string)(_type=="quanti" ? "the allelic effects, " : "") +
                  "the mutation frequency, and/or the initial frequency. " \
                  "Note, no '%' is assumed in front of the filename.",1);
    
    add_parameter(_type+"_ini_allele_model" + trait, INT2, false, 0, 1, "0", false,
                  "How should the initial genotypes be initialized:\n" \
                  "  0: polymorphic (alleles are randomly drawn (uniform distribution)\n" \
                  "  1: monomorphic (all individuals are fixed for a single allele with " \
                  "index " + _type + "_all/2'." ,1);
    

}

// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
void
TTraitProto::ini_base(TMetapop* pMetapop, mut_model_t* mutation_models, unsigned int nbModels)
{
	_popPtr = pMetapop;
	_protoGenome = _popPtr->get_protoGenome();
	string trait = get_trait_indexStr_();

	if(get_trait_index() <= 1){    	// only once for all types of this trait together
		_protoGenome->set_genome_positions(_type, this);  	// set the positions of all loci of this type of trait
		_ini_fstat_file = get_parameter(_type+"_ini_genotypes")->get_arg();
	}

	// get the parameters (they are validated if (and only if) they are used)
	_nb_locus           = (unsigned int)get_parameter_value(_type+"_loci"+trait);
	_ini_allele_model   = (ini_model_t) get_parameter_value(_type+"_ini_allele_model"+trait);   //0: uniform, 1: monomorph


	//----------------------------------------------------------------------------
	// set the number of alleles
	Param* pParam = get_parameter(_type+"_all"+trait);
	if(_nb_allele) delete _nb_allele;
	_nb_allele = new unsigned int[_nb_locus];
	if(pParam->is_matrix()){
		TMatrix* m = pParam->get_matrix();
		unsigned int size = m->get_dims(NULL);
		double* array = m->get();       // a matrix is always double
		// check the matrix sizes with the number of patches
		if(size>_nb_locus){
			string name = _type+"_all"+trait;
			warning("Parameter %s: The matrix length exceeds the number of loci! Not all values will be used!\n", name.c_str());
		}
		else if(_nb_locus % size){
			string name = _type+"_all"+trait;
			warning("Parameter %s: The matrix length is not an entire subset of the number of loci!\n", name.c_str());
		}
		for(unsigned int l=0; l<_nb_locus; ++l){
			_nb_allele[l] = (unsigned int)array[l % size];
		}
		delete m;
	}
	else{
		unsigned int nbAllele = (unsigned int)pParam->get_value();
		for(unsigned int l=0; l<_nb_locus; ++l){
			_nb_allele[l] = nbAllele;
		}
	}
	//----------------------------------------------------------------------------
	// create the loci
	if(_aLocus) delete[] _aLocus;
	_aLocus = new TLocus[_nb_locus];
	for(unsigned int l=0; l<_nb_locus; ++l){
		_aLocus[l].set_nb_allele(_nb_allele[l]);
		_aLocus[l].set_traitPointer(this);
		_aLocus[l].set_locus_id_trait(l);
	}

	//----------------------------------------------------------------------------
	// genetic map
	ARRAY::create_1D(_locus_index, _nb_locus, (unsigned int)my_NAN);
	add_loci2genome(trait);                        // set the positions of the loci

	//----------------------------------------------------------------------------
	// set mutation function
	set_mutation_model(trait, mutation_models, nbModels);
	set_mutation_rates(trait);
}


//------------------------------------------------------------------------------
/** for the derived classes to compare directly the paramters in this base calss */
TTraitProto&
TTraitProto::setTTraitProto(const TTraitProto& T){             // for operator=
	_copyTraitPrototypeParameters(T);
	return *this;
}

//------------------------------------------------------------------------------
bool
TTraitProto::isEqualTTraitProto(const TTraitProto& T){         // for operator==
	if(_type              != T._type              )  return false;
	if(_nb_locus          != T._nb_locus          )  return false;
	if(_absolute_index    != T._absolute_index    )  return false;
	if(_trait_index       != T._trait_index       )  return false;

	for(unsigned int i=0; i<_nb_locus; ++i){
		if(_nb_allele[i] != T._nb_allele[i])  return false;
	}

	return true;
}

//------------------------------------------------------------------------------
bool
TTraitProto::isUnequalTTraitProto(const TTraitProto& T){       // for operator!=
	return !(*this == T);
}

//------------------------------------------------------------------------------
/** destructor */
TTrait::~TTrait ( ){
	if(sequence) delete[] sequence;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void*
TTrait::get_allele(const unsigned int& loc, const unsigned int& all)  const  {
	if(loc<pTraitProto->_nb_locus && all<pTraitProto->_nb_allele[loc]) return (void*)&sequence[loc][all];
	return 0;
}

//------------------------------------------------------------------------------
/** same as the copy constructor */
void
TTrait::_copyTTraitParameters(const TTrait& T)
{
	assert(!sequence);
	pTraitProto   = T.pTraitProto;
}

// ----------------------------------------------------------------------------------------
// ini_sequence
// ----------------------------------------------------------------------------------------
/** initialization of the sequence */
void
TTrait::ini_sequence (TPatch* patch)
{
	for(unsigned int l=0; l<pTraitProto->_nb_locus; ++l){
		pTraitProto->_aLocus[l].ini_sequence(sequence[l], patch);
	}
} // ini_sequence

// ----------------------------------------------------------------------------------------
// ini_base
// ----------------------------------------------------------------------------------------
/** set pointers of the trait sequence to the genome sequence */
void
TTrait::ini(TIndividual* ind)
{
	assert(ind->genome.get_sequence());
	if(sequence){delete[] sequence; sequence=NULL;}

	unsigned int nbLocus = 	pTraitProto->_nb_locus;
	TLocus* aLocus       = pTraitProto->_aLocus;
	sequence             = new ALLELE*[nbLocus];
	ALLELE** cur_seq = ind->genome.get_sequence();
	for(unsigned int l=0; l<nbLocus; ++l){
		sequence[l] = cur_seq[aLocus[l].get_locus_id_tot()];
	}
}





