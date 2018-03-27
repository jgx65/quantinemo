/** @file coal_deme.h
*
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


#ifndef coal_demeH
#define coal_demeH

#include <map>
#include <vector>
#include "genealogy.h"
//#include "functions.h"
//#include "simcomponent.h"
//#include "stathandler.h"
//#include "filehandler.h"
#include "types.h"

//#include "lifecycleevent.h"
using namespace std;


/** class to store for a GIVEN PATCH (and generation) the the populated patches and its immigrants
	* the population size is set after the immigrants have been specified!!!
	*
class dbPopTemp{
private:
	POP_SIZE                          _popSize;       // population size
	vector<pair<PATCH_ID, POP_SIZE> > _immigrants;    // vector of the immigrants (first: from id; secons: # immigrants)

public:
	dbPopTemp(){_popSize=0;}
	~dbPopTemp(){}

	inline void set_popSize(const POP_SIZE& size){assert(!_popSize);  _popSize = size;}
	inline void add_immigrants(const PATCH_ID& from, const POP_SIZE& size){
		_immigrants.push_back(pair<PATCH_ID, POP_SIZE>(from, size));
	}

	inline POP_SIZE&                           get_popSize()   {return _popSize;}
	inline vector<pair<PATCH_ID, POP_SIZE> >&  get_immigrants(){return _immigrants;}
};
*/

/** this class contains the same information as the class dbPop, but consumes less memory,
	* thus it is more efficient to use this class for the database. The drawback is that this
	* class cannot be directly feed by the information. One has to feed the class dbPop and
	* before changing to the next generation one has to dump the class dbPop to dbPopII.
	* In other words this is time consuming.
	*/


class dbPopII{
public:
	PATCH_ID    curID;           // current patch id
	POP_SIZE    popSize;         // population size

	PATCH_ID    nbNeighbour;     // number of neighbours (size of the array)
	MIGR_SIZE*  nbImmigrants;    // array of number of immigrants per neighbour
	PATCH_ID*   neighbourID;     // array of the neighbours-id

	dbPopII(): curID(0), popSize(0), nbNeighbour(0), nbImmigrants(0), neighbourID(0){}
	~dbPopII(){
		if(nbImmigrants) delete[] nbImmigrants;
		if(neighbourID)   delete[] neighbourID;
	}
};

/** class to store for a GIVEN PATCH (and generation) the the populated patches and its immigrants
	* the population size is set after the immigrants have been specified!!!
	*/
class dbPop{
private:
	POP_SIZE                          _popSize;       // population size
	vector<pair<PATCH_ID, POP_SIZE> > _immigrants;    // vector of the immigrants (first: from id; secons: # immigrants)

public:
	dbPop(){_popSize=0;}
	~dbPop(){}

	inline void set_popSize(const POP_SIZE& size){assert(!_popSize);  _popSize = size;}
	inline void add_immigrants(const PATCH_ID& from, const POP_SIZE& size){ _immigrants.push_back(pair<PATCH_ID, POP_SIZE>(from, size));}

	inline POP_SIZE&                           get_popSize()   {return _popSize;}
	inline vector<pair<PATCH_ID, POP_SIZE> >&  get_immigrants(){return _immigrants;}
};


#ifndef _DB_COALESCENCE
/** this class allows to store for a GIVEN GENERATION the pop sizes and immigrants. Note, that the
	* numbers in the db are diploid individuals, whereas the coalescence simulations
	* are haploid!
	*/
class dbCoalescence{
private:
	map<PATCH_ID, dbPop> _dbPatches;

public:
	dbCoalescence(){}
	~dbCoalescence(){}

	inline void add_immigrants(const unsigned int& to, const unsigned int& from, const unsigned int& size){
		_dbPatches[to].add_immigrants((PATCH_ID)from, (POP_SIZE) size);
	}

	inline void set_popSize(const unsigned int& pop, const unsigned int& size){
		_dbPatches[pop].set_popSize((POP_SIZE) size);
	}

	bool check_db();

	inline map<PATCH_ID, dbPop>&              get_patches(){return _dbPatches;}
	inline vector<pair<PATCH_ID, POP_SIZE> >& get_immigrants(const PATCH_ID& pop){return _dbPatches[pop].get_immigrants();}
	inline POP_SIZE&                          get_popSize(const PATCH_ID& pop) {return _dbPatches[pop].get_popSize();}
};
#endif

////////////////////////////////////////////////////////////////////////////////
//template <class <NODE> class LCE_Coalescence;
//class LCE_CoalescenceRecomb;
class TMetapop;
class LCE_Disperse;
class LCE_Coalescence;
class LCE_CoalescenceRecomb;

class TDeme {
protected:
	unsigned int _deme_size;     // in diploids
	unsigned int _lineages;      // in haploids (it is ChainNodeList.size())
	unsigned int _id;            // id of the corresponding patch

	vector<TTNode*> ChainNodeList;      // chained list of Nodes
	vector<TTNode*> ChainNodeList_temp; // temporal chained list of Nodes used for recombination purposes

public:
    TMetapop* _popPtr;
    
	TDeme() : _deme_size(0), _lineages(0), _id(my_NAN), _popPtr(0) {
	}
	TDeme(TMetapop* p, unsigned int id) : _deme_size(0), _lineages(0), _id(id), _popPtr(p) {
	}
	~TDeme(){
        if(!ChainNodeList.empty()){
            TTNode* ee = ChainNodeList.front();
        }
		assert(ChainNodeList.empty());
		assert(!_lineages);
	}

	inline void set_deme_size(unsigned int size) {assert(_id!=my_NAN);_deme_size = size;}
	inline void set_deme_size(unsigned int size, const unsigned int& id) {_deme_size = size; _id=id;}
	inline void set_id(unsigned int i) {_id=i;}
	inline void update_lineages(){_lineages=(unsigned int)ChainNodeList.size();}
	inline unsigned int update_get_lineages(){_lineages=(unsigned int)ChainNodeList.size(); return _lineages;}

	inline unsigned int get_deme_size() 	const	{return _deme_size;}
	inline unsigned int get_id()          const {return _id;}
	inline vector<TTNode*>& get_chainNodeList() {return ChainNodeList;}  // output may be modified!!!
    bool isNodeInChainNodeList(TTNode* pNode);
    vector<TTNode*>::iterator getNodePosInChainNodeList(TTNode* pNode);
	inline vector<TTNode*>& get_chainNodeList_temp() {return ChainNodeList_temp;}  // output may be modified!!!
	inline unsigned int merge_chainNodeLists(){
		if(!ChainNodeList_temp.empty()){
			ChainNodeList.insert(ChainNodeList.end(), ChainNodeList_temp.begin(), ChainNodeList_temp.end());
			ChainNodeList_temp.clear();
		}
		_lineages =(unsigned int) ChainNodeList.size();
		return _lineages;
	}
	inline unsigned int& get_lineages()         {return _lineages;}      // output may be modified!!!

	inline double get_probCoal(){assert(_deme_size); return ((double)_lineages*(_lineages-1))/(4*_deme_size);} // deme size is in diploid individuals
	inline double get_probCoal(unsigned int size){assert(_deme_size); return ((double)(size+_lineages)*(size+_lineages-1))/(4*_deme_size);} // deme size is in diploid individuals

	inline void clear_chainNodeList(){
		ChainNodeList.clear();
		_lineages=0;
	}

	void setNodes(const unsigned int& n);
	void set_sample_size(const unsigned int& size, TTNode** sampleNodes, unsigned int& nbNodes, unsigned int demeID);


	void migrate(const PATCH_ID& sinkID, unsigned int nbLineages, LCE_Coalescence& coal);
	void migrate(const PATCH_ID& sinkID, unsigned int nbLineages, LCE_Coalescence& coal, vector<TDeme*>& newDemes);
	void migrate(const PATCH_ID& sinkID, unsigned int nbLineages, LCE_Coalescence& coal, map<PATCH_ID, TDeme*>& demes, map<PATCH_ID, TDeme*>& newDemes);
	void migrate(TDeme* sink , unsigned int geneID, LCE_Coalescence& coal);
	void migrate_recomb(TDeme* sink , unsigned int geneID, LCE_Coalescence& coal);
	unsigned int post_migration_lineages();

	bool evacuate_lineages_disp_matrix(LCE_Coalescence* pCoal, bool all=false);
	bool evacuate_lineages_island(LCE_Coalescence* pCoal, bool all=false);
	bool evacuate_lineages_stepping_stone(LCE_Coalescence* pCoal, bool all=false);
	bool evacuate_stepping_stone(unsigned int lineages, unsigned int steps, LCE_Coalescence* pCoal);

	// used when all markers are unlinked
	unsigned int singleCoalescentEvent(const unsigned long& time, LCE_Coalescence& coal);
	unsigned int multipleCoalescentEvents(const unsigned long& time, LCE_Coalescence& coal);
	unsigned int mixedCoalescentEvents(const unsigned long& time, LCE_Coalescence& coal);
	TTNode* coalesce(TTNode* desc1, TTNode* desc2, const unsigned long& time, LCE_Coalescence& coal);

//	// used for the first tree for linked markers
//	unsigned int singleCoalescentEventII(const unsigned long& time, LCE_CoalescenceRecomb& coal);
//	unsigned int multipleCoalescentEventsII(const unsigned long& time, LCE_CoalescenceRecomb& coal);
//	unsigned int mixedCoalescentEventsII(const unsigned long& time, LCE_CoalescenceRecomb& coal);
//	//TTNode* coalesce(TTNode* desc1, TTNode* desc2, const long& time, LCE_Coalescence& coal);
//
//	// used for subsequent trees for linked markers (with multiple recombining nodes)
//	unsigned int singleCoalescentEventIII(const unsigned long& time, LCE_CoalescenceRecomb& coal, vector<TTNode*>& activeNodes);
//	unsigned int multipleCoalescentEventsIII(const unsigned long& time, LCE_CoalescenceRecomb& coal, vector<TTNode*>& activeNodes);
//	unsigned int mixedCoalescentEventsIII(const unsigned long& time, LCE_CoalescenceRecomb& coal, vector<TTNode*>& activeNodes);
//	TTNode* coalesce_recomb(TTNode* desc1, TTNode* desc2, const unsigned long& time, LCE_CoalescenceRecomb& coal, unsigned int demeID);
//	TTNode* coalesce_recomb1(TTNode* desc1, TTNode* desc2, const unsigned long& time, LCE_CoalescenceRecomb& coal);
//	TTNode* coalesce_recomb2(TTNode* desc1, TTNode* desc2, const unsigned long& time, LCE_CoalescenceRecomb& coal);
//
//	// used for subsequent trees for linked markers (with a single recombining node)
//	unsigned int singleCoalescentEventIV(const unsigned long& time, LCE_CoalescenceRecomb& coal, TTNode*& recombNode);
//	unsigned int multipleCoalescentEventsIV(const unsigned long& time, LCE_CoalescenceRecomb& coal, TTNode*& recombNode);
//	unsigned int mixedCoalescentEventsIV(const unsigned long& time, LCE_CoalescenceRecomb& coal, TTNode*& recombNode);

//	TTNode* recombine_back(TTNode* desc1, TTNode* desc2, const unsigned long& time, LCE_CoalescenceRecomb& coal);
//	TTNode* recombining_branch(TTNode* node, const unsigned long& time, LCE_CoalescenceRecomb& coal);

	bool ChainNodeList_remove(TTNode* node);
};



#endif /* COAL_DEME_H_ */
