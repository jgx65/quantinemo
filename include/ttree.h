/** @file ttree.h
*
*   Copyright (C) 2008 Samuel Neuenschwander <samuel.neuenschwander@unil.ch>
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
//---------------------------------------------------------------------------

#ifndef ttreeH
#define ttreeH

#include <assert.h>
#include <map>
#include <sstream>
using namespace std;
//---------------------------------------------------------------------------
/** the tree implemented is a n-tree where all leaves have the same depth.
	* The implementation contains two classes: the TTree and the TNode class.
	* The TTree class is derived form the TNode class and contains some additional
	* functions to deal with the entire tree. In theory it would be feasible to have
	* a single class, however I prefer for clarity reasons to have a tree (i.e. root) and nodes.
	*
	* The tree is implementation as a template where the key and the element may be any
	* data type. However, if the tree is read by a stream (i.e. no keys are specified)
	* the default is to use a unsigned int as key data type.
	*/
//---------------------------------------------------------------------------
/** template K: key for the tree
	* template E: data type of the element to store
	*/

const unsigned int NOT_SET=99999;

/** input and output declarations */
template <class K, class E> class TNode;
template <class K, class E> class TTree;
template <class K, class E> ostream &operator<< (ostream &os, const TNode<K,E>& o);
template <class K, class E> istream &operator>> (istream &is, TNode<K,E>& o);

/******************************************************************************/
// class TNode /////////////////////////////////////////////////////////////////
template <class K=unsigned int, class E=double>
class TNode{
protected:
	map<K, TNode<K,E>* >           _children;  // container for the childern
	TNode<K,E>*                    _pParent;   // pointer to the parent node
	K                              _key;       // key of the node (repeat of the _childern container)
	E*					 		               _value;     // value of the node if it is a leave
	unsigned int                   _depth;     // depth of the current node (0: leave; >0: branch)

public:
	// input output functions
	friend ostream &operator<< <>(ostream &os, const TNode<K,E>& o);
	friend istream &operator>> <>(istream &is, TNode<K,E>& o);

	// cunstructors
	TNode();                                          // default root constructor
	TNode(unsigned int size);                         // root constructor with known depth
	TNode(unsigned int size, E def);                  // root constructor with default value and known depth
	TNode(K key, TNode<K,E>* p);                      // constructor for a branch node
	TNode(K key, TNode<K,E>* p, E elem);              // constructor for a leave
	TNode(K* key, E val, K cur_key, TNode<K,E>* p);   // recursive constructor until the leave

	// destructor
	~TNode(){reset();}
	void reset();

	// setter
	void        set(K* key, E val);  // set the value recursively
	void        set(E val);          // it is a leave set the value
	TNode<K,E>* add(K key);          // add node at the specified key
	TNode<K,E>* push_back();	       // add a new node (may be also a leave)
	TNode<K,E>* push_back(E val);	   // add a new leave with element val

	// getter
	TNode<K,E>* get(K key);          // get the node with key K (if not present create a new node with key K)
	TNode<K,E>* get_static(K key);   // get the node with key K and NULL if not present
	TNode<K,E>* first();             // get the first node
	TNode<K,E>* last();              // get the after last node
	TNode<K,E>* get_parent(){return _pParent;}        // get parent node
	TNode<K,E>& operator[] (K i);
	E* get_static(K* key);           // when a branch -> get the value recursive (and NULL if not present)
	E& get(K* key);                  // when a branch -> get the value recursive (and create one if not present)
	E& get(){return *_value;}        // when a leave  -> get the value
	K  get_key(){return _key;}       // get the key

	E* next(K* key, bool first=false); // get the next element and change the key accordingly
	E* first(K* key){return next(key, true);} // get the first element and change the key accordingly

	void        change_key(K newKey);

	unsigned int size() {return _children.size();}
	bool         empty(){return _children.empty();}


	// default value
	void set_default(E def);
	E get_default();

	bool is_root(TNode<K,E>* p){return (bool)p->_pParent;}         	// the node in question
	bool is_root(){return (bool)_pParent;}                         	// this node
	TNode<K,E>* get_root(){return _pParent ? _pParent->get_root() : this;}    // get the node recursively

	bool is_leave(TNode<K,E>* p){return p->_pParent && p->_value;}	// the node in question
	bool is_leave(){return _pParent && _value;}                     // the node in question

	void set_depth_backward();
	unsigned int get_depth(){return _depth;}

	void print_tabbed(ostream &os, unsigned int tab=0) const;
	void print_bracket(ostream& os, bool first=true, const K* index = NULL) const;

	void read(string& ss);
	void read(istream &is);

	// functions to retrieve the keys, nodes, or elments (only for the lowest level)  as an array
	void get_keys_as_array(K* &aKey, unsigned int& size);
	void get_nodes_as_array(TNode<K,E>** &aNodes, unsigned int& size);
	void get_nodes_as_array_static(TNode<K,E>** &aNodes, unsigned int& size);
	void get_elements_as_array(E* &array, unsigned int& size);
};

/******************************************************************************/
// class TTree /////////////////////////////////////////////////////////////////
template <class K=unsigned int, class E=double>
class TTree : public TNode<K,E>{
public:
	TTree():TNode<K,E>(){}                                                      // default constructor
	TTree(unsigned int size): TNode<K,E>(size){}
	TTree(unsigned int size, E default_value): TNode<K,E>(size, default_value){}
	TTree(K* key, E val): TNode<K,E>(key, val, this){}  // constructor together with the first element
	TTree(string text){read(text);}                     // create a tree based on the string
	TTree(istream& text){read(text);}                   // create a tree based on the stream

	~TTree(){}                                          // destructor

	// setter
	void set(K* key, const unsigned int& size, E val){check_depth(size); TNode<K,E>::set(key, val);}

	// getter
	E&   get_i(K* key, const unsigned int& size){check_depth(size); return TNode<K,E>::get(key);}
	E*   get_static_i(K* key, const unsigned int& size) const;
	E*   next_i(K* key, const unsigned int& size){check_depth(size); return TNode<K,E>::next(key);}
	E*   first_i(K* key, const unsigned int& size){check_depth(size); return TNode<K,E>::first(key);}

	// read tree
	void read(string& ss);      // read an entire tree
	void read(istream &is);     // read an entire tree

	void check_depth(unsigned int size) const;
};

//---------------------------------------------------------------------------
template <class K, class E>
E* TTree<K,E>::get_static_i(K* key, const unsigned int& size) const {
	check_depth(size);
	return TNode<K,E>::get_static(key);
}

//---------------------------------------------------------------------------
/** input output functions */
template <class K, class E>
ostream &operator<< (ostream &os, const TNode<K,E>& o){
	o.print_bracket(os, 0);
	return os;
}

template <class K, class E>
istream &operator>> (istream &is, TNode<K,E>& o){
	o.read(is);
	return is;
}

//---------------------------------------------------------------------------
/** cunstructors */
template <class K, class E>
TNode<K,E>::TNode()                        // default root constructor
{
	_pParent = NULL;
	_value   = NULL;
	_depth   = NOT_SET;
}
template <class K, class E>
TNode<K,E>::TNode(unsigned int size)        // root constructor with known depth
{
	_pParent = NULL;
	_value   = NULL;
	_depth   = size;
}
template <class K, class E>
TNode<K,E>::TNode(unsigned int size, E def) // root constructor with default value and known depth
{
	_pParent = NULL;
	_value   = new E[1];
	*_value  = def;
	_depth   = size;
}
template <class K, class E>
TNode<K,E>::TNode(K key, TNode<K,E>* p)                // constructor for a branch node
{
	_pParent = p;
	_value   = NULL;
	_key     = key;
	p->_depth == NOT_SET ? _depth = NOT_SET : _depth = p->_depth-1;
}
template <class K, class E>
TNode<K,E>::TNode(K key, TNode<K,E>* p, E elem)        // constructor for a leave
{
	_pParent = p;
	_key     = key;
	set(elem);
}
template <class K, class E>
TNode<K,E>::TNode(K* key, E val, K cur_key, TNode<K,E>* p)  // recursive constructor until the leave
{
	_pParent = p;
	_value   = NULL;
	_depth   = p->_depth-1;
	_key     = cur_key;
	if(!_depth) set(val);         // it is a leave
	else        set(key, val);
}


//---------------------------------------------------------------------------
/** set the value */
template <class K, class E>
void TNode<K,E>::set(E val){                 // it is a leave set the value
	_value   = new E[1];
	*_value  = val;
	_depth   = 0;
	set_depth_backward();          // set the depth of the tree
}

//---------------------------------------------------------------------------
/** getter and setter of the default value */
template <class K, class E>
void TNode<K,E>:: set_default(E def){
	if(_pParent) return _pParent->set_default(def);    // search recusively the root
	if(!_value) _value = new E[1];
	*_value = def;
}

template <class K, class E>
E TNode<K,E>:: get_default(){
	if(_pParent) return _pParent->get_default();       // search recusively the root
	if(!_value) throw("The default value has not yet been specified!");
	return *_value;
}

//---------------------------------------------------------------------------
/** resets the node and all lower nodes */
template <class K, class E>
void TNode<K,E>::reset()
{
	typename map<K, TNode<K,E>* >::iterator _pos = _children.begin();
	for(; _pos != _children.end(); ++_pos){
		if(_pos->second) delete _pos->second;
	}
	_children.clear();

	if(_value) {delete[] _value; _value=NULL;}
	_depth = NOT_SET;
	_pParent = NULL;
}

//---------------------------------------------------------------------------
/** set the value recursively following the array */
template <class K, class E>
void TNode<K,E>::set(K* key, E val)
{
	assert(_depth != NOT_SET);
	if(_depth){           							 // it is a branch
		typename map<K, TNode<K,E>* >::iterator _pos = _children.find(*key);
		if(_pos != _children.end()){       // child is already present
			_pos->second.set(key+1, val);
		}
		else{														   // child has to be created
			_children.insert(pair<K, TNode<K,E>* >(*key, new TNode<K,E>(key+1, val, *key, this)));
		}
	}
	else set(val);										   // it is a leave
}

//---------------------------------------------------------------------------
template <class K, class E>
E* TNode<K,E>::get_static(K* key)
{
	if(_depth){           							 // it is a branch
		typename map<K, TNode<K,E>* >::iterator pos = _children.find(*key);
		if(pos != _children.end()) return pos->second->get(key+1); // child is already present
		else                       return NULL;                    // not present
	}
	else{															 // it is a leave
		if(_value) return &_value;
		else       return NULL;          // not present (should a error be returned?)
	}
}

//---------------------------------------------------------------------------
/** returns a reference to the node given by the key */
template <class K, class E>
TNode<K,E>& TNode<K,E>::operator[] (K key)
{
	assert(_depth);
	typename map<K, TNode<K,E>* >::iterator _pos = _children.find(key);
	if(_pos == _children.end()){
		_pos = _children.insert(pair<K, TNode<K,E>* >(key, new TNode<K,E>(key, this))).first;
	}
	return *_pos->second;
}

//---------------------------------------------------------------------------
/** returns a reference to an element pointed to by a key array */
template <class K, class E>
E& TNode<K,E>::get(K* key)
{
	if(_depth){           						 // it is a branch
		typename map<K, TNode<K,E>* >::iterator _pos = _children.find(*key);
		if(_pos == _children.end()){     // not yet present
			_pos = _children.insert(pair<K, TNode<K,E>* >(*key, new TNode<K,E>(*key, this))).first;
		}
		return _pos->second->get(key+1); // recursive
	}
	else{															 // it is a leave
		if(!_value){                     // not yet present (initialize it)
			_value = new E[1];
			*_value = get_default();
		}
		return *_value;
	}
}

//---------------------------------------------------------------------------
/** returns a pointer to the next element and changes the key accordingly.
	* returns NULL if there are no more elements
	* if first is true, take the first child, without considering the key (first is by default false)
	*/
template <class K, class E>
E* TNode<K,E>::next(K* key, bool first)
{
	typename map<K, TNode<K,E>* >::iterator pos;

	// if the first node has to be taken (ignoring the key)
	if(first){
		pos = _children.begin();                 // get the first child
		*key = pos->first;                       // adapt the key
		if(_depth==1) return &pos->second->get();// get the first element
		return pos->second->next(++key, true);   // go on recursively with the first element
	}

	// get the child and go one with the search
	pos = _children.find(*key);                // get the current child node depending on the key
	assert(pos != _children.end());            // this child must be available

	// go on with recursive search
	if(_depth==1){
		++pos;                                   // get the next element
		if(pos == _children.end()) return NULL;  // that was the last key => return NULL
		*key = pos->first;                       // adapt the key
		return &pos->second->get();              // return the element
	}
	E* val = pos->second->next(key+1);         // recursive search
	if(val) return val;                        // Yeap, we have the element => return the value

	// the specified child has no morde elments => go to the next node and retry
	++pos;
	if(pos == _children.end()) return NULL;    // that was the last key => return NULL

	// there is a next child
	*key = pos->first;                         // adapt the key
	return pos->second->next(++key, true);     // go on recursively for the first element
}

//---------------------------------------------------------------------------
/** add an leave node to the end of the children */
template <class K, class E>
TNode<K,E>* TNode<K,E>::push_back(E val)
{
	assert(_depth == NOT_SET || _depth == 1);
	TNode<K,E>* newNode;
	if(_children.empty()){         // it is the first element of this node
		newNode = new TNode<K,E>((K)0, this, val);
		_children.insert(pair<K, TNode<K,E>* >((K)0, newNode));
	}
	else{
		K key = _children.rbegin()->first;
		newNode = new TNode<K,E>(++key, this, val);
		_children.insert(pair<K, TNode<K,E>* >(key, newNode));
	}
	return newNode;
}

//---------------------------------------------------------------------------
/** add a node to the end of the children (this may also be a leave)*/
template <class K, class E>
TNode<K,E>* TNode<K,E>::push_back()
{
	K key;
	if(_children.empty()) key = (K)0;                         // it is the first element of this node
	else{key = _children.rbegin()->first; ++key;}             // interate the last key

	TNode<K,E>* newNode = new TNode<K,E>(key, this);          // create the new node
	_children.insert(pair<K, TNode<K,E>* >(key, newNode));    // inser the node
	return newNode;
}

//---------------------------------------------------------------------------
/** add a node at key K and return the node (if the key is already present an error is returned) */
template <class K, class E>
TNode<K,E>* TNode<K,E>::add(K key)
{
	assert(_depth == NOT_SET || _depth);    // it can't be a leave
	TNode<K,E>*  newNode = new TNode<K,E>(key, this);
	typename map<K, TNode<K,E>* >::iterator pos = _children.find(key);
	if(pos != _children.end()) throw("A node of this key is already present!");
	_children.insert(pair<K, TNode<K,E>* >(key, newNode));
	return newNode;
}

//---------------------------------------------------------------------------
// get the node with key K and NULL if not present
template <class K, class E>
TNode<K,E>* TNode<K,E>::get_static(K key)
{
	typename map<K, TNode<K,E>* >::const_iterator pos = _children.find(key);
	if(pos != _children.end()) return pos->second;
	return NULL;
}

//---------------------------------------------------------------------------
// get the node with key K  and create a new node with key K if not present
template <class K, class E>
TNode<K,E>* TNode<K,E>::get(K key)
{
	typename map<K, TNode<K,E>* >::iterator pos = _children.find(key);
	if(pos != _children.end()) return pos->second;

	// it is not present create it
	TNode<K,E> newNode = new TNode<K,E>(key, this);
	_children.insert(pair<K, TNode<K,E>* >(key, newNode));
	return newNode;
}

//---------------------------------------------------------------------------
// get the first node and NULL if not present
template <class K, class E>
TNode<K,E>* TNode<K,E>::first()
{
	if(_children.empty()) return NULL;
	return _children.begin()->second;
}

//---------------------------------------------------------------------------
// get the last node and NULL if not prsent
template <class K, class E>
TNode<K,E>* TNode<K,E>::last()
{
	if(_children.empty()) return NULL;
	return _children.rbegin()->second;
}

//---------------------------------------------------------------------------
/** an array is returned containing the keys of the children.
	* If the array is not yet generated it is done or adjuster to the number of nodes
	*/
template <class K, class E>
void TNode<K,E>::get_keys_as_array(K* &aKeys, unsigned int& size)
{
	if(size != _children.size()){                // adjust the arrays if needed
		size = _children.size();
		if(aKeys) delete[] aKeys;
		aKeys = size ? new K[size] : NULL;
	}
	typename map<K, TNode<K,E>* >::iterator pos = _children.begin();
	for(unsigned int i=0; pos != _children.end(); ++pos, ++i){
		aKeys[i]  = pos->first;
	}
}

//---------------------------------------------------------------------------
/** an array is returned containing the pointers to the nodes of the children.
	* If the array is not yet generated it is done or adjuster to the number of nodes
	*/
template <class K, class E>
void TNode<K,E>::get_nodes_as_array(TNode<K,E>** &aNodes, unsigned int& size)
{
	if(size != _children.size()){                // adjust the arrays if needed
		size = _children.size();
		if(aNodes) delete[] aNodes;
		aNodes = size ? new TNode<K,E>*[size] : NULL;
	}
	typename map<K, TNode<K,E>* >::iterator pos = _children.begin();
	for(unsigned int i=0; pos != _children.end(); ++pos, ++i){
		aNodes[i] = pos->second;
	}
}

//---------------------------------------------------------------------------
/** an array is returned containing the pointers to the nodes of the children.
	* in contrast to the previous function this array is adjusted to missing children,
	* i.e. a missing child will be visible as a NULL-pointer. The key thus reflects the position in the array.
	*/
template <class K, class E>
void TNode<K,E>::get_nodes_as_array_static(TNode<K,E>** &aNodes, unsigned int& size)
{
	if(size != _children.rbegin()->first+1){                // adjust the arrays if needed
		size = _children.rbegin()->first+1;
		if(aNodes) delete[] aNodes;
		aNodes = size ? new TNode<K,E>*[size] : NULL;
	}
	typename map<K, TNode<K,E>* >::iterator pos = _children.begin();
	for(unsigned int i=0; i<size; ++i){
		if(i!=pos->first) aNodes[i] = NULL;            // the element is jumped
		else{
			 aNodes[i] = pos->second;
    	++pos;
		}
	}
}
//---------------------------------------------------------------------------
/** an array is returned containing a copy of the elements of the children.
	* If the array is not yet generated it is done or adjuster to the number of nodes
	*/
template <class K, class E>
void TNode<K,E>::get_elements_as_array(E* &array, unsigned int& size)
{
	if(_depth!=1) throw("Elements may only be retrieved from leaves!");
	if(size != (unsigned int)_children.size()){                // adjust the arrays if needed
		size = (unsigned int)_children.size();
		if(array) delete[] array;
		array = size ? new E[size] : NULL;
	}
	typename map<K, TNode<K,E>* >::iterator pos = _children.begin();
	for(unsigned int i=0; pos != _children.end(); ++pos, ++i){
		array[i] = *pos->second->_value;
	}
}

//---------------------------------------------------------------------------
// change the key of this node
template <class K, class E>
void TNode<K,E>::change_key(K newKey)
{
	// find the corresponding node in the parent
	typename map<K, TNode<K,E>* >::iterator pos = _pParent->_children.find(_key);
	_pParent->_children.erase(pos);     // erase the entry in the children container of the parents
	_pParent->_children.insert(pair<K, TNode<K,E>* >(newKey, this));    // insert the new key
	_key = newKey;
}

//---------------------------------------------------------------------------
/* print the tree in a tabbed version*/
template <class K, class E>
void TNode<K,E>::print_tabbed(ostream& os, unsigned int tab) const
{
	if(!_depth) os  << " " << *_value;         // it is a leave
	else{
		typename map<K, TNode<K,E>* >::const_iterator pos =_children.begin();
		for(; pos != _children.end(); ++pos){
			os << "\n";
			os.width(tab);
			os << "";     // insert the tab empty characters
			os << pos->first << ":";
			pos->second->print_tabbed(os, tab+1);
		}
	}
}

//---------------------------------------------------------------------------
/* print the tree in a bracketed version*/
template <class K, class E>
void TNode<K,E>::print_bracket(ostream& os, bool first, const K* index) const
{
	if(!_depth){                           // it is a leave
		if(!first) os << " ";
		os << *_value;
	}
	else{
		os << "{";
		if(index) os << *index << ": ";      // if there is a jump
		typename map<K, TNode<K,E>* >::const_iterator pos =_children.begin();
		for(unsigned int i = 0, j=0; pos != _children.end(); ++pos, ++i, ++j){
			if(i != pos->first){      // if there is a jump
				pos->second->print_bracket(os, (bool)!j, &(pos->first));
				i = pos->first;
			}
			else pos->second->print_bracket(os, (bool)!j, NULL);
		}
		os << "}";
	}
}

//---------------------------------------------------------------------------
/** set the depht of the parents recursively */
template <class K, class E>
void TNode<K,E>::set_depth_backward()
{
		if(_pParent && _pParent->_depth == NOT_SET){
		_pParent->_depth = _depth+1;   // set the depth
		_pParent->set_depth_backward();   // iterate to the node
	};
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
template <class K, class E>
void TTree<K,E>::check_depth(unsigned int size) const{
	if(TNode<K,E>::_depth == NOT_SET) TNode<K,E>::_depth=size;
	else if(TNode<K,E>::_depth != size) throw("Tree: wrong depth of the tree!");
}


// ----------------------------------------------------------------------------------------
// read
// ----------------------------------------------------------------------------------------
/** reading of a tree with identical depth
	* the tree must be clean, ie.e no comments,...
	* format {{{ }{ }}{4: {3: }{ }}}
	* reads until the next bracket shows up and then passes the stream either
	* to a new lower node or back to the parent node
	*/

template <class K, class E>        // no text after the string allowed
void TTree<K,E>::read(string& ss)
{
	istringstream is(ss);
	read(is);                    // read the tree
	is >> ws;                        // check if the string is empty, i.e. the initial string contained only the tree!
	if(!is.eof()) throw("There is text at after the end of the tree!");
}

/** function to initialize the readin of the entire tree */
template <class K, class E>
void TTree<K,E>::read(istream &is)
{
	// reset the tree
	TNode<K,E>::reset();

	// remove the first bracket
	char c;
	is >> ws;
	is.get(c);
	if(c != '{') throw("The tree has to start with a bracket '{'!");

	// read the tree
	TNode<K,E>::read(is);
}

// ----------------------------------------------------------------------------------------
// read_tree
// ----------------------------------------------------------------------------------------
template <class K, class E>        // no text after the string allowed
void TNode<K,E>::read(string& ss)
{
	istringstream is(ss);
	read(is);
	is >> ws;
	if(!is.eof()) throw("There is text after the end of the tree!");
}

template <class K, class E>     	// reads until the tree is finsihed and returnes the rest of the stream
void TNode<K,E>::read(istream &is)
{
	string rowStr;
	E elmnt;
	char c;

	// get the text until a new bracket shows up
	while(is.get(c) && c!='}' && c!='{'){
		rowStr += c;
	}

	// if there is any text
	if(!rowStr.empty()){
		// if there are alredy children no text is allowed
		if(!_children.empty()) throw("Text between brackets not allowed!");

		// check if a row specification is present
		string::size_type pos = rowStr.find_first_of(':');
		if(pos != string::npos){
			istringstream ss(rowStr.substr(0, pos));
			unsigned int key;
			char b;
			if ((ss >> key) && !ss.get(b)){ // it is effectively a key specification: change the key
				if(_pParent->get_static(key)) throw("Key specification already used!");
				change_key(key); // change the key of the present node
				rowStr = rowStr.substr(pos+1, rowStr.size());	// remove the key
			}
		}

		// text is only allowed at leaves and not at nodes!
		if(c=='{' && !rowStr.empty()) throw("Text between brackets not allowed!");

		// read the elements elment by element and add them to the node
		istringstream ROW(rowStr);
		ROW >> ws;
		while(ROW >> elmnt) {
			push_back(elmnt);                      // add the element
			ROW >> ws;                             // for security reasons remove any space
		}
		if(!ROW.eof()) throw("The elements of the bracket could not be read!");
		rowStr.clear();                          // empty the string
	}

	// how to proceed?
	if(is.eof()) throw("The number of brackets is wrong!"); // end of the stream: not possible
	if(c=='{') return push_back()->read(is);	 // add a new node at the end and continue reading


	assert(c == '}');
	is >> ws;
	if(_pParent) return _pParent->read(is);    // go one up and continue reading

	// that's the root
	if(!is.eof())  throw("The matrix could not be read correctly!");
	return;                       // uff the end of the tree is reached
}
#endif


