/** @file tmatrix.h
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
//---------------------------------------------------------------------------

#ifndef tmatrixH
#define tmatrixH
//---------------------------------------------------------------------------

#include "functions.h"
#include "tstring.h"

using namespace std;



/** This matrix is variable in both dimensions and allows to read any matrix input
 * The variable matrix may be of any kind. Values are stored in a double vector.
 * If the number of column is not identical for all rows, the function get_nbCols returns NaN.
 **/

template <class T>
class TMatrixVar {
private:

	unsigned int _rows, _cols;   // if the rows have not the same number of elements: _cols = my_NAN
	bool         _oneDims;       // was the input matrix just 1D?

	vector< vector<T> > _vals;


	map<unsigned int, map<unsigned int, string> > _strMatrix;     // _strMatrix[row][col]: for text in number matrixes

	static T matrix_NAN;
	static T matrix_STR;

public:

	TMatrixVar() : _rows(0), _cols(0), _oneDims(0){ }
	TMatrixVar(string text): _rows(0), _cols(0), _oneDims(0) {read(text);}

	TMatrixVar  (const TMatrixVar<T>& mat) {
		_rows   = mat._rows;
		_cols   = mat._cols;
		_vals   = mat._vals;
		_oneDims = mat._oneDims;
		_strMatrix = mat._strMatrix;
	}

    ~TMatrixVar () {	}

	/**@brief Accessor to element at row i and column j**/
	T operator()(const int& i, const int& j) {return get(i,j);}
	T get (unsigned int i, unsigned int j) {
		assert(i<_rows && j<_vals[i].size());
		return _vals[i][j];
	}

	/**@brief Accessor to row i**/
	vector<T>* operator()(const int& i) {return get(i);}
	vector<T>* get (unsigned int i) {
		assert(i < _rows);
		return &_vals[i];
	}
	/**@brief Accessor to the whole double vector**/
	vector<vector<T> >* get() {return &_vals;}

	/**@brief Accessor to the matrix dimensions
	 *@param dims an array of at least 2 elements to store the row [0] and column [1] numbers. May be NULL.
	 *@return the total size of the matrix
	 **/
	unsigned int get_dims (unsigned int* dims) {
		if(dims) {
			dims[0] = _rows;
			dims[1] = _cols;
		}
		return _rows*_cols;
	}

	/** add a new element into row i of the matrix **/
	void add(const unsigned int& i, T elem){
		while(i >= _vals.size()){
			_vals.push_back(vector<T>());                    // add a new row
		}
		_rows = _vals.size();
		_vals[i].push_back(elem);
	}

	void set_row(const unsigned int& i, vector<T> row){
		while(i >= _vals.size()){
			_vals.push_back(vector<T>());                    // add a new row
		}
		_rows = _vals.size();
		_vals[i] = row;
	}

	/** set the dimension of the matrix. If not all the rows have the same number
	 * of elements, _cols will be set to my_NAN.
	 */
	void set_dimensions(){
		// get the number of rows
		_rows = (unsigned int)_vals.size();
		if(!_rows) {_cols = my_NAN;  return;}    // check if the matrix is empty

		// check if all rows have the same number of elements
		_cols = (unsigned int)_vals[0].size();
		for(unsigned int i=1; i<_rows; ++i){
			if(_vals[i].size() != _cols){
				_cols = my_NAN;
				return;
			}
		}
	}

	bool         get_oneDims ( ){return _oneDims;}         /**is it a 1D matrix.*/
	unsigned int getNbRows ( ) {return _rows;}             /**Gives the number of rows.*/
	unsigned int getNbCols ( ) {return _cols;}             /**Gives the number of columns. (my_NAN if the columns have different sizes)*/

	/**Gives the number of columns of the row i */
	unsigned int getNbCols (unsigned int i) {
		assert(i<_rows);
		return (unsigned int)_vals[i].size();
	}
    
    
    /** get the totalnumber of elements */
    unsigned int getNbElements(){
        unsigned int nb=0;
        typename vector<vector<T> >::iterator cur, end;
        for(cur=_vals.begin(), end=_vals.end(); cur!=end; ++cur){
            nb += cur->size();
        }
        return nb;
    }

    map<unsigned int, map<unsigned int, string> >* get_strMatrix(){return _strMatrix.empty() ? NULL : &_strMatrix;}

	/**Reset members to zero state.*/
	void reset ( ) {
		_rows = 0;
		_cols = 0;
		_vals.clear();
	}

	/** extracts from a string the matrix (1D or 2D) */
	void read(string text, bool allowStr=false){
		reset();
		istringstream FILE;
		FILE.str(text);
		string rowStr;

		T elmnt;
		char c;
		_oneDims = false;
		unsigned int line=0;
        
		//remove the first enclosing bracket
		FILE >> c;
		assert(c=='{');

		//then read the rows
		unsigned int row=0, col;
		while(FILE.peek() != '}' && !_oneDims){
			FILE >> ws;
			if(FILE.peek()!='{'){
				if(!_vals.empty()) throw "Could not read matrix: text between rows!";

				// it is a one dimensional array
				FILE.putback('{');
				_oneDims = true;
			}

			_vals.push_back(vector<T>());                    // add a new row

			//read a row enclosed by {...}:
			rowStr = STRING::readUntilCharacter(FILE, line, '}', false, false);

			// check if it has a row indicator and add the missing rows if necessary
			string::size_type pos = rowStr.find_first_of(':');
			if(pos != string::npos){
				unsigned int nb = strTo<unsigned int>(rowStr.substr(0, pos));
				rowStr = rowStr.substr(pos+1, rowStr.size());
				if(nb <_vals.size()) throw "Could not read matrix: row indicator must exceed the number of the previous row!";
				while(_vals.size() < nb){
					_vals.push_back(vector<T>());               // add a new row
				}
			}

			// read the row
			istringstream ROW;
			ROW.str(rowStr);
			ROW >> ws;
			col=0;		// column counter
			while(!ROW.eof()) {
				if(ROW >> elmnt) _vals.back().push_back(elmnt);   // add the element
				else {
					ROW.clear();
					if(ROW >> rowStr && rowStr==my_NANstr) _vals.back().push_back(matrix_NAN);      // if a NaN is read
					else if(allowStr){
						_vals.back().push_back(matrix_STR);
						_strMatrix[row][col]=rowStr;

					}
					else throw("Could not read the matrix!");
				}
				ROW >> ws;
				++col; // column counter
			}
			++row;	// row counter
		}
		set_dimensions();
	}

	friend ostream& operator<<(ostream &os, const TMatrixVar<T>& m){
		unsigned int c, r;
		for(c=0; c<m._vals.size(); ++c){
			os << '{';
			for(r=0; r<m._vals[c].size(); ++r){
				if(r) os << ' ';
				os << m._vals[c][r];
			}
			os << '}';
		}
		return os;
	}
};


////////////////////////////////////////////////////////////////////////////////
/**A class to handle matrix in params, coerces matrix into a vector of same total size**/
class TMatrix {
private:

	unsigned int _rows, _cols;

	bool _oneDims;

	double* _val;

	map<unsigned int, map<unsigned int, string> > _strMatrix;     // _strMatrix[line][col]: for text in number matrixes

public:

	TMatrix() : _rows(0), _cols(0), _oneDims(0), _val(0) { }

	TMatrix(string tt, bool allowStr=false) : _rows(0), _cols(0), _val(0){read(tt,allowStr);}

	/**copy constructor.**/
	TMatrix(const TMatrix& mat)
	{
		_rows = mat._rows;
		_cols = mat._cols;
		_oneDims = mat._oneDims;
		_val = new double [_rows*_cols];
		ARRAY::copy<double>(_val, mat._val, _rows*_cols);

		_strMatrix = mat._strMatrix;
	}

	/**@brief Creates an array of doubles of size = rows*cols**/
	TMatrix ( unsigned int rows, unsigned int cols )
	{
		_val = new double [rows*cols];
		_rows = rows;
		_cols = cols;
		_oneDims = false;
	}

	~TMatrix () {
		if(_val) delete [] _val;
	}

	bool operator==(const TMatrix& m){
		if(_rows != m._rows) return false;
		if(_cols != m._cols) return false;
		for(unsigned int i=0; i<_rows*_cols; ++i){
			if(_val[i] != m._val[i]) return false;
		}
		return true;
	}

	bool operator!=(const TMatrix& m) {return !(*this==m);}

	map<unsigned int, map<unsigned int, string> >* get_strMatrix(){return _strMatrix.empty() ? NULL : &_strMatrix;}

	/**Sets element at row i and column j to value val**/
	inline void set (unsigned int i, unsigned int j, double val) {
		assert(i<_rows);
		assert(j<_cols);
		_val[i*_cols + j] = val;
	}

	inline void assign (int val){
        memset(_val, val, _rows * _cols * sizeof(double));
	}

	/**Re-allocate the existing matrix with assigned rows and cols dimensions**/
	inline void reset (unsigned int rows, unsigned int cols) {
		if(_val) delete [] _val;
		_val = new double [rows*cols];
		_rows = rows;
		_cols = cols;
	}

	/**Reset the existing matrix to the new dimensions and copies the array**/
	inline void reset (unsigned int rows, unsigned int cols, double* array) {
		reset(rows, cols);
		ARRAY::copy(_val, array, rows*cols);
	}
	/**Reset members to zero state.*/
	inline void reset ( ){
		_rows = 0;
		_cols = 0;
		if(_val) delete [] _val;
		_val = NULL;
	}

	/**@brief Accessor to element at row i and column j**/
	inline double get (unsigned int i, unsigned int j) const {
		assert(i<_rows);
		assert(j<_cols);
		return _val[i*_cols + j];
	}
	/**Accessor to the whole array.*/
	inline double* get () const {return _val;}

	/**Accessor to the matrix dimensions.
	 *@param dims an array of at least 2 elements to store the row [0] and column [1] numbers. May be NULL.
	 *@return the total size of the matrix
	 **/
	inline unsigned int get_dims (unsigned int* dims) {
		if(dims) { dims[0] = _rows; dims[1] = _cols; }
		return _rows*_cols;
	}
	inline unsigned int get_dims (unsigned int& xMax, unsigned int& yMax) {
		xMax = _rows;
		yMax = _cols;
		return _rows*_cols;
	}


	inline bool         get_oneDims ( ) const {return _oneDims;}         /**is it a 1D matrix.*/
	inline unsigned int getNbRows ( ) const {return _rows;}             /**Gives the number of rows.*/
	inline unsigned int getNbCols ( ) const {return _cols;}             /**Gives the number of columns.*/
	inline unsigned int length    ( ) const {return _rows*_cols;}       /**Returns the number of elements in the matrix.*/

	/**Gives access to a column of the matrix.*/
	void getColumnView (unsigned int col, unsigned int n, double* array){
		assert(col < _cols);
		assert(n == _rows);
		for(unsigned int i = 0; i < _rows; ++i){
			array[i] = _val[i*_cols + col];
		}
	}

	/**Gives access to a row of the matrix.*/
	void getRowView (unsigned int row, unsigned int n, double* array){
		assert(row < _rows);
		assert(n == _cols);
		for(unsigned int i = 0, stride = row*_cols; i < _cols; ++i){
			array[i] = _val[stride + i];
		}
	}

	/**Adds a value to an element of the matrix.*/
	inline void plus (unsigned int i, unsigned int j, double value){
		assert(i<_rows);
		assert(j<_cols);
		_val[i*_cols + j] += value;
	}

	/**Substracts a value from an element of the matrix.*/
	inline void minus (unsigned int i, unsigned int j, double value){
		assert(i<_rows);
		assert(j<_cols);
		_val[i*_cols + j] -= value;
	}

	/**Multiply an element of the matrix by a value.*/
	inline void multi (unsigned int i, unsigned int j, double value){
		assert(i<_rows);
		assert(j<_cols);
		_val[i*_cols + j] *= value;
	}

	/**Divide an element of the matrix by a value.*/
	inline void divide (unsigned int i, unsigned int j, double value){
		assert(i<_rows);
		assert(j<_cols);
		if(value) _val[i*_cols + j] /= value;
		else      _val[i*_cols + j] = my_NAN;
	}

	void show_up(){
		message("TMatrix dimensions: rows = %i, columns = %i\n",_rows,_cols);
		for(unsigned int i = 0; i < _rows; i++) {
			for(unsigned int j = 0; j < _cols; j++){
				message("%.3f ",_val[i*_cols + j]);
			}
			message("\n");
		}
	}

	inline string toStr(){
		ostringstream o;
		o << *this;
		return o.str();
	}

	// functions to read matrixes from files or streams
	map<string, int> read_matrix(string filename, unsigned int with_bracket=0);
	map<string, int> read_matrix(string filename, unsigned int& line);
	void             read_matrix_bracket (istream & FILE);
	void             read_matrix_bracket (istream & FILE, unsigned int& line);
	void             read_matrix_none (istream & FILE);
	void             read_matrix_none (istream & FILE, unsigned int& line);

	void read(string text, bool allowStr=false);

	friend ostream& operator<<(ostream &os, const TMatrix& m)
	{
		unsigned int c, r;
		if(!m._oneDims) os << '{';
		for(r=0; r<m._rows; ++r){
			os << '{';
			for(c=0; c<m._cols; ++c){
				if(c) os << ' ';
				os << m.get(r, c);
			}
			os << '}';
		}
		if(!m._oneDims) os << '}';
		return os;
	}
};

#endif

