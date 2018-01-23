/** @file tmatrix.cpp
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

#include "tmatrix.h"
#include <fstream>
using namespace std;



template<class T> T TMatrixVar<T>::matrix_NAN = my_NAN;
template<> string TMatrixVar<string>::matrix_NAN = my_NANstr;
template<class T> T TMatrixVar<T>::matrix_STR = my_STR;
template<> string TMatrixVar<string>::matrix_STR = my_STRstr;

// ----------------------------------------------------------------------------------------
// read
// ----------------------------------------------------------------------------------------
void TMatrix::read(string text, bool allowStr)
{
	reset();
	unsigned int i, j, size;

	// read it a as a TMatrixVar object
	TMatrixVar<double> m;
	try {
		m.read(text, allowStr);
	}
	catch(const char* err) {
		error("ParameterMatrix '%s': %s\n", text.c_str(), err);
	}

	_oneDims = m.get_oneDims();
	_rows = m.getNbRows();

	// find the number of columns (it could have empty columns)
	for (i = 0; i < _rows; i++) {
		_cols = m.getNbCols(i);
		if (_cols) break;
	}

	//check for matrix coherence:
	for (i = 0; i < _rows; i++) {
		size = m.getNbCols(i);
		if (!size) {            // if it is a empty column, fill it with my_NANs
			for (j = 0; j < _cols; j++) {
				m(i)->push_back(my_NAN);
			}
		}
		else if (size != _cols) {
			throw "Could not read matrix: not same number of columns per row!";
		}
	}

	//copy to input TMatrix:
	if (_val) delete _val;
	_val = new double[_cols * _rows];
	for (i = 0; i < _rows; ++i) {
		for (j = 0; j < _cols; ++j) {
			_val[i * _cols + j] = m.get(i, j);
		}
	}
	if (m.get_strMatrix()) _strMatrix = *m.get_strMatrix();
}

// ----------------------------------------------------------------------------------------
// read_matrix
// ----------------------------------------------------------------------------------------
/** read a matrix file and return the file information if available
 * with_bracket defines if the matrix is defined with brackets or just by new lines:
 * 0: the program decides (default)
 * 1: no bracket
 * 2: bracket
 */
map<string, int> TMatrix::read_matrix(string filename, unsigned int with_bracket)
{
	ifstream IN(filename.c_str());
	if (!IN) error("File '%s' could not be opened!", filename.c_str());

	unsigned int line = 0;
	map<string, int> fileInfo;

	// read the file info if available
	STRING::removeCommentAndSpace(IN, line);
	fileInfo = STRING::readFileInfo(IN, line);
	if (fileInfo.empty())
		error("The file info box of file '%s' could not be read!\n",
				filename.c_str());

	STRING::removeCommentAndSpace(IN, line);

	// read the matrix (check if the matrix has brackets)
	try {
		switch(with_bracket){
			default:
			case 0:	if (IN.peek() == '{') read_matrix_bracket(IN, line);
			else read_matrix_none(IN, line);
			break;
			case 1: read_matrix_none(IN, line);
			break;
			case 2: read_matrix_bracket(IN, line);
			break;
		}
	}
	catch(const char* err) {
		error("File '%s': %s\n", filename.c_str(), err);
	}

	IN.close();

	return fileInfo;
}

// check if the matrix is
// ----------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------
// read_matrix_bracket
// ----------------------------------------------------------------------------------------
/** reading a matrix from a stream. The matrix has to be specified by brackets */
void TMatrix::read_matrix_bracket(istream & IN) {
	unsigned int line = 0;
	read_matrix_bracket(IN, line);
}

void TMatrix::read_matrix_bracket(istream & IN, unsigned int& line) {
	vector<vector<double> > tmpMat;
	double elmnt;
	char c;
	int openBrack = 1,  // the one which was just removed
			closeBrack = 0;

	//remove the first enclosing bracket
	STRING::removeCommentAndSpace(IN, line);
	IN.get(c);
	assert(c=='{');

	//then read the rows
	for (int i = 0; openBrack != closeBrack && IN; ++i) {
		tmpMat.push_back(vector<double>());  // add a new row

		// go to next row
		STRING::removeCommentAndSpace(IN, line);
		IN.get(c);             //first character:
		if (c == '{') ++openBrack;
		else throw "Could not read matrix: brackets are misplaced!";

		//read a row enclosed by {...}:
		while (IN) {
			STRING::removeCommentAndSpace(IN, line);
			c = IN.peek();
			if (c == ',' || c == ';') continue;
			if (c == '}') {
				++closeBrack;
				IN.get(c);
				break;
			}    //go to next row

			//read a row element:
			if (IN >> elmnt) tmpMat[i].push_back(elmnt);
			else throw("Could not read matrix!");

		}
	}

	//check for matrix coherence:
	unsigned int rows = (unsigned int)tmpMat.size();      // get the number of rows
	unsigned int cols = (unsigned int)tmpMat[0].size();   // get the size of the first row
	unsigned int i, j;
	for (i = 1; i < rows; i++) {
		if (tmpMat[i].size() != cols) {
			throw "Could not read matrix: not same number of elements in all rows!";
		}
	}

	//copy to input TMatrix:
	reset(rows, cols);
	for (i = 0; i < rows; ++i) {
		for (j = 0; j < cols; ++j) {
			set(i, j, tmpMat[i][j]);
		}
	}
}

// ----------------------------------------------------------------------------------------
// read_matrix_none
// ----------------------------------------------------------------------------------------
/** reading a matrix from a stream. The matrix has to be specified WITHOUT brackets,
 * each line corresponds to a row
 */
void TMatrix::read_matrix_none(istream & IN) {
	unsigned int line = 0;
	read_matrix_none(IN, line);

}

void TMatrix::read_matrix_none(istream & IN, unsigned int& line) {

	vector<vector<double>*> tmpMat;
	unsigned int lineTemp, i=0;
	double elmnt;
	string lineStr, text;

	// read line by line
	while (STRING::readUntilCharacter(IN, lineStr, line, '\n', true, false)) {
		if (lineStr.empty()) continue;	// if empty line
		istringstream LINE;             // make a new stream
		LINE.str(lineStr);              // allocate the line to the new stream
		LINE >> ws;

		//read a row (no more comments present):
		tmpMat.push_back(new vector<double>);  // add a new row
		while (LINE.good() && !LINE.eof()) {
			//read a row element:
			if (LINE.peek() == '{') {
				_strMatrix[i][(unsigned int)tmpMat[i]->size()] = STRING::readBrackets2String(LINE, lineTemp, '}');
				tmpMat[i]->push_back(my_NAN);
			}
			else {
				if (LINE >> elmnt) tmpMat[i]->push_back(elmnt);
				else if (LINE >> text && text == my_NANstr) {
					tmpMat[i]->push_back(my_NAN);
				}
				else throw("Could not read matrix!");
			}
			LINE >> ws;
		}
		++i;
	}

	//check for matrix coherence:
	unsigned int rows = (unsigned int)tmpMat.size();       // get the number of rows
	unsigned int cols = (unsigned int)tmpMat[0]->size();   // get the size of the first row
	unsigned int j;
	for (i = 1; i < rows; i++) {
		if (tmpMat[i]->size() != cols) {
			throw "Could not read matrix: not same number of elements in all rows of the matrix!";
		}
	}

	//copy to input TMatrix:
	reset(rows, cols);
	for (i = 0; i < rows; ++i) {
		for (j = 0; j < cols; ++j) {
			set(i, j, (*tmpMat[i])[j]);
		}
		delete tmpMat[i];
	}
}

