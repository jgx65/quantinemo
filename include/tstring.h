/** @file tstring.h
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

#ifndef tstringH
#define tstringH

#include "types.h"
#include <string>
#include <sstream>
#include <map>
#include <vector>
#include <cmath>

using namespace std;
//---------------------------------------------------------------------------
class STRING {
private:


public:

	/** file reading routines */
	static bool removeCommentAndSpace(istream & IN, unsigned int& line, int comment=0, char end='\0', bool removed = false);
	static string readBrackets2String(istream & IN, unsigned int& linecnt, char dim=')');

	static string find_first_of(istream & IN, string chars, bool stopAtSpace=false, unsigned int include=0);
	static bool readUntilCharacter(istream & IN, string& out, unsigned int& linecnt, char item, bool includeFirst, bool includeLast);
	static string readUntilCharacter(istream & IN, unsigned int& linecnt, char item, bool includeFirst, bool includeLast);
	static map<string, int> readFileInfo(istream & IN, unsigned int& line);

	static bool file_exists(const string& name);
	static bool file_exists(string& name, const string& dir);
    
    static string remove_file_extension(string file);
    static string remove_file_path(string file);
    static string get_file_path(string file);
    static string get_file_basename(string file);
    static string get_file_name(string file);


	// -----------------------------------------------------------------------------
	template <typename T>
	static vector<T> strMatrix2vector(string m){
		vector<T> vec;
		T val;

		istringstream LINE;                   // make a new stream
		// assert(m[0] == '{' && m[m.length()-1] == '}');
		string s = m.substr(1, m.length()-2); // removing the brackets
		LINE.str(s);  // allocate the matrix to the new stream

		while(LINE >> val){
			vec.push_back(val);
		}
		return vec;
	}

	//------------------------------------------------------------------------------
	/** transformes the string matrix (1 dimensional) to an array of doubles, the
     array has to be created before and be passed  */
	template <typename T>
	static void strMatrix2array(const string& m, T* array, const int& size){
		istringstream LINE;                   // make a new stream
		assert(m[0] == '{' && m[m.length()-1] == '}');
		LINE.str(m.substr(1, m.length()-2));  // allocate the matrix to the new stream (removing the brackets)
        
		for(int i=0; i<size; ++i){
			LINE >> array[i];
		}
	}
    

	//------------------------------------------------------------------------------
};
#endif
