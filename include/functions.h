/** @file functions.h
 *
 *   Copyright (C) 2006 Frederic Guillaume    <guillaum@zoology.ubc.ca>
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

#ifndef functionsH
#define functionsH

#include <vector>
#include <cstdlib>
#include <assert.h>
#include <time.h>
#include <string>
#include <cmath>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <map>
#include "tarray.h"
#include "types.h"
#include "ttree.h"
#include <iomanip>
#include <typeinfo>
#include <stdexcept>


#ifdef __BCPLUSPLUS__                 // if Borland is used
#include <direct.h>
#else
#include <unistd.h>
#endif




using namespace std;

class Param;

extern void message(const char* message, ...);
extern void warning(const char* message, ...);
extern void error(const char* message, ...);

extern string getElapsedTime(double time);
extern string getElapsedTime(time_t end, time_t start);

extern unsigned int nb_warning;
extern unsigned int nb_error;
extern unsigned int verbose_message;
extern unsigned int verbose_warning;
extern unsigned int verbose_error;

//------------------------------------------------------------------------------
/** check if a file exitis /
bool file_exist(string file)
{
	ifstream IN(file.c_str());
	if (!IN) return false;
    else IN.close();
    return true;
}

*/
 //------------------------------------------------------------------------------
inline string get_cwd()
{
	char temp[1024];
	if (!getcwd(temp, 1024)) return std::string("");
	string path(temp);
	char lastChar = path[path.size() - 1]; // add the directoy seperator if not present
	if (lastChar != '/' && lastChar != '\\') path += SEP;
	return path;
}

//------------------------------------------------------------------------------
/** returns the rounded value with the specified number of digits */
inline int my_round(double value)
{
	return (int) (value < 0 ? value - 0.5 : value + 0.5);
}

inline double my_round(double value, int digits)
{
	if (!digits) return my_round(value);
	return floor(value * pow((double) 10, digits) + 0.5)
			* pow((double) 10, -digits);
}

//------------------------------------------------------------------------------
/** remove trailing/forgoing/edge space of a string */
inline void rem_trailing_space(string& s)
{
	s.erase(s.find_last_not_of(" \t\r\n\v\f") + 1);
}

/** remove forgoing space of a string */
inline void rem_forgoing_space(string& s)
{
	s.erase(0, s.find_first_not_of(" \t\r\n\v\f"));
}

/** remove edge space of a string */
inline void rem_edge_space(string& s)
{
	rem_forgoing_space(s);
	rem_trailing_space(s);
}

//------------------------------------------------------------------------------
/** similar to sample in R:
 * returns nb UNIQUE integer numbers between 0 and max-1
 * the returned array must be deleted
 * the returned array is sorted!!!
 */
//unsigned int* sample(const unsigned int& nb, const unsigned int& max);

/** similar to sample in R WITHOUT replacement:
 * input:  -array to sample from (array is changed)
 *         -size: total size of the array
 *         -nb: points to sample
 *         -order: does the order of the samples matters? (true=yes)
 * output: -nothing, but the array is modified
 the sampled items are the "nb" first elements of the array
 *
 template <typename T>
 void sample(T* array, const unsigned int& size, const unsigned int& nb, bool order=true){
 if(size<nb) nb=size;
 if(!nb) return;                                 // as it could cause loop problems

 if(order || nb<size/2){                         // draw the "good" ones (start left)
 unsigned int pos, i;
 for(i = 0; i< nb; ++i){
 pos = TReplicate::r.Uniform(size-i);         // get the pos
 if(pos+i != i) ARRAY::swap(array, i, pos+i);
 }
 }
 else{                                           // draw the bad ones (start right)
 unsigned int pos, i;
 for(i = size-1; i>=nb; --i){
 pos = TReplicate::r.Uniform(i+1);         		// get the pos
 if(pos != i) ARRAY::swap(array, i, pos);
 }
 }
 }
 */
//------------------------------------------------------------------------------
/** returns the new population size assuming logistic growth limited by the carrying capacity
 * r of 0 is stable!!!!!!!!!!!!!!!!!!!!!!!!!!
 */
inline unsigned int logisticGrowth(const double& r, const unsigned int& K,
		const unsigned int& N)
{
	if (!K) return 0;
	return my_round(N + r * N * (1.0 - (double) N / K));
}

//------------------------------------------------------------------------------
/** returns the new population size assuming logistic growth limited by the carrying capacity
 * The Beverton-Hold model is implemented:  newN[i]<-N*r*K/(r*N-N+K)
 * Beverton RJH and Holt SJ (1957). On the Dynamics of Exploited Fish Populations.
 * The model (respectivly r) was modified such as r=0 results in a stable population size
 */
inline double beverton_hold(const double& r, const unsigned int& K,
		const unsigned int& N)
{
	if (!K) return 0;
	return N * K * (1 + r) / (N * (1 + r) - N + K);
}

//------------------------------------------------------------------------------
/** that is an implementation of the general logistic curve
 * (Richards, F.J. 1959 A flexible growth function for empirical use. J. Exp. Bot. 10: 290--300.)
 * The curve is defined by 5 parameters.
 * Do to size problems a slope of more than (+/-)1e4 is regarded as instantenious change from min to max
 */
inline double generalLogisticCurve(const double& x,     // time
		const double& min,   // the lower asymptote
		const double& max,   // the upper asymptote
		const double& max_r, // the time of maximum growth
		const double& r,     // the growth rate
		const double& s)   // affects near which asymptote maximum growth occurs
{

	// check if the slope exceeds the limits
	if (r >= 1e4) {   // slope is too big            => use instantenious change
		if (x < max_r) return min;
		else if (x > max_r) return max;
		else return min + ((max - min) / pow(1 + s, 1.0 / s));  // if x == max_r
	}
	else if (r <= -1e4) { // slope is too big (negative) => use instantenious change
		if (x < max_r) return max;
		else if (x > max_r) return min;
		else return min + ((max - min) / pow(1 + s, 1.0 / s));  // if x == max_r
	}

	// the slope is small enough to use the generalized logistic function
	assert(s);
	return min
			+ ((max - min)
					/ pow(1 + (s * (double) exp((long double) (max_r - x) * r)),
							1.0 / s));
}

// ----------------------------------------------------------------------------------------
// read_Fstat_file
// ----------------------------------------------------------------------------------------
/** this function reads any FSTAT file, a vector is returned which contains each individual.
 * Each individual is determined by a double unsigned int array:
 *       v[ind][0][0] = patch             // current patch (compulsory, index starting at 0)
 *       v[ind][0][1] = NaN               // currently not used
 *       v[ind][1][0] = nb_loci           // number of locus found in the file (compulsory)
 *       v[ind][1][1] = max_allele        // maximal index of the alleles (compulsory)
 *       v[ind][2][0] = age               // off: 0, adlt: 1
 *       v[ind][2][1] = sex               // mal: 0, fem: 1
 *       v[ind][3][0] = id                // id of the individual
 *       v[ind][3][1] = natal patch       // id of the individual's natal patch (index starting at 0)
 *       v[ind][4][0] = id_mom            // id of mother
 *       v[ind][4][1] = mom's patch       // id of mother's natal patch (index starting at 0)
 *       v[ind][5][0] = id_dad            // id of father
 *       v[ind][5][1] = dad's patch       // id of father's natal patch (index starting at 0)
 *       v[ind][loc+6][all] (compulsory)  // allele starting at 0
 * Note, that the arrays of the vector has to be deleted
 */
vector<unsigned int**>* read_Fstat_file(string filename);
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename T>     // for any data: make a string and then get the length
inline unsigned int getNbDigits(T n)
{
	string s = toStr(n);
	return (unsigned int)s.size();
}

template<>
inline unsigned int getNbDigits<unsigned int>(unsigned int n)
{
	unsigned int count = 0;
	if (!n) return 1;
	while (n > 0) {
		n /= 10;
		count++;
	}
	return count;
}

template<>
inline unsigned int getNbDigits<int>(int n)
{
	if (n < 0) return getNbDigits<unsigned int>((unsigned int) -n);
	return getNbDigits((unsigned int) n);
}

template<>
inline unsigned int getNbDigits<string>(string n)
{
	return (unsigned int)n.size();
}

template<>
inline unsigned int getNbDigits<char*>(char* n)
{
	string s(n);
	return (unsigned int)s.size();
}

double process_mem_usage();
double max_mem_usage();

unsigned int binomial_coefficient(unsigned int n, unsigned int m);
double rarefaction(unsigned int n, unsigned int k, unsigned int i);
double factln(const unsigned int& n);
double gammln(double xx);

// number - string conversion //////////////////////////////////////////////////

class BadConversion: public std::runtime_error {
public:
	BadConversion(const string& s) :
			std::runtime_error(s)
	{
	}
};

template<typename T>
inline string toStr(const T& x)
{
	ostringstream o;
	if (!(o << std::setprecision(10)<< x))
		throw BadConversion(string("toStr(") + typeid(x).name() + ")");
	return o.str();
}

template<typename T, typename M>
inline string toStr(const T& x, const M& maxTerm, char fillingChar = '0')
{
	ostringstream o;
    if (!(o << std::setfill(fillingChar) << setw(getNbDigits(maxTerm)) << right
			<< x)) {
		throw BadConversion(string("toStr(") + typeid(x).name() + ")");
	}
	return o.str();
}

template<typename T>
inline T strTo(const string& s, bool failIfLeftoverChars = true)
{
	T x;
	istringstream i(s);
	char c;
	if (!(i >> x)) {
		throw "This is not a number!"; //BadConversion(s);
	}
	i >> ws;                             // remove with space
	if (failIfLeftoverChars && i.get(c)) {
		throw "This is not a number!"; //BadConversion(s);
	}
	return x;
}

template<typename T>
inline T strTo_withNAN(const string& s, bool failIfLeftoverChars = true)
{
	T x;
	istringstream i(s);
	char c;
	if (!(i >> x)) return my_NAN;
	i >> ws;                             // remove with space
	if (failIfLeftoverChars && i.get(c)) return my_NAN;
	return x;
}

template<typename T>
inline bool isNumber(const string& s, bool failIfLeftoverChars = true)
{
	T x;
	istringstream i(s);
	char c;
	if (!(i >> x)) return false;
	i >> ws;                             // remove with space
	if (failIfLeftoverChars && i.get(c)) return false;
	return true;
}

inline bool isNumber(const char& c)
{
	return (c >= '0' && c <= '9');
}

inline bool isChar(const char& c)
{
	if (c <  'A') return false;
	if (c <= 'Z') return true;
	if (c <  'a') return false;
	if (c <= 'z') return true;
	return false;
}

/*
 int getTotalRAM();
 int getAvailRAM();
 int getTotalMemory();
 int getAvailMemory();
 */

extern int check_expiration(string exeDate, int exp, string& expDate);
// extern time_t to_time_t(char const *time);
extern tm to_tm(char const *time);

/** deleting sex specific arrays where both may point to the same array */
template<typename T>
void delete_sex_specific_array(T* &array1, T* &array2)
{
	if (array1) {
		if (array1 != array2 && array2) {
			delete[] array2;
		}
		delete[] array1;
	}
	array1 = array2 = NULL;
}

template<typename T>
void delete_sex_specific_array(T** &array)
{
	if (array[1]) {
		if (array[1] != array[0] && array[0]) {
			delete[] array[0];
		}
		delete[] array[1];
	}
	array[1] = array[0] = NULL;
}

/** is c the end of a line?
 * Windows: "\r\n"
 * Mac: 	"\r"
 * Linux: 	"\n"
 */
extern bool isEOL(char c);
extern bool isEOL(char c, istream& is);

extern istream& safeGetline(istream& is, string& t);
extern istream& safeIgnore(istream& is);


// ----------------------------------------------------------------------------------------
// convertion between centiMorgan (cM) and r (recombiantion rate)
// ----------------------------------------------------------------------------------------
extern double cM2r(double cM);
extern double r2cM(double r);

// ----------------------------------------------------------------------------------------
// set_path if it does not exist
// ----------------------------------------------------------------------------------------
void check_path(string& path);

// -----------------------------------------------------------------------------
inline double ProbDensNorm(double x, double mu = 0, double sd = 1)
{
    double m = (x - mu) / sd;
    return exp(-m * m / 2.) / (sqrt(2 * PI) * sd);
}

// -----------------------------------------------------------------------------
vector<unsigned int> get_logtime_occurences(Param* pParam, unsigned int totGen, bool isCoal);
extern void runScript(string scriptname, string filename);


#endif
