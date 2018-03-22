/** @file functions.cpp
 *
 *   Copyright (C) 2006 Frederic Guillaume    <guillaum@zoology.ubc.ca>
 *   Copyright (C) 2008 Samuel Neuenschwander <samuel.neuenschwander@unil.ch>
 *   Copyright (C) 2018 Frederic Michaud <frederic.michaud@unil.ch>
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


#include "functions.h"

#include <sys/stat.h>


#ifdef __BCPLUSPLUS__                 // if Borland is used
#include <windows.h>
#include <psapi.h>
#include <stdlib.h>
#include <stdio.h>
#include <dir.h>
#else
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#endif

#include <stdarg.h>
#include <fstream>

#include <dirent.h>

#include "param.h"

using namespace std;


void message (const char* message, ...)
{
    if(!verbose_message) return;
    va_list ap;
	va_start(ap, message);
	vprintf(message,ap);
	va_end(ap);
}

void warning (const char* message, ...)
{
    if(!verbose_warning) return;
    ++nb_warning;
    va_list ap;
    va_start(ap, message);
    printf("\n***WARNING*** ");
    vprintf(message,ap);
    va_end(ap);
}

void error (const char* message, ...)
{
    if(!verbose_error) return;
	++nb_error;
	va_list ap;
	va_start(ap, message);
	fprintf(stderr,"\n***ERROR*** ");
	vfprintf(stderr,message,ap);
	va_end(ap);
	throw 1;
}

//------------------------------------------------------------------------------
string getElapsedTime(time_t end, time_t start)
{
    return getElapsedTime(difftime(end, start));
}

//------------------------------------------------------------------------------
string getElapsedTime(double time)
{
    unsigned int e_time = (unsigned int) time;
    unsigned int hour = e_time / 3600;
	unsigned int min  = ((e_time % 3600) / 60);
	unsigned int sec  = (e_time % 3600) % 60;

	ostringstream o;
	if(hour<10) o << setfill('0') << setw(2) << hour;
	else o << hour;
	o << ":" << setfill('0') << setw(2) << min;
	o << ":" << setfill('0') << setw(2) << sec;

	return o.str();
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
 * Note, that the arrays of the vector has to be deleted
 */
vector<unsigned int**>*
read_Fstat_file(string filename){
	vector<unsigned int**>* v = new vector<unsigned int**>();     // vector containing the individuals
	try{
		ifstream FILE(filename.c_str());
		if(!FILE) error("FSTAT file '%s' could not be opened!\n", filename.c_str());

		// read the heading line
		unsigned int nbPop, nbLoci, maxAllele, nbDigit, i;
		FILE >> nbPop >> nbLoci >> maxAllele >> nbDigit >> ws;
		nbDigit = (unsigned int)pow((double)10, (int)nbDigit);

		for(i = 0; i < nbLoci; ++i) {
			safeIgnore(FILE);  // remove line by line
			if(FILE.fail()) {
					error("FSTAT file: The number of loci does not match the heading number (%i)!\n", nbLoci);
			}
		}
		// read individual by individual
		unsigned int** curInd;
		unsigned int val, val1;
		string text;
		string text1;
		string::size_type pos;
		while(FILE.good() && !FILE.eof()){      // for each individual
			// get the line of an individual
			getline(FILE, text);           				// read the line
			istringstream LINE;             			// make a new stream
			LINE.str(text);              					// allocate the line to the new stream
			LINE >> ws;
			if(LINE.eof()) continue;              // the line is empty

			// create the array for the individual
			curInd = ARRAY::new_2D(nbLoci+6, 2, (unsigned int) my_NAN);

			LINE >> val;                               // get the pop (index starts at 1)
			curInd[0][0] = val-1;                      // store the pop (index starts at 0)
			curInd[1][0] = nbLoci;                     // store the number of loci
			curInd[1][1] = maxAllele;                  // store the maximal index of the alleles

			// get the alleles
			for(i=0; i<nbLoci; ++i){                   // read locus by locus
				try{
					// read the locus genotype
					if(!(LINE >> val)) error("FSTAT file: The number of loci (%i) does not match the heading number (%i, line %i)!\n",
							i+1, nbLoci, v->size()+1);
					if(!val)error("FSTAT file: quantiNemo2 does not support missing genotypes (line %i)!\n", v->size()+1);

					// first allele of the locus
					val1 = val/nbDigit;  // get first allele
					if(!val1 || val1>maxAllele) error("FSTAT file: The allele index (%i) must be within the range 1 and %i (line %i, locus %i)!\n",
							val1, maxAllele, v->size()+1, i+1);
					curInd[i+6][0] = val1-1;

					// second allele of the locus
					val1 = val % nbDigit; // get second allele
					if(!val1 || val1>maxAllele) error("FSTAT file: The allele index (%i) must be within the range 1 and %i (line %i, locus %i)!\n",
							val1, maxAllele, v->size()+1, i+1);
					curInd[i+6][1] = val1-1;
				}
				catch(...) {
					error("FSTAT file: The allele index is not a number (line %i, locus %i)!\n", v->size()+1, i+1);
				}
			}

			// check if the line is finsihed or has suplement information
			LINE >> ws;
			if(!LINE.eof()){       // there is suplement information => read it
				LINE >> val;  curInd[2][0] = (val-1);   // age
				LINE >>       curInd[2][1];             // sex

				LINE >> text;                           // id of individual
				pos = text.find('_');
				if(pos != string::npos){
					curInd[3][0] = strTo<unsigned int>(text.substr(0,pos));   // index
					curInd[3][1] = strTo<unsigned int>(text.substr(pos+1));   // patch
					if(curInd[3][1] != my_NAN) --curInd[3][1];
				}
				else  error("FSTAT file: The suplement individual information or the number of loci is incorrect(line %i)!\n",v->size()+1);


				LINE >> text;                           // id of mother
				pos = text.find('_');
				if(pos != string::npos){
					curInd[4][0] = strTo<unsigned int>(text.substr(0,pos));   // index
					curInd[4][1] = strTo<unsigned int>(text.substr(pos+1));   // patch
					if(curInd[4][1] != my_NAN) --curInd[4][1];
				}
				else  error("FSTAT file: The suplement individual information or the number of loci is incorrect(line %i)!\n",v->size()+1);

				LINE >> text;                           // id of father
				pos = text.find('_');
				if(pos != string::npos){
					curInd[5][0] = strTo<unsigned int>(text.substr(0,pos));   // index
					curInd[5][1] = strTo<unsigned int>(text.substr(pos+1));   // patch
					if(curInd[5][1] != my_NAN) --curInd[5][1];
				}
				else  error("FSTAT file: The suplement individual information or the number of loci is incorrect(line %i)!\n",v->size()+1);

				LINE >> text >> ws;											// remove the fitness: fitness is not considered
				if(!LINE.eof())  error("FSTAT file: The suplement individual information is incorrect(line %i)!\n",v->size()+1);
			}
			v->push_back(curInd);      // add the individual to the vector
		}

		FILE.close();
	}catch(...) {throw("Problems reading the FSTAT file!");}

	return v;
}

// -----------------------------------------------------------------------------
/** From the Numerical Recieps.
 * Returns the binomial coefficient
 * The floor function cleans up roundoff error for smaller values of n and k.
 */
unsigned int binomial_coefficient(unsigned int n, unsigned int k)
{
	return (unsigned int) floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
}

// -----------------------------------------------------------------------------
/** Computes the rarefaction  r = C(n-i, k) / C(n, k)
 *  n: current sample size
 *  k: smallest sample size
 *  i: number of allele i present in this sample
 */
double rarefaction(unsigned int n, unsigned int k, unsigned int i)
{
	return exp(factln(n-i)+factln(n-k)-factln(n)-factln(n-i-k));
}

// -----------------------------------------------------------------------------
/**From the Numerical Recieps.
 * Returns ln(n!).
 */
double factln(const unsigned int& n)
{
	if (n <= 1) return 0.0;
	static float a[101]; // A static array is automatically initialized to zero.
	if (n <= 100) return a[n] ? a[n] : (a[n]=gammln(n+1.0)); // In range of table.
	else return gammln(n+1.0); // Out of range of table.
}

// -----------------------------------------------------------------------------
/**From the Numerical Recieps.
 A function to return the natural log of the Gamma function of x.
 Adapted from the algorithm described in the book Numerical Recipes by
 Press et al.
 It can be used to compute factorials since ln(n!) = lnGamma(n + 1)
 */
double gammln(double xx)
{
    double x, y, tmp, ser = 1.000000000190015;
    static double cof[6] = { 76.18009172947146, -86.50532032941677,
        24.01409824083091, -1.231739572450155, 0.1208650973866179e-2,
        -0.5395239384953e-5 };
    int j;
    y = x = xx;
    tmp = x + 5.5;
    tmp -= (x + 0.5) * log(tmp);
    for (j = 0; j < 6; ++j)
        ser += cof[j] / ++y;
    
    return -tmp + log(2.5066282746310005 * ser / x);
}

//////////////////////////////////////////////////////////////////////////////
/* the function returns the memory consumption in MBs
 * On failure, returns 0.0
 */
double process_mem_usage()
{
	double vsize=0;
#ifdef  __BCPLUSPLUS__                 // if Borland is used
	/*
		DWORD aProcesses[1024], cbNeeded;
		unsigned int i;

		if(!EnumProcesses(aProcesses, sizeof(aProcesses), &cbNeeded ) ) return 0.0;
		PROCESS_MEMORY_COUNTERS pmc;
		HANDLE hProcess = OpenProcess(  PROCESS_QUERY_INFORMATION |
																		PROCESS_VM_READ,
																		FALSE, *aProcesses);
		if(NULL == hProcess) return 0.0;
		if(!GetProcessMemoryInfo(hProcess, &pmc, sizeof(pmc))) return 0.0;
		CloseHandle(hProcess);
		return pmc.PagefileUsage/1024;
	 */
#else                    // LINUX
	// 'file' stat seems to give the most reliable results
	ifstream stat_stream("/proc/self/stat",ios_base::in);

	// dummy vars for leading entries in stat that we don't care about
	string pid, comm, state, ppid, pgrp, session, tty_nr,
	tpgid, flags, minflt, cminflt, majflt, cmajflt,
	utime, stime, cutime, cstime, priority, nice,
	O, itrealvalue, starttime;

	stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
	>> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
	>> utime >> stime >> cutime >> cstime >> priority >> nice
	>> O >> itrealvalue >> starttime >> vsize; // don't care about the rest

	vsize /= 1024.0 * 1024.0;
#endif
	return vsize; // code should never pass here
}

//////////////////////////////////////////////////////////////////////////////
/* the function returns the memory consumption in MBs
 * On failure, returns 0.0
 */
double max_mem_usage()
{
	double vsize=0;
#ifdef  __BCPLUSPLUS__                 // if Borland is used
	/*
		DWORD aProcesses[1024], cbNeeded;
		unsigned int i;

		if(!EnumProcesses(aProcesses, sizeof(aProcesses), &cbNeeded ) ) return 0.0;
		PROCESS_MEMORY_COUNTERS pmc;
		HANDLE hProcess = OpenProcess(  PROCESS_QUERY_INFORMATION |
																		PROCESS_VM_READ,
																		FALSE, *aProcesses);
		if(NULL == hProcess) return 0.0;
		if(!GetProcessMemoryInfo(hProcess, &pmc, sizeof(pmc))) return 0.0;
		CloseHandle(hProcess);
		return pmc.PagefileUsage/1024;
	 */
#else                    // LINUX
	// 'file' stat seems to give the most reliable results
	ifstream stat_stream("/proc/self/stat",ios_base::in);

	// dummy vars for leading entries in stat that we don't care about
	string pid, comm, state, ppid, pgrp, session, tty_nr,
	tpgid, flags, minflt, cminflt, majflt, cmajflt,
	utime, stime, cutime, cstime, priority, nice,
	O, itrealvalue, starttime;

	stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
	>> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
	>> utime >> stime >> cutime >> cstime >> priority >> nice
	>> O >> itrealvalue >> starttime >> vsize; // don't care about the rest

	vsize /= 1024.0 * 1024.0;
#endif
	return vsize; // code should never pass here
}


// ----------------------------------------------------------------------------------------
// check_expiration
// ----------------------------------------------------------------------------------------
/** stops the program if the expiration time has reached
 * exeDate has the format May 19 2011
 * exp is in days
 * returns the expiration date if still valid
 */
int check_expiration(string exeDate, int exp, string& expDate)
{
	if(exp == my_NAN) return my_NAN;       // exe expires never

	// compilation time of exe
	tm tm = to_tm(exeDate.c_str());                              // string to struct

	// time of expiration
	tm.tm_mday += exp;                                           // add the number of days until expiration
	time_t t_exp = mktime(&tm);                                  // struct to time
	char datetime[20];
	strftime(datetime, 20, "%b %d %Y", localtime(&t_exp));
	expDate = datetime;

	return (int)difftime(t_exp,  time(NULL))/60/60/24;           // return current days until expiration
}

tm to_tm(char const *time)
{
	char s_month[5];
	int month, day, year;
	struct tm t = {0};
	static const char month_names[] = "JanFebMarAprMayJunJulAugSepOctNovDec";

	sscanf(time, "%s %d %d", s_month, &day, &year);
	month = (int)(strstr(month_names, s_month)-month_names)/3;

	t.tm_mon = month;
	t.tm_mday = day;
	t.tm_year = year - 1900;
	t.tm_isdst = -1;
	return t;
}

// time_t to_time_t(char const *time){return mktime(&to_tm(time));}

// ----------------------------------------------------------------------------------------
// isEOL
// ----------------------------------------------------------------------------------------
/** is c the end of a line?
 * Windows: "\r\n"
 * Mac: 	"\r"
 * Linux: 	"\n"
 */
bool isEOL(char c)
{
	return (c=='\r' || c=='\n');
}

/** same as above, but removes the second character from the stream if present
 */
bool isEOL(char c, istream& is)
{
	if(c=='\r'){
		if(is.peek()=='\n') is.get(c);
		return true;
	}
	return (c=='\n');
}

// ----------------------------------------------------------------------------------------
// safeGetline
// ----------------------------------------------------------------------------------------
/* returns a line with whatever the end of line is characterized
 */
istream& safeGetline(istream& is, string& t)
{
    t.clear();

    // The characters in the stream are read one-by-one using a std::streambuf.
    // That is faster than reading them one-by-one using the std::istream.
    // Code that uses streambuf this way must be guarded by a sentry object.
    // The sentry object performs various tasks,
    // such as thread synchronization and updating the stream state.

    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();

    for(;;) {
        int c = sb->sbumpc();
        switch (c) {
        case '\r':
            c = sb->sgetc();
            if(c == '\n')
                sb->sbumpc();
            return is;
        case '\n':
        case EOF:
            return is;
        default:
            t += (char)c;
        }
    }
    return is; // not used but code analysis desires it
}

// ----------------------------------------------------------------------------------------
// safeIgnore
// ----------------------------------------------------------------------------------------
/* returns a line with whatever the end of line is characterized
 */
istream& safeIgnore(istream& is)
{
    // The characters in the stream are read one-by-one using a std::streambuf.
    // That is faster than reading them one-by-one using the std::istream.
    // Code that uses streambuf this way must be guarded by a sentry object.
    // The sentry object performs various tasks,
    // such as thread synchronization and updating the stream state.

    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();

    for(;;) {
        int c = sb->sbumpc();
        switch (c) {
        case '\r':
            c = sb->sgetc();
            if(c == '\n')
                sb->sbumpc();
            return is;
        case '\n':
        case EOF:
            return is;
        }
	}
    return is; // not used but code analysis desires it
}

// ----------------------------------------------------------------------------------------
// convertion between centiMorgan (cM) and r (recombination rate)
// ----------------------------------------------------------------------------------------
double cM2r(double cM){
    return (1-exp(-2*cM/100))/2;
}
double r2cM(double r){
    return -log(1-2*r)/2*100;
}

// ----------------------------------------------------------------------------------------
// set_logtime_occurences
// ----------------------------------------------------------------------------------------
/** return a vector of the loged time for the given logtime */
vector<unsigned int>
get_logtime_occurences(Param* pParam, unsigned int totGen, bool isCoal)
{
    vector<unsigned int> vIndex;
    // if there is not logtime parameter (coalescence...)
    if(!pParam){
        vIndex.push_back(totGen);
        return vIndex;
    }

    if(pParam && pParam->isTemporalParam()){
        map<int, string>* pArgs = pParam->get_temporal_args();
        map<int, string>::iterator pos = pArgs->begin();  // set the iterator
        
        // loop through the generations and sum up the occurences
        for(unsigned int occ, i=1; i<=totGen; ++i){
            assert(i!=1 || (i==1 && i==(unsigned int) pos->first)); // first one has to be set
            if(i == (unsigned int) pos->first){                      // new temporal parameter
                occ = strTo<unsigned int>(pos->second);     // change occurence
                ++pos;                                  // go the the next position
            }
            if(!((i) % occ)) vIndex.push_back(i); //++ _tot_occurrence;
        }
    }
    else {     // no temporal change
        unsigned int occ = (unsigned int) pParam->get_value();
        for(unsigned int i=1; i<=totGen; ++i){
            if(!((i) % occ)) vIndex.push_back(i); //++ _tot_occurrence;
        }
    }
    
    // coalescence simulations: stats
    if (isCoal && pParam->get_name()=="statistics" && (vIndex.size() != 1 || vIndex[0] != totGen)) {
        warning("Coaelscence: %s may only be computed for the last generation: logtime adjusted!\n",
                pParam->get_name().c_str());
        vIndex.clear(); // simplest way to achieve it, although not the most effective...
        vIndex.push_back(totGen);
    }
    
    return vIndex;
}

// ----------------------------------------------------------------------------------------
// runScript
// ----------------------------------------------------------------------------------------
/** launch a script with the filename as parameter */
void runScript(string scriptname, string filename)
{
    // check if the script exists
    ifstream script(scriptname.c_str(), ios::in);
    if (!script.is_open()) {
        error("could not open script %s!\n", scriptname.c_str());
        return;
    }
    script.close();
    
    string cmd = scriptname + " " + filename;
    
#ifdef _DEBUG
    message("Executing shell script \"%s\" ",cmd.c_str());
#endif
    fflush(stdout);
    
    if (system(cmd.c_str()) < 0) {
        error("Execution of `sh %s' failed: %s\n", cmd.c_str(), strerror(errno));
        return;
    }
    
#ifdef _DEBUG
    message("...done\n");
#endif
}

// ----------------------------------------------------------------------------------------
// set_path if it does not exist
// ----------------------------------------------------------------------------------------
void check_path(string& path)
{
#if defined(__BCPLUSPLUS__) && defined(_DEBUG)
    cout << "Warning C++Builder: path is not checked!!!" << endl;
#else
        if (path.size()) {
        if (path[path.length() - 1] != SEP) path += SEP;
        
        // for each folder
        string::size_type cur, next;
        string curFolder;
        cur = 0;
        next = path.find_first_of("/\\", cur);
        while (next != string::npos) {
            curFolder = path.substr(0, next + 1);
            if(next<=1){
                // go to next folder
                cur = next + 1;
                next = path.find_first_of("/\\", cur);
                continue;
            }
            DIR *dirname = opendir(curFolder.c_str());
            if (!dirname) {      // if folder does not exist
#ifdef  __BCPLUSPLUS__                 // if Borland is used
                if((mkdir(curFolder.c_str())) == -1) {
#else
                    string cmd = "mkdir -p \"" + curFolder + "\""; // "" are needed if a folder name contains spaces
                    if (system(cmd.c_str()) < 0) {
#endif
                        
                        error("could not create directory \"%s\", saving in wd.\n",
                              path.c_str());
                        path = "";
                        return;
                    }
                }
                else closedir(dirname);
                
                // go to next folder
                cur = next + 1;
                next = path.find_first_of("/\\", cur);
            }
        }
#endif
    }




