/** @file main.cpp
 *
 *   Copyright (C) 2006 Frederic Guillaume    <guillaum@zoology.ubc.ca>
 *   Copyright (C) 2008 Samuel Neuenschwander <samuel.neuenschwander@unil.ch>
 *   Copyright (C) 2018 Frederic Michaud <frederic.a.michaud@gmail.com>
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
unsigned int nb_error;
unsigned int nb_warning;
unsigned int verbose_message;
unsigned int verbose_warning;
unsigned int verbose_error;

#include "version.h"
#include "tsim_manager.h"


#include "functions.h"

using namespace std;





//--------------------------------------------------------------------
int main (int argc, char **argv)
{
    int returnVal = 0;  // value to return by the program
    verbose_message = 1;
    verbose_warning = 1;
    verbose_error = 1;
    
	try{
		try{
			

#ifdef _DEBUG
            message("\n***** DEBUG MODE *****\n\n");

#endif
            // read settings file(s) ///////////////////////////////////////////
            TSimManager simManager(argc, argv);
            
            // run simulations /////////////////////////////////////////////////
            simManager.run();
        
            if(nb_warning) message("\nCaution: %i WARNINGs were raised. See above for details.", nb_warning);
			if(nb_error)   message("\nCaution: %i ERRORs were raised. See above for details.", nb_error);
			message("\nquantiNemo terminated successfully!\n");
        }
		catch(const bad_alloc& x){
			error("Out of memory (%s)!\n", x.what());
		}
		catch(const char* text){
			error("%s!\n", text);
		}
	}
	catch(const int value){
		if(nb_warning) message("\nCaution: %i WARNINGs were raised.", nb_warning);
		if(nb_error)   message("\nCaution: %i ERRORs were raised.", nb_error);
		if(value!=1111) message("\nquantiNemo was not able to run successfully!\n");  // error re-throws
		//else message("\nquantiNemo terminated!\n");
        returnVal=1;
        
	}
	catch(...){
		if(nb_warning) message("\nCaution: %i WARNINGs were raised.", nb_warning);
		if(nb_error)   message("\nCaution: %i ERRORs were raised.", nb_error);
		message("\nquantiNemo was not able to run successfully!\n");  // error re-throws
        returnVal=1;
	}
    
#ifdef _DEBUG
	message("\nquantiNemo debug: memory cleanup succeeded!\n");
#ifdef  __BCPLUSPLUS__	// this allows to stop the program before the console window is closed
	message("\ndebug C++Builder: type any character to finish and close the window:!\n");
    int w;
 	 cin >> w;
#endif
#else
    //int z;
    //cin >> z;
#endif
	return returnVal;
}










