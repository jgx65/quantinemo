/** @file types.h
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

#ifndef typesH
#define typesH

//#define _RANDOM_C11


// use the assert macro only when debugging
// #define _DEBUG
#ifndef _DEBUG
#define NDEBUG  // deactivate assert
#endif


typedef unsigned short int POP_SIZE;
typedef unsigned short int PATCH_ID;
typedef unsigned short int MIGR_SIZE;
typedef unsigned short ALLELE;
#define PI 3.1415926535897932384626433832795028841972

#include <string>
#include <assert.h>

using namespace std;

/**Sex types, males are always 0 and females 1!!**/
typedef enum {
    MAL=0, FEM=1
} sex_t;

/**Array index of the age classes in the patch sizes and containers arrays.**/
typedef enum {
    OFFSx=0, ADLTx=1
}age_idx;

/**Age class flags.*/  // binary representation
typedef unsigned int age_t;
#define NONE 0        // 00000000   /** No age flag.*/
#define OFFSPRG 1     // 00000001   /** Offspring age class flag.*/
#define ADULTS 2      // 00000010   /** Adults age class flag (breeders).*/
#define ALL 3         // 00000011   /** Adults and offspring flag */

inline age_t age_idx2t(const age_idx& AGE){
    return (AGE == ADLTx ? ADULTS : OFFSPRG);
}

inline age_idx age_t2idx(const age_t& AGE){
    assert(AGE==ADULTS || AGE==OFFSPRG);
    return (AGE == ADULTS ? ADLTx : OFFSx);
}



#ifndef LONG_MAX
#define LONG_MAX 2147483647L
#endif
#ifndef ULONG_MAX
#define ULONG_MAX (LONG_MAX * 2UL + 1)
#endif


#define ploidy 2


#define my_NAN 99999
#define my_NANstr (string)("NaN")

#define my_STR 22222
#define my_STRstr (string)("STR")

#define SAMPLED 55555
#define RECOMB 66666

#define NB_AGE_CLASSES 2

#ifdef __BCPLUSPLUS__                 // if Borland is used
#define SEP '\\'
#else
#define SEP '/'
#endif


/**Ordering type used to record statistics in the StatRecorders.**/
typedef enum {
    FLAT = 2,       // values are stored for each generation and replicate
    GEN  = 4,       // values are stored for each generation: replicate values are added
    RPL  = 6,       // values are stored for each replicate: generation values are added
    PARAM = 8,
}st_order;

/**mutation models.**/
typedef enum {
    KAM,
    SSM,
    RMM,
    IMM,
    NO
}mut_model_t;

/**sequence initialization models.**/
typedef enum {
    INI_UNIF,
    INI_MONO,
    INI_DIST,
}ini_model_t;

/**Trait types**/
typedef string trait_t;
/**Max number of characters in the trait's type descriptor.*/
#define TRAIT_T_MAX 5
#define DELE "delet"
#define DISP "disp"
#define FDISP "fdisp"
#define MDISP "mdisp"
#define NTRL "ntrl"
#define DQUANT "quanti"

/**Param's types**/
typedef enum {
    DBL,INT2,STR,MAT,DIST,MAT_VAR,
    INT_MAT,DBL_MAT,STR_MAT
}param_t;



#endif

