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
typedef unsigned char ALLELE;
#define PI 3.1415926535897932384626433832795028841972

#include <string>
#include <assert.h>
#include <stdint.h>

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

class AlleleRef;

// ----------------------------------------------------------------------------------------
// AlleleContainer owns one individual's diploid allele block and is the single place
// that knows the storage layout. Two layouts exist and the choice is made once per run:
//
//   byte mode   (default)  one ALLELE per allele at _data[l*ploidy + c].
//                          The general case: supports up to 255 allelic states.
//   packed mode            one *bit* per allele at bit (l*ploidy + c) of _bits.
//                          Only valid when every locus is diallelic, which is checked
//                          once at initialisation (see TGenomeProto::ini_packed_mode).
//                          Uses 8x less memory for the allele block.
//
// The mode is a static because it is a property of the run and not of the individual:
// keeping it out of the object avoids growing the per-individual footprint, and lets
// the hot inheritance kernels be selected once (see the templated kernels in tgenome.cpp)
// instead of branching per allele access.
class AlleleContainer {
    static bool  _packed;      // run-wide layout choice; false = byte mode
    ALLELE*      _data;        // byte mode storage   (0 in packed mode)
    uint64_t*    _bits;        // packed mode storage (0 in byte mode)
    unsigned int _nb_locus;

    static size_t nb_words(unsigned int nb_locus){
        return ((size_t)nb_locus*ploidy + 63) / 64;
    }
public:
    // -- run-wide layout selection (must be set before any genome is allocated) --
    static void set_packed(bool p){ _packed = p; }
    static bool packed()          { return _packed; }

    AlleleContainer()                        : _data(0), _bits(0), _nb_locus(0) {}
    AlleleContainer(const AlleleContainer& o): _data(0), _bits(0), _nb_locus(0) { *this = o; }
    ~AlleleContainer()                       { delete[] _data; delete[] _bits; }
    AlleleContainer& operator=(const AlleleContainer& o){
        if(this != &o){
            if(_nb_locus != o._nb_locus){
                clear();
                allocate(o._nb_locus);
            }
            if(_packed){
                for(size_t w = 0, n = nb_words(_nb_locus); w < n; ++w) _bits[w] = o._bits[w];
            }
            else {
                for(size_t i = 0, n = (size_t)_nb_locus*ploidy; i < n; ++i) _data[i] = o._data[i];
            }
        }
        return *this;
    }
    void allocate(unsigned int nb_locus){
        delete[] _data; _data = 0;
        delete[] _bits; _bits = 0;
        _nb_locus = nb_locus;
        if(_packed) _bits = new uint64_t[nb_words(nb_locus)]();
        else        _data = new ALLELE[(size_t)nb_locus*ploidy];
    }
    void clear(){ delete[] _data; _data = 0; delete[] _bits; _bits = 0; _nb_locus = 0; }
    bool         allocated() const { return _packed ? _bits != 0 : _data != 0; }
    unsigned int nb_locus()  const { return _nb_locus; }

    // -- layout-agnostic value access (used everywhere outside the hot kernels) --
    inline ALLELE get(size_t locus, size_t copy) const {
        size_t i = locus*ploidy + copy;
        return _packed ? (ALLELE)((_bits[i>>6] >> (i&63)) & 1ULL) : _data[i];
    }
    inline void set(size_t locus, size_t copy, ALLELE v){
        size_t i = locus*ploidy + copy;
        if(_packed){
            uint64_t m = 1ULL << (i&63);
            if(v) _bits[i>>6] |=  m;
            else  _bits[i>>6] &= ~m;
        }
        else _data[i] = v;
    }

    // allele() keeps the original `seq.allele(l,c)` read/write syntax working by
    // returning a proxy; see AlleleRef below.
    inline AlleleRef allele(size_t locus, size_t copy);
    inline ALLELE    allele(size_t locus, size_t copy) const { return get(locus, copy); }

    // -- direct, unbranched accessors for the mode-specialised hot kernels --
    inline ALLELE byte_at(size_t i) const        { return _data[i]; }
    inline void   set_byte_at(size_t i, ALLELE v){ _data[i] = v; }
    inline ALLELE bit_at(size_t i) const         { return (ALLELE)((_bits[i>>6] >> (i&63)) & 1ULL); }
    inline void   set_bit_at(size_t i, ALLELE v){
        uint64_t m = 1ULL << (i&63);
        if(v) _bits[i>>6] |=  m;
        else  _bits[i>>6] &= ~m;
    }

    // -- per-locus kernels (mutation / initialisation) operate on a small scratch copy --
    // In byte mode locus_ptr() still hands out the real storage so those kernels work in
    // place exactly as before; in packed mode callers must round-trip through
    // read_locus()/write_locus() because a bit has no address.
    inline ALLELE*       locus_ptr(size_t locus)       { assert(!_packed); return _data + locus*ploidy; }
    inline const ALLELE* locus_ptr(size_t locus) const { assert(!_packed); return _data + locus*ploidy; }
    inline void read_locus (size_t locus, ALLELE* out) const {
        for(size_t c = 0; c < ploidy; ++c) out[c] = get(locus, c);
    }
    inline void write_locus(size_t locus, const ALLELE* in){
        for(size_t c = 0; c < ploidy; ++c) set(locus, c, in[c]);
    }
};

// Proxy standing in for an ALLELE lvalue. A packed allele is a single bit and therefore
// has no address, so allele() cannot return ALLELE&; this restores assignment and
// implicit read at the call sites that used the reference.
class AlleleRef {
    AlleleContainer* _c;
    size_t           _locus, _copy;
public:
    AlleleRef(AlleleContainer* c, size_t l, size_t cp): _c(c), _locus(l), _copy(cp) {}
    operator ALLELE() const { return _c->get(_locus, _copy); }
    AlleleRef& operator=(ALLELE v)          { _c->set(_locus, _copy, v); return *this; }
    AlleleRef& operator=(const AlleleRef& o){ _c->set(_locus, _copy, (ALLELE)o); return *this; }
    AlleleRef& operator++()                 { _c->set(_locus, _copy,
                                                      (ALLELE)(_c->get(_locus, _copy) + 1)); return *this; }
};

inline AlleleRef AlleleContainer::allele(size_t locus, size_t copy){
    return AlleleRef(this, locus, copy);
}

// Access policies for the hot inheritance/recombination kernels. Those kernels are
// templated on one of these and the right instantiation is bound once per run to the
// existing _inherit_func_ptr/_recombine_func_ptr, so the inner loop never tests the
// layout: the byte instantiation compiles to exactly the code that ran before.
struct ByteAccess {
    static inline ALLELE get(const AlleleContainer& c, size_t i)        { return c.byte_at(i); }
    static inline void   set(AlleleContainer& c, size_t i, ALLELE v)    { c.set_byte_at(i, v); }
};
struct PackedAccess {
    static inline ALLELE get(const AlleleContainer& c, size_t i)        { return c.bit_at(i); }
    static inline void   set(AlleleContainer& c, size_t i, ALLELE v)    { c.set_bit_at(i, v); }
};


#define my_NAN 9999999
#define my_NANstr (string)("NaN")

#define my_STR 2222222
#define my_STRstr (string)("STR")

#define SAMPLED 5555555
#define RECOMB 6666666

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

