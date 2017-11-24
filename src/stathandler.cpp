/** @file stathandler.cpp
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

#include "stathandler.h"
#include "tmetapop.h"
#include <limits>

// ------------------------------------------------------------------------------
/** add the stat for each pair of patches */
template<class SH>
bool StatHandler<SH>::add_pairwisePatch(const string & end,
                                        const string & Title, const string & Name, const st_order & Order,
                                        const age_t & AGE, const unsigned int&ARG, double(SH::*getStat)(void),
                                        double(SH::*getStatBoolArg)(bool), double(SH::*getStatUintArg)(unsigned int),
                                        double(SH::*getStatAGE)(const age_t&), double(SH::*getStatBoolArgAGE)(bool,
                                                                                                              const age_t&), double(SH::*getStatUintArgAGE)(unsigned int, const age_t&))
{
    assert(!ARG && (getStatUintArg || getStatUintArgAGE));
    
    // check if the pair is specified
    if (end != "pair")
        return add_pairwisePatch_ij(end, Title, Name, Order, AGE, ARG, getStat,
                                    getStatBoolArg, getStatUintArg, getStatAGE, getStatBoolArgAGE,
                                    getStatUintArgAGE);
    
    // all pairs have to be computed
    unsigned int nbPatch = get_nbPatch();
    if (nbPatch == 1)
        return add(Title, Name, Order, AGE, ARG, getStat, getStatBoolArg,
                   getStatUintArg, getStatAGE, getStatBoolArgAGE, getStatUintArgAGE);
    
    // several patches are available
    string t1, t2, n1, n2;
    unsigned int i, j;
    
    for (i = 0; i < nbPatch - 1; ++i) {
        if (get_vPatch(i)->get_sampleID() == my_NAN) continue;
        n1 = "_p" + toStr(i + 1, nbPatch) + "-"; // "_p1-2"   first patch (name)
        t1 = " between patch " + toStr(i + 1, (unsigned int)1) + " and ";
        // first patch (description)
        for (j = i + 1; j < nbPatch; ++j) {
            if (get_vPatch(j)->get_sampleID() == my_NAN) continue;
            n2 = toStr(j + 1, nbPatch); // second patch (name)
            t2 = toStr(j + 1); // second patch (description)
            add(Title + t1 + t2, Name + n1 + n2, Order, AGE, toID(i, j, nbPatch, nbPatch),
                getStat, getStatBoolArg, getStatUintArg, getStatAGE,
                getStatBoolArgAGE, getStatUintArgAGE);
        }
    }
    
    return true;
}

// ------------------------------------------------------------------------------
/** add the stat for the specified pair of patches */
template<class SH>
bool StatHandler<SH>::add_pairwisePatch_ij(const string & end, const string & Title, const string & Name,
                                           const st_order & Order, const age_t & AGE, const unsigned int&ARG,
                                           double(SH::*getStat)(void), double(SH::*getStatBoolArg)(bool),
                                           double(SH::*getStatUintArg)(unsigned int),
                                           double(SH::*getStatAGE)(const age_t&),
                                           double(SH::*getStatBoolArgAGE)(bool, const age_t&),
                                           double(SH::*getStatUintArgAGE)(unsigned int, const age_t&))
{
    assert(!ARG && (getStatUintArg || getStatUintArgAGE));
    if (end[0] != 'p')return false;
    
    // extract the patch ids
    string::size_type pos = end.find('-'); // get the separation of the two loci
    if (pos == string::npos) return false; // no "-" found
    unsigned int p1, p2, nbPatch = _popPtr->getPatchNbr();
    try {
        p1 = strTo<unsigned int>(end.substr(1, pos - 1)); // first patch id
        p2 = strTo<unsigned int>(end.substr(pos + 1)); // second patch id
    }
    catch(...) {
        return false;
    } // this parameter may not be set for each trait separately
    if (p1 > nbPatch) return false; // id of patch exeeds the range
    if (p2 > nbPatch) return false; // id of patch exeeds the range
    if (p1 == p2)     return false; // that is not a pair of patches
    
    if (p1 > p2) swap(p1, p2); // if the order of the patches is wrong, swap it
    
    // check if the patches are really sampled
    if (get_vPatch(p1 - 1)->get_sampleID() == my_NAN) {
        warning("Stat '%s' cannot be computed, since the patch %i is not sampled!\n", Name.c_str(), p1);
        return false;
    }
    if (get_vPatch(p2 - 1)->get_sampleID() == my_NAN) {
        warning("Stat '%s' cannot be computed, since the patch %i is not sampled!\n", Name.c_str(), p2);
        return false;
    }
    
    add(Title + " between patches " + toStr(p1) + " and " + toStr(p2), // description
        Name + "_p" + toStr(p1, nbPatch) + "-" + toStr(p2, nbPatch), // name
        Order, AGE, (p1 - 1) + (p2 - 1)*nbPatch, getStat, getStatBoolArg,
        getStatUintArg, getStatAGE, getStatBoolArgAGE, getStatUintArgAGE);
    
    return true;
}

// ------------------------------------------------------------------------------
/** add the stat for each pop */
template<class SH>
bool StatHandler<SH>::add_perPatch(const string & end, const string & Title, const string & Name, const st_order & Order,
                                   const age_t & AGE, const unsigned int&ARG, double(SH::*getStat)(void),
                                   double(SH::*getStatBoolArg)(bool), double(SH::*getStatUintArg)(unsigned int),
                                   double(SH::*getStatAGE)(const age_t&), double(SH::*getStatBoolArgAGE)(bool,const age_t&),
                                   double(SH::*getStatUintArgAGE)(unsigned int, const age_t&))
{
    assert(!ARG && (getStatUintArg || getStatUintArgAGE));
    
    // check if the patch is specified (just an "p")
    if (end != "p"){
        return add_perPatch_i(end, Title, Name, Order, AGE, ARG, getStat,
                              getStatBoolArg, getStatUintArg, getStatAGE, getStatBoolArgAGE,
                              getStatUintArgAGE);
    }
    
    // the stat will be computed for each pop
    unsigned int nbPatch = get_nbPatch();
    if (nbPatch == 1){
        return add(Title, Name, Order, AGE, ARG, getStat, getStatBoolArg,
                   getStatUintArg, getStatAGE, getStatBoolArgAGE, getStatUintArgAGE);
    }
    // several patches are available
    string t, n;
    
    for (unsigned int i = 0; i < nbPatch; ++i) {
        if (get_vPatch(i)->get_sampleID() == my_NAN) continue;
        n = "_p" + toStr(i + 1, nbPatch); // name: "_p1"
        t = " of population " + toStr(i + 1); // description
        add(Title + t, Name + n, Order, AGE, i, getStat, getStatBoolArg,
            getStatUintArg, // pass the patch id
            getStatAGE, getStatBoolArgAGE, getStatUintArgAGE);
    }

    return true;
}

// ------------------------------------------------------------------------------
/** add the stat for the specified pop */
template<class SH>
bool StatHandler<SH>::add_perPatch_i(const string & end,
                                     const string & Title, const string & Name, const st_order & Order,
                                     const age_t & AGE, const unsigned int&ARG, double(SH::*getStat)(void),
                                     double(SH::*getStatBoolArg)(bool), double(SH::*getStatUintArg)(unsigned int),
                                     double(SH::*getStatAGE)(const age_t&), double(SH::*getStatBoolArgAGE)(bool,const age_t&),
                                     double(SH::*getStatUintArgAGE)(unsigned int, const age_t&))
{
 if (end[0] != 'p') return false;
    
    // extract the locus
    unsigned int p, nbPatch = _popPtr->getPatchNbr();
    try {
        p = strTo<unsigned int>(end.substr(1));
    }
    catch(...) {
        return false;
    } // this parameter may not be set for each patch separately
    if (p > nbPatch) return false; // id of patch exeeds the range
    
    // check if the pop is really sampled
    if (get_vPatch(p - 1)->get_sampleID() == my_NAN) {
        warning("Stat '%s' cannot be computed, since the patch %i is not sampled!\n", Name.c_str(), p);
        return false;
    }
    
    add(Title + " of population " + toStr(p), // description
        Name + "_p" + toStr(p, nbPatch), // name
        Order, AGE, p - 1, getStat, getStatBoolArg, getStatUintArg, // patch id
        getStatAGE, getStatBoolArgAGE, getStatUintArgAGE);

    return true;
}

// ------------------------------------------------------------------------------
/** add the stat for each pop */
template<class SH>
bool StatHandler<SH>::add_perPatch_andLocus(const string & end, const string & Title, const string & Name,
                                            const st_order & Order, const age_t & AGE, const unsigned int&ARG,
                                            double(SH::*getStat)(void), double(SH::*getStatBoolArg)(bool),
                                            double(SH::*getStatUintArg)(unsigned int),
                                            double(SH::*getStatAGE)(const age_t&),
                                            double(SH::*getStatBoolArgAGE)(bool, const age_t&),
                                            double(SH::*getStatUintArgAGE)(unsigned int, const age_t&))
{
 assert(!ARG && (getStatUintArg || getStatUintArgAGE));
    assert(end == "p_l" || end == "l_p");
    
    // check if therer are several patches and loci
    unsigned int nbPatch = get_nbPatch();
    if (nbPatch == 1) {
        if (_nb_locus == 1){
            return add(Title, Name, Order, AGE, ARG, getStat, getStatBoolArg,
                       getStatUintArg, // single patch and sinlge locus
                       getStatAGE, getStatBoolArgAGE, getStatUintArgAGE);
        }
        else{
            return add_perLocus("l", Title, Name, Order, AGE, ARG, getStat,
                                getStatBoolArg, getStatUintArg, // single patch
                                getStatAGE, getStatBoolArgAGE, getStatUintArgAGE);
        }
    }
    else if (_nb_locus == 1){
        return add_perPatch("p", Title, Name, Order, AGE, ARG, getStat,
                            getStatBoolArg, getStatUintArg, // single locus
                            getStatAGE, getStatBoolArgAGE, getStatUintArgAGE);
    }
    
    string t1, t2, n1, n2;
    unsigned int i, l;
    
    for (i = 0; i < nbPatch; ++i) {
        if (get_vPatch(i)->get_sampleID() == my_NAN) continue;
        n1 = "_p" + toStr(i + 1, nbPatch) + "_l"; // "_p1_l1"   patch (name)
        t1 = " of patch " + toStr(i + 1, (unsigned int)1) + " and locus ";
        // locus (description)
        for (l = 0; l < _nb_locus; ++l) {
            n2 = toStr(l + 1, _nb_locus); // locus (name)
            t2 = toStr(l + 1); // locus (description)
            add(Title + t1 + t2, Name + n1 + n2, Order, AGE, toID(i, l, _nb_patch, _nb_locus),
                getStat, getStatBoolArg, getStatUintArg, getStatAGE,
                getStatBoolArgAGE, getStatUintArgAGE);
        }
    }
    
    return true;
}

// ------------------------------------------------------------------------------
/** add the stat for the each locus */
template<class SH>
bool StatHandler<SH>::add_perLocus(const string & end,const string & Title, const string & Name, const st_order & Order,
                                   const age_t & AGE, const unsigned int&ARG, double(SH::*getStat)(void),
                                   double(SH::*getStatBoolArg)(bool), double(SH::*getStatUintArg)(unsigned int),
                                   double(SH::*getStatAGE)(const age_t&), double(SH::*getStatBoolArgAGE)(bool,
                                                                                                         const age_t&), double(SH::*getStatUintArgAGE)(unsigned int, const age_t&))
{

 assert(!ARG && (getStatUintArg || getStatUintArgAGE));
    
    // check if the locus is specified (just an "l")
    if (end != "l"){
        return add_perLocus_i(end, Title, Name, Order, AGE, ARG, getStat,
                              getStatBoolArg, getStatUintArg, getStatAGE, getStatBoolArgAGE,
                              getStatUintArgAGE);
    }
    
    // the stat will be computed for each locus
    if (_nb_locus == 1){
        return add(Title, Name, Order, AGE, ARG, getStat, getStatBoolArg,
                   getStatUintArg, getStatAGE, getStatBoolArgAGE, getStatUintArgAGE);
    }
    
    // several patches are available
    string t, n;
    
    for (unsigned int i = 0; i < _nb_locus; ++i) {
        n = "_l" + toStr(i + 1, _nb_locus); // name: "_l1"
        t = " of locus " + toStr(i + 1); // description
        add(Title + t, Name + n, Order, AGE, i, getStat, getStatBoolArg,
            getStatUintArg, getStatAGE, getStatBoolArgAGE, getStatUintArgAGE);
    }
    
    return true;
}

// ------------------------------------------------------------------------------
/** add the stat for the specified locus */
template<class SH>
bool StatHandler<SH>::add_perLocus_i(const string & end,
                                     const string & Title, const string & Name, const st_order & Order,
                                     const age_t & AGE, const unsigned int&ARG, double(SH::*getStat)(void),
                                     double(SH::*getStatBoolArg)(bool), double(SH::*getStatUintArg)(unsigned int),
                                     double(SH::*getStatAGE)(const age_t&), double(SH::*getStatBoolArgAGE)(bool,
                                                                                                           const age_t&), double(SH::*getStatUintArgAGE)(unsigned int, const age_t&))
{
 assert(!ARG && (getStatUintArg || getStatUintArgAGE));
    
    if (end[0] != 'l') return false;
    
    // extract the locus
    unsigned int l;
    
    try {
        l = strTo<unsigned int>(end.substr(1));
    }
    catch(...) {
        return false;
    } // this parameter may not be set for each locus separately
    if (l > this->_nb_locus) return false; // id of locus exeeds the range
    
    add(Title + " of locus " + toStr(l), // description
        Name + "_l" + toStr(l, _nb_locus), // name
        Order, AGE, (l - 1), getStat, getStatBoolArg, getStatUintArg, getStatAGE,
        getStatBoolArgAGE, getStatUintArgAGE);

    return true;
}

// ------------------------------------------------------------------------------
/** add the stat for each pair of loci */
template<class SH> bool
StatHandler<SH>::add_pairwiseLocus(const string & end,
                                   const string & Title, const string & Name, const st_order & Order,
                                   const age_t & AGE, const unsigned int&ARG, double(SH::*getStat)(void),
                                   double(SH::*getStatBoolArg)(bool), double(SH::*getStatUintArg)(unsigned int),
                                   double(SH::*getStatAGE)(const age_t&), double(SH::*getStatBoolArgAGE)(bool,const age_t&),
                                   double(SH::*getStatUintArgAGE)(unsigned int, const age_t&))
{
 assert(!ARG && (getStatUintArg || getStatUintArgAGE));
    
    // check if the pair is specified
    if (end[0] == 'l') return add_pairwiseLocus_ij(end, Title, Name, Order, AGE, ARG, getStat,
                                                   getStatBoolArg, getStatUintArg, getStatAGE, getStatBoolArgAGE, getStatUintArgAGE);
    
    // all pairs have to be computed
    if (this->_nb_locus == 1) return add(Title, Name, Order, AGE, ARG, getStat, getStatBoolArg,
                                         getStatUintArg, getStatAGE, getStatBoolArgAGE, getStatUintArgAGE);
    
    // several loci are available
    string t1, t2, n1, n2;
    unsigned int i, j;
    
    for (i = 0; i < _nb_locus; ++i) {
        n1 = "_l" + toStr(i + 1, _nb_locus) + "-"; // "_l1-2"   first locus (name)
        t1 = " between loci " + toStr(i + 1, (unsigned int)1) + " and ";
        // first locus (description)
        for (j = i + 1; j < _nb_locus; ++j) {
            n2 = toStr(j + 1, _nb_locus); // second locus (name)
            t2 = toStr(i + 1); // second locus (description)
            add(Title + t1 + t2, Name + n1 + n2, Order, AGE, toID(i, j, _nb_locus, _nb_locus),
                getStat, getStatBoolArg, getStatUintArg, getStatAGE, getStatBoolArgAGE,
                getStatUintArgAGE);
        }
    }

    return true;
}

// ------------------------------------------------------------------------------
/** add the stat for the specified pair of loci */
template<class SH> bool
StatHandler<SH>::add_pairwiseLocus_ij(const string & end, const string & Title, const string & Name,
                                      const st_order & Order, const age_t & AGE, const unsigned int&ARG,
                                      double(SH::*getStat)(void), double(SH::*getStatBoolArg)(bool),
                                      double(SH::*getStatUintArg)(unsigned int),
                                      double(SH::*getStatAGE)(const age_t&),
                                      double(SH::*getStatBoolArgAGE)(bool, const age_t&),
                                      double(SH::*getStatUintArgAGE)(unsigned int, const age_t&))
{
 assert(!ARG && (getStatUintArg || getStatUintArgAGE));
    if (end[0] != 'l')
        return false;
    
    // extract the loci
    string::size_type pos = end.find("-"); // get the separation of the two loci
    if (pos == string::npos) return false;
    // no "-" found
    unsigned int l1, l2;
    
    try {
        l1 = strTo<unsigned int>(end.substr(1, pos - 1)); // first locus
        l2 = strTo<unsigned int>(end.substr(pos + 1)); // second locus
    }
    catch(...) {
        return false;
    } // this parameter may not be set for each trait separately
    if (l1 > this->_nb_locus) return false; // id of locus exeeds the range
    if (l2 > this->_nb_locus) return false; // id of locus exeeds the range
    if (l1 == l2)             return false; // that is not a pair of loci
    
    if (l1 > l2) swap(l1, l2); // if the order of the loci is wrong, swap it
    
    return add(Title + " between loci " + toStr(l1) + " and " + toStr(l2), // description
               Name + "_l" + toStr(l1, _nb_locus) + "-" + toStr(l2, _nb_locus), // name
               Order, AGE, toID(l1-1, l2-1, _nb_locus, _nb_locus), getStat, getStatBoolArg,
               getStatUintArg, getStatAGE, getStatBoolArgAGE, getStatUintArgAGE);

    return true;
 }

// ------------------------------------------------------------------------------
/** add the stat for each pair of loci for each patch separately */
template<class SH> bool
StatHandler<SH>::add_pairwiseLocus_perPatch(const string & end,
                                            const string & Title, const string & Name, const st_order & Order,
                                            const age_t & AGE, const unsigned int&ARG, double(SH::*getStat)(void),
                                            double(SH::*getStatBoolArg)(bool), double(SH::*getStatUintArg)(unsigned int),
                                            double(SH::*getStatAGE)(const age_t&), double(SH::*getStatBoolArgAGE)(bool,const age_t&),
                                            double(SH::*getStatUintArgAGE)(unsigned int, const age_t&))
{
    assert(!ARG && (getStatUintArg || getStatUintArgAGE));
    
    // check if the pair is specified
    if (end.empty()) return add_pairwiseLocus_ij(end, Title, Name, Order, AGE, ARG, getStat,
                                                 getStatBoolArg, getStatUintArg, getStatAGE, getStatBoolArgAGE, getStatUintArgAGE);
    
    // all pairs have to be computed
    if (this->_nb_locus == 1) return false; // two loci are needed
    
    // the stat will be computed for each pop
    unsigned int nbPatch = get_nbPatch();
    if (nbPatch == 1) return add_pairwiseLocus(end,Title, Name, Order, AGE, ARG, getStat, getStatBoolArg,
                                               getStatUintArg, getStatAGE, getStatBoolArgAGE, getStatUintArgAGE);
    
    // several loci are available
    string t, t1, t2, n, n1, n2;
    unsigned int p, i, j;
    
    for (p = 0; p < nbPatch; ++p) {
        if (get_vPatch(p)->get_sampleID() == my_NAN) continue;
        n = "_p" + toStr(p + 1, nbPatch); // name: "_p1"
        t = " of population " + toStr(p + 1); // description
        for (i = 0; i < _nb_locus; ++i) {
            n1 = "_l" + toStr(i + 1, _nb_locus) + "-"; // "_l1-2"   first locus (name)
            t1 = " between loci " + toStr(i + 1, (unsigned int)1) + " and "; // first locus (description)
            for (j = i + 1; j < _nb_locus; ++j) {
                n2 = toStr(j + 1, _nb_locus); // second locus (name)
                t2 = toStr(i + 1);            // second locus (description)
                add(Title + t + t1 + t2, Name + n + n1 + n2, Order, AGE,
                    toID(p, i, j, _nb_patch, _nb_locus, _nb_locus),
                    getStat, getStatBoolArg, getStatUintArg, getStatAGE, getStatBoolArgAGE,
                    getStatUintArgAGE);
            }
        }
    }
    
    return true;
}

// ------------------------------------------------------------------------------
/** add the stat for each pair of loci for each patch separately */
template<class SH> bool
StatHandler<SH>::add_pairwiseLocus_perPatch_p(const string & end,
                                              const string & Title, const string & Name, const st_order & Order,
                                              const age_t & AGE, const unsigned int&ARG, double(SH::*getStat)(void),
                                              double(SH::*getStatBoolArg)(bool), double(SH::*getStatUintArg)(unsigned int),
                                              double(SH::*getStatAGE)(const age_t&), double(SH::*getStatBoolArgAGE)(bool,const age_t&),
                                              double(SH::*getStatUintArgAGE)(unsigned int, const age_t&))
{
    assert(!ARG && (getStatUintArg || getStatUintArgAGE));
    
    if (end[0] != 'p') return false;
    
    // extract the patch
    unsigned int p, nbPatch = _popPtr->getPatchNbr();
    try {
        p = strTo<unsigned int>(end.substr(1));
    }
    catch(...) {
        return false;
    } // this parameter may not be set for each patch separately
    if (p > nbPatch) return false; // id of patch exeeds the range
    
    // check if the pop is really sampled
    if (get_vPatch(p - 1)->get_sampleID() == my_NAN) {
        warning("Stat '%s' cannot be computed, since the patch %i is not sampled!\n", Name.c_str(), p);
        return false;
    }
    
    // all pairs have to be computed
    if (this->_nb_locus == 1) return false; // two loci are needed
    
    // several loci are available
    string t, t1, t2, n, n1, n2;
    unsigned int i, j;
    
    n = "_p" + toStr(p, nbPatch); // name: "_p1"
    t = " of population " + toStr(p); // description
    for (i = 0; i < _nb_locus; ++i) {
        n1 = "_l" + toStr(i + 1, _nb_locus) + "-"; // "_l1-2"   first locus (name)
        t1 = " between loci " + toStr(i + 1, (unsigned int)1) + " and "; // first locus (description)
        for (j = i + 1; j < _nb_locus; ++j) {
            n2 = toStr(j + 1, _nb_locus); // second locus (name)
            t2 = toStr(i + 1);            // second locus (description)
            add(Title + t1 + t2, Name + n1 + n2, Order, AGE,
                toID(p-1, i, j, _nb_patch, _nb_locus, _nb_locus),
                getStat, getStatBoolArg, getStatUintArg, getStatAGE, getStatBoolArgAGE,
                getStatUintArgAGE);
        }
    }
    
    return true;
}

// ------------------------------------------------------------------------------
/** add the stat for each pair of loci for each patch separately */
template<class SH> bool
StatHandler<SH>::add_pairwisePatch_perLocus(const string & end,
                                            const string & Title, const string & Name, const st_order & Order,
                                            const age_t & AGE, const unsigned int&ARG, double(SH::*getStat)(void),
                                            double(SH::*getStatBoolArg)(bool), double(SH::*getStatUintArg)(unsigned int),
                                            double(SH::*getStatAGE)(const age_t&), double(SH::*getStatBoolArgAGE)(bool,const age_t&),
                                            double(SH::*getStatUintArgAGE)(unsigned int, const age_t&))
{
 assert(!ARG && (getStatUintArg || getStatUintArgAGE));
    
    // check if the pair is specified
    if (end.empty()) return add_pairwisePatch_ij(end, Title, Name, Order, AGE, ARG, getStat,
                                                 getStatBoolArg, getStatUintArg, getStatAGE, getStatBoolArgAGE, getStatUintArgAGE);
    
    // all pairs have to be computed
    unsigned int nbPatch = get_nbPatch();
    if (nbPatch == 1) return false; // two patches are needed
    
    // the stat will be computed for each locus
    if (this->_nb_locus == 1) return add_pairwisePatch(end,Title, Name, Order, AGE, ARG, getStat, getStatBoolArg,
                                                       getStatUintArg, getStatAGE, getStatBoolArgAGE, getStatUintArgAGE);
    
    // several loci are available
    string t, t1, t2, n, n1, n2;
    unsigned int p1, p2, l;
    
    for (p1 = 0; p1 < nbPatch; ++p1) {
        if (get_vPatch(p1)->get_sampleID() == my_NAN) continue;
        n = "_p" + toStr(p1 + 1, nbPatch); // name: "_p1"
        t = " between populations " + toStr(p1 + 1) + " and "; // description
        for(p2 = p1+1; p2 < nbPatch; ++p2) {
            if (get_vPatch(p2)->get_sampleID() == my_NAN) continue;
            n1 = "-p" + toStr(p2 + 1, nbPatch) + "_l"; // "_p1-2"   second patch (name)
            t1 =  toStr(p2 + 1) + " of locus "; // first locus (description)
            for (l = 0; l < _nb_locus; ++l) {
                n2 = toStr(l + 1, _nb_locus); // locus (name)
                t2 = toStr(l + 1);            // locus (description)
                add(Title + t1 + t1 + t2, Name + n + n1 + n2, Order, AGE,
                    toID(p1, p2, l, _nb_patch, _nb_patch, _nb_locus),
                    getStat, getStatBoolArg, getStatUintArg, getStatAGE, getStatBoolArgAGE,
                    getStatUintArgAGE);
            }
        }
    }
    
    return true;
}

// ------------------------------------------------------------------------------
/** add the stat for each pair of loci for each patch separately */
template<class SH> bool
StatHandler<SH>::add_pairwisePatch_perLocus_l(const string & end,
                                              const string & Title, const string & Name, const st_order & Order,
                                              const age_t & AGE, const unsigned int&ARG, double(SH::*getStat)(void),
                                              double(SH::*getStatBoolArg)(bool), double(SH::*getStatUintArg)(unsigned int),
                                              double(SH::*getStatAGE)(const age_t&), double(SH::*getStatBoolArgAGE)(bool,const age_t&),
                                              double(SH::*getStatUintArgAGE)(unsigned int, const age_t&))
{
    assert(!ARG && (getStatUintArg || getStatUintArgAGE));
    
    if (end[0] != 'l') return false;
    
    // extract the patch
    unsigned int l, nbPatch = _popPtr->getPatchNbr();
    try {
        l = strTo<unsigned int>(end.substr(1));
    }
    catch(...) {
        return false;
    } // this parameter may not be set for each patch separately
    if (l > this->_nb_locus) return false; // id of locus exeeds the range
    
    // all pairs have to be computed
    if (nbPatch == 1) return false; // two patches are needed
    
    // several loci are available
    string t, t1, n, n1;
    unsigned int p1, p2;
    
    for (p1 = 0; p1 < nbPatch; ++p1) {
        n = "_p" + toStr(p1+1, nbPatch); // name: "_p1"
        t = " between populations " + toStr(p1+1) + " and "; // description
        for (p2 = p1+1; p2 < nbPatch; ++p2) {
            n1 = "-p" + toStr(p2 + 1, nbPatch) + "_l" + toStr(l, _nb_locus); // "_p1-2_l1"
            t1 = toStr(p2 + 1) + " of locus " + toStr(l); // second patch (description)
            add(Title + t + t1, Name + n + n1, Order, AGE,
                toID(p1, p2, l-1, _nb_patch, _nb_patch, _nb_locus),
                getStat, getStatBoolArg, getStatUintArg, getStatAGE, getStatBoolArgAGE,
                getStatUintArgAGE);
        }
    }
    
    return true;
}

// ------------------------------------------------------------------------------
/** Adds a StatRecorder to the list, it is also added to the StatHandlerBase::_stats list.
 Two types of function variables are passed to this function. The "getter" and "setter".
 A "getter" returns a double value that will be stored in the StatRecorder structure. It
 may or may not take an argument. Only one getter should be passed to the stat recorder.
 The setter is used to set variables in the SH class which are then read by the getter.
 A setter and a getter may be given together but a setter alone will issue an error at
 runtime as no getter is present in the stat recorder.
 * @param Title 						the stat title
 * @param Name 							the stat name (headers in the output file)
 * @param Order 						stat table ordering flag
 * @param AGE 							age on which the stat should be processed
 * @param ARG 							the argument to pass to the S function
 * @param getSt 						function ptr to a S getter
 * @param getStBoolArg 			function ptr to a S getter with boolean argument
 * @param getStUintArg 			function ptr to a S getter with unsigned int argument
 * @param getStAGE 					function ptr to a S getter with age argument
 * @param getStBoolArgAGE 	function ptr to a S getter with boolean and age argument
 * @param getStUintArgAGE		function ptr to a S getter with unsigned int and age argument
 */
template<class SH>
bool StatHandler<SH>::add(const string & Title, string Name,
                          const st_order & Order, const age_t & AGE, const unsigned int&ARG,
                          double(SH::*getStat)(void), double(SH::*getStatBoolArg)(bool),
                          double(SH::*getStatUintArg)(unsigned int),
                          double(SH::*getStatAGE)(const age_t&),
                          double(SH::*getStatBoolArgAGE)(bool, const age_t&),
                          double(SH::*getStatUintArgAGE)(unsigned int, const age_t&))
{
    StatRecorder<SH> *new_rec= new StatRecorder<SH>(_popPtr->get_pStat_db());
    if(_trait_index) Name += "_t" + toStr(_trait_index); // "_t2"
    new_rec->set(Title, Name, Order, AGE, ARG, getStat, getStatBoolArg,
                 getStatUintArg, getStatAGE, getStatBoolArgAGE, getStatUintArgAGE);
    
    StatHandlerBase::add(new_rec);
    _recorders.push_back(new_rec);
    new_rec->init_stat_rec(_popPtr);
    
#ifdef _DEBUG
    if (_recorders.size() > 1) message(" ");
    message("%s", Name.c_str());
#endif
    
    return true;
}

// ------------------------------------------------------------------------------
template<class SH>
StatHandler<SH>::~StatHandler()
{
    ARRAY::delete_2D(_hsnei_locus, NB_AGE_CLASSES);
    ARRAY::delete_2D(_htnei_locus, NB_AGE_CLASSES);
    ARRAY::delete_2D(_fst_locus, NB_AGE_CLASSES);
    ARRAY::delete_2D(_fis_locus, NB_AGE_CLASSES);
    ARRAY::delete_2D(_fit_locus, NB_AGE_CLASSES);
    ARRAY::delete_2D(_fst_WC_locus, NB_AGE_CLASSES);
    ARRAY::delete_2D(_fis_WC_locus, NB_AGE_CLASSES);
    ARRAY::delete_2D(_fit_WC_locus, NB_AGE_CLASSES);
    ARRAY::delete_2D(_ho_locus, NB_AGE_CLASSES);
    ARRAY::delete_2D(_hs_locus, NB_AGE_CLASSES);
    ARRAY::delete_2D(_ht_locus, NB_AGE_CLASSES);
    ARRAY::delete_2D(_het0_locus, NB_AGE_CLASSES);
    ARRAY::delete_2D(_het1_locus, NB_AGE_CLASSES);
    ARRAY::delete_2D(_fst_bn_locus, NB_AGE_CLASSES);
    
    if(_popPtr){
        if(!get_current_nbSamplePatch()) _popPtr->_current_nbSamplePatch = get_last_nbSamplePatch(); // since the stats were not made at the last generation
        ARRAY::delete_3D(_fst_matrix_wc, NB_AGE_CLASSES, get_current_nbSamplePatch());
        ARRAY::delete_3D(_fst_matrix, NB_AGE_CLASSES, get_current_nbSamplePatch());
        ARRAY::delete_3D(_coa_matrix, NB_AGE_CLASSES, get_current_nbSamplePatch());
        ARRAY::delete_3D(_alleleFreq_local, NB_AGE_CLASSES, get_current_nbSamplePatch());
        ARRAY::delete_2D(_alleleFreq_global, NB_AGE_CLASSES);
        ARRAY::delete_3D(_locusFreq_local, NB_AGE_CLASSES, get_current_nbSamplePatch());
        ARRAY::delete_2D(_locusFreq_global, NB_AGE_CLASSES);
    }
    
    for (REC_IT pos = _recorders.begin(); pos != _recorders.end(); ++pos) {
        delete *pos;
    }
    _recorders.clear();
    
    ARRAY::delete_2D(_computed, _computed_size);
}

// ------------------------------------------------------------------------------
/** stats are stored in the db	*/
template<class SH> void
StatHandler<SH>::execute()
{
    for (REC_IT rec = _recorders.begin(); rec != _recorders.end(); rec++) {
        (*rec)->setStats(_popPtr->getCurrentAge(), _popPtr->getCurrentGeneration(),
                         _popPtr->getCurrentReplicate(), dynamic_cast<SH*>(this), _popPtr);
    }
}

// ------------------------------------------------------------------------------
/** stats are directly written to the file */
template<class SH> void
StatHandler<SH>::execute(ostream& FH)
{
    if(_recorders.empty()) return;
    
    unsigned int curGen = _popPtr->getCurrentGeneration();
    unsigned int curRep = _popPtr->getCurrentReplicate();
    
    for (REC_IT rec = _recorders.begin(); rec != _recorders.end(); ++rec) {
        (*rec)->setStats(_popPtr->getCurrentAge(), curGen, curRep, dynamic_cast<SH*>(this), FH);
    }
}

// ------------------------------------------------------------------------------
/** stats are stored in the db	*/
template<class SH> void
StatHandler<SH>::execute_param()
{
    for (REC_IT rec = _recorders.begin(); rec != _recorders.end(); rec++) {
        (*rec)->setParams(_popPtr->getCurrentAge(), _popPtr->getCurrentGeneration(),
                         _popPtr->getCurrentReplicate(), dynamic_cast<SH*>(this), _popPtr);
    }
}

// ------------------------------------------------------------------------------
/** stats are directly written to the file */
template<class SH> void
StatHandler<SH>::execute_param(ostream& FH)
{
    if(_recorders.empty()) return;
    
    unsigned int curGen = _popPtr->getCurrentGeneration();
    unsigned int curRep = _popPtr->getCurrentReplicate();
    
    for (REC_IT rec = _recorders.begin(); rec != _recorders.end(); ++rec) {
        (*rec)->setParams(_popPtr->getCurrentAge(), curGen, curRep, dynamic_cast<SH*>(this), FH);
    }
}

// ------------------------------------------------------------------------------
template<class SH>void
StatHandler<SH>::set(TTraitProto * TT)
{
    _trait_index = TT->get_trait_index();
    _SHLinkedTrait = TT;
    _SHLinkedTraitIndex = TT->get_absolute_index();
    _popPtr = TT->get_popPtr();
    _nb_locus = _SHLinkedTrait->get_nb_locus();
    _nb_allele = _SHLinkedTrait->get_nb_allele();     // redirect the pointer
    _nb_patch = _SHLinkedTrait->get_popPtr()->get_nbPatch();
    _nb_allele_max = _SHLinkedTrait->get_nb_allele_max();
}

// ------------------------------------------------------------------------------
template<class SH> bool
StatHandler<SH>::init()
{
    StatHandlerBase::init();
    return true;
}

// ------------------------------------------------------------------------------
/** check if it has already been computed: return true if it was already computed
 * the array has to have a length of three
 * if the age is not used it is set to NONE (0)
 */
template<class SH>
bool StatHandler<SH>::already_computed(unsigned int*array, const age_idx & AGE)
{
    if (array[0] == get_current_generation() && array[1] == get_current_replicate() &&
        ((AGE + 1) & array[2])) { // check if this age class has already been computed
        return true;
    }
    
    array[0] = get_current_generation();
    array[1] = get_current_replicate();
    array[2] |= AGE + 1;
    // add this age class: Caution: 0: none; 1: offs; 2: adlts; 3 both
    return false;
}

// ------------------------------------------------------------------------------
/** check if it has already been computed: return true if it was already computed
 * the array has to have a length of three
 * if the age is not used it is set to NONE (0)
 */
template<class SH>
bool StatHandler<SH>::already_computed(unsigned int*array, age_t AGE)
{
    if (array[0] == get_current_generation() && array[1] == get_current_replicate() &&
        (AGE & array[2])) { // check if these age classes have already been computed
        return true;
    }
    
    array[0] = get_current_generation();
    array[1] = get_current_replicate();
    array[2] |= AGE;
    // add these age classes
    return false;
}

// ----------------------------------------------------------------------------------------
//////////////////////// FSTAT functions //////////////////////////////////////////////////

// ----------------------------------------------------------------------------------------
// get_allele_freq
// ----------------------------------------------------------------------------------------
/** function to retrieve the local allele frequencies with only a single index
 * i is the patch id, stored are the frequencies however under the sampleID!
 */
template<class SH>
double StatHandler<SH>::get_allele_freq_local(unsigned int i, const age_t & AGE)
{
    age_idx age_pos = age_t2idx(AGE);
    if (!set_alleleFreq(age_pos)) return my_NAN; // pop empty
    
    unsigned int p,l,a;
    fromID(i, p, l, a, _nb_patch, _nb_locus, _nb_allele_max);
    unsigned int sampleID = get_vPatch(p)->get_sampleID();
    if(sampleID==SAMPLED) return my_NAN;
    
    map<unsigned char, double>::iterator pos = _alleleFreq_local[age_pos][sampleID][l].find(a);
    if (pos == _alleleFreq_local[age_pos][sampleID][l].end()) return 0; // allele not present
    return pos->second;
}

// ----------------------------------------------------------------------------------------
// get_allele_freq_global
// ----------------------------------------------------------------------------------------
/** function to retrieve the global allele frequencies with only a single index
 * i is the patch id, stored are the frequencies however under the sampleID!
 */
template<class SH>
double StatHandler<SH>::get_allele_freq_global(unsigned int i, const age_t & AGE)
{
    age_idx age_pos = age_t2idx(AGE);
    if (!set_alleleFreq(age_pos)) return my_NAN; // pop empty
    
    unsigned int l,a;
    fromID(i, l, a, _nb_locus, _nb_allele_max);
    
    map<unsigned char, double>::iterator pos = _alleleFreq_global[age_pos][l].find(a);
    if (pos == _alleleFreq_global[age_pos][l].end()) return 0; // allele not present
    return pos->second;
}

// ----------------------------------------------------------------------------------------
// get_locus_freq_local
// ----------------------------------------------------------------------------------------
/** function to retrieve the local locus genotype frequencies with only a single index
 * i is the patch id, stored are the frequencies however under the sampleID!
 */
template<class SH>
double StatHandler<SH>::get_locus_freq_local(unsigned int i, const age_t & AGE)
{
    age_idx age_pos = age_t2idx(AGE);
    if (!set_locusFreq(age_pos)) return my_NAN; // pop empty
    
    unsigned int p,l,a1,a2;
    fromID(i, p, l, a1, a2, _nb_patch, _nb_locus, _nb_allele_max, _nb_allele_max);
    unsigned int sampleID = get_vPatch(p)->get_sampleID();
    if(sampleID==SAMPLED) return my_NAN;
    
    map<unsigned char, map<unsigned char, double> >::iterator pos = _locusFreq_local[age_pos][sampleID][l].find(a1);
    if (pos == _locusFreq_local[age_pos][sampleID][l].end()) return 0; // allele not present
    
    map<unsigned char, double>::iterator pos2 = pos->second.find(a2);
    if (pos2 == pos->second.end()) return 0; // allele not present
    return pos2->second;
}

// ----------------------------------------------------------------------------------------
// get_locus_freq_global
// ----------------------------------------------------------------------------------------
/** function to retrieve the global locus genotype frequencies with only a single index
 * i is the patch id, stored are the frequencies however under the sampleID!
 */
template<class SH>
double StatHandler<SH>::get_locus_freq_global(unsigned int i, const age_t & AGE)
{
    age_idx age_pos = age_t2idx(AGE);
    if (!set_locusFreq(age_pos)) return my_NAN; // pop empty
    
    unsigned int l, a1, a2;
    fromID(i, l, a1, a2, _nb_locus, _nb_allele_max, _nb_allele_max);
    assert(a1<=a2);
    
    map<unsigned char, map<unsigned char, double> >::iterator pos = _locusFreq_global[age_pos][l].find(a1);
    if (pos == _locusFreq_global[age_pos][l].end()) return 0; // allele not present
    
    map<unsigned char, double>::iterator pos2 = pos->second.find(a2);
    if (pos2 == pos->second.end()) return 0; // allele not present
    return pos2->second;
}

// ----------------------------------------------------------------------------------------
// set_alleleFreqs_local
// ----------------------------------------------------------------------------------------
/** computation of the local and global allele frequencies
 * (have to compute it for all CURRENTLY sampled patches!)
 * _alleleFreq_global[age][locus][allele]
 * _alleleFreq_local[age][sampleID][locus][allele]
 * returns true if the freqs could be computed, i.e. samples were present
 */
template<class SH> bool
StatHandler<SH>::set_alleleFreq(const age_idx & AGE)
{
    // check if the table has already been  created
    if (already_computed(_computed[4], AGE)) return(!_alleleFreq_global[AGE][0].empty());
    
    unsigned int l, p;
    
    if (!_alleleFreq_global) { // first time             <map<unsigned char, double>*>
        ARRAY::create_1D<map<unsigned char, double> *>(_alleleFreq_global, NB_AGE_CLASSES, NULL); // [age][locus][allele]
        ARRAY::create_1D<map<unsigned char, double> **>(_alleleFreq_local, NB_AGE_CLASSES, NULL); // [age][sampleID][locus][allele]
    }
    
    if (!_alleleFreq_global[AGE]) { // first time this age class
        ARRAY::create_1D<map<unsigned char, double> >(_alleleFreq_global[AGE], _nb_locus); // [age][locus][allele]
        ARRAY::create_2D<map<unsigned char, double> >(_alleleFreq_local[AGE],get_current_nbSamplePatch(), _nb_locus); // [age][sampleID][locus][allele]
    }
    else { // reset
        // global freqs
        for (l = 0; l < _nb_locus; ++l) {
            _alleleFreq_global[AGE][l].clear();
        }
        
        // local freqs
        if (get_current_nbSamplePatch() != get_last_nbSamplePatch()) { // number of sampled and populated patches change over time
            ARRAY::delete_2D<map<unsigned char, double> >(_alleleFreq_local[AGE],get_last_nbSamplePatch());
            ARRAY::create_2D<map<unsigned char, double> >(_alleleFreq_local[AGE],get_current_nbSamplePatch(), _nb_locus); // [age][sampleID][locus][allele]
        }
        else { // the number of sampled patches has not changed
            for (p = 0; p < get_current_nbSamplePatch(); ++p) {
                assert(_alleleFreq_local[AGE][p]);
                for (l = 0; l < _nb_locus; ++l) {
                    _alleleFreq_local[AGE][p][l].clear();
                }
            }
        }
    }
    
    
    
    // get the local frequencies and global allele counts for each patch
    unsigned int nbAllele = 0; // total number of alleles
    vector<Patch*>::iterator curPop, endPop;
    for (curPop = get_vSamplePatch().begin(), endPop = get_vSamplePatch().end(); curPop != endPop; ++curPop) {
        assert((*curPop)->get_sampleID()<get_current_nbSamplePatch());
        nbAllele += set_alleleFreq_ofPatch(*curPop, AGE,_alleleFreq_local[AGE][(*curPop)->get_sampleID()], _alleleFreq_global[AGE]);
    }
    
    // compute global allele frequencies
    map<unsigned char, double>::iterator curAllele, endAllele;
    for (l = 0; l < _nb_locus; ++l) {
        curAllele = _alleleFreq_global[AGE][l].begin();
        endAllele = _alleleFreq_global[AGE][l].end();
        for (; curAllele != endAllele; ++curAllele) {
            curAllele->second /= nbAllele; // from counts to frequencies
        }
    }
    
    return(nbAllele != 0);
}

// ----------------------------------------------------------------------------------------
// set_locusFreqs
// ----------------------------------------------------------------------------------------
/** computation of the local and global locus genotype frequencies
 * (have to compute it for all CURRENTLY sampled patches!)
 * _locusFreq_global[age][locus][allele]
 * _locusFreq_local[age][sampleID][locus][allele]
 * returns true if the freqs could be computed, i.e. samples were present
 */
template<class SH> bool
StatHandler<SH>::set_locusFreq(const age_idx & AGE)
{
    // check if the table has already been  created
    if (already_computed(_computed[34], AGE)) return(!_locusFreq_global[AGE][0].empty());
    
    unsigned int l, p;
    
    if (!_locusFreq_global) { // first time             <map<unsigned char, double>*>
        ARRAY::create_1D<map<unsigned char, map<unsigned char, double> >*>(_locusFreq_global, NB_AGE_CLASSES, NULL); // [age][locus][allele]
        ARRAY::create_1D<map<unsigned char, map<unsigned char, double> >**>(_locusFreq_local, NB_AGE_CLASSES, NULL); // [age][sampleID][locus][allele]
    }
    
    if (!_locusFreq_global[AGE]) { // first time this age class
        ARRAY::create_1D<map<unsigned char, map<unsigned char, double> > >(_locusFreq_global[AGE], _nb_locus); // [age][locus][allele]
        ARRAY::create_2D<map<unsigned char, map<unsigned char, double> > >(_locusFreq_local[AGE],	get_current_nbSamplePatch(), _nb_locus); // [age][sampleID][locus][allele]
    }
    else { // reset
        // global freqs
        for (l = 0; l < _nb_locus; ++l) {
            _locusFreq_global[AGE][l].clear();
        }
        
        // local freqs
        if (get_current_nbSamplePatch() != get_last_nbSamplePatch()) { // number of sampled and populated patches change over time
            ARRAY::delete_2D<map<unsigned char, map<unsigned char, double> > >(_locusFreq_local[AGE],get_last_nbSamplePatch());
            ARRAY::create_2D<map<unsigned char, map<unsigned char, double> > >(_locusFreq_local[AGE],get_current_nbSamplePatch(), _nb_locus); // [age][sampleID][locus][allele]
        }
        else { // the number of sampled patches has not changed
            for (p = 0; p < get_current_nbSamplePatch(); ++p) {
                assert(_locusFreq_local[AGE][p]);
                for (l = 0; l < _nb_locus; ++l) {
                    _locusFreq_local[AGE][p][l].clear();
                }
            }
        }
    }
    
    // get the local frequencies and global allele counts for each patch
    unsigned int nbInds = 0; // total number of alleles
    vector<Patch*>::iterator curPop, endPop;
    for (curPop = get_vSamplePatch().begin(), endPop = get_vSamplePatch().end(); curPop != endPop; ++curPop) {
        assert((*curPop)->get_sampleID()<get_current_nbSamplePatch());
        nbInds += set_locusFreq_ofPatch(*curPop, AGE,_locusFreq_local[AGE][(*curPop)->get_sampleID()], _locusFreq_global[AGE]);
    }
    
    // compute global locus genotype frequencies
    map<unsigned char, map<unsigned char, double> >::iterator curAllele1, endAllele1;
    map<unsigned char, double>::iterator curAllele2, endAllele2;
    for (l = 0; l < _nb_locus; ++l) {
        curAllele1 = _locusFreq_global[AGE][l].begin();
        endAllele1 = _locusFreq_global[AGE][l].end();
        for (; curAllele1 != endAllele1; ++curAllele1) {
            curAllele2 = curAllele1->second.begin();
            endAllele2 = curAllele1->second.end();
            for (; curAllele2 != endAllele2; ++curAllele2) {
                curAllele2->second /= nbInds; // from counts to frequencies
            }
        }
    }
    
    return(nbInds != 0);
}

// ----------------------------------------------------------------------------------------
// set_alleleFreq_ofPatch (per patch)
// ----------------------------------------------------------------------------------------
/** computation of the allele frequencies for each patch separately
 * Caution: only the sampled pops are available
 * the number of allels is returned
 */
template<class SH> unsigned int
StatHandler<SH>::set_alleleFreq_ofPatch(Patch * curPatch,
                                        const age_idx & AGE,
                                        map<unsigned char, double>*&freqs,
                                        map<unsigned char, double>*&global_freqs)
{
    unsigned int l, p, nbAllele;
    unsigned char** genes;
    
    // if patch is empty all frequencies are zero
    vector<TIndividual*>&curFem = curPatch->get_sampled_inds(FEM, AGE);
    vector<TIndividual*>&curMal = curPatch->get_sampled_inds(MAL, AGE);
    nbAllele = ((unsigned int)curFem.size() + (unsigned int)curMal.size()) * ploidy;
    if (!nbAllele) return 0; // patch is empty
    
    // get female allele counts
    map<unsigned char, double>::iterator poss;
    vector<TIndividual*>::iterator curInd, endInd;
    
    for (curInd = curFem.begin(), endInd = curFem.end(); curInd != endInd;
         ++curInd) {
        genes = (unsigned char**)(*curInd)->getTrait(_SHLinkedTraitIndex)->get_sequence();
        for (l = 0; l < _nb_locus; ++l) {
            for (p = 0; p < ploidy; ++p) {
                ++freqs[l][genes[l][p]]; // local freqs (it is initialized with zero)
            }
        }
    }
    
    // get male allele counts
    for (curInd = curMal.begin(), endInd = curMal.end(); curInd != endInd;
         ++curInd) {
        genes = (unsigned char**)(*curInd)->getTrait(_SHLinkedTraitIndex)
        ->get_sequence();
        for (l = 0; l < _nb_locus; ++l) {
            for (p = 0; p < ploidy; ++p) {
                ++freqs[l][genes[l][p]]; // local freqs (it is initialized with zero)
            }
        }
    }
    
    // compute local allele frequencies
    map<unsigned char, double>::iterator pos, end;
    map<unsigned char, double>* cur_globalFreqs;
    for (l = 0; l < _nb_locus; ++l) {
        cur_globalFreqs = &(global_freqs[l]);
        for (pos = freqs[l].begin(), end = freqs[l].end(); pos != end; ++pos) {
            (*cur_globalFreqs)[pos->first] += pos->second;
            pos->second /= nbAllele;
        }
    }
    
    return nbAllele;
}

// ----------------------------------------------------------------------------------------
// set_alleleFreq_ofPatch (per patch)
// ----------------------------------------------------------------------------------------
/** computation of the allele frequencies for each patch separately
 * Caution: only the sampled pops are available
 * the number of allels is returned
 */
template<class SH> unsigned int
StatHandler<SH>::set_locusFreq_ofPatch(Patch * curPatch,
                                       const age_idx & AGE,
                                       map<unsigned char, map<unsigned char, double> >*&freqs,
                                       map<unsigned char, map<unsigned char, double> >*&global_freqs)
{
    unsigned int l, nbInds;
    
    // if patch is empty all frequencies are zero
    vector<TIndividual*>&curFem = curPatch->get_sampled_inds(FEM, AGE);
    vector<TIndividual*>&curMal = curPatch->get_sampled_inds(MAL, AGE);
    nbInds = ((unsigned int)curFem.size() + (unsigned int)curMal.size());
    if (!nbInds) return 0; // patch is empty
    
    // get female allele counts
    map<unsigned char, double>::iterator poss;
    vector<TIndividual*>::iterator curInd, endInd;
    unsigned char a1, a2;
    unsigned char** genes;
    
    for (curInd = curFem.begin(), endInd = curFem.end(); curInd != endInd;
         ++curInd) {
        genes = (unsigned char**)(*curInd)->getTrait(_SHLinkedTraitIndex)->get_sequence();
        for (l = 0; l < _nb_locus; ++l) {
            a1 = genes[l][0];
            a2 = genes[l][1];
            if(a1>a2) swap(a1,a2);
            ++freqs[l][a1][a2]; // local freqs (it is initialized with zero)
        }
    }
    
    // get male allele counts
    for (curInd = curMal.begin(), endInd = curMal.end(); curInd != endInd;
         ++curInd) {
        genes = (unsigned char**)(*curInd)->getTrait(_SHLinkedTraitIndex)
        ->get_sequence();
        for (l = 0; l < _nb_locus; ++l) {
            a1 = genes[l][0];
            a2 = genes[l][1];
            if(a1>a2) swap(a1,a2);
            ++freqs[l][a1][a2]; // local freqs (it is initialized with zero)
        }
    }
    
    // compute local locus gentoype frequencies
    map<unsigned char, map<unsigned char, double> >::iterator pos1, end1;
    map<unsigned char, double>::iterator pos2, end2;
    map<unsigned char, map<unsigned char, double> >* cur_globalFreqs;
    map<unsigned char, double>* cur_globalFreqs2;
    for (l = 0; l < _nb_locus; ++l) {
        cur_globalFreqs = &(global_freqs[l]);
        for (pos1 = freqs[l].begin(), end1 = freqs[l].end(); pos1 != end1; ++pos1) {
            cur_globalFreqs2 = &(*cur_globalFreqs)[pos1->first];
            for (pos2 = pos1->second.begin(), end2 = pos1->second.end(); pos2 != end2; ++pos2) {
                global_freqs[l][pos1->first][pos2->first] += pos2->second;
                (*cur_globalFreqs2)[pos2->first] += pos2->second;           // sum up for gloabl counts
                pos2->second /= nbInds;                                     // counts to frwqs
            }
        }
    }
    
    return nbInds;
}

// ----------------------------------------------------------------------------------------
// set_alleleFreq_ofPatch_allInds (per patch)
// ----------------------------------------------------------------------------------------
/** computation of the allele frequencies for each patch separately
 * Caution: ALL indivuals of the patch are considered
 * the number of allels is returned
 */
template<class SH> unsigned int
StatHandler<SH>::set_alleleFreq_ofPatch_allInds(Patch * curPatch,
                                                const age_idx & AGE,
                                                map<unsigned char, double>*&freqs,
                                                map<unsigned char, double>*&global_freqs)
{
    unsigned int l, p, nbAllele;
    unsigned char** genes;
    
    // if patch is empty all frequencies are zero
    vector<TIndividual*>&curFem = curPatch->get_all_inds(FEM, AGE);
    vector<TIndividual*>&curMal = curPatch->get_all_inds(MAL, AGE);
    nbAllele = ((unsigned int)curFem.size() + (unsigned int)curMal.size()) * ploidy;
    if (!nbAllele) return 0; // patch is empty
    
    // get female allele counts
    map<unsigned char, double>::iterator poss;
    vector<TIndividual*>::iterator curInd, endInd;
    
    for (curInd = curFem.begin(), endInd = curFem.end(); curInd != endInd;
         ++curInd) {
        genes = (unsigned char**)(*curInd)->getTrait(_SHLinkedTraitIndex)->get_sequence();
        for (l = 0; l < _nb_locus; ++l) {
            for (p = 0; p < ploidy; ++p) {
                ++freqs[l][genes[l][p]]; // local freqs (it is initialized with zero)
            }
        }
    }
    
    // get male allele counts
    for (curInd = curMal.begin(), endInd = curMal.end(); curInd != endInd;
         ++curInd) {
        genes = (unsigned char**)(*curInd)->getTrait(_SHLinkedTraitIndex)
        ->get_sequence();
        for (l = 0; l < _nb_locus; ++l) {
            for (p = 0; p < ploidy; ++p) {
                ++freqs[l][genes[l][p]]; // local freqs (it is initialized with zero)
            }
        }
    }
    
    // compute local allele frequencies
    map<unsigned char, double>::iterator pos, end;
    map<unsigned char, double>* cur_globalFreqs;
    for (l = 0; l < _nb_locus; ++l) {
        cur_globalFreqs = &(global_freqs[l]);
        for (pos = freqs[l].begin(), end = freqs[l].end(); pos != end; ++pos) {
            (*cur_globalFreqs)[pos->first] += pos->second;
            pos->second /= nbAllele;
        }
    }
    
    return nbAllele;
}

// ----------------------------------------------------------------------------------------
// get_locusGenotypeCounts_ofPatch_andSex (per patch & sex)
// ----------------------------------------------------------------------------------------
/** computation of the locus genotype count for each patch and sex separately
 * the number of indivduals is returned
 * Note, freqs is not checked if it is empty!
 */
template<class SH> unsigned int
StatHandler<SH>::get_locusGenotypeCounts_ofPatch_andSex(Patch * curPatch,
                                                        const age_idx & AGE, sex_t SEX,
                                                        unsigned int traitID,
                                                        map<unsigned char, map<unsigned char, double> >*& freqs)
{
    unsigned int l;
    unsigned char** genes;
    unsigned char a1, a2;
    assert(freqs);
    
    // if patch is empty all frequencies are zero
    vector<TIndividual*>&curFem = curPatch->get_sampled_inds(SEX, AGE);
    if (curFem.empty()) return 0; // patch is empty
    
    // get locus genotype counts
    vector<TIndividual*>::iterator curInd, endInd;
    for (curInd = curFem.begin(), endInd = curFem.end(); curInd != endInd; ++curInd) {
        genes = (unsigned char**)(*curInd)->getTrait(traitID)->get_sequence();
        for (l = 0; l < _nb_locus; ++l) {
            a1 = genes[l][0];
            a2 = genes[l][1];
            if(a1<a2) ++freqs[l][a1][a2];
            else      ++freqs[l][a2][a1];
        }
    }
    
    return (unsigned int)curFem.size();
}

// ----------------------------------------------------------------------------------------
// get_locusGenotypeCounts_ofPatch_andSex (per patch)
// ----------------------------------------------------------------------------------------
/** computation of the locus genotype counts for each patch separately
 * the number of indivduals is returned
 * Note, freqs is not checked if it is empty!
 */
template<class SH> unsigned int
StatHandler<SH>::get_locusGenotypeCounts_ofPatch(Patch * curPatch, const age_idx & AGE,
                                                 unsigned int traitID,
                                                 map<unsigned char, map<unsigned char, double> >*& freqs)
{
    assert(freqs);
    
    return (get_locusGenotypeCounts_ofPatch_andSex(curPatch, AGE, FEM, traitID, freqs)
            + get_locusGenotypeCounts_ofPatch_andSex(curPatch, AGE, MAL, traitID, freqs));
}

// ----------------------------------------------------------------------------------------
// get_locusGenotypeFreqs_ofPatch (per patch)
// ----------------------------------------------------------------------------------------
/** computation of the locus genotype frequencies for each patch separately
 * the number of indivduals is returned
 * here the _SHLinkedTraitIndex is used!
 */
template<class SH> unsigned int
StatHandler<SH>::get_locusGenotypeFreqs_ofPatch(Patch * curPatch, const age_idx & AGE,
                                                map<unsigned char, map<unsigned char, double> >*& freqs)
{
    get_locusGenotypeFreqs_ofPatch(curPatch, AGE, _SHLinkedTraitIndex, freqs);
}

// ----------------------------------------------------------------------------------------
// get_locusGenotypeFreqs_ofPatch (per patch)
// ----------------------------------------------------------------------------------------
/** computation of the locus genotype frequencies for each patch separately
 * the number of indivduals is returned
 */
template<class SH> unsigned int
StatHandler<SH>::get_locusGenotypeFreqs_ofPatch(Patch * curPatch, const age_idx & AGE,
                                                unsigned int traitID,
                                                map<unsigned char, map<unsigned char, double> >*& freqs)
{
    unsigned int l;
    if(!freqs) freqs = ARRAY::new_1D<map<unsigned char, map<unsigned char, double> > >(_nb_locus);
    else{
        for (l = 0; l < _nb_locus; ++l) {
            freqs[l].clear();
        }
    }
    
    unsigned int nbInd = get_locusGenotypeCounts_ofPatch(curPatch, AGE, traitID, freqs);
    
    // compute local frequencies
    map<unsigned char, map<unsigned char, double> >::iterator pos1, end1;
    map<unsigned char, double>::iterator pos2, end2;
    for (l = 0; l < _nb_locus; ++l) {
        for (pos1 = freqs[l].begin(), end1 = freqs[l].end(); pos1 != end1; ++pos1) {
            for (pos2 = pos1->second.begin(), end2 = pos1->second.end(); pos2 != end2; ++pos2) {
                pos2->second /= nbInd;
            }
        }
    }
    return nbInd;
}

// ----------------------------------------------------------------------------------------
// get_locusGenotypeFreqs_ofPatch_allInds (per patch)
// ----------------------------------------------------------------------------------------
/** computation of the locus genotype frequencies for each patch separately
 * Caution: all individuals  are considered
 * the number of indivduals is returned
 */
template<class SH> unsigned int
StatHandler<SH>::get_locusGenotypeFreqs_ofPatch_allInds(Patch * curPatch, const age_idx & AGE,
                                                        unsigned int traitID,
                                                        map<unsigned char, map<unsigned char, double> >*& freqs)
{
    unsigned int l;
    if(!freqs) freqs = ARRAY::new_1D<map<unsigned char, map<unsigned char, double> > >(_nb_locus);
    else{
        for (l = 0; l < _nb_locus; ++l) {
            freqs[l].clear();
        }
    }
    
    // if patch is empty all frequencies are zero
    vector<TIndividual*>&curFem = curPatch->get_all_inds(FEM, AGE);
    vector<TIndividual*>&curMal = curPatch->get_all_inds(MAL, AGE);
    unsigned int nbInd = (unsigned int) (curFem.size() + curMal.size());
    if (!nbInd) return 0; // patch is empty
    
    // get locus genotype counts FEM
    vector<TIndividual*>::iterator curInd, endInd;
    unsigned char** genes;
    unsigned char a1, a2;
    for (curInd = curFem.begin(), endInd = curFem.end(); curInd != endInd; ++curInd) {
        genes = (unsigned char**)(*curInd)->getTrait(traitID)->get_sequence();
        for (l = 0; l < _nb_locus; ++l) {
            a1 = genes[l][0];
            a2 = genes[l][1];
            if(a1<a2) ++freqs[l][a1][a2];
            else      ++freqs[l][a2][a1];
        }
    }
    
    // get locus genotype counts MAL
    for (curInd = curMal.begin(), endInd = curMal.end(); curInd != endInd; ++curInd) {
        genes = (unsigned char**)(*curInd)->getTrait(traitID)->get_sequence();
        for (l = 0; l < _nb_locus; ++l) {
            a1 = genes[l][0];
            a2 = genes[l][1];
            if(a1<a2) ++freqs[l][a1][a2];
            else      ++freqs[l][a2][a1];
        }
    }
    
    // compute local genotype frequencies
    map<unsigned char, map<unsigned char, double> >::iterator pos1, end1;
    map<unsigned char, double>::iterator pos2, end2;
    for (l = 0; l < _nb_locus; ++l) {
        for (pos1 = freqs[l].begin(), end1 = freqs[l].end(); pos1 != end1; ++pos1) {
            for (pos2 = pos1->second.begin(), end2 = pos1->second.end(); pos2 != end2; ++pos2) {
                pos2->second /= nbInd;
            }
        }
    }
    return nbInd;
}

// ----------------------------------------------------------------------------------------
// set_alleleFreq_global (globally)
// ----------------------------------------------------------------------------------------
/** computation of the global allele frequencies of the given patches using the local allele frerqs
 */
template<class SH> void
StatHandler<SH>::set_alleleFreq_global(map<unsigned char, double> **localFreqs,
                                       map<unsigned char, double> *globalFreqs,
                                       unsigned int*popSizes, const unsigned int&tot_popSize,
                                       const unsigned int&nbPatch)
{
    // get the global allele frequencies for each locus
    map<unsigned char, double>::iterator pos, end;
    for (unsigned int l = 0; l < _nb_locus; ++l) { // for each locus
        set_alleleFreq_global_ofLocus(l, localFreqs, &globalFreqs[l], popSizes, tot_popSize, nbPatch);
    }
}

// ----------------------------------------------------------------------------------------
// set_alleleFreq_global (globally)
// ----------------------------------------------------------------------------------------
/** computation of the global allele frequencies for the locus and given patches
 * popSizes has to be the sample size!!!
 */
template<class SH> void
StatHandler<SH>::set_alleleFreq_global_ofLocus(const unsigned int&l, map<unsigned char, double> **localFreqs,
                                               map<unsigned char, double> *globalFreqs, unsigned int*popSizes,
                                               const unsigned int&tot_popSize, const unsigned int&nbPatch)
{
    assert(tot_popSize == ARRAY::sum(popSizes, nbPatch));
    
    globalFreqs->clear(); // clear the global freqs of this locus
    
    // get the global allele count
    map<unsigned char, double>::iterator pos, end;
    for (unsigned int p = 0; p < nbPatch; ++p) { // for each patch
        for (pos = localFreqs[p][l].begin(), end = localFreqs[p][l].end();
             pos != end; ++pos) { // for each present allele
            (*globalFreqs)[pos->first] += (pos->second * popSizes[p]);
            // initialized with zero
        }
    }
    
    // get the global allele frequency
    pos = globalFreqs->begin();
    end = globalFreqs->end();
    for (; pos != end; ++pos) { // for each present allele
        pos->second /= tot_popSize; // devide by the total size
    } // "the factor 2 of diploidy can be ignored..."
}

// ----------------------------------------------------------------------------------------
// get_genotypeFreq
// ----------------------------------------------------------------------------------------
/** get the genotype frequencies of the genotypes at the given patch and loci
 * the returned map has to be deleted!
 * returned map: [allele1][allele2]
 */
template<class SH>map<unsigned char, map<unsigned char,double> > *
StatHandler<SH>::get_genotypeFreq(const age_idx & AGE, Patch* curPop, const unsigned int&l1, const unsigned int&l2)
{
    map<unsigned char, map<unsigned char, double> > *freqs = new map<unsigned char, map<unsigned char, double> >;
    unsigned int a;
    unsigned char ** genes;
    vector<TIndividual*>::iterator curInd, endInd;
    
    // females
    vector<TIndividual*>& curFem = curPop->get_sampled_inds(FEM, AGE);
    for (curInd = curFem.begin(), endInd = curFem.end(); curInd != endInd; ++curInd) {
        genes = (unsigned char**)(*curInd)->getTrait(_SHLinkedTraitIndex)->get_sequence();
        for (a = 0; a < ploidy; ++a) {
            ++(*freqs)[genes[l1][a]][genes[l2][a]]; // it is initialized with zero
        }
    }
    
    // males
    vector<TIndividual*>&curMal = curPop->get_sampled_inds(MAL, AGE);
    for (curInd = curMal.begin(), endInd = curMal.end(); curInd != endInd; ++curInd) {
        genes = (unsigned char**)(*curInd)->getTrait(_SHLinkedTraitIndex)->get_sequence();
        for (a = 0; a < ploidy; ++a) {
            ++(*freqs)[genes[l1][a]][genes[l2][a]]; // it is initialized with zero
        }
    }
    
    // compute genotype frequencies
    unsigned int totSize = 2*((unsigned int)curFem.size()+(unsigned int)curMal.size()); // from inds to number of alleles
    map<unsigned char, map<unsigned char, double> >::iterator pos1, end1;
    map<unsigned char, double>::iterator pos2, end2;
    for (pos1 = freqs->begin(), end1 = freqs->end(); pos1 != end1; ++pos1) {
        for (pos2 = pos1->second.begin(), end2 = pos1->second.end(); pos2 != end2; ++pos2) {
            pos2->second /= totSize;
        }
    }
    
    return freqs;
}

// ----------------------------------------------------------------------------------------
// setFstat
// ----------------------------------------------------------------------------------------
/** compute global F-satistics averaged across loci and pops following Nei & Chesser (1983)
 * Hs  = Mean_patch(Mean_loci(1-SUM_allele(p^2))  // p: allele freq of locus l and patch p NaN if metapop empty
 * Ht  = Mean_loci(1-SUM_allele(mean(p)^2))       // p: allele freq of locus l NaN if metapop empty
 * Fst = 1-Hs/Ht NaN if metapop empty, 0 if Ht=0 (all individuals are monomorph for a SINGLE allele)
 * Fis = Ho/Hs NaN if metapop empty, 0 if Hs=0 (all individual of a patch are monomorph for an allele (may differ among patches))
 * Fit = Ho/Ht NaN if metapop empty, 0 if Ht=0 (all individuals are monomorph for a SINGLE allele)
 */
template<class SH>void StatHandler<SH>::setFstat_Nei_Chesser(const age_idx & AGE)
{
    // check if the table has already been computed
    if (already_computed(_computed[6], AGE)) return;
    
    vector<Patch*>&aPatch = get_vSamplePatch();
    unsigned int nbPops;
    double H = getHarmonicMean_ofPopSize(AGE, aPatch, nbPops); // nbPops (full) is set by the function
    
    // if metapop is empty
    if (!nbPops) {
        _hsnei[AGE] = _htnei[AGE] = _fis[AGE] = _fit[AGE] = _fst[AGE] = my_NAN;
        return;
    }
    
    getHo(AGE); // mean observed heterozygostiy
    getHs(AGE); // mean expected heterozyogsity
    
    // Fis
    _hsnei[AGE] = H != 1 ? H / (H - 1.0) * (_hs[AGE] - (_ho[AGE] / (2.0 * H))) : _hs[AGE]; // Nei's corrections:
    _fis[AGE] = _hsnei[AGE] ? 1.0 - (_ho[AGE] / _hsnei[AGE]) : 0; // monomorphic = ht=0: after definition 0 and not NaN!
    
    if (nbPops > 1) { // more than one colonized patch
        getHt(AGE);
        _htnei[AGE] = _ht[AGE] + (_hsnei[AGE] / (H * nbPops)) - (_ho[AGE] / (2.0 * H * nbPops)); // Nei's corrections:
        _fit[AGE] = _htnei[AGE] ? 1.0 - (_ho[AGE] / _htnei[AGE]) : 0;     // monomorphic = ht=0: after definition 0 and not NaN!
        _fst[AGE] = _htnei[AGE] ? 1.0 - (_hsnei[AGE] / _htnei[AGE]) : 0;  // monomorphic = ht=0: after definition 0 and not NaN!
    }
    else _ht[AGE] = _htnei[AGE] = _fit[AGE] = _fst[AGE] = my_NAN;
}

// ----------------------------------------------------------------------------------------
// setFstatMatrix
// ----------------------------------------------------------------------------------------
/** computes the pairwise Fst following Nei & Chesser (1983) and stores them in a matrix.
 * All pairs are computed!
 * Values are stored in a matrix.
 */
template<class SH>void StatHandler<SH>::setFstat_Nei_Chesser_perPatchPair(const age_idx & AGE)
{
    // check if the table has already been computed
    if (already_computed(_computed[7], AGE)) return;
    
    double H, hsnei, htnei, ho, ht, hs;
    unsigned int i, j;
    
    if (!_fst_matrix) ARRAY::create_1D<double**>(_fst_matrix, NB_AGE_CLASSES, NULL); // first time
    if (!_fst_matrix[AGE]) 	ARRAY::create_2D(_fst_matrix[AGE], get_current_nbSamplePatch(), get_current_nbSamplePatch(), (double)my_NAN);
    else if (get_current_nbSamplePatch() != get_last_nbSamplePatch()) { // number of sampeld patches may change over time
        ARRAY::delete_2D(_fst_matrix[AGE], get_last_nbSamplePatch());
        ARRAY::create_2D(_fst_matrix[AGE], get_current_nbSamplePatch(), get_current_nbSamplePatch(), (double)my_NAN);
    }
    else ARRAY::reset_2D(_fst_matrix[AGE], get_current_nbSamplePatch(), get_current_nbSamplePatch(), (double)my_NAN);
    
    double * hs_pop = new double[get_current_nbSamplePatch()];
    double * ho_pop = new double[get_current_nbSamplePatch()];
    vector<Patch*>::iterator curPop1, curPop2, endPop = get_vSamplePatch().end();
    for (i = 0, curPop1 = get_vSamplePatch().begin(); curPop1 != endPop; ++i, ++curPop1) { // for each patch
        ho_pop[i] = getHo_ofPatch(*curPop1, AGE);
        hs_pop[i] = getHs_ofPatch(*curPop1, AGE);
        
    }
    
    // for each pair of pops
    // we have to check all pops, since the Fst after Nei is also defined for a
    // single populated pop, but not if both pops are empty
    for (i = 0, curPop1 = get_vSamplePatch().begin(); curPop1 != endPop;
         ++i, ++curPop1) { // for each patch
        if (hs_pop[i] == my_NAN) continue; // if the patch is not populated of the given AGE
        
        curPop2 = curPop1;
        for (j = i + 1, ++curPop2; curPop2 != endPop; ++j, ++curPop2) {
            // for each other patch
            assert(j == (*curPop2)->get_sampleID());
            if (hs_pop[j] == my_NAN) continue; // if the patch is not populated
            
            H = 2.0 / ((1.0 / (*curPop1)->sampleSize(AGE)) + (1.0 / (*curPop2)->sampleSize(AGE))); // harmonic mean of N
            hs = (hs_pop[i] + hs_pop[j]) / 2.0;
            ho = (ho_pop[i] + ho_pop[j]) / 2.0;
            ht = getHt(AGE, *curPop1, *curPop2);
            
            // Nei's corrections:
            hsnei = H != 1 ? H / (H - 1.0) * (hs - (ho / (2.0 * H))) : hs;
            htnei = ht + (hsnei / (H * 2)) - (ho / (2.0 * H * 2));
            _fst_matrix[AGE][i][j] = htnei ? 1.0 - (hsnei / htnei) : 0; // monomorphic = ht=0: after definition 0 and not NaN!
        }
    }
    
    delete[]hs_pop;
    delete[]ho_pop;
}

// ----------------------------------------------------------------------------------------
// setFstatMatrix
// ----------------------------------------------------------------------------------------
/** computes the Fst between the specified patches following Nei & Chesser (1983)
 * p1 and p2 are the patchID's
 * nothing is stored in a matrix
 */
template<class SH> double
StatHandler<SH>::getFstat_Nei_Chesser_perPatchPair(const age_idx & AGE, Patch * curPatch1, Patch * curPatch2)
{
    // get population sizes
    unsigned int N1 = curPatch1->sampleSize(AGE);
    if (!N1) return my_NAN; // if not populated
    
    unsigned int N2 = curPatch2->sampleSize(AGE);
    if (!N2) return my_NAN; // if not populated
    
    double ho = (getHo_ofPatch(curPatch1, AGE) + getHo_ofPatch(curPatch2, AGE)) / 2;
    double hs = (getHs_ofPatch(curPatch1, AGE) + getHs_ofPatch(curPatch2, AGE)) / 2;
    double ht = getHt(AGE, curPatch1, curPatch2);
    double H = 2.0 / ((1.0 / N1) + (1.0 / N2)); // harmonic mean
    
    // Nei's corrections:
    double hsnei = H != 1 ? H / (H - 1.0) * (hs - (ho / (2.0 * H))) : hs;
    double htnei = ht + (hsnei / (H * 2)) - (ho / (2.0 * H * 2));
    if (!htnei) return 0; // monomorphic = ht=0: after definition 0 and not NaN!
    return 1.0 - (hsnei / htnei);
}

// ----------------------------------------------------------------------------------------
// setFstatMatrix
// ----------------------------------------------------------------------------------------
/** computes the Fst between the specified patches following Nei & Chesser (1983)
 * p1 and p2 are the patchID's
 * nothing is stored in a matrix
 */
template<class SH> double
StatHandler<SH>::getFstat_Nei_Chesser_perPatchPair_andLocus(const age_idx & AGE,
                                                            Patch * curPatch1, Patch * curPatch2, unsigned int locus)
{
    // get population sizes
    unsigned int N1 = curPatch1->sampleSize(AGE);
    if (!N1) return my_NAN; // if not populated
    
    unsigned int N2 = curPatch2->sampleSize(AGE);
    if (!N2) return my_NAN; // if not populated
    
    double ho = (getHo_ofPatch_andLocus(AGE, curPatch1, locus) + getHo_ofPatch_andLocus(AGE, curPatch2, locus)) / 2;
    double hs = (getHs_ofPatch_andLocus(AGE, curPatch1, locus) + getHs_ofPatch_andLocus(AGE, curPatch2, locus)) / 2;
    double ht = getHt_ofLocus(AGE, locus, curPatch1, curPatch2);
    double H = 2.0 / ((1.0 / N1) + (1.0 / N2)); // harmonic mean
    
    // Nei's corrections:
    double hsnei = H != 1 ? H / (H - 1.0) * (hs - (ho / (2.0 * H))) : hs;
    double htnei = ht + (hsnei / (H * 2)) - (ho / (2.0 * H * 2));
    if (!htnei) return 0; // monomorphic = ht=0: after definition 0 and not NaN!
    return 1.0 - (hsnei / htnei);
}

// ----------------------------------------------------------------------------------------
// setFstat
// ----------------------------------------------------------------------------------------
// following Nei & Chesser (1983)
template<class SH>void StatHandler<SH>::setFstat_Nei_Chesser_perLocus
(const age_idx & AGE) {
    // check if the table has already been computed
    if (already_computed(_computed[8], AGE))
        return;
    
    if (!_hsnei_locus) { // first time
        ARRAY::create_1D<double*>(_hsnei_locus, NB_AGE_CLASSES, NULL);
        ARRAY::create_1D<double*>(_htnei_locus, NB_AGE_CLASSES, NULL);
        ARRAY::create_1D<double*>(_fis_locus, NB_AGE_CLASSES, NULL);
        ARRAY::create_1D<double*>(_fit_locus, NB_AGE_CLASSES, NULL);
        ARRAY::create_1D<double*>(_fst_locus, NB_AGE_CLASSES, NULL);
        ARRAY::create_1D<double*>(_hs_locus, NB_AGE_CLASSES, NULL);
        ARRAY::create_1D<double*>(_ht_locus, NB_AGE_CLASSES, NULL);
    }
    
    if (!_hsnei_locus[AGE]) {
        _hsnei_locus[AGE] = new double[_nb_locus];
        _htnei_locus[AGE] = new double[_nb_locus];
        _fis_locus[AGE] = new double[_nb_locus];
        _fit_locus[AGE] = new double[_nb_locus];
        _fst_locus[AGE] = new double[_nb_locus];
        _hs_locus[AGE] = new double[_nb_locus];
        _ht_locus[AGE] = new double[_nb_locus];
    }
    
    // compute harmonic mean of patch sizes:
    unsigned int i, nbPatch;
    
    double H = getHarmonicMean_ofPopSize(AGE, get_vSamplePatch(), nbPatch);
    // nbPatch is set by the function
    
    // if the patches are empty
    if (!nbPatch) {
        for (i = 0; i < _nb_locus; ++i) {
            _hsnei_locus[AGE][i] = _htnei_locus[AGE][i] = _fis_locus[AGE][i] = _fit_locus[AGE][i] = _fst_locus[AGE][i] = my_NAN;
        }
        return;
    }
    
    getHo_perLocus(AGE); // _ho: computed for all loci at a time
    
    for (i = 0; i < _nb_locus; ++i) { // for each locus
        _hs_locus[AGE][i] = getHs_ofLocus(i, AGE); // compute Hs
        _hsnei_locus[AGE][i] = H != 1 ? H / (H - 1.0) * (_hs_locus[AGE][i] - (_ho_locus[AGE][i] / (2.0 * H))) : _hs_locus[AGE][i]; // Nei's corrections:
        _fis_locus[AGE][i] = _hsnei_locus[AGE][i] ? 1.0 - (_ho_locus[AGE][i] / _hsnei_locus[AGE][i]) : 0; // monomorphic = ht=0: after definition 0 and not NaN!
        
        if (nbPatch > 1) { // more than one colonized patch
            _ht_locus[AGE][i] = getHt_ofLocus(AGE, i); // compute Ht
            _htnei_locus[AGE][i] = _ht_locus[AGE][i] + (_hsnei_locus[AGE][i] / (H * nbPatch)) - (_ho_locus[AGE][i] / (2.0 * H * nbPatch)); // Nei's corrections:
            _fit_locus[AGE][i] = _htnei_locus[AGE][i] ? 1.0 - (_ho_locus[AGE][i] / _htnei_locus[AGE][i]) : 0;    // monomorphic = ht=0: after definition 0 and not NaN!
            _fst_locus[AGE][i] = _htnei_locus[AGE][i] ? 1.0 - (_hsnei_locus[AGE][i] / _htnei_locus[AGE][i]) : 0; // monomorphic = ht=0: after definition 0 and not NaN!
        }
        else _ht_locus[AGE][i] = _htnei_locus[AGE][i] = _fit_locus[AGE][i] = _fst_locus[AGE][i] = my_NAN;
    }
}

// ----------------------------------------------------------------------------------------
// setFstat
// ----------------------------------------------------------------------------------------
// following fdist (Beaumont & Nichols (1996)
template<class SH>void StatHandler<SH>::setFstat_Beaumont_Nichols_perLocus
(const age_idx & AGE) {
    // check if the table has already been computed
    if (already_computed(_computed[31], AGE))
        return;
    
    if (!_het0_locus) { // first time
        ARRAY::create_1D<double*>(_het0_locus, NB_AGE_CLASSES, NULL);
        ARRAY::create_1D<double*>(_het1_locus, NB_AGE_CLASSES, NULL);
        ARRAY::create_1D<double*>(_fst_bn_locus, NB_AGE_CLASSES, NULL);
    }
    if (!_het0_locus[AGE]) {
        _het0_locus[AGE] = new double[_nb_locus];
        _het1_locus[AGE] = new double[_nb_locus];
        _fst_bn_locus[AGE] = new double[_nb_locus];
    }
    
    unsigned int i;
    
    if (!set_alleleFreq(AGE)) { // pop empty
        for (i = 0; i < _nb_locus; ++i) {
            _het0_locus[AGE][i] = _het1_locus[AGE][i] = _fst_bn_locus[AGE][i] = my_NAN;
        }
        return;
    }
    
    unsigned int l, size1, nbPatch = 0;
    double x0, x2, yy;
    map<unsigned char, double>::iterator pos1, end1, pos2, end2;
    map<unsigned char, double> *curAllele1, *curAllele2;
    vector<Patch*>::iterator curPop1, curPop2, endPop = get_vSamplePatch().end();
    
    for (l = 0; l < _nb_locus; ++l) { // for each locus
        x0 = 0.0;
        // for each pop
        for (curPop1 = get_vSamplePatch().begin(); curPop1 != endPop; ++curPop1) {
            // for each pop
            size1 = 2 * (*curPop1)->sampleSize(AGE);
            if (!size1) continue; // if the patch is empty for these age classes
            ++nbPatch;
            x2 = 0.0;
            curAllele1 = _alleleFreq_local[(*curPop1)->get_sampleID()][l];
            for (pos1 = curAllele1->begin(), end1 = curAllele1->end(); pos1 != end1;
                 ++pos1) { // for each allele
                x2 += pos1->second * pos1->second;
            }
            x0 += (size1 * x2 - 1) / (size1 - 1);
        }
        
        // for each pair of pops
        yy = 0.0;
        
        for (curPop1 = get_vSamplePatch().begin(); curPop1 != endPop; ++curPop1) {
            // for each pop
            size1 = 2 * (*curPop1)->sampleSize(AGE);
            if (!size1) continue; // if the patch is empty for these age classes
            curAllele1 = _alleleFreq_local[(*curPop1)->get_sampleID()][l];
            curPop2 = curPop1;
            for (++curPop2; curPop2 != endPop; ++curPop2) { // for each pair of pops
                curAllele2 = _alleleFreq_local[(*curPop2)->get_sampleID()][l];
                for (pos1 = curAllele1->begin(), end1 = curAllele1->end(),
                     // for each allele in pop1...
                     end2 = curAllele2->end(); pos1 != end1; ++pos1, ++pos2) {
                    pos2 = curAllele2->find(pos1->first);
                    // ...check if it is present in pop2
                    if (pos2 == end2) continue;
                    yy += pos1->second * pos2->second;
                }
            }
        }
        
        _het0_locus[AGE][l] = 1.0 - x0 / nbPatch;
        _het1_locus[AGE][l] = 1.0 - 2 * yy / (nbPatch * (nbPatch - 1));
        if (_het1_locus[AGE][l] < 1.0e-10)
            _fst_bn_locus[AGE][l] = my_NAN;
        else
            _fst_bn_locus[AGE][l] = 1.0 - _het0_locus[AGE][l] / _het1_locus[AGE][l];
    }
}

// ----------------------------------------------------------------------------------------
// getNbAllele_ofLocus
// ----------------------------------------------------------------------------------------
/** get the mean number of alleles across patches of the given locus
 */
template<class SH>
double StatHandler<SH>::getNbAllele_ofLocus(unsigned int l,	const age_t & AGE)
{
    // get the local allele frequencies if necessary
    age_idx age_pos = age_t2idx(AGE);
    double curNb, sum = 0;
    unsigned int nbFullPops = 0;
    
    vector<Patch*>::iterator curPop, endPop = get_vSamplePatch().end();
    for (curPop = get_vSamplePatch().begin(); curPop != endPop; ++curPop) {
        // for each patch
        curNb = getNbAllele_ofPatch_andLocus(age_pos, *curPop, l);
        if (curNb == my_NAN) 	continue; // if patch is empty
        sum += curNb;
        ++nbFullPops;
    }
    return nbFullPops ? sum / nbFullPops : my_NAN;
}

// ----------------------------------------------------------------------------------------
// getNbAllele_ofPatch
// ----------------------------------------------------------------------------------------
/** get the mean number of alleles per locus of this patch specified by the patchID
 */
template<class SH>
double StatHandler<SH>::getNbAllele_ofPatch(unsigned int p,	const age_t & AGE)
{
    return getNbAllele_ofPatch(get_vPatch(p), age_t2idx(AGE));
}

// ----------------------------------------------------------------------------------------
// getNbAllele_ofPatch
// ----------------------------------------------------------------------------------------
/** get the mean number of alleles per locus of the given patch
 * !! the allele frequencies have to be set before (not anymore tested) !!
 */
template<class SH>
double StatHandler<SH>::getNbAllele_ofPatch(Patch * p,	const age_idx & AGE)
{
    if (p->get_sampleID() == SAMPLED) return my_NAN;
    
    double sum = 0;
    for (unsigned int l = 0; l < _nb_locus; ++l) {
        sum += getNbAllele_ofPatch_andLocus(AGE, p,	l);
    }
    return sum / _nb_locus;
}

// ----------------------------------------------------------------------------------------
// getNbAllele_ofPatch_andLocus
// ----------------------------------------------------------------------------------------
/** get the number of alleles of the given locus and patch */
template<class SH>
unsigned int StatHandler<SH>::getNbAllele_ofPatch_andLocus(const age_idx & AGE, Patch * p,	const unsigned int& l)
{
    if (p->get_sampleID() == SAMPLED) return my_NAN;
    set_alleleFreq(AGE);
    return (unsigned int)_alleleFreq_local[AGE][p->get_sampleID()][l].size();
}

// ----------------------------------------------------------------------------------------
// getNbAllele
// ----------------------------------------------------------------------------------------
/** get the mean number of alleles per locus averaged across all sampled populated patches
 */
template<class SH>double StatHandler<SH>::getNbAllele(const age_t & AGE) {
    // get the local allele frequencies if necessary
    age_idx age_pos = age_t2idx(AGE);
    if (!set_alleleFreq(age_pos)) return my_NAN;
    
    double curNb, sum = 0;
    unsigned int nbFullPops = 0;
    
    vector<Patch*>::iterator curPop, endPop = get_vSamplePatch().end();
    for (curPop = get_vSamplePatch().begin(); curPop != endPop; ++curPop) {
        // for each patch
        curNb = getNbAllele_ofPatch(*curPop, age_pos);
        if (curNb == my_NAN) 	continue; // if patch is empty
        sum += curNb;
        ++nbFullPops;
    }
    return sum / nbFullPops;
}

// ----------------------------------------------------------------------------------------
// getNbAlleleTot
// ----------------------------------------------------------------------------------------
/** return mean number of alleles per locus of the entire metapopulation
 */
template<class SH>
double StatHandler<SH>::getNbAlleleTot(const age_idx & AGE) {
    if (!set_alleleFreq(AGE)) return my_NAN;
    unsigned int l, sum = 0;
    for (l = 0; l < _nb_locus; ++l) {
        sum += _alleleFreq_global[AGE][l].size();
    }
    return(double)sum / _nb_locus;
}

// ----------------------------------------------------------------------------------------
// getNbAlleleTot
// ----------------------------------------------------------------------------------------
/** return total number of alleles of locus l */
template<class SH>
unsigned int StatHandler<SH>::getNbAlleleTot_ofLocus(const age_idx & AGE, const unsigned int& l) {
    if (!set_alleleFreq(AGE)) return my_NAN;
    return (unsigned int)_alleleFreq_global[AGE][l].size();
}

// ----------------------------------------------------------------------------------------
// getNbFixedLocus_ofPatch
// ----------------------------------------------------------------------------------------
/** get the number of fixed locus of this patch specified by the patchID
 * doubel since it is a stat
 */
template<class SH>
double
StatHandler<SH>::getNbFixedLocus_ofPatch(unsigned int p, const age_t & AGE) {
    return (double)getNbFixedLocus_ofPatch(get_vPatch(p), age_t2idx(AGE));
}

// ----------------------------------------------------------------------------------------
// getNbFixedLocus_ofPatch
// ----------------------------------------------------------------------------------------
/** get the number of fixed locus of this patch specified by the patchID
 * !! the allele frequencies have to be set before (not anymore tested) !!
 */
template<class SH>
unsigned int StatHandler<SH>::getNbFixedLocus_ofPatch(Patch * p, const age_idx & AGE) {
    
    assert(p->get_sampleID()!=my_NAN);
    if(p->get_sampleID()==SAMPLED) return my_NAN;
    set_alleleFreq(AGE);
    
    map<unsigned char, double> *curPatchFreqs = _alleleFreq_local[AGE][p->get_sampleID()];
    if (curPatchFreqs[0].empty()) return my_NAN; // if pop is empty
    
    unsigned int l, sum = 0;
    for (l = 0; l < _nb_locus; ++l) {
        if (curPatchFreqs[l][0] == 1) ++sum;
    }
    
    return sum;
}

// ----------------------------------------------------------------------------------------
// getNbFixedLocus_ofPatch
// ----------------------------------------------------------------------------------------
/** get the number of fixed locus of this patch specified by the patchID
 * !! the allele frequencies have to be set before (not anymore tested) !!
 */
template<class SH>
unsigned int StatHandler<SH>::getNbFixedLocus_ofPatch_andLocus(const age_idx & AGE, Patch * p, const unsigned int& l) {
    
    assert(p->get_sampleID()==my_NAN);
    if(p->get_sampleID()==SAMPLED) return my_NAN;
    set_alleleFreq(AGE);
    return (unsigned int)(_alleleFreq_local[AGE][p->get_sampleID()][l][0]==1);
}

// ----------------------------------------------------------------------------------------
// getNbFixedLocus
// ----------------------------------------------------------------------------------------
/** get the number of fixed locus averaged across all sampled populated patches
 */
template<class SH>
double StatHandler<SH>::getNbFixedLocus(const age_t & AGE)
{
    // get the local allele frequencies if necessary
    age_idx age_pos = age_t2idx(AGE);
    if (!set_alleleFreq(age_pos))	return my_NAN;
    
    unsigned int curNb, sum = 0;
    unsigned int nbFullPops = 0;
    
    vector<Patch*>::iterator curPop = get_vSamplePatch().begin();
    vector<Patch*>::iterator endPop = get_vSamplePatch().end();
    for (; curPop != endPop; ++curPop) { // for each patch
        curNb = getNbFixedLocus_ofPatch(*curPop, age_pos);
        if (curNb == my_NAN) continue; // if patch is empty
        sum += curNb;
        ++nbFullPops;
    }
    return(double)sum / nbFullPops;
}

// ----------------------------------------------------------------------------------------
// getNbFixedLocus
// ----------------------------------------------------------------------------------------
/** get the number of fixed locus averaged across all sampled populated patches
 */
template<class SH>
double StatHandler<SH>::getNbFixedLocus_ofLocus(unsigned int l, const age_t & AGE)
{
    // get the local allele frequencies if necessary
    age_idx age_pos = age_t2idx(AGE);
    unsigned int curNb, sum = 0;
    unsigned int nbFullPops = 0;
    
    vector<Patch*>::iterator curPop = get_vSamplePatch().begin();
    vector<Patch*>::iterator endPop = get_vSamplePatch().end();
    for (; curPop != endPop; ++curPop) { // for each patch
        curNb = getNbFixedLocus_ofPatch_andLocus(age_pos, *curPop, l);
        if (curNb == my_NAN) continue; // if patch is empty
        sum += curNb;
        ++nbFullPops;
    }
    return(double)sum / nbFullPops;
}

// ----------------------------------------------------------------------------------------
// getNbFixedLocusTot
// ----------------------------------------------------------------------------------------
/** get the number of fixed loci of the entire metapopulation
 * double is returned since all stats have to return a double
 */
template<class SH>
unsigned int StatHandler<SH>::getNbFixedLocusTot(const age_idx & AGE)
{
    if (!set_alleleFreq(AGE)) return my_NAN;
    
    unsigned int l, fix_loc_global = 0;
    for (l = 0; l < _nb_locus; ++l) {
        if (_alleleFreq_global[AGE][l].size() == 1) ++fix_loc_global;
    }
    return fix_loc_global;
}

// ----------------------------------------------------------------------------------------
// getNbFixedLocusTot
// ----------------------------------------------------------------------------------------
/** get the number of fixed loci of the entire metapopulation
 * double is returned since all stats have to return a double
 */
template<class SH>
unsigned int StatHandler<SH>::getNbFixedLocusTot_ofLocus(const age_idx & AGE, const unsigned int& l)
{
    if (!set_alleleFreq(AGE)) return my_NAN;
    return (unsigned int) (_alleleFreq_global[AGE][l].size() == 1);
}

// ----------------------------------------------------------------------------------------
// getHarmonicMean_ofPopSize
// ----------------------------------------------------------------------------------------
/** compute harmonic mean of the popualtion sizes of the given patches
 * the number of populated patches is adapted (last parameter "nbPopFull")!
 */
template<class SH>
double StatHandler<SH>::getHarmonicMean_ofPopSize(const age_idx & AGE, const vector<Patch*>&aPatch,
                                                  unsigned int&nbPopFull) {
    double harmonic = 0;
    
    unsigned int nbind;
    
    nbPopFull = 0;
    vector<Patch*>::const_iterator curPop, endPop;
    for (curPop = aPatch.begin(), endPop = aPatch.end(); curPop != endPop; ++curPop) {
        nbind = (*curPop)->sampleSize(AGE);
        if (!nbind) continue; // if the patch is empty
        harmonic += 1.0 / nbind;
        ++nbPopFull;
    }
    if (!nbPopFull) return(double)my_NAN;
    return(double)nbPopFull/harmonic;
}

// ----------------------------------------------------------------------------------------
// setHo
// ----------------------------------------------------------------------------------------
/** get mean observed heterozygosity */
template<class SH>double StatHandler<SH>::getHo(const age_idx & AGE) {
    // compute the observed heterozygosity for each locus
    getHo_perLocus(AGE);
    
    // if the metapop is empty
    if (_ho_locus[AGE][0] == my_NAN) {
        _ho[AGE] = my_NAN;
        return _ho[AGE];
    }
    
    // compute the mean observed heterozygosity
    _ho[AGE] = 0;
    for (unsigned int k = 0; k < _nb_locus; ++k) {
        _ho[AGE] += _ho_locus[AGE][k];
    }
    
    _ho[AGE] /= _nb_locus;
    
    return _ho[AGE];
}

// ----------------------------------------------------------------------------------------
// setHo
// ----------------------------------------------------------------------------------------
/** get observed heterozygoxity for each loci individually */
template<class SH>double* StatHandler<SH>::getHo_perLocus(const age_idx & AGE)
{
    // check if the table has already been computed
    if (already_computed(_computed[10], AGE)) return _ho_locus[AGE];
    
    if (!_ho_locus) ARRAY::create_1D<double*>(_ho_locus, NB_AGE_CLASSES, NULL); // first time
    if (_ho_locus[AGE]) ARRAY::reset_1D(_ho_locus[AGE], _nb_locus, (double)0);
    else ARRAY::create_1D(_ho_locus[AGE], _nb_locus, (double)0);
    
    double* ho = new double[_nb_locus]; // temp array
    
    // compute Ho per patch
    vector<Patch*>::iterator curPop, endPop = get_vSamplePatch().end();
    unsigned int l, nbFullPatch = 0;
    for (curPop = get_vSamplePatch().begin(); curPop != endPop; ++curPop) {
        getHo_ofPatch_andLocus(AGE, *curPop, ho); // compute ho for each locus of the patch i
        if (ho[0] == my_NAN) continue; // check if patch is empty
        for (l = 0; l < _nb_locus; ++l) {
            _ho_locus[AGE][l] += ho[l]; // sum up the ho of each patch
        }
        ++nbFullPatch;
    }
    
    if (!nbFullPatch) ARRAY::reset_1D(_ho_locus[AGE], _nb_locus, (double)my_NAN);
    // if metapop is empty
    else {
        for (l = 0; l < _nb_locus; ++l) {
            _ho_locus[AGE][l] /= nbFullPatch; // compute mean across patches
        }
    }
    
    delete[]ho;
    
    return _ho_locus[AGE];
}

// ----------------------------------------------------------------------------------------
// setHo
// ----------------------------------------------------------------------------------------
/** get observed heterozygosity of the given patch */
template<class SH>double StatHandler<SH>::getHo_ofPatch(Patch * curPatch,
                                                        const age_idx & AGE) {
    assert(curPatch->get_sampleID() != my_NAN);
    if (curPatch->get_sampleID() == SAMPLED)
        return my_NAN; // if currently not populated
    
    // compute the observed heterozygosity for each locus of that patch
    double* ho_locus = getHo_ofPatch_andLocus(AGE, curPatch);
    
    // if the metapop is empty
    if (ho_locus[0] == my_NAN) {
        delete[]ho_locus;
        return my_NAN;
    }
    
    // compute the mean observed heterozygosity
    double ho = 0;
    for (unsigned int l = 0; l < _nb_locus; ++l) {
        ho += ho_locus[l];
    }
    
    delete[]ho_locus;
    
    return ho /= _nb_locus;
}

// ----------------------------------------------------------------------------------------
// setHo
// ----------------------------------------------------------------------------------------
/** get observed heterozygoxity for each locus of the selected patch individually
 * if no array is passed then the returned array has to be deleted !
 */
template<class SH>
double* StatHandler<SH>::getHo_ofPatch_andLocus(const age_idx & AGE, Patch * cur_patch, double*array)
{
    unsigned int l;
    unsigned char** genes;
    
    vector<TIndividual*>&curFem = cur_patch->get_sampled_inds(FEM, AGE);
    vector<TIndividual*>&curMal = cur_patch->get_sampled_inds(MAL, AGE);
    unsigned int size = (unsigned int)curFem.size() + (unsigned int)curMal.size();
    
    // if patch is empty
    if (size){
        if (!array) 	ARRAY::create_1D(array, _nb_locus, (double)0); // create the array if not passed
        else ARRAY::reset_1D(array, _nb_locus, (double)0);
    }
    else{                                                      // patch is empty
        if (!array) 	ARRAY::create_1D(array, _nb_locus, (double)my_NAN); // create the array if not passed
        else ARRAY::reset_1D(array, _nb_locus, (double)my_NAN);
        return array;
    }
    
    // count heterozygote females
    vector<TIndividual*>::iterator curInd, endInd;
    for (curInd = curFem.begin(), endInd = curFem.end(); curInd != endInd; ++curInd) {
        genes = (unsigned char**)(*curInd)->getTrait(_SHLinkedTraitIndex)->get_sequence();
        for (l = 0; l < _nb_locus; ++l) {
            array[l] += (genes[l][0] != genes[l][1]);
        }
    }
    
    // count heterozygote males
    for (curInd = curMal.begin(), endInd = curMal.end(); curInd != endInd;++curInd) {
        genes = (unsigned char**)(*curInd)->getTrait(_SHLinkedTraitIndex)->get_sequence();
        for (l = 0; l < _nb_locus; ++l) {
            array[l] += (genes[l][0] != genes[l][1]);
        }
    }
    
    // compute the mean heterozygosity for each locus
    for (l = 0; l < _nb_locus; ++l) {
        array[l] /= size;
    }
    
    return array;
}

// ----------------------------------------------------------------------------------------
// getHo_ofPatchperLocus
// ----------------------------------------------------------------------------------------
/** get observed heterozygoxity for the fgiven locus and patch
 */
template<class SH> double
StatHandler<SH>::getHo_ofPatch_andLocus(const age_idx & AGE, Patch * cur_patch, const unsigned int&l)
{
    unsigned char** genes;
    
    vector<TIndividual*>&curFem = cur_patch->get_sampled_inds(FEM, AGE);
    vector<TIndividual*>&curMal = cur_patch->get_sampled_inds(MAL, AGE);
    unsigned int ho = 0, size = (unsigned int)curFem.size() + (unsigned int)curMal.size();
    
    // if patch is empty
    if (!size)
        return my_NAN; // if patch is empty
    
    // count heterozygote females
    vector<TIndividual*>::iterator curInd, endInd;
    for (curInd = curFem.begin(), endInd = curFem.end(); curInd != endInd; ++curInd) {
        genes = (unsigned char**)(*curInd)->getTrait(_SHLinkedTraitIndex)->get_sequence();
        ho += (genes[l][0] != genes[l][1]);
    }
    
    // count heterozygote males
    for (curInd = curMal.begin(), endInd = curMal.end(); curInd != endInd; ++curInd) {
        genes = (unsigned char**)(*curInd)->getTrait(_SHLinkedTraitIndex)->get_sequence();
        ho += (genes[l][0] != genes[l][1]);
    }
    
    return(double)ho / size;
}

// ----------------------------------------------------------------------------------------
/** get observed heterozygosity of the patch for each allele separately array[locus][allele]
 * if no array is passed then the returned array has to be deleted !
 */
template<class SH>map<unsigned char,
double> *StatHandler<SH>::getHo_ofPatchperAllele(const age_idx & AGE, Patch * cur_patch,
                                                 map<unsigned char, double> *array)
{
    if (array) delete[]array;
    array = new map<unsigned char, double>[_nb_locus];
    
    unsigned int l;
    unsigned char** genes;
    
    vector<TIndividual*>&curFem = cur_patch->get_sampled_inds(FEM, AGE);
    vector<TIndividual*>&curMal = cur_patch->get_sampled_inds(MAL, AGE);
    unsigned int size = (unsigned int)curFem.size() + (unsigned int)curMal.size();
    
    // if patch is empty
    if (!size) return array; // if patch is empty
    
    // count heterozygote females
    vector<TIndividual*>::iterator curInd, endInd;
    for (curInd = curFem.begin(), endInd = curFem.end(); curInd != endInd; ++curInd) {
        genes = (unsigned char**)(*curInd)->getTrait(_SHLinkedTraitIndex)->get_sequence();
        for (l = 0; l < _nb_locus; ++l) {
            if (genes[l][0] != genes[l][1]) {
                ++array[l][genes[l][0]]; // auto-initialization with zero
                ++array[l][genes[l][1]]; // auto-initialization with zero
            }
        }
    }
    
    // count heterozygote males
    for (curInd = curMal.begin(), endInd = curMal.end(); curInd != endInd; ++curInd) {
        genes = (unsigned char**)(*curInd)->getTrait(_SHLinkedTraitIndex)->get_sequence();
        for (l = 0; l < _nb_locus; ++l) {
            if (genes[l][0] != genes[l][1]) {
                ++array[l][genes[l][0]]; // auto-initialization with zero
                ++array[l][genes[l][1]]; // auto-initialization with zero
            }
        }
    }
    
    // compute the mean heterozygosity for each locus
    map<unsigned char, double>::iterator pos, end;
    for (l = 0; l < _nb_locus; ++l) { // for each locus
        for (pos = array[l].begin(), end = array[l].end(); pos != end; ++pos) {
            // for each allele
            pos->second /= size;
        }
    }
    
    return array;
}

// ----------------------------------------------------------------------------------------
// setHo
// ----------------------------------------------------------------------------------------
/** get observed heterozygoxity for each locus individually
 * Ho = mean(sum_patch(a1 != a2)/nb_ind_of_patch)
 * the returned array must be deleted!
 */
template<class SH>double* StatHandler<SH>::getHo_perLocus(const age_idx & AGE,
                                                          vector<Patch*>&vPatch) {
    // check if the table has already been computed
    if (already_computed(_computed[11], AGE)) return _ho_locus;
    
    unsigned int size, nbFullPatch = 0, k;
    unsigned int* ho = new unsigned int[_nb_locus]; // temp parameter
    double* ho_locus = new double[_nb_locus];
    // is returned and must be deleted in the calling function
    unsigned char** genes;
    
    vector<TIndividual*>::iterator curInd, endInd;
    
    // compute Ho per patch
    vector<Patch*>::iterator curPop, endPop;
    for (curPop = vPatch.begin(), endPop = vPatch.end();
         curPop != endPop; ++curPop) {
        ARRAY::reset_1D(ho, _nb_locus, (unsigned int)0); // reset all ho to zero
        vector<TIndividual*>&curFem = (*curPop)->get_sampled_inds(FEM, AGE);
        vector<TIndividual*>&curMal = (*curPop)->get_sampled_inds(MAL, AGE);
        size = (unsigned int)(curFem.size() + curMal.size());
        if (!size)
            continue; // if patch is empty
        ++nbFullPatch;
        
        // females
        for (curInd = curFem.begin(), endInd = curFem.end(); curInd != endInd;
             ++curInd) {
            genes = (unsigned char**)(*curInd)->getTrait(_SHLinkedTraitIndex)
            ->get_sequence();
            for (k = 0; k < _nb_locus; ++k) {
                ho[k] += (genes[k][0] != genes[k][1]);
            }
        }
        
        // males
        for (curInd = curMal.begin(), endInd = curMal.end(); curInd != endInd;
             ++curInd) {
            genes = (unsigned char**)(*curInd)->getTrait(_SHLinkedTraitIndex)
            ->get_sequence();
            for (k = 0; k < _nb_locus; ++k) {
                ho[k] += (genes[k][0] != genes[k][1]);
            }
        }
        
        // sum up the mean ho per locus for each patch
        for (k = 0; k < _nb_locus; ++k) {
            _ho_locus[k] += (double)ho[k] / size;
        }
    }
    
    // compute mean across patches
    for (k = 0; k < _nb_locus; ++k) {
        _ho_locus[k] /= nbFullPatch;
    }
    
    delete[]ho;
    return ho_locus;
}

// ----------------------------------------------------------------------------------------
// getHs
// ----------------------------------------------------------------------------------------
/** return mean expected heterozygosity (across loci and patches) hs = 1 - sum(p^2) */
template<class SH>double StatHandler<SH>::getHs(const age_idx & AGE) {
    // check if the table has already been computed
    if (already_computed(_computed[12], AGE))
        return _hs[AGE];
    
    // compute mean expected heterozygosity
    _hs[AGE] = 0;
    for (unsigned int l = 0; l < _nb_locus; ++l) {
        _hs[AGE] += getHs_ofLocus(l, AGE);
    }
    _hs[AGE] /= _nb_locus;
    
    return _hs[AGE];
}

// ----------------------------------------------------------------------------------------
// getHs_perPatch_perLocus
// ----------------------------------------------------------------------------------------
/** return mean expected heterozygosity of locus l across patches: hs = 1 - sum(p^2) */
template<class SH>double StatHandler<SH>::getHs_ofLocus(const unsigned int&l, const age_idx & AGE)
{
    if (!set_alleleFreq(AGE)) return my_NAN;
    
    unsigned int nbFullPatch = 0;
    double val, hs = 0;
    vector<Patch*>::iterator curPop, endPop = get_vSamplePatch().end();
    for (curPop = get_vSamplePatch().begin(); curPop != endPop; ++curPop) {
        val = getHs_ofPatch_andLocus(AGE, *curPop, l);
        if (val == my_NAN) continue;         // patch is empty
        hs += val;
        ++nbFullPatch;
    }
    
    return nbFullPatch ? hs / nbFullPatch : my_NAN;
}

// ----------------------------------------------------------------------------------------
// getHs_perPatch_perPatch
// ----------------------------------------------------------------------------------------
/** return mean expected heterozygosity of patch p across loci: hs = 1 - sum(p^2) */
template<class SH> double
StatHandler<SH>::getHs_ofPatch(Patch * p,const age_idx & AGE)
{
    if (!set_alleleFreq(AGE)) return my_NAN;
    
    double hs = 0;
    for (unsigned int l = 0; l < _nb_locus; ++l) {
        hs += getHs_ofPatch_andLocus(AGE, p, l);
    }
    return hs / _nb_locus;
}

// ----------------------------------------------------------------------------------------
// getHs_perPatch
// ----------------------------------------------------------------------------------------
/** return mean UNBIASED expected heterozygosity of patch p across loci (Nei 1987, eq. 7.39, p. 164) */
template<class SH> double
StatHandler<SH>::getHsUnbiased_ofPatch(unsigned int p, const age_t & AGE)
{
    age_idx age_pos = age_t2idx(AGE);
    Patch* curPatch = get_vPatch(p);
    unsigned int size = curPatch->sampleSize(age_pos);
    if (size < 2) return my_NAN;
    return (double)size / (size - 1) * (getHs_ofPatch(curPatch, age_pos) - (getHo_ofPatch(curPatch, age_pos) / (2 * size)));
}

// ----------------------------------------------------------------------------------------
// getHs_perPatch
// ----------------------------------------------------------------------------------------
/** return mean UNBIASED expected heterozygosity of patch p across loci (Nei 1987, eq. 7.39, p. 164) */
template<class SH> double
StatHandler<SH>::getHsUnbiased_ofPatch_andLocus(const age_idx & AGE, Patch * p, const unsigned int&l)
{
    unsigned int size = p->sampleSize(AGE);
    if (size < 2) return my_NAN;
    return(double)size / (size - 1) * (getHs_ofPatch_andLocus(AGE, p, l) - (getHo_ofPatch_andLocus(AGE, p, l) / (2 * size)));
}

// ----------------------------------------------------------------------------------------
// getHs_perPatchandLocus
// ----------------------------------------------------------------------------------------
/** return expected heterozygosity of patch p and locus l: hs = 1 - sum(p^2)
 * returns NaN if the pop is empty
 */
template<class SH> double
StatHandler<SH>::getHs_ofPatch_andLocus(const age_idx & AGE, Patch * p, const unsigned int&l)
{
    assert(set_alleleFreq(AGE));
    assert(p->get_sampleID()!=my_NAN);
    if(p->get_sampleID()==SAMPLED) return my_NAN;
    
    double hs = 1;
    map<unsigned char, double>::iterator pos, end;
    pos = _alleleFreq_local[AGE][p->get_sampleID()][l].begin();
    end = _alleleFreq_local[AGE][p->get_sampleID()][l].end();
    for (; pos != end; ++pos) {
        hs -= pos->second * pos->second;
    }
    return hs==1 ? (double) my_NAN : hs;
}

// ----------------------------------------------------------------------------------------
// getHt
// ----------------------------------------------------------------------------------------
/** return mean total expected heterozygosity (across loci) ht = 1 - sum(p^2) */
template<class SH> double
StatHandler<SH>::getHt(const age_idx & AGE)
{
    // check if the table has already been computed
    if (already_computed(_computed[13], AGE)) return _ht[AGE];
    
    // Ht is NaN if only a single pop is populated
    if (get_nbSamplePatch(AGE) < 2) {
        _ht[AGE] = my_NAN;
        return _ht[AGE];
    } // 2 or more pops are needed
    
    _ht[AGE] = 0;
    for (unsigned int l = 0; l < _nb_locus; ++l) {
        _ht[AGE] += getHt_ofLocus(AGE, l);
    }
    _ht[AGE] /= _nb_locus;
    
    return _ht[AGE];
}

// ----------------------------------------------------------------------------------------
// getHt
// ----------------------------------------------------------------------------------------
/** return mean total expected heterozygosity (across loci) ht = 1 - sum(mean(p)^2)
 *  across 2 patches
 */
template<class SH>
double StatHandler<SH>::getHt(const age_idx & AGE, Patch * p1, Patch * p2)
{
    double ht = 0;
    assert(p1->sampleSize(AGE) && p2->sampleSize(AGE));
    for (unsigned int l = 0; l < _nb_locus; ++l) {
        ht += getHt_ofLocus(AGE, l, p1, p2);
    }
    return ht /= _nb_locus;
}

// ----------------------------------------------------------------------------------------
// getHt_perPatch_perLocus
// ----------------------------------------------------------------------------------------
/** return mean expected heterozygosity of locus l across all patches: hs = 1 - sum(mean(p)^2)
 * used by Nei
 */
template<class SH> double
StatHandler<SH>::getHt_ofLocus(const age_idx & AGE,	const unsigned int&l)
{
    // get the local allele frequencies if necessary
    if (!set_alleleFreq(AGE)) return my_NAN;
    
    double ht = 1;
    map<unsigned char, double>::iterator pos, end;
    pos = _alleleFreq_global[AGE][l].begin();
    end = _alleleFreq_global[AGE][l].end();
    for (; pos != end; ++pos) {
        ht -= pos->second * pos->second;
    }
    
    return ht;
}

// ----------------------------------------------------------------------------------------
/** return mean expected heterozygosity of locus l across patches: ht = 1 - sum(mean(p)^2)
 * used by Nei for the given patches
 */
template<class SH>
double StatHandler<SH>::getHt_ofLocus(const age_idx & AGE,const unsigned int&l, Patch * p1, Patch * p2)
{
    set_alleleFreq(AGE);
    
    // get the pop sizes
    unsigned int size[2];
    
    size[0] = p1->sampleSize(AGE);
    if (!size[0]) return my_NAN;
    size[1] = p2->sampleSize(AGE);
    if (!size[1]) return my_NAN;
    unsigned int tot_size = size[0] + size[1];
    
    // compute the global allele frequencies for the two patches
    map<unsigned char, double>allFreqGlobal;
    map<unsigned char, double> *allFreqLocal[2];
    
    allFreqLocal[0] = _alleleFreq_local[AGE][p1->get_sampleID()];
    allFreqLocal[1] = _alleleFreq_local[AGE][p2->get_sampleID()];
    set_alleleFreq_global_ofLocus(l, allFreqLocal, &allFreqGlobal, size, tot_size, 2);
    
    // compute the Ht
    double ht = 1;
    map<unsigned char, double>::iterator pos, end;
    for (pos = allFreqGlobal.begin(), end = allFreqGlobal.end(); pos != end;
         ++pos) {
        ht -= pos->second * pos->second;
    }
    
    return ht;
}

// ----------------------------------------------------------------------------------------
// setFstat_Weir_Cockerham
// ----------------------------------------------------------------------------------------
/** computes the global F-statistics following Weir and Cockerham (1984) */
template<class SH> void
StatHandler<SH>::setFstat_Weir_Cockerham(const age_idx & AGE)
{
    // check if the table has already been computed
    if (already_computed(_computed[14], AGE)) return;
    
    if (!set_alleleFreq(AGE)) { // all pops are emtpy
        _fst_wc[AGE] = _fis_wc[AGE] = _fit_wc[AGE] = my_NAN;
    }
    
    setFstat_Weir_Cockerham(AGE, get_vSamplePatch(), _alleleFreq_local[AGE],
                            _alleleFreq_global[AGE], _fst_wc[AGE], _fis_wc[AGE], _fit_wc[AGE]);
}

// ----------------------------------------------------------------------------------------
// setFstat_Weir_Cockerham_perLocus
// ----------------------------------------------------------------------------------------
/** computes the global F-statistics following Weir and Cockerham for each locus (1984) */
template<class SH> void
StatHandler<SH>::setFstat_Weir_Cockerham_perLocus(const age_idx & AGE)
{
    // check if the table has already been computed
    if (already_computed(_computed[15], AGE)) return;
    
    if (!_fis_WC_locus) { // first time
        ARRAY::create_1D<double*>(_fis_WC_locus, NB_AGE_CLASSES, NULL);
        ARRAY::create_1D<double*>(_fit_WC_locus, NB_AGE_CLASSES, NULL);
        ARRAY::create_1D<double*>(_fst_WC_locus, NB_AGE_CLASSES, NULL);
    }
    if (!_fis_WC_locus[AGE]) {
        _fis_WC_locus[AGE] = new double[_nb_locus];
        _fit_WC_locus[AGE] = new double[_nb_locus];
        _fst_WC_locus[AGE] = new double[_nb_locus];
    }
    
    set_alleleFreq(AGE);
    
    unsigned int* pop_sizes = new unsigned int[get_current_nbSamplePatch()];
    map<unsigned char, double> **ho_patch_allele = new map<unsigned char,
    double> *[get_current_nbSamplePatch()]; // ho_patch_allele[p][l][a]
    unsigned int l, tot_size = 0, cur_size, p, nbFullPatch = 0;
    double a2, b2, w2, tot_square_size = 0, nc;
    
    // get pop sizes
    vector<Patch*>::iterator curPop, endPop = get_vSamplePatch().end();
    for (p = 0, curPop = get_vSamplePatch().begin(); curPop != endPop;
         ++curPop, ++p) {
        pop_sizes[p] = cur_size = (*curPop)->sampleSize(AGE);
        if (!cur_size) {
            ho_patch_allele[p] = NULL;
            continue; // if pop is empty
        }
        tot_size += cur_size;
        ho_patch_allele[p] = getHo_ofPatchperAllele(AGE, *curPop);
        tot_square_size += cur_size * cur_size;
        ++nbFullPatch;
    }
    
    if (nbFullPatch < 2) {
        for (l = 0; l < _nb_locus; ++l) {
            _fis_WC_locus[AGE][l] = _fit_WC_locus[AGE][l] = _fst_WC_locus[AGE][l]
            = my_NAN;
        }
        return;
    }
    
    nc = (tot_size - (tot_square_size / tot_size)) / (nbFullPatch - 1);
    
    for (l = 0; l < _nb_locus; ++l) { // for each locus
        a2 = b2 = w2 = 0;
        
        get_sigma_of_locus_Weir_Cockerham(a2, b2, w2, get_current_nbSamplePatch(),
                                          nbFullPatch, _alleleFreq_local[AGE], _alleleFreq_global[AGE], l,
                                          tot_size, pop_sizes, ho_patch_allele, nc);
        
        if (a2 || b2 || w2) {
            _fst_WC_locus[AGE][l] = a2 / (a2 + b2 + w2);
            _fit_WC_locus[AGE][l] = (a2 + b2) / (a2 + b2 + w2);
        }
        else _fst_WC_locus[AGE][l] = _fit_WC_locus[AGE][l] = my_NAN;
        
        if (b2 || w2) _fis_WC_locus[AGE][l] = b2 / (b2 + w2);
        else          _fis_WC_locus[AGE][l] = my_NAN;
    }
    
    delete[]pop_sizes;
    
    ARRAY::delete_2D(ho_patch_allele, get_current_nbSamplePatch());
}

// ----------------------------------------------------------------------------------------
// setFstat_Weir_Cockerham_perLocus
// ----------------------------------------------------------------------------------------
/** computes the pairwise Fst following Weir and Cockerham (1984).
 * all pairs are computed
 * Values are stored in a matrix
 */
template<class SH> void
StatHandler<SH>::setFstat_Weir_Cockerham_perPatchPair(const age_idx & AGE)
{
    // check if the table has already been computed
    if (already_computed(_computed[23], AGE)) return;
    set_alleleFreq(AGE);
    if (!_fst_matrix_wc) ARRAY::create_1D<double**>(_fst_matrix_wc, NB_AGE_CLASSES, NULL);
    
    // first time
    if (!_fst_matrix_wc[AGE]) ARRAY::create_2D(_fst_matrix_wc[AGE], get_current_nbSamplePatch(), get_current_nbSamplePatch(), (double)my_NAN);
    else if (get_current_nbSamplePatch() != get_last_nbSamplePatch()) {
        ARRAY::delete_2D(_fst_matrix_wc[AGE], get_last_nbSamplePatch());
        ARRAY::create_2D(_fst_matrix_wc[AGE], get_current_nbSamplePatch(),	get_current_nbSamplePatch(), (double)my_NAN);
    }
    else ARRAY::reset_2D(_fst_matrix_wc[AGE], get_current_nbSamplePatch(), get_current_nbSamplePatch(), (double)my_NAN);
    
    unsigned int* pop_sizes = new unsigned int[get_current_nbSamplePatch()];
    map<unsigned char, double> **ho_patch_allele = new map<unsigned char, double> *[get_current_nbSamplePatch()]; // ho_patch_allele[p][l][a]
    unsigned int tot_size, l, i, j, cur_pop_sizes[2];
    double a2, b2, w2, tot_square_size, nc;
    map<unsigned char, double> *cur_ho_patch_allele[2];
    map<unsigned char, double> *cur_allele_freq_local[2];
    map<unsigned char, double> *cur_allele_freq_global = new map<unsigned char, double>[_nb_locus];
    
    // get pop sizes and Ho
    vector<Patch*>::iterator curPop1, curPop2, endPop = get_vSamplePatch().end();
    for (i = 0, curPop1 = get_vSamplePatch().begin(); curPop1 != endPop; ++curPop1, ++i) {
        assert((*curPop1)->get_sampleID() == i);
        pop_sizes[i] = (*curPop1)->sampleSize(AGE);
        if (pop_sizes[i]) ho_patch_allele[i] = getHo_ofPatchperAllele(AGE, *curPop1);
        else              ho_patch_allele[i] = NULL;
    }
    
    // for each pair of pops
    for (i = 0, curPop1 = get_vSamplePatch().begin(); curPop1 != endPop; ++curPop1, ++i) { // first patch
        assert((*curPop1)->get_sampleID() == i);
        if (!pop_sizes[i]) continue; // both pops have to be populated
        cur_ho_patch_allele[0] = ho_patch_allele[i];
        cur_pop_sizes[0] = pop_sizes[i];
        cur_allele_freq_local[0] = _alleleFreq_local[AGE][(*curPop1)->get_sampleID()];
        
        curPop2 = curPop1;
        for (j = i + 1, ++curPop2; curPop2 != endPop; ++curPop2, ++j) { // second patch
            assert(j == (*curPop2)->get_sampleID());
            if (!pop_sizes[j]) continue; // both pops have to be populated!
            cur_ho_patch_allele[1] = ho_patch_allele[j];
            cur_pop_sizes[1] = pop_sizes[j];
            cur_allele_freq_local[1] = _alleleFreq_local[AGE][(*curPop2)->get_sampleID()];
            
            tot_size = pop_sizes[i] + pop_sizes[j];
            tot_square_size = pop_sizes[i] * pop_sizes[i] + pop_sizes[j] * pop_sizes[j];
            
            // compute the global allele frequencies
            set_alleleFreq_global(cur_allele_freq_local, cur_allele_freq_global, cur_pop_sizes, tot_size, 2);
            
            nc = (tot_size - (tot_square_size / tot_size)); // /(nbFullPatch-1);  = /1
            
            // loop over each locus
            a2 = b2 = w2 = 0;
            for (l = 0; l < _nb_locus; ++l) { // for each locus
                get_sigma_of_locus_Weir_Cockerham(a2, b2, w2, 2, 2,
                                                  cur_allele_freq_local, cur_allele_freq_global, l, tot_size,
                                                  cur_pop_sizes, cur_ho_patch_allele, nc);
            }
            
            if (a2 || b2 || w2) _fst_matrix_wc[AGE][i][j] = a2 / (a2 + b2 + w2);
        }
    }
    
    delete[]pop_sizes;
    ARRAY::delete_2D(ho_patch_allele, get_current_nbSamplePatch());
    delete[]cur_allele_freq_global;
}

// ----------------------------------------------------------------------------------------
// setFstat_Weir_Cockerham_perLocus
// ----------------------------------------------------------------------------------------
/** computes the Fst between the two specified pops following Weir and Cockerham (1984)
 * The Fst is computed just for this pair of patches.
 * Values are not stored in a matrix
 */
template<class SH>
double StatHandler<SH>::getFstat_Weir_Cockerham_perPatchPair(const age_idx & AGE, Patch * pop1, Patch * pop2)
{
    if (!set_alleleFreq(AGE)) return my_NAN;
    
    // get pop sizes
    unsigned int nbPop = 2;
    
    unsigned int size[2];
    
    assert(pop1->get_sampleID() != my_NAN);
    size[0] = pop1->sampleSize(AGE);
    if (!size[0]) return my_NAN;
    
    assert(pop1->get_sampleID() != my_NAN);
    size[1] = pop2->sampleSize(AGE);
    if (!size[1]) return my_NAN;
    unsigned int tot_size = size[0] + size[1];
    unsigned int tot_square_size = size[0] * size[0] + size[1] * size[1];
    
    // get Ho
    map<unsigned char, double> *ho[2];
    
    ho[0] = getHo_ofPatchperAllele(AGE, pop1);
    ho[1] = getHo_ofPatchperAllele(AGE, pop2);
    
    // allele freqs
    map<unsigned char, double> *allFreqLocal[2];
    
    allFreqLocal[0] = _alleleFreq_local[AGE][pop1->get_sampleID()];
    allFreqLocal[1] = _alleleFreq_local[AGE][pop2->get_sampleID()];
    
    // compute the global allele frequencies
    map<unsigned char, double> *allFreqGlobal = new map<unsigned char,
    double>[_nb_locus];
    set_alleleFreq_global(allFreqLocal, allFreqGlobal, size, tot_size, nbPop);
    
    double nc = (tot_size - (tot_square_size / tot_size)) / (nbPop - 1);
    
    // loop over each locus
    double a2 = 0, b2 = 0, w2 = 0;
    for (unsigned int l = 0; l < _nb_locus; ++l) { // for each locus
        get_sigma_of_locus_Weir_Cockerham(a2, b2, w2, 2, nbPop, allFreqLocal,
                                          allFreqGlobal, l, tot_size, size, ho, nc);
    }
    
    delete[]allFreqGlobal;
    delete[]ho[0];
    delete[]ho[1];
    
    if (a2 || b2 || w2)
        return a2 / (a2 + b2 + w2);
    
    return my_NAN;
}

// ----------------------------------------------------------------------------------------
// getFstat_Weir_Cockerham_perPatchPair_andLocus
// ----------------------------------------------------------------------------------------
/** computes the Fst between the two specified pops and locus following Weir and Cockerham (1984)
 * The Fst is computed just for this pair of patches and the specified locus.
 * Values are not stored in a matrix
 */
template<class SH>
double StatHandler<SH>::getFstat_Weir_Cockerham_perPatchPair_andLocus(const age_idx & AGE,
                                                                      Patch * pop1, Patch * pop2, unsigned int locus)
{
    if (!set_alleleFreq(AGE)) return my_NAN;
    
    // get pop sizes
    unsigned int nbPop = 2;
    
    unsigned int size[2];
    
    assert(pop1->get_sampleID() != my_NAN);
    size[0] = pop1->sampleSize(AGE);
    if (!size[0]) return my_NAN;
    
    assert(pop1->get_sampleID() != my_NAN);
    size[1] = pop2->sampleSize(AGE);
    if (!size[1]) return my_NAN;
    unsigned int tot_size = size[0] + size[1];
    unsigned int tot_square_size = size[0] * size[0] + size[1] * size[1];
    
    // get Ho
    map<unsigned char, double> *ho[2];
    
    ho[0] = getHo_ofPatchperAllele(AGE, pop1);
    ho[1] = getHo_ofPatchperAllele(AGE, pop2);
    
    // allele freqs
    map<unsigned char, double> *allFreqLocal[2];
    
    allFreqLocal[0] = _alleleFreq_local[AGE][pop1->get_sampleID()];
    allFreqLocal[1] = _alleleFreq_local[AGE][pop2->get_sampleID()];
    
    // compute the global allele frequencies
    map<unsigned char, double> *allFreqGlobal = new map<unsigned char,double>[_nb_locus];
    set_alleleFreq_global(allFreqLocal, allFreqGlobal, size, tot_size, nbPop);
    
    double nc = (tot_size - (tot_square_size / tot_size)) / (nbPop - 1);
    
    double a2 = 0, b2 = 0, w2 = 0;
    get_sigma_of_locus_Weir_Cockerham(a2, b2, w2, 2, nbPop, allFreqLocal,
                                      allFreqGlobal, locus, tot_size, size, ho, nc);
    
    delete[]allFreqGlobal;
    delete[]ho[0];
    delete[]ho[1];
    
    if (a2 || b2 || w2)
        return a2 / (a2 + b2 + w2);
    
    return my_NAN;
}

// ----------------------------------------------------------------------------------------
// setFstat_Weir_Cockerham
// ----------------------------------------------------------------------------------------
/** computes the global F-statistics following Weir and Cockerham (1984) for the given patches */
template<class SH>
void StatHandler<SH>::setFstat_Weir_Cockerham(const age_idx & AGE, vector<Patch*>&aPatch,
                                              map<unsigned char, double> **alleleFreqs,
                                              map<unsigned char, double> *alleleFreqsGlobal,
                                              double&fst, double&fis, double&fit)
{
    unsigned int* pop_sizes = new unsigned int[get_current_nbSamplePatch()];
    map<unsigned char, double> **ho_patch_allele = new map<unsigned char, double> *[get_current_nbSamplePatch()]; // ho_patch_allele[p][l][a]
    unsigned int tot_size = 0, cur_size, p, l, nbFullPatch = 0;
    double a2, b2, w2, tot_square_size = 0, nc;
    
    // get pop sizes
    vector<Patch*>::iterator curPop, endPop = get_vSamplePatch().end();
    for (p = 0, curPop = get_vSamplePatch().begin(); curPop != endPop; ++curPop, ++p) {
        pop_sizes[p] = cur_size = (*curPop)->sampleSize(AGE);
        if (!cur_size) {
            ho_patch_allele[p] = NULL;
            continue; // if pop is empty
        }
        ++nbFullPatch; // get the number of populated pops
        tot_size += cur_size;
        ho_patch_allele[p] = getHo_ofPatchperAllele(AGE, *curPop);
        tot_square_size += cur_size * cur_size;
    }
    
    if (nbFullPatch < 2) {
        fst = fis = fit = my_NAN;
        delete[]pop_sizes;
        ARRAY::delete_2D(ho_patch_allele, get_current_nbSamplePatch());
        return;
    }
    
    nc = (tot_size - (tot_square_size / tot_size)) / (nbFullPatch - 1);
    
    // compute sigma(a2), sigma(b2) and sigma(w2) for each locus and sum them up
    a2 = b2 = w2 = 0;
    for (l = 0; l < _nb_locus; ++l) {
        get_sigma_of_locus_Weir_Cockerham(a2, b2, w2, get_current_nbSamplePatch(),
                                          nbFullPatch, alleleFreqs, alleleFreqsGlobal, l, tot_size, pop_sizes,
                                          ho_patch_allele, nc);
    }
    
    if (a2 || b2 || w2) {
        fst = a2 / (a2 + b2 + w2);
        fit = (a2 + b2) / (a2 + b2 + w2);
    }
    else fst = _fit_wc[AGE] = my_NAN;
    
    if (b2 || w2) fis = b2 / (b2 + w2);
    else          fis = my_NAN;
    
    delete[]pop_sizes;
    
    ARRAY::delete_2D(ho_patch_allele, get_current_nbSamplePatch());
}

// ----------------------------------------------------------------------------------------
// setFstat_Weir_Cockerham
// ----------------------------------------------------------------------------------------
/** computes the three sigma for each locus separately (F-statistics following Weir and Cockerham (1984)) */
template<class SH>
void StatHandler<SH>::get_sigma_of_locus_Weir_Cockerham(double&sigma_a2, double&sigma_b2, double&sigma_w2,
                                                        const unsigned int&nbPatch,
                                                        const unsigned int&nbFullPatch,
                                                        map<unsigned char, double> **alleleFreqs,
                                                        map<unsigned char,double> *alleleFreqsGlobal,
                                                        const unsigned int&l,
                                                        const unsigned int&tot_size,
                                                        unsigned int*pop_sizes,
                                                        map<unsigned char, double> **ho_patch_allele,
                                                        const double&nc)
{
    double p_freq, sum1, sum2, nom, bloc, val;
    unsigned int p;
    
    map<unsigned char, double>::iterator pos_global, end_global, pos_local,
    end_local;
    pos_global = alleleFreqsGlobal[l].begin();
    end_global = alleleFreqsGlobal[l].end();
    for (; pos_global != end_global; ++pos_global) { // for each allele
        p_freq = pos_global->second; // the global allele frequency of allele a
        
        sum1 = tot_size * p_freq * (1 - p_freq);
        
        sum2 = 0;
        nom = 0;
        for (p = 0; p < nbPatch; ++p) { // for each patch
            if (!pop_sizes[p]) continue; // if population is emtpy
            pos_local = alleleFreqs[p][l].find(pos_global->first);
            
            // get local allele frequency of allele a in patch p
            if (pos_local == alleleFreqs[p][l].end()) val = -p_freq; // if not present in this population
            else{                                                    // if present in this population
                val = pos_local->second - p_freq;
                nom  += ho_patch_allele[p][l][pos_local->first] * pop_sizes[p];
            }
            sum2 += pop_sizes[p] * val * val;
        }
        
        if ((tot_size != nbFullPatch)) { // check if it can be computed!!!!
            bloc = (sum1 - sum2 - (nom / 4)) / (tot_size - nbFullPatch);
        }
        else {
            sigma_a2 = sigma_b2 = sigma_w2 = 0;
            cout << "\nFST-problems\n";
            return;
        }
        
        sigma_a2 += (sum2 / (nbFullPatch - 1) - bloc) / nc; // sigma(a)^2
        sigma_b2 += bloc - (nom / (4 * tot_size)); // sigma(b)^2
        sigma_w2 += nom / (2 * tot_size); // sigma(w)^2
    }
}

// ----------------------------------------------------------------------------------------
// getChordDist_perPatchPair
// ----------------------------------------------------------------------------------------
/** returns the genetic distance following Cavalli_Sforza & Bodmer (1971)
 * Dc=D^2=Sum_loci(1-Sum_allele(sqrt(p1*p2)))/(sum_loci(numAllele))
 *       =(numLoci - Sum_loci(Sum_allele(sqrt(p1*p2))))/(sum_loci(numAllele) - numLoci)
 */
template<class SH>
double StatHandler<SH>::getChordDist2_perPatchPair(const age_idx& AGE, Patch* pop1, Patch* pop2)
{
    if (!set_alleleFreq(AGE)) return my_NAN;
    
    assert(pop1->get_sampleID() != my_NAN);          // they have to be sampled at any time
    unsigned int sampleID1 = pop1->get_sampleID();
    if(sampleID1 == SAMPLED) return my_NAN;          // not sampled now
    
    assert(pop2->get_sampleID() != my_NAN);          // they have to be sampled at any time
    unsigned int sampleID2 = pop2->get_sampleID();
    if(sampleID2 == SAMPLED) return my_NAN;          // not sampled now
    assert(pop2->sampleSize(AGE));
    
    // allele freqs
    map<unsigned char, double> *allFreqs1 = _alleleFreq_local[AGE][sampleID1];
    map<unsigned char, double> *allFreqs2 = _alleleFreq_local[AGE][sampleID2];
    map<unsigned char, double>::iterator curAll1, curAll2, endAll1, endAll2;
    
    unsigned int l, nbAll=0;
    double product=0;
    for(l=0; l<_nb_locus; ++l){
        curAll1 = allFreqs1[l].begin();
        endAll1 = allFreqs1[l].end();
        curAll2 = allFreqs2[l].begin();
        endAll2 = allFreqs2[l].end();
        while(curAll1!=endAll1 && curAll2!=endAll2){
            if(curAll1->first==curAll2->first){		               // the allele is present in both pops
                product += sqrt(curAll1->second*curAll2->second);  // the product changes only if the allele occurs in both pops
                ++curAll1;
                ++curAll2;
            }
            else if(curAll1->first<curAll2->first) ++curAll1;    // allele is only present in pop1 (product of pi*pj is 0)
            else                                   ++curAll2;    // allele is only present in pop2 (product of pi*pj is 0)
            ++nbAll;                               // increment the number of alleles
        }
        // maybe the last few alleles just occur in one of the pops:
        for(; curAll1!=endAll1; ++curAll1) ++nbAll;// for all alleles just occuring in pop1 (product of pi*pj is 0)
        for(; curAll2!=endAll2; ++curAll2) ++nbAll;// for all alleles just occuring in pop2 (product of pi*pj is 0)
    }
    
    nbAll -= _nb_locus;			                   // correct the number of alleles
    if(!nbAll) return my_NAN;
    return (4.0*(_nb_locus - product))/nbAll;
}

// ----------------------------------------------------------------------------------------
// getChordDist_perPatchPair
// ----------------------------------------------------------------------------------------
/** returns the genetic distance following Cavalli_Sforza & Bodmer (1971)
 * Dc=2/(Pi*r)*Sum_loci(sqrt(2*(1-Sum_allele(sqrt(p1*p2)))))
 */
template<class SH>
double StatHandler<SH>::getChordDist_perPatchPair(const age_idx& AGE, Patch* pop1, Patch* pop2)
{
    if (!set_alleleFreq(AGE)) return my_NAN;
    
    assert(pop1->get_sampleID() != my_NAN);          // they have to be sampled at any time
    unsigned int sampleID1 = pop1->get_sampleID();
    if(sampleID1 == SAMPLED) return my_NAN;          // not sampled now
    
    assert(pop2->get_sampleID() != my_NAN);          // they have to be sampled at any time
    unsigned int sampleID2 = pop2->get_sampleID();
    if(sampleID2 == SAMPLED) return my_NAN;          // not sampled now
    assert(pop2->sampleSize(AGE));
    
    // allele freqs
    map<unsigned char, double> *allFreqs1 = _alleleFreq_local[AGE][sampleID1];
    map<unsigned char, double> *allFreqs2 = _alleleFreq_local[AGE][sampleID2];
    map<unsigned char, double>::iterator curAll1, curAll2, endAll1, endAll2;
    
    unsigned int l;
    double sum_all, sum_loc=0;
    for(l=0; l<_nb_locus; ++l){
        curAll1 = allFreqs1[l].begin();
        endAll1 = allFreqs1[l].end();
        curAll2 = allFreqs2[l].begin();
        endAll2 = allFreqs2[l].end();
        sum_all=0;
        while(curAll1!=endAll1 && curAll2!=endAll2){
            if(curAll1->first==curAll2->first){		               // the allele is present in both pops
                sum_all += sqrt(curAll1->second*curAll2->second);  // the product changes only if the allele occurs in both pops
                ++curAll1;
                ++curAll2;
            }
            else if(curAll1->first<curAll2->first) ++curAll1;    // allele is only present in pop1 (product of pi*pj is 0)
            else                                   ++curAll2;    // allele is only present in pop2 (product of pi*pj is 0)
        }
        sum_loc += sqrt(2*(1-sum_all));
    }
    
    return 2/PI/_nb_locus*sum_loc;
}

// ----------------------------------------------------------------------------------------
// getChordDist
// ----------------------------------------------------------------------------------------
/** returns the average chord distance following Cavalli_Sforza & Bodmer (1971)
 * Dc=D^2=4*Sum_loci(1-Sum_allele(sqrt(p1*p2)))/(sum_loci(numAllele - 1))
 *       =(4*numLoci - 4*Sum_loci(Sum_allele(sqrt(p1*p2)))/(sum_loci(numAllele) - numLoci))
 */
template<class SH>
double StatHandler<SH>::getChordDist(const age_t& AGE)
{
    age_idx curAge= age_t2idx(AGE);
    
    unsigned int nbStats=0;
    double sumChor=0, value;
    vector<Patch*>::iterator curPop1, curPop2, endPop = get_vSamplePatch().end();
    for(curPop1 = get_vSamplePatch().begin(); curPop1 != endPop; ++curPop1) { // for each pair of pops
        for(curPop2 = curPop1; curPop2 != endPop; ++curPop2) {
            value = getChordDist_perPatchPair(curAge, *curPop1, *curPop2);
            if(value==my_NAN) continue;
            sumChor += value;
            ++nbStats;
        }
    }
    if(!nbStats) return my_NAN;
    return sumChor/nbStats;
}

/* _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/ */

// ******** coancestries analysis ********
// ----------------------------------------------------------------------------------------
// coancestry
// ----------------------------------------------------------------------------------------
template<class SH>double StatHandler<SH>::Coancestry(unsigned char**seq1,
                                                     unsigned char**seq2) {
    unsigned int k, p = 0;
    for (k = 0; k < _nb_locus; ++k) {
        p += !(seq1[k][0] ^ seq2[k][0]) + !(seq1[k][0] ^ seq2[k][1]) + !
        (seq1[k][1] ^ seq2[k][0]) + !(seq1[k][1] ^ seq2[k][1]);
    }
    
    return(double)p / (4. * _nb_locus);
}

// ----------------------------------------------------------------------------------------
// setCoaMatrixTheta
// ----------------------------------------------------------------------------------------
/** computes the pairwise coancestry matrix: diagonal (theta) */
template<class SH> void
StatHandler<SH>::setCoaMatrixTheta(const age_idx & AGE)
{
    // check if the table has already been computed
    if (already_computed(_computed[16], AGE)) return;
    
    double coa;
    unsigned int i, tot_size, wt = 0;
    
    // create the matrix if not present
    if(!_coa_matrix) ARRAY::create_1D<double**>(_coa_matrix, NB_AGE_CLASSES, NULL); // first time
    if(!_coa_matrix[AGE]) 	ARRAY::create_2D(_coa_matrix[AGE], get_current_nbSamplePatch(), get_current_nbSamplePatch(), (double)my_NAN);
    else if(get_last_nbSamplePatch() != get_current_nbSamplePatch() && !already_computed(_computed[17], AGE)) { // the number of sampled patches may change over time
        ARRAY::delete_2D(_coa_matrix[AGE], get_last_nbSamplePatch());
        ARRAY::create_2D(_coa_matrix[AGE], get_current_nbSamplePatch(),
                         get_current_nbSamplePatch(), (double)my_NAN);
    }
    else { // reset
        for (i = 0; i < get_current_nbSamplePatch(); ++i) {
            _coa_matrix[AGE][i][i] = my_NAN; // just reset the diagonal
        }
    }
    
    // first fill the diagonale: within deme coancestry (theta)
    _mean_theta[AGE] = 0;
    vector<TIndividual*>::iterator curInd1, curInd2, endInd;
    vector<Patch*>::iterator curPop, endPop = get_vSamplePatch().end();
    for (i = 0, curPop = get_vSamplePatch().begin(); curPop != endPop; ++i, ++curPop) {
        vector<TIndividual*>&curFem = (*curPop)->get_sampled_inds(FEM, AGE);
        vector<TIndividual*>&curMal = (*curPop)->get_sampled_inds(MAL, AGE);
        coa = 0;
        tot_size = (unsigned int)curFem.size() + (unsigned int)curMal.size();
        if (tot_size) {
            // fem-fem coa
            for (curInd1 = curFem.begin(), endInd = curFem.end(); curInd1 != endInd; ++curInd1) {
                for (curInd2 = curInd1 + 1; curInd2 != endInd; ++curInd2) {
                    coa += Coancestry((unsigned char**)(*curInd1)->getTrait
                                      (_SHLinkedTraitIndex)->get_sequence(),
                                      (unsigned char**)(*curInd2)->getTrait(_SHLinkedTraitIndex)->get_sequence());
                }
            }
            
            // mal-mal coa
            for (curInd1 = curMal.begin(), endInd = curMal.end(); curInd1 != endInd;++curInd1) {
                for (curInd2 = curInd1 + 1; curInd2 != endInd; ++curInd2) {
                    coa += Coancestry((unsigned char**)(*curInd1)->getTrait
                                      (_SHLinkedTraitIndex)->get_sequence(),
                                      (unsigned char**)(*curInd2)->getTrait(_SHLinkedTraitIndex)->get_sequence());
                }
            }
            // fem-mal coa
            coa += get_coancestry(*curPop, FEM, *curPop, MAL, AGE);
            
            if (tot_size == 1) coa = my_NAN;
            else               coa /= tot_size * (tot_size - 1) / 2.0;
        } // end if
        
        _coa_matrix[AGE][i][i] = coa;
        if (coa == my_NAN) _mean_theta[AGE] = my_NAN;
        else               _mean_theta[AGE] += tot_size * coa;
        wt += tot_size;
    } // end for patchNbr
    
    if (wt) {
        if (_mean_theta[AGE] != my_NAN) _mean_theta[AGE] /= wt; // weighted average
    }
    else _mean_theta[AGE] = my_NAN;
}

// ----------------------------------------------------------------------------------------
// setCoaMatrixAlpha
// ----------------------------------------------------------------------------------------
/** computes the pairwise coancestry matrix: upper half part (alpha) */
template<class SH>void StatHandler<SH>::setCoaMatrixAlpha(const age_idx & AGE) {
    // fill the first upper half of the matrix: between deme coancestry (alpha)
    // check if the table has already been computed
    if (already_computed(_computed[17], AGE)) return;
    
    vector<Patch*>::iterator curPop1, curPop2, endPop = get_vSamplePatch().end();
    unsigned int i, j, Fsize1, Msize1, Fsize2, Msize2, tot_size, wt = 0;
    
    double coa;
    
    // create the matrix if not present
    if (!_coa_matrix)
        ARRAY::create_1D<double**>(_coa_matrix, NB_AGE_CLASSES, NULL); // first time
    if (!_coa_matrix[AGE])
        ARRAY::create_2D(_coa_matrix[AGE], get_current_nbSamplePatch(),
                         get_current_nbSamplePatch(), (double)my_NAN);
    else if (get_last_nbSamplePatch() != get_current_nbSamplePatch() && !already_computed
             (_computed[16], AGE)) { // the number of sampled patches may change over time
        ARRAY::delete_2D(_coa_matrix[AGE], get_last_nbSamplePatch());
        ARRAY::create_2D(_coa_matrix[AGE], get_current_nbSamplePatch(),
                         get_current_nbSamplePatch(), (double)my_NAN);
    }
    else { // just reset upper half (without diagonal)
        for (i = 0; i < get_current_nbSamplePatch(); ++i) {
            for (j = i + 1; j < get_current_nbSamplePatch(); ++j) {
                _coa_matrix[AGE][i][j] = my_NAN;
            }
        }
    }
    
    _mean_alpha[AGE] = 0;
    for (i = 0, curPop1 = get_vSamplePatch().begin(); curPop1 != endPop;
         ++i, ++curPop1) { // first patch
        Fsize1 = (*curPop1)->sampleSize(FEM, AGE); // size of patch 1
        Msize1 = (*curPop1)->sampleSize(MAL, AGE); // size of patch 1
        if (!Fsize1 || !Msize1) { // if patch is empty continue with the next patch
            for (j = i + 1; j < get_current_nbSamplePatch(); ++j) {
                // set the corresponding elements to zero
                _coa_matrix[AGE][i][j] = my_NAN;
            }
            continue;
        }
        curPop2 = curPop1;
        for (j = i + 1, ++curPop2; curPop2 != endPop; ++j, ++curPop2) {
            // second patch
            assert(j == (*curPop2)->get_sampleID());
            Fsize2 = (*curPop2)->sampleSize(FEM, AGE); // size of patch 1
            Msize2 = (*curPop2)->sampleSize(MAL, AGE); // size of patch 1
            tot_size = Fsize1 + Msize1 + Fsize2 + Msize2;
            coa = 0;
            if (Fsize2 || Msize2) { // if patch is not empty
                coa += get_coancestry((*curPop1), FEM, (*curPop2), FEM, AGE);
                // fem-fem coa
                coa += get_coancestry((*curPop1), FEM, (*curPop2), MAL, AGE);
                // fem-mal coa
                coa += get_coancestry((*curPop1), MAL, (*curPop2), FEM, AGE);
                // mal-fem coa
                coa += get_coancestry((*curPop1), MAL, (*curPop2), MAL, AGE);
                // mal-mal coa
                coa /= (Fsize1 * Fsize2) + (Fsize1 * Msize2) + (Msize1 * Fsize2) +
                (Msize1 * Msize2);
            } // endif
            
            _coa_matrix[AGE][i][j] = coa;
            _mean_alpha[AGE] += tot_size * coa;
            wt += tot_size;
        } // end for P2
    } // end for P1
    if (wt)
        _mean_alpha[AGE] /= wt; // weighted average
    else
        _mean_alpha[AGE] = my_NAN;
}

// ----------------------------------------------------------------------------------------
// setSexspecific_Theta
// ----------------------------------------------------------------------------------------
template<class SH>double StatHandler<SH>::get_coancestry(Patch * P1,
                                                         const sex_t & SEX1, Patch * P2, const sex_t & SEX2, const age_idx & AGE)
{
    double sum = 0;
    vector<TIndividual*>::iterator cur1, cur2, end1, end2;
    // for each individual of patch 1
    for (cur1 = P1->get_sampled_inds(SEX1, AGE).begin(),
         end1 = P1->get_sampled_inds(SEX1, AGE).end(); cur1 != end1; ++cur1) {
        for (cur2 = P2->get_sampled_inds(SEX2, AGE).begin(),
             // for each individual of patch 2
             end2 = P2->get_sampled_inds(SEX2, AGE).end(); cur2 != end2; ++cur2) {
            sum += Coancestry((unsigned char**)(*cur1)->getTrait(_SHLinkedTraitIndex)
                              ->get_sequence(), (unsigned char**)(*cur2)->getTrait
                              (_SHLinkedTraitIndex)->get_sequence());
        }
    }
    
    return sum;
}

// ----------------------------------------------------------------------------------------
// setSexspecific_Theta
// ----------------------------------------------------------------------------------------
template<class SH>void StatHandler<SH>::setSexspecific_Theta
(const age_idx & AGE) {
    // check if the table has already been computed
    if (already_computed(_computed[18], AGE))
        return;
    
    vector<Patch*>::iterator curPop, endPop;
    vector<TIndividual*>::iterator curInd1, curInd2, endInd1, endInd2;
    unsigned int Fsize, Msize, FFsize, MMsize, FMsize, nbFullPatch = 0;
    
    double mean, grand_mean;
    
    Theta_FF[AGE] = Theta_MM[AGE] = Theta_FM[AGE] = 0;
    _mean_theta[AGE] = 0;
    
    for (curPop = get_vSamplePatch().begin(), endPop = get_vSamplePatch().end();
         curPop != endPop; ++curPop) {
        vector<TIndividual*>&curFem = (*curPop)->get_sampled_inds(FEM, AGE);
        vector<TIndividual*>&curMal = (*curPop)->get_sampled_inds(MAL, AGE);
        Fsize = (unsigned int)curFem.size();
        Msize = (unsigned int)curMal.size();
        
        FFsize = Fsize > 1 ? Fsize * (Fsize - 1) / 2 : 0;
        MMsize = Msize > 1 ? Msize * (Msize - 1) / 2 : 0;
        FMsize = Fsize * Msize;
        
        if (!FMsize)
            continue; // if patch is empty
        ++nbFullPatch;
        
        grand_mean = mean = 0;
        if (FFsize) { // femal - female
            for (curInd1 = curFem.begin(), endInd1 = curFem.end();
                 curInd1 != endInd1; ++curInd1) {
                for (curInd2 = curInd1 + 1; curInd2 != endInd1; ++curInd2) {
                    mean += Coancestry((unsigned char**)(*curInd1)->getTrait
                                       (_SHLinkedTraitIndex)->get_sequence(),
                                       (unsigned char**)(*curInd2)->getTrait(_SHLinkedTraitIndex)
                                       ->get_sequence());
                }
            }
            Theta_FF[AGE] += mean / FFsize;
        }
        grand_mean += mean;
        
        mean = 0;
        if (MMsize) { // male - male
            for (curInd1 = curMal.begin(), endInd1 = curMal.end();
                 curInd1 != endInd1; ++curInd1) {
                for (curInd2 = curInd1 + 1; curInd2 != endInd1; ++curInd2) {
                    mean += Coancestry((unsigned char**)(*curInd1)->getTrait
                                       (_SHLinkedTraitIndex)->get_sequence(),
                                       (unsigned char**)(*curInd2)->getTrait(_SHLinkedTraitIndex)
                                       ->get_sequence());
                }
            }
            Theta_MM[AGE] += mean / MMsize;
        }
        grand_mean += mean;
        
        mean = 0;
        if (FMsize) { // female - male
            for (curInd1 = curFem.begin(), endInd1 = curFem.end();
                 curInd1 != endInd1; ++curInd1) {
                for (curInd2 = curMal.begin(), endInd2 = curMal.end();
                     curInd2 != endInd2; ++curInd2) {
                    mean += Coancestry((unsigned char**)(*curInd1)->getTrait
                                       (_SHLinkedTraitIndex)->get_sequence(),
                                       (unsigned char**)(*curInd2)->getTrait(_SHLinkedTraitIndex)
                                       ->get_sequence());
                }
            }
            Theta_FM[AGE] += mean / FMsize;
        }
        grand_mean += mean;
        _mean_theta[AGE] += grand_mean / (FFsize + MMsize + FMsize);
    }
    
    if (nbFullPatch) {
        _mean_theta[AGE] /= nbFullPatch;
        Theta_FF[AGE] /= nbFullPatch;
        Theta_MM[AGE] /= nbFullPatch;
        Theta_FM[AGE] /= nbFullPatch;
    }
    else
        _mean_theta[AGE] = Theta_FF[AGE] = Theta_MM[AGE] = Theta_FM[AGE] = my_NAN;
}

/* _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/ */

// ******** kinship analysis ********
// ----------------------------------------------------------------------------------------
// setSibsStats
// ----------------------------------------------------------------------------------------
/** sets the sib coancestry and at the same time the kinship as it does not cost a lot of extra effort
 * (the kinship can separately been set by the function setKinship()
 * only computable for individuals  with a previous generation
 */
template<class SH>void StatHandler<SH>::setSibStats(const age_idx & AGE) {
    // check if the table has already been computed
    if (already_computed(_computed[19], AGE))
        return;
    
    unsigned int i, tot_size = 0;
    
    // if it is the first generation stop here
    if (get_current_generation() == 1 || _popPtr->getCoalescence()) {
        for (i = 0; i < 5; ++i) {
            _sib_prop[AGE][i] = _sib_coa[AGE][i] = my_NAN;
        }
        return;
    }
    
    // counters initialization
    for (i = 0; i < 5; ++i) {
        _sib_prop[AGE][i] = _sib_coa[AGE][i] = 0.0;
    }
    
    vector<Patch*>::iterator curPop, endPop = get_vSamplePatch().end();
    vector<TIndividual*>::iterator curInd1, curInd2, endInd1, endInd2;
    for (curPop = get_vSamplePatch().begin(); curPop != endPop; ++curPop) {
        // for each patch
        vector<TIndividual*>&curFem = (*curPop)->get_sampled_inds(FEM, AGE);
        vector<TIndividual*>&curMal = (*curPop)->get_sampled_inds(MAL, AGE);
        tot_size = (unsigned int)curFem.size() + (unsigned int)curMal.size();
        
        // female - female
        for (curInd1 = curFem.begin(), endInd1 = curFem.end(); curInd1 != endInd1;
             ++curInd1) { // for each female
            for (curInd2 = curInd1 + 1; curInd2 != endInd1; ++curInd2) {
                // and another one
                setSibCoa(*curInd1, *curInd2, AGE);
            }
            // selfed offspring counter:
            if ((*curInd1)->getIsSelfed())
                _sib_prop[AGE][4]++;
        }
        
        // male - male
        for (curInd1 = curMal.begin(), endInd1 = curMal.end(); curInd1 != endInd1;
             ++curInd1) { // for each male
            for (curInd2 = curInd1 + 1; curInd2 != endInd1; ++curInd2) {
                // and another one
                setSibCoa(*curInd1, *curInd2, AGE);
            }
            // selfed offspring counter:
            if ((*curInd1)->getIsSelfed())
                _sib_prop[AGE][4]++;
        }
        
        // female - male
        for (curInd1 = curFem.begin(), endInd1 = curFem.end(); curInd1 != endInd1;
             ++curInd1) { // for each female
            for (curInd2 = curMal.begin(), endInd2 = curMal.end();
                 curInd2 != endInd2; ++curInd2) { // for each male
                setSibCoa(*curInd1, *curInd2, AGE);
            }
        }
    }
    
    double tot = _sib_prop[AGE][0] + _sib_prop[AGE][1] + _sib_prop[AGE][2] + _sib_prop[AGE][3];
    
    for (i = 0; i < 4; ++i) {
        _sib_coa[AGE][i] = ((_sib_prop[AGE][i]) ? _sib_coa[AGE][i] / _sib_prop[AGE][i] : my_NAN);
        _sib_prop[AGE][i] = tot ? _sib_prop[AGE][i] / tot : my_NAN;
    }
    _sib_prop[AGE][4] = tot_size ? _sib_prop[AGE][4] / tot_size : my_NAN;
}

// ----------------------------------------------------------------------------------------
// setSibCoa
// ----------------------------------------------------------------------------------------
/** sets the sib coancestry and at the same time the kinship as it does n ot cost a lot of extra effort
 * (the kinship can separately been set by the function setKinship()
 */
template<class SH>void StatHandler<SH>::setSibCoa(TIndividual * I1,
                                                  TIndividual * I2, const age_idx & AGE) {
    double coa = Coancestry((unsigned char**)I1->getTrait(_SHLinkedTraitIndex)
                            ->get_sequence(), (unsigned char**)I2->getTrait(_SHLinkedTraitIndex)
                            ->get_sequence());
    if (I1->getMotherID() == I2->getMotherID()) {
        if (I1->getFatherID() == I2->getFatherID()) {
            _sib_prop[AGE][3]++;
            _sib_coa[AGE][3] += coa;
        } // full sibs
        else {
            _sib_prop[AGE][1]++;
            _sib_coa[AGE][1] += coa;
        } // maternal half sibs
    }
    else {
        if (I1->getFatherID() == I2->getFatherID()) {
            _sib_prop[AGE][2]++;
            _sib_coa[AGE][2] += coa;
        } // paternal half sibs
        else {
            _sib_prop[AGE][0]++;
            _sib_coa[AGE][0] += coa;
        } // non sibs
    }
}

// ----------------------------------------------------------------------------------------
// overall getDprime
// ----------------------------------------------------------------------------------------
/** computes overall locus pairwise linkage disequilibrium measuerment D' for known gametic phase
 * as done in Arlequin:
 * Dprime = sum(n*Dprime_i)/sum(n)
 */
template<class SH> double
StatHandler<SH>::getDprime(const age_idx & AGE,  const unsigned int&l1, const unsigned int&l2)
{
    double Dprime = 0, dp;
    unsigned int size, totSize = 0;
    vector<Patch*>::iterator curPop, endPop = get_vSamplePatch().end();
    for (curPop = get_vSamplePatch().begin(); curPop != endPop; ++curPop) {
        size = (*curPop)->sampleSize(AGE);
        if(!size) continue;
        dp = getDprime_ofPatch(AGE, *curPop, l1, l2);
        if(dp == my_NAN) continue;
        Dprime += dp;
        totSize += size;
    }
    return totSize ? Dprime/totSize : my_NAN;
}

// ----------------------------------------------------------------------------------------
// getDprime
// ----------------------------------------------------------------------------------------
/** computes locus pairwise linkage disequilibrium measurement D' for known gametic phase
 * as done in Arlequin for patch p:
 * if(Dab<0) Dab'=Dab/min[p(a)*p(b); (1-p(a))*(1-p(b))]  Hedrick (1987)
 * else      Dab'=Dab/min[(1-p(a))*p(b); p(a)*(1-p(b))]
 * D' = Sum(p(a)*p(b)*|Dab'|)
 */
template<class SH> double
StatHandler<SH>::getDprime_ofPatch(const age_idx & AGE, Patch* p,
                                   const unsigned int&l1, const unsigned int&l2)
{
    set_alleleFreq(AGE); // allele freuqencies are used
    
    map<unsigned char, double>* freqs = _alleleFreq_local[AGE][p->get_sampleID()];
    if (freqs[l1].size() == 1) return my_NAN; // locus 1 is monomorph
    if (freqs[l2].size() == 1) return my_NAN; // locus 2 is monomorph
    
    map<unsigned char, map<unsigned char, double> > *pAB = get_genotypeFreq(AGE, p, l1, l2);
    map<unsigned char, map<unsigned char, double> > geno1;
    map<unsigned char, double>::iterator geno, a1, a2, end1=freqs[l1].end(), end2=freqs[l2].end();
    
    // for each pariwise allele combination (also if the given genotype does not exist)
    double D, pA, pB, Dprime=0;
    for(a1=freqs[l1].begin(); a1!=end1; ++a1){
        pA=a1->second;
        map<unsigned char, double>& geno1=(*pAB)[a1->first];
        for(a2=freqs[l2].begin(); a2!=end2; ++a2){
            pB=a2->second;
            geno = geno1.find(a2->first);
            if(geno==geno1.end()) D =  - pA*pB;
            else                  D = geno->second - pA*pB;
            
            // compute D' for each genotype (Dab')
            if (D < 0) D /= min(pA*pB, (1-pA)*(1-pB));
            else       D /= min((1-pA)*pB, pA*(1-pB));
            
            // compute D' across all genotypes (D')
            if (D > 0) Dprime += pA * pB * D; // make it absolute
            else       Dprime -= pA * pB * D; // make it absolute
        }
    }
    delete pAB;
    return Dprime;
}

// ----------------------------------------------------------------------------------------
// overall getR2
// ----------------------------------------------------------------------------------------
/** computes overall locus pairwise linkage disequilibrium measurement R2 for known gametic phase
 * as done in Arlequin:
 * R2 = sum(n*R2_i)/sum(n)
 */
template<class SH> double
StatHandler<SH>::getR2(const age_idx & AGE,  const unsigned int&l1, const unsigned int&l2)
{
    double R2 = 0, r2;
    unsigned int size, totSize = 0;
    vector<Patch*>::iterator curPop, endPop = get_vSamplePatch().end();
    for (curPop = get_vSamplePatch().begin(); curPop != endPop; ++curPop) {
        size = (*curPop)->sampleSize(AGE);
        if(!size) continue;
        r2 = getR2_ofPatch(AGE, *curPop, l1, l2);
        if(r2 == my_NAN) continue;
        R2 += size * r2;
        totSize += size;
    }
    return totSize ? R2/totSize : my_NAN;
}

// ----------------------------------------------------------------------------------------
// getR2
// ----------------------------------------------------------------------------------------
/** computes locus pairwise linkage disequilibrium measurement R2 for known gametic phase
 * as done in Arlequin:
 * R2ab=Dab^2/[p(a)*(1-p(a))*p(b)*(1-p(b))]
 * R2 = Sum(p(a)*p(b)*R2ab)  (FSTAT); Zhao (2005)
 *    = Sum(Dab^2/[(1-p(a))*(1-p(b))]    // p(a) and p(b) are reduced by the devision
 */
template<class SH> double
StatHandler<SH>::getR2_ofPatch(const age_idx & AGE, Patch* p, const unsigned int&l1, const unsigned int&l2)
{
    set_alleleFreq(AGE); // allele frequencies are used
    
    map<unsigned char, double>* freqs = _alleleFreq_local[AGE][p->get_sampleID()];
    if (freqs[l1].size() == 1) return my_NAN; // locus 1 is monomorph
    if (freqs[l2].size() == 1) return my_NAN; // locus 2 is monomorph
    
    map<unsigned char, map<unsigned char, double> > *pAB = get_genotypeFreq(AGE, p, l1, l2);
    map<unsigned char, map<unsigned char, double> > geno1;
    map<unsigned char, double>::iterator geno, a1, a2, end1=freqs[l1].end(), end2=freqs[l2].end();
    
    // for each pariwise allele combination (also if the given genotype does not exist)
    double D, pA, pB, R2=0;
    for(a1=freqs[l1].begin(); a1!=end1; ++a1){
        pA=a1->second;
        map<unsigned char, double>& geno1=(*pAB)[a1->first];
        for(a2=freqs[l2].begin(); a2!=end2; ++a2){
            pB=a2->second;
            geno = geno1.find(a2->first);
            if(geno==geno1.end()) D =  - pA*pB;
            else                  D = geno->second - pA*pB;
            R2 += D*D / ((1-pA)*(1-pB));
        }
    }
    delete pAB;
    return R2;
}

// ----------------------------------------------------------------------------------------
// overall getDstar
// ----------------------------------------------------------------------------------------
/** computes overall locus pairwise linkage disequilibrium measuerment D* for known gametic phase
 * Dstar = sum(n*Dstar_i)/sum(n)
 */
template<class SH> double
StatHandler<SH>::getDstar(const age_idx & AGE,  const unsigned int&l1, const unsigned int&l2)
{
    double Dstar = 0, d;
    unsigned int size, totSize = 0;
    vector<Patch*>::iterator curPop, endPop = get_vSamplePatch().end();
    for (curPop = get_vSamplePatch().begin(); curPop != endPop; ++curPop) {
        size = (*curPop)->sampleSize(AGE);
        if(!size) continue;
        d = getDstar_ofPatch(AGE, *curPop, l1, l2);
        if(d==my_NAN) continue;
        Dstar += size * d;
        totSize += size;
    }
    return totSize ? Dstar/totSize : my_NAN;
}

// ----------------------------------------------------------------------------------------
// getDstar
// ----------------------------------------------------------------------------------------
/** computes locus pairwise linkage disequilibrium measurement D* for known gametic phase
 * D* = Sum(Dab^2)/[(1-sum(p(a)^2))(1-sum(p(b)^2))]  (Marujama, 1982)
 *    = Sum(Dab^2)/[Hs(a)*Hs(b)]
 */
template<class SH> double
StatHandler<SH>::getDstar_ofPatch(const age_idx & AGE, Patch* p, const unsigned int&l1, const unsigned int&l2)
{
    set_alleleFreq(AGE); // allele freuqencies are used
    
    map<unsigned char, double>* freqs = _alleleFreq_local[AGE][p->get_sampleID()];
    if (freqs[l1].size() == 1) return my_NAN; // locus 1 is monomorph
    if (freqs[l2].size() == 1) return my_NAN; // locus 2 is monomorph
    
    map<unsigned char, map<unsigned char, double> > *pAB = get_genotypeFreq(AGE, p, l1, l2);
    map<unsigned char, map<unsigned char, double> > geno1;
    map<unsigned char, double>::iterator geno, a1, a2, end1=freqs[l1].end(), end2=freqs[l2].end();
    
    // for each pariwise allele combination (also if the given genotype does not exist)
    double D, pA, pB, sumD=0;
    for(a1=freqs[l1].begin(); a1!=end1; ++a1){
        pA=a1->second;
        map<unsigned char, double>& geno1=(*pAB)[a1->first];
        for(a2=freqs[l2].begin(); a2!=end2; ++a2){
            pB=a2->second;
            geno = geno1.find(a2->first);
            if(geno==geno1.end()) D =  - pA*pB;
            else                  D = geno->second - pA*pB;
            sumD += D*D;
        }
    }
    
    delete pAB;
    
    assert(getHs_ofPatch_andLocus(AGE, p, l1) && getHs_ofPatch_andLocus(AGE, p, l2));
    return sumD / (getHs_ofPatch_andLocus(AGE, p, l1) * getHs_ofPatch_andLocus(AGE, p, l2));
}


// ----------------------------------------------------------------------------------------
// overall getChi2
// ----------------------------------------------------------------------------------------
/** computes overall locus pairwise linkage disequilibrium measurment Chi2 for known gametic phase
 * Chi2 = sum(n*Chi2_i)/sum(n)
 */
template<class SH> double
StatHandler<SH>::getChi2(const age_idx & AGE,  const unsigned int&l1, const unsigned int&l2)
{
    double Chi2 = 0, chi2;
    unsigned int size, totSize = 0;
    vector<Patch*>::iterator curPop, endPop = get_vSamplePatch().end();
    for (curPop = get_vSamplePatch().begin(); curPop != endPop; ++curPop) {
        size = (*curPop)->sampleSize(AGE);
        if(!size) continue;
        chi2 = getChi2_ofPatch(AGE, *curPop, l1, l2);
        if(chi2==my_NAN) continue;
        Chi2 += size * chi2;
        totSize += size;
    }
    return totSize ? Chi2/totSize : my_NAN;
}
// ----------------------------------------------------------------------------------------
// getChi2
// ----------------------------------------------------------------------------------------
/** computes locus pairwise linkage disequilibrium measurement Chi2 for known gametic phase
 * Chi2 = 2*N*Sum(Sum(Dab^2/[p(a)*p(b)]))      Hedrick (1987); Hill (1975)
 * normalized:
 * Chi2' = Chi2/[2N(l-1)]                      Yamazaki (1977)
 *       = Sum(Sum(Dab^2/[p(a)*p(b)]))/(l-1); l: number of alleles of the locus with less alleles
 */
template<class SH> double
StatHandler<SH>::getChi2_ofPatch(const age_idx & AGE, Patch* p, const unsigned int&l1, const unsigned int&l2)
{
    set_alleleFreq(AGE); // allele frequencies are used
    
    map<unsigned char, double>* freqs = _alleleFreq_local[AGE][p->get_sampleID()];
    if (freqs[l1].size() == 1) return my_NAN; // locus 1 is monomorph
    if (freqs[l2].size() == 1) return my_NAN; // locus 2 is monomorph
    
    map<unsigned char, map<unsigned char, double> > *pAB = get_genotypeFreq(AGE, p, l1, l2);
    map<unsigned char, map<unsigned char, double> > geno1;
    map<unsigned char, double>::iterator geno, a1, a2, end1=freqs[l1].end(), end2=freqs[l2].end();
    
    // for each pairwise allele combination (also if the given genotype does not exist)
    double D, pA, pB, Chi2=0;
    for(a1=freqs[l1].begin(); a1!=end1; ++a1){
        pA=a1->second;
        map<unsigned char, double>& geno1=(*pAB)[a1->first];
        for(a2=freqs[l2].begin(); a2!=end2; ++a2){
            pB=a2->second;
            geno = geno1.find(a2->first);
            if(geno==geno1.end()) D =  - pA*pB;
            else                  D = geno->second - pA*pB;
            Chi2 += D*D / (pA*pB);              // the 2N* are reduced
        }
    }
    
    delete pAB;
    
    return Chi2 / (min(freqs[l1].size(), freqs[l2].size()) - 1);
}

/* _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/ */

// ******** stat options ********
// ----------------------------------------------------------------------------------------
// set_stat_coancestry
// ----------------------------------------------------------------------------------------
/** coancestry stat options */
template<class SH> bool
StatHandler<SH>::set_stat_coancestry(string t, string i, string trait,
                                     string token, string end,
                                     age_t AGE, string ageStr)
{
    if (t == "theta")
        return add("Within patch coancestry (" + ageStr + ", " + trait + ")",
                   i + "." + token + ".theta", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getMeanTheta);
    if (t == "alpha")
        return add("Between patch coancestry (" + ageStr + ", " + trait + ")",
                   i + "." + token + ".alpha", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getMeanAlpha);
    if (t == "thetaFF")
        return add("Mean within patch, within females coancestry (" + ageStr +
                   ", " + trait + ")", i + "." + token + ".thetaFF", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getTheta_FF);
    if (t == "thetaMM")
        return add("Mean within patch, within males coancestry (" + ageStr + ", " +
                   trait + ")", i + "." + token + ".thetaMM", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getTheta_MM);
    if (t == "thetaFM")
        return add("Mean within patch, between sexes coancestry (" + ageStr +
                   ", " + trait + ")", i + "." + token + ".thetaFM", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getTheta_FM);
    if (t == "coa.fsib")
        return add("Coancestry of full-sib (" + ageStr + ", " + trait + ")",
                   i + "." + token + ".coa.fsib", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getSibCoaMeans);
    if (t == "coa.phsib")
        return add("Coancestry of paternal half-sib (" + ageStr + ", " + trait +
                   ")", i + "." + token + ".coa.phsib", FLAT, AGE, 1, 0, 0, 0, 0, 0, &StatHandler<SH>::getSibCoaMeans);
    if (t == "coa.mhsib")
        return add("Coancestry of maternal half-sib (" + ageStr + ", " + trait +
                   ")", i + "." + token + ".coa.mhsib", FLAT, AGE, 2, 0, 0, 0, 0, 0, &StatHandler<SH>::getSibCoaMeans);
    if (t == "coa.nsib")
        return add("Coancestry of non-sib (" + ageStr + ", " + trait + ")",
                   i + "." + token + ".coa.nsib", FLAT, AGE, 3, 0, 0, 0, 0, 0, &StatHandler<SH>::getSibCoaMeans);
    if (t == "theta_")
        return add_perPatch(end, "Within coancestry (" + ageStr + ", " + trait + ")",
                            i + "." + token + ".theta", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getCoaTheta);
    if (t == "alpha_")
        return add_pairwisePatch(end, "Mean coancestry (" + ageStr + ", " + trait + ")",
                                 i + "." + token + ".alpha", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getCoaAlpha);
    
    if (t == "coa") {
        add("Within patch coancestry (" + ageStr + ", " + trait + ")",
            i + "." + token + ".theta", FLAT, AGE, 0, 0, 0, 0, 	&StatHandler<SH>::getMeanTheta);
        add("Between patch coancestry (" + ageStr + ", " + trait + ")",
            i + "." + token + ".alpha", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getMeanAlpha);
        add("Mean within patch, within females coancestry (" + ageStr + ", " +
            trait + ")", i + "." + token + ".thetaFF", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getTheta_FF);
        add("Mean within patch, within males coancestry (" + ageStr + ", " +
            trait + ")", i + "." + token + ".thetaMM", FLAT, AGE, 0, 0, 0, 0,	&StatHandler<SH>::getTheta_MM);
        add("Mean within patch, between sexes coancestry (" + ageStr + ", " +
            trait + ")", i + "." + token + ".thetaFM", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getTheta_FM);
        add("Coancestry of full-sib (" + ageStr + ", " + trait + ")",
            i + "." + token + ".coa.fsib", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getSibCoaMeans);
        add("Coancestry of paternal half-sib (" + ageStr + ", " + trait + ")",
            i + "." + token + ".coa.phsib", FLAT, AGE, 1, 0, 0, 0, 0, 0, &StatHandler<SH>::getSibCoaMeans);
        add("Coancestry of maternal half-sib (" + ageStr + ", " + trait + ")",
            i + "." + token + ".coa.mhsib", FLAT, AGE, 2, 0, 0, 0, 0, 0, &StatHandler<SH>::getSibCoaMeans);
        add("Coancestry of non-sib (" + ageStr + ", " + trait + ")",
            i + "." + token + ".coa.nsib", FLAT, AGE, 3, 0, 0, 0, 0, 0, &StatHandler<SH>::getSibCoaMeans);
        return true;
    }
    
    return false;
}

// ----------------------------------------------------------------------------------------
// set_stat_LD
// ----------------------------------------------------------------------------------------
/** linkage disequilibrium stat options */
template<class SH> bool
StatHandler<SH>::set_stat_LD(string t, string i, string trait,
                             string token, string end, age_t AGE, string ageStr)
{
    if (t == "Dprime") return add_pairwiseLocus(end, "Overall linkage disequilibrium D'(" + ageStr + ", " + trait + ")",
                                                i + "." + token + ".Dprime", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getDprime_ij);
    if (t == "Dprime_") return add_pairwiseLocus_perPatch_p(end, "Linkage disequilibirum D' (" + ageStr + ", " + trait + ")",
                                                            i + "." + token + ".Dprime", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getDprime_ij_ofPatch);
    
    if (t == "Dstar") return add_pairwiseLocus(end, "Linkage disequilibrium D*(" + ageStr + ", " + trait + ")",
                                               i + "." + token + ".Dstar", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getDstar_ij);
    if (t == "Dstar_") return add_pairwiseLocus_perPatch(end, "Linkage disequilibrium D*(" + ageStr + ", " + trait + ")",
                                                         i + "." + token + ".Dstar", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getDstar_ij_ofPatch);
    if (t == "R2") return add_pairwiseLocus(end, "Linkage disequilibrium r2(" + ageStr + ", " + trait + ")",
                                            i + "." + token + ".R2", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getR2_ij);
    if (t == "R2_") return add_pairwiseLocus_perPatch(end, "Linkage disequilibrium r2(" + ageStr + ", " + trait + ")",
                                                      i + "." + token + ".R2", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getR2_ij_ofPatch);
    if (t == "Chi2") return add_pairwiseLocus(end, "Linkage disequilibrium Chi2(" + ageStr + ", " + trait + ")",
                                              i + "." + token + ".Chi2", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getChi2_ij);
    if (t == "Chi2_") return add_pairwiseLocus_perPatch(end, "Linkage disequilibrium Chi2(" + ageStr + ", " + trait + ")",
                                                        i + "." + token + ".Chi2", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getChi2_ij_ofPatch);
    
    return false;
}

// ----------------------------------------------------------------------------------------
// set_stat_fstat
// ----------------------------------------------------------------------------------------
/** fstat stat options */
template<class SH> bool
StatHandler<SH>::set_stat_fstat(string t, string i, string trait,
                                string token, string end, age_t AGE, string ageStr)
{
    // number of fixed loci
    if (t == "nbFixLoc")
        return add("Mean number of fixed loci/patch/locus (" + ageStr + ", " + trait + ")",
                   i + "." + token + ".nbFixLoc", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getNbFixedLocus);
    if (t == "nbFixLoc_"){
        if (end == "p_l" || end == "l_p")
            return add_perPatch_andLocus(end, "Number of fixed loci (" + ageStr + ", " + trait + ")",
                                         i + "." + token + ".nbFixLoc", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getNbFixedLocus_ofPatch_andLocus);
        if (end[0] == 'p')
            return add_perPatch(end,
                                "Mean number of fixed loci/locus (" + ageStr + ", " + trait + ")",
                                i + "." + token + ".nbFixLoc", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getNbFixedLocus_ofPatch);
        if (end[0] == 'l')
            return add_perLocus(end,
                                "Mean number of fixed loci/patch (" + ageStr + ", " + trait + ")",
                                i + "." + token + ".nbFixLoc", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getNbFixedLocus_ofLocus);
    }
    if (t == "nbFixLoc.tot")
        return add("Mean total number of fixed loci/locus (" + ageStr + ", " + trait + ")",
                   i + "." + token + ".nbFixLoc.tot", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getNbFixedLocusTot);
    if (t == "nbFixLoc.tot_")
        return add_perLocus(end, "Total number of fixed loci (" + ageStr + ", " + trait + ")",
                            i + "." + token + ".nbFixLoc.tot", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getNbFixedLocusTot_ofLocus);
    
    // number of alleles
    if (t == "nbAll")
        return add("Mean number of alleles/patch/locus (" + ageStr + ", " + trait + ")",
                   i + "." + token + ".nbAll", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getNbAllele);
    if (t == "nbAll_"){
        if (end == "p_l" || end == "l_p")
            return add_perPatch_andLocus(end, "Number of alleles (" + ageStr + ", " + trait + ")",
                                         i + "." + token + ".nbAll", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getNbAllele_ofPatch_andLocus);
        if (end[0] == 'p')
            return add_perPatch(end, "Mean number of alleles/locus (" + ageStr + ", " + trait + ")",
                                i + "." + token + ".nbAll", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getNbAllele_ofPatch);
        if (end[0] == 'l')
            return add_perLocus(end, "Mean number of alleles/patch (" + ageStr + ", " + trait + ")",
                                i + "." + token + ".nbAll", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getNbAllele_ofLocus);
    }
    if (t == "nbAll.tot")
        return add("Mean total number of alleles/locus (" + ageStr + ", " + trait + ")",
                   i + "." + token + ".nbAll.tot", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getNbAlleleTot);
    if (t == "nbAll.tot_")
        return add_perLocus(end, "Total number of alleles (" + ageStr + ", " + trait + ")",
                            i + "." + token + ".nbAll.tot", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getNbAlleleTot_ofLocus);
    
    // observed heterozygosity
    if (t == "ho")
        return add("Ho (Nei & Chesser, 1983)(" + ageStr + ", " + trait + ")",
                   i + "." + token + ".ho", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getHo);
    if (t == "ho_") {
        if (end == "p_l" || end == "l_p")
            return add_perPatch_andLocus(end, "Ho (Nei & Chesser, 1983)(" + ageStr + ", " + trait + ")",
                                         i + "." + token + ".ho", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getHo_ofPatch_andLocus);
        if (end[0] == 'p')
            return add_perPatch(end, "Ho (Nei & Chesser, 1983)(" + ageStr + ", " + trait + ")",
                                i + "." + token + ".ho", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getHo_ofPatch);
        if (end[0] == 'l')
            return add_perLocus(end, "Ho (Nei & Chesser, 1983)(" + ageStr + ", " + trait + ")",
                                i + "." + token + ".ho", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getHo_ofLocus);
    }
    
    // expected heterozygosity
    if (t == "hs")
        return add("Hs (Nei & Chesser, 1983)(" + ageStr + ", " + trait + ")",
                   i + "." + token + ".hs", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getHsnei);
    if (t == "hs_") {
        if (end == "p_l" || end == "l_p")
            return add_perPatch_andLocus(end, "Hs (Nei & Chesser, 1983)(" + ageStr + ", " + trait + ")",
                                         i + "." + token + ".hs", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getHsUnbiased_ofPatch_andLocus);
        if (end[0] == 'p')
            return add_perPatch(end, "Hs (Nei & Chesser, 1983)(" + ageStr + ", " + trait + ")",
                                i + "." + token + ".hs", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getHsUnbiased_ofPatch);
        if (end[0] == 'l')
            return add_perLocus(end, "Hs (Nei & Chesser, 1983)(" + ageStr + ", " + trait + ")",
                                i + "." + token + ".hs", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getHsnei_ofLocus);
    }
    
    // expected heterozygosity (simple)
    if (t == "hs.p2")
        return add("Hs (1-sum(p2))(" + ageStr + ", " + trait + ")",
                   i + "." + token + ".hs.p2", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getHs);
    if (t == "hs.p2_") {
        if (end == "p_l" || end == "l_p")
            return add_perPatch_andLocus(end, "Hs (1-sum(p2))(" + ageStr + ", " + trait + ")",
                                         i + "." + token + ".hs.p2", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getHs_ofPatch_andLocus);
        if (end[0] == 'p')
            return add_perPatch(end, "Hs (1-sum(p2))(" + ageStr + ", " + trait + ")",
                                i + "." + token + ".hs.p2", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getHs_ofPatch);
        if (end[0] == 'l')
            return add_perLocus(end, "Hs (1-sum(p2))(" + ageStr + ", " + trait + ")",
                                i + "." + token + ".hs.p2", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getHs_ofLocus);
    }
    
    // total expected heterozygosity
    if (t == "ht")
        return add("Ht (Nei & Chesser, 1983)(" + ageStr + ", " + trait + ")",
                   i + "." + token + ".ht", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getHtnei);
    if (t == "ht_")
        return add_perLocus(end, "Ht (Nei & Chesser, 1983)(" + ageStr + ", " + trait + ")",
                            i + "." + token + ".ht.p2", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getHtnei_ofLocus);
    
    // total expected heterozygosity (simple)
    if (t == "ht.p2")
        return add("Ht (1-sum(p2))(" + ageStr + ", " + trait + ")",
                   i + "." + token + ".ht.p2", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getHt);
    if (t == "ht.p2_")
        return add_perLocus(end, "Ht (1-sum(p2))(" + ageStr + ", " + trait + ")",
                            i + "." + token + ".ht.p2", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getHt_ofLocus);
    
    // Allelic richness
    if (t == "rs")
        return add("Rs (Mousadik & Petit, 1996)(" + ageStr + ", " + trait + ")",
                   i + "." + token + ".rs", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getRs);
    if (t == "rs_") {
        if (end == "p_l" || end == "l_p")
            return add_perPatch_andLocus(end, "Rs (Mousadik & Petit, 1996)(" + ageStr + ", " + trait + ")",
                                         i + "." + token + ".rs", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getRs_ofPatch_andLocus);
        if (end[0] == 'p')
            return add_perPatch(end, "Rs (Mousadik & Petit, 1996)(" + ageStr + ", " + trait + ")",
                                i + "." + token + ".rs", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getRs_ofPatch);
        if (end[0] == 'l')
            return add_perLocus(end, "Rs (Mousadik & Petit, 1996)(" + ageStr + ", " + trait + ")",
                                i + "." + token + ".rs", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getRs_ofLocus);
    }
    if (t == "rt")
        return add("Rt (Mousadik & Petit, 1996)(" + ageStr + ", " + trait + ")",
                   i + "." + token + ".rt", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getRt);
    if (t == "rt_")
        return add_perLocus(end, "Rt (Mousadik & Petit, 1996)(" + ageStr + ", " + trait + ")",
                            i + "." + token + ".rt", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getRt_ofLocus);
    
    // allelic range
    if (t == "r")
        return add("Allelic range (" + ageStr + ", " + trait + ")",
                   i + "." + token + ".r", FLAT, AGE, 0, 0, 0, 0,
                   &StatHandler<SH>::getRange);
    if (t == "r_") {
        if (end == "p_l" || end == "l_p")
            return add_perPatch_andLocus(end, "Allelic range (" + ageStr + ", " + trait + ")", i + "." + token + ".r",
                                         FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getRange_ofPatch_andLocus);
        if (end[0] == 'p')
            return add_perPatch(end, "Allelic range (" + ageStr + ", " + trait + ")",
                                i + "." + token + ".r", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getRange_ofPatch);
        if (end[0] == 'l')
            return add_perLocus(end, "Allelic range (" + ageStr + ", " + trait + ")",
                                i + "." + token + ".r", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getRange_ofLocus);
    }
    if (t == "r.tot")
        return add("overall allelic range (" + ageStr + ", " + trait + ")",
                   i + "." + token + ".r.tot", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getRangeTot);
    if (t == "r.tot_")
        return add_perLocus(end, "overall allelic range (" + ageStr + ", " + trait + ")",
                            i + "." + token + ".r.tot", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getRangeTot_ofLocus);
    
    // garza williamson
    if (t == "gw")
        return add("GW (Garza & Williamson 2001)(" + ageStr + ", " + trait + ")",
                   i + "." + token + ".gw", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getGW);
    if (t == "gw_") {
        if (end == "p_l" || end == "l_p")
            return add_perPatch_andLocus(end, "GW (Garza & Williamson 2001)(" + ageStr + ", " + trait + ")",
                                         i + "." + token + ".gw", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getGW_ofPatch_andLocus);
        if (end[0] == 'p')
            return add_perPatch(end, "GW (Garza & Williamson 2001)(" + ageStr + ", " + trait + ")",
                                i + "." + token + ".gw", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getGW_ofPatch);
        if (end[0] == 'l')
            return add_perLocus(end, "GW (Garza & Williamson 2001)(" + ageStr + ", " + trait + ")",
                                i + "." + token + ".gw", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getGW_ofLocus);
    }
    if (t == "gw.tot")
        return add("overall GW (Garza & Williamson 2001)(" + ageStr + ", " +
                   trait + ")", i + "." + token + ".gw.tot", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getGWTot);
    if (t == "gw.tot_")
        return add_perLocus(end, "overall GW (Garza & Williamson 2001)(" + ageStr + ", " + trait + ")",
                            i + "." + token + ".gw.tot", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getGWTot_ofLocus);
    
    if (t == "gendiv") {
        add("Mean number of alleles/patch/locus (" + ageStr + ", " + trait + ")",
            i + "." + token + ".nbAll", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getNbAllele);
        add("Mean total number of alleles/locus (" + ageStr + ", " + trait + ")",
            i + "." + token + ".nbAll.tot", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getNbAlleleTot);
        add("Mean number of fixed loci/patch/locus (" + ageStr + ", " + trait + ")",
            i + "." + token + ".nbFixLoc", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getNbFixedLocus);
        add("Mean total number of fixed loci/locus (" + ageStr + ", " + trait + ")",
            i + "." + token + ".nbFixLoc.tot", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getNbFixedLocusTot);
        add("Ho (Nei & Chesser, 1983)(" + ageStr + ", " + trait + ")",
            i + "." + token + ".ho", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getHo);
        add("Hs (Nei & Chesser, 1983)(" + ageStr + ", " + trait + ")",
            i + "." + token + ".hs", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getHsnei);
        add("Ht (Nei & Chesser, 1983)(" + ageStr + ", " + trait + ")",
            i + "." + token + ".ht", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getHtnei);
        return true;
    }
    
    // F-stats following Nei and Chesser
    if (t == "fst")
        return add("Fst (Nei & Chesser, 1983)(" + ageStr + ", " + trait + ")",
                   i + "." + token + ".fst", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getFst);
    if (t == "fis_")
        return add_perLocus(end, "Fis (Nei & Chesser, 1983)(" + ageStr + ", " + trait + ")",
                            i + "." + token + ".fis", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getFis_ofLocus);
    if (t == "fis")
        return add("Fis (Nei & Chesser, 1983)(" + ageStr + ", " + trait + ")",
                   i + "." + token + ".fis", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getFis);
    if (t == "fit")
        return add("Fit (Nei & Chesser, 1983)(" + ageStr + ", " + trait + ")",
                   i + "." + token + ".fit", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getFit);
    if (t == "fit_")
        return add_perLocus(end, "Fit (Nei & Chesser, 1983)(" + ageStr + ", " + trait + ")",
                            i + "." + token + ".fit", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getFit_ofLocus);
    if (t == "fst_") {
        if (end == "pair")
            return add_pairwisePatch(end, "Fst (Nei & Chesser, 1983)(" + ageStr + ", " + trait + ")",
                                     i + "." + token + ".fst", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getFst_ij);
        if (end == "pair_l" || end == "l_pair")
            return add_pairwisePatch_perLocus(end, "Fst (Nei & Chesser, 1983)(" + ageStr + ", " + trait + ")",
                                              i + "." + token + ".fst", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getFst_ij_l_single);
        if (end[0] == 'p')
            return add_pairwisePatch_ij(end, "Fst (Nei & Chesser, 1983)(" + ageStr + ", " + trait + ")",
                                        i + "." + token + ".fst", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getFst_ij_single);
        if (end[0] == 'l')
            return add_perLocus(end, "Fst (Nei & Chesser, 1983)(" + ageStr + ", " + trait + ")",
                                i + "." + token + ".fst", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getFst_ofLocus);
    }
    
    if (t == "fstat") {
        add("Fst (Nei & Chesser, 1983)(" + ageStr + ", " + trait + ")",
            i + "." + token + ".fst", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getFst);
        add("Fis (Nei & Chesser, 1983)(" + ageStr + ", " + trait + ")",
            i + "." + token + ".fis", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getFis);
        add("Fit (Nei & Chesser, 1983)(" + ageStr + ", " + trait + ")",
            i + "." + token + ".fit", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getFit);
        return true;
    }
    
    // F-stats following Weir and Cockerham
    if (t == "fst.wc")
        return add("Fst (Weir & Cockerham, 1984)(" + ageStr + ", " + trait + ")",
                   i + "." + token + ".fst.wc", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getFst_WC);
    if (t == "fis.wc")
        return add("Fis (Weir & Cockerham, 1984)(" + ageStr + ", " + trait + ")",
                   i + "." + token + ".fis.wc", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getFis_WC);
    if (t == "fis.wc_")
        return add_perLocus(end, "Fis (Weir & Cockerham, 1984)(" + ageStr + ", " + trait + ")",
                            i + "." + token + ".fis.wc", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getFis_WC_ofLocus);
    if (t == "fit.wc")
        return add("Fit (Weir & Cockerham, 1984)(" + ageStr + ", " + trait + ")",
                   i + "." + token + ".fit.wc", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getFit_WC);
    if (t == "fit.wc_")
        return add_perLocus(end, "Fit (Weir & Cockerham, 1984)(" + ageStr + ", " + trait + ")",
                            i + "." + token + ".fit.wc", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getFit_WC_ofLocus);
    if (t == "fst.wc_") {
        if (end == "pair")
            return add_pairwisePatch(end, "Fst(Weir & Cockerham, 1984)(" + ageStr + ", " + trait + ")",
                                     i + "." + token + ".fst.wc", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getFst_WC_ij);
        if (end == "pair_l" || end == "l_pair")
            return add_pairwisePatch_perLocus(end, "Fst(Weir & Cockerham, 1984)(" + ageStr + ", " + trait + ")",
                                              i + "." + token + ".fst.wc", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getFst_WC_ij_l_single);
        if (end[0] == 'p')
            return add_pairwisePatch_ij(end, "Fst(Weir & Cockerham, 1984)(" + ageStr + ", " + trait + ")",
                                        i + "." + token + ".fst.wc", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getFst_WC_ij_single);
        if (end[0] == 'l')
            return add_perLocus(end, "Fst (Weir & Cockerham, 1984)(" + ageStr + ", " + trait + ")",
                                i + "." + token + ".fst.wc", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getFst_WC_ofLocus);
    }
    if (t == "fstat.wc") {
        add("Fst (Weir & Cockerham, 1984)(" + ageStr + ", " + trait + ")",
            i + "." + token + ".fst.wc", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getFst_WC);
        add("Fis (Weir & Cockerham, 1984)(" + ageStr + ", " + trait + ")",
            i + "." + token + ".fis.wc", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getFis_WC);
        add("Fit (Weir & Cockerham, 1984)(" + ageStr + ", " + trait + ")",
            i + "." + token + ".fit.wc", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getFit_WC);
        return true;
    }
    
    // fdist stats per locus (Beaumont et al 2002)
    if (t == "het0_")
        return add_perLocus(end, "Het0 (Beaumont & Nichols 1996)(" + ageStr + ", " + trait + ")",
                            i + "." + token + ".het0", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getHet0_ofLocus);
    if (t == "het1_")
        return add_perLocus(end, "Het1 (Beaumont & Nichols 1996)(" + ageStr + ", " + trait + ")",
                            i + "." + token + ".het1", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getHet1_ofLocus);
    if (t == "fst.fdist_")
        return add_perLocus(end, "Fst (Beaumont & Nichols 1996)(" + ageStr + ", " + trait + ")",
                            i + "." + token + ".fst.fdist", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getFst_fdist_ofLocus);
    
    // chrod distance following Cavalli-Sforza and Bodmer (1971)
    if (t == "chord")
        return add("average chord distance (Cavalli-Sforza & Bodmer, 1971)(" + ageStr + ", " + trait + ")",
                   i + "." + token + ".chord", FLAT, AGE, 0, 0, 0, 0, &StatHandler<SH>::getChordDist);
    if (t == "chord_"){
        if (end == "pair")
            return add_pairwisePatch(end, "chord distance (Cavalli-Sforza & Bodmer, 1971)(" + ageStr + ", " + trait + ")",
                                     i + "." + token + ".chord", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getChordDist_ij);
        if (end[0] == 'p')
            return add_pairwisePatch_ij(end, "chord distance (Cavalli-Sforza & Bodmer, 1971)(" + ageStr + ", " + trait + ")",
                                        i + "." + token + ".chord", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getChordDist_ij);
        //		if (end[0] == 'l')
        //			return add_perLocus(end, "average chord distance (Cavalli-Sforza & Bodmer, 1971)(" + ageStr + ", " + trait + ")",
        //			i + "." + token + ".chord", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getChordDist_ofLocus);
    }
    
    return false;
}

// ----------------------------------------------------------------------------------------
// set_stat_all_freq_local
// ----------------------------------------------------------------------------------------
/** local allele freuqencies for each patch stat options "adlt_a.freq_p1_l1_a1"
 */
template<class SH> bool
StatHandler<SH>::set_stat_all_freq_local(string t, string i, string trait,
                                         string token, string end, age_t AGE,
                                         string ageStr)
{
    // group allele frequencies
    if (t == "a.freq") { // for all pops, loci and alleles
        unsigned int a, l, p, nbPatch = get_nbPatch();
        string _i, n;
        
        for (p = 0; p < nbPatch; ++p) {
            if (get_vPatch(p)->get_sampleID() == my_NAN) continue;
            for (l = 0; l < _nb_locus; ++l) {
                for (a = 0; a < _nb_allele[l]; ++a) {
                    set_stat_all_freq_local_text(t,i,trait,ageStr,p,l,a,_i,n);
                    add(n, _i, FLAT, AGE, toID(p,l,a,_nb_patch,_nb_locus,_nb_allele_max),
                        0, 0, 0, 0, 0,&StatHandler<SH>::get_allele_freq_local);
                }
            }
        }
        return true;
    }
    
    if (t == "a.freq_") { // for the specified patch,  and/or locus, and/or allele
        unsigned int p, l, a;
        string _i, n;
        
        // get the patch id
        p = extract_id(end, 'p', get_nbPatch());
        if(p == my_NAN) return false;
        
        // check if the patch is really sampled
        if (get_vPatch(p)->get_sampleID() == my_NAN) {
            warning("Stat '%s' cannot be computed, since the patch %i is not sampled!\n",
                    (t + end).c_str(), p+1);
            return false;
        }
        
        // if just the patch is specified
        if (end.empty()) {
            for (l = 0; l < _nb_locus; ++l) { // for each locus
                for (a = 0; a < _nb_allele[l]; ++a) { // for each allele
                    set_stat_all_freq_local_text(t,i,trait,ageStr,p,l,a,_i,n);
                    add(n, _i, FLAT, AGE, toID(l,a,_nb_locus,_nb_allele_max),
                        0, 0, 0, 0, 0, &StatHandler<SH>::get_allele_freq_local);
                    
                }
            }
            return true;
        }
        
        // get the locus id
        l = extract_id(end, 'l', _nb_locus);
        if(l == my_NAN) return false;
        
        // if the patch and locus are specified
        if (end.empty()) {
            for (a = 0; a < _nb_allele[l]; ++a) { // for each allele
                set_stat_all_freq_local_text(t,i,trait,ageStr,p,l,a,_i,n);
                add(n, _i, FLAT, AGE, toID(l,a,_nb_locus,_nb_allele_max), 0, 0, 0, 0, 0,
                    &StatHandler<SH>::get_allele_freq_local);
            }
            return true;
        }
        
        // get the allele
        a = extract_id(end, 'a', _nb_allele[l]);
        if(a == my_NAN) return false;
        
        set_stat_all_freq_local_text(t,i,trait,ageStr,p,l,a,_i,n);
        add(n, _i, FLAT, AGE, toID(l,a,_nb_locus,_nb_allele_max), 0, 0, 0, 0, 0,
            &StatHandler<SH>::get_allele_freq_local);
        return true;
    }
    
    return false;
}

// ----------------------------------------------------------------------------------------
// set_stat_all_freq_local_text
// ----------------------------------------------------------------------------------------
/** local allele frequencies text generator */
template<class SH> void
StatHandler<SH>::set_stat_all_freq_local_text(string t, string i, string trait,
                                              string ageStr, unsigned int p,
                                              unsigned int l, unsigned int a,
                                              string& stat, string & text)
{
    stat = i + "." + (ageStr=="adult" ? "adlt" : "off") + ".a.freq";
    text = "Local allele frequency of";
    
    // patch
    if (get_nbTotSamplePatch() > 1) {
        stat += "_p" + toStr(p+1, get_nbPatch()); // "_p1"
        text += " population " + toStr(p+1);
    }
    
    // locus
    if (_nb_locus > 1) {
        stat += "_l" + toStr(l + 1, _nb_locus); // "_l1"
        if (get_nbTotSamplePatch() > 1) text += ",";
        text += " locus " + toStr(l + 1);
    }
    
    // allele
    stat += "_a" + toStr(a + 1, _nb_allele_max); // "_a1"
    if (get_nbTotSamplePatch() > 1 || _nb_locus > 1) text += ", and";
    text += " allele " + toStr(a + 1) + " (" + ageStr + ", " + trait + ")";
}

// ----------------------------------------------------------------------------------------
// set_stat_all_freq_global
// ----------------------------------------------------------------------------------------
/** globel allele frequencies for each locus stat options */
template<class SH> bool
StatHandler<SH>::set_stat_all_freq_global(string t, string i, string trait,
                                          string token, string end, age_t AGE,
                                          string ageStr)
{
    // group allele frequencies
    if (t == "a.freq.global") { // for all loci and alleles
        unsigned int a, l;
        string _i, n;
        
        for (l = 0; l < _nb_locus; ++l) {
            for (a = 0; a < _nb_allele[l]; ++a) {
                set_stat_all_freq_global_text(t,i,trait,ageStr,l,a,_i,n);
                add(n, _i, FLAT, AGE, toID(l,a,_nb_locus,_nb_allele_max), 0, 0, 0, 0, 0,
                    &StatHandler<SH>::get_allele_freq_global);
            }
        }
        return true;
    }
    if (t == "a.freq.global_") { // for the specified locus and/or allele
        unsigned int l, a;
        string _i, n;
        
        // get the locus
        l = extract_id(end, 'l', _nb_locus);
        if(l == my_NAN) return false;
        
        // the allele is not specified
        if (end.empty()) { // for all allele of the specified locus
            for (a = 0; a < _nb_allele[l]; ++a) {
                set_stat_all_freq_global_text(t,i,trait,ageStr,l,a,_i,n);
                add(n, _i, FLAT, AGE, toID(l,a,_nb_locus,_nb_allele_max), 0, 0, 0, 0, 0,
                    &StatHandler<SH>::get_allele_freq_global);
            }
            return true;
        }
        
        // get the allele
        a = extract_id(end, 'a', _nb_allele[l]);
        if(a == my_NAN) return false;
        
        set_stat_all_freq_global_text(t,i,trait,ageStr,l,a,_i,n);
        add(n, _i, FLAT, AGE, toID(l,a,_nb_locus,_nb_allele_max), 0, 0, 0, 0, 0,
            &StatHandler<SH>::get_allele_freq_global);
        return true;
    }
    
    return false;
}

// ----------------------------------------------------------------------------------------
// set_stat_all_freq_gloabl_text
// ----------------------------------------------------------------------------------------
/** global allele frequencies text generator */
template<class SH> void
StatHandler<SH>::set_stat_all_freq_global_text(string t, string i, string trait,
                                               string ageStr,
                                               unsigned int l, unsigned int a,
                                               string& stat, string & text)
{
    stat = i + "." + (ageStr=="adult" ? "adlt" : "off") + ".a.freq_global";
    text = "Global allele frequency of";
    
    // locus
    if (_nb_locus > 1) {
        stat += "_l" + toStr(l + 1, _nb_locus); // "_l1"
        if (get_nbTotSamplePatch() > 1) text += ",";
        text += " locus " + toStr(l + 1);
    }
    
    // allele
    stat += "_a" + toStr(a + 1, _nb_allele_max); // "_a1"
    if (get_nbTotSamplePatch() > 1 || _nb_locus > 1) text += ", and";
    text += " allele " + toStr(a + 1) + " (" + ageStr + ", " + trait + ")";
}

// ----------------------------------------------------------------------------------------
// set_stat_locus_freq_local
// ----------------------------------------------------------------------------------------
/** local locus freuqencies for each patch stat options "a.freq_p1_l1_a1"
 */
template<class SH> bool
StatHandler<SH>::set_stat_locus_freq_local(string t, string i, string trait,
                                           string token, string end, age_t AGE,
                                           string ageStr)
{
    // group allele frequencies
    if (t == "l.freq") { // for all pops, loci and alleles
        unsigned int a1, a2, l, p, nbPatch = get_nbPatch();
        string n, _i;
        
        for (p = 0; p < nbPatch; ++p) {
            if (get_vPatch(p)->get_sampleID() == my_NAN) continue;
            for (l = 0; l < _nb_locus; ++l) {
                for (a1 = 0; a1 < _nb_allele[l]; ++a1) {
                    for (a2 = a1; a2 < _nb_allele[l]; ++a2) {
                        set_stat_locus_freq_local_text(t,i,trait,ageStr,p,l,a1,a2,_i,n);
                        add(n, _i, FLAT, AGE, toID(p,l,a1,a2,_nb_patch,_nb_locus,_nb_allele_max,_nb_allele_max),
                            0, 0, 0, 0, 0,  &StatHandler<SH>::get_locus_freq_local);
                    }
                }
            }
        }
        return true;
    }
    
    if (t == "l.freq_") { // for the specified patch,  and/or locus, and/or allele
        unsigned int p, l, a1, a2;
        string _i, n;
        
        // get the patch id
        p = extract_id(end, 'p', get_nbPatch());
        if(p == my_NAN) return false;
        
        // check if the patch is really sampled
        if (get_vPatch(p)->get_sampleID() == my_NAN) {
            warning("Stat '%s' cannot be computed, since the patch %i is not sampled!\n",
                    (t + end).c_str(), p+1);
            return false;
        }
        
        // if just the patch is specified
        if (end.empty()) {
            for (l = 0; l < _nb_locus; ++l) {
                for (a1 = 0; a1 < _nb_allele[l]; ++a1) {
                    for (a2 = a1; a2 < _nb_allele[l]; ++a2) {
                        set_stat_locus_freq_local_text(t,i,trait,ageStr,p,l,a1,a2,_i,n);
                        add(n, _i, FLAT, AGE, toID(p,l,a1,a2,_nb_patch,_nb_locus,_nb_allele_max,_nb_allele_max),
                            0, 0, 0, 0, 0,  &StatHandler<SH>::get_locus_freq_local);
                    }
                }
            }
            return true;
        }
        
        // get the locus id
        l = extract_id(end, 'l', _nb_locus);
        if(l == my_NAN) return false;
        
        // if the patch and locus are specified
        if (end.empty()) {
            for (a1 = 0; a1 < _nb_allele[l]; ++a1) {
                for (a2 = a1; a2 < _nb_allele[l]; ++a2) {
                    set_stat_locus_freq_local_text(t,i,trait,ageStr,p,l,a1,a2,_i,n);
                    add(n, _i, FLAT, AGE, toID(p,l,a1,a2,_nb_patch,_nb_locus,_nb_allele_max,_nb_allele_max),
                        0, 0, 0, 0, 0,  &StatHandler<SH>::get_locus_freq_local);
                }
            }
            return true;
        }
        
        // get the first allele of the genotype
        unsigned int nbDigits = getNbDigits(_nb_allele_max);
        try {
            a1 = strTo<unsigned int>(end.substr(0, nbDigits)) - 1;    // get the first allele
            end.erase(0, nbDigits);
        }
        catch(...) {
            return false;
        }
        if (a1 >= _nb_allele[l]) return false;
        
        // if only a sinlge allele is specified
        if (end.empty()) { // for all locus gentoypes with the specified allele
            unsigned int aa1, aa2;
            for (a2 = 0; a2 < _nb_allele[l]; ++a2) {
                if(a1>a2) {aa1 = a2; aa2 = a1;} // a1 has to be smaller than a2
                else {aa1 = a1; aa2 = a2;}
                set_stat_locus_freq_local_text(t,i,trait,ageStr,p,l,aa1,aa2,_i,n);
                add(n, _i, FLAT, AGE, toID(l,aa1,aa2,_nb_locus,_nb_allele_max,_nb_allele_max), 0, 0, 0, 0, 0,
                    &StatHandler<SH>::get_locus_freq_local);
            }
            return true;
        }
        
        // get the second allele
        if((unsigned int)end.size()!= nbDigits) return false;
        try {
            a2 = strTo<unsigned int>(end.substr(0, nbDigits)) - 1;    // get the second allele
        }
        catch(...) {
            return false;
        }
        if (a2 >= _nb_allele[l]) return false;
        if(a1>a2) swap(a1, a2);     // a1 has to be smaller than a2
        
        set_stat_locus_freq_local_text(t,i,trait,ageStr,p,l,a1,a2,_i,n);
        add(n, _i, FLAT, AGE, toID(l,a1,a2,_nb_locus,_nb_allele_max,_nb_allele_max), 0, 0, 0, 0, 0,
            &StatHandler<SH>::get_locus_freq_local);
        return true;
    }
    
    return false;
}

// ----------------------------------------------------------------------------------------
// set_stat_locus_freq_local_text
// ----------------------------------------------------------------------------------------
/** local locus frequencies text generator */
template<class SH> void
StatHandler<SH>::set_stat_locus_freq_local_text(string t, string i, string trait,
                                                string ageStr, unsigned int p,
                                                unsigned int l, unsigned int a1,
                                                unsigned char a2,
                                                string& stat, string & text)
{
    stat = i + "." + (ageStr=="adult" ? "adlt" : "off") + ".l.freq";
    text = "Local locus frequency of";
    
    // patch
    if (get_nbTotSamplePatch() > 1) {
        stat += "_p" + toStr(p+1, get_nbPatch()); // "_p1"
        text += " population " + toStr(p+1);
    }
    
    // locus
    if (_nb_locus > 1) {
        stat += "_l" + toStr(l + 1, _nb_locus); // "_l1"
        if (get_nbTotSamplePatch() > 1) text += ", ";
        text += " locus " + toStr(l + 1);
    }
    
    // alleles
    assert(a1<=a2);
    stat += "_" + toStr(a1 + 1, _nb_allele_max) + toStr(a2 + 1, _nb_allele_max); // "_0101"
    if (get_nbTotSamplePatch() > 1 || _nb_locus > 1) text += ", and";
    text += " genotype " + toStr(a1 + 1, _nb_allele_max) + toStr(a2 + 1, _nb_allele_max)
    + " (" + ageStr + ", " + trait + ")";
}

// ----------------------------------------------------------------------------------------
// set_stat_locus_freq_global
// ----------------------------------------------------------------------------------------
/** globel locus frequencies for each locus genotype stat options "l.freq.global_l1_0101"
 */
template<class SH> bool
StatHandler<SH>::set_stat_locus_freq_global(string t, string i, string trait,
                                            string token, string end, age_t AGE,
                                            string ageStr)
{
    // group allele frequencies
    if (t == "l.freq.global") { // for all loci and alleles
        unsigned int a1, a2, l;
        string _i, n;
        
        for (l = 0; l < _nb_locus; ++l) {
            for (a1 = 0; a1 < _nb_allele[l]; ++a1) {
                for (a2 = a1; a2 < _nb_allele[l]; ++a2) {
                    set_stat_locus_freq_global_text(t,i,trait,ageStr,l,a1,a2,_i,n);
                    add(n, _i, FLAT, AGE, toID(l,a1,a2,_nb_patch,_nb_allele_max,_nb_allele_max),
                        0, 0, 0, 0, 0, &StatHandler<SH>::get_locus_freq_global);
                }
            }
        }
        return true;
    }
    if (t == "l.freq.global_") { // for the specified locus genotype
        unsigned int l, a1, a2;
        string n, _i;
        
        // get the locus id
        l = extract_id(end, 'l', _nb_locus);
        if(l == my_NAN) return false;
        
        // the genotype is not specified
        if (end.empty()) { // for all allele of the specified locus
            for (a1 = 0; a1 < _nb_allele[l]; ++a1) {
                for (a2 = a1; a2 < _nb_allele[l]; ++a2) {
                    set_stat_locus_freq_global_text(t,i,trait,ageStr,l,a1,a2,_i,n);
                    add(n, _i, FLAT, AGE, toID(l,a1,a2,_nb_patch,_nb_allele_max,_nb_allele_max),
                        0, 0, 0, 0, 0, &StatHandler<SH>::get_locus_freq_global);
                }
            }
            return true;
        }
        
        // get the first allele of the genotype
        unsigned int nbDigits = getNbDigits(_nb_allele_max);
        try {
            a1 = strTo<unsigned int>(end.substr(0, nbDigits)) - 1;    // get the first allele
            end.erase(0, nbDigits);
        }
        catch(...) {
            return false;
        }
        if (a1 >= _nb_allele[l]) return false;
        
        // if only a sinlge allele is specified
        if (end.empty()) { // for all locus gentoypes with the specified allele
            unsigned int aa1, aa2;
            for (a2 = 0; a2 < _nb_allele[l]; ++a2) {
                if(a1>a2) {aa1 = a2; aa2 = a1;} // a1 has to be smaller than a2
                else {aa1 = a1; aa2 = a2;}
                set_stat_locus_freq_global_text(t,i,trait,ageStr,l,aa1,aa2,_i,n);
                add(n, _i, FLAT, AGE, toID(l,aa1,aa2,_nb_locus,_nb_allele_max,_nb_allele_max), 0, 0, 0, 0, 0,
                    &StatHandler<SH>::get_locus_freq_global);
            }
            return true;
        }
        
        // get the second allele
        if((unsigned int)end.size()!= nbDigits) return false;
        try {
            a2 = strTo<unsigned int>(end.substr(0, nbDigits)) - 1;    // get the second allele
        }
        catch(...) {
            return false;
        }
        if (a2 >= _nb_allele[l]) return false;
        if(a1>a2) swap(a1, a2);     // a1 has to be smaller than a2
        
        set_stat_locus_freq_global_text(t,i,trait,ageStr,l,a1,a2,_i,n);
        add(n, _i, FLAT, AGE, toID(l,a1,a2,_nb_locus,_nb_allele_max,_nb_allele_max), 0, 0, 0, 0, 0,
            &StatHandler<SH>::get_locus_freq_global);
        return true;
    }
    
    return false;
}

// ----------------------------------------------------------------------------------------
// set_stat_locus_global_local_text
// ----------------------------------------------------------------------------------------
/** global locus frequencies text generator */
template<class SH> void
StatHandler<SH>::set_stat_locus_freq_global_text(string t, string i, string trait,
                                                 string ageStr,
                                                 unsigned int l, unsigned int a1,
                                                 unsigned char a2,
                                                 string& stat, string & text)
{
    stat = i + "." + (ageStr=="adult" ? "adlt" : "off") + "l.freq_global";
    text = "Global locus frequency of";
    
    // locus
    if (_nb_locus > 1) {
        stat += "_l" + toStr(l + 1, _nb_locus); // "_l1"
        text += " locus " + toStr(l + 1);
    }
    
    // alleles
    assert(a1<=a2);
    stat += "_" + toStr(a1 + 1, _nb_allele_max) + toStr(a2 + 1, _nb_allele_max); // "_0101"
    if (get_nbTotSamplePatch() > 1 || _nb_locus > 1) text += ", and";
    text += " genotype " + toStr(a1 + 1, _nb_allele_max) + toStr(a2 + 1, _nb_allele_max)
    + " (" + ageStr + ", " + trait + ")";
}


// ----------------------------------------------------------------------------------------
// extract_patchID
// ----------------------------------------------------------------------------------------
/** extract the ID from the text and remove the corresponding part of the text
 * returned is the patch id and NaN if not correct
 * text: the input text from which the id is extracted (the id starts at 1, an dthe return starts at 0!!!)
 * first: the first character to extract, respectively the specification of the id
 * max: the max posssible id
 * charOblig: is the character "first" needed or not
 * Note: if max is 1 the return is always 0;
 */
template<class SH>
unsigned int StatHandler<SH>::extract_id(string& text, char first, unsigned max, bool charOblig)
{
    if(text[0]!=first){
        if(charOblig){
            if(get_nbPatch()>1) return my_NAN;              // if more then one patch
            else return 0;                                  // if single patch is present
        }
        else text.erase(0,1);                               // remove "first"
        
    }
    
    string::size_type pos = text.find('_');
    unsigned int id;                                        // id starts at 0
    try {
        if (pos == string::npos){                           // there is no text behind the index
            id = strTo<unsigned int>(text) - 1;
            text.clear();
        }
        else{                                           // there is text behind the index
            id = strTo<unsigned int>(text.substr(0, pos)) - 1;
            text.erase(0,pos+1);                        // remove also the underscore
        }
        if (id >= max) return my_NAN;
    }
    catch(...) {
        return false;
    }
    return id;
}

// ----------------------------------------------------------------------------------------
// get_nbInds
// ----------------------------------------------------------------------------------------
/** get the total number of individuals  of the given age class and sex */
template<class SH>
unsigned int StatHandler<SH>::get_nbInds(const age_idx & AGE, const sex_t & SEX)
{
    unsigned int size = 0;
    vector<Patch*>::iterator curPop = get_vSamplePatch().begin();
    vector<Patch*>::iterator endPop = get_vSamplePatch().end();
    for (; curPop != endPop; ++curPop) {
        size += (*curPop)->size(SEX, AGE);
    }
    
    return size;
}

// ----------------------------------------------------------------------------------------
// get_nbInds
// ----------------------------------------------------------------------------------------
/** get the total number of individuals  of the given age class and sex */
template<class SH>
unsigned int StatHandler<SH>::get_nbInds(const age_t & AGE, const sex_t & SEX)
{
    unsigned int size = 0;
    vector<Patch*>::iterator curPop = get_vSamplePatch().begin();
    vector<Patch*>::iterator endPop = get_vSamplePatch().end();
    for (; curPop != endPop; ++curPop) {
        size += (*curPop)->size(SEX, AGE);
    }
    
    return size;
}

// ----------------------------------------------------------------------------------------
// get_nbInds
// ----------------------------------------------------------------------------------------
/** get the total number of individuals  of the given age class and sex */
template<class SH>
double StatHandler<SH>::get_meanNbInds(const age_idx & AGE, const sex_t & SEX)
{
    unsigned int size = 0, curSize, nbPatch = 0;
    vector<Patch*>::iterator curPop = get_vSamplePatch().begin();
    vector<Patch*>::iterator endPop = get_vSamplePatch().end();
    for (; curPop != endPop; ++curPop) {
        curSize = (*curPop)->size(SEX, AGE);
        if (!curSize)
            continue;
        size += curSize;
        ++nbPatch;
    }
    return nbPatch ? (double)size / nbPatch : my_NAN;
}

// ----------------------------------------------------------------------------------------
// get_nbInds
// ----------------------------------------------------------------------------------------
/** get the total number of individuals  of the given age class and sex */
template<class SH>double StatHandler<SH>::get_meanNbInds(const age_t & AGE,
                                                         const sex_t & SEX) {
    unsigned int size = 0, curSize, nbPatch = 0;
    vector<Patch*>::iterator curPop = get_vSamplePatch().begin();
    vector<Patch*>::iterator endPop = get_vSamplePatch().end();
    for (; curPop != endPop; ++curPop) {
        curSize = (*curPop)->size(SEX, AGE);
        if (!curSize)
            continue;
        size += curSize;
        ++nbPatch;
    }
    return nbPatch ? (double)size / nbPatch : my_NAN;
}

// ----------------------------------------------------------------------------------------
// get_nbInds
// ----------------------------------------------------------------------------------------
/** get the total number of individuals  of the given age class */
template<class SH>double StatHandler<SH>::get_meanNbInds(const age_idx & AGE) {
    unsigned int size = 0, curSize, nbPatch = 0;
    vector<Patch*>::iterator curPop = get_vSamplePatch().begin();
    vector<Patch*>::iterator endPop = get_vSamplePatch().end();
    for (; curPop != endPop; ++curPop) {
        curSize = (*curPop)->size(AGE);
        if (!curSize)
            continue;
        size += curSize;
        ++nbPatch;
    }
    return nbPatch ? (double)size / nbPatch : my_NAN;
}

// ----------------------------------------------------------------------------------------
// get_nbInds
// ----------------------------------------------------------------------------------------
/** get the total number of individuals  of the given age class */
template<class SH>
double StatHandler<SH>::get_meanNbInds(const age_t & AGE)
{
    unsigned int size = 0, curSize, nbPatch = 0;
    vector<Patch*>::iterator curPop = get_vSamplePatch().begin();
    vector<Patch*>::iterator endPop = get_vSamplePatch().end();
    for (; curPop != endPop; ++curPop) {
        curSize = (*curPop)->size(AGE);
        if (!curSize)
            continue;
        size += curSize;
        ++nbPatch;
    }
    return nbPatch ? (double)size / nbPatch : my_NAN;
}

// ----------------------------------------------------------------------------------------
// get_nbSamples
// ----------------------------------------------------------------------------------------
/** get the total number of samples of the given age class and sex */
template<class SH>
unsigned int StatHandler<SH>::get_nbSamples(const age_idx & AGE, const sex_t & SEX)
{
    unsigned int size = 0;
    vector<Patch*>::iterator cur, end;
    for (cur = get_vSamplePatch().begin(), end = get_vSamplePatch().end();
         cur != end; ++cur) {
        size += (*cur)->sampleSize(SEX, AGE);
    }
    
    return size;
}

// ----------------------------------------------------------------------------------------
// get_nbSamples
// ----------------------------------------------------------------------------------------
/** get the total number of samples of the given age class and sex */
template<class SH>
unsigned int StatHandler<SH>::get_nbSamples(const age_t & AGE, const sex_t & SEX)
{
    unsigned int size = 0;
    vector<Patch*>::iterator cur, end;
    for (cur = get_vSamplePatch().begin(), end = get_vSamplePatch().end();
         cur != end; ++cur) {
        size += (*cur)->sampleSize(SEX, AGE);
    }
    
    return size;
}

// ----------------------------------------------------------------------------------------
// get_nbSamples
// ----------------------------------------------------------------------------------------
/** get the mean number of samples of the given age class and sex */
template<class SH>
double StatHandler<SH>::get_meanNbSamples(const age_idx & AGE, const sex_t & SEX)
{
    unsigned int size = 0, curSize, nbPatch = 0;
    vector<Patch*>::iterator cur, end;
    for (cur = get_vSamplePatch().begin(), end = get_vSamplePatch().end();
         cur != end; ++cur) {
        curSize = (*cur)->sampleSize(SEX, AGE);
        if (!curSize)
            continue;
        size += curSize;
        ++nbPatch;
    }
    return nbPatch ? (double)size / nbPatch : (double)my_NAN;
}

// ----------------------------------------------------------------------------------------
// get_nbSamples
// ----------------------------------------------------------------------------------------
/** get the total number of samples of the given age class and sex */
template<class SH>
double StatHandler<SH>::get_meanNbSamples(const age_t & AGE, const sex_t & SEX)
{
    unsigned int size = 0, curSize, nbPatch = 0;
    vector<Patch*>::iterator cur, end;
    for (cur = get_vSamplePatch().begin(), end = get_vSamplePatch().end();
         cur != end; ++cur) {
        curSize = (*cur)->sampleSize(SEX, AGE);
        if (!curSize)
            continue;
        size += curSize;
        ++nbPatch;
    }
    return nbPatch ? (double)size / nbPatch : (double)my_NAN;
}

// ----------------------------------------------------------------------------------------
// get_nbSamples
// ----------------------------------------------------------------------------------------
/** get the mean number of samples of the given age class and sex */
template<class SH>double StatHandler<SH>::get_meanNbSamples
(const age_idx & AGE) {
    unsigned int size = 0, curSize, nbPatch = 0;
    vector<Patch*>::iterator cur, end;
    for (cur = get_vSamplePatch().begin(), end = get_vSamplePatch().end();
         cur != end; ++cur) {
        curSize = (*cur)->sampleSize(AGE);
        if (!curSize)
            continue;
        size += curSize;
        ++nbPatch;
    }
    return nbPatch ? (double)size / nbPatch : (double)my_NAN;
}

// ----------------------------------------------------------------------------------------
// get_nbSamples
// ----------------------------------------------------------------------------------------
/** get the total number of samples of the given age class and sex */
template<class SH>double StatHandler<SH>::get_meanNbSamples(const age_t & AGE) {
    unsigned int size = 0, curSize, nbPatch = 0;
    vector<Patch*>::iterator cur, end;
    for (cur = get_vSamplePatch().begin(), end = get_vSamplePatch().end();
         cur != end; ++cur) {
        curSize = (*cur)->sampleSize(AGE);
        if (!curSize)
            continue;
        size += curSize;
        ++nbPatch;
    }
    return nbPatch ? (double)size / nbPatch : (double)my_NAN;
}

// ----------------------------------------------------------------------------------------
// get_minSampleSize
// ----------------------------------------------------------------------------------------
/** return the minimal sample size */
template<class SH>
unsigned int StatHandler<SH>::get_minSampleSize(const age_idx & AGE)
{
    //unsigned int size = UINT_MAX, curSize;
    unsigned int size = std::numeric_limits<unsigned int>::max(), curSize;
    vector<Patch*>::iterator cur, end;
    for (cur = get_vSamplePatch().begin(), end = get_vSamplePatch().end(); cur != end; ++cur) {
        curSize = (*cur)->sampleSize(AGE);
        if (!curSize) 	continue; // empty
        if (curSize < size) size = curSize;
    }
    return size != std::numeric_limits<unsigned int>::max() ? size : (unsigned int)my_NAN;
}

// ----------------------------------------------------------------------------------------
// get_Rs
// ----------------------------------------------------------------------------------------
// allelic richness
template<class SH>double StatHandler<SH>::getRs(const age_idx & AGE) {
    unsigned int minN = get_minSampleSize(AGE);
    double rs = 0;
    for (unsigned int l = 0; l < _nb_locus; ++l) {
        rs += getRs_ofLocus(l, AGE, minN);
    }
    return rs / _nb_locus;
}

// ----------------------------------------------------------------------------------------
// get_Rs
// ----------------------------------------------------------------------------------------
template<class SH>
double StatHandler<SH>::getRs_ofPatch(Patch * p, const age_idx & AGE, unsigned int minN)
{
    if(!p->sampleSize(AGE)) return my_NAN;
    if (minN == my_NAN) minN = get_minSampleSize(AGE);
    double rs = 0;
    for (unsigned int l = 0; l < _nb_locus; ++l) {
        rs += getRs_ofPatch_andLocus(AGE, p, l, minN);
    }
    return rs / _nb_locus;
}

// ----------------------------------------------------------------------------------------
// get_Rs
// ----------------------------------------------------------------------------------------
template<class SH>
double StatHandler<SH>::getRs_ofLocus(const unsigned int&l,	const age_idx & AGE, unsigned int minN)
{
    if (minN == my_NAN) minN = get_minSampleSize(AGE);
    unsigned int nbFullPatch = 0;
    double val, rs = 0;
    vector<Patch*>::iterator curPop, endPop = get_vSamplePatch().end();
    for (curPop = get_vSamplePatch().begin(); curPop != endPop; ++curPop) {
        val = getRs_ofPatch_andLocus(AGE, *curPop, l, minN);
        if (val == my_NAN) continue;
        rs += val;
        ++nbFullPatch;
    }
    return nbFullPatch ? rs / nbFullPatch : my_NAN;
}

// ----------------------------------------------------------------------------------------
// get_Rs
// ----------------------------------------------------------------------------------------
/** allelic richness with rarefaction following El Mousadik and Petit [1996]
 *  r = C(curSize-nbAllele, minSize) / C(curSize, minSize)
 * minSize is the samllest pop size
 * due to range problems computed as
 *  r = exp(factln(n-i)+factln(n-k)-factln(n)-factln(n-i-k))
 */
template<class SH> double
StatHandler<SH>::getRs_ofPatch_andLocus(const age_idx & AGE, Patch * p,
                                        const unsigned int&l, unsigned int minN)
{
    if (!set_alleleFreq(AGE)) return my_NAN;
    if (p->get_sampleID()==SAMPLED) return my_NAN; // pop not populated
    assert(minN != my_NAN);
    double rs = 0;
    unsigned int nbAll, curSize = 2 * p->sampleSize(AGE); // number of alleles
    minN *= 2; // number of alleles
    double const_part = factln(curSize-minN)-factln(curSize);
    map<unsigned char, double>::iterator pos, end;
    pos = _alleleFreq_local[AGE][p->get_sampleID()][l].begin();
    end = _alleleFreq_local[AGE][p->get_sampleID()][l].end();
    for (; pos != end; ++pos) { // for each allele
        nbAll = my_round(pos->second * curSize);   // current number of this allele
        rs += 1 - exp(factln(curSize-nbAll)-factln(curSize-nbAll-minN)+const_part);
    }
    
    return rs;
}

// ----------------------------------------------------------------------------------------
// get_Rt
// ----------------------------------------------------------------------------------------
template<class SH>double StatHandler<SH>::getRt(const age_idx & AGE) {
    // Ht is NaN if only a single pop is populated
    if (get_nbSamplePatch(AGE) < 2)
        return my_NAN; // 2 or more pops are needed
    unsigned int minN = get_minSampleSize(AGE);
    double rt = 0;
    for (unsigned int l = 0; l < _nb_locus; ++l) {
        rt += getRt_ofLocus(l, AGE, minN);
    }
    return rt / _nb_locus;
}

// ----------------------------------------------------------------------------------------
// get_Rt
// ----------------------------------------------------------------------------------------
/** total allelic richness with rarefaction following El Mousadik and Petit [1996]
 *  r = C(curSize-nbAllele, minSize) / C(curSize, minSize)
 * minSize is the samllest pop size
 * due to range problems computed as
 *  r = exp(factln(n-i)+factln(n-k)-factln(n)-factln(n-i-k))
 */
template<class SH>
double StatHandler<SH>::getRt_ofLocus(const unsigned int&l, const age_idx & AGE, unsigned int minN)
{
    if (!set_alleleFreq(AGE)) return my_NAN;
    if (minN == my_NAN) minN = get_minSampleSize(AGE);
    double rt = 0;
    unsigned int nbAll, curSize = 2 * _popPtr->sampleSize(AGE); // number of alleles
    minN *= 2; // minimal number of alleles
    double const_part = factln(curSize-minN)-factln(curSize);
    map<unsigned char, double>::iterator pos, end;
    pos = _alleleFreq_global[AGE][l].begin();
    end = _alleleFreq_global[AGE][l].end();
    for (; pos != end; ++pos) { // for each allele of this locus
        nbAll = my_round(pos->second * curSize);   // current number of this allele
        rt += 1 - exp(factln(curSize-nbAll)-factln(curSize-nbAll-minN)+const_part);
    }
    return rt;
}

// ----------------------------------------------------------------------------------------
// get_Rs
// ----------------------------------------------------------------------------------------
// allelic richness
template<class SH>
double StatHandler<SH>::getRange(const age_idx & AGE)
{
    double range = 0;
    for (unsigned int l = 0; l < _nb_locus; ++l) {
        range += getRange_ofLocus(l, AGE);
    }
    return range / _nb_locus;
}

// ----------------------------------------------------------------------------------------
// get_Rs
// ----------------------------------------------------------------------------------------
template<class SH>
double StatHandler<SH>::getRange_ofPatch(Patch * p, const age_idx & AGE)
{
    if(!p->sampleSize(AGE)) return my_NAN;
    unsigned int range = 0;
    for (unsigned int l = 0; l < _nb_locus; ++l) {
        range += getRange_ofPatch_andLocus(AGE, p, l);
    }
    return(double)range / _nb_locus;
}

// ----------------------------------------------------------------------------------------
// get_Rs
// ----------------------------------------------------------------------------------------
template<class SH>
double StatHandler<SH>::getRange_ofLocus(const unsigned int&l, const age_idx & AGE)
{
    unsigned int nbFullPatch = 0, range = 0, val;
    vector<Patch*>::iterator curPop, endPop = get_vSamplePatch().end();
    for (curPop = get_vSamplePatch().begin(); curPop != endPop; ++curPop) {
        val = getRange_ofPatch_andLocus(AGE, *curPop, l);
        if (val == my_NAN) continue;
        range += val;
        ++nbFullPatch;
    }
    
    return nbFullPatch ? (double)range / nbFullPatch : my_NAN;
}

// ----------------------------------------------------------------------------------------
// get_Rs
// ----------------------------------------------------------------------------------------
template<class SH>
unsigned int StatHandler<SH>::getRange_ofPatch_andLocus(const age_idx & AGE, Patch * p, const unsigned int&l)
{
    if (p->get_sampleID()==SAMPLED) return my_NAN; // pop not populated
    if (!set_alleleFreq(AGE)) return my_NAN; // get the local allele frequencies if necessary
    return(unsigned int)(_alleleFreq_local[AGE][p->get_sampleID()][l].rbegin()->first       // last
                         - _alleleFreq_local[AGE][p->get_sampleID()][l].begin()->first);     // first
}

// ----------------------------------------------------------------------------------------
// get_Range
// ----------------------------------------------------------------------------------------
template<class SH>
double StatHandler<SH>::getRangeTot(const age_idx & AGE)
{
    unsigned int range = 0;
    for (unsigned int l = 0; l < _nb_locus; ++l) {
        range += getRangeTot_ofLocus(l, AGE);
    }
    return(double)range / _nb_locus;
}

// ----------------------------------------------------------------------------------------
// get_Range
// ----------------------------------------------------------------------------------------
template<class SH>unsigned int StatHandler<SH>::getRangeTot_ofLocus
(const unsigned int&l, const age_idx & AGE) {
    if (!set_alleleFreq(AGE))
        return my_NAN; // get the local allele frequencies if necessary
    return(unsigned int)(_alleleFreq_global[AGE][l].rbegin()->first - _alleleFreq_global[AGE][l].begin()->first);
}

// ----------------------------------------------------------------------------------------
// get_GW
// ----------------------------------------------------------------------------------------
// get gw global
template<class SH>double StatHandler<SH>::getGW(const age_idx & AGE)
{
    unsigned int nbFullPatch=0;
    double gw = 0, curGW;
    vector<Patch*>::iterator curPop, endPop = get_vSamplePatch().end();
    for (curPop = get_vSamplePatch().begin(); curPop != endPop; ++curPop) {
        curGW = getGW_ofPatch(*curPop, AGE);
        if(curGW == my_NAN) continue;
        gw += curGW;
        ++nbFullPatch;
    }
    return nbFullPatch ? gw / nbFullPatch : my_NAN;
}

// ----------------------------------------------------------------------------------------
// get_GW
// ----------------------------------------------------------------------------------------
template<class SH>
double StatHandler<SH>::getGW_ofPatch(Patch * p, const age_idx & AGE)
{
    if(!p->sampleSize(AGE)) return my_NAN;
    
    unsigned int sumAllele = 0, sumRange = 0;
    for(unsigned int l = 0; l < _nb_locus; ++l) {
        sumAllele += _alleleFreq_local[AGE][p->get_sampleID()][l].size();
        sumRange  += getRange_ofPatch_andLocus(AGE, p, l);
    }
    return (double)sumAllele / (_nb_locus+sumRange);      // modified after Excoffier et al (2005)
}

// ----------------------------------------------------------------------------------------
// get_GW
// ----------------------------------------------------------------------------------------
template<class SH>
double StatHandler<SH>::getGW_ofLocus(const unsigned int&l, const age_idx & AGE)
{
    unsigned int nbFullPatch=0;
    double gw = 0, curGW;
    vector<Patch*>::iterator curPop, endPop = get_vSamplePatch().end();
    for (curPop = get_vSamplePatch().begin(); curPop != endPop; ++curPop) {
        curGW = getGW_ofPatch_andLocus(AGE, *curPop, l);
        if(curGW == my_NAN) continue;
        gw += curGW;
        ++nbFullPatch;
    }
    return nbFullPatch ? gw / nbFullPatch : my_NAN;
}

// ----------------------------------------------------------------------------------------
// get_GW
// ----------------------------------------------------------------------------------------
template<class SH>
double StatHandler<SH>::getGW_ofPatch_andLocus(const age_idx & AGE, Patch * p, const unsigned int&l)
{
    if (!set_alleleFreq(AGE)) return my_NAN; // get the local allele frequencies if necessary
    if (p->get_sampleID()==SAMPLED) return my_NAN; // pop not populated
    return(double)_alleleFreq_local[AGE][p->get_sampleID()][l].size() /
    (1 + getRange_ofPatch_andLocus(AGE, p, l));
}

// ----------------------------------------------------------------------------------------
// get_GW
// ----------------------------------------------------------------------------------------
template<class SH>double StatHandler<SH>::getGWTot_ofLocus
(const unsigned int&l, const age_idx & AGE) {
    if (!set_alleleFreq(AGE))
        return my_NAN; // get the local allele frequencies if necessary
    return(double)_alleleFreq_global[AGE][l].size() / (1 + getRangeTot_ofLocus(l,
                                                                               AGE));
}
