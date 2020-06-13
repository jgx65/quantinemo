/** @file tarray.h
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

#ifndef tarrayH
#define tarrayH
//---------------------------------------------------------------------------

#include "functions.h"
#include "types.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <sstream>
using namespace std;


class ARRAY{
private:
    
public:
    //------------------------------------------------------------------------------
    // array functions
    template <typename T>
    static void delete_3D(T*** &array, unsigned int l, unsigned int a){
        try{
            if(array){
                unsigned int i, j;
                for(i=0; i<l; ++i){
                    if(array[i]){
                        for(j=0; j<a; ++j){
                            if(array[i][j]) delete[] array[i][j];
                        }
                        delete[] array[i];
                    }
                }
                delete[] array;
                array = NULL;
            }
        }catch(...) {throw ("Out of memory!\n");}
    }
    
    template <typename T>
    static void delete_3D(T*** &array, unsigned int l, unsigned int* a){
        try{
            if(array){
                unsigned int i, j;
                for(i=0; i<l; ++i, ++a){
                    if(array[i]){
                        for(j=0; j<*a; ++j){
                            if(array[i][j]) delete[] array[i][j];
                        }
                        delete[] array[i];
                    }
                }
                delete[] array;
                array = NULL;
            }
        }catch(...) {throw ("Out of memory!\n");}
    }
    
    template <typename T>
    static void create_3D(T*** &array, unsigned int l, unsigned int a, unsigned int t){    // with initialization
        try{
            if(array) delete_3D<T>(array, l, a);
            unsigned int i, j;
            array = new T** [l];
            for(i=0; i < l; ++i){
                array[i] = new T*[a];
                for(j=0; j < a; ++j){
                    array[i][j] = new T[t];
                }
            }
        }catch(...) {throw ("Out of memory!\n");}
    }
    
    template <typename T>
    static void create_3D(T*** &array, unsigned int l, unsigned int a, unsigned int* t){    // without initialization
        try{
            if(array) delete_3D<T>(array, l, a);
            unsigned int i, j;
            array = new T** [l];
            for(i=0; i < l; ++i){
                array[i] = new T*[a];
                for(j=0; j < a; ++j, ++t){
                    array[i][j] = new T[*t];
                }
            }
        }catch(...) {throw ("Out of memory!\n");}
    }
    
    template <typename T>
    static void create_3D(T*** &array, unsigned int l, unsigned int a, unsigned int t, T init){    // with initialization
        try{
            if(array) delete_3D<T>(array, l, a);
            unsigned int i, j, k;
            array = new T** [l];
            for(i=0; i < l; ++i){
                array[i] = new T*[a];
                for(j=0; j < a; ++j){
                    array[i][j] = new T[t];
                    for(k = 0; k < t; ++k){
                        array[i][j][k] = init;       //init values
                    }
                }
            }
        }catch(...) {throw ("Out of memory!\n");}
    }
    
    template <typename T>
    static void create_3D(T*** &array, unsigned int l, unsigned int a1, unsigned int* a2, T init){    // with initialization
        try{
            if(array) delete_3D<T>(array, l, a1);
            unsigned int i, j, k;
            unsigned int* cur_a2;
            array = new T** [l];
            for(i=0; i < l; ++i){
                array[i] = new T*[a1];
                for(j=0, cur_a2=a2; j < a1; ++j, ++cur_a2){
                    array[i][j] = new T[*cur_a2];
                    for(k = 0; k < *cur_a2; ++k){
                        array[i][j][k] = init;       //init values
                    }
                }
            }
        }catch(...) {throw ("Out of memory!\n");}
    }
    
    // a1 and a2 are in relation to l, thus a1 and a1 have size l!!!
    template <typename T>
    static void create_3D(T*** &array, unsigned int l, unsigned int* a1, unsigned int* a2, T init){    // with initialization
        try{
            if(array) delete_3D<T>(array, l, a1);
            unsigned int i, j, k;
            array = new T** [l];
            for(i=0; i < l; ++i, ++a1, ++a2){
                array[i] = new T*[*a1];
                for(j=0; j < *a1; ++j){
                    array[i][j] = new T[*a2];
                    for(k = 0; k < *a2; ++k){
                        array[i][j][k] = init;       //init values
                    }
                }
            }
        }catch(...) {throw ("Out of memory!\n");}
    }
    
    template <typename T>
    static T*** new_3D(unsigned int l, unsigned int a, unsigned int t){
        T*** array = NULL;
        create_3D(array, l, a, t);
        return array;
    }
    
    template <typename T>
    static T*** new_3D(unsigned int l, unsigned int a, unsigned int t, T init){
        T*** array = NULL;
        create_3D(array, l, a, t, init);
        return array;
    }
    
    template <typename T>
    static void reset_3D(T*** &array, unsigned int l, unsigned int a, unsigned int t, T init){
        assert(array);
        unsigned int i, j, k;
        for(i=0; i < l; ++i){
            assert(array[i]);
            for(j=0; j < a; ++j){
                assert(array[i][j]);
                for(k = 0; k < t; ++k){
                    array[i][j][k] = init;       //init values
                }
            }
        }
    }
    
    // 2Darray functions
    template <typename T>
    static void delete_2D(T** &array, unsigned int l){
        if(array){
            unsigned int i;
            for(i=0; i<l; ++i){
                if(array[i]) delete[] array[i];
            }
            delete[] array;
            array = NULL;
        }
    }
    
    template <typename T>
    static void create_2D(T** &array, unsigned int l, unsigned int a){
        try{
            if(array) delete_2D<T>(array, l);
            unsigned int i;
            array = new T* [l];
            for(i=0; i < l; ++i){
                array[i] = new T[a];
            }
        }catch(...) {throw ("Out of memory!\n");}
    }
    
    template <typename T>                          // if the last arrays differ in length (specified by an array of lengths)
    static void create_2D(T** &array, unsigned int l, unsigned int* a){
        try{
            if(array) delete_2D<T>(array, l);
            unsigned int i;
            array = new T* [l];
            for(i=0; i < l; ++i, ++a){
                array[i] = new T[*a];
            }
        }catch(...) {throw ("Out of memory!\n");}
    }
    
    template <typename T>
    static void create_2D(T** &array, unsigned int l, unsigned int a, T init){   // with initialization
        try{
            if(array) delete_2D<T>(array, l);
            unsigned int i, j;
            array = new T* [l];
            for(i=0; i < l; ++i){
                array[i] = new T[a];
                for(j=0; j < a; ++j){
                    array[i][j] = init;       //init values
                }
            }
        }catch(...) {throw ("Out of memory!\n");}
    }
    
    template <typename T>                 // if the last arrays differ in length (specified by an array of lengths)
    static void create_2D(T** &array, unsigned int l, unsigned int* a, T init){   // with initialization
        try{
            if(array) delete_2D<T>(array, l);
            unsigned int i, j;
            array = new T* [l];
            for(i=0; i < l; ++i, ++a){
                array[i] = new T[*a];
                for(j=0; j < *a; ++j){
                    array[i][j] = init;       //init values
                }
            }
        }catch(...) {throw ("Out of memory!\n");}
    }
    
    template <typename T>
    static T** new_2D(unsigned int l, unsigned int a){
        T** array = NULL;
        create_2D(array, l, a);
        return array;
    }
    template <typename T>
    static T** new_2D(unsigned int l, unsigned int a, T init){
        T** array = NULL;
        create_2D(array, l, a, init);
        return array;
    }
    
    
    template <typename T>
    static void reset_2D(T** &array, unsigned int l, unsigned int a, T init){
        assert(array);
        int i, j;
        for(i=0; i < l; ++i){
            assert(array[i]);
            for(j=0; j < a; ++j){
                array[i][j] = init;       //init values
            }
        }
    }
    
    // 1Darray functions
    template <typename T>
    static void delete_1D(T* &array){
        if(array){
            delete[] array;
            array = NULL;
        }
    }
    
    template <typename T>
    static void create_1D(T* &array, unsigned int l){
        try{
            if(array) delete_1D<T>(array);
            array = new T [l];
        }catch(...) {throw ("Out of memory!\n");}
    }
    
    template <typename T>
    static void create_1D(T* &array, unsigned int l, T init){    // with initialization
        try{
            if(array) delete_1D<T>(array);
            unsigned int i;
            if(l) array = new T [l];
            for(i=0; i < l; ++i){
                array[i] = init;       //init values
            }
        }catch(...) {throw ("Out of memory!\n");}
    }
    
    template <typename T>
    static T* new_1D(unsigned int l){
        T* array = NULL;
        create_1D(array, l);
        return array;
    }
    
    template <typename T>
    static T* new_1D(unsigned int l, T init){
        T* array = NULL;
        create_1D(array, l, init);
        return array;
    }
    
    template <typename T>
    static void reset_1D(T* &array, unsigned int l, T init){
        assert(array);
        unsigned int i;
        for(i=0; i < l; ++i){
            array[i] = init;       //init values
        }
    }
    
    
    // swap two elements in an array
    template <typename T>
    static void swap(T* a, unsigned int i, unsigned int j){
        T temp = a[j];
        a[j] = a[i];
        a[i] = temp;
    }
    
private:
    template <typename T>
    static int increment(const T& i1, const T& i2){
        if( i1 < i2) return 1;
        else return -1;
    }
    template <typename T>
    static int decrement(const T& i1, const T& i2){
        if( i1 < i2) return -1;
        else return 1;
    }
public:
    
    
    //******************************************************************************
    // quicksort for one dimensional array
    //******************************************************************************
private:
    template <typename T>
    static void quicksort( T* a, unsigned int left, unsigned int right, int (*comp)(const T&, const T&)){
        if (left < right) {
            unsigned int lp, rp, i = left + 1, j = left + 1;
            T x = a[left];
            while (j <= right) {
                if((*comp)(a[j], x)<0){
                    //if (a[j] < x) {
                    swap(a, i, j);
                    i++;
                }
                j++;
            }
            swap(a, left, i-1);
            lp = i - 2;
            rp = i;
            
            if (i>1) quicksort(a, left, lp, comp);
            quicksort(a, rp, right, comp);
        }
    }
    
    //------------------------------------------------------------------------------
public:
    template <typename T>
    static void quicksort( T* a, unsigned int n, bool incremental=true){
        if(incremental) quicksort<T>(a, 0, n-1, &decrement<T>);
        else            quicksort<T>(a, 0, n-1, &increment<T>);
    }
    
    //******************************************************************************
    // quicksortIndex for one dimensional array where the index of the sorted array is returned
    //******************************************************************************
    /** get the sorted index
     * it is possible to pass an array for the output (aIndex), otherwise the output
     * array is created and has to be deleted afterwards.
     * The original unsorted array will not be affected by the sort
     */
private:
    template <typename T>
    static void quicksortIndex( T* a, unsigned int left, unsigned int right, unsigned int* index, int (*comp)(const T&, const T&)){
        if (left < right) {
            unsigned int lp, rp, i = left + 1, j = left + 1;
            T x = a[index[left]];
            while (j <= right) {
                if((*comp)(a[index[j]], x)<0){
                    swap(index, i, j);
                    i++;
                }
                j++;
            }
            swap(index, left, i-1);
            lp = i - 2;
            rp = i;
            
            if (i>1) quicksortIndex(a, left, lp, index, comp);
            quicksortIndex(a, rp, right, index, comp);
        }
    }
    
    
    //------------------------------------------------------------------------------
public:
    template <typename T>
    static unsigned int* quicksortIndex( T* a, unsigned int n, unsigned int* index, bool incremental=true){
        // set the indexes (1,2,3,4,...size-1)
        if(!index) index = new unsigned int[n];
        for(unsigned int i = 0; i<n; ++i){
            index[i] = i;
        }
        
        if(incremental) quicksortIndex(a, 0, n-1, index, &decrement<T>);
        else            quicksortIndex(a, 0, n-1, index, &increment<T>);
        return index;
    }
    
    //******************************************************************************
    // quicksortBoth for one dimensional array where the index AND the array are sorted
    //******************************************************************************
    /** sort the array a and change the elements of array index at the same time
     * Both arrays have to have the same size
     */
private:
    template <typename T, typename S>
    static void quicksortBoth( T* a, unsigned int left, unsigned int right, S* index, int (*comp)(const T&, const T&)){
        if (left < right) {
            unsigned int lp, rp, i = left + 1, j = left + 1;
            T x = a[left];
            while (j <= right) {
                if((*comp)(a[j], x)<0){
                    swap(a, i, j);
                    swap(index, i, j);
                    i++;
                }
                j++;
            }
            swap(a, left, i-1);
            swap(index, left, i-1);
            lp = i - 2;
            rp = i;
            
            if (i>1) quicksortBoth(a, left, lp, index, comp);
            quicksortBoth(a, rp, right, index, comp);
        }
    }
    
    
    //------------------------------------------------------------------------------
public:
    template <typename T, typename S>
    static void quicksortBoth( T* a, unsigned int n, S* index, bool incremental=true){
        if(incremental) quicksortBoth<T,S>(a, 0, n-1, index, &decrement<T>);
        else            quicksortBoth<T,S>(a, 0, n-1, index, &increment<T>);
    }
    
    //------------------------------------------------------------------------------
    /** makes the array a cumulative array */
    template <typename T>
    static void cumulative(T* a, unsigned int n){
        for(unsigned int i = 1; i<n; ++i){
            a[i] += a[i-1];
        }
    }
    //------------------------------------------------------------------------------
    /** return a cumulative array to a non cumulative array */
    template <typename T>
    static void undo_cumulative(T* a, unsigned int n){
        for(unsigned int i = n-1; i>0; --i){
            a[i] -= a[i-1];
        }
    }
    //------------------------------------------------------------------------------
    /** the smallest values become the biggest ones
     * a value of zero is approximated by 1e-9
     */
    template <typename T>
    static void reverse_cumulative(T* a, unsigned int n){
        a[0] = a[0] ? (1.0/a[0]) : 1e9; // first element
        for(unsigned int i = 1; i<n; ++i){
            a[i] = a[i] ? (1.0/a[i] + a[i-1]) : (1e9 + a[i-1]);
        }
    }
    
    //------------------------------------------------------------------------------
    /** return a cumulative array to a non cumulative array
     * it is not possible to return to the absolut values, but the relation between the values is maintained
     */
    template <typename T>
    static void undo_reverse_cumulative(T* a, unsigned int n){
        for(unsigned int i = n-1; i>0; --i){
            a[i] = 1.0/(a[i]-a[i-1]);
        }
        a[0] = 1.0/a[0]; // last element
    }
    
    //******************************************************************************
    // compute mean, variance, standard variance
    //******************************************************************************
    /** compute the sum of an array */
    template <typename T>
    static T sum(T* array, unsigned int size){
        if(!size) return my_NAN;
        
        T sum=0;
        T* end = &array[size];          // get the end
        for(; array != end; ++array){
            if(*array == my_NAN) continue;   // continue if the value is missing
            sum += *array;
        }
        return  sum;
    }
    
    //------------------------------------------------------------------------------
    /** compute the sum of an array */
    template <typename T>
    static T sum2D(T** array, unsigned int size, unsigned int pos){
        if(!size) return my_NAN;
        
        T sum=0;
        T elem;
        for(unsigned int i=0; i<size; ++i){
            elem = array[i][pos];
            if(elem == my_NAN) continue;   // continue if the value is missing
            sum += elem
            ;
        }
        return  sum;
    }
    
    //------------------------------------------------------------------------------
    /** compute the mean of an array */
    template <typename T>
    static double mean(T* array, unsigned int size){
        if(!size) return my_NAN;
        
        double sum=0;
        unsigned int nb=0;
        T* end = array + size;          // get the end
        
        for(; array != end; ++array){   // for each element
            if(*array == my_NAN) continue;// if it is not a missing value
            sum += *array;
            ++nb;
        }
        return  nb ? sum/nb : my_NAN;
    }
    
    //------------------------------------------------------------------------------
    /** compute the mean of an array */
    template <typename T>
    static double mean(vector<T> vec)
    {
        if(vec.empty()) return my_NAN;
        
        double sum=0;
        
        typename vector< T >::const_iterator cur, end;
        for(cur=vec.begin(), end=vec.end(); cur!=end; ++cur){
            sum += *cur;
        }
        return  sum/vec.size();
    }
    
    
    //------------------------------------------------------------------------------
    /** compute the mean of a column of a 2D matrix */
    template <typename T>
    static double mean2D(T** array, unsigned int size, unsigned int pos){
        if(!size) return my_NAN;
        
        double sum=0;
        T elem;
        unsigned int nb=0;
        
        for(unsigned int i=0; i<size; ++i){  // for each element
            elem = array[i][pos];
            if(elem == my_NAN) continue;// if it is not a missing value
            sum += elem;
            ++nb;
        }
        return  nb ? sum/nb : my_NAN;
    }
    
    //------------------------------------------------------------------------------
    /** compute the median of an array */
    template <typename T>
    static double median(T* array, unsigned int size){
        if(!size) return my_NAN;
        
        quicksort(array, size);
        if(size % 2) return array[(size-1)/2];  	// odd number
        unsigned int idx = size/2;                // even number
        return (array[idx] + array[idx-1])/2;
    }
    
    //------------------------------------------------------------------------------
    /** compute the mean of a column of a 2D matrix */
    template <typename T>
    static double median2D(T** array, unsigned int size, unsigned int pos){
        if(!size) return my_NAN;
        T* new_array = new T[size];
        for(unsigned int i=0; i<size; ++i){  // for each element
            new_array[i] = array[i][pos];
        }
        
        double val =  median(new_array, size);
        delete[] new_array;
        return val;
    }
    
    //------------------------------------------------------------------------------
    /** compute the sum of an array */
    template <typename T>
    static string print_2D(T** array, unsigned int size1, unsigned int size2){
        ostringstream out;
        for(unsigned int j, i=0; i<size1; ++i){
            if(j) out << "\n";
            for(j=0; j<size2; ++j){
                if(j) out << "\t";
                out << array[i][j];
            }
        }
        return out.str();
    }
    
    
    //------------------------------------------------------------------------------
    /** compute the unbiased variance of an array */
    template <typename T>
    static double varUnbiased(T* array, unsigned  int size, double mu=my_NAN){
        if (size == 1) return my_NAN;
        double bias_var = var(array, size, mu);
        double unbias_var = bias_var*size/(size-1);
        return unbias_var;
    }
    
    //------------------------------------------------------------------------------
    /** compute the biased variance of an vector */
    template <typename T>
    static double varUnbiased(vector<T>& vec, double mu=my_NAN){
        unsigned int size = vec.size();
        if ( size == 1) return my_NAN;
        double bias_var = var(vec, mu);
        double unbias_var = bias_var*(size)/(size - 1);
        return unbias_var;
    }
    
    //------------------------------------------------------------------------------
    /** compute the biased variance of a vector */
    template <typename T>
    static double var(vector<T>& vec, double mu=my_NAN){
        double bias_var = var(&vec[0], vec.size(), mu);
        return bias_var;
    }
    
    
    //------------------------------------------------------------------------------
    /** compute the biased variance of an array */
    template <typename T>
    static double var(T* array, unsigned  int size, double mu=my_NAN){
        if(mu == my_NAN) mu = mean(array, size);  // if the mean has not be computed
        if(mu == my_NAN) return my_NAN;                    // the array has no valid values
        
        double diff, sum=0;
        unsigned int nb=0;
        T* end = array + size;  // get the end
        
        for(; array != end; ++array) {                // for each element
            if(*array == my_NAN) continue;              // if it is not a valid element
            diff = (*array)-mu;
            sum  += diff*diff;
            ++nb;
        }
        
        sum /=nb;
        return (abs(sum) > 1e-30) ? sum : 0;                                 // no control as it was alredy checked
    }
    

    
    //------------------------------------------------------------------------------
    /** compute the biased variance of a column of a 2D matrix */
    template <typename T>
    static double var2D(T** array, unsigned int size, unsigned int pos, double mu=my_NAN){
        if(mu == my_NAN) mu = mean2D(array, size, pos);  // if the mean has not be computed
        if(mu == my_NAN) return my_NAN;             // the array has no valid values
        
        double diff, sum=0;
        unsigned int nb=0;
        T elem;
        
        for(unsigned int i=0; i<size; ++i) {                // for each element
            elem = array[i][pos];
            if(elem == my_NAN) continue;              // if it is not a valid element
            diff = elem-mu;
            sum  += diff*diff;
            ++nb;
        }
        
        sum /=nb;
        return (abs(sum) > 1e-30) ? sum : 0;                                 // no control as it was alredy checked
    }
    
    //------------------------------------------------------------------------------
    /** compute the standard deviation (sd) of an array */
    template <typename T>
    static double sd(T* array, unsigned int size){
        
        // special cases
        if(size<=1) return my_NAN;
        
        return sqrt(var(array, size));
    }
    
    //------------------------------------------------------------------------------
    /** transforms the array that the sum of the array is 1
     * returns true, if this was already the case, and false if the array had to be adjusted
     */
    template <typename T>
    static bool make_frequency(T* array, unsigned int size){
        
        double Sum = (double) sum(array, size);
        if(abs(1-Sum) < 1e-4) return true;
        
        // if all values are zero
        if(!Sum) Sum = 1.0/size;
        
        T* end = array + size;          	// get the end
        
        for(; array != end; ++array) {
            (*array) = (*array)/Sum;
            if(*array<1e-13) *array = 0;		// due to accuracy problems
        }
        return false;
    }
    
    //------------------------------------------------------------------------------
    /** copy an array (with its contents) from the source to the destination
     * both arrays have to have the same size!!!
     */
    template <typename T>
    static void copy(T* destination, T* source, unsigned int size){
        std::copy(source, source+size, destination);
        //memcpy(destination, source, size*sizeof(source));
        //for(unsigned int i=0; i<size; ++i, ++destination, ++source){
        //	*destination = *source;
        //}
    }
    
    //------------------------------------------------------------------------------
    /** transpose a 2D array: array[i][j] => array[j][i]
     * a new array is created which has to be deleted
     * the old array still exists
     */
    template <typename T>
    static T** transpose(T** in, unsigned int s1, unsigned int s2)
    {
        T** array = new T*[s2];
        unsigned int i, j;
        for(i=0; i<s2; ++i){
            array[i] = new T[s1];
            for(j=0; j<s1; ++j){
                array[i][j] = in[j][i];
            }
        }
        return array;
    }
    
    //------------------------------------------------------------------------------
    /** quicksearch within a sorted array
     * if the element is not found the element right to it (larger) will be returned
     * if the element is outside the range of the array NaN is returned
     */
    template <typename T>
    static unsigned int searchBest(T* array, unsigned int size, const T& elem)
    {
        if(!size)                return my_NAN;		// empty array
        if(elem < array[0])      return my_NAN;   // elem too small
        if(elem > array[size-1]) return my_NAN;   // elem too large
        if(elem == array[0])     return 0;         // first element (is not checked below...
        
        double test;
        unsigned int start, end, middle;
        
        start = 0;
        end = size - 1;
        while (end - start > 1) {
            middle = (start + end) / 2;
            test = array[middle];
            if (elem < test) 		    end = middle;
            else if (elem  > test)	start = middle;
            else					          return middle;
        }
        return end;
    }
    
    //------------------------------------------------------------------------------
    /** quicksearch within a sorted array
     * if the element is not found NaN is returned
     */
    template <typename T>
    static unsigned int search(T* array, unsigned int size, const T& elem)
    {
        if(!size)                return my_NAN;		// empty array
        if(elem < array[0])      return my_NAN;   // elem too small
        if(elem > array[size-1]) return my_NAN;   // elem too large
        
        double test;
        unsigned int start, end, middle;
        
        start = 0;
        end = size - 1;
        while (end - start > 1) {
            middle = (start + end) / 2;
            test = array[middle];
            if (elem < test) 		    end = middle;
            else if (elem  > test)	start = middle;
            else					          return middle;
        }
        return array[end]==elem ? end : my_NAN;
    }
    
    
    //------------------------------------------------------------------------------
    /** quicksearch within a sorted array
     * if the element is not found return the next smaller one and
     * if there is no smaller one return NaN
     */
    template <typename T>
    static unsigned int search_smallerEqual(T* array, unsigned int size, const T& elem)
    {
        if(!size)                return my_NAN;		// empty array
        if(elem < array[0])      return my_NAN;   // elem too small
        if(elem > array[size-1]) return size-1;   // elem too large
        
        T test;
        unsigned int start, end, middle;
        
        start = 0;
        end = size - 1;
        while (end - start > 1) {
            middle = (start + end) / 2;
            test = array[middle];
            if (elem < test)        end = middle;
            else if (elem  > test)	start = middle;
            else                    return middle;
        }
        
        return array[end]==elem ? end : start;
    }
    
    
    
    //------------------------------------------------------------------------------
    /** search of an element within an UNsorted array
     * if the element is not found NaN is returned
     */
    template <typename T>
    static unsigned int searchUnsorted(T* array, unsigned int size, const T& elem)
    {
        
        unsigned int i=0,j=size-1;
        unsigned int k=(i+j)/2;
        unsigned int l=k+1;
        unsigned int index=my_NAN;
        
        while((i<=k)&&(l<=j)){
            if(array[i]==elem) index=i;
            if(array[k]==elem) index=k;
            if(array[l]==elem) index=l;
            if(array[j]==elem) index=j;
            i++;
            k--;
            l++;
            j--;
        }
        
        if(size%2!=0 && array[i]==elem) index=i;
        return index;
    }
    
    
    //------------------------------------------------------------------------------
    /** copies the elements from a vector to an array
     * if the array is generated and of the correct size, the array is overtaken, otherwise resized
     */
    template <typename T>
    static void vector2array(vector<T> vec, T* & array, unsigned int& size)
    {        
        if(size != vec.size() || !array){
            if(array) delete[] array;
            size = (unsigned int)vec.size();
            array = new T[size];
        }
        std::copy(vec.begin(), vec.end(), array);
        /*
         for(unsigned int i=0; i<size; ++i){
         array[i] = vec[i];
         }
         */
    }
    
    //------------------------------------------------------------------------------
    /** copies the elements from an array to a vector	*/
    template <typename T>
    static void array2vector(T* & array, vector<T> vec, unsigned int& size)
    {
        vec.resize(size);
        std::copy(array, array + size, vec.begin() );
    }
    
};

#endif
