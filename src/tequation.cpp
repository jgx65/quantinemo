/** @file tequation.cpp
 *
 *   Copyright (C) 2008 Samuel Neuenschwander <samuel.neuenschwander@unil.ch>
 *   Copyright (C) 2018 Frederic Michaud <frederic.michaud@unil.ch>
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


#include "tequation.h"


//-----------------------------------------------------------------------------
// TEquation::read_input
//-----------------------------------------------------------------------------
void
TEquation::read_input(string input, vector<int> qtrait)
{
    string cur_str;
    size_t pos1=0, pos2=input.find_first_of("+-*/");
    while(pos2!=string::npos){
        // get the sign
        switch(input[pos2]){
            case '+': arithmetic_vec.push_back(&TEquation::plus);  break;
            case '-': arithmetic_vec.push_back(&TEquation::minus); break;
            case '*': arithmetic_vec.push_back(&TEquation::mult);  break;
            case '/': arithmetic_vec.push_back(&TEquation::div);   break;
            default: error("Could not read the equation sign!\n");
        }
        
        // get first agument
        cur_str = input.substr(pos1,pos2-pos1);
        rem_forgoing_space(cur_str);
        switch(cur_str[0]){
            case 'Z': {
                values_vec.push_back(&TEquation::getPheno);
                add_index(strTo<int>(cur_str.substr(1)), input, qtrait);
                break;
            }
            case 'G': {
                values_vec.push_back(&TEquation::getGeno);
                add_index(strTo<int>(cur_str.substr(1)), input, qtrait);
                break;
            }
            default:{
                values_vec.push_back(&TEquation::getConst);
                try{
                    index_vec.push_back(strTo<double>(cur_str));
                } catch(...){
                    error("Could not read %s", input.c_str());
                }
                break;
            }
        }
        
        
        // prepare for next search
        pos1 = pos2+1;
        pos2 = input.find_first_of("+-*/", pos1);
    }
    
    // get last argument
    cur_str = input.substr(pos1);
    rem_edge_space(cur_str);
    switch(cur_str[0]){
        case 'Z': {
            values_vec.push_back(&TEquation::getPheno);
            add_index(strTo<int>(cur_str.substr(1)), input, qtrait);
            break;
       }
        case 'G': {
            values_vec.push_back(&TEquation::getGeno);
            add_index(strTo<int>(cur_str.substr(1)), input, qtrait);
            break;
        }
        default:{
            values_vec.push_back(&TEquation::getConst);
            try{
                index_vec.push_back(strTo<double>(cur_str));
            } catch(...){
                error("Could not read %s", input.c_str());
            }
            break;
        }
    }
    end_arithmetic = arithmetic_vec.end();
}

//-----------------------------------------------------------------------------
// TEquation::getValue
//-----------------------------------------------------------------------------
/** return the the value of the equation */
double
TEquation::getValue(TIndividual* ind)
{
    cur_arithmetic = arithmetic_vec.begin();
    cur_value = values_vec.begin();
    cur_index = index_vec.begin();
    
    double val = (this->**cur_value)(ind, *cur_index);  // get the first element
    for(++cur_value, ++cur_index; cur_arithmetic!=end_arithmetic; ++cur_arithmetic, ++cur_index, ++cur_value){
        val = (this->**cur_arithmetic)(val, (this->**cur_value)(ind, *cur_index));
    }
    return val;
}

//-----------------------------------------------------------------------------
// TEquation::add_index
//-----------------------------------------------------------------------------
/** add the trait index to the inde_vector. Index is qtrait specific, change to global */
void
TEquation::add_index(int cur_index, string input, vector<int> qtrait)
{
    try{
        if(cur_index<1) error("Quantitative trait index is out of range!\n");
        if(cur_index>(unsigned int)qtrait.size()) error("Quantitative trait index is out of range!\n");
        index_vec.push_back(qtrait[cur_index-1]);
    } catch(...){
        error("Could not read %s", input.c_str());
    }
}

