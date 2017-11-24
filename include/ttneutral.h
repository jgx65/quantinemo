/** @file ttneutral.h
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



#ifndef ttneutralH
#define ttneutralH

#include "ttrait.h"
#include "filehandler.h"
#include "stathandler.h"

class TTNeutralFH;
class TTNeutralSH;
class TTNeutralProto;

/**Microsatellites genome.*/
class TTNeutral : public TTrait
{
    
    TTNeutralProto* pProto;
    
private:
    
public:
    virtual void set_from_prototype(TTraitProto* T);
    
    TTNeutral () : pProto(0){}
    
    TTNeutral(const TTNeutral& T) : pProto(T.pProto){
        _copyTTraitParameters(T);  // copy the parameters of TTrait
    }
    
    virtual ~TTNeutral ();
    
    ///@}
    ///@name Implementations
    ///@{
    virtual TTNeutral& operator= (const TTrait& T);
    virtual bool    operator== (const TTrait& T);
    virtual bool    operator!= (const TTrait& T);
    virtual void    reset                ( ){}
    virtual void*   set_trait            (void* value)           {return NULL;}
    inline virtual void**  get_sequence  ( )  const              {return (void**)sequence;}
    inline virtual void*   get_allele    (const unsigned int& loc, const unsigned int& all)  const;
    virtual void    set_sequence         (void** seq)            {reset();sequence = (unsigned char**)seq;}
    virtual void    set_value            ( )                     { }
    virtual void    set_value            (double value)          {return;}
    virtual double  get_value            ( )					           {return my_NAN;}
    virtual void    show_up              ( );
    virtual TTNeutral*  clone       ( )                      {return new TTNeutral(*this);}
    ///@}
};

/**Prototype class for the TTNeutral trait class.**/
class TTNeutralProto : public TTraitProto {
    friend class TTNeutral; // we allow to access these parameters from TTNeutral directly
protected:
    void ini_paramset();
    
public:
    TTNeutralFH* _writer;
    TTNeutralSH* _stats;
    
    TTNeutralProto ( );
    TTNeutralProto (int i);
    TTNeutralProto(const TTNeutralProto& T);
    
    virtual ~TTNeutralProto ( );
    
    //implementation of TTraitProto:
    virtual void                     init (TMetapop* pMetapop);
    
    virtual TTNeutral*          hatch ();
    
    virtual TTNeutralProto*      clone () {return new TTNeutralProto(*this);}
    
    //implementation of SimComponent:
    virtual void loadFileServices ( FileServices* loader );
    
    virtual void loadStatServices ( StatServices* loader );
    
    string get_info();
    
    void temporal_change(const unsigned int& gen);
};

////////////////////////////////////////////////////////////////////////////////
/**A file handler to save the neutral markers genotypes in the FSTAT format*/
class TTNeutralFH: public FileHandler {
    
public:
    
    TTNeutralFH (){}
    
    virtual ~TTNeutralFH ( ) { }
    
};

////////////////////////////////////////////////////////////////////////////////
/**The stat handler for neutral markers. */
class TTNeutralSH: public StatHandler<TTNeutralSH> {
    
public:
    
    TTNeutralSH (TTNeutralProto* TT){
        set(TT);
    }
    
    virtual ~TTNeutralSH ( ){}
    
    virtual bool init ( ) ;
    
    virtual bool setStatRecorders (const string& token);
    virtual string getName() {return "NeutralSH";}
    
};



#endif //TTNEUTRALGENES_H

