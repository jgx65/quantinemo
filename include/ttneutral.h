/** @file ttneutral.h
 *
 *   Copyright (C) 2006 Frederic Guillaume    <guillaum@zoology.ubc.ca>
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



#ifndef ttneutralH
#define ttneutralH

#include "ttrait.h"
#include "filehandler.h"
#include "stathandler.h"

class TTNeutralFH;
class TTNeutralSH;
class TTraitNeutralProto;

/**Microsatellites genome.*/
class TTraitNeutral : public TTrait
{
    
    TTraitNeutralProto* pProto;
    
private:
    
public:
    virtual void set_from_prototype(TTraitProto* T);
    
    TTraitNeutral () : pProto(0){}
    
    TTraitNeutral(const TTraitNeutral& T) : pProto(T.pProto){
        _copyTTraitParameters(T);  // copy the parameters of TTrait
    }
    
    virtual ~TTraitNeutral ();
    
    ///@}
    ///@name Implementations
    ///@{
    virtual TTraitNeutral& operator= (const TTrait& T);
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
    virtual TTraitNeutral*  clone       ( )                      {return new TTraitNeutral(*this);}
    ///@}
};

/**Prototype class for the TTNeutral trait class.**/
class TTraitNeutralProto : public TTraitProto {
    friend class TTraitNeutral; // we allow to access these parameters from TTNeutral directly
protected:
    void ini_paramset();
    
public:
    TTNeutralFH* _writer;
    TTNeutralSH* _stats;
    
    TTraitNeutralProto ( );
    TTraitNeutralProto (int i);
    TTraitNeutralProto(const TTraitNeutralProto& T);
    
    virtual ~TTraitNeutralProto ( );
    
    //implementation of TTraitProto:
    virtual void                     init (TMetapop* pMetapop);
    
    virtual TTraitNeutral*          hatch ();
    
    virtual TTraitNeutralProto*      clone () {return new TTraitNeutralProto(*this);}
    
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
    
    TTNeutralSH (TTraitNeutralProto* TT){
        set(TT);
    }
    
    virtual ~TTNeutralSH ( ){}
    
    virtual bool init ( ) ;
    
    virtual bool setStatRecorders (const string& token);
    virtual string getName() {return "NeutralSH";}
    
};



#endif //TTNEUTRALGENES_H

