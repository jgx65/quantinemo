/** @file filehandler.h
 *
 *   Copyright (C) 2006 Frederic Guillaume    <guillaum@zoology.ubc.ca>
 *   Copyright (C) 2008 Samuel Neuenschwander <samuel.neuenschwander@unil.ch>
 *   Copyright (C) 2018 Frederic Michaud <frederic.a.michaud@gmail.com>

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

#ifndef filehandlerH
#define filehandlerH

#include "ttrait.h"

class TPatch;
class TIndividual;
class TMetapop;

/** Interface to handle file input/output for any SimComponent.
 *  Stores the periodicity parameters and the file path and extension. The replicate file name
 *  is given by the FileServices. A file handler might be set to write output at specific generation
 *  of a specific replicate or at some periodic time during the simulation.
 */
class FileHandler : public Handler {
    
protected:
    string _name;
    
    FileServices* _service;
    
    Param*       _GenerationParam;

    unsigned int       _GenerationOccurrence;  /**Occurence of the stats recording as set by the user (see parameter "stat_log_time").*/
    
    string  _path;
    string  _short_path;
    string  _extension;
    string  _current_filename;         // current file name with the replicate and generation counter
    string _filename;                  // base filename
    string _filename_ped;              // filename for plink .ped file
    string _filename_pheno;            // filename for plink .pheno file
    string _rep_filename;
    unsigned int last_generation;       // when was the last output created
    
    vector<int> _quantiTraitIndexes;   // used for PLINK
    
    void write_Fstat (const age_idx& cur_age, const sex_t& cur_sex, ostream& FILE,
                      TPatch* current_patch, const int& nbPatchDigit, const int& position,
                      bool extened);
    void write_Arlequin (const age_idx& cur_age, const sex_t& cur_sex, ostream& FILE,
                         TPatch* current_patch, const unsigned int& nbIndDigit,
                         const int& position, bool extened);
    void write_Plink_ped (const age_idx& cur_age, const sex_t& cur_sex, ostream& FILE,
                          TPatch* current_patch, const int& nbPatchDigit, const int& position,
                          char sep, bool extended, bool appended);
    void write_Plink_pheno (const age_idx& cur_age, const sex_t& cur_sex, ostream& FILE,
                          TPatch* current_patch, const int& nbPatchDigit, const int& position,
                          char sep);
    
protected:
    string _script;       // a script name which should be executed after the write of a file
    vector<TTraitProto*> _trait;
    vector<int> _TTidx;
    unsigned int _nb_trait;
    
    age_t _age;           // which age should be outputted (ALL; ADULTS, OFFSPRG)?
    int   _sex;           // which sex should be outputet (0: both; 1: only female; 2: only male)?
    unsigned int _format; // 0: none
                          // 1: FSTAT
                          // 2: FSTAT extended
                          // 3: Arlequin
                          // 4: Arlequin extended
                          // 5: PLINK
                          // 5: PLINK extended
    
public:
    
    FileHandler () : _service(0), _GenerationParam(0),
    _GenerationOccurrence(0), _age(0), _sex(0), last_generation(my_NAN) {
    }
    
    FileHandler (Param* gen_occ, string path, string script, int sex, int age,
                 int format, TTraitProto* trait, string name, string ext, TMetapop* ptr)
    : _service(0), _GenerationParam(0), _GenerationOccurrence(0), last_generation(my_NAN)
    {
        set(gen_occ, 0, path, script, sex, age, format, trait, name, ext, NULL, ptr);
    }
    
    FileHandler (unsigned int gen_occ_int, string path, string filename, string script,
                 int sex, int age,
                 int format, TTraitProto* trait, string name, string ext, TMetapop* ptr)
    : _service(0), _GenerationParam(0), _GenerationOccurrence(0), last_generation(my_NAN)
    {
        set(NULL, gen_occ_int, path, filename, script, sex, age, format, trait, name, ext, NULL, ptr);
    }
    
    
    virtual ~FileHandler ( ) { }
    /**Called by notifier during simulation setup, performs file checking.
     *@return false if filename already exists on disk, true otherwise.
     */
    virtual bool  init ( );
    ///@name Accessors
    ///@{
    
    FileServices* get_service () {return _service;}
    
    void set_service (FileServices* srv) {
        if(!srv) return;
        _service = srv;
    }
    
    string&       get_path () {return _path;}
    string&       get_short_path () {return _short_path;}
    
    string get_name()           {return _name;}
    void   set_name(string s)   {_name=s;}
    
    age_t  get_age()            {return _age;}
    
    vector<TTraitProto*>& get_trait(){ return _trait;}
    
    void set_path (string& path);
    void   set_filename(string s);
    
    string&       get_extension ( )  {return _extension;}    
    void          set_extension (string ext) {_extension = ext;}
    
    string&      get_filename ();
    string       getGenerationReplicateFileName ();
    string&      getReplicateFileName ();
    
    virtual int  get_GenerationOccurrence  ()                {return _GenerationOccurrence;}
    virtual void set_GenerationOccurrence  (unsigned int val){_GenerationOccurrence = val;}
    virtual void set_GenerationParam       (Param* p)        {_GenerationParam = p;}
    
    void write_individual_info_to_stream(ostream& FILE, TIndividual* ind, const age_idx& cur_age, const sex_t& cur_sex, char sep=' ');
    
    virtual string get_script ()            {return _script;}
    virtual void set_script (string val)    {_script = val;}
    virtual void set_sex(int i)             {_sex = i;}
    virtual void set_age(int i)             {
        switch(i){
            case 0: _age=ADULTS;    break;
            case 1: _age=OFFSPRG;   break;
            case 2: _age=ALL; 	    break;
        }
    }
    virtual void set_format(int f) {_format = f;}
    unsigned int get_format() {return _format;}
    
    ///@}
    
    /**Sets the handler parameters.
     @param gen_occ generation occurrence
     @param path the file path
     */
    virtual bool set (Param* gen_occ, string path, string filename, string script, int sex, int age,
                      int format, TTraitProto* trait, string name, string ext, FileServices* loader,
                      TMetapop* ptr) {
        return set(gen_occ, 0, path, filename, script, sex, age, format, trait, name, ext, loader, ptr);
    }
    virtual bool set (unsigned int gen_occ_int, string path, string filename, string script,
                      int sex, int age, int format, TTraitProto* trait, string name,
                      string ext, FileServices* loader, TMetapop* ptr) {
        return set(NULL, gen_occ_int, path, filename, script, sex, age, format, trait, name, ext, loader, ptr);
    }
    virtual bool set (Param* gen_occ, unsigned int gen_occ_int, string path, string filename, string script, int sex, int age, int format, TTraitProto* trait, string name, string ext, FileServices* loader, TMetapop* ptr);
    
    virtual void set (TTraitProto* trait) {
        if(trait){
            _trait.push_back(trait);
            _TTidx.push_back(trait->get_absolute_index());
            _nb_trait = (unsigned int)_TTidx.size();
            _popPtr = trait->get_popPtr();
        }
    }
    
    /**Default behaviour of the class, called by Handler::update().**/
    virtual void  FHwrite();
    void (FileHandler::*writeGenotype_func_ptr)(bool);
    virtual void  write_Fstat (bool extened=false);         // FSTAT format
    virtual void  write_Arlequin (bool extened=false);      // Arlequin format
    virtual void  write_Plink_ped (bool extended);          // Plink format
    virtual void  write_Plink_pheno (bool extended, bool appended);        // Plink format
    virtual void  write_Plink_map (bool extended);          // Plink format (computed once per replicate
    virtual void  write_Plink_map(string file, sex_t SEX, bool extended);
    
    virtual void update();
    
    virtual void execute_script(string script, string filename);
    
    virtual unsigned int get_tot_occurrence();
    
    virtual age_t get_age_class(){return _age;}
};
#endif //FILEHANDLER_H


