/** @file lce_disperse.h
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

#ifndef lce_disperseH
#define lce_disperseH

#include "lifecycleevent.h"


/**The base class of the dispersal LCEs, all events move offspring to the post-dispersal patch containers.
 * Stores the dispersal matrices and dispersal model parameters and interface.
 **/
class LCE_Disperse: public LCE
{
protected:
    int    _disp_model;             // 0: island; 1: propagule island, 2: 1D SS; 3: 2D SS
    int    _border_model;           // 0: torus; 1: reflecting; 2: absorbing
    int    _lattice_range;          // 0: 4 neighbours; 1: 8 neighbours
    unsigned int _rel_abs_disp_rate[2];   // 0: only disp rates; 1: only absolute numbers;
    double _disp_propagule_prob;
    unsigned int _nb_patch;
    unsigned int _x_size, _y_size;  // dimensions of the 2D lattice 0: nbRows(x); 1: nbColumns(y)
    
    // used when the migration rate differs depending on the popualtion density (for each sex)
    double* _disp_factor[2];        //    0: min           (default: 0)
    //    1: max           (default: 1)
    //    2: max growth    (default: 0)
    //    3: growth rate   (default: 1e9)
    //    4: symmetry      (default: 1)
    
    /**The sex-specific dispersal matrices, [0] for males, [1] for females, might be used as connectivity matrix as well*/
    TMatrix* _dispMatrix[2];        // diagonal contains the number/proportion of emigrants!
    
    double _migr_rate[2];           // migration rate of males [0], and females [1]
    double _tot_emigRate[2];        // total emigration rate
    double _migr_rateIn[2];         // 2D-SS: migration from the edge to the middle (negative numbers mean that the ind get lost (absorbing boundaries))
    double _migr_rateOut[2];        // 2D-SS: migration from the edge to the outter
    double _migr_rateCorner[2];     // 2D-SS: migration from the corner to the middle
    
    double _migr_rate_propagule[2];
    double _disp_long_range_coef[2];// long range disperssal coefficient (currently for geometric function)
    
    virtual	bool _computeTotEmigrants(Patch* curPatch, unsigned int* nbEmigrants,
                                      double* migrTotRate, double* sum_m, double* factor);
    virtual unsigned int _computeTotEmigrants(Patch* curPatch, const double& migrTotRate,
                                              double& sum_m, const double& factor,
                                              const sex_t& SEX);
    
    virtual bool _sendEmigrants(Patch* curPatch,const unsigned int& home,
                                const unsigned int& target, unsigned int* totEmigr,
                                double* migrRate, double* sum_m, double* factor);
    virtual unsigned int _sendEmigrants(Patch* curPatch, const unsigned int& home,
                                        const unsigned int& target, unsigned int& totEmigr,
                                        const double& migrRate, double& sum_m,
                                        const double& factor, const sex_t& SEX);
    
    virtual void _sendEmigrant(Patch*, const unsigned int&,const unsigned int&,
                               const unsigned int&,const sex_t&,const unsigned int&);
    virtual void _sendEmigrant(Patch*, const unsigned int&, const unsigned int&,
                               const unsigned int&,const sex_t&);
    
    virtual void _absorbEmigrant(Patch* curPatch, const unsigned int& nbEmigr,
                                 const sex_t& SEX, const unsigned int& popSize);
    virtual void _absorbEmigrant(Patch* curPatch, const unsigned int& nbEmigr,
                                 const sex_t& SEX);
    
    virtual void migrate_2_neighbours(Patch* curPatch, const unsigned int& home,
                                      const unsigned int& n1, const unsigned int& n2,
                                      double* m1, double* m2, double* tot_m);
    virtual void migrate_4_neighbours(Patch* curPatch, const unsigned int& home,
                                      const unsigned int& n1, const unsigned int& n2,
                                      const unsigned int& n3, const unsigned int& n4,
                                      double* m1, double* m2, double* m3, double* m4,
                                      double* tot_m);
    virtual void migrate_8_neighbours(Patch* curPatch, const unsigned int& home,
                                      const unsigned int& n1, const unsigned int& n2,
                                      const unsigned int& n3, const unsigned int& n4,
                                      const unsigned int& n5, const unsigned int& n6,
                                      const unsigned int& n7, const unsigned int& n8,
                                      double* m1, double* m2, double* m3, double* m4,
                                      double* m5, double* m6, double* m7, double* m8,
                                      double* tot_m);
    
    unsigned int get_nbEmigrant(const double& m, double& sum_m, unsigned int& totEmigrants);
    
    void _get_lattice_dims();
    inline unsigned int Coord2ID(const unsigned int& x, const unsigned int& y){
        assert(_disp_model>=3 && x<_x_size && y<_y_size);
        return _x_size*y + x;
    }
    inline void ID2Coord(const unsigned int& id, unsigned int& x, unsigned int& y){
        assert(_disp_model>=3 && id<_nb_patch);
        y = id/_x_size;
        x = id - y*_x_size;
    }
    
    /** get the factor for the migration rate:
     * the _tot_emigraton_rate is passed such that the factor could be adjusted
     * if the total emigration rate would exceed 1.
     * This test is performed within the called function since get_migr_factor_one()
     * is most often called and here the test is not needed (perfomance...)
     */
    void _setDispersalFactor();
    void _setDispersalFactor(const sex_t& SEX);
    
    // function pointer to switch beteen immigraten and emigration
    void emigrate(sex_t SEX, age_idx from_age, unsigned int from_deme, age_idx to_age, unsigned int to_deme, unsigned int nbMigr){
        _popPtr->move_random(SEX, from_age, from_deme, to_age, to_deme, nbMigr);
    }
    void immigrate(sex_t SEX, age_idx from_age, unsigned int from_deme, age_idx to_age, unsigned int to_deme, unsigned int nbMigr){
        _popPtr->move_random(SEX, to_age, to_deme, from_age, from_deme, nbMigr);
    }
    
public:
    void _setDispersalFactor_friction(const sex_t& SEX);
    void _setDispersal_direction();

    void immigrate(Patch* curPatch, vector<Patch*> vNeighbours, unsigned int totSize, unsigned int nb_mgir, sex_t SEX, age_idx fromAge=OFFSx, age_idx toAge=ADLTx);
    
protected:
    inline void get_migr_factor(Patch* p, double* factor){
        factor[FEM] = (this->*get_migr_factor_funcPtr[FEM])(p, FEM); // female factor
        factor[MAL] = (this->*get_migr_factor_funcPtr[MAL])(p, MAL); // male factor
    }
    double (LCE_Disperse::* get_migr_factor_funcPtr[2])(Patch* p, sex_t s);
    inline double get_migr_factor_one(Patch* p, sex_t s){return 1.0;}

    inline double get_migr_factor_min(Patch* p, sex_t s){
        return _disp_factor[s][0];
    }
    inline double get_migr_factor_max(Patch* p, sex_t s){
        return _disp_factor[s][1];
    }
    inline double get_migr_factor_k_threshold(Patch* p, sex_t s){
        return (p->get_density(OFFSx)<_disp_factor[s][2]) ? get_migr_factor_min(p,s) : get_migr_factor_max(p,s);
    }
    inline double get_migr_factor_k_logistic(Patch* p, sex_t s){
        return generalLogisticCurve(p->get_density(OFFSx), _disp_factor[s][0],
                                    _disp_factor[s][1], _disp_factor[s][2],
                                    _disp_factor[s][3], _disp_factor[s][4]);
    }
    inline double get_migr_factor_saturation(Patch* p, sex_t s){
		double K(double(p->get_K(s)));
		double nb_ind(double(p->size(s, OFFSx)));
		double factor(0);
		if(nb_ind > 0) factor  = (nb_ind-K*(1-exp(-nb_ind/K)))/nb_ind;
		return factor ;}
    
    
    /** same as above but combined with the friction */
    inline double get_migr_factor_one_friction(Patch* p, sex_t s){
        return p->get_friction(s);
    }
    inline double get_migr_factor_min_friction(Patch* p, sex_t s){
        return _disp_factor[s][0] * p->get_friction(s);
    }
    inline double get_migr_factor_max_friction(Patch* p, sex_t s){
        return _disp_factor[s][1] * p->get_friction(s);
    }
    inline double get_migr_factor_k_threshold_friction(Patch* p, sex_t s){
        return (p->get_friction(s)*p->get_density(OFFSx)<_disp_factor[s][2]) ? get_migr_factor_min_friction(p,s) : get_migr_factor_min_friction(p,s);
    }
    inline double get_migr_factor_k_logistic_friction(Patch* p, sex_t s){
        return p->get_friction(s)*generalLogisticCurve(p->get_density(OFFSx), _disp_factor[s][0],
                                                       _disp_factor[s][1], _disp_factor[s][2],
                                                       _disp_factor[s][3], _disp_factor[s][4]);
    }
    
public:
    void (LCE_Disperse::* migration_func_ptr) ();
    
    virtual void migrate_zeroMigration();
    
    virtual void migrate_island           ( );
    virtual void migrate_island_propagule ( );
    
    virtual void migrate_matrix           ( );
    
    virtual void migrate_1D_ss( );
    virtual void migrate_1D_ss_all( );
    virtual void migrate_1D_ss_full( );
    
    virtual void migrate_2D_ss_4Neighbour();
    virtual void migrate_2D_ss_4Neighbour_all();
    virtual void migrate_2D_ss_4Neighbour_full();
    
    virtual void migrate_2D_ss_8Neighbour();
    virtual void migrate_2D_ss_8Neighbour_all();
    virtual void migrate_2D_ss_8Neighbour_full();
    
    virtual void migrate_2D_ss_geometric();
    virtual void migrate_2D_ss_geometric(const sex_t& SEX);
    
    virtual void immigrate_island( );
    
    vector<unsigned int> get_neighbours(const unsigned int& curPatch);
    vector<Patch*>       get_neighbours(Patch* curPatch);
    
    LCE_Disperse(int rank = my_NAN);
    
    virtual ~LCE_Disperse();
    
    virtual void execute ();
    
    ///@name Implementations
    ///@{
    virtual bool init(Metapop* popPtr);
    virtual LCE_Disperse* clone () {return new LCE_Disperse();}
    virtual void loadFileServices ( FileServices* loader ) {}
    virtual void loadStatServices ( StatServices* loader ) {}
    virtual age_t removeAgeClass ( ) {return OFFSPRG;}
    virtual age_t addAgeClass ( ) {return ADULTS;}
    virtual age_t requiredAgeClass () {return OFFSPRG;}
    
    ///@}
    
    void preDispersal();
    
    void (LCE_Disperse::*func_ptr_postDispersal)();
    virtual void postDispersal_noSample_withFull();
    virtual void postDispersal_noSample_noFull();
    virtual void postDispersal_withSample_withFull();
    virtual void postDispersal_withSample_noFull();
    virtual void set_func_ptr_withSample_withFull()   {func_ptr_postDispersal = &LCE_Disperse::postDispersal_withSample_withFull;}
    virtual void set_func_ptr_withSample_noFull()     {func_ptr_postDispersal = &LCE_Disperse::postDispersal_withSample_noFull;}
    virtual void set_func_ptr_noSample_withFull()     {func_ptr_postDispersal = &LCE_Disperse::postDispersal_noSample_withFull;}
    virtual void set_func_ptr_noSample_noFull()       {func_ptr_postDispersal = &LCE_Disperse::postDispersal_noSample_noFull;}
    
    ///@name Dispersal Matrix
    ///@{
    void setDispersalMatrix  ( );
    void setDispersalRate    ( );
    bool checkDispMatrix     (TMatrix* mat, sex_t sex);
    bool setDispMatrix       (TMatrix* mat);
    bool setDispMatrix       (sex_t sex, TMatrix* mat);
    ///@}
    
    int  get_disp_model      ( )  {return _disp_model;}
    string get_disp_model_str( );
    
    unsigned int get_border_model(){return _border_model;}
    unsigned int get_lattice_range(){return _lattice_range;}
    
    
    virtual void executeBeforeEachGeneration(const unsigned int& gen){}
    virtual void temporal_change(const unsigned int& gen);
    
    ///@}
};

class LCE_DisperseCoalescence: public LCE_Disperse {
public:
    LCE_DisperseCoalescence(int rank = my_NAN) : LCE_Disperse(rank){
    };
    ~LCE_DisperseCoalescence(){}
    
    void _sendEmigrant(Patch* curPatch, const unsigned int& home, const unsigned int& target,
                       const unsigned int& nbEmigr);
    void _sendEmigrant(Patch* curPatch, const unsigned int& home, const unsigned int& target,
                       const unsigned int& nbEmigr, const sex_t& SEX, const unsigned int& popSize);
    void migrate_zeroMigration();
    void migrate_island_propagule();
    void migrate_matrix();
    unsigned int _sendEmigrants(Patch* curPatch, const unsigned int& home, const unsigned int& target,
                                unsigned int& totEmigr, const double& migrRate, double& sum_m,
                                const double& factor);
    void migrate_island();
    bool _computeTotEmigrants(Patch* curPatch, unsigned int& totEmigr, const double& migrTotRate,
                              double& sum_m, const double& factor);
    void postDispersal_noSample_withFull();
    void postDispersal_noSample_noFull();
    void migrate_8_neighbours(Patch* curPatch, const unsigned int& home,
                              const unsigned int& n1, const unsigned int& n2, const unsigned int& n3, const unsigned int& n4,
                              const unsigned int& n5, const unsigned int& n6, const unsigned int& n7, const unsigned int& n8,
                              double* m1, double* m2, double* m3, double* m4,
                              double* m5, double* m6, double* m7, double* m8,
                              double* tot_m);
    void migrate_4_neighbours(Patch* curPatch, const unsigned int& home,
                              const unsigned int& n1, const unsigned int& n2, const unsigned int& n3, const unsigned int& n4,
                              double* m1, double* m2, double* m3, double* m4,
                              double* tot_m);
    void migrate_2_neighbours(Patch* curPatch, const unsigned int& home,
                              const unsigned int& n1, const unsigned int& n2,
                              double* m1, double* m2, double* tot_m);
    
    void execute();
};

#endif //LCEDISPERSE_H

