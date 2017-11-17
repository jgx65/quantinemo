/** @file random.h
 *
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

#ifndef uniformH
#define uniformH

#include <random>
#include <algorithm>
using namespace std;

/**Random number generation class, uses various types of random generation depending on the implementation.*/
class RAND {
private:
    
    
    std::mt19937 engine; // random number generator
    
    void initialize(unsigned int seed){engine.seed(seed);}

public:
    
    RAND() {set_seed((unsigned int)time(0));}
    RAND(unsigned int seed){set_seed(seed);}
    RAND(vector<unsigned int> seed){set_seed(seed);}
    ~RAND(){}
    
    void set_seed(unsigned int seed) {engine.seed(seed);}
    void set_seed(vector<unsigned int> seed){
        std::seed_seq seq(seed.begin(), seed.end());
        engine.seed(seq);
    }
    
    
    //------------------------------------------------------------------------------
    /** Returns a uniformly distributed random number from [0.0, 1.0[. **/
    double Uniform()
    {
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        return distribution(engine);
    }
    
    //------------------------------------------------------------------------------
    /** Returns a uniformly distributed entire random number from [0, max[.*/
    unsigned int Uniform(unsigned int max)
    {
        std::uniform_int_distribution<unsigned int> distribution(0, max-1);
        return distribution(engine);
    }
    
    //------------------------------------------------------------------------------
    /** Returns a uniformly distributed entire random number from [min, max[.*/
    unsigned int Uniform(unsigned int min, unsigned int max)
    {
        assert(max>min);
        std::uniform_int_distribution<unsigned int> distribution(min, max-1);
        return distribution(engine);
    }
    
    //------------------------------------------------------------------------------
    /**Returns a random boolean number.*/
    bool Bool()
    {
        std::bernoulli_distribution distribution(0.5);
        return distribution(engine);
    }
    
    //------------------------------------------------------------------------------
    /*Produces random non-negative integer values i, distributed according to 
     discrete probability function. The value obtained is the number of successes 
     in a sequence of t yes/no experiments, each of which succeeds with probability p.
     */
    inline unsigned int Binomial(double p, unsigned int N)
    {
        if(p>=1) return N; // due to rounding problems p may exceed slightly 1

        assert(N>0);
        assert(p>=0 && p<=1);
        std::binomial_distribution<> distribution(N, p);
        return distribution(engine);
        
    }
    
    //------------------------------------------------------------------------------
    // some standard distributions
    inline double Normal(double mean=0, double std_dev=1)
    {
        std::normal_distribution<double> distribution(mean, std_dev);
        return distribution(engine);
    }
    
    //------------------------------------------------------------------------------
    /* Produces random non-negative integer values i, distributed according to 
     discrete probability function. The value obtained is the probability of 
     exactly i occurrences of a random event if the expected, mean number of 
     its occurrence under the same conditions (on the same time/space interval) is μ.
     */
    inline unsigned int Poisson(double mean)
    {
        assert(mean>0.0);
        std::poisson_distribution<unsigned int> distribution(mean); // default std::poisson_distribution<int>
        return distribution(engine);
    }
    
    
    //------------------------------------------------------------------------------
    /* exponential_distribution expects the rate parameter lambda as input
     variable which is 1.0 divided by the desired mean
     */
    inline double random_exponential(double mean)
    {
        assert(mean>0.0);
        std::exponential_distribution<double> distribution(1.0 / mean);
        return distribution(engine);
    }
    
    
    //------------------------------------------------------------------------------
    /* Returns a geometrically distributed random variate between zero and infinity
     */
    unsigned int Geometric(double p)
    {
        if (p < 1.e-10 || p == 1.0) return 0;
        geometric_distribution<int> distribution(p);
        return distribution(engine);
    }
    
    //------------------------------------------------------------------------------
    /* Returns a exponentially distributed random variate between zero and infinity
     */
    double Exponential(double lamda)
    {
        if (lamda <= 0) return 0;
        std::exponential_distribution<double> distribution(lamda);
        return distribution(engine);
    }
    
    //------------------------------------------------------------------------------
    //Returns a exponentially distributed random variate between zero and infinity
    double Exponential_withMean(double mean)
    {
        if (mean <= 0) return 0;
        std::exponential_distribution<double> distribution(1.0/mean);
        return distribution(engine);
    }
    
    //------------------------------------------------------------------------------
    double Gamma(double a, double b)
    {
        if (b <= 0) throw("Gamma: scale must be positive!");
        std::gamma_distribution<double> distribution(a, b);
        return distribution(engine);
    }
    
    /* -----------------------------------------------------------------------------
     BetaRandom
     
     returns a variate that is Beta distributed on the interval [a,b]
     with shape parameters alpha and beta.
     
     The Beta function has two shaping parameters, alpha and beta.
     Setting these parameters to 1.5 and 1.5 yields a normal-like
     distribution, but without tails. If alpha and beta are equal to
     1 it is a uniform distribution.
     
     If alpha and beta are less than 1, use a rejection algorithm;
     Otherwise use the fact that if x is distributed Gamma(alpha) and y
     Gamma(beta) then x/(x+y) is Beta(alpha, beta).
     
     The rejection algorithm first a Beta variate is found over the
     interval [0, 1] with not the most efficient algorithm.  This is then
     scaled at the end to desired range.
     
     It may be tempting to re-use the second number drawn as the first
     random number of the next iteration, and simply draw one more.
     *** Don't do it.  You will produce an incorrect distribution.  You
     must draw two new numbers for the rejection sampling to be correct.
     
     References:
     - Ripley, Stochastic Simulations, John Wiley and Sons, 1987, p 90.
     - J.H.Maindonald, Statistical Computation, John Wiley and Sons,
     1984, p 370.
     */
    double Beta(double alpha, double beta, double a, double b)
    {
        if (b <= a) throw("Beta distribution: the ranges are mixed up!");
        return (a + Beta(alpha, beta) * (b - a)); /* Scale to interval [a, b] */
    }
    
    double Beta(double alpha, double beta)
    {
        double u1, u2, w;
        
        if (alpha <= 0 || beta <= 0)
        throw("Beta distribution: bad shape or scale!");
        
        if ((alpha < 1) && (beta < 1)) {
            /* use rejection */
            do {
                u1 = Uniform(); /* Draw two numbers */
                u2 = Uniform();
                
                u1 = pow(u1, 1 / alpha); /* alpha and beta are > 0 */
                u2 = pow(u2, 1 / beta);
                
                w = u1 + u2;
                
            } while (w > 1.0);
        }
        else {
            /* use relation to Gamma */
            u1 = Gamma(alpha);
            u2 = Gamma(beta);
            w = u1 + u2;
        }
        
        return u1 / w;
    } /* Beta*/
    
    double Gamma(double a)
    {
        int ia;
        double u, b, p, x, y = 0.0, recip_a;
        
        if (a <= 0) throw("Gamma: shape must be positive!");
        
        ia = (int) floor(a); /* integer part */
        a -= ia; /* fractional part */
        if (ia > 0) {
            y = igamma_dev(ia); /* gamma deviate w/ integer argument ia */
            if (a == 0.0) return (y);
        }
        
        /* get gamma deviate with fractional argument "a" */
        b = (M_E + a) / M_E;
        recip_a = 1.0 / a;
        for (;;) {
            u = Uniform();
            p = b * u;
            if (p > 1) {
                x = -log((b - p) / a);
                if (Uniform() > pow(x, a - 1)) continue;
                break;
            }
            else {
                x = pow(p, recip_a);
                if (Uniform() > exp(-x)) continue;
                break;
            }
        }
        return (x + y);
    }
    
    /****************************************************************
     gamma deviate for integer shape argument.  Code modified from pp
     292-293 of Numerical Recipes in C, 2nd edition.
     ****************************************************************/
    double igamma_dev(int ia)
    {
        int j;
        double am, e, s, v1, v2, x, y;
        
        if (ia < 1) throw("Gamma: argument must be >=1!");
        
        if (ia < 6) {
            x = 1.0;
            for (j = 0; j < ia; j++)
            x *= Uniform();
            x = -log(x);
        }
        else {
            do {
                do {
                    do { /* next 4 lines are equivalent */
                        v1 = 2.0 * Uniform() - 1.0; /* to y = tan(Pi * uni()).     */
                        v2 = 2.0 * Uniform() - 1.0;
                    } while (v1 * v1 + v2 * v2 > 1.0);
                    y = v2 / v1;
                    am = ia - 1;
                    s = sqrt(2.0 * am + 1.0);
                    x = s * y + am;
                } while (x <= 0.0);
                e = (1.0 + y * y) * exp(am * log(x / am) - s * y);
            } while (Uniform() > e);
        }
        return (x);
    }
    
    //----------------------------------------------------------------------------
    template<typename T>
    int AfterDistribution(T* array, int size, int isCumul)
    {
        std::discrete_distribution<> distribution(array, array+size);
        return distribution(engine);
    }
    
    //----------------------------------------------------------------------------
    /** array must be cumulative */
    template<typename T>
    int AfterDistribution(T* array, unsigned int size)
    {
        std::discrete_distribution<> distribution(array, array+size);
        return distribution(engine);
    }
    
     //----------------------------------------------------------------------------
    /* returns a variate such that the log of the variate is normally distributed.*/
    double LogNormal(const double& meanlog, const double& sdlog)
    {
        return exp(Normal(meanlog, sdlog));
        
        double mean=exp(meanlog + (sdlog*sdlog/2));
        double sd=sqrt(exp(sdlog*sdlog-1)*exp((2*meanlog)+(sdlog*sdlog)));
        lognormal_distribution<> distribution(mean, sd);
        return distribution(engine);
    }
    
    // ----------------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------------
    // Multivariate Normal random deviate (MNR)
    // ----------------------------------------------------------------------------------------
private:
    unsigned int _size;             // dimension of the matrix
    double* _covar;                 // variance-covariance matrix (size = _size^2)
    double* _corre;                 // pairwise correlation matrix (size = _size*(_size-1)/2)
    double* _meanv;                 // Mean vector of multivariate normal distribution (size = _size, values all set to 0)
    double* _e;                     // returned erandom deviates (size = _size)
    double* _tempMNR;               // temp array used to compute the deviates (size = _size)
    double* _param;                 // array filled by a call to setgmn()
    //     and used to to call genmn()    (size = (_size*(_size+3)/2 + 1)
    
    // ----------------------------------------------------------------------------------------
    // clear the arrays
    // ----------------------------------------------------------------------------------------
    /** clrear the arrays */
private:
    void _clear_MNR()
    {
        if (_covar) delete[] _covar;
        if (_corre) delete[] _corre;
        if (_meanv) delete[] _meanv;
        if (_e) delete[] _e;
        if (_tempMNR) delete[] _tempMNR;
        if (_param) delete[] _param;
    }
    
    // ----------------------------------------------------------------------------------------
    // initialize the arrays
    // ----------------------------------------------------------------------------------------
    /** initialize the arrays */
private:
    void _ini_MNR()
    {
        _size = 0;
        _covar = NULL;
        _corre = NULL;
        _meanv = NULL;
        _e = NULL;
        _tempMNR = NULL;
        _param = NULL;
    }
    // ----------------------------------------------------------------------------------------
    /** initialize the arrays */
private:
    void _ini_MNR(const unsigned int& s)
    {
        _size = s;                                // number fo dimensions
        
        if (_e) delete[] _e;                      // random deviate output array
        _e = new double[_size];
        
        if (_param) delete[] _param; // matrix filled by setgmn and used by genmn
        _param = new double[_size * (_size + 3) / 2 + 1];
        
        if (_meanv) delete[] _meanv; // Mean vector of multivariate normal distribution
        _meanv = new double[_size];               // here always 0
        for (unsigned int i = 0; i < _size; ++i) {
            _meanv[i] = 0;
        }
        
        if (_tempMNR) delete[] _tempMNR;         // working matrix used in genmn
        _tempMNR = new double[_size];
    }
    
    // ----------------------------------------------------------------------------------------
    // get random deviates
    // ----------------------------------------------------------------------------------------
    /** random deviates are returned. The last set variance-covaraince matrix is used */
public:
    double* get_MNR_deviates()
    {
        assert(_param);
        _genmn(_param, _e, _tempMNR);
        return _e;
    }
    
    // ----------------------------------------------------------------------------------------
    /** the MNR is initialized by passing an entire covariance matrix as a 2D array
     * and random deviates are returned
     */
public:
    double* get_MNR_deviates(double** matrix, const unsigned int& size)
    {
        double* covar = get_MNR_matrix2covar(matrix, size);
        _set_MNR(covar, size);
        delete[] covar;
        return get_MNR_deviates();
    }
    
    // ----------------------------------------------------------------------------------------
    /** the MNR is initialized by passing an entire covariance matrix as an 1D array
     * size is the dimension of the matrix and not the size of the array which is passed
     * the passed array should be of size = _size^2
     * and random deviates are returned
     */
public:
    double* get_MNR_deviates(double* covar, const unsigned int& size)
    {
        _set_MNR(covar, size);
        return get_MNR_deviates();
    }
    
    // ----------------------------------------------------------------------------------------
    /** the MNR is initialized by passing an entire G-matrix as an 1D array
     * size is the dimension of the G-matrix and not the size of the array which is passed
     * the passed array should be of size = _size^2
     * and random deviates are returned
     */
public:
    double* get_MNR_deviates(double* var, double* corr, const unsigned int& size)
    {
        double* covar = get_MNR_correl2covar(var, corr, size);
        _set_MNR(covar, size);
        delete[] covar;
        return get_MNR_deviates();
    }
    
    // ----------------------------------------------------------------------------------------
    /** random deviates are returned, following the "param" passed.
     * "Param" may be optained by calling set_MNR_getParam.
     * The returned array should not be deleted!
     */
public:
    double* get_MNR_deviates_fromParam(double* param, const unsigned int& size)
    {
        if (size != _size) _ini_MNR(size);
        _genmn(param, _e, _tempMNR);
        return _e;
    }
    
    // ----------------------------------------------------------------------------------------
    // generate and get the params used to call the get_MNR_deviate function efficently
    // ----------------------------------------------------------------------------------------
    /** returnes a pointer to the param array (array still belongs to this class)
     */
public:
    double* get_MNR_params()
    {
        return _param;
    }
    
    // ----------------------------------------------------------------------------------------
    /** returnes a copy of to the param array  (array belongs to the caller)
     */
public:
    double* get_MNR_params_copy()
    {
        unsigned int tot = _size * (_size + 3) / 2 + 1;
        double* a = new double[tot];
        for (unsigned int i = 0; i < tot; ++i) {
            a[i] = _param[i];
        }
        return a;
    }
    
    // ----------------------------------------------------------------------------------------
    // set the variance-covariance matrix
    // ----------------------------------------------------------------------------------------
    /** the MNR is initialized by passing an entire covaraince matrix as a single array
     * size is the dimension of the covariance and not the size of the array which is passed
     * the passed array should be of size = _size^2
     */
private:
    void _set_MNR(double* covar, const unsigned int& size) // params remain in object
    {
        if (!_size) _ini_MNR(size);
        _setgmn(_meanv, covar, _size, _param);   // _param is filled
    }
    
    // ----------------------------------------------------------------------------------------
    // Creation of the covaraince matrix
    // ----------------------------------------------------------------------------------------
    /** The covariance matrix is created by passing a 2D matrix
     * the returned array has to be deleted
     */
public:
    double* get_MNR_matrix2covar(double** matrix, const unsigned int& size)
    {
        
        double* m = new double[size * size];
        unsigned int t, c, l;    // 2D array has to be transformed to a 1D array
        for (t = 0, l = 0; t < size; ++t) {
            for (c = t + 1; c < size; ++c, ++l) {
                m[l] = matrix[t][c];
            }
        }
        return m;
    }
    
    // ----------------------------------------------------------------------------------------
    /** The covariance matrix is created and returned.
     * Input: - number of dimensions              (size)
     *        - variances                         (size)
     *        - pairwise correlation coefficients (size = _size*(_size-1)/2).
     * the covariances are approximated as follows: CoVar_1-2 = correlation_1-2*sqrt(Var_1*Var_2)
     * Example for 2 dimensions:
     *    [           var1           ] [correl_1_2*sqrt(var1*var2)]
     *    [correl_1_2*sqrt(var1*var2)] [          var2            ]
     * the returned array has to be deleted!
     */
public:
    double* get_MNR_correl2covar(double* var, double* corr, const unsigned int& size)
    {
        int tt, cc;
        unsigned int t, c, l;
        double* m = new double[size * size];
        
        // set the variances (diagonal)
        for (l = 0; l < size; ++l) {
            m[l + l * size] = var[l];
        }
        
        // set the covariances following the variances such
        for (t = 0, l = 0; t < size; ++t) {
            for (c = t + 1; c < size; ++c, ++l) {
                cc = c * size;
                tt = t * size;
                m[c * size + t] = m[t * size + c] = corr[l]
                * sqrt(var[t] * var[c]);
            }
        }
        return m;
    }
    
private:
    // ---------------------------------------------------------------------------
    /** GENerate Multivariate Normal random deviate
     
     Arguments
     parm --> Parameters needed to generate multivariate normal
     deviates (MEANV and Cholesky decomposition of
     COVM). Set by a previous call to SETGMN.
     1 : 1                - size of deviate, P
     2 : P + 1            - mean vector
     P+2 : P*(P+3)/2 + 1  - upper half of cholesky
     decomposition of cov matrix
     x    <-- Vector deviate generated.
     work <--> Scratch array
     
     Method
     1) Generate P independent standard normal deviates - Ei ~ N(0,1)
     2) Using Cholesky decomposition find A s.t. trans(A)*A = COVM
     3) trans(A)E + MEANV ~ N(MEANV,COVM)
     */
    void _genmn(double *parm, double *x, double *work)
    {
        static long i, icount, j, p, D1, D2, D3, D4;
        static double ae;
        p = (long) (*parm);
        
        // Generate P independent normal deviates - WORK ~ N(0,1)
        for (i = 1; i <= p; i++) {
            *(work + i - 1) = Normal();
        }
        for (i = 1, D3 = 1, D4 = (p - i + D3) / D3; D4 > 0; D4--, i += D3) {
            /*
             PARM (P+2 : P*(P+3)/2 + 1) contains A, the Cholesky
             decomposition of the desired covariance matrix.
             trans(A)(1,1) = PARM(P+2)
             trans(A)(2,1) = PARM(P+3)
             trans(A)(2,2) = PARM(P+2+P)
             trans(A)(3,1) = PARM(P+4)
             trans(A)(3,2) = PARM(P+3+P)
             trans(A)(3,3) = PARM(P+2-1+2P)  ...
             trans(A)*WORK + MEANV ~ N(MEANV,COVM)
             */
            icount = 0;
            ae = 0.0;
            for (j = 1, D1 = 1, D2 = (i - j + D1) / D1; D2 > 0; D2--, j += D1) {
                icount += (j - 1);
                ae +=
                (*(parm + i + (j - 1) * p - icount + p)
                 * *(work + j - 1));
            }
            *(x + i - 1) = ae + *(parm + i);
        }
        
    }
    
    // ---------------------------------------------------------------------------
    /** SET Generate Multivariate Normal random deviate Function
     Places P, MEANV, and the Cholesky factoriztion of COVM
     in GENMN.
     
     Arguments
     meanv  --> Mean vector of multivariate normal distribution.
     covm  <--> (Input) Covariance   matrix    of  the  multivariate
     normal distribution
     (Output) Destroyed on output
     p      --> Dimension of the normal, or length of MEANV.
     parm  <--  Array of parameters needed to generate multivariate norma
     deviates (P, MEANV and Cholesky decomposition of COVM).
     1 : 1                - P
     2 : P + 1            - MEANV
     P+2 : P*(P+3)/2 + 1  - Cholesky decomposition of COVM
     Needed dimension is (p*(p+3)/2 + 1)
     */
    void _setgmn(double *meanv, double *covm, long p, double *parm)
    {
        //extern void spofa(double *a,long lda,long n,long *info);
        static long T1;
        static long i, icount, info, j, D2, D3, D4, D5;
        T1 = p * (p + 3) / 2 + 1;
        
        // TEST THE INPUT
        if (p <= 0) throw("P nonpositive in SETGMN!");
        
        *parm = p;
        
        // PUT P AND MEANV INTO PARM
        for (i = 2, D2 = 1, D3 = (p + 1 - i + D2) / D2; D3 > 0; D3--, i += D2) {
            *(parm + i - 1) = *(meanv + i - 2);
        }
        
        // Cholesky decomposition to find A s.t. trans(A)*(A) = COVM
        _spofa(covm, p, p, &info);
        if (info != 0) throw("COVM not positive definite in SETGMN!");
        
        icount = p + 1;
        /*
         PUT UPPER HALF OF A, WHICH IS NOW THE CHOLESKY FACTOR, INTO PARM
         COVM(1,1) = PARM(P+2)
         COVM(1,2) = PARM(P+3)
         :
         COVM(1,P) = PARM(2P+1)
         COVM(2,2) = PARM(2P+2)  ...
         */
        for (i = 1, D4 = 1, D5 = (p - i + D4) / D4; D5 > 0; D5--, i += D4) {
            for (j = i - 1; j < p; j++) {
                icount += 1;
                *(parm + icount - 1) = *(covm + i - 1 + j * p);
            }
        }
        
    }
    
    // ---------------------------------------------------------------------------
    double _sdot(long n, double *sx, long incx, double *sy, long incy)
    {
        static long i, ix, iy, m, mp1;
        static double sdot, stemp;
        stemp = sdot = 0.0;
        if (n <= 0) return sdot;
        if (incx != 1 && incy != 1) {
            ix = iy = 1;
            if (incx < 0) ix = (-n + 1) * incx + 1;
            if (incy < 0) iy = (-n + 1) * incy + 1;
            for (i = 1; i <= n; i++) {
                stemp += (*(sx + ix - 1) * *(sy + iy - 1));
                ix += incx;
                iy += incy;
            }
            sdot = stemp;
            return sdot;
        }
        
        m = n % 5L;
        if (m != 0) {
            for (i = 0; i < m; i++) {
                stemp += (*(sx + i) * *(sy + i));
            }
            if (n < 5) {
                sdot = stemp;
                return sdot;
            }
        }
        
        mp1 = m + 1;
        for (i = mp1; i <= n; i += 5) {
            stemp += (*(sx + i - 1) * *(sy + i - 1) + *(sx + i) * *(sy + i)
                      + *(sx + i + 1) * *(sy + i + 1)
                      + *(sx + i + 2) * *(sy + i + 2)
                      + *(sx + i + 3) * *(sy + i + 3));
        }
        
        sdot = stemp;
        return sdot;
    }
    
    // ---------------------------------------------------------------------------
    /*
     SPOFA FACTORS A REAL SYMMETRIC POSITIVE DEFINITE MATRIX.
     SPOFA IS USUALLY CALLED BY SPOCO, BUT IT CAN BE CALLED
     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
     (TIME FOR SPOCO) = (1 + 18/N)*(TIME FOR SPOFA) .
     ON ENTRY
     A       REAL(LDA, N)
     THE SYMMETRIC MATRIX TO BE FACTORED.  ONLY THE
     DIAGONAL AND UPPER TRIANGLE ARE USED.
     LDA     INTEGER
     THE LEADING DIMENSION OF THE ARRAY  A .
     N       INTEGER
     THE ORDER OF THE MATRIX  A .
     ON RETURN
     A       AN UPPER TRIANGULAR MATRIX  R  SO THAT  A = TRANS(R)*R
     WHERE  TRANS(R)  IS THE TRANSPOSE.
     THE STRICT LOWER TRIANGLE IS UNALTERED.
     IF  INFO .NE. 0 , THE FACTORIZATION IS NOT COMPLETE.
     INFO    INTEGER
     = 0  FOR NORMAL RETURN.
     = K  SIGNALS AN ERROR CONDITION.  THE LEADING MINOR
     OF ORDER  K  IS NOT POSITIVE DEFINITE.
     LINPACK.  THIS VERSION DATED 08/14/78 .
     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
     SUBROUTINES AND FUNCTIONS
     BLAS SDOT
     FORTRAN SQRT
     INTERNAL VARIABLES
     */
    void _spofa(double *a, long lda, long n, long *info)
    {
        static long j, jm1, k;
        static double t, s;
        
        // BEGIN BLOCK WITH ...EXITS
        for (j = 1; j <= n; j++) {
            *info = j;
            s = 0.0;
            jm1 = j - 1;
            if (jm1 >= 1) {
                for (k = 0; k < jm1; k++) {
                    t = *(a + k + (j - 1) * lda)
                    - _sdot(k, (a + k * lda), 1L, (a + (j - 1) * lda),
                            1L);
                    t /= *(a + k + k * lda);
                    *(a + k + (j - 1) * lda) = t;
                    s += (t * t);
                }
            }
            
            s = *(a + j - 1 + (j - 1) * lda) - s;
            
            // ......EXIT
            if (s <= 0.0) return;
            *(a + j - 1 + (j - 1) * lda) = sqrt(s);
        }
        *info = 0;
        return;
        
    }
    
    //------------------------------------------------------------------------------
public:
    unsigned int operator()(unsigned int x)
    {
        return Uniform(x);
    }
    double operator()()
    {
        return Uniform();
    }
    
    //------------------------------------------------------------------------------
    /** randomize the order of the array
     */
    template<typename T>
    void randomize(T* & array, const unsigned int& size)
    {
        random_shuffle(array, array + size, *this);
    }
    
    //------------------------------------------------------------------------------
    /** randomize the order of the vector
     */
    template<typename T>
    void randomize(vector<T> & vec)
    {
        random_shuffle(vec.begin(), vec.end(), *this);
    }
    
    //------------------------------------------------------------------------------
    /** return a new array of size "sampleSize" sampled randomly from the array "array"
     * with or without replacement
     * "array" is not modified
     * the returned array has to be deleted
     */
    template<typename T>
    T* sample(T* & array, const unsigned int& size,
              const unsigned int& sampleSize, bool replacement = false)
    {
        T* array_out;
        
        if (replacement) {        // sampling with replacement
            // draw randomly the indexes
            array_out = new T[sampleSize];
            for (unsigned int i = 0; i < sampleSize; ++i) {
                array_out[i] = array[Uniform(size)];
            }
        }
        else {                   // sampling without replacement
            // copy the array, randomize the order an return the entire array although it may be too long (should not cause a problem...)
            assert(sampleSize<=size);
            array_out = new T[size];
            ARRAY::copy(array_out, array, size);
            randomize(array_out, size);     // randomize its order
        }
        
        // return the array
        return array_out;
    }
    
    //------------------------------------------------------------------------------
    /** return a new vector of size "sampleSize" sampled randomly from the array "array"
     * with or without replacement
     * "vec" is not modified
     */
    template<typename T>
    vector<T> sample(vector<T> & vec, const unsigned int& sampleSize,
                     bool replacement = false)
    {
        vector<T> out;
        unsigned int i, size = (unsigned int)vec.size();
        out.reserve(sampleSize);
        
        if (replacement) {          // sampling with replacement
            for (i = 0; i < sampleSize; ++i) {
                out.push_back(vec[Uniform(size)]);  // draw randomly the indexes
            }
        }
        else {                     // sampling without replacement
            // copy the array, randomize the order and return the entire array although it may be too long (should not cause a problem...)
            if (sampleSize > size)
            throw("sampling size exceeds number of elements!");
            if (size / 2 < sampleSize) {   // more efficient to draw the loosers
                vector<unsigned int> vIndex = sample<unsigned int>(0, size,
                                                                   size - sampleSize, false, true);	// the vector has to be sorted
                vIndex.push_back(size);  // add a looser after the range (never reached) to allow to simply loop
                vector<unsigned int>::iterator curIndex = vIndex.begin();
                for (i = 0; i < size; ++i) {
                    if (i == *curIndex) ++curIndex; 	// index is present: this is a looser
                    else out.push_back(vec[i]); // this is a winner: copy it
                }
            }
            else {                     // more efficient to draw the winners
                vector<unsigned int> vIndex = sample<unsigned int>(0, size,
                                                                   sampleSize);
                vector<unsigned int>::iterator curIndex = vIndex.begin(),
                endIndex = vIndex.end();
                for (; curIndex != endIndex; ++curIndex) {
                    out.push_back(vec[*curIndex]);  // this is a winner: copy it
                }
            }
        }
        
        // return the array
        return out;
    }
    
    //------------------------------------------------------------------------------
    /** return a new vector of size "size" with randomly sampled values between min(included) and max (excluded)
     * with or without replacement
     */
    template<typename T>
    vector<T> sample(T min, T max, unsigned int sampleSize, bool replacement =
                     false, bool sortedOutput= false)
    {
        vector<T> out;
        T range = max - min;
        unsigned int i;
        
        if (replacement) {          // sampling with replacement
            // draw randomly the indexes
            out.reserve(sampleSize);
            for (i = 0; i < sampleSize; ++i) {
                out.push_back(min + (T) Uniform() * range);
            }
        }
        else {                     // sampling without replacement
            if (sampleSize > range)
            throw("sampling size exceeds number of elements!");
            
            // if few indexes have to be drawn
            if (range > 10 * sampleSize) { // method draw and check if they were already drawn
                unsigned int j;
                T value;
                out.reserve(sampleSize);
                for (i = 0; i < sampleSize;) {
                    value = min + (T) (Uniform() * range); // draw a value
                    for (j = 0; j < i; ++j) { // check if the value is already drawn
                        if (out[j] == value) break;
                    }
                    if (j != i) continue;        // is the value already drawn
                    out.push_back(value);     // it is a new value
                    ++i;
                }
            }
            else { // the entire index is generated and then the values are randomly swaped
                out.reserve((unsigned int) range);
                for (T j = min; j < max; ++j) {    // generate the ordered index
                    out.push_back(j);
                }
                for (i = 0; i < sampleSize; ++i) { // sample the first "sampleSize" elements randomly
                    swap(out[i], out[i + Uniform((unsigned int) range - i)]);
                }
                out.resize(sampleSize);
            }
        }
        
        // return the array
        if(sortedOutput) sort(out.begin(), out.end()); // sort the vector if needed
        return out;
    }
    
    //------------------------------------------------------------------------------
    /** similar to sample in R:
     * returns nb UNIQUE integer numbers between 0 and max-1
     * the returned array must be deleted
     * the returned array is sorted!!!
     */
    unsigned int* sample(const unsigned int& nb, const unsigned int& max)
    {
        assert(nb<=max);
        
        // create a temp vector (0, 1, 2, ... max-1)
        vector<unsigned int> vec;
        for(unsigned int i = 0; i< max; ++i){
            vec.push_back(i);
        }
        
        unsigned int* array = new unsigned int[nb];     // array to return
        if(nb<max/2){                                   // draw the "good" ones
            unsigned int pos;
            for(unsigned int i = 0; i< nb; ++i){
                pos = Uniform(max-i);          // get the pos
                array[i] = vec[pos];                        // store the value to the array
                vec.erase(vec.begin()+pos);                 // remove the element
            }
        }
        else{                                           // remove the bad ones
            for(unsigned int i = 0; i< max-nb; ++i){
                vec.erase(vec.begin()+Uniform(max-i));
            }
            for(unsigned int i = 0; i< nb; ++i){          // copy the vector to the array
                array[i] = vec[i];
            }
        }
        return array;
    }
    
    
};
#endif
