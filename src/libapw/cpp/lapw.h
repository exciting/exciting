#ifndef __LAPW_H__
#define __LAPW_H__

/*! \file lapw.h
    \brief Defines basic classes, variables and inline functions.
    
    \defgroup functions global functions
*/

#include <vector>
#include <string>
#include <map>
#include <complex>
#include <algorithm>
#include <iomanip>
#include "typedefs.h"
#include "config.h"
#include "mdarray.h"
#include "linalg.h"
#include "timer.h"

void lapw_fft(int32_t direction, complex16 *data);
double lapw_spline_integrate(int32_t n, double *x, double *f);

/// describes single atomic level
struct atomic_level 
{
    /// principal quantum number
    int n;
    
    /// angular momentum quantum number
    int l;
    
    /// quantum number k
    int k;
    
    /// level occupancy
    int occupancy;
};

/// describes radial solution
struct radial_solution_descriptor
{
    /// principal quantum number
    int n;

    /// angular momentum quantum number
    int l;

    /// order of energy derivative
    int dme;

    /// energy of the solution
    double enu;

    /// automatically determine energy
    bool auto_enu;
};

/// collection of radial solutions for a given angular momentum quantum number l
class radial_l_descriptor
{
    public:
        radial_l_descriptor()
        {
        }
        
        radial_l_descriptor(int l) : l(l)
        {
        }
    
        /// angular momentum quantum number
        int l;

        /// list of radial solution descriptors
        std::vector<radial_solution_descriptor> radial_solution_descriptors;
};

struct mt_radial_function_descriptor
{
    /// l-quantum number
    int l;
    
    /// order of radial function for a given l (\f$ u_l,\dot{u_l},... \f$ from the APW set, then local orbitals)
    int order;
    
    /// index of local orbital in the list of local orbitals for a given species
    int idxlo;
};

class mt_radial_index
{
    public:

        mt_radial_index() : lmax(0), nrflmax(0)
        {
        }
        
        /// add radial function for an angular momentum quantum number l to the list of radial functions
        void add(int l, int idxlo = -1)
        {
            lmax = std::max(lmax, l);

            nrfl.resize(lmax + 1, 0);

            mt_radial_function_descriptor rfd;
            rfd.l = l;
            rfd.order = nrfl[l]++;
            rfd.idxlo = idxlo;
            
            radial_function_descriptors.push_back(rfd);
            
            nrflmax = std::max(nrflmax, nrfl[l]);
        }

        /// initialize backward mapping
        void init()
        {
            index_by_lo.set_dimensions(lmax + 1, nrflmax);
            index_by_lo.allocate();
            index_by_lo.zero();

            for (int i = 0; i < (int)radial_function_descriptors.size(); i++)
                index_by_lo(radial_function_descriptors[i].l, radial_function_descriptors[i].order) = i;
        }

        inline const mt_radial_function_descriptor& operator [](int i)
        {
            return radial_function_descriptors[i];
        }
        
        inline int operator()(int l, int order)
        {
            return index_by_lo(l, order);
        }
         
        /// number of radial functions for a given l
        inline int nrf(int l)
        {
            return nrfl[l];
        }
        
        /// total number of radial functions
        inline int size()
        {
            return radial_function_descriptors.size();
        }

        inline int getlmax()
        {
            return lmax;
        }
        
        inline int getnrflmax()
        {
            return nrflmax;
        }

    private:

        /// maximum l for which indices are built
        int lmax;

        /// list of descriptors for radial functions \f$ f_{\ell \lambda}^{\alpha}(r) \f$
        std::vector<mt_radial_function_descriptor> radial_function_descriptors;
        
        /// number of radial functions for a given l
        std::vector<int> nrfl;

        /// maximum number of radial functions over all values of l
        int nrflmax;

        mdarray<int,2> index_by_lo;
};

struct mt_function_descriptor
{
    /// l-quantum number
    int l;

    /// magnetic quantum number
    int m;

    /// composite lm index
    int lm;
    
    /// order of radial function for a given l (\f$ u_l,\dot{u_l},... \f$ from the APW set, then local orbitals)
    int order;

    /// index of radial function
    int idxrf;
    
    /// index of local orbital in the list of local orbitals for a given species
    int idxlo;
};

/*! \brief muffin-tin combined indices

    Inside each muffin-tin sphere \f$ \alpha \f$ wave-functions are expanded in radial functions and spherical harmonics:
    \f[
        \psi_{i {\bf k}}({\bf r})= \sum_{L} \sum_{\lambda=1}^{N_{\ell}^{\alpha}} 
           F_{L \lambda}^{i {\bf k},\alpha}f_{\ell \lambda}^{\alpha}(r) 
           Y_{\ell m}(\hat {\bf r}) 
    \f]
    where functions \f$ f_{\ell \lambda}^{\alpha}(r) \f$ label both APW and local orbitals.
*/
class mt_index
{
       
    public:
        
        mt_index() : apw_function_descriptors_size(0), lo_function_descriptors_size(0)
        {
        }
        
        /// composite lm index by l and m
        inline static int idxlm(int l, int m)
        {
            return l * l + l + m;
        }

        /// initialize muffin-tin indices and backward mapping
        void init(mt_radial_index& ri)
        {
            for (int i = 0; i < ri.size(); i++)
            {
                int l = ri[i].l;
                for (int m = -l; m <= l; m++)
                {
                    mt_function_descriptor fd;
                    fd.l = l;
                    fd.m = m;
                    fd.lm = idxlm(l, m);
                    fd.order = ri[i].order;
                    fd.idxrf = i;
                    fd.idxlo = ri[i].idxlo;
                    function_descriptors.push_back(fd);
                }
                if (ri[i].idxlo == -1)
                    apw_function_descriptors_size = function_descriptors.size();
            }
            
            lo_function_descriptors_size = function_descriptors.size() - apw_function_descriptors_size;

            int lmmax = pow(ri.getlmax() + 1, 2);
            index_by_lmo.set_dimensions(lmmax, ri.getnrflmax());
            index_by_lmo.allocate();
            index_by_lmo.zero();

            for (int i = 0; i < (int)function_descriptors.size(); i++)
                index_by_lmo(function_descriptors[i].lm, function_descriptors[i].order) = i;
        }
        
        /// total number of basis functions
        inline int size()
        {
            return function_descriptors.size();
        }
        
        /// number of APW basis functions
        inline int apw_size()
        {
            return apw_function_descriptors_size;
        }

        /// number of local orbital basis functions        
        inline int lo_size()
        {
            return lo_function_descriptors_size;
        }
        
        inline const mt_function_descriptor& operator[](int i)
        {
            return function_descriptors[i];
        }

        inline int operator()(int lm, int order)
        {
           return index_by_lmo(lm, order);
        }

        inline int operator()(int l, int m, int order)
        {
           return index_by_lmo(idxlm(l, m), order);
        }

    private:

        /// size of the APW part of basis function descriptors
        int apw_function_descriptors_size;

        /// size of the local orbital part of basis function descriptors
        int lo_function_descriptors_size;
        
        /// list of descriptors for basis basis functions \f$ f_{\ell \lambda}^{\alpha}(r)Y_{\ell m}(\hat {\bf r}) \f$
        std::vector<mt_function_descriptor> function_descriptors;

        /// mapping from l,m,order to a global index of basis function
        mdarray<int,2> index_by_lmo;
};

/// species related variables
class Species 
{
    public:

        Species() 
        {
        }
    
        Species(const std::string& symbol) : symbol(symbol) 
        {
        }
        
        std::string name;
        std::string symbol;
        int number;
        double mass;
        double rmin;
        double rmax;
        double rmt;
        int nrmt;
        
        /// l for U correction
        int lu;
  
        /// list of core levels 
        std::vector<atomic_level> core;

        /// list of radial descriptors used to construct local orbitals
        std::vector<radial_l_descriptor> lo_descriptors;

        /// list of radial descriptors used to construct apw functions 
        std::vector<radial_l_descriptor> apw_descriptors;

        mt_radial_index radial_index;
        mt_index index;
        
        mdarray<double,1> radial_mesh;
};

class Atom 
{
    public:
        
        Atom()
        {
        }

        Atom(Species *species, int ic) : species(species), idxclass(ic)
        {
        }
        
        /// atom position in lattice coordinates
        double posl[3];
        
        /// atom position in Cartesian coordinates
        double posc[3];
        
        /// muffin-tin magnetic field in Cartesian coordinates
        double bfcmt[3];
        
        /// pointer to species class
        Species *species;
        
        /// offset of APW matching coefficients in the array apwalm 
        int offset_apw;
        int offset_lo;
        int offset_wfmt;
        
        /// index of symmetry class
        int idxclass;
};

/// global read-only variables initialized once at the beginning
struct lapw_global_variables
{
    /// maximum number of G+k vectors
    int ngkmax;
    
    /// maximum order of APW functions 
    int apwordmax;
    
    /// maximum l for APW functions
    int lmaxapw;
    
    /// (lmaxapw+1)^2
    int lmmaxapw;
    
    /// maximum l for potential and density
    int lmaxvr;
    
    /// (lmaxvr+1)^2
    int lmmaxvr;
    
    int lmaxlu;
    
    int lmmaxlu;
    
    bool ldapu;
    
    /// number of G-vectors
    int ngvec;
    
    /// FFT grid dimensions
    int ngrid[3];
    
    /// number of FFT grid points
    int ngrtot;
    
    /// maximum number of radial points
    int nrmtmax;
    
    /// number of first-variational states
    int nstfv;
    
    /// number of second-variational states
    int nstsv;
    
    /// maximum number of radial functions over all species
    int nrfmtmax;
    
    /// maximum number (order) order of radial functions over all species and l-channels 
    int ordrfmtmax;
    
    /// tolerance for generalized eigen-value solver
    double evaltol;
    
    /// spin-polarized calculation (T/F)
    bool spinpol;
    
    bool spinorb;
    
    /// number of dimensions (0,1 or 3) of the magnetization field
    int ndmag;
    
    /// number of spinor components (1 or 2)
    int nspinor;

    /// mapping from G-vector index to FFT grid
    std::vector<int> igfft;

    /// plane-wave coefficients of the characteristic function 
    std::vector<complex16> cfunig;

    /// characteristic function on the real-space grid
    std::vector<double> cfunir;
    
    /// FFT grid range
    mdarray<int,2> intgv;
    
    /// G-vector integer coordinates
    mdarray<int,2> ivg;
    
    /// mapping from G-vector integer coordinates to global G-vector index
    mdarray<int,3> ivgig;
    
    /// Gaunt coefficients <Y|R|Y>
    mdarray<complex16,3> gntyry;

    /// indices of non-zero gntyry elements
    mdarray<std::vector<int>,2> L3_gntyry;

    /// non-zero gntyry elements
    //mdarray<std::vector<complex16>,2> L3_gntyry_data;
    
    int size_wfmt_apw;
    int size_wfmt_lo;
    int size_wfmt;
    std::vector<int> l_by_lm;
    
    int natmcls;
    std::vector<int> ic2ias;
    std::vector<int> natoms_in_class;
    
    /// direct lattice vectors
    std::vector< std::vector<double> > avec;
    
    /// 3x3 matrix of lattice vectors 
    double avec_m[9];
    
    /// inverse matrix of lattice vectors
    double ainv_m[9];
    
    /// volume of the unit cell
    double omega;
    
    /// reciprocal lattice vectors
    std::vector< std::vector<double> > bvec;
    
    /// 3x3 matrix of reciprocal vectors
    double bvec_m[9];
    
    /// inverse matrix of reciprocal vectors
    double binv_m[9];
    
    /// list of species
    std::vector<Species*> species;
    
    /// list of atoms
    std::vector<Atom*> atoms;

    int max_mt_index_size;

};
extern lapw_global_variables lapw_global;

/// describes the Bloch-states for a specific k-point 
class bloch_states_k
{
    public:

        bloch_states_k(int ngk);

        void copy_apwalm(complex16 *apwalm_);
        
        void generate_scalar_wave_functions(); 
        
        void generate_spinor_wave_functions(diagonalization mode);
        
        void test_scalar_wave_functions(int use_fft);
        
        void test_spinor_wave_functions(int use_fft);
        
        /// number of G+k vectors
        int ngk;

        /// number of LAPW basis functions (number of plane-waves + number of local orbitals)
        int lapw_basis_size;

        /// number of wave-function expansion coefficients (number of plane-waves + number of MT lm- coefficients)
        int wave_function_size;

        /// G+k vectors in Cartesian coordinates
        mdarray<double,2> vgkc;
        
        /// index of G-vector 
        std::vector<int> idxg;
        
        /// position of the G-vector in FFR arrays
        std::vector<int> idxgfft;
        
        /// weight of k-point
        double weight;
        
        /// APW matching coefficients
        mdarray<complex16,2> apwalm;

        /// first-variational eigen vectors
        mdarray<complex16,2> evecfv;

        /// second-variational eigen vectors
        mdarray<complex16,2> evecsv;

        /// full-diagonalization eigen vectors
        mdarray<complex16,2> evecfd;

        /// scalar wave functions
        mdarray<complex16,2> scalar_wave_functions;

        /// spinor wave functions
        mdarray<complex16,3> spinor_wave_functions;

        /// eigen-values of the first-variational equation
        std::vector<double> evalfv;

        /// eigen-values of the second-variational equation
        std::vector<double> evalsv;

        /// occupation of the second-variational states
        std::vector<double> occsv;
};


/// global variables that change during the run time
struct lapw_runtime_variables
{
    mdarray<double,4> hmltrad;
    mdarray<double,4> ovlprad;
    mdarray<double,4> socrad;
    mdarray<double,5> beffrad;
    mdarray<double,5> apwfr;
    mdarray<double,3> apwdfr;
    mdarray<double,2> beffir;
    mdarray<complex16,2> beffig;
    mdarray<complex16,1> veffig;

    /// collection of bloch states for a local fraction of k-points
    std::vector<bloch_states_k*> bloch_states;
    
    mdarray<complex16,5> dmatu;
    mdarray<complex16,5> vmatu;
    
    mdarray<double,3> rfmt;
};
extern lapw_runtime_variables lapw_runtime;

extern complex16 zone;
extern complex16 zzero;
extern complex16 zi;
extern double y00;

inline int idxG12(int ig1, int ig2)
{
    return lapw_global.ivgig(lapw_global.ivg(0, ig1) - lapw_global.ivg(0, ig2), 
                             lapw_global.ivg(1, ig1) - lapw_global.ivg(1, ig2),
                             lapw_global.ivg(2, ig1) - lapw_global.ivg(2, ig2));
}

template <typename T>
inline void L3_sum_gntyry(int lm1, int lm2, T *v, std::complex<double>& zsum)
{
    for (int k = 0; k < (int)lapw_global.L3_gntyry(lm1, lm2).size(); k++)
    {
        int lm3 = lapw_global.L3_gntyry(lm1, lm2)[k];
        zsum += lapw_global.gntyry(lm3, lm1, lm2) * v[lm3];
    }
}

inline bool use_spin_block(int ispn1, int ispn2)
{
    return ((lapw_global.ndmag == 1 && ispn1 == ispn2) || (lapw_global.ndmag != 1));
}

#endif // __LAPW_H__

