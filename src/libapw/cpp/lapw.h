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

/// composite lm index by l and m
inline int idxlm(int l, int m)
{
    return l * l + l + m;
}

/// describes single atomic level
struct atomic_level 
{
    /// principal quantum number
    unsigned int n;
    
    /// angular quantum number
    unsigned int l;
    
    /// quantum number k
    unsigned int k;
    
    /// level occupancy
    int occupancy;
};

/// describes radial solution
struct radial_solution_descriptor
{
    /// principal quantum number
    unsigned int n;

    /// angular quantum number
    unsigned int l;

    /// order of energy derivative
    unsigned int dme;

    /// energy of the solution
    double enu;

    /// automatically determine energy
    bool auto_enu;
};

/// collection of radial solutions for a given angular quantum number l
class radial_l_descriptor
{
    public:
        radial_l_descriptor()
        {
        }
        
        radial_l_descriptor(unsigned int l) : l(l)
        {
        }
    
        /// angular quantum number
        unsigned int l;

        /// list of radial solution descriptors
        std::vector<radial_solution_descriptor> radial_solution_descriptors;
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

class mtci
{
    public:
        mtci(unsigned int l, int m, unsigned int order, unsigned int idxrf) : l(l), m(m), order(order), idxrf(idxrf), idxlo(-1)
        {
            lm = idxlm(l, m);
        }
        mtci(unsigned int l, int m, unsigned int order, unsigned int idxrf, unsigned int idxlo) : l(l), m(m), order(order), idxrf(idxrf), idxlo(idxlo)
        {
            lm = idxlm(l, m);
        }
        
        unsigned int l;
        int m;
        unsigned int lm;
        unsigned int order;
        unsigned int idxrf;
        unsigned int idxlo;
};

/*struct mt_composite_index_descriptor
{
    unsigned int l;
    int m;
    unsigned int lm;

    /// order of radial function for a givel l (u_l, \dot{u_l}... from the APW set, then local orbitals)
    unsigned int ordrf;

    /// index of radial function in the array of radial functions
    unsigned int idxrf;

    /// index of local orbital in the list of local orbitals for a given species
    unsigned int idxlo;
};*/

class index_mapping
{
    private:

        struct radial_function_descriptor
        {
            /// l-quantum number
            unsigned int l;
            
            /// order of radial function for a given l (\f$ u_l,\dot{u_l},... \f$ from the APW set, then local orbitals)
            unsigned int ordrf;
            
            /// index of radial function in the array of radial functions
            unsigned int idxrf;

            /// index of local orbital in the list of local orbitals for a given species
            unsigned int idxlo;
        };

        class basis_function_descriptor
        {
            public:
                basis_function_descriptor(const radial_function_descriptor& rfd) : rfd(rfd)
                {
                }

                int m;
                unsigned int lm;
                radial_function_descriptor rfd;
        };

        /// maximum l for which indices are build
        unsigned int lmax;

        /// (lmax+1)^2
        unsigned int lmmax;
       
        /// size of the APW part of basis function descriptors
        unsigned int apw_basis_function_descriptors_size;

        /// size of the local orbital part of basis function descriptors
        unsigned int lo_basis_function_descriptors_size;
        
        /// list of descriptors for basis basis functions \f$ f_{\ell \lambda}^{\alpha}(r)Y_{\ell m}(\hat {\bf r}) \f$
        std::vector<basis_function_descriptor> basis_function_descriptors;

        /// local-orbital part of basis function descriptors
        basis_function_descriptor *lo_basis_function_descriptors;

        /// list of descriptors for radial functions \f$ f_{\ell \lambda}^{\alpha}(r) \f$
        std::vector<radial_function_descriptor> radial_function_descriptors;
        
        /// number of radial functions for a given l
        std::vector<unsigned int> nrfl;

        /// maximum number of radial functions over all values of l
        unsigned int nrflmax;

        /// mapping from l,m,order to a global index of basis function
        mdarray<unsigned int,2> idxbf_by_lmo;
        
        mdarray<unsigned int,2> idxrf_by_lo;

        /// l for a given radial function
        std::vector<int> l_by_idxrf;
    
    public:
        
        index_mapping() : lmax(0), 
                          apw_basis_function_descriptors_size(0),
                          lo_basis_function_descriptors_size(0),
                          lo_basis_function_descriptors(0), 
                          nrflmax(0)
        {
            radial_function_descriptors.clear();
            basis_function_descriptors.clear();
        }

        static inline unsigned int idxlm(unsigned int l, int m)
        {
            return l * l + l + m;
        }

        void add(unsigned int l, int idxlo = -1)
        {
            lmax = std::max(lmax, l);

            nrfl.resize(lmax + 1, 0);

            radial_function_descriptor rfd;
            rfd.l = l;
            rfd.ordrf = nrfl[l]++;
            nrflmax = std::max(nrflmax, nrfl[l]);

            rfd.idxrf = radial_function_descriptors.size();
            rfd.idxlo = idxlo;
            
            radial_function_descriptors.push_back(rfd);

            basis_function_descriptor bfd(rfd);

            for (int m = -l; m <= (int)l; m++)
            {
                bfd.m = m;
                bfd.lm = idxlm(l, m);
                basis_function_descriptors.push_back(bfd);
                
                if (idxlo != -1 && lo_basis_function_descriptors == 0)
                {
                    lo_basis_function_descriptors = &basis_function_descriptors.back();
                    apw_basis_function_descriptors_size = basis_function_descriptors.size() - 1;
                } 
            }
        }
        
        void init()
        {
            lmmax = pow(lmax + 1, 2);

            if (apw_basis_function_descriptors_size == 0)
            {
                apw_basis_function_descriptors_size = basis_function_descriptors.size();
            }
            else
            {
                lo_basis_function_descriptors_size = basis_function_descriptors.size() - apw_basis_function_descriptors_size;
            }

            idxbf_by_lmo.set_dimensions(lmmax, nrflmax);
            idxbf_by_lmo.allocate();
            idxbf_by_lmo.zero();

            for (unsigned int i = 0; i < basis_function_descriptors.size(); i++)
                idxbf_by_lmo(basis_function_descriptors[i].lm, basis_function_descriptors[i].rfd.ordrf) = i;
            
            idxrf_by_lo.set_dimensions(lmax + 1, nrflmax);
            idxrf_by_lo.allocate();
            idxrf_by_lo.zero();

            for (unsigned int i = 0; i < radial_function_descriptors.size(); i++)
                idxrf_by_lo(radial_function_descriptors[i].l, radial_function_descriptors[i].ordrf) = i;
            
            
            /*for (unsigned int i = 0; i < basis_function_descriptors.size(); i++)
            {
                std::cout << "idxbf=" << i 
                          << "  l=" << basis_function_descriptors[i].rfd.l
                          << "  ord=" << basis_function_descriptors[i].rfd.ordrf
                          << "  m=" << basis_function_descriptors[i].m << std::endl;


            }*/
        }
        
        inline unsigned int getidxrf(unsigned int l, unsigned int ordrf)
        {
            return idxrf_by_lo(l, ordrf);
        } 
        
        inline unsigned int getnrf(unsigned int l)
        {
            return nrfl[l];
        }

        inline unsigned int getnbf()
        {
            return basis_function_descriptors.size();
        }

        inline unsigned int getbfl(unsigned int idxbf)
        {
            return basis_function_descriptors[idxbf].rfd.l;
        }
        
        inline unsigned int getbfm(unsigned int idxbf)
        {
            return basis_function_descriptors[idxbf].m;
        }

        inline unsigned int getbflm(unsigned int idxbf)
        {
            return basis_function_descriptors[idxbf].lm;
        }
        
        inline unsigned int getbfordrf(unsigned int idxbf)
        {
            return basis_function_descriptors[idxbf].rfd.ordrf;
        }

        inline unsigned int getidxbf(unsigned int lm, unsigned int ordrf)
        {
            return idxbf_by_lmo(lm, ordrf);
        }
        
        inline unsigned int getidxbf(unsigned int l, int m, unsigned int ordrf)
        {
            return idxbf_by_lmo(idxlm(l, m), ordrf);
        }


};

/// species related variables
class Species 
{
    public:

        Species() : size_ci_lo(0), size_ci_apw(0), nrfmt(0)
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

        /// list of combined indices for the muffin-tin representation of wave-functions    
        std::vector<mtci> ci;

        /// local-orbital part of combined indices
        mtci *ci_lo;
        unsigned int size_ci_lo;
        unsigned int size_ci_apw;
        mdarray<int,2> ci_by_lmo;
        std::vector<int> ci_by_idxrf;
        std::vector<int> l_by_idxrf;
        std::vector<int> rfmt_order;
        unsigned int nrfmt;
        
        index_mapping idxmap;
        
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
        
        double posl[3];
        double posc[3];
        double bfcmt[3];
        int symclass;
        Species *species;
        unsigned int offset_apw;
        unsigned int offset_lo;
        unsigned int offset_wfmt;
        int idxclass;
};

/// global read-only variables initialized once at the beginning
struct lapw_global_variables
{
    /// maximum number of G+k vectors
    unsigned int ngkmax;
    
    /// maximum order of APW functions 
    unsigned int apwordmax;
    
    /// maximum l for APW functions
    unsigned int lmaxapw;
    
    /// (lmaxapw+1)^2
    unsigned int lmmaxapw;
    
    /// maximum l for potential and density
    unsigned int lmaxvr;
    
    /// (lmaxvr+1)^2
    unsigned int lmmaxvr;
    
    unsigned int lmaxlu;
    
    unsigned int lmmaxlu;
    
    bool ldapu;
    
    /// number of G-vectors
    unsigned int ngvec;
    
    /// FFT grid dimensions
    int ngrid[3];
    
    /// number of FFT grid points
    unsigned int ngrtot;
    
    /// maximum number of radial points
    unsigned int nrmtmax;
    
    /// number of first-variational states
    unsigned int nstfv;
    
    /// number of second-variational states
    unsigned int nstsv;
    
    /// maximum number of radial functions over all species
    unsigned int nrfmtmax;
    
    /// maximum number (order) order of radial functions over all species and l-channels 
    unsigned int ordrfmtmax;
    
    /// tolerance for generalized eigen-value solver
    double evaltol;
    
    /// spin-polarized calculation (T/F)
    bool spinpol;
    
    bool spinorb;
    
    /// number of dimensions (0,1 or 3) of the magnetization field
    unsigned int ndmag;
    
    /// number of spinor components (1 or 2)
    unsigned int nspinor;

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
    
    unsigned int size_wfmt_apw;
    unsigned int size_apwalm; 
    unsigned int size_wfmt_lo;
    unsigned int size_wfmt;
    std::vector<int> l_by_lm;
    
    unsigned int natmcls;
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

};
extern lapw_global_variables lapw_global;

/// describes the Bloch-states for a specific k-point 
class bloch_states_k
{
    public:

        bloch_states_k(unsigned int ngk);

        void copy_apwalm(complex16 *apwalm_);
        
        void generate_scalar_wave_functions(); 
        
        void generate_spinor_wave_functions(diagonalization mode);
        
        void test_scalar_wave_functions(int use_fft);
        
        void test_spinor_wave_functions(int use_fft);
        
        /// returns global index of G1-G2 vector, where G1 and G2 vectors are specified by their local indices
        inline unsigned int idxG12(unsigned int ig1, unsigned int ig2)
        {
            return lapw_global.ivgig(lapw_global.ivg(0, idxg[ig1]) - lapw_global.ivg(0, idxg[ig2]), 
                                     lapw_global.ivg(1, idxg[ig1]) - lapw_global.ivg(1, idxg[ig2]),
                                     lapw_global.ivg(2, idxg[ig1]) - lapw_global.ivg(2, idxg[ig2]));

        }

        /// number of G+k vectors
        unsigned int ngk;

        /// number of LAPW basis functions (number of plane-waves + number of local orbitals)
        unsigned int lapw_basis_size;

        /// number of wave-function expansion coefficients (number of plane-waves + number of MT lm- coefficients)
        unsigned int wave_function_size;

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

/*class lapw_eigen_states
{
    public:

        lapw_eigen_states(kpoint *kp) : kp(kp)
        {
        }
                
        void pack_apwalm(complex16 *apwalm_);
        
        void generate_scalar_wf(); 
        
        void generate_spinor_wf(diagonalization mode);
        
        void test_scalar_wf(int use_fft);
        
        void test_spinor_wf(int use_fft);
        
        kpoint *kp;
        mdarray<complex16,2> apwalm;
        mdarray<complex16,2> evecfv;
        mdarray<complex16,2> evecsv;
        mdarray<complex16,2> evecfd;
        mdarray<complex16,2> scalar_wf;
        mdarray<complex16,3> spinor_wf;
        std::vector<double> evalfv;
        std::vector<double> evalsv;
        std::vector<double> occsv;
}; */

/*class Geometry 
{
    public:
    
        std::vector< std::vector<double> > avec;
        double avec_m[9];
        double ainv_m[9];
        std::vector< std::vector<double> > bvec;
        double bvec_m[9];
        double binv_m[9];
        std::vector<Species*> species;
        std::vector<Atom> atoms;
        double omega;
};*/

/*class Parameters
{
    public:
        unsigned int ngkmax;
        unsigned int apwordmax;
        unsigned int lmmaxapw;
        unsigned int natmtot;
        unsigned int nspecies;
        unsigned int lmaxvr;
        unsigned int lmmaxvr;
        unsigned int lmaxapw;
        unsigned int ngvec;
        unsigned int ngrtot;
        unsigned int nlomax;
        unsigned int nrmtmax;
        unsigned int nstfv;
        unsigned int nstsv;
        unsigned int nmatmax;
        unsigned int nrfmtmax;
        unsigned int ordrfmtmax;
        double evaltol;
        int ngrid[3];
        bool spinpol;
        unsigned int ndmag;
        unsigned int nspinor;
        std::vector<int> igfft;
        std::vector< std::complex<double> > cfunig;
        std::vector<double> cfunir;
        mdarray<int,2> intgv;
        mdarray<int,2> ivg;
        mdarray<int,3> ivgig;
        mdarray<std::complex<double>,3> gntyry;
        mdarray<std::vector<int>,2> L3_gntyry;
        mdarray<std::vector<complex16>,2> L3_gntyry_data;
        
        /// collection of bloch states for a local fraction of k-points
        std::vector<bloch_states_k*> bloch_states;
        
        unsigned int size_wfmt_apw;
        unsigned int size_apwalm; 
        unsigned int size_wfmt_lo;
        unsigned int size_wfmt;
        mdarray<double,4> hmltrad;
        mdarray<double,4> ovlprad;
        mdarray<double,5> beffrad;
        mdarray<double,5> apwfr;
        mdarray<double,3> apwdfr;
        mdarray<double,2> beffir;
        mdarray<complex16,2> beffig;
        mdarray<complex16,1> veffig;
        std::vector<int> l_by_lm;
        
        unsigned int natmcls;
        std::vector<int> ic2ias;
        std::vector<int> natoms_in_class;
};*/


/*
    global variables
*/
//extern Geometry geometry;
//extern Parameters p;
extern complex16 zone;
extern complex16 zzero;
extern complex16 zi;
extern double y00;

/*
    inline functions
*/
/*inline int idxlm(int l, int m)
{
    return l * l + l + m;
}*/

/*inline int idxG12(const kpoint& kp, int ig1, int ig2)
{
    return p.ivgig(p.ivg(0, kp.idxg[ig1]) - p.ivg(0, kp.idxg[ig2]), 
                   p.ivg(1, kp.idxg[ig1]) - p.ivg(1, kp.idxg[ig2]),
                   p.ivg(2, kp.idxg[ig1]) - p.ivg(2, kp.idxg[ig2]));
}

inline int idxG12(const kpoint *kp, int ig1, int ig2)
{
    return p.ivgig(p.ivg(0, kp->idxg[ig1]) - p.ivg(0, kp->idxg[ig2]), 
                   p.ivg(1, kp->idxg[ig1]) - p.ivg(1, kp->idxg[ig2]),
                   p.ivg(2, kp->idxg[ig1]) - p.ivg(2, kp->idxg[ig2]));
}*/

template <typename T>
inline void L3_sum_gntyry(int lm1, int lm2, T *v, std::complex<double>& zsum)
{
    for (unsigned int k = 0; k < lapw_global.L3_gntyry(lm1, lm2).size(); k++)
    {
        int lm3 = lapw_global.L3_gntyry(lm1, lm2)[k];
        zsum += lapw_global.gntyry(lm3, lm1, lm2) * v[lm3];
    }
}
        

#endif // __LAPW_H__
