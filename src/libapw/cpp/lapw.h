#ifndef __LAPW_H__
#define __LAPW_H__

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

/*
    forward class declarations
*/
struct atomic_level;
struct radial_solution_descriptor;
class radial_l_channel_descriptor;
class mtci;
class Species;
class Atom;
class kpoint;
class Geometry;
class Parameters;
class lapw_eigen_states;

/*
    function declarations
*/

void lapw_set_sv(lapw_eigen_states& eigen_data);

void lapw_fft(int32_t direction, complex16 *data);

void lapw_density(lapw_eigen_states& eigen_states, mdarray<double,6>& densmt, mdarray<double,3>& densir);

inline int idxlm(int l, int m);

/*
    actual class definitions
*/
struct atomic_level 
{
    int n;
    int l;
    int k;
    int occupancy;
};

struct radial_solution_descriptor
{
    int n;
    int l;
    int dme;
    double enu;
    bool auto_enu;
};

class radial_l_channel_descriptor
{
    public:
        radial_l_channel_descriptor()
        {
        }
        
        radial_l_channel_descriptor(int l) : l(l)
        {
        }
    
        int l;
        std::vector<radial_solution_descriptor> radial_solution_descriptors;
};

// muffin-tin combined indices
class mtci
{
    public:
        mtci(int l, int m, int order, int idxrf) : l(l), m(m), order(order), idxrf(idxrf), idxlo(-1)
        {
            lm = idxlm(l, m);
        }
        mtci(int l, int m, int order, int idxrf, int idxlo) : l(l), m(m), order(order), idxrf(idxrf), idxlo(idxlo)
        {
            lm = idxlm(l, m);
        }
        
        int l;
        int m;
        int lm;
        int order;
        int idxrf;
        int idxlo;
};

class Species 
{
    public:

        Species() : size_ci_lo(0), size_ci_apw(0), nrfmt(0)
        {
        }
    
        Species(const std::string& symbol) : symbol(symbol) 
        {
        };
        
        std::string name;
        std::string symbol;
        int number;
        double mass;
        double rmin;
        double rmax;
        double rmt;
        int nrmt;
  
        std::vector<atomic_level> core;
        std::vector<radial_l_channel_descriptor> lo_descriptors;
        std::vector<radial_l_channel_descriptor> apw_descriptors;
        std::vector<mtci> ci;
        mtci *ci_lo;
        unsigned int size_ci_lo;
        unsigned int size_ci_apw;
        mdarray<int,2> ci_by_lmo;
        std::vector<int> ci_by_idxrf;
        std::vector<int> l_by_idxrf;
        std::vector<int> rfmt_order;
        unsigned int nrfmt;
};

class Atom 
{
    public:
        
        Atom()
        {
        }

        Atom(Species *species) : species(species)
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
};

class kpoint
{
    public:

        kpoint(unsigned int ngk) : ngk(ngk)
        {
            idxg.resize(ngk);
            idxgfft.resize(ngk);
        }
        
        unsigned int ngk;
    
        mdarray<double,2> vgkc;
        std::vector<int> idxg;
        std::vector<int> idxgfft;
        double weight;
};

class Geometry 
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
};

class Parameters
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
        std::vector<kpoint*> kpoints;
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

        //inline void L3_sum_gntyry(int lm1, int lm2, double *v, std::complex<double>& zsum)
        //{
        //    for (unsigned int k = 0; k < L3_gntyry(lm1, lm2).size(); k++)
        //    {
        //        int lm3 = L3_gntyry(lm1, lm2)[k];
        //        zsum += gntyry(lm3, lm1, lm2) * v[lm3];
        //    }
        //}
        //
        //inline void L3_sum_gntyry(int lm1, int lm2, complex16 *v, complex16& zsum)
        //{
        //    for (unsigned int k = 0; k < L3_gntyry(lm1, lm2).size(); k++)
        //    {
        //        int lm3 = L3_gntyry(lm1, lm2)[k];
        //        zsum += gntyry(lm3, lm1, lm2) * v[lm3];
        //    }
        //}
};

class lapw_eigen_states
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
        //mdarray<complex16,2> apwalm1;
        mdarray<complex16,2> evecfv;
        mdarray<complex16,2> evecsv;
        mdarray<complex16,2> evecfd;
        mdarray<complex16,2> scalar_wf;
        mdarray<complex16,3> spinor_wf;
        std::vector<double> evalfv;
        std::vector<double> evalsv;
        std::vector<double> occsv;
};

/*
    global variables
*/
extern Geometry geometry;
extern Parameters p;
extern complex16 zone;
extern complex16 zzero;
extern complex16 zi;
extern double y00;

/*
    inline functions
*/
inline int idxlm(int l, int m)
{
    return l * l + l + m;
}

inline int idxG12(const kpoint& kp, int ig1, int ig2)
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
}

template <typename T>
inline void L3_sum_gntyry(int lm1, int lm2, T *v, std::complex<double>& zsum)
{
    for (unsigned int k = 0; k < p.L3_gntyry(lm1, lm2).size(); k++)
    {
        int lm3 = p.L3_gntyry(lm1, lm2)[k];
        zsum += p.gntyry(lm3, lm1, lm2) * v[lm3];
    }
}
        

#endif // __LAPW_H__
