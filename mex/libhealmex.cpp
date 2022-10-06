#include <mex.hpp>
#include <mexAdapter.hpp>

#include <complex>
#include <utility> // tuple
#include <vector>

#include <libsharp/sharp_almhelpers.h>
#include <libsharp/sharp_geomhelpers.h>
#include <libsharp/sharp_cxx.h>

#include <healpix_map.h>
#include <alm.h>
#include <sharp_almhelpers.h>
#include <sharp_geomhelpers.h>
#include <sharp_cxx.h>
#include <alm_healpix_tools.h>
#include <alm_powspec_tools.h>
#include <error_handling.h>
#include <powspec.h>
#include <rangeset.h>
#include <rotmatrix.h>
#include <xcomplex.h>

#include "libmcm.cpp"

using namespace std;
using namespace matlab::data;
using matlab::mex::ArgumentList;

using complex64 = complex<double>;
using healpix = Healpix_Base2;
using healmap = Healpix_Map<double>;
using healalm = Alm<complex64>;

/* Useful type predicates */

inline bool ischar(Array& a)      { return a.getType() == ArrayType::CHAR; }
inline bool isbool(Array& a)      { return a.getType() == ArrayType::LOGICAL; }
inline bool isint32(Array& a)     { return a.getType() == ArrayType::INT32; }
inline bool isint64(Array& a)     { return a.getType() == ArrayType::INT64; }
inline bool isdouble(Array& a)    { return a.getType() == ArrayType::DOUBLE; }
inline bool iscomplex64(Array& a) { return a.getType() == ArrayType::COMPLEX_DOUBLE; }
inline bool isint(Array& a) {
    switch (sizeof(int)) {
        case sizeof(int32_t): return a.getType() == ArrayType::INT32; break;
        case sizeof(int64_t): return a.getType() == ArrayType::INT64; break;
        default: return false;
    }
}
inline bool islong(Array& a) {
    switch (sizeof(long)) {
        case sizeof(int32_t): return a.getType() == ArrayType::INT32; break;
        case sizeof(int64_t): return a.getType() == ArrayType::INT64; break;
        default: return false;
    }
}

inline bool isscalar(Array& a) { return a.getNumberOfElements() == 1; }
inline bool isvector(Array& a) {
    auto dims = a.getDimensions();
    if (dims.size() > 2)
        return false;
    if (dims[0] == 1 || dims[1] == 1)
        return true;
    return false;
}
inline bool ismatrix(Array& a) {
    auto dims = a.getDimensions();
    if (dims.size() > 2)
        return false;
    return true;
}

inline double angdist(pointing p1, pointing p2){
  /*double gamma = p2.phi - p1.phi; 
    return acos(cos(p1.theta)*cos(p2.theta) + sin(p1.theta)*sin(p2.theta)*cos(gamma));*/
  vec3 v1 = p1.to_vec3();
  vec3 v2 = p2.to_vec3();
  return acos(dotprod(v1,v2));
}


template <typename T, typename A> buffer_ptr_t<T> buffer(A& array);

/*
 * MEX entry point
 */

enum libhealpix_mex_calls {
    id_heartbeat        = -1,
    id_nest2ring        =  1,
    id_ring2nest        =  2,
    id_pix2ring         = 10,
    id_pix2vec          = 11,
    id_pix2zphi         = 12,
    id_pix2ang          = 13,
    id_vec2pix          = 14,
    id_zphi2pix         = 15,
    id_ang2pix          = 16,
    id_pix2xyf          = 17,
    id_xyf2pix          = 18,
    id_query_disc       = 19,
    id_neighbors        = 20,
    id_map2alm_iter     = 53,
    id_map2alm_spin_iter= 54,
    id_alm2map          = 55,
    id_alm2map_spin     = 56,
    id_alm2map_der1     = 57,
    id_map2alm_pure     = 58,
    id_alm2cl           = 61,
    id_almxfl           = 62,
    id_rotate_alm_coord = 65,
    id_rotate_alm_euler = 66,
    id_rotate_alm_matrix= 67,
    id_apodize_mask     = 68,
    id_shrink_mask      = 69,
    id_smooth_mask      = 70,
    id_alm2map_polonly  = 71,
    id_smoothing_pol    = 72,
    id_smoothing        = 73,
    id_scan_rings_observed = 74,
    id_compute_mcm      = 75,
    id_map2alm_polonly_iter = 76
};

class MexFunction : public matlab::mex::Function {
public:
    MexFunction()
    {
        // Sanity check of complex data representation. This assumption must
        // hold for us to pass data on to the HEALPix routines.
        //  - Generate complex doubles in a Matlab array
        //  - Extract a pointer to the data buffer backing the Matlab array
        //  - Reinterpret cast the complex pointer as a double pointer
        //  - Check that the values match the expected packed real+imag
        //    layout
        auto ml_arr = factory.createArray({2, 1}, {
                complex<double>{1, 2}, complex<double>{-3, -4}});
        auto ml_buf = buffer<complex<double>>(ml_arr);
        double* packed = reinterpret_cast<double*>(ml_buf.get());
        if (packed[0] != 1 || packed[1] != 2 || packed[2] != -3 || packed[3] != -4) {
            error("Unexpected complex format!");
        }
    }

    void operator ()(ArgumentList outputs, ArgumentList inputs) {
        if (inputs.size() < 1) {
            error("Missing dispatch command.");
        }
        auto ml_dispatch = inputs[0];
        if (!isint64(ml_dispatch) || !isscalar(ml_dispatch)) {
            error("Dispatch argument must be a scalar of type int64.");
        }
        int64_t dispatch = ml_dispatch[0];

        auto n_in  = inputs.size() - 1;
        auto n_out = outputs.size();

        #define CHECK_WRAP(expr) \
            do { \
                expr \
            } while (false)
        #define CHECK_NINOUT(funcname, nin, nout) CHECK_WRAP( \
            if (n_in != (nin)) { \
                error(funcname ": Expected " #nin " inputs, got %d.", n_in); \
            } \
            if (n_out != (nout)) { \
                error(funcname ": Expected " #nout " outputs, got %d.", n_out); \
            } )
        #define CHECK_INPUT_SCALAR(funcname, argname, num) CHECK_WRAP( \
            if (!isscalar(inputs[num])) { \
                error(funcname ": argument " #num " (" argname ") must be a scalar.", num); \
            } )
        #define CHECK_INPUT_VECTOR(funcname, argname, num) CHECK_WRAP( \
            if (!isvector(inputs[num])) { \
                error(funcname ": argument " #num " (" argname ") must be a vector.", num); \
            } )
        #define CHECK_INPUT_MATRIX(funcname, argname, num) CHECK_WRAP( \
            if (!ismatrix(inputs[num])) { \
                error(funcname ": argument " #num " (" argname ") must be a matrix.", num); \
            } )
        #define CHECK_INPUT_BOOL(funcname, argname, num) CHECK_WRAP( \
            if (!isbool(inputs[num])) { \
                error(funcname ": argument " #num " (" argname ") must of type int.", num); \
            } )
        #define CHECK_INPUT_INT(funcname, argname, num) CHECK_WRAP( \
            if (!isint(inputs[num])) { \
                error(funcname ": argument " #num " (" argname ") must of type int.", num); \
            } )
        #define CHECK_INPUT_INT32(funcname, argname, num) CHECK_WRAP( \
            if (!isint32(inputs[num])) { \
                error(funcname ": argument " #num " (" argname ") must of type int32.", num); \
            } )
        #define CHECK_INPUT_INT64(funcname, argname, num) CHECK_WRAP( \
            if (!isint64(inputs[num])) { \
                error(funcname ": argument " #num " (" argname ") must of type int64.", num); \
            } )
        #define CHECK_INPUT_DOUBLE(funcname, argname, num) CHECK_WRAP( \
            if (!isdouble(inputs[num])) { \
                error(funcname ": argument " #num " (" argname ") must of type double.", num); \
            } )
        #define CHECK_INPUT_COMPLEX64(funcname, argname, num) CHECK_WRAP( \
            if (!iscomplex64(inputs[num])) { \
                error(funcname ": argument " #num " (" argname ") must of type complex double.", num); \
            } )
        #define CHECK_INPUT_CHAR(funcname, argname, num) CHECK_WRAP( \
            if (!ischar(inputs[num])) { \
                error(funcname ": argument " #num " (" argname ") must of type char.", num); \
            } )

        switch (dispatch) {
            case id_heartbeat:
                outputs[0] = factory.createScalar(true);
                break;

            case id_nest2ring:
                CHECK_NINOUT("nest2ring", 2, 1);
                CHECK_INPUT_SCALAR("nest2ring", "nside", 1);
                CHECK_INPUT_INT64("nest2ring", "nside", 1);
                CHECK_INPUT_INT64("nest2ring", "ipix", 2);
                mex_nest2ring(outputs, inputs);
                break;

            case id_ring2nest:
                CHECK_NINOUT("ring2nest", 2, 1);
                CHECK_INPUT_SCALAR("ring2nest", "nside", 1);
                CHECK_INPUT_INT64("ring2nest", "nside", 1);
                CHECK_INPUT_INT64("ring2nest", "ipix", 2);
                mex_ring2nest(outputs, inputs);
                break;

            case id_pix2ring:
                CHECK_NINOUT("pix2ring", 3, 3);
                CHECK_INPUT_SCALAR("pix2ring", "nside", 1);
                CHECK_INPUT_INT64("pix2ring", "nside", 1);
                CHECK_INPUT_SCALAR("pix2ring", "nest", 2);
                CHECK_INPUT_BOOL("pix2ring", "nest", 2);
                CHECK_INPUT_INT64("pix2ring", "ipix", 3);
                mex_pix2ring(outputs, inputs);
                break;

            case id_pix2vec:
                CHECK_NINOUT("pix2vec", 3, 3);
                CHECK_INPUT_SCALAR("pix2vec", "nside", 1);
                CHECK_INPUT_INT64("pix2vec", "nside", 1);
                CHECK_INPUT_SCALAR("pix2vec", "nest", 2);
                CHECK_INPUT_BOOL("pix2vec", "nest", 2);
                CHECK_INPUT_INT64("pix2vec", "ipix", 3);
                mex_pix2vec(outputs, inputs);
                break;

            case id_pix2zphi:
                CHECK_NINOUT("pix2zphi", 3, 2);
                CHECK_INPUT_SCALAR("pix2zphi", "nside", 1);
                CHECK_INPUT_INT64("pix2zphi", "nside", 1);
                CHECK_INPUT_SCALAR("pix2zphi", "nest", 2);
                CHECK_INPUT_BOOL("pix2zphi", "nest", 2);
                CHECK_INPUT_INT64("pix2zphi", "ipix", 3);
                mex_pix2zphi(outputs, inputs);
                break;

            case id_pix2ang:
                CHECK_NINOUT("pix2ang", 3, 2);
                CHECK_INPUT_SCALAR("pix2ang", "nside", 1);
                CHECK_INPUT_INT64("pix2ang", "nside", 1);
                CHECK_INPUT_SCALAR("pix2ang", "nest", 2);
                CHECK_INPUT_BOOL("pix2ang", "nest", 2);
                CHECK_INPUT_INT64("pix2ang", "ipix", 3);
                mex_pix2ang(outputs, inputs);
                break;

            case id_vec2pix:
                CHECK_NINOUT("vec2pix", 5, 1);
                CHECK_INPUT_SCALAR("vec2pix", "nside", 1);
                CHECK_INPUT_INT64("vec2pix", "nside", 1);
                CHECK_INPUT_SCALAR("vec2pix", "nest", 2);
                CHECK_INPUT_BOOL("vec2pix", "nest", 2);
                CHECK_INPUT_DOUBLE("vec2pix", "x", 3);
                CHECK_INPUT_DOUBLE("vec2pix", "y", 3);
                CHECK_INPUT_DOUBLE("vec2pix", "z", 3);
                mex_vec2pix(outputs, inputs);
                break;

            case id_zphi2pix:
                CHECK_NINOUT("zphi2pix", 4, 1);
                CHECK_INPUT_SCALAR("zphi2pix", "nside", 1);
                CHECK_INPUT_INT64("zphi2pix", "nside", 1);
                CHECK_INPUT_SCALAR("zphi2pix", "nest", 2);
                CHECK_INPUT_BOOL("zphi2pix", "nest", 2);
                CHECK_INPUT_DOUBLE("zphi2pix", "z", 3);
                CHECK_INPUT_DOUBLE("zphi2pix", "phi", 4);
                mex_zphi2pix(outputs, inputs);
                break;

            case id_ang2pix:
                CHECK_NINOUT("ang2pix", 4, 1);
                CHECK_INPUT_SCALAR("ang2pix", "nside", 1);
                CHECK_INPUT_INT64("ang2pix", "nside", 1);
                CHECK_INPUT_SCALAR("ang2pix", "nest", 2);
                CHECK_INPUT_BOOL("ang2pix", "nest", 2);
                CHECK_INPUT_DOUBLE("ang2pix", "theta", 3);
                CHECK_INPUT_DOUBLE("ang2pix", "phi", 4);
                mex_ang2pix(outputs, inputs);
                break;

            case id_pix2xyf:
                CHECK_NINOUT("pix2xyf", 3, 3);
                CHECK_INPUT_SCALAR("pix2xyf", "nside", 1);
                CHECK_INPUT_INT64("pix2xyf", "nside", 1);
                CHECK_INPUT_SCALAR("pix2xyf", "nest", 2);
                CHECK_INPUT_BOOL("pix2xyf", "nest", 2);
                CHECK_INPUT_INT64("pix2xyf", "ipix", 3);
                mex_pix2xyf(outputs, inputs);
                break;

            case id_xyf2pix:
                CHECK_NINOUT("xyf2pix", 5, 1);
                CHECK_INPUT_SCALAR("xyf2pix", "nside", 1);
                CHECK_INPUT_INT64("xyf2pix", "nside", 1);
                CHECK_INPUT_SCALAR("xyf2pix", "nest", 2);
                CHECK_INPUT_BOOL("xyf2pix", "nest", 2);
                CHECK_INPUT_INT("xyf2pix", "x", 3);
                CHECK_INPUT_INT("xyf2pix", "y", 4);
                CHECK_INPUT_INT("xyf2pix", "f", 5);
                mex_xyf2pix(outputs, inputs);
                break;

            case id_query_disc:
                CHECK_NINOUT("query_disc", 5, 1);
                CHECK_INPUT_SCALAR("query_disc", "nside", 1);
                CHECK_INPUT_INT64("query_disc", "nside", 1);
                CHECK_INPUT_SCALAR("query_disc", "nest", 2);
                CHECK_INPUT_BOOL("query_disc", "nest", 2);
                CHECK_INPUT_VECTOR("query_disc", "rvec", 3);
                CHECK_INPUT_DOUBLE("query_disc", "rvec", 3);
                CHECK_INPUT_SCALAR("query_disc", "radius", 4);
                CHECK_INPUT_DOUBLE("query_disc", "radius", 4);
                CHECK_INPUT_SCALAR("query_disc", "inclusive", 5);
                CHECK_INPUT_BOOL("query_disc", "inclusive", 5);
                mex_query_disc(outputs, inputs);
                break;

            case id_neighbors:
                CHECK_NINOUT("neighbors", 3, 1);
                CHECK_INPUT_SCALAR("neighbors", "nside", 1);
                CHECK_INPUT_INT64("neighbors", "nside", 1);
                CHECK_INPUT_SCALAR("neighbors", "nest", 2);
                CHECK_INPUT_BOOL("neighbors", "nest", 2);
                CHECK_INPUT_VECTOR("neighbors", "ipix", 3);
                CHECK_INPUT_INT64("neighbors", "ipix", 3);
                mex_neighbors(outputs, inputs);
                break;

            case id_map2alm_iter:
                CHECK_NINOUT("map2alm_iter", 8, 3);
                CHECK_INPUT_SCALAR("map2alm_iter", "nside", 1);
                CHECK_INPUT_INT64("map2alm_iter", "nside", 1);
                CHECK_INPUT_DOUBLE("map2alm_iter", "mapT", 2);
                CHECK_INPUT_DOUBLE("map2alm_iter", "mapQ", 3);
                CHECK_INPUT_DOUBLE("map2alm_iter", "mapU", 4);
                CHECK_INPUT_SCALAR("map2alm_iter", "lmax", 5);
                CHECK_INPUT_INT32("map2alm_iter", "lmax", 5);
                CHECK_INPUT_SCALAR("map2alm_iter", "mmax", 6);
                CHECK_INPUT_INT32("map2alm_iter", "mmax", 6);
                CHECK_INPUT_DOUBLE("map2alm_iter", "rwghts", 7);
                CHECK_INPUT_SCALAR("map2alm_iter", "iter", 8);
                CHECK_INPUT_INT32("map2alm_iter", "iter", 8);
                mex_map2alm_iter(outputs, inputs);
                break;

            case id_map2alm_spin_iter:
                CHECK_NINOUT("map2alm_spin_iter", 8, 2);
                CHECK_INPUT_SCALAR("map2alm_spin_iter", "nside", 1);
                CHECK_INPUT_INT64("map2alm_spin_iter", "nside", 1);
                CHECK_INPUT_DOUBLE("map2alm_spin_iter", "map1", 2);
                CHECK_INPUT_DOUBLE("map2alm_spin_iter", "map2", 3);
                CHECK_INPUT_SCALAR("map2alm_spin_iter", "spin", 4);
                CHECK_INPUT_INT32("map2alm_spin_iter", "spin", 4);
                CHECK_INPUT_SCALAR("map2alm_spin_iter", "lmax", 5);
                CHECK_INPUT_INT32("map2alm_spin_iter", "lmax", 5);
                CHECK_INPUT_SCALAR("map2alm_spin_iter", "mmax", 6);
                CHECK_INPUT_INT32("map2alm_spin_iter", "mmax", 6);
                CHECK_INPUT_DOUBLE("map2alm_spin_iter", "rwghts", 7);
                CHECK_INPUT_SCALAR("map2alm_spin_iter", "iter", 8);
                CHECK_INPUT_INT32("map2alm_spin_iter", "iter", 8);
                mex_map2alm_spin_iter(outputs, inputs);
                break;

            case id_map2alm_polonly_iter:
                CHECK_NINOUT("map2alm_polonly_iter", 8, 2);
                CHECK_INPUT_SCALAR("map2alm_polonly_iter", "nside", 1);
                CHECK_INPUT_INT64("map2alm_polonly_iter", "nside", 1);
                CHECK_INPUT_CHAR("map2alm_polonly_iter", "order", 2);
                CHECK_INPUT_DOUBLE("map2alm_polonly_iter", "mapQ", 3);
                CHECK_INPUT_DOUBLE("map2alm_polonly_iter", "mapU", 4);
                CHECK_INPUT_SCALAR("map2alm_polonly_iter", "lmax", 5);
                CHECK_INPUT_INT32("map2alm_polonly_iter", "lmax", 5);
                CHECK_INPUT_SCALAR("map2alm_polonly_iter", "mmax", 6);
                CHECK_INPUT_INT32("map2alm_polonly_iter", "mmax", 6);
                CHECK_INPUT_DOUBLE("map2alm_polonly_iter", "rwghts", 7);
                CHECK_INPUT_SCALAR("map2alm_polonly_iter", "iter", 8);
                CHECK_INPUT_INT32("map2alm_polonly_iter", "iter", 8);
                mex_map2alm_polonly_iter(outputs, inputs);
                break;
				
            case id_map2alm_pure:
                CHECK_NINOUT("map2alm_pure", 10, 2);
                CHECK_INPUT_SCALAR("map2alm_pure", "nside", 1);
                CHECK_INPUT_INT64("map2alm_pure", "nside", 1);
                CHECK_INPUT_CHAR("map2alm_pure", "order", 2);
                CHECK_INPUT_DOUBLE("map2alm_pure", "mapQ", 3);
                CHECK_INPUT_DOUBLE("map2alm_pure", "mapU", 4);
                CHECK_INPUT_DOUBLE("map2alm_pure", "mapW", 5);
                CHECK_INPUT_SCALAR("map2alm_pure", "lmax", 6);
                CHECK_INPUT_INT32("map2alm_pure", "lmax", 6);
                CHECK_INPUT_SCALAR("map2alm_pure", "mmax", 7);
                CHECK_INPUT_INT32("map2alm_pure", "mmax", 7);
                CHECK_INPUT_DOUBLE("map2alm_pure", "rwghts", 8);
                CHECK_INPUT_SCALAR("map2alm_pure", "iter", 9);
                CHECK_INPUT_INT32("map2alm_pure", "iter", 9);
                mex_map2alm_pure(outputs, inputs);
                break;

            case id_alm2map:
                CHECK_NINOUT("alm2map", 6, 3);
                CHECK_INPUT_SCALAR("alm2map", "lmax", 1);
                CHECK_INPUT_INT32("alm2map", "lmax", 1);
                CHECK_INPUT_SCALAR("alm2map", "mmax", 2);
                CHECK_INPUT_INT32("alm2map", "mmax", 2);
                CHECK_INPUT_COMPLEX64("alm2map", "almsT", 3);
                CHECK_INPUT_COMPLEX64("alm2map", "almsG", 4);
                CHECK_INPUT_COMPLEX64("alm2map", "almsC", 5);
                CHECK_INPUT_SCALAR("alm2map", "nside", 6);
                CHECK_INPUT_INT64("alm2map", "nside", 6);
                mex_alm2map(outputs, inputs);
                break;

            case id_alm2map_polonly:
                CHECK_NINOUT("alm2map_polonly", 7, 2);
                CHECK_INPUT_SCALAR("alm2map_polonly", "lmax", 1);
                CHECK_INPUT_INT32("alm2map_polonly", "lmax", 1);
                CHECK_INPUT_SCALAR("alm2map_polonly", "mmax", 2);
                CHECK_INPUT_INT32("alm2map_polonly", "mmax", 2);
                CHECK_INPUT_COMPLEX64("alm2map_polonly", "almsG", 3);
                CHECK_INPUT_COMPLEX64("alm2map_polonly", "almsC", 4);
                CHECK_INPUT_SCALAR("alm2map_polonly", "nside", 5);
                CHECK_INPUT_INT64("alm2map_polonly", "nside", 5);
                CHECK_INPUT_CHAR("alm2map_polonly", "order", 6);
                CHECK_INPUT_DOUBLE("alm2map_polonly", "rwghts", 7);
                mex_alm2map_polonly(outputs, inputs);
                break;

            case id_alm2map_spin:
                CHECK_NINOUT("alm2map_spin", 6, 2);
                CHECK_INPUT_SCALAR("alm2map_spin", "lmax", 1);
                CHECK_INPUT_INT32("alm2map_spin", "lmax", 1);
                CHECK_INPUT_SCALAR("alm2map_spin", "mmax", 2);
                CHECK_INPUT_INT32("alm2map_spin", "mmax", 2);
                CHECK_INPUT_COMPLEX64("alm2map_spin", "alms1", 3);
                CHECK_INPUT_COMPLEX64("alm2map_spin", "alms2", 4);
                CHECK_INPUT_SCALAR("alm2map_spin", "spin", 5);
                CHECK_INPUT_INT32("alm2map_spin", "spin", 5);
                CHECK_INPUT_SCALAR("alm2map_spin", "nside", 6);
                CHECK_INPUT_INT64("alm2map_spin", "nside", 6);
                mex_alm2map_spin(outputs, inputs);
                break;

            case id_alm2map_der1:
                CHECK_NINOUT("alm2map_der1", 4, 3);
                CHECK_INPUT_SCALAR("alm2map_der1", "lmax", 1);
                CHECK_INPUT_INT32("alm2map_der1", "lmax", 1);
                CHECK_INPUT_SCALAR("alm2map_der1", "mmax", 2);
                CHECK_INPUT_INT32("alm2map_der1", "mmax", 2);
                CHECK_INPUT_COMPLEX64("alm2map_der1", "alms", 3);
                CHECK_INPUT_SCALAR("alm2map_der1", "nside", 4);
                CHECK_INPUT_INT64("alm2map_der1", "nside", 4);
                mex_alm2map_der1(outputs, inputs);
                break;

            case id_alm2cl:
                CHECK_NINOUT("alm2cl", 4, 1);
                CHECK_INPUT_SCALAR("alm2cl", "lmax", 1);
                CHECK_INPUT_INT32("alm2cl", "lmax", 1);
                CHECK_INPUT_SCALAR("alm2cl", "mmax", 2);
                CHECK_INPUT_INT32("alm2cl", "mmax", 2);
                CHECK_INPUT_COMPLEX64("alm2cl", "alms1", 3);
                CHECK_INPUT_COMPLEX64("alm2cl", "alms2", 4);
                mex_alm2cl(outputs, inputs);
                break;

            case id_almxfl:
                CHECK_NINOUT("almxfl", 4, 1);
                CHECK_INPUT_SCALAR("almxfl", "lmax", 1);
                CHECK_INPUT_INT32("almxfl", "lmax", 1);
                CHECK_INPUT_SCALAR("almxfl", "mmax", 2);
                CHECK_INPUT_INT32("almxfl", "mmax", 2);
                CHECK_INPUT_COMPLEX64("almxfl", "alms", 3);
                CHECK_INPUT_VECTOR("almxfl", "alms", 3);
                CHECK_INPUT_COMPLEX64("almxfl", "fl", 4);
                CHECK_INPUT_VECTOR("almxfl", "fl", 4);
                mex_almxfl(outputs, inputs);
                break;

            case id_smoothing_pol:
                CHECK_NINOUT("smoothing_pol", 11, 2);
                CHECK_INPUT_SCALAR("smoothing_pol", "nside", 1);
                CHECK_INPUT_INT64("smoothing_pol", "nside", 1);
                CHECK_INPUT_CHAR("smoothing_pol", "order", 2);
                CHECK_INPUT_DOUBLE("smoothing_pol", "mapQ", 3);
                CHECK_INPUT_DOUBLE("smoothing_pol", "mapU", 4);
                CHECK_INPUT_DOUBLE("smoothing_pol", "fle", 5);
                CHECK_INPUT_VECTOR("smoothing_pol", "fle", 5);
                CHECK_INPUT_DOUBLE("smoothing_pol", "flb", 6);
                CHECK_INPUT_VECTOR("smoothing_pol", "flb", 6);
                CHECK_INPUT_SCALAR("smoothing_pol", "lmax", 7);
                CHECK_INPUT_INT32("smoothing_pol", "lmax", 7);
                CHECK_INPUT_SCALAR("smoothing_pol", "mmax", 8);
                CHECK_INPUT_INT32("smoothing_pol", "mmax", 8);
                CHECK_INPUT_SCALAR("smoothing_pol", "mmin", 9);
                CHECK_INPUT_INT32("smoothing_pol", "mmin", 9);
                CHECK_INPUT_DOUBLE("smoothing_pol", "rwghts", 10);
                CHECK_INPUT_SCALAR("smoothing_pol", "iter", 11);
                CHECK_INPUT_INT32("smoothing_pol", "iter", 11);
                mex_smoothing_pol(outputs, inputs);
                break;

            case id_smoothing:
                CHECK_NINOUT("smoothing", 9, 1);
                CHECK_INPUT_SCALAR("smoothing", "nside", 1);
                CHECK_INPUT_INT64("smoothing", "nside", 1);
                CHECK_INPUT_CHAR("smoothing", "order", 2);
                CHECK_INPUT_DOUBLE("smoothing", "map", 3);
                CHECK_INPUT_DOUBLE("smoothing", "fl", 4);
                CHECK_INPUT_VECTOR("smoothing", "fl", 4);
                CHECK_INPUT_SCALAR("smoothing", "lmax", 5);
                CHECK_INPUT_INT32("smoothing", "lmax", 5);
                CHECK_INPUT_SCALAR("smoothing", "mmax", 6);
                CHECK_INPUT_INT32("smoothing", "mmax", 6);
                CHECK_INPUT_SCALAR("smoothing", "mmin", 7);
                CHECK_INPUT_INT32("smoothing", "mmin", 7);
                CHECK_INPUT_DOUBLE("smoothing", "rwghts", 8);
                CHECK_INPUT_SCALAR("smoothing", "iter", 9);
                CHECK_INPUT_INT32("smoothing", "iter", 9);
                mex_smoothing(outputs, inputs);
                break;

            case id_rotate_alm_coord:
                CHECK_NINOUT("rotate_alm_coord", 6, 3);
                CHECK_INPUT_SCALAR("rotate_alm_coord", "itransform", 1);
                CHECK_INPUT_INT32("rotate_alm_coord", "itransform", 1);
                CHECK_INPUT_SCALAR("rotate_alm_coord", "lmax", 2);
                CHECK_INPUT_INT32("rotate_alm_coord", "lmax", 2);
                CHECK_INPUT_SCALAR("rotate_alm_coord", "mmax", 3);
                CHECK_INPUT_INT32("rotate_alm_coord", "mmax", 3);
                CHECK_INPUT_COMPLEX64("rotate_alm_coord", "almsT", 4);
                CHECK_INPUT_COMPLEX64("rotate_alm_coord", "almsG", 5);
                CHECK_INPUT_COMPLEX64("rotate_alm_coord", "almsC", 6);
                mex_rotate_alm_coord(outputs, inputs);
                break;

            case id_rotate_alm_euler:
                CHECK_NINOUT("rotate_alm_euler", 6, 3);
                CHECK_INPUT_VECTOR("rotate_alm_euler", "euler", 1);
                CHECK_INPUT_DOUBLE("rotate_alm_euler", "euler", 1);
                CHECK_INPUT_SCALAR("rotate_alm_euler", "lmax", 2);
                CHECK_INPUT_INT32("rotate_alm_euler", "lmax", 2);
                CHECK_INPUT_SCALAR("rotate_alm_euler", "mmax", 3);
                CHECK_INPUT_INT32("rotate_alm_euler", "mmax", 3);
                CHECK_INPUT_COMPLEX64("rotate_alm_euler", "almsT", 4);
                CHECK_INPUT_COMPLEX64("rotate_alm_euler", "almsG", 5);
                CHECK_INPUT_COMPLEX64("rotate_alm_euler", "almsC", 6);
                mex_rotate_alm_euler(outputs, inputs);
                break;

            case id_rotate_alm_matrix:
                CHECK_NINOUT("rotate_alm_matrix", 6, 3);
                CHECK_INPUT_MATRIX("rotate_alm_matrix", "matrix", 1);
                CHECK_INPUT_DOUBLE("rotate_alm_matrix", "matrix", 1);
                CHECK_INPUT_SCALAR("rotate_alm_matrix", "lmax", 2);
                CHECK_INPUT_INT32("rotate_alm_matrix", "lmax", 2);
                CHECK_INPUT_SCALAR("rotate_alm_matrix", "mmax", 3);
                CHECK_INPUT_INT32("rotate_alm_matrix", "mmax", 3);
                CHECK_INPUT_COMPLEX64("rotate_alm_matrix", "almsT", 4);
                CHECK_INPUT_COMPLEX64("rotate_alm_matrix", "almsG", 5);
                CHECK_INPUT_COMPLEX64("rotate_alm_matrix", "almsC", 6);
                mex_rotate_alm_matrix(outputs, inputs);
                break;

            case id_apodize_mask:
                CHECK_NINOUT("apodize_mask", 4, 1);
                CHECK_INPUT_SCALAR("apodize_mask", "nside", 1);
                CHECK_INPUT_INT64("apodize_mask", "nside", 1);
                CHECK_INPUT_CHAR("apodize_mask", "order", 2);
                CHECK_INPUT_DOUBLE("apodize_mask", "map", 3);
                CHECK_INPUT_SCALAR("apodize_mask", "radius", 4);
                CHECK_INPUT_DOUBLE("apodize_mask", "radius", 4);
                mex_apodize_mask(outputs, inputs);
                break;

            case id_shrink_mask:
                CHECK_NINOUT("shrink_mask", 4, 1);
                CHECK_INPUT_SCALAR("shrink_mask", "nside", 1);
                CHECK_INPUT_INT64("shrink_mask", "nside", 1);
                CHECK_INPUT_CHAR("shrink_mask", "order", 2);
                CHECK_INPUT_DOUBLE("shrink_mask", "map", 3);
                CHECK_INPUT_SCALAR("shrink_mask", "radius", 4);
                CHECK_INPUT_DOUBLE("shrink_mask", "radius", 4);
                mex_shrink_mask(outputs, inputs);
                break;

            case id_smooth_mask:
                CHECK_NINOUT("smooth_mask", 5, 1);
                CHECK_INPUT_SCALAR("smooth_mask", "nside", 1);
                CHECK_INPUT_INT64("smooth_mask", "nside", 1);
                CHECK_INPUT_CHAR("smooth_mask", "order", 2);
                CHECK_INPUT_DOUBLE("smooth_mask", "map", 3);
                CHECK_INPUT_SCALAR("smooth_mask", "radius", 4);
                CHECK_INPUT_DOUBLE("smooth_mask", "radius", 4);
                CHECK_INPUT_DOUBLE("smooth_mask", "rwghts", 5);
                mex_smooth_mask(outputs, inputs);
				break;
            
			case id_scan_rings_observed:
				CHECK_NINOUT("scan_rings_observed", 1, 1);
				CHECK_INPUT_DOUBLE("scan_rings_observed", "map", 1);
				CHECK_INPUT_VECTOR("scan_rings_observed", "map", 1);
				mex_scan_rings_observed(outputs, inputs);
				break;
	
			case id_compute_mcm:
                CHECK_NINOUT("compute_mcm", 11, 1);
                CHECK_INPUT_SCALAR("compute_mcm", "lmax", 1);
                CHECK_INPUT_INT32("compute_mcm", "lmax", 1);
                CHECK_INPUT_SCALAR("compute_mcm", "lmax_mask", 2);
                CHECK_INPUT_INT32("compute_mcm", "lmax_mask", 2);
                CHECK_INPUT_SCALAR("compute_mcm", "nside", 3);
                CHECK_INPUT_INT64("compute_mcm", "nside", 3);
                CHECK_INPUT_CHAR("compute_mcm", "order", 4);
                CHECK_INPUT_DOUBLE("compute_mcm", "mask", 5);
                CHECK_INPUT_DOUBLE("compute_mcm", "rwghts", 6);
                CHECK_INPUT_SCALAR("compute_mcm", "iter", 7);
                CHECK_INPUT_INT32("compute_mcm", "iter", 7);
                CHECK_INPUT_SCALAR("compute_mcm", "pe1", 8);
                CHECK_INPUT_INT32("compute_mcm", "pe1", 8);
                CHECK_INPUT_SCALAR("compute_mcm", "pe2", 9);
                CHECK_INPUT_INT32("compute_mcm", "pe2", 9);
                CHECK_INPUT_SCALAR("compute_mcm", "pb1", 10);
                CHECK_INPUT_INT32("compute_mcm", "pb1", 10);
                CHECK_INPUT_SCALAR("compute_mcm", "pb2", 11);
                CHECK_INPUT_INT32("compute_mcm", "pb2", 11);
                mex_compute_mcm(outputs, inputs);
				break;
				
            default:
                error("Unhandled dispatch type %d", dispatch);
        }
    }

private:
    std::shared_ptr<matlab::engine::MATLABEngine> engine = getEngine();
    ArrayFactory factory;

    #define DISPATCH_FN(name) \
        void mex_##name(ArgumentList& outputs, ArgumentList& inputs)

    DISPATCH_FN(nest2ring);
    DISPATCH_FN(ring2nest);

    DISPATCH_FN(pix2ring);
    DISPATCH_FN(pix2vec);
    DISPATCH_FN(pix2zphi);
    DISPATCH_FN(pix2ang);
    DISPATCH_FN(vec2pix);
    DISPATCH_FN(zphi2pix);
    DISPATCH_FN(ang2pix);
    DISPATCH_FN(pix2xyf);
    DISPATCH_FN(xyf2pix);
    DISPATCH_FN(query_disc);
    DISPATCH_FN(neighbors);

    DISPATCH_FN(map2alm_iter);
    DISPATCH_FN(map2alm_polonly_iter);
    DISPATCH_FN(map2alm_pure);
    DISPATCH_FN(alm2map);
    DISPATCH_FN(alm2map_polonly);

    DISPATCH_FN(alm2cl);
    DISPATCH_FN(almxfl);
    DISPATCH_FN(apodize_mask);
    DISPATCH_FN(shrink_mask);
    DISPATCH_FN(smooth_mask);
    DISPATCH_FN(smoothing_pol);
    DISPATCH_FN(smoothing);

    /* Utility functions */
    DISPATCH_FN(scan_rings_observed);
    DISPATCH_FN(compute_mcm);
    DISPATCH_FN(map2alm_spin_iter);
    DISPATCH_FN(alm2map);
    DISPATCH_FN(alm2map_spin);
    DISPATCH_FN(alm2map_der1);

    DISPATCH_FN(alm2cl);
    DISPATCH_FN(almxfl);
    DISPATCH_FN(rotate_alm_coord);
    DISPATCH_FN(rotate_alm_euler);
    DISPATCH_FN(rotate_alm_matrix);

    #undef DISPATCH_FN

    /*
     * Convenient wrappers to translate variadic arguments into a
     * vector<Array> list required by the Matlab API. i.e.
     *
     *   auto ml_args = scalar_args(args...);
     */

    // One argument base case
    template <typename T1>
    inline vector<Array>&& _scalar_args(vector<Array>&& args, T1 arg1) {
        args.emplace_back(factory.createScalar(arg1));
        return std::move(args);
    }
    // Multiple argument (recursive) case
    template <typename T1, typename... TN>
    inline vector<Array>&& _scalar_args(vector<Array>&& args, T1 arg1, TN... argN) {
        args.emplace_back(factory.createScalar(arg1));
        return std::move(_scalar_args(std::move(args), argN...));
    }

    template <typename... TN>
    inline vector<Array> scalar_args(TN... argN) {
        return std::move(_scalar_args(vector<Array>(), argN...));
    }

    /* call's Matlab's error(msg, ...) with [optional] arguments */
    template <typename... TN>
    void error(const string& msg, TN... argN) {
        auto&& args = scalar_args("libhealpix:error", msg, argN...);
        engine->feval(u"error", 0, args);
    }

};

/* Semi-generic helper functions */

// Extracts an explicitly-typed memory buffer from a Matlab array
template <typename T, typename A>
buffer_ptr_t<T> buffer(A& array) {
    TypedArray<T> array_t = array;
    return array_t.release();
}

template <typename T, typename A>
tuple<buffer_ptr_t<T>, size_t> bufferlen(A& array) {
    TypedArray<T> array_t = array;
    size_t len = array.getNumberOfElements();
    return tuple{array_t.release(), len};
}

// Extracts an explicitly-typed scalar from a Matlab array
template <typename T, typename A>
T scalar(A array) {
    TypedArray<T> array_t = array;
    return array_t[0];
}

// Turns a Matlab arrays giving the Nside into a HEALPix object, assuming
// RING ordering.
inline
healpix nsideorder(const TypedArray<int64_t> ml_nside)
{
    int64_t nside = ml_nside[0];
    return healpix(nside, RING, SET_NSIDE);
}

// Turns a pair of Matlab arrays giving the Nside and ordering scheme into
// a HEALPix object.
inline
healpix nsideorder(const TypedArray<int64_t> ml_nside,
                   const TypedArray<bool> ml_nested)
{
    int64_t nside = ml_nside[0];
    auto order = ml_nested[0] ? NEST : RING;
    return healpix(nside, order, SET_NSIDE);
}

// Computes the length of the complete alm vector for a given lmax & mmax
inline
int64_t alm_getn(int64_t lmax, int64_t mmax)
{
    if (lmax < 0)
        return 0;
    if (mmax < 0 || lmax < mmax)
        return 0;
    return ((mmax + 1) * (mmax + 2)) / 2 + (mmax + 1) * (lmax - mmax);
}

// Scans a map vector to determine observed rings and indicates as observed
// or not as boolean mask in output return array ring.
//   - map must be a full-sky map vector (i.e. length 12*nside*nside)
//   - ring must be a vector of length (4*nside-1)
void _scan_rings_observed(int64_t nside, const double* map, bool* ring) {
    const int64_t npix  = 12 * nside * nside;
    const int64_t nring = 4 * nside - 1;
    // Pixel offset to beginning of each ring (in northern hemisphere)
    ssize_t base = 0;
    // N.B.: Ring numbering traditionally 1-based in HEALPix
    for (ssize_t rr = 1; rr <= 2*nside; ++rr) {
        bool n_occupied = false;
        bool s_occupied = false;

        // Polar cap rings are the first (and last) nside rings.
        // Equitorial belt is (2*nside-1) rings; the equator is considered
        //   part of the northern hemisphere's nside rings, leaving (nside-1)
        //   for the southern hemisphere.

        // Within polar cap, each ring is 4*rr long; equitorial belt is 4*nside.
        ssize_t ringlen = (rr <= nside) ? 4*rr : 4*nside;
        for (ssize_t pp = 0; pp < ringlen; ++pp) {
            double nval = map[base + pp];        // Northern hemisphere
            n_occupied |= nval != 0 && !isnan(nval);

            double sval = map[npix-1-base - pp]; // Southern hemisphere
            s_occupied |= sval != 0 && !isnan(sval);
        }
        ring[rr-1] = n_occupied || s_occupied;
        base += ringlen;

/* libhealpix doesn't provide a map2alm_spin_iter interface which takes a number
 * of iterations - it only has the *_iter2 interface variant which works until
 * some abs/rel convergence is achieved. Write our own map2alm_spin_iter variant
 * for consistency with the other iterative options.
 */
template<typename T> void map2alm_spin_iter(
        const Healpix_Map<T> &map1, const Healpix_Map<T> &map2,
        Alm<xcomplex<T> > &alm1, Alm<xcomplex<T> > &alm2,
        int spin, int num_iter, const arr<double> &weight)
{
    planck_assert(map1.Scheme()==RING, "map2alm_spin_iter: maps must be in RING scheme");
    planck_assert(map1.conformable(map2), "map2alm_spin_iter: maps are not conformable");
    planck_assert(alm1.conformable(alm1), "map2alm_spin_iter: a_lm are not conformable");
    planck_assert(map1.fullyDefined() && map2.fullyDefined(), "map contains undefined pixels");

    sharp_cxxjob<T> job;
    job.set_weighted_Healpix_geometry(map1.Nside(), &weight[0]);
    job.set_triangular_alm_info(alm1.Lmax(), alm1.Mmax());
    job.map2alm_spin(&map1[0], &map2[0], &alm1(0,0), &alm2(0,0), spin, false);

    if (num_iter == 0)
        return;

    Healpix_Map<T> map1b(map1.Nside(), map1.Scheme(), SET_NSIDE);
    Healpix_Map<T> map2b(map1.Nside(), map1.Scheme(), SET_NSIDE);
    for (int iter = 0; iter < num_iter; ++iter) {
        job.alm2map_spin(&alm1(0,0), &alm2(0,0), &map1b[0], &map2b[0], spin, false);
        #pragma omp parallel for
        for (size_t ii = 0; ii < map1.Npix(); ++ii) {
            map1b[ii] = map1[ii] - map1b[ii];
            map2b[ii] = map2[ii] - map2b[ii];
        }
        job.map2alm_spin(&map1[0], &map2[0], &alm1(0,0), &alm2(0,0), spin, true);
    }
}

/* Externally callable function implementations */

#define DISPATCH_FN(name) \
    void MexFunction::mex_##name(ArgumentList& outputs, ArgumentList& inputs)

DISPATCH_FN(nest2ring) {
    healpix base = nsideorder(inputs[1]);
    TypedArray<int64_t> ml_ipix  = inputs[2];

    auto ipix  = buffer<int64_t>(ml_ipix);
    auto dims  = ml_ipix.getDimensions();
    auto npix  = ml_ipix.getNumberOfElements();
    auto rpix = factory.createBuffer<int64_t>(npix);

    #pragma omp parallel for
    for (size_t ii = 0; ii < npix; ++ii) {
        rpix[ii] = base.nest2ring(ipix[ii]);
    }

    outputs[0] = factory.createArrayFromBuffer(dims, move(rpix));
}

DISPATCH_FN(ring2nest) {
    healpix base = nsideorder(inputs[1]);
    TypedArray<int64_t> ml_ipix  = inputs[2];

    auto ipix  = buffer<int64_t>(ml_ipix);
    auto dims  = ml_ipix.getDimensions();
    auto npix  = ml_ipix.getNumberOfElements();
    auto rpix = factory.createBuffer<int64_t>(npix);

    #pragma omp parallel for
    for (size_t ii = 0; ii < npix; ++ii) {
        rpix[ii] = base.ring2nest(ipix[ii]);
    }

    outputs[0] = factory.createArrayFromBuffer(dims, move(rpix));
}

DISPATCH_FN(pix2ring) {
    healpix base = nsideorder(inputs[1], inputs[2]);
    TypedArray<int64_t> ml_ipix  = inputs[3];

    auto ipix  = buffer<int64_t>(ml_ipix);
    auto dims  = ml_ipix.getDimensions();
    auto npix  = ml_ipix.getNumberOfElements();
    auto ring = factory.createBuffer<int64_t>(npix);

    #pragma omp parallel for
    for (size_t ii = 0; ii < npix; ++ii) {
        ring[ii] = base.pix2ring(ipix[ii]);
    }

    outputs[0] = factory.createArrayFromBuffer(dims, move(ring));
}

DISPATCH_FN(pix2vec) {
    healpix base = nsideorder(inputs[1], inputs[2]);
    TypedArray<int64_t> ml_ipix  = inputs[3];

    auto ipix = buffer<int64_t>(ml_ipix);
    auto npix = ml_ipix.getNumberOfElements();
    auto x = factory.createBuffer<double>(npix);
    auto y = factory.createBuffer<double>(npix);
    auto z = factory.createBuffer<double>(npix);

    #pragma omp parallel for
    for (size_t ii = 0; ii < npix; ++ii) {
        auto r = base.pix2vec(ipix[ii]);
        x[ii] = r.x;
        y[ii] = r.y;
        z[ii] = r.z;
    }

    outputs[0] = factory.createArrayFromBuffer({npix}, move(x));
    outputs[1] = factory.createArrayFromBuffer({npix}, move(y));
    outputs[2] = factory.createArrayFromBuffer({npix}, move(z));
}

DISPATCH_FN(pix2zphi) {
    healpix base = nsideorder(inputs[1], inputs[2]);
    TypedArray<int64_t> ml_ipix  = inputs[3];

    auto ipix  = buffer<int64_t>(ml_ipix);
    auto dims  = ml_ipix.getDimensions();
    auto npix  = ml_ipix.getNumberOfElements();
    auto z     = factory.createBuffer<double>(npix);
    auto phi   = factory.createBuffer<double>(npix);

    #pragma omp parallel for
    for (size_t ii = 0; ii < npix; ++ii) {
        base.pix2zphi(ipix[ii], z[ii], phi[ii]);
    }

    outputs[0] = factory.createArrayFromBuffer(dims, move(z));
    outputs[1] = factory.createArrayFromBuffer(dims, move(phi));
}

DISPATCH_FN(pix2ang) {
    healpix base = nsideorder(inputs[1], inputs[2]);
    TypedArray<int64_t> ml_ipix  = inputs[3];

    auto ipix  = buffer<int64_t>(ml_ipix);
    auto dims  = ml_ipix.getDimensions();
    auto npix  = ml_ipix.getNumberOfElements();
    auto theta = factory.createBuffer<double>(npix);
    auto phi   = factory.createBuffer<double>(npix);

    #pragma omp parallel for
    for (size_t ii = 0; ii < npix; ++ii) {
        auto pointing = base.pix2ang(ipix[ii]);
        theta[ii] = pointing.theta;
        phi[ii]   = pointing.phi;
    }

    outputs[0] = factory.createArrayFromBuffer(dims, move(theta));
    outputs[1] = factory.createArrayFromBuffer(dims, move(phi));
}

DISPATCH_FN(vec2pix) {
    healpix base = nsideorder(inputs[1], inputs[2]);
    TypedArray<double> ml_x = inputs[3];
    TypedArray<double> ml_y = inputs[4];
    TypedArray<double> ml_z = inputs[5];

    auto dims = ml_x.getDimensions();
    auto npix = ml_x.getNumberOfElements();
    if (dims != ml_y.getDimensions() || dims != ml_z.getDimensions()) {
        error("vec2pix: x, y, and z must be same sizes");
    }
    auto x = buffer<double>(ml_x);
    auto y = buffer<double>(ml_y);
    auto z = buffer<double>(ml_z);
    auto ipix  = factory.createBuffer<int64_t>(npix);

    #pragma omp parallel for
    for (size_t ii = 0; ii < npix; ++ii) {
        auto v = vec3(x[ii], y[ii], z[ii]);
        ipix[ii] = base.vec2pix(v);
    }

    outputs[0] = factory.createArrayFromBuffer({npix}, move(ipix));
}

DISPATCH_FN(zphi2pix) {
    healpix base = nsideorder(inputs[1], inputs[2]);
    TypedArray<double>  ml_z     = inputs[3];
    TypedArray<double>  ml_phi   = inputs[4];

    auto dims = ml_z.getDimensions();
    auto npix = ml_z.getNumberOfElements();
    if (dims != ml_phi.getDimensions()) {
        error("zphi2pix: theta and phi must be same sizes");
    }
    auto z    = buffer<double>(ml_z);
    auto phi  = buffer<double>(ml_phi);
    auto ipix = factory.createBuffer<int64_t>(npix);

    #pragma omp parallel for
    for (size_t ii = 0; ii < npix; ++ii) {
        ipix[ii] = base.zphi2pix(z[ii], phi[ii]);
    }

    outputs[0] = factory.createArrayFromBuffer({npix}, move(ipix));
}

DISPATCH_FN(ang2pix) {
    healpix base = nsideorder(inputs[1], inputs[2]);
    TypedArray<double>  ml_theta = inputs[3];
    TypedArray<double>  ml_phi   = inputs[4];

    auto dims = ml_theta.getDimensions();
    auto npix = ml_theta.getNumberOfElements();
    if (dims != ml_phi.getDimensions()) {
        error("ang2pix: theta and phi must be same sizes");
    }
    auto theta = buffer<double>(ml_theta);
    auto phi   = buffer<double>(ml_phi);
    auto ipix  = factory.createBuffer<int64_t>(npix);

    #pragma omp parallel for
    for (size_t ii = 0; ii < npix; ++ii) {
        auto point = pointing(theta[ii], phi[ii]);
        ipix[ii] = base.ang2pix(point);
    }

    outputs[0] = factory.createArrayFromBuffer({npix}, move(ipix));
}

DISPATCH_FN(pix2xyf) {
    healpix base = nsideorder(inputs[1], inputs[2]);
    TypedArray<int64_t> ml_ipix = inputs[3];

    auto ipix  = buffer<int64_t>(ml_ipix);
    auto dims  = ml_ipix.getDimensions();
    auto npix  = ml_ipix.getNumberOfElements();
    auto x = factory.createBuffer<int>(npix);
    auto y = factory.createBuffer<int>(npix);
    auto f = factory.createBuffer<int>(npix);

    #pragma omp parallel for
    for (size_t ii = 0; ii < npix; ++ii) {
        base.pix2xyf(ipix[ii], x[ii], y[ii], f[ii]);
    }

    outputs[0] = factory.createArrayFromBuffer({npix}, move(x));
    outputs[1] = factory.createArrayFromBuffer({npix}, move(y));
    outputs[2] = factory.createArrayFromBuffer({npix}, move(f));
}

DISPATCH_FN(xyf2pix) {
    healpix base = nsideorder(inputs[1], inputs[2]);
    TypedArray<int> ml_x = inputs[3];
    TypedArray<int> ml_y = inputs[4];
    TypedArray<int> ml_f = inputs[5];

    auto npix = ml_x.getNumberOfElements();
    if (npix != ml_y.getNumberOfElements()) {
        error("xyf2pix: x and y must have the same number of elements");
    }
    if (npix != ml_f.getNumberOfElements()) {
        error("xyf2pix: x and f must have the same number of elements");
    }
    auto x = buffer<int>(ml_x);
    auto y = buffer<int>(ml_y);
    auto f = buffer<int>(ml_f);
    auto ipix = factory.createBuffer<int64_t>(npix);

    #pragma omp parallel for
    for (size_t ii = 0; ii < npix; ++ii) {
        ipix[ii] = base.xyf2pix(x[ii], y[ii], f[ii]);
    }

    outputs[0] = factory.createArrayFromBuffer({npix}, move(ipix));
}

DISPATCH_FN(query_disc) {
    healpix base = nsideorder(inputs[1], inputs[2]);
    auto [buf_rvec, len_rvec] = bufferlen<double>(inputs[3]);
    auto radius = scalar<double>(inputs[4]);
    auto inclusive = scalar<bool>(inputs[5]);

    if (len_rvec != 3) {
        error("query_disc: rvec must be length 3 vector");
    }
    if (radius < 0) {
        error("query_disc: radius must be positive value");
    }

    vec3 rvec = vec3(buf_rvec[0], buf_rvec[1], buf_rvec[2]);
    pointing point = pointing(rvec);
    rangeset<int64_t> pixset;

    point.normalize();
    if (inclusive)
        base.query_disc_inclusive(point, radius, pixset);
    else
        base.query_disc(point, radius, pixset);

    auto npix = pixset.nval();
    auto buf_ipix = factory.createBuffer<uint64_t>(npix);
    // based on _query_disc.pyx pixset_to_array()
    for (size_t ii = 0, jj = 0; ii < pixset.size(); ++ii) {
        for (auto ip = pixset.ivbegin(ii); ip < pixset.ivend(ii); ++ip) {
            buf_ipix[jj++] = ip;
        }
    }

    outputs[0] = factory.createArrayFromBuffer({(size_t)npix}, move(buf_ipix));
}

DISPATCH_FN(neighbors) {
    healpix base = nsideorder(inputs[1], inputs[2]);
    auto [ipix, npix] = bufferlen<int64_t>(inputs[3]);

    auto nei = factory.createBuffer<int64_t>(8 * npix);
    #pragma omp parallel for
    for (size_t ii = 0; ii < npix; ++ii) {
        fix_arr<int64_t,8> res;
        base.neighbors(ipix[ii], res);
        nei[0 + 8*ii] = res[0];
        nei[1 + 8*ii] = res[1];
        nei[2 + 8*ii] = res[2];
        nei[3 + 8*ii] = res[3];
        nei[4 + 8*ii] = res[4];
        nei[5 + 8*ii] = res[5];
        nei[6 + 8*ii] = res[6];
        nei[7 + 8*ii] = res[7];
    }

    outputs[0] = factory.createArrayFromBuffer({8, npix}, move(nei));
}

DISPATCH_FN(map2alm_iter) {
    healpix base = nsideorder(inputs[1]);
    auto [buf_mapT, len_mapT] = bufferlen<double>(inputs[2]);
    auto [buf_mapQ, len_mapQ] = bufferlen<double>(inputs[3]);
    auto [buf_mapU, len_mapU] = bufferlen<double>(inputs[4]);
    auto lmax = scalar<int32_t>(inputs[5]);
    auto mmax = scalar<int32_t>(inputs[6]);
    auto [buf_wght, len_wght] = bufferlen<double>(inputs[7]);
    auto iter = scalar<int32_t>(inputs[8]);

    if (len_mapQ != len_mapU) {
        error("map2alm_iter: mapQ and mapU must have the same length");
    }
    if (len_mapQ > 0 && len_mapT != len_mapQ) {
        error("map2alm_iter: mapT, mapQ, and mapU have mismatched lengths");
    }
    auto nalms = Alm_Base::Num_Alms(lmax, mmax);
    bool do_pol = len_mapQ > 0;

    auto mapT = healmap();
    auto mapQ = healmap();
    auto mapU = healmap();
    auto almsT = healalm();
    auto almsG = healalm();
    auto almsC = healalm();

    auto buf_almsT = factory.createBuffer<complex64>(nalms);
    auto buf_almsG = factory.createBuffer<complex64>(do_pol ? nalms : (size_t)0);
    auto buf_almsC = factory.createBuffer<complex64>(do_pol ? nalms : (size_t)0);
    {
        arr<complex64> almtmp(buf_almsT.get(), nalms);
        almsT.Set(almtmp, lmax, mmax);
        arr<double> maptmp(buf_mapT.get(), len_mapT);
        mapT.Set(maptmp, base.Scheme());
    }
    if (do_pol) {
        {
            arr<complex64> almtmp(buf_almsG.get(), nalms);
            almsG.Set(almtmp, lmax, mmax);
            arr<double> maptmp(buf_mapQ.get(), len_mapQ);
            mapQ.Set(maptmp, base.Scheme());
        }
        {
            arr<complex64> almtmp(buf_almsC.get(), nalms);
            almsC.Set(almtmp, lmax, mmax);
            arr<double> maptmp(buf_mapU.get(), len_mapU);
            mapU.Set(maptmp, base.Scheme());
        }
    }

    arr<double> rwghts(buf_wght.get(), len_wght);

    if (do_pol) {
        map2alm_pol_iter(mapT, mapQ, mapU,
                almsT, almsG, almsC, iter, rwghts);
    } else {
        map2alm_iter(mapT, almsT, iter, rwghts);
    }

    outputs[0] = factory.createArrayFromBuffer({nalms}, move(buf_almsT));
    outputs[1] = factory.createArrayFromBuffer({do_pol ? nalms : (size_t)0}, move(buf_almsG));
    outputs[2] = factory.createArrayFromBuffer({do_pol ? nalms : (size_t)0}, move(buf_almsC));
}


DISPATCH_FN(map2alm_pure) {
    healpix base = nsideorder(inputs[1], inputs[2]);
    auto [buf_mapQ, len_mapQ] = bufferlen<double>(inputs[3]);
    auto [buf_mapU, len_mapU] = bufferlen<double>(inputs[4]);
    auto [buf_mapPw, len_mapPw] = bufferlen<double>(inputs[5]);
    auto lmax = scalar<int32_t>(inputs[6]);
    auto mmax = scalar<int32_t>(inputs[7]);
    auto [buf_wght, len_wght] = bufferlen<double>(inputs[8]);
    auto iter = scalar<int32_t>(inputs[9]);
    auto pureE = scalar<bool>(inputs[10]);

    auto mapQ = healmap();
    auto mapU = healmap();
    auto papQ = healmap();
    auto papU = healmap();
    auto wapQ = healmap();
    auto wapU = healmap();
    auto mask = healmap();
    {
        arr<double> tmp(buf_mapQ.get(), len_mapQ);
        mapQ.Set(tmp, base.Scheme());
		papQ.SetNside(base.Nside(),base.Scheme());
    }
    {
        arr<double> tmp(buf_mapU.get(), len_mapU);
        mapU.Set(tmp, base.Scheme());
		papU.SetNside(base.Nside(),base.Scheme());
    }
    {
        arr<double> tmp(buf_mapPw.get(), len_mapPw);
		wapQ.Set(tmp, base.Scheme());
		wapU.SetNside(base.Nside(),base.Scheme());
        mask.SetNside(base.Nside(),base.Scheme());
    }

    arr<double> rwghts(buf_wght.get(), len_wght);

    auto nalms = Alm_Base::Num_Alms(lmax, mmax);
    auto almsG = healalm();
    auto almsC = healalm();
    auto plmsG = healalm();
    auto plmsC = healalm();
    auto wlmsG = healalm();
    auto wlmsC = healalm();

    auto buf_almsG = factory.createBuffer<complex64>(nalms);
    auto buf_almsC = factory.createBuffer<complex64>(nalms);
	{
        arr<complex64> tmp(buf_almsG.get(), nalms);
        almsG.Set(tmp, lmax, mmax);
    }
    {
        arr<complex64> tmp(buf_almsC.get(), nalms);
        almsC.Set(tmp, lmax, mmax);
    }
	{
		plmsG.Set(lmax, mmax);
		plmsC.Set(lmax, mmax);
	}
	{
		wlmsG.Set(lmax, mmax);
		wlmsC.Set(lmax, mmax);
	}
	
	arr<double> f_l;
	f_l.alloc(lmax+1);
	
	sharp_cxxjob<double> job;
	job.set_weighted_Healpix_geometry (base.Nside(), &rwghts[0]);
	job.set_triangular_alm_info (lmax, mmax);
	
	#pragma omp parallel for
	for (int ip=0; ip< wapQ.Npix(); ip++) {
		if (wapQ[ip]!=0.) mask[ip]=1.; 
		else mask[ip]=0.;
		wapU[ip]=0.;
	}	
	map2alm_spin_iter(job,wapQ,wapU,wlmsG,wlmsC,0,iter);

	#pragma omp parallel for
	for (int ip=0; ip< papQ.Npix(); ip++) {
		papQ[ip]=mapQ[ip]*wapQ[ip];
		papU[ip]=mapU[ip]*wapQ[ip];
	}
	
	map2alm_spin_iter(job,papQ,papU,almsG,almsC,2,iter);
	
	// Compute spin-1 mask
	for(int l=0;l<=lmax;l++) //The minus sign is because of the definition of E-modes
		f_l[l]=-sqrt(((double)l+1.)*(double)l);
		
	wlmsG.ScaleL(f_l);
	wlmsC.ScaleL(f_l);
	job.alm2map_spin(&wlmsG(0,0),&wlmsC(0,0),&wapQ[0],&wapU[0],1,false);
	
	// Product with spin-1 mask
	#pragma omp parallel for
	for(int ip=0;ip<papQ.Npix();ip++) {
		papQ[ip]=(wapQ[ip]*mapQ[ip]+wapU[ip]*mapU[ip])*mask[ip];
		papU[ip]=(wapQ[ip]*mapU[ip]-wapU[ip]*mapQ[ip])*mask[ip];
	}
	
	// Compute SHT, multiply by 2*sqrt((l+1)!(l-2)!/((l-1)!(l+2)!)) and add to alm_out
	map2alm_spin_iter(job,papQ,papU,plmsG,plmsC,1,iter);
	
	for(int l=0;l<=lmax;l++) {
		if(l>1)
		  f_l[l]=2./sqrt(((double)l+2.)*((double)l-1.));
		else
		  f_l[l]=0;
	}	
	
	#pragma omp parallel for
	for (int m=0; m<=mmax; m++) {
		for (int l=m; l<=lmax; l++) {
			if (pureE==true) almsG(l,m)+=plmsG(l,m)*f_l[l];
			almsC(l,m)+=plmsC(l,m)*f_l[l];
		}
	}
	
	// Compute spin-2 mask
	for(int l=0;l<=lmax;l++) { //The extra minus sign is because of the scalar SHT below (E-mode def for s=0)
		if(l>1)
			f_l[l]=-sqrt(((double)l-1.)*((double)l+2.));
		else
			f_l[l]=0;
	}
	
	wlmsG.ScaleL(f_l);
	wlmsC.ScaleL(f_l);
	job.alm2map_spin(&wlmsG(0,0),&wlmsC(0,0),&wapQ[0],&wapU[0],2,false);

	// Product with spin-2 mask
	#pragma omp parallel for
	for(int ip=0;ip<papQ.Npix();ip++) {
		papQ[ip]=(wapQ[ip]*mapQ[ip]+wapU[ip]*mapU[ip])*mask[ip];
		papU[ip]=(wapQ[ip]*mapU[ip]-wapU[ip]*mapQ[ip])*mask[ip];
	}
	
	// Compute SHT, multiply by sqrt((l-2)!/(l+2)!) and add to alm_out
	map2alm_iter(papQ, plmsG, iter, rwghts);
	map2alm_iter(papU, plmsC, iter, rwghts);
	
	for(int l=0;l<=lmax;l++) {
		if(l>1)
			f_l[l]=1./sqrt(((double)l+2.)*((double)l+1.)*(double)l*((double)l-1.));
		else
			f_l[l]=0;
	}
	#pragma omp parallel for
	for (int m=0; m<=mmax; m++) {
		for (int l=m; l<=lmax; l++) {
			if (pureE==true) almsG(l,m)+=plmsG(l,m)*f_l[l];
			almsC(l,m)+=plmsC(l,m)*f_l[l];
		}
	}
	
    outputs[0] = factory.createArrayFromBuffer({nalms}, move(buf_almsG));
    outputs[1] = factory.createArrayFromBuffer({nalms}, move(buf_almsC));
}


void apodize(const Healpix_Map<double> & distmap, Healpix_Map<double> & mask, double radius, bool pixbool,bool inside){
	double omega = M_PI/radius;
	int    sign_coef ;
	
	if (pixbool) sign_coef = 1;  // Transition from 0 to 1 if pixbool is 1 
	else         sign_coef = -1; // Transition from 1 to 0 if pixbool is 0
	if (inside) sign_coef = - sign_coef; // invert the transition to apodize inside
	double oradius =  radius;//60/60*degr2rad; 
	#pragma omp parallel for  
	for(int i = 0; i<distmap.Npix(); i++)
	{
		if ((distmap[i] >= oradius) || (approx(distmap[i],Healpix_undef))) mask[i] = 1-pixbool;
		//if (approx(distmap[i],Healpix_undef)) mask[i] = 1-pixbool;
		else if ((inside) && (distmap[i]>=radius)) mask[i] = pixbool;
		else if ((!inside) && (distmap[i]>=radius)) mask[i] = 1-pixbool;
		else if (distmap[i] ==0) mask[i] = pixbool;
		else{
			mask[i] = 0.5 + sign_coef*(0.5* cos(omega * distmap[i]));
		}
	}
}

// Distance map : starting value is 2*pi outside of the mask and 0 inside
// To apodize inside : for each point inside the mask, get all pixels inside a disk of radius "radius"
// If the disk cross the border of the mask, update the distance value of all pixel inside of the mask 
// which are close to the border with the minimum distance value to the border 
// Credits: Marc Betoule, Xpure
void computeDistance_in(const Healpix_Map<double> & mask, Healpix_Map<double> & distmap, double radius,bool pixbool){
  vector<int> listpix;
  double mind;
  for(int i = 0; i<mask.Npix();i++){
    // for each point in the mask
    if(mask[i] == pixbool){
      pointing pixcenter = mask.pix2ang(i);
      // get a disk
      mask.query_disc(pixcenter ,radius, listpix);
      // maximum value of the distance map is radius
      mind = 2*M_PI;//radius;
      for(vector<int>::iterator p = listpix.begin(); p!=listpix.end(); p++){
		double v = distmap[*p];
		// for each point of the disk outside of the mask
		if(v==2*M_PI){
			pointing p2 = mask.pix2ang(*p);
		    // get the distance
			double da = angdist(pixcenter,p2);
		    // get the minimum distance
		    if (da<mind) mind = da;
		}
      }
      // update the distance value if a minimum is found
      //      distmap[i] = mind==radius?0:mind;
      distmap[i] = mind>radius?0:mind;
    }
  }
}


// Distance map : starting value is 2*pi outside of the mask and 0 inside
// To apodize outside : for each point inside the mask, get all pixels inside a disk of radius "radius"
// If the disk cross the border of the mask, update the distance value of all pixel outside of the mask 
// which are close to the border with the minimum distance value to the border 
void computeDistance_out(const Healpix_Map<double> & mask, Healpix_Map<double> & distmap, double radius,bool pixbool){
  vector<int> listpix;
  for(int i = 0; i<mask.Npix();i++){
    // for each point in the mask
    if(mask[i] == pixbool){
      pointing pixcenter = mask.pix2ang(i);
      // get a disk
      mask.query_disc(pixcenter ,radius, listpix);
      for(vector<int>::iterator p = listpix.begin(); p!=listpix.end(); p++){
	double v = distmap[*p];
	// for each point of the disk outside of the mask
      	if(v != 0){
      	  pointing p2 = mask.pix2ang(*p);
	  // get the distance
      	  double da = angdist(pixcenter,p2);
	  // get the minimum distance
      	  distmap[*p] = v<da?v:da; 
	}
      }
    }
  }
}


DISPATCH_FN(shrink_mask) {
    healpix base = nsideorder(inputs[1], inputs[2]);
    auto [buf_map, len_map] = bufferlen<double>(inputs[3]);
    TypedArray<double> radius_ = inputs[4];
	double radius=radius_[0];
	
	double radiusrad=radius*M_PI/180;

    auto npix = 12 * base.Nside() * base.Nside();

    auto map = healmap();
	auto amap = healmap();
    auto buf_amap = factory.createBuffer<double>(npix);
    auto dum_map = healmap();
    {
        arr<double> tmp(buf_map.get(), len_map);
        map.Set(tmp, base.Scheme());
		dum_map.SetNside(base.Nside(),base.Scheme());
    }
    {
        arr<double> tmp(buf_amap.get(), npix);
        amap.Set(tmp, base.Scheme());
    }
	#pragma omp for schedule(dynamic)
	for(int ip=0;ip<npix;ip++) {
		dum_map[ip]=1;
	}
	
	#pragma omp parallel default(none) shared(npix,map,dum_map,radiusrad,base)
	{
    rangeset<int64_t> listpix;

	#pragma omp for schedule(dynamic)
	for(int ip=0;ip<npix;ip++) {
		if(map[ip]<=0) {
			base.query_disc(pointing(base.pix2vec(ip)),radiusrad,listpix);
			int n = listpix.size();
			for(int i=0;i<n;i++) {
				for(int j=listpix.ivbegin(i);j<listpix.ivend(i);j++) {
					dum_map[j]=0;
				}
			}	
		}
	} //end omp for
	} //end omp parallel
	
	#pragma omp parallel for
	for(size_t ip=0;ip<amap.Npix();ip++) {
		amap[ip]=map[ip]*dum_map[ip];
	}

    outputs[0] = factory.createArrayFromBuffer({npix}, move(buf_amap));
}

DISPATCH_FN(smooth_mask) {
    healpix base = nsideorder(inputs[1], inputs[2]);
    auto [buf_map, len_map] = bufferlen<double>(inputs[3]);
    TypedArray<double> radius_ = inputs[4];
	double radius=radius_[0];
    auto [buf_wght, len_wght] = bufferlen<double>(inputs[5]);
	
	double radiusrad=radius*M_PI/180;

    auto npix = 12 * base.Nside() * base.Nside();

    auto map = healmap();
	auto amap = healmap();
    auto buf_amap = factory.createBuffer<double>(npix);
    {
        arr<double> tmp(buf_map.get(), len_map);
        map.Set(tmp, base.Scheme());
    }
    {
        arr<double> tmp(buf_amap.get(), npix);
        amap.Set(tmp, base.Scheme());
    }
	int32_t lmax = 3 * base.Nside() - 1;
	int32_t mmax = lmax ;
	
	auto dum_alms = healalm();
	{
		dum_alms.Set(lmax, mmax);
	}

    arr<double> rwghts(buf_wght.get(), len_wght);

	map2alm_iter(map, dum_alms, 3, rwghts);

	arr<double> f_l;
	f_l.alloc(lmax+1);
	
	double sigma=0.00012352884853326381*radiusrad*180*60*2.355/M_PI;
	for(size_t l=0;l<=lmax;l++) 
		f_l[l]=exp(-0.5*l*(l+1)*sigma*sigma);
		
	#pragma omp parallel for
	for (size_t m=0; m<=mmax; ++m) {
		for (size_t l=m; l<=lmax; ++l) {
			dum_alms(l,m)*=f_l[l];
		}
	}
	
	alm2map(dum_alms, amap);

    outputs[0] = factory.createArrayFromBuffer({npix}, move(buf_amap));
}

DISPATCH_FN(apodize_mask) {
    healpix base = nsideorder(inputs[1], inputs[2]);
    auto [buf_map, len_map] = bufferlen<double>(inputs[3]);
    TypedArray<double> radius_ = inputs[4];
	double radius=radius_[0]*M_PI/180;

    auto npix = 12 * base.Nside() * base.Nside();

    auto map = healmap();
	auto amap = healmap();
    auto buf_amap = factory.createBuffer<double>(npix);
    auto distmap = healmap();
    {
        arr<double> tmp(buf_map.get(), len_map);
        map.Set(tmp, base.Scheme());
		distmap.SetNside(base.Nside(),base.Scheme());
    }
    {
        arr<double> tmp(buf_amap.get(), npix);
        amap.Set(tmp, base.Scheme());
    }
	distmap.fill(2*M_PI);
	
	bool pixbool  = 0;
	
	#pragma omp parallel for
	for(int i = 0; i<map.Npix(); i++){
		if(map[i] == pixbool) distmap[i] = 0;}
	
	computeDistance_out(map, distmap, radius, pixbool);
	
	// check that the radius of the input distance map is compatible
    double min, max;
    distmap.minmax(min,max);
    double pixsize = sqrt(M_PI/(3*distmap.Nside()*distmap.Nside()));
    if(max < radius-pixsize){
      cout << "radius too large for the input dist map\n";
      exit(-1);
    }
	
    // call the right function to apodize
    apodize(distmap, amap, radius, pixbool, false);
	
	#pragma omp for schedule(dynamic)
	for(int ip=0;ip<npix;ip++) {
		amap[ip]=amap[ip];
	}

    outputs[0] = factory.createArrayFromBuffer({npix}, move(buf_amap));
}

DISPATCH_FN(map2alm_spin_iter) {
    healpix base = nsideorder(inputs[1]);
    auto [buf_map1, len_map1] = bufferlen<double>(inputs[2]);
    auto [buf_map2, len_map2] = bufferlen<double>(inputs[3]);
    auto spin = scalar<int32_t>(inputs[4]);
    auto lmax = scalar<int32_t>(inputs[5]);
    auto mmax = scalar<int32_t>(inputs[6]);
    auto [buf_wght, len_wght] = bufferlen<double>(inputs[7]);
    auto iter = scalar<int32_t>(inputs[8]);

    if (len_map1 != len_map2) {
        error("map2alm_spin_iter: map1 and map2 must have the same length");
    }
    auto nalms = Alm_Base::Num_Alms(lmax, mmax);

    auto map1 = healmap();
    auto map2 = healmap();
    auto alms1 = healalm();
    auto alms2 = healalm();

    auto buf_alms1 = factory.createBuffer<complex64>(nalms);
    auto buf_alms2 = factory.createBuffer<complex64>(nalms);
    {
        arr<complex64> almtmp(buf_alms1.get(), nalms);
        alms1.Set(almtmp, lmax, mmax);
        arr<double> maptmp(buf_map1.get(), len_map1);
        map1.Set(maptmp, base.Scheme());
    }
    {
        arr<complex64> almtmp(buf_alms2.get(), nalms);
        alms2.Set(almtmp, lmax, mmax);
        arr<double> maptmp(buf_map2.get(), len_map2);
        map2.Set(maptmp, base.Scheme());
    }

    arr<double> rwghts(buf_wght.get(), len_wght);
    map2alm_spin_iter(map1, map2, alms1, alms2, spin, iter, rwghts);

    outputs[0] = factory.createArrayFromBuffer({nalms}, move(buf_alms1));
    outputs[1] = factory.createArrayFromBuffer({nalms}, move(buf_alms2));
}

DISPATCH_FN(alm2map_polonly) {
    auto lmax = scalar<int32_t>(inputs[1]);
    auto mmax = scalar<int32_t>(inputs[2]);
    auto [buf_almsG, len_almsG] = bufferlen<complex64>(inputs[3]);
    auto [buf_almsC, len_almsC] = bufferlen<complex64>(inputs[4]);
    healpix base = nsideorder(inputs[5], inputs[6]);
    auto [buf_wght, len_wght] = bufferlen<double>(inputs[7]);

    arr<double> rwghts(buf_wght.get(), len_wght);

    auto almsG = healalm();
    auto almsC = healalm();
    {
        arr<complex64> tmp(buf_almsG.get(), len_almsG);
        almsG.Set(tmp, lmax, mmax);
    }
    {
        arr<complex64> tmp(buf_almsC.get(), len_almsC);
        almsC.Set(tmp, lmax, mmax);
    }

    auto npix = 12 * base.Nside() * base.Nside();
    auto mapQ = healmap();
    auto mapU = healmap();
    auto buf_mapQ = factory.createBuffer<double>(npix);
    auto buf_mapU = factory.createBuffer<double>(npix);
    {
        arr<double> tmp(buf_mapQ.get(), npix);
        mapQ.Set(tmp, base.Scheme());
    }
    {
        arr<double> tmp(buf_mapU.get(), npix);
        mapU.Set(tmp, base.Scheme());
    }

    sharp_cxxjob<double> job;

	job.set_weighted_Healpix_geometry (base.Nside(), &rwghts[0]);
	job.set_triangular_alm_info (lmax, mmax);
	
	job.alm2map_spin(&almsG(0,0),&almsC(0,0),&mapQ[0],&mapU[0],2,false);

    outputs[0] = factory.createArrayFromBuffer({npix}, move(buf_mapQ));
    outputs[1] = factory.createArrayFromBuffer({npix}, move(buf_mapU));
}

DISPATCH_FN(alm2map) {
    auto lmax = scalar<int32_t>(inputs[1]);
    auto mmax = scalar<int32_t>(inputs[2]);
    auto [buf_almsT, len_almsT] = bufferlen<complex64>(inputs[3]);
    auto [buf_almsG, len_almsG] = bufferlen<complex64>(inputs[4]);
    auto [buf_almsC, len_almsC] = bufferlen<complex64>(inputs[5]);
    healpix base = nsideorder(inputs[6]);

    if (len_almsG != len_almsC) {
        error("alm2map: almsG and almsC must have the same size");
    }
    if (len_almsG > 0 && len_almsT != len_almsT) {
        error("alm2map: almsT, almsG, and almsC have mismatched sizes");
    }
    auto npix = (size_t)12 * base.Nside() * base.Nside();
    bool do_pol = len_almsG > 0;

    auto almsT = healalm();
    auto almsG = healalm();
    auto almsC = healalm();
    auto mapT = healmap();
    auto mapQ = healmap();
    auto mapU = healmap();

    auto buf_mapT = factory.createBuffer<double>(npix);
    auto buf_mapQ = factory.createBuffer<double>(do_pol ? npix : (size_t)0);
    auto buf_mapU = factory.createBuffer<double>(do_pol ? npix : (size_t)0);

    {
        arr<complex64> almtmp(buf_almsT.get(), len_almsT);
        almsT.Set(almtmp, lmax, mmax);
        arr<double> maptmp(buf_mapT.get(), npix);
        mapT.Set(maptmp, base.Scheme());
    }
    if (do_pol) {
        {
            arr<complex64> almtmp(buf_almsG.get(), len_almsG);
            almsG.Set(almtmp, lmax, mmax);
            arr<double> maptmp(buf_mapQ.get(), npix);
            mapQ.Set(maptmp, base.Scheme());
        }
        {
            arr<complex64> almtmp(buf_almsC.get(), len_almsC);
            almsC.Set(almtmp, lmax, mmax);
            arr<double> maptmp(buf_mapU.get(), npix);
            mapU.Set(maptmp, base.Scheme());
        }
    }

    if (do_pol) {
        alm2map_pol(almsT, almsG, almsC, mapT, mapQ, mapU);
    } else {
        alm2map(almsT, mapT);
    }

    outputs[0] = factory.createArrayFromBuffer({npix}, move(buf_mapT));
    outputs[1] = factory.createArrayFromBuffer({do_pol ? npix : (size_t)0}, move(buf_mapQ));
    outputs[2] = factory.createArrayFromBuffer({do_pol ? npix : (size_t)0}, move(buf_mapU));
}

DISPATCH_FN(alm2map_spin) {
    auto lmax = scalar<int32_t>(inputs[1]);
    auto mmax = scalar<int32_t>(inputs[2]);
    auto [buf_alms1, len_alms1] = bufferlen<complex64>(inputs[3]);
    auto [buf_alms2, len_alms2] = bufferlen<complex64>(inputs[4]);
    auto spin = scalar<int32_t>(inputs[5]);
    healpix base = nsideorder(inputs[6]);

    if (len_alms1 != len_alms2) {
        error("alm2map_spin: almsG and almsC must have the same size");
    }
    auto npix = (size_t)12 * base.Nside() * base.Nside();

    auto alms1 = healalm();
    auto alms2 = healalm();
    auto map1 = healmap();
    auto map2 = healmap();

    auto buf_map1 = factory.createBuffer<double>(npix);
    auto buf_map2 = factory.createBuffer<double>(npix);
    {
        arr<complex64> almtmp(buf_alms1.get(), len_alms1);
        alms1.Set(almtmp, lmax, mmax);
        arr<double> maptmp(buf_map1.get(), npix);
        map1.Set(maptmp, base.Scheme());
    }
    {
        arr<complex64> almtmp(buf_alms2.get(), len_alms2);
        alms2.Set(almtmp, lmax, mmax);
        arr<double> maptmp(buf_map2.get(), npix);
        map2.Set(maptmp, base.Scheme());
    }

    alm2map_spin(alms1, alms2, map1, map2, spin, false);

    outputs[0] = factory.createArrayFromBuffer({npix}, move(buf_map1));
    outputs[1] = factory.createArrayFromBuffer({npix}, move(buf_map2));
}

DISPATCH_FN(alm2map_der1) {
    auto lmax = scalar<int32_t>(inputs[1]);
    auto mmax = scalar<int32_t>(inputs[2]);
    auto [buf_alms, len_alms] = bufferlen<complex64>(inputs[3]);
    healpix base = nsideorder(inputs[4]);

    auto npix = (size_t)12 * base.Nside() * base.Nside();
    auto alms = healalm();
    {
        arr<complex64> tmp(buf_alms.get(), len_alms);
        alms.Set(tmp, lmax, mmax);
    }

    auto map    = healmap();
    auto mapdth = healmap();
    auto mapdph = healmap();

    auto buf_map    = factory.createBuffer<double>(npix);
    auto buf_mapdth = factory.createBuffer<double>(npix);
    auto buf_mapdph = factory.createBuffer<double>(npix);
    {
        arr<double> tmp(buf_map.get(), npix);
        map.Set(tmp, base.Scheme());
    }
    {
        arr<double> tmp(buf_mapdth.get(), npix);
        mapdth.Set(tmp, base.Scheme());
    }
    {
        arr<double> tmp(buf_mapdph.get(), npix);
        mapdph.Set(tmp, base.Scheme());
    }

    alm2map_der1(alms, map, mapdth, mapdph);

    outputs[0] = factory.createArrayFromBuffer({npix}, move(buf_map));
    outputs[1] = factory.createArrayFromBuffer({npix}, move(buf_mapdth));
    outputs[2] = factory.createArrayFromBuffer({npix}, move(buf_mapdph));
}

DISPATCH_FN(alm2cl) {
    auto lmax = scalar<int32_t>(inputs[1]);
    auto mmax = scalar<int32_t>(inputs[2]);
    if (lmax < 0 || mmax < 0) {
        error("lmax and mmax must be non-negative values");
    }
    auto [buf_alms1, len_alms1] = bufferlen<complex64>(inputs[3]);
    auto [buf_alms2, len_alms2] = bufferlen<complex64>(inputs[4]);

    auto alms1 = healalm();
    auto alms2 = healalm();
    {
        arr<complex64> tmp(buf_alms1.get(), len_alms1);
        alms1.Set(tmp, lmax, mmax);
    }
    {
        arr<complex64> tmp(buf_alms2.get(), len_alms2);
        alms2.Set(tmp, lmax, mmax);
    }

    auto powspec = factory.createBuffer<double>((size_t)lmax + 1);

    // Essentially extract_crosspowspec from alm_powspec_tools.{h,cc}, but
    // avoids PowSpec object which can't be given a pre-allocated buffer to
    // make use of.
    for (int ll = 0; ll <= lmax; ++ll) {
        powspec[ll] = alms1(ll,0).real() * alms2(ll,0).real();
        int endm = ll < mmax ? ll : mmax;
        for (int mm = 0; mm <= endm; ++mm) {
            powspec[ll] += 2 * (alms1(ll,mm).real() * alms2(ll,mm).real()
                              + alms1(ll,mm).imag() * alms2(ll,mm).imag());
        }
        powspec[ll] /= (2*ll + 1);
    }

    outputs[0] = factory.createArrayFromBuffer({(size_t)lmax+1}, move(powspec));
}

DISPATCH_FN(almxfl) {
    auto lmax = scalar<int32_t>(inputs[1]);
    auto mmax = scalar<int32_t>(inputs[2]);
    auto [buf_alms, len_alms] = bufferlen<complex64>(inputs[3]);
    auto [buf_fl,   len_fl]   = bufferlen<complex64>(inputs[4]);

    if (len_alms != alm_getn(lmax, mmax)) {
        error("almxfl: Length of alms(=%d) is incompatible with given"
              " lmax(=%d) and mmax(=%d)",
                len_alms, lmax, mmax);
    }
    size_t ii = 0;
    for (size_t mm = 0; mm <= mmax; ++mm) {
        for (size_t ll = mm; ll <= lmax; ++ll) {
            complex64 fl = ll < len_fl ? buf_fl[ll] : 0.0;
            buf_alms[ii] *= fl;
            ++ii;
        }
    }

    outputs[0] = factory.createArrayFromBuffer({len_alms}, move(buf_alms));
}


DISPATCH_FN(smoothing_pol) {
    healpix base = nsideorder(inputs[1], inputs[2]);
    auto [buf_mapQ, len_mapQ] = bufferlen<double>(inputs[3]);
    auto [buf_mapU, len_mapU] = bufferlen<double>(inputs[4]);
    auto [buf_fle,   len_fle]   = bufferlen<double>(inputs[5]);
    auto [buf_flb,   len_flb]   = bufferlen<double>(inputs[6]);
    auto lmax = scalar<int32_t>(inputs[7]);
    auto mmax = scalar<int32_t>(inputs[8]);
    auto mmin = scalar<int32_t>(inputs[9]);
    auto [buf_wght, len_wght] = bufferlen<double>(inputs[10]);
    auto iter = scalar<int32_t>(inputs[11]);

    auto mapQ = healmap();
    auto mapU = healmap();
    {
        arr<double> tmp(buf_mapQ.get(), len_mapQ);
        mapQ.Set(tmp, base.Scheme());
    }
    {
        arr<double> tmp(buf_mapU.get(), len_mapU);
        mapU.Set(tmp, base.Scheme());
    }

    arr<double> rwghts(buf_wght.get(), len_wght);

    auto nalms = Alm_Base::Num_Alms(lmax, mmax);
	
    auto almsG = healalm();
    auto almsC = healalm();
	almsG.Set(lmax, mmax);
	almsC.Set(lmax, mmax);
		
	sharp_cxxjob<double> job;
	job.set_weighted_Healpix_geometry (base.Nside(), &rwghts[0]);
	job.set_triangular_alm_info (lmax, mmax);

	map2alm_spin_iter(job,mapQ,mapU,almsG,almsC,2,iter);
	
	double fle,flb;
	
	#pragma omp parallel for
	for (int m=0; m<=mmax; m++) {
		for (int l=m; l<=lmax; l++) {
			if ( l < len_fle && m >= mmin ) fle=buf_fle[l];
			else fle=0.0;
			if ( l < len_flb && m >= mmin ) flb=buf_flb[l];
			else flb=0.0;
			almsG(l,m)*=fle;
			almsC(l,m)*=flb;
		}
	}
	
	job.alm2map_spin(&almsG(0,0),&almsC(0,0),&mapQ[0],&mapU[0],2,false);

    outputs[0] = factory.createArrayFromBuffer({len_mapQ}, move(buf_mapQ));
    outputs[1] = factory.createArrayFromBuffer({len_mapU}, move(buf_mapU));
}

DISPATCH_FN(smoothing) {
    healpix base = nsideorder(inputs[1], inputs[2]);
    auto [buf_map, len_map] = bufferlen<double>(inputs[3]);
    auto [buf_fl,   len_fl]   = bufferlen<double>(inputs[4]);
    auto lmax = scalar<int32_t>(inputs[5]);
    auto mmax = scalar<int32_t>(inputs[6]);
    auto mmin = scalar<int32_t>(inputs[7]);
    auto [buf_wght, len_wght] = bufferlen<double>(inputs[8]);
    auto iter = scalar<int32_t>(inputs[9]);

    auto map = healmap();
    {
        arr<double> tmp(buf_map.get(), len_map);
        map.Set(tmp, base.Scheme());
    }

    arr<double> rwghts(buf_wght.get(), len_wght);

    auto nalms = Alm_Base::Num_Alms(lmax, mmax);
	
    auto alms = healalm();
	alms.Set(lmax, mmax);

	map2alm_iter(map,alms,iter,rwghts);
	
	double fl;
	
	#pragma omp parallel for
	for (int m=0; m<=mmax; m++) {
		for (int l=m; l<=lmax; l++) {
			if ( l < len_fl && m >= mmin ) fl=buf_fl[l];
			else fl=0.0;
			alms(l,m)*=fl;
		}
	}
	
	alm2map(alms, map);

    outputs[0] = factory.createArrayFromBuffer({len_map}, move(buf_map));
}

// Forward declaration of HEALPix's Trafo object. Header isn't installed, so
// declare directly. Also cheat and declare a static private member as public.
enum coordsys { Ecliptic, Equatorial, Galactic };
class Trafo {
public:
    static void coordsys2matrix (double iepoch, double oepoch,
                                 coordsys isys, coordsys osys,
                                 rotmatrix &matrix);
};
inline void coordsys2matrix (int num, rotmatrix& rm) {
    #define c2m Trafo::coordsys2matrix
    switch (num) {
        case  1: c2m(2000, 2000, Equatorial, Galactic,   rm); break;
        case  2: c2m(2000, 2000, Galactic,   Equatorial, rm); break;
        case  3: c2m(2000, 2000, Equatorial, Ecliptic,   rm); break;
        case  4: c2m(2000, 2000, Ecliptic,   Equatorial, rm); break;
        case  5: c2m(2000, 2000, Ecliptic,   Galactic,   rm); break;
        case  6: c2m(2000, 2000, Galactic,   Ecliptic,   rm); break;
        case  7: c2m(1950, 1950, Equatorial, Galactic,   rm); break;
        case  8: c2m(1950, 1950, Galactic,   Equatorial, rm); break;
        case  9: c2m(1950, 1950, Equatorial, Ecliptic,   rm); break;
        case 10: c2m(1950, 1950, Ecliptic,   Equatorial, rm); break;
        case 11: c2m(1950, 1950, Ecliptic,   Galactic,   rm); break;
        case 12: c2m(1950, 1950, Galactic,   Ecliptic,   rm); break;
        default: c2m(2000, 2000, Galactic,   Galactic,   rm); break;
    }
    #undef c2m
}

DISPATCH_FN(rotate_alm_coord) {
    auto itransform = scalar<int32_t>(inputs[1]);
    auto lmax = scalar<int32_t>(inputs[2]);
    auto mmax = scalar<int32_t>(inputs[3]);
    auto [buf_almsT, len_almsT] = bufferlen<complex64>(inputs[4]);
    auto [buf_almsG, len_almsG] = bufferlen<complex64>(inputs[5]);
    auto [buf_almsC, len_almsC] = bufferlen<complex64>(inputs[6]);

    if (len_almsG != len_almsC) {
        error("rotate_alm: almsG and almsC must have the same size");
    }
    if (len_almsG > 0 && len_almsT != len_almsT) {
        error("rotate_alm: almsT, almsG, and almsC have mismatched sizes");
    }
    auto almsT = healalm();
    auto almsG = healalm();
    auto almsC = healalm();
    {
        arr<complex64> tmp(buf_almsT.get(), len_almsT);
        almsT.Set(tmp, lmax, mmax);
    }
    if (len_almsG > 0) {
        {
            arr<complex64> tmp(buf_almsG.get(), len_almsG);
            almsG.Set(tmp, lmax, mmax);
        }
        {
            arr<complex64> tmp(buf_almsC.get(), len_almsC);
            almsC.Set(tmp, lmax, mmax);
        }
    }

    rotmatrix rm;
    coordsys2matrix(itransform, rm);
    if (len_almsG > 0) {
        rotate_alm(almsT, almsG, almsC, rm);
    } else {
        rotate_alm(almsT, rm);
    }

    outputs[0] = factory.createArrayFromBuffer({len_almsT}, move(buf_almsT));
    outputs[1] = factory.createArrayFromBuffer({len_almsG}, move(buf_almsG));
    outputs[2] = factory.createArrayFromBuffer({len_almsC}, move(buf_almsC));
}

DISPATCH_FN(rotate_alm_euler) {
    auto [euler, len_euler] = bufferlen<double>(inputs[1]);
    auto lmax = scalar<int32_t>(inputs[2]);
    auto mmax = scalar<int32_t>(inputs[3]);
    auto [buf_almsT, len_almsT] = bufferlen<complex64>(inputs[4]);
    auto [buf_almsG, len_almsG] = bufferlen<complex64>(inputs[5]);
    auto [buf_almsC, len_almsC] = bufferlen<complex64>(inputs[6]);

    if (len_euler != 3) {
        error("rotate_alm: euler must be 3-vector of angles");
    }
    if (len_almsG != len_almsC) {
        error("rotate_alm: almsG and almsC must have the same size");
    }
    if (len_almsG > 0 && len_almsT != len_almsT) {
        error("rotate_alm: almsT, almsG, and almsC have mismatched sizes");
    }
    auto almsT = healalm();
    auto almsG = healalm();
    auto almsC = healalm();
    {
        arr<complex64> tmp(buf_almsT.get(), len_almsT);
        almsT.Set(tmp, lmax, mmax);
    }
    if (len_almsG > 0) {
        {
            arr<complex64> tmp(buf_almsG.get(), len_almsG);
            almsG.Set(tmp, lmax, mmax);
        }
        {
            arr<complex64> tmp(buf_almsC.get(), len_almsC);
            almsC.Set(tmp, lmax, mmax);
        }
    }

    double psi = euler[0];
    double theta = euler[1];
    double phi = euler[2];
    if (len_almsG > 0) {
        rotate_alm(almsT, almsG, almsC, psi, theta, phi);
    } else {
        rotate_alm(almsT, psi, theta, phi);
    }

    outputs[0] = factory.createArrayFromBuffer({len_almsT}, move(buf_almsT));
    outputs[1] = factory.createArrayFromBuffer({len_almsG}, move(buf_almsG));
    outputs[2] = factory.createArrayFromBuffer({len_almsC}, move(buf_almsC));
}

DISPATCH_FN(rotate_alm_matrix) {
    auto [buf_rm, len_rm] = bufferlen<double>(inputs[1]);
    auto lmax = scalar<int32_t>(inputs[2]);
    auto mmax = scalar<int32_t>(inputs[3]);
    auto [buf_almsT, len_almsT] = bufferlen<complex64>(inputs[4]);
    auto [buf_almsG, len_almsG] = bufferlen<complex64>(inputs[5]);
    auto [buf_almsC, len_almsC] = bufferlen<complex64>(inputs[6]);

    if (len_rm != 9) {
        error("rotate_alm: matrix must be a 3-by-3 rotation matrix");
    }
    if (len_almsG != len_almsC) {
        error("rotate_alm: almsG and almsC must have the same size");
    }
    if (len_almsG > 0 && len_almsT != len_almsT) {
        error("rotate_alm: almsT, almsG, and almsC have mismatched sizes");
    }
    auto almsT = healalm();
    auto almsG = healalm();
    auto almsC = healalm();
    {
        arr<complex64> tmp(buf_almsT.get(), len_almsT);
        almsT.Set(tmp, lmax, mmax);
    }
    if (len_almsG > 0) {
        {
            arr<complex64> tmp(buf_almsG.get(), len_almsG);
            almsG.Set(tmp, lmax, mmax);
        }
        {
            arr<complex64> tmp(buf_almsC.get(), len_almsC);
            almsC.Set(tmp, lmax, mmax);
        }
    }

    // Matlab is column major, arguments in row-major order
    double a00 = buf_rm[0], a01 = buf_rm[3], a02 = buf_rm[6]; // row 1
    double a10 = buf_rm[1], a11 = buf_rm[4], a12 = buf_rm[7]; // row 2
    double a20 = buf_rm[2], a21 = buf_rm[5], a22 = buf_rm[8]; // row 3
    rotmatrix rm = rotmatrix(a00, a01, a02, a10, a11, a12, a20, a21, a22);
    if (len_almsG > 0) {
        rotate_alm(almsT, almsG, almsC, rm);
    } else {
        rotate_alm(almsT, rm);
    }

    outputs[0] = factory.createArrayFromBuffer({len_almsT}, move(buf_almsT));
    outputs[1] = factory.createArrayFromBuffer({len_almsG}, move(buf_almsG));
    outputs[2] = factory.createArrayFromBuffer({len_almsC}, move(buf_almsC));
}

DISPATCH_FN(scan_rings_observed) {
    auto [buf_map, len_map] = bufferlen<double>(inputs[1]);
    auto nside = static_cast<int>(sqrt(len_map / 12));

    auto buf_rings = factory.createBuffer<bool>(2*nside);

    _scan_rings_observed(nside, buf_map.get(),  buf_rings.get());
    outputs[0] = factory.createArrayFromBuffer({2*nside}, move(buf_rings));
}

DISPATCH_FN(compute_mcm) {
    auto lmax = scalar<int32_t>(inputs[1]);
    auto lmax_mask = scalar<int32_t>(inputs[2]);

    healpix base = nsideorder(inputs[3], inputs[4]);
    auto [buf_map, len_map] = bufferlen<double>(inputs[5]);
    auto [buf_wght, len_wght] = bufferlen<double>(inputs[6]);
    auto iter = scalar<int32_t>(inputs[7]);

    auto pe1 = scalar<int32_t>(inputs[8]);
    auto pe2 = scalar<int32_t>(inputs[9]);
    auto pb1 = scalar<int32_t>(inputs[10]);
    auto pb2 = scalar<int32_t>(inputs[11]);

    auto map = healmap();
    {
        arr<double> tmp(buf_map.get(), len_map);
        map.Set(tmp, base.Scheme());
    }

    arr<double> rwghts(buf_wght.get(), len_wght);

    auto alms = healalm();
	alms.Set(lmax_mask,lmax_mask);

    map2alm_iter(map, alms, iter, rwghts);
	
	std::vector<double> pcl_masks(lmax_mask+1);

    for (int ll = 0; ll <= lmax; ++ll) {
        pcl_masks[ll] = alms(ll,0).real() * alms(ll,0).real();
        int endm = ll < lmax ? ll : lmax;
        for (int mm = 0; mm <= endm; ++mm) {
            pcl_masks[ll] += 2 * (alms(ll,mm).real() * alms(ll,mm).real()
                              + alms(ll,mm).imag() * alms(ll,mm).imag());
        }
        pcl_masks[ll] /= (4*M_PI);
    }
	
	int s1=2;
	int s2=2;	
	
	int ncls=7;
	
    int pure_any;
    int npure_0s;
    int npure_ss;
    pure_any=pe1 || pe2 || pb1 || pb2;    
    
    if(pure_any) {
      npure_0s=2;
      npure_ss=3;
    }
    else {
      npure_0s=1;
      npure_ss=1;
    }
	
    std::vector<std::vector<std::vector<double>>> xi_00;
    std::vector<std::vector<std::vector<std::vector<double>>>> xi_0s;
    std::vector<std::vector<std::vector<std::vector<double>>>> xi_pp;
    std::vector<std::vector<std::vector<std::vector<double>>>> xi_mm;

	
    xi_00=std::vector<std::vector<std::vector<double>>>(1,std::vector<std::vector<double>>(lmax+1,std::vector<double>(lmax+1,0.0)));
    xi_0s=std::vector<std::vector<std::vector<std::vector<double>>>>(1,std::vector<std::vector<std::vector<double>>>(npure_0s,std::vector<std::vector<double>>(lmax+1,std::vector<double>(lmax+1,0.0))));
    xi_pp=std::vector<std::vector<std::vector<std::vector<double>>>>(1,std::vector<std::vector<std::vector<double>>>(npure_ss,std::vector<std::vector<double>>(lmax+1,std::vector<double>(lmax+1,0.0))));
    xi_mm=std::vector<std::vector<std::vector<std::vector<double>>>>(1,std::vector<std::vector<std::vector<double>>>(npure_ss,std::vector<std::vector<double>>(lmax+1,std::vector<double>(lmax+1,0.0))));

    compute_mcm_namaster(lmax, lmax_mask, 1, pcl_masks,s1,s2,pe1,pb1,pe2,pb2,1,-1,-1,-1,xi_00,xi_0s,xi_pp,xi_mm);

	std::vector<size_t> dims = {9,9,lmax+1,lmax+1};
	
	TypedArray<double> cm=factory.createArray<double>(dims);

    // #pragma omp parallel default(none) shared(lmax,s1,s2,pe1,pb1,pe2,pb2,c,cm)
    // {
      int ll2,ll3;
      int sign_overall=1;
      if((s1+s2) & 1)
        sign_overall=-1;
    
      // #pragma omp for schedule(dynamic)
      for(ll2=0;ll2<=lmax;ll2++) {
        for(ll3=0;ll3<=lmax;ll3++) {
          double fac=(2*ll3+1.)*sign_overall;
          cm[0][0][ll2][ll3]=fac*xi_00[0][ll2][ll3]; //TT,TT
          cm[1][1][ll2][ll3]=fac*xi_0s[0][pe2][ll2][ll3]; //TE,TE
          cm[4][4][ll2][ll3]=fac*xi_0s[0][pb2][ll2][ll3]; //TB,TB
          cm[2][3][ll2][ll3]=fac*xi_mm[0][pe2+pe2][ll2][ll3]; //EE,BB
          cm[5][8][ll2][ll3]=-fac*xi_mm[0][pe2+pb2][ll2][ll3]; //EB,BE
          cm[8][5][ll2][ll3]=-fac*xi_mm[0][pb2+pe2][ll2][ll3]; //BE,EB
          cm[3][2][ll2][ll3]=fac*xi_mm[0][pb2+pb2][ll2][ll3]; //BB,EE
          cm[2][2][ll2][ll3]=fac*xi_pp[0][pe2+pe2][ll2][ll3]; //EE,EE
          cm[5][5][ll2][ll3]=fac*xi_pp[0][pe2+pb2][ll2][ll3]; //EB,EB
          cm[8][8][ll2][ll3]=fac*xi_pp[0][pb2+pe2][ll2][ll3]; //BE,BE
          cm[3][3][ll2][ll3]=fac*xi_pp[0][pb2+pb2][ll2][ll3]; //BB,BB

        }
      } //end omp for
    // } //end omp parallel

    outputs[0] = cm;
}
