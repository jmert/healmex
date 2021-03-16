#include <mex.hpp>
#include <mexAdapter.hpp>

#include <complex>
#include <utility> // tuple
#include <vector>

#include <healpix_map.h>
#include <alm.h>
#include <alm_healpix_tools.h>
#include <alm_powspec_tools.h>
#include <powspec.h>
#include <rotmatrix.h>

using namespace std;
using namespace matlab::data;
using matlab::mex::ArgumentList;

using complex64 = complex<double>;
using healpix = Healpix_Base2;
using healmap = Healpix_Map<double>;
using healalm = Alm<complex64>;

/* Useful type predicates */

inline bool ischar(Array& a)      { return a.getType() == ArrayType::CHAR; }
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

template <typename T, typename A> buffer_ptr_t<T> buffer(A& array);

/*
 * MEX entry point
 */

enum libhealpix_mex_calls {
    id_heartbeat        = -1,
    id_nest2ring        =  1,
    id_ring2nest        =  2,
    id_pix2vec          = 11,
    id_pix2zphi         = 12,
    id_pix2ang          = 13,
    id_vec2pix          = 14,
    id_zphi2pix         = 15,
    id_ang2pix          = 16,
    id_pix2xyf          = 17,
    id_xyf2pix          = 18,
    id_map2alm_iter     = 53,
    id_alm2map          = 55,
    id_alm2map_der1     = 56,
    id_alm2cl           = 61,
    id_almxfl           = 62,
    id_rotate_alm_coord = 65,
    id_rotate_alm_euler = 66,
    id_rotate_alm_matrix= 67,
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

            case id_pix2vec:
                CHECK_NINOUT("pix2vec", 3, 3);
                CHECK_INPUT_SCALAR("pix2vec", "nside", 1);
                CHECK_INPUT_INT64("pix2vec", "nside", 1);
                CHECK_INPUT_CHAR("pix2vec", "order", 2);
                CHECK_INPUT_INT64("pix2vec", "ipix", 3);
                mex_pix2vec(outputs, inputs);
                break;

            case id_pix2zphi:
                CHECK_NINOUT("pix2zphi", 3, 2);
                CHECK_INPUT_SCALAR("pix2zphi", "nside", 1);
                CHECK_INPUT_INT64("pix2zphi", "nside", 1);
                CHECK_INPUT_CHAR("pix2zphi", "order", 2);
                CHECK_INPUT_INT64("pix2zphi", "ipix", 3);
                mex_pix2zphi(outputs, inputs);
                break;

            case id_pix2ang:
                CHECK_NINOUT("pix2ang", 3, 2);
                CHECK_INPUT_SCALAR("pix2ang", "nside", 1);
                CHECK_INPUT_INT64("pix2ang", "nside", 1);
                CHECK_INPUT_CHAR("pix2ang", "order", 2);
                CHECK_INPUT_INT64("pix2ang", "ipix", 3);
                mex_pix2ang(outputs, inputs);
                break;

            case id_vec2pix:
                CHECK_NINOUT("vec2pix", 5, 1);
                CHECK_INPUT_SCALAR("vec2pix", "nside", 1);
                CHECK_INPUT_INT64("vec2pix", "nside", 1);
                CHECK_INPUT_CHAR("vec2pix", "order", 2);
                CHECK_INPUT_DOUBLE("vec2pix", "x", 3);
                CHECK_INPUT_DOUBLE("vec2pix", "y", 3);
                CHECK_INPUT_DOUBLE("vec2pix", "z", 3);
                mex_vec2pix(outputs, inputs);
                break;

            case id_zphi2pix:
                CHECK_NINOUT("zphi2pix", 4, 1);
                CHECK_INPUT_SCALAR("zphi2pix", "nside", 1);
                CHECK_INPUT_INT64("zphi2pix", "nside", 1);
                CHECK_INPUT_CHAR("zphi2pix", "order", 2);
                CHECK_INPUT_DOUBLE("zphi2pix", "z", 3);
                CHECK_INPUT_DOUBLE("zphi2pix", "phi", 4);
                mex_zphi2pix(outputs, inputs);
                break;

            case id_ang2pix:
                CHECK_NINOUT("ang2pix", 4, 1);
                CHECK_INPUT_SCALAR("ang2pix", "nside", 1);
                CHECK_INPUT_INT64("ang2pix", "nside", 1);
                CHECK_INPUT_CHAR("ang2pix", "order", 2);
                CHECK_INPUT_DOUBLE("ang2pix", "theta", 3);
                CHECK_INPUT_DOUBLE("ang2pix", "phi", 4);
                mex_ang2pix(outputs, inputs);
                break;

            case id_pix2xyf:
                CHECK_NINOUT("pix2xyf", 3, 3);
                CHECK_INPUT_SCALAR("pix2xyf", "nside", 1);
                CHECK_INPUT_INT64("pix2xyf", "nside", 1);
                CHECK_INPUT_CHAR("pix2xyf", "order", 2);
                CHECK_INPUT_INT64("pix2xyf", "ipix", 3);
                mex_pix2xyf(outputs, inputs);
                break;

            case id_xyf2pix:
                CHECK_NINOUT("xyf2pix", 5, 1);
                CHECK_INPUT_SCALAR("xyf2pix", "nside", 1);
                CHECK_INPUT_INT64("xyf2pix", "nside", 1);
                CHECK_INPUT_CHAR("xyf2pix", "order", 2);
                CHECK_INPUT_INT("xyf2pix", "x", 3);
                CHECK_INPUT_INT("xyf2pix", "y", 4);
                CHECK_INPUT_INT("xyf2pix", "f", 5);
                mex_xyf2pix(outputs, inputs);
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

    DISPATCH_FN(pix2vec);
    DISPATCH_FN(pix2zphi);
    DISPATCH_FN(pix2ang);
    DISPATCH_FN(vec2pix);
    DISPATCH_FN(zphi2pix);
    DISPATCH_FN(ang2pix);
    DISPATCH_FN(pix2xyf);
    DISPATCH_FN(xyf2pix);

    DISPATCH_FN(map2alm_iter);
    DISPATCH_FN(alm2map);
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
                   const CharArray ml_order)
{
    int64_t nside = ml_nside[0];
    auto order = string2HealpixScheme(ml_order.toAscii());
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

DISPATCH_FN(alm2map) {
    auto lmax = scalar<int32_t>(inputs[1]);
    auto mmax = scalar<int32_t>(inputs[2]);
    auto [buf_almsT, len_almsT] = bufferlen<complex64>(inputs[3]);
    auto [buf_almsG, len_almsG] = bufferlen<complex64>(inputs[4]);
    auto [buf_almsC, len_almsC] = bufferlen<complex64>(inputs[5]);
    healpix base = nsideorder(inputs[6]);

    if (len_almsG != len_almsC) {
        error("rotate_alm: almsG and almsC must have the same size");
    }
    if (len_almsG > 0 && len_almsT != len_almsT) {
        error("rotate_alm: almsT, almsG, and almsC have mismatched sizes");
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
