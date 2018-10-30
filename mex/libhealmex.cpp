#include <mex.hpp>
#include <mexAdapter.hpp>

#include <complex>
#include <utility> // tuple
#include <vector>

#include <healpix_map.h>
#include <alm.h>
#include <alm_healpix_tools.h>
#include <powspec.h>

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

inline bool isscalar(Array& a) { return a.getNumberOfElements() == 1; }

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
    id_map2alm_iter     = 53,
    id_map2alm_pol_iter = 54,
    id_alm2map          = 55,
    id_alm2map_pol      = 56,
    id_alm2cl           = 61,
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
                CHECK_NINOUT("pix2vec", 3, 1);
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
                CHECK_NINOUT("vec2pix", 3, 1);
                CHECK_INPUT_SCALAR("vec2pix", "nside", 1);
                CHECK_INPUT_INT64("vec2pix", "nside", 1);
                CHECK_INPUT_CHAR("vec2pix", "order", 2);
                CHECK_INPUT_DOUBLE("vec2pix", "vec", 3);
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

            case id_map2alm_iter:
                CHECK_NINOUT("map2alm_iter", 7, 1);
                CHECK_INPUT_SCALAR("map2alm_iter", "nside", 1);
                CHECK_INPUT_INT64("map2alm_iter", "nside", 1);
                CHECK_INPUT_CHAR("map2alm_iter", "order", 2);
                CHECK_INPUT_DOUBLE("map2alm_iter", "map", 3);
                CHECK_INPUT_SCALAR("map2alm_iter", "lmax", 4);
                CHECK_INPUT_INT32("map2alm_iter", "lmax", 4);
                CHECK_INPUT_SCALAR("map2alm_iter", "mmax", 5);
                CHECK_INPUT_INT32("map2alm_iter", "mmax", 5);
                CHECK_INPUT_DOUBLE("map2alm_iter", "rwghts", 6);
                CHECK_INPUT_SCALAR("map2alm_iter", "iter", 7);
                CHECK_INPUT_INT32("map2alm_iter", "iter", 7);
                mex_map2alm_iter(outputs, inputs);
                break;

            case id_map2alm_pol_iter:
                CHECK_NINOUT("map2alm_pol_iter", 9, 3);
                CHECK_INPUT_SCALAR("map2alm_pol_iter", "nside", 1);
                CHECK_INPUT_INT64("map2alm_pol_iter", "nside", 1);
                CHECK_INPUT_CHAR("map2alm_pol_iter", "order", 2);
                CHECK_INPUT_DOUBLE("map2alm_pol_iter", "mapT", 3);
                CHECK_INPUT_DOUBLE("map2alm_pol_iter", "mapQ", 4);
                CHECK_INPUT_DOUBLE("map2alm_pol_iter", "mapU", 5);
                CHECK_INPUT_SCALAR("map2alm_pol_iter", "lmax", 6);
                CHECK_INPUT_INT32("map2alm_pol_iter", "lmax", 6);
                CHECK_INPUT_SCALAR("map2alm_pol_iter", "mmax", 7);
                CHECK_INPUT_INT32("map2alm_pol_iter", "mmax", 7);
                CHECK_INPUT_DOUBLE("map2alm_pol_iter", "rwghts", 8);
                CHECK_INPUT_SCALAR("map2alm_pol_iter", "iter", 9);
                CHECK_INPUT_INT32("map2alm_pol_iter", "iter", 9);
                mex_map2alm_pol_iter(outputs, inputs);
                break;

            case id_alm2map:
                CHECK_NINOUT("alm2map", 5, 1);
                CHECK_INPUT_SCALAR("alm2map", "lmax", 1);
                CHECK_INPUT_INT32("alm2map", "lmax", 1);
                CHECK_INPUT_SCALAR("alm2map", "mmax", 2);
                CHECK_INPUT_INT32("alm2map", "mmax", 2);
                CHECK_INPUT_COMPLEX64("alm2map", "alms", 3);
                CHECK_INPUT_SCALAR("alm2map", "nside", 4);
                CHECK_INPUT_INT64("alm2map", "nside", 4);
                CHECK_INPUT_CHAR("alm2map", "order", 5);
                mex_alm2map(outputs, inputs);
                break;

            case id_alm2map_pol:
                CHECK_NINOUT("alm2map_pol", 7, 3);
                CHECK_INPUT_SCALAR("alm2map_pol", "lmax", 1);
                CHECK_INPUT_INT32("alm2map_pol", "lmax", 1);
                CHECK_INPUT_SCALAR("alm2map_pol", "mmax", 2);
                CHECK_INPUT_INT32("alm2map_pol", "mmax", 2);
                CHECK_INPUT_COMPLEX64("alm2map_pol", "almsT", 3);
                CHECK_INPUT_COMPLEX64("alm2map_pol", "almsG", 4);
                CHECK_INPUT_COMPLEX64("alm2map_pol", "almsC", 5);
                CHECK_INPUT_SCALAR("alm2map_pol", "nside", 6);
                CHECK_INPUT_INT64("alm2map_pol", "nside", 6);
                CHECK_INPUT_CHAR("alm2map_pol", "order", 7);
                mex_alm2map_pol(outputs, inputs);
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

    DISPATCH_FN(map2alm_iter);
    DISPATCH_FN(map2alm_pol_iter);
    DISPATCH_FN(alm2map);
    DISPATCH_FN(alm2map_pol);

    DISPATCH_FN(alm2cl);

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
    auto vec  = factory.createBuffer<double>(3 * npix);
    auto x = &vec[0];
    auto y = &vec[npix];
    auto z = &vec[2*npix];

    #pragma omp parallel for
    for (size_t ii = 0; ii < npix; ++ii) {
        auto r = base.pix2vec(ipix[ii]);
        x[ii] = r.x;
        y[ii] = r.y;
        z[ii] = r.z;
    }

    outputs[0] = factory.createArrayFromBuffer({npix, 3}, move(vec));
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
    TypedArray<double>  ml_vec   = inputs[3];

    auto dims  = ml_vec.getDimensions();
    if (dims.size() != 2 || dims[1] != 3) {
        error("vec2pix: vec must be a N-by-3 matrix of Cartesian coordinates");
    }
    auto npix  = dims[0];
    auto vec   = buffer<double>(ml_vec);
    auto ipix  = factory.createBuffer<int64_t>(npix);
    auto x = &vec[0];
    auto y = &vec[npix];
    auto z = &vec[2*npix];

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

DISPATCH_FN(map2alm_iter) {
    healpix base = nsideorder(inputs[1], inputs[2]);
    auto [buf_map, len_map] = bufferlen<double>(inputs[3]);
    auto lmax = scalar<int32_t>(inputs[4]);
    auto mmax = scalar<int32_t>(inputs[5]);
    auto [buf_wght, len_wght] = bufferlen<double>(inputs[6]);
    auto iter = scalar<int32_t>(inputs[7]);

    auto map = healmap();
    {
        arr<double> tmp(buf_map.get(), len_map);
        map.Set(tmp, base.Scheme());
    }

    arr<double> rwghts(buf_wght.get(), len_wght);

    auto nalms = Alm_Base::Num_Alms(lmax, mmax);
    auto alms = healalm();
    auto buf_alms = factory.createBuffer<complex64>(nalms);
    {
        arr<complex64> tmp(buf_alms.get(), nalms);
        alms.Set(tmp, lmax, mmax);
    }

    map2alm_iter(map, alms, iter, rwghts);

    outputs[0] = factory.createArrayFromBuffer({nalms}, move(buf_alms));
}

DISPATCH_FN(map2alm_pol_iter) {
    healpix base = nsideorder(inputs[1], inputs[2]);
    auto [buf_mapT, len_mapT] = bufferlen<double>(inputs[3]);
    auto [buf_mapQ, len_mapQ] = bufferlen<double>(inputs[4]);
    auto [buf_mapU, len_mapU] = bufferlen<double>(inputs[5]);
    auto lmax = scalar<int32_t>(inputs[6]);
    auto mmax = scalar<int32_t>(inputs[7]);
    auto [buf_wght, len_wght] = bufferlen<double>(inputs[8]);
    auto iter = scalar<int32_t>(inputs[9]);

    auto mapT = healmap();
    auto mapQ = healmap();
    auto mapU = healmap();
    {
        arr<double> tmp(buf_mapT.get(), len_mapT);
        mapT.Set(tmp, base.Scheme());
    }
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
    auto almsT = healalm();
    auto almsG = healalm();
    auto almsC = healalm();

    auto buf_almsT = factory.createBuffer<complex64>(nalms);
    auto buf_almsG = factory.createBuffer<complex64>(nalms);
    auto buf_almsC = factory.createBuffer<complex64>(nalms);
    {
        arr<complex64> tmp(buf_almsT.get(), nalms);
        almsT.Set(tmp, lmax, mmax);
    }
    {
        arr<complex64> tmp(buf_almsG.get(), nalms);
        almsG.Set(tmp, lmax, mmax);
    }
    {
        arr<complex64> tmp(buf_almsC.get(), nalms);
        almsC.Set(tmp, lmax, mmax);
    }

    map2alm_pol_iter(mapT, mapQ, mapU,
            almsT, almsG, almsC, iter, rwghts);

    outputs[0] = factory.createArrayFromBuffer({nalms}, move(buf_almsT));
    outputs[1] = factory.createArrayFromBuffer({nalms}, move(buf_almsG));
    outputs[2] = factory.createArrayFromBuffer({nalms}, move(buf_almsC));
}

DISPATCH_FN(alm2map) {
    auto lmax = scalar<int32_t>(inputs[1]);
    auto mmax = scalar<int32_t>(inputs[2]);
    auto [buf_alms, len_alms] = bufferlen<complex64>(inputs[3]);
    healpix base = nsideorder(inputs[4], inputs[5]);

    auto alms = healalm();
    {
        arr<complex64> tmp(buf_alms.get(), len_alms);
        alms.Set(tmp, lmax, mmax);
    }

    auto npix = 12 * base.Nside() * base.Nside();
    auto map = healmap();
    auto buf_map = factory.createBuffer<double>(npix);
    {
        arr<double> tmp(buf_map.get(), npix);
        map.Set(tmp, base.Scheme());
    }

    alm2map(alms, map);

    outputs[0] = factory.createArrayFromBuffer({npix}, move(buf_map));
}

DISPATCH_FN(alm2map_pol) {
    auto lmax = scalar<int32_t>(inputs[1]);
    auto mmax = scalar<int32_t>(inputs[2]);
    auto [buf_almsT, len_almsT] = bufferlen<complex64>(inputs[3]);
    auto [buf_almsG, len_almsG] = bufferlen<complex64>(inputs[4]);
    auto [buf_almsC, len_almsC] = bufferlen<complex64>(inputs[5]);
    healpix base = nsideorder(inputs[6], inputs[7]);

    auto almsT = healalm();
    auto almsG = healalm();
    auto almsC = healalm();
    {
        arr<complex64> tmp(buf_almsT.get(), len_almsT);
        almsT.Set(tmp, lmax, mmax);
    }
    {
        arr<complex64> tmp(buf_almsG.get(), len_almsG);
        almsG.Set(tmp, lmax, mmax);
    }
    {
        arr<complex64> tmp(buf_almsC.get(), len_almsC);
        almsC.Set(tmp, lmax, mmax);
    }

    auto npix = 12 * base.Nside() * base.Nside();
    auto mapT = healmap();
    auto mapQ = healmap();
    auto mapU = healmap();
    auto buf_mapT = factory.createBuffer<double>(npix);
    auto buf_mapQ = factory.createBuffer<double>(npix);
    auto buf_mapU = factory.createBuffer<double>(npix);
    {
        arr<double> tmp(buf_mapT.get(), npix);
        mapT.Set(tmp, base.Scheme());
    }
    {
        arr<double> tmp(buf_mapQ.get(), npix);
        mapQ.Set(tmp, base.Scheme());
    }
    {
        arr<double> tmp(buf_mapU.get(), npix);
        mapU.Set(tmp, base.Scheme());
    }

    alm2map_pol(almsT, almsG, almsC, mapT, mapQ, mapU);

    outputs[0] = factory.createArrayFromBuffer({npix}, move(buf_mapT));
    outputs[1] = factory.createArrayFromBuffer({npix}, move(buf_mapQ));
    outputs[2] = factory.createArrayFromBuffer({npix}, move(buf_mapU));
}

DISPATCH_FN(alm2cl) {
    auto lmax = scalar<int32_t>(inputs[1]);
    auto mmax = scalar<int32_t>(inputs[2]);
    auto [buf_alms1, len_alms1] = bufferlen<complex64>(inputs[3]);
    auto [buf_alms2, len_alms2] = bufferlen<complex64>(inputs[3]);

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

    auto powspec = factory.createBuffer<double>(lmax + 1);

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

    outputs[0] = factory.createArrayFromBuffer({lmax+1}, move(powspec));
}
