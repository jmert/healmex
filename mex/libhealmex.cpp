#include <mex.hpp>
#include <mexAdapter.hpp>

#include <complex>
#include <utility> // tuple
#include <vector>

#include <healpix_map.h>

using namespace std;
using namespace matlab::data;
using matlab::mex::ArgumentList;

using healpix = Healpix_Base2;
using healmap = Healpix_Map<double>;

/* Useful type predicates */

inline bool ischar(Array& a)      { return a.getType() == ArrayType::CHAR; }
inline bool isint64(Array& a)     { return a.getType() == ArrayType::INT64; }
inline bool isdouble(Array& a)    { return a.getType() == ArrayType::DOUBLE; }
inline bool iscomplex64(Array& a) { return a.getType() == ArrayType::COMPLEX_DOUBLE; }

inline bool isscalar(Array& a) { return a.getNumberOfElements() == 1; }

template <typename T, typename A> buffer_ptr_t<T> buffer(A array);

/*
 * MEX entry point
 */

enum libhealpix_mex_calls {
    id_heartbeat    = -1,
    id_nest2ring    =  1,
    id_ring2nest    =  2,
    id_pix2vec      = 11,
    id_pix2zphi     = 12,
    id_pix2ang      = 13,
    id_vec2pix      = 14,
    id_zphi2pix     = 15,
    id_ang2pix      = 16,
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
        auto ml_buf = buffer<complex<double>>(move(ml_arr));
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
buffer_ptr_t<T> buffer(A array) {
    TypedArray<T> array_t = array;
    return array_t.release();
}

// Turns a Matlab array into a healpix arr<T> object
template <typename T>
tuple<arr<T>, buffer_ptr_t<T>> array2arrbuf(Array& ml_array)
{
    auto nel = ml_array.getNumberOfElements();
    auto buf = buffer<T>(ml_array);
    auto vec = arr<T>(buf.get(), nel);
    return tuple{move(vec), move(buf)};
}

// Turns a Matlab array into a healpix map object.
tuple<healmap, buffer_ptr_t<double>> array2mapbuf(Array& ml_map)
{
    auto [hpx, buf] = array2arrbuf<double>(ml_map);
    auto map = healmap(hpx, RING);
    return tuple{move(map), move(buf)};
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
