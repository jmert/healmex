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
    id_pix2vec      = 1,
    id_pix2zphi     = 2,
    id_pix2ang      = 3,
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

        switch (dispatch) {
            case id_heartbeat:
                outputs[0] = factory.createScalar(true);
                break;

            case id_pix2ang:
                CHECK_NINOUT("pix2ang", 2, 2);
                CHECK_INPUT_SCALAR("pix2ang", "nside", 1);
                CHECK_INPUT_INT64("pix2ang", "nside", 1);
                CHECK_INPUT_INT64("pix2ang", "ipix", 2);
                mex_pix2ang(outputs, inputs);
                break;

            default:
                error("Unhandled dispatch type %d", dispatch);
        }
    }

private:
    std::shared_ptr<matlab::engine::MATLABEngine> engine = getEngine();
    ArrayFactory factory;

    void mex_pix2ang(ArgumentList& outputs, ArgumentList& inputs);

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

/* Externally callable function implementations */

void MexFunction::mex_pix2ang(ArgumentList& outputs, ArgumentList& inputs) {
    TypedArray<int64_t> ml_nside = inputs[1];
    TypedArray<int64_t> ml_ipix  = inputs[2];
    healpix base;

    int64_t nside = ml_nside[0];
    base.SetNside(nside, RING);

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

