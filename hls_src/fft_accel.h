#include "ap_axi_sdata.h"
#include "hls_stream.h"
#include <stdint.h>
#include <math.h>
#include <cstdlib>
#include <exception>
#include <memory>
#include "coef_init.h"
/* ****************************** DEFINES ************************************** */

#define REVERSAL
#define NPOINTS  		32
#define FFTRADIX     	5
#define FFTRAD_1     	6

#define POW2(casc)      ((1) << casc)
#define W_IDX(idx, casc)(((idx) % POW2(casc)) * (POW2(FFTRADIX - 1 - casc)))

typedef ap_axiu<32, 4, 1, 1> stream_1ch;

/* ****************************** STRUCTURES ************************************ */
template <typename T>
struct cmpx_t
{
	T re_;
	T im_;
};

union
{
	uint32_t 		uint;
	cmpx_t<int16_t> cmpx;

}uint2cmpx;

/* ****************************** FUNCTIONS ************************************** */


template <typename T, typename U> 		        U sum_pair(cmpx_t<T> x_0){return ((U)x_0.im_ + (U)x_0.re_);}
template <typename T, typename U> 		        U sub_pair(cmpx_t<T> x_0){return ((U)x_0.im_ - (U)x_0.re_);}
template <typename T, typename U> 		 cmpx_t<T>cnv_pair(cmpx_t<U> x_0){return {(T)x_0.re_ , (T)x_0.im_};}
template <typename T, typename U, int F> cmpx_t<T>scl_pair(cmpx_t<U> x_0){return {(T)(x_0.re_ >> F), (T)(x_0.im_ >> F)};}
template <typename T>cmpx_t<T> operator+(const cmpx_t<T> & l,const cmpx_t<T> & r) {return {l.re_ + r.re_,l.im_ + r.im_};}
template <typename T>cmpx_t<T> operator-(const cmpx_t<T> & l,const cmpx_t<T> & r) {return {l.re_ - r.re_,l.im_ - r.im_};}

using namespace hls;

/* ****************************** FUNCTIONS DECLARATION *************************** */

void FFT_TOP  	(stream<stream_1ch> &in_stream, stream<stream_1ch> &out_stream);

/* ****************************** C++ TEMPLATES ************************************ */

/**
 *  butterfly dit (decimation-in-time) implementation
 *
 *  @param x_0 point in complex value converted in word cmpx_t<int16_t> -> uint32_t
 *  @param y_0 point in complex value converted in word cmpx_t<int16_t> -> uint32_t
 *  @param w_0 coeff in complex value converted in word cmpx_t<int16_t> -> uint32_t
 *  @return x_1, y_1 points in complex value converted in word cmpx_t<int16_t> -> uint32_t
 *           ------             -----
 *  x_0---->| z^-2 |--+-----+->|  +  |------> x_1
 *           ------    \   /    -----
 *                      \ /
 *                       X
 *           ------    /  \     -----
 *  y_0---->| cmpx |--+----+-->|  -  |------> y_1
 *           -+----             -----
 *           /
 *  w_0---->
 *
 *   x_1 = x_0 - y_0 * w_0
 *   y_1 = x_0 + y_0 * w_0
 */
template <typename T, typename U, typename V, uint8_t F>
void butter_dit(T x0, T y0, T w0, T *x1, T *y1)
{
#pragma HLS INLINE
	cmpx_t<V>   x_0 = {0, 0}, y_0 = {0, 0}, w_0 = {0, 0};

	uint2cmpx.uint = x0;
	x_0 = uint2cmpx.cmpx;

	uint2cmpx.uint = y0;
	y_0 = uint2cmpx.cmpx;

	uint2cmpx.uint = w0;
	w_0 = uint2cmpx.cmpx;

	cmpx_t<V>  x_1 = {0, 0};
	cmpx_t<V>  y_1 = {0, 0};

	cmpx_t<U> cmpx_mlt = {0, 0};
	cmpx_mlt.re_  = sum_pair<V, U>(y_0) * (U)w_0.re_ - sum_pair<V, U>(w_0) * (U)y_0.im_;
	cmpx_mlt.im_  = sum_pair<V, U>(y_0) * (U)w_0.re_ + sub_pair<V, U>(w_0) * (U)y_0.re_;

	cmpx_t<V> scaled_mlt = scl_pair<V, U, F>(cmpx_mlt);

	cmpx_t<U> dout_0 = cnv_pair<U, V>(x_0) + cnv_pair<U, V>(scaled_mlt);
	cmpx_t<U> dout_1 = cnv_pair<U, V>(x_0) - cnv_pair<U, V>(scaled_mlt);

	x_1 = scl_pair<V, U, 1>(dout_0);
	y_1 = scl_pair<V, U, 1>(dout_1);

	*x1 = (T &) x_1;
	*y1 = (T &) y_1;
} // butter_dit

/**
 * Bit reversal operation
 *
 * @param 	Addr 		U-bits address for reversal
 * @return 	reversal 	address
 *
 */
template <typename T>
T revBits(T Addr)
{
#pragma HLS INLINE
	T revAddr = 0;
	rb_L:for(uint8_t idx = 0; idx < FFTRADIX; idx ++)
	{
		revAddr <<= 1;
		revAddr |= Addr & 1;
		Addr >>= 1;
	} // for

	return revAddr;
} // revBits

template <typename T>
void reverse_stage(T x[NPOINTS], T y[NPOINTS])
{
	T temp = 0;
	uint16_t idx_r = 0;

	revst_L:for(uint16_t idx_d = 0;  idx_d < NPOINTS; idx_d ++)
	{
		idx_r = revBits<uint16_t>(idx_d);
		y[idx_r] = x[idx_d];
	} // for
} // reverse_stage



/*
 * 	reading a sample from axis interface
 */
template <typename T, int U, int TI, int TD>
T read_stream(ap_axiu <sizeof(T)*8,U,TI,TD> const &e)
{
#pragma HLS INLINE
	union
	{
		int ival;
		T oval;
	} converter;
	converter.ival = e.data;
	T ret = converter.oval;

	volatile ap_uint<sizeof(T)> strb = e.strb;
	volatile ap_uint<sizeof(T)> keep = e.keep;
	volatile ap_uint<U> user = e.user;
	volatile ap_uint<1> last = e.last;
	volatile ap_uint<TI> id = e.id;
	volatile ap_uint<TD> dest = e.dest;

	return ret;
} // read_stream

template <typename T, int U, int TI, int TD, int NPTS>
void pop_input(stream<stream_1ch> &in_stream,T y[NPTS])
{
#pragma HLS INLINE
	pin_L:for(uint16_t idx = 0; idx <  NPTS; idx ++)
		y[idx] = read_stream<T, U, TI, TD>(in_stream.read());
} // pop_input


/*
 * 	writing a sample to axis interface
 */
template <typename T, int U, int TI, int TD>
ap_axiu <sizeof(T)*8,U,TI,TD> write_stream(T const &v, bool last = false)
{
#pragma HLS INLINE
	ap_axiu<sizeof(T)*8,U,TI,TD> e;

	union
	{
		int oval;
		T ival;
	} converter;
	converter.ival = v;
	e.data = converter.oval;

	// set it to sizeof(T) ones
	e.strb = -1;
	e.keep = 15; //e.strb;
	e.user = 0;
	e.last = last ? 1 : 0;
	e.id = 0;
	e.dest = 0;
	return e;
}

template <typename T, int U, int TI, int TD, int NPTS>
void push_output(stream<stream_1ch> &out_stream,T y[NPTS])
{
#pragma HLS INLINE
	pout_L:for(uint16_t idx = 0; idx <  NPTS; idx ++)
		out_stream.write(write_stream<T, U, TI, TD>(y[idx], (idx == NPTS - 1)));
} // push_output

/**
 *
 * FFT stage (reading array + butterfly processing)
 *
 * @param x[] 	input array
 * @param w[] 	input coefficients
 * @param casc 	stage number
 *
 * @return y[] 	output array
 *
 */

template <typename T, typename U, typename V>
void n_stage(T x[NPOINTS], T y[NPOINTS], uint8_t casc)
{
	T x0 = 0, y0 = 0;
	T x1 = 0, y1 = 0;
	nstage_L:for(uint16_t idx = 0; idx < NPOINTS / 2; idx ++)
	{
#pragma HLS PIPELINE
		uint16_t  d 	 = ((idx % POW2(casc)) == 0)  ? 2*idx : d + 1;
		uint16_t _idx1   = d + 0;
		uint16_t _idx2   = d + POW2(casc);
		T w0 = wcoe[W_IDX(idx, casc)];

		x0 		 = x[_idx1];
		y0 		 = x[_idx2];
		butter_dit<T, U, V, 15>(x0, y0, w0, &x1, &y1);
		y[_idx1] = x1;
		y[_idx2] = y1;
	} // nstage_L
} // n_stage



/**
 *  FFT CORE
 *
 */
template <typename T, typename U, typename V, int TU, int TI, int TD>
void wrapped_fft_hw (stream<stream_1ch> &in_stream, stream<stream_1ch> &out_stream)
{
#pragma HLS DATAFLOW
	T     mem_bram[FFTRAD_1][NPOINTS];
#pragma HLS ARRAY_PARTITION variable=mem_bram dim=1 type=block factor=6
#pragma HLS BIND_STORAGE variable=mem_bram type=ram_t2p impl=uram
	T x[NPOINTS];
#pragma HLS BIND_STORAGE variable=x type=ram_t2p impl=uram


	pop_input  <T, TU, TI, TD, NPOINTS>( in_stream, x);
	reverse_stage<T>(x, mem_bram[0]);
	ffthw_L:for(uint8_t casc = 0; casc < FFTRADIX; casc ++)
#pragma HLS UNROLL
			n_stage    <T,U,V>(  mem_bram[casc], mem_bram[casc + 1], casc);

	push_output<T, TU, TI, TD, NPOINTS>(out_stream, mem_bram[FFTRADIX]);
} // wrapped_fft_hw
