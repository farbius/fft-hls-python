#include "fft_accel.h"





void FFT_TOP(stream<stream_1ch> &in_stream, stream<stream_1ch> &out_stream)
{
#pragma HLS DATAFLOW
#pragma HLS INTERFACE ap_ctrl_none port=return
#pragma HLS INTERFACE axis port=in_stream
#pragma HLS INTERFACE axis port=out_stream

	wrapped_fft_hw <uint32_t, int32_t, int16_t, 4, 1, 1>(in_stream, out_stream);
}

/*

void NSTAGE_TOP(stream<stream_1ch> &in_stream, stream<stream_1ch> &out_stream, uint8_t casc)
{
#pragma HLS DATAFLOW
#pragma HLS INTERFACE ap_ctrl_none port=return
#pragma HLS INTERFACE axis port=in_stream
#pragma HLS INTERFACE axis port=out_stream

	wrapped_n_stage <uint32_t, int32_t, int16_t, 4, 1, 1>(in_stream, out_stream);
}
*/
