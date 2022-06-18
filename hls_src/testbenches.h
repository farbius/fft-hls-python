#include "fft_accel.h"


/* ****************************** DEFINES ************************************** */

// #define CONSOLE // output to console (for few FFT points)
#define DIN_RE 	"..\\..\\..\\..\\..\\sim_files\\data_re.txt"
#define DIN_IM 	"..\\..\\..\\..\\..\\sim_files\\data_im.txt"
#define DOUT    "..\\..\\..\\..\\..\\sim_files\\dataOUT.txt"

using namespace std;

/* ****************************** C++ TEMPLATES ************************************** */

template <typename T, typename U, uint8_t F>
void but_dif_sw(cmpx_t<T> x_0, cmpx_t<T> y_0, cmpx_t<T> w_0, cmpx_t<T>* x_1, cmpx_t<T>* y_1)
{
	cmpx_t<U> cmpx_mlt = {0, 0};
	cmpx_mlt.re_  = sum_pair<T, U>(y_0) * (U)w_0.re_ - sum_pair<T, U>(w_0) * (U)y_0.im_;
	cmpx_mlt.im_  = sum_pair<T, U>(y_0) * (U)w_0.re_ + sub_pair<T, U>(w_0) * (U)y_0.re_;

	cmpx_t<T> scaled_mlt = scl_pair<T, U, F>(cmpx_mlt);

	cmpx_t<U> dout_0 = cnv_pair<U, T>(x_0) + cnv_pair<U, T>(scaled_mlt);
	cmpx_t<U> dout_1 = cnv_pair<U, T>(x_0) - cnv_pair<U, T>(scaled_mlt);


	*x_1 = scl_pair<T, U, 1>(dout_0);
	*y_1 = scl_pair<T, U, 1>(dout_1);

}

void read_txt(uint32_t *dout)
{
	// array
	std::cout << "READING INPUT" << std::endl;
	cmpx_t<int16_t> _dout[NPOINTS];
	uint16_t idx = 0;

	ifstream fpReIn (DIN_RE);
	ifstream fpImIn (DIN_IM);

	if (fpReIn.is_open())
	{
		idx = 0;
		while ( !fpReIn.eof())
		{
			fpReIn >> _dout[idx].re_;
			idx++;
		}
		fpReIn.close();
	}

	else cout << "Unable to open file" << endl;


	if (fpImIn.is_open())
	{
		idx = 0;
		while ( !fpImIn.eof())
		{
			fpImIn >> _dout[idx].im_;
			idx++;
		}
		fpImIn.close();
	}

	else std::cout << "Unable to open file" << std::endl;

	memcpy(&dout[0], &_dout, NPOINTS*sizeof(uint32_t));

	std::cout  << "DONE!" << std::endl;

}

void test_top()
{
	cout  << endl << "START SIMULATION" << endl;

	uint32_t din_uint[NPOINTS];
	stream<stream_1ch> 	src, dst;

	read_txt(&din_uint[0]);

#ifdef CONSOLE
	cmpx_t<int16_t> din_points = {0, 0};
	cout << endl << "Input Points " << endl;
	cout << setw(3) << "NUM" << "\t";
	cout << setw(6) << "REAL" << "\t";
	cout << setw(6) << "IMAG" << "\n";
	for(int idx = 0; idx < NPOINTS; idx++)
	{
		uint2cmpx.uint = din_uint[idx];
		din_points     = uint2cmpx.cmpx;

		cout << setw(3) << idx << ":" << "\t";
		cout << setw(6) << din_points.re_ << "\t";
		cout << setw(6) << din_points.im_ << "\n";

	}
#endif

	for(unsigned int idx = 0; idx < NPOINTS; idx ++)
	{
	   stream_1ch VallIn;
	   VallIn.data = din_uint[idx];
	   VallIn.last = (idx == NPOINTS - 1);
	   src << VallIn;
	}


	cout  << endl << "RUN HARDWARE TEST FOR " << NPOINTS << " FFT points"  << endl;
	FFT_TOP(src, dst);
	cout  << "DONE!" << endl;

	cmpx_t<int16_t> dout_points[NPOINTS];

	for(unsigned int idx = 0; idx < NPOINTS; idx ++)
	{
		stream_1ch VallOut;
		dst >> VallOut;
		uint2cmpx.uint = VallOut.data;
		dout_points[idx].re_ = uint2cmpx.cmpx.re_;
		dout_points[idx].im_ = uint2cmpx.cmpx.im_;
	}

#ifdef CONSOLE
	std::cout << endl << "Output Points " << endl;
	cout << setw(3) << "NUM" << "\t";
	cout << setw(6) << "REAL" << "\t";
	cout << setw(6) << "IMAG" << "\n";
	for(int idx = 0; idx < NPOINTS; idx++)
	{
		cout << setw(3) << idx << ":" << "\t";
		cout << setw(6) << dout_points[idx].re_ << "\t";
		cout << setw(6) << dout_points[idx].im_ << "\n";
	}


#else
	cout << endl << "WRITING RESULT TO TXT FILE" << endl;

	ofstream myfile;
	myfile.open (DOUT);
	for(int k = 0; k < NPOINTS; k++)
	{
		myfile << dout_points[k].re_  << "\n";
		myfile << dout_points[k].im_  << "\n";
	}
	myfile.close();
	cout  << "DONE!" << endl;

#endif


	cout << endl << "END SIMULATION" << endl << endl;
}

/*
 *
 *  TESTING BUTTERFLY
 *
 */
/*void test_but()
{

	cout << "START SIMULATION" << endl;

	pair<int16_t, int16_t> x_0[NPAIRS], y_0[NPAIRS];
	pair<int16_t, int16_t> x_1[NPAIRS], y_1[NPAIRS];
	pair<int16_t, int16_t> w_0;

	w_0  = make_pair((int16_t)rand() % 32767, (int16_t)rand() % 32767);

	cout << " INITIALIZATION " << endl;
	for (unsigned int k = 0; k < NPAIRS; k ++ )
	{
		x_0[k] = make_pair((int16_t)rand() % 32767, (int16_t)rand() % 32767);
		y_0[k] = make_pair((int16_t)rand() % 32767, (int16_t)rand() % 32767);

		x_1[k] = make_pair(0, 0);
		y_1[k] = make_pair(0, 0);
	}

	cout << " SOFTWARE " << endl;
	for (unsigned int k = 0; k < NPAIRS; k ++ )
		but_dif_sw<int16_t, int32_t, 15>(x_0[k], y_0[k], w_0, &x_1[k], &y_1[k]);


	uint32_t *x0_32[NPAIRS];
	uint32_t *y0_32[NPAIRS];
	uint32_t *w0_32 = (uint32_t *) &w_0;

	uint32_t x1_32[NPAIRS];
    uint32_t y1_32[NPAIRS];

    pair<int16_t, int16_t> x_2[NPAIRS], y_2[NPAIRS];

    for (unsigned int k = 0; k < NPAIRS; k ++ )
    {
    	 x0_32[k] = (uint32_t *) &x_0[k];
    	 y0_32[k] = (uint32_t *) &y_0[k];
    	 x1_32[k] = 0;
    	 y1_32[k] = 0;
    }

    cmpx16_t *x0_cmpx[NPAIRS];
    cmpx16_t *y0_cmpx[NPAIRS];
    cmpx16_t *w0_cmpx = (cmpx16_t *) &w_0;

    cmpx16_t x1_cmpx[NPAIRS];
    cmpx16_t y1_cmpx[NPAIRS];

	for (unsigned int k = 0; k < NPAIRS; k ++ )
	{
		 x0_cmpx[k] = (cmpx16_t *) &x_0[k];
		 y0_cmpx[k] = (cmpx16_t *) &y_0[k];
		 x1_cmpx[k] = {0, 0};
		 y1_cmpx[k] = {0, 0};
	}

    cout << " HARDWARE " << endl;
    for (unsigned int k = 0; k < NPAIRS; k ++ )
    	butter_dit(x0_32[k], y0_32[k], w0_32, &x1_32[k], &y1_32[k]);
    	// butter_dit_0<int16_t, int32_t, 15>(x_0[k], y_0[k], w_0[k], &x_2[k], &y_2[k]);

    uint32_t *x0_G[NPAIRS];
	uint32_t *y0_G[NPAIRS];
	for (unsigned int k = 0; k < NPAIRS; k ++ )
	{
		 x0_G[k] = (uint32_t *) &x_1[k];
		 y0_G[k] = (uint32_t *) &y_1[k];
	}



    for (unsigned int k = 0; k < NPAIRS; k ++ )
	{
		if(*x0_G[k] != x1_32[k])
		{
			cout << " x error " << endl;
			cout << "NOT PASSED " << endl;
			exit(0);
		}

		if(*y0_G[k] != y1_32[k])
		{
			cout << " y error " << endl;
			cout << "NOT PASSED " << endl;
			exit(0);
		}

	}
		cout << "PASSED " << endl;





	//
	cout << "END SIMULATION" << endl;

}*/
