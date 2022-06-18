// #include "fft_accel.h"
#include <iostream>
#include <fstream>
#include "testbenches.h"

// #define TEST_BUT
#define TEST_TOP

using namespace std;


int main()
{

#ifdef TEST_BUT
	test_but();
#endif


#ifdef TEST_TOP
	test_top();
#endif

#ifdef MAIN

	stream<stream_1ch> 	src, dst;
	cout << "Initialization input" << endl;
	string line;

	cmpx16_t din_cmpx[NPOINTS];
	uint32_t din_uint[NPOINTS];

	uint16_t idx = 0;

	ifstream fpReIn (DIN_RE);
	if (fpReIn.is_open())
	{
		idx = 0;
		while ( !fpReIn.eof())
		{
			fpReIn >> din_cmpx[idx].re_;
			idx++;
		}
		fpReIn.close();
	}

	else cout << "Unable to open file" << endl;

	ifstream fpImIn (DIN_IM);
	if (fpImIn.is_open())
	{
		idx = 0;
		while ( !fpImIn.eof())
		{
			fpImIn >> din_cmpx[idx].im_;
			idx++;
		}
		fpImIn.close();
	}

	else cout << "Unable to open file" << endl;

	memcpy(&din_uint, &din_cmpx, NPOINTS*sizeof(uint32_t));
	for(uint16_t p = 0; p < NPOINTS; p ++)
	{
		cout << "d_in_re["<< p << "] = " << din_cmpx[p].re_  << " d_in_im["<< p << "] = " << din_cmpx[p].im_ << endl;
	}

	for(uint16_t points = 0; points < NPOINTS; points ++)
	{
		stream_1ch VallIn;
		VallIn.data = din_uint[points % NPOINTS];
		src << VallIn;
		// cout << "d_in["<< points << "] = 0x" << (uint32_t)(VallIn.data >> 0) << endl;
	}
	cout << "Initialization input done" << endl;

	cout << "Start Function" << endl;
	fft_accel(src, dst);
	cout << "End Function" << endl;

	cout << "Reading output" << endl;
	for(uint16_t points = 0; points < NPOINTS; points ++)
	{
		stream_1ch tmp_1;
		uint32_t tmp_0 = 0;
		dst.read(tmp_1);

		tmp_0 = tmp_1.data;
		// cout << "tmp_1 = " << (uint32_t)(tmp_0) << endl;
		// cout << "d_out["<< points << "] = 0x" << (uint32_t)(tmp_0) << endl;
		cout << "x_re = " << (int16_t)(tmp_0 >> 0) << " x_im = " << (int16_t)(tmp_0 >> 16) <<endl;
	}

#endif



return 0;
}




