// #include "fft_accel.h"
#include <iostream>
#include <fstream>
#include "testbenches.h"

#define TEST_BUT
#define TEST_TOP

using namespace std;


int main()
{

#ifdef TEST_BUT
	test_butterfly();
#endif


#ifdef TEST_TOP
	test_fft();
#endif


return 0;
}




