#include "fft_sw.h"

/*

void butter_dit(cmpx16_t x_0, cmpx16_t y_0, cmpx16_t w_0, cmpx16_t *x_1, cmpx16_t *y_1);
void  coef_init(cmpx16_t *w);
void 	n_stage(cmpx16_t *x, cmpx16_t *y, cmpx16_t *w, uint8_t casc);

uint16_t revBits(uint16_t Addr, uint8_t Nbytes);
*/




int main()
{
    printf("<- Start App \n");
    struct timeval t_0, t_1, t_diff;
    size_t         t_usec = 0;

    cmpx16_t w[NFFT];

    cmpx16_t x0[NFFT];
    cmpx16_t x1[NFFT];




    FILE *fp_re, *fp_im, *fp_out;


    fp_re = fopen("data_re.txt", "r");
    fp_im = fopen("data_im.txt", "r");
    if(fp_re == NULL)    	printf("Error \n");
    if(fp_im == NULL)    	printf("Error \n");

    for(uint16_t idx = 0; idx < NFFT; idx ++)
    {
    	fscanf(fp_re, "%d", &x0[idx].re_);
    	fscanf(fp_im, "%d", &x0[idx].im_);
    }

    fclose(fp_re);
    fclose(fp_im);
   /* for(uint16_t idx = 0; idx < NFFT; idx ++)
   {
	printf("x.re = %d x.im = %d \n", x0[idx].re_, x0[idx].im_);
   }*/
    printf("<- Initialization of coefficients \n");
    coef_init(&w[0]);


    printf("<- Software FFT \n");
    gettimeofday(&t_0, NULL);
    /* Function for evaluation */

    fft_top(&x0[0], &x1[0], &w[0]);

    /*for(uint16_t idx = 0; idx < NFFT; idx ++)
   {
	printf("y.re = %3d    y.im = %3d \n", x1[idx].re_, x1[idx].im_);
   }*/

    /*                         */
    gettimeofday(&t_1, NULL);
    timersub(&t_1, &t_0, &t_diff);
    t_usec	= t_diff.tv_sec * 1000000 + t_diff.tv_usec;
    printf("<- Software FFT run time: %ld us\n", t_usec );

    fp_out = fopen("dataOUT.txt", "w");
    for(uint16_t idx = 0; idx < NFFT; idx ++)
	{
    	fprintf(fp_out, "%d \n", x1[idx].re_);
    	fprintf(fp_out, "%d \n", x1[idx].im_);
	}

    fclose(fp_out);
    printf("<- End App \n");
    return 0;
}


