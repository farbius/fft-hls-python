#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <math.h>

#define AMPL 				(32768)
#define NFFT 				(8)
#define NRDX 				(uint8_t)(log2(NFFT))
#define POW2(casc)      	((1) << (casc))
#define W_IDX(idx, casc)	(((idx) % POW2(casc)) * (POW2(NRDX - 1 - casc)))

typedef struct
{
	int16_t re_;
	int16_t im_;
} cmpx16_t;

typedef struct
{
	int32_t re_;
	int32_t im_;
} cmpx32_t;

inline int32_t 	sum_cmpx(cmpx16_t x)	{return ((int32_t)x.im_ + (int32_t)x.re_);};
inline int32_t 	sub_cmpx(cmpx16_t x)	{return ((int32_t)x.im_ - (int32_t)x.re_);};
inline cmpx32_t cnv_cmpx(cmpx16_t x)	{
	cmpx32_t y;
	y.re_ = (int32_t)x.re_;
	y.im_ = (int32_t)x.im_;
	return y;
};

void butter_dit(cmpx16_t x_0, cmpx16_t y_0, cmpx16_t w_0, cmpx16_t *x_1, cmpx16_t *y_1);
void  coef_init(cmpx16_t *w);
void 	n_stage(cmpx16_t *x, cmpx16_t *y, cmpx16_t *w, uint8_t casc);

uint16_t revBits(uint16_t Addr, uint8_t Nbytes);




int main()
{
    printf("<- Start App \n");
    struct timeval t_0, t_1, t_diff;
    size_t         t_usec = 0;

    cmpx16_t w[NFFT];
    coef_init(&w[0]);
   /* for(uint16_t idx = 0; idx < NFFT; idx ++)
    {
    	printf("w.re = %d w.im = %d \n", w[idx].re_, w[idx].im_);
    }*/

    cmpx16_t x0[NFFT];
    cmpx16_t x1[NFFT];
    cmpx16_t x2[NFFT];


    FILE *fp_re, *fp_im;


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
    for(uint16_t idx = 0; idx < NFFT; idx ++)
   {
	printf("x.re = %d x.im = %d \n", x0[idx].re_, x0[idx].im_);
   }



    printf("\n");
    printf("\n");

    for(uint16_t casc = 0; casc < NRDX; casc ++)
        for(uint16_t idx = 0; idx < NFFT / 2; idx ++)
        	printf("W_IDX = %d \n", W_IDX(idx, casc));


    printf("\n");
    printf("\n");

    gettimeofday(&t_0, NULL);
    /* Function for evaluation */
    n_stage(&x0[0], &x1[0], &w[0], 0);
    n_stage(&x1[0], &x0[0], &w[0], 1);
    n_stage(&x0[0], &x1[0], &w[0], 2);


    for(uint16_t idx = 0; idx < NFFT; idx ++)
   {
	printf("y.re = %3d    y.im = %3d \n", x1[idx].re_, x1[idx].im_);
   }

    /*                         */
    gettimeofday(&t_1, NULL);
    timersub(&t_1, &t_0, &t_diff);
    t_usec	= t_diff.tv_sec * 1000000 + t_diff.tv_usec;
    printf("<- run time: %ld us\n", t_usec );

    printf("<- End App \n");
    return 0;
}


uint16_t revBits(uint16_t Addr, uint8_t Nbytes)
{
	uint16_t revAddr = 0;
	for(uint8_t idx = 0; idx < Nbytes; idx ++)
	{
		revAddr <<= 1;
		revAddr |= Addr & 1;
		Addr >>= 1;
	}
	return revAddr;
}

void n_stage(cmpx16_t *x, cmpx16_t *y, cmpx16_t *w, uint8_t casc)
{
	cmpx16_t x0, y0;
	cmpx16_t x1, y1;
	uint16_t 	d = 0;
	for(uint16_t idx = 0; idx < NFFT / 2; idx ++)
	{
		uint16_t _widx 	= W_IDX(idx, casc);
		        	d 	= ((idx % POW2(casc)) == 0)  ? 2*idx : d + 1;
		uint16_t _idx1  = revBits((d + 0), 			(uint8_t)NRDX);
		uint16_t _idx2  = revBits((d + POW2(casc)),	(uint8_t)NRDX);
		printf("_wdx = %d, idx1 = %d, idx2 = %d \n", _widx, _idx1, _idx2);

		x0 = x[_idx1];
		y0 = x[_idx2];
		butter_dit(x0, y0, w[_widx], &x1, &y1);

		y[_idx1] = x1;
		y[_idx2] = y1;
	}

}

void butter_dit(cmpx16_t x_0, cmpx16_t y_0, cmpx16_t w_0, cmpx16_t *x_1, cmpx16_t *y_1)
{

	int32_t  ywr 	= sum_cmpx(y_0)*(int32_t)w_0.re_ - sum_cmpx(w_0)*(int32_t)y_0.im_;
	int32_t  ywi 	= sum_cmpx(y_0)*(int32_t)w_0.re_ + sub_cmpx(w_0)*(int32_t)y_0.re_;

	cmpx32_t d_x 	= {(int32_t)x_0.re_ + (ywr >> 15), (int32_t)x_0.im_ + (ywi >> 15)};
	cmpx32_t d_y 	= {(int32_t)x_0.re_ - (ywr >> 15), (int32_t)x_0.im_ - (ywi >> 15)};
	*x_1 			= (cmpx16_t){(int16_t)(d_x.re_ >> 1) , (int16_t)(d_x.im_ >> 1)};
	*y_1 			= (cmpx16_t){(int16_t)(d_y.re_ >> 1) , (int16_t)(d_y.im_ >> 1)};
}

void coef_init(cmpx16_t *w)
{
	for(uint16_t idx = 0; idx < NFFT; idx ++)
	{
		w[idx].re_ = (int16_t)round(AMPL * cos(- 2 * M_PI * idx / NFFT));
		w[idx].im_ = (int16_t)round(AMPL * sin(- 2 * M_PI * idx / NFFT));

	}
}
