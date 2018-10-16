#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#define PI 3.141592653589793
#define SQRTPI 1.77245385091
#define GCONST 1.8867186677527936 //= (2*PI*)/4/sqrt(ln(2))
#define GCONST_SQUARED 3.5597073312468765 //= ( (2*PI*)/4/sqrt(ln(2)) )**2
#define SQRT2 1.4142135623730951

//Compile with: gcc -fPIC -shared -o lineshapes_cfuncs.so -O3 -ffast-math lineshapes_c.c

int sl_t(double amp, double wid, double wid0, double center, double* p2, double* weights, double* t, int npoints, int ntheta, double* fid_re, double * fid_im){
/* This function calculates a Super-Lorentzian FID made up of Gaussians*/    
    int i,j;
    double W[ntheta], W_squared[ntheta];
    double a,b,c;
    double a_squared;

    //Calculate W(theta)
    for(i=0;i<ntheta;i++){
        W_squared[i] = wid0*wid0 + wid*wid*p2[i]*p2[i];
        //W[i] = sqrt(W_squared[i]);
    }

    for(i=0;i<npoints;i++){
        fid_re[i]=0.0;
        for(j=0;j<ntheta;j++){
            //calculate Gaussian
            fid_re[i]+=exp(-W_squared[j]*GCONST_SQUARED*t[i]*t[i])*weights[j]; //*2*W[j]*GCONST/SQRTPI
        }

        fid_im[i]=amp*fid_re[i]*sin(-2*PI*center*t[i]);
        fid_re[i]*=amp*cos(-2*PI*center*t[i]);
    }

    return 1;
}

int sl_t_lorentz(double amp, double wid, double wid0, double center, double* p2, double* weights, double* t, int npoints, int ntheta, double* fid_re, double * fid_im){
/* This function calculates a Super-Lorentzian FID made up of Lorentzians.
   This is unphysical but might be a useful approximation? */    
    int i,j;
    double W[ntheta];
    double a,b,c;
    double a_squared;

    //Calculate W(theta)
    for(i=0;i<ntheta;i++){
        W[i] = sqrt(wid0*wid0 + wid*wid*p2[i]*p2[i]);
        printf("%.3e\n",W[i]);
    }

    for(i=0;i<npoints;i++){
        fid_re[i]=0.0;
        for(j=0;j<ntheta;j++){
            //calculate Gaussian
            fid_re[i]+=exp(-W[j]*t[i])*weights[j]; //*2*W[j]*GCONST/SQRTPI
        }

        fid_im[i]=amp*fid_re[i]*sin(-2*PI*center*t[i]);
        fid_re[i]*=amp*cos(-2*PI*center*t[i]);
    }

    return 1;
}

int sl_f(double amp, double wid, double wid0, double center, double* p2, double* weights, double* freqs, int npoints, int ntheta, double* spect_re){
/* This function calculates a Super-Lorentzian lineshape (spectrum) made up of Gaussians*/        
    int i,j;
    double W[ntheta];
    double a,b,c;

    //Calculate W(theta)
    for(i=0;i<ntheta;i++){
        W[i] = sqrt(wid0*wid0 + wid*wid*p2[i]*p2[i]);
    }

    for(i=0;i<npoints;i++){
        spect_re[i]=0.0;
        for(j=0;j<ntheta;j++){
            a = W[j];
            b = weights[j]/SQRTPI/SQRT2/a;
            c = pow((freqs[i]-center)/a,2.0);
            spect_re[i] += b*exp(-0.5*c);
        }
    }

    return 1;
}


/// TESTING CODE BELOW
/*
int main(int argc,char *argv[]){
    int i,j;
    int n_t = 131000;
    int n_ang = 500;
    double dwell = 1e-6;
    double wid = 20e3;
    double wid0 = 1.0; //can't be 0
    double center = 0;
    double sw = 500e3;

    double t[n_t], freq[n_t], fid_re[n_t], fid_im[n_t], thetas[n_ang], weights[n_ang], spect_re[n_t];

    for(i=0;i<n_ang;i++){
        thetas[i] = ((double)i)/n_ang*PI/2;
        weights[i] = 1.0;
    }

    for(i=0;i<n_t;i++){
        t[i]=i*dwell;
        freq[i] = sw/n_t*(i - n_t/2.0);
    }


    sl_t(wid, wid0, center, thetas, weights, t, n_t, n_ang,fid_re,fid_im);
    sl_f(wid, wid0, center, thetas, weights, freq, n_t, n_ang,spect_re);


    FILE *fstream=NULL;
    fstream = fopen("sl_c.txt","w");
    for(i=0;i<n_t;i++){
//        fprintf(fstream,"%f %f\n",fid_re[i],fid_im[i]);
        fprintf(fstream,"%f %f\n",freq[i],spect_re[i]);
    }
    fclose(fstream);


    return 0;
}



*/
