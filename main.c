//
//  main.c
//  kappasampler
//
//  Created by Jordy on 13/07/2017.
//  Copyright © 2017 Jordy. All rights reserved.
//

#include <stdio.h>
#include "gsl/gsl_sf_hyperg.h"
#include "gsl/gsl_sf_gamma.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "math.h"
#include <stdlib.h>
#include <time.h>

gsl_rng *r;

#define KAPPA 1
#define PWL 0
#define TH 0

double kappa=4.0;
double p=3.0;

struct f_params {
  double u;
};

double find_y(double u, double (*df)( double),double (*f)( double,void *),double Thetae){
// Get maximum for window
  int status, steps = 0, max_steps = 10;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double solution;
  double low = 1e-2;//sqrt(0.001/Thetae);
  double high = 1e3;//sqrt(1e5/Thetae);
//  fprintf(stderr,"%e %e\n",f(low,u),f(high,u));
  gsl_function F;
  struct f_params params = {u};
  F.function = f;
  F.params = &params;
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc(T);
  gsl_root_fsolver_set(s, &F, low, high);
  //printf("trying to find %e\n",u);
  do {
    steps++;
    status = gsl_root_fsolver_iterate(s);
    solution = gsl_root_fsolver_root(s);
    low = gsl_root_fsolver_x_lower(s);
    high = gsl_root_fsolver_x_upper(s);
    status = gsl_root_test_interval(low, high, 0, 1e-4);
  } while (status=GSL_CONTINUE && steps < max_steps);
  //fprintf(stderr,"steps %d solution %e\n",steps,solution);
return solution;
}

double dF3(double y){
    double value,denom,num;
    double y2=y*y;
    
    num=4.*y2*sqrt(kappa)*pow((kappa+y2)/(kappa),-kappa)*gsl_sf_gamma(kappa);
    denom=sqrt(M_PI)*(y2+kappa)*gsl_sf_gamma(kappa-1./2.);
    value=num/denom;
    
    return value;
}
double dF4(double y){
    double value,denom,num;
    double y2 = y*y;
    
    num=2.*(kappa-1)*y2*y*pow((kappa+y2)/(kappa),-kappa-1.);
    denom=kappa;
    value=num/denom;
    
    return value;
}

double dF5(double y){
    double value,denom,num;
    double y2 = y*y;

    num=8*y2*y2*pow((kappa+y2)/(kappa),-kappa-1) *  gsl_sf_gamma(kappa);
    denom=3*sqrt(M_PI)*pow(kappa,3./2.)*  gsl_sf_gamma(kappa-3./2.);
    value=num/denom;
    
    return value;
}
double dF6(double y){
    double value,denom,num;
    double y2 = y*y;

    num=(kappa*kappa-3.*kappa+2)*pow(y,5.)*pow((kappa+y2)/(kappa),-kappa-1);
    denom=kappa*kappa;
    value=num/denom;
    
    return value;
}


double F3(double y,void *params){
    struct f_params *p = (struct f_params*) params;
    double u = p->u;
    double value,denom,num;
    double y2=y*y;
    double a = 1;
    double b = -kappa - 1./2.;
    double c = 1./2.;
    double z = -y2/kappa;
    double hyp2F1;

//   if(abs(z)<1e-8){
//	hyp2F1=1+a*b*z/c + a*(1+a)*b*(1+b)/(2*c*(1+c))*z*z;
  //   }
   // else{
	hyp2F1 = pow(1.-z, -a) * tgamma(c) * tgamma(b-a)
    / (tgamma(b)*tgamma(c-a))
    * gsl_sf_hyperg_2F1(a, c-b, a-b+1., 1./(1.-z))
    + pow(1.-z, -b) * tgamma(c) * tgamma(a-b)
    / (tgamma(a) * tgamma(c-b))
    * gsl_sf_hyperg_2F1(b, c-a, b-a+1., 1./(1.-z));
    //}

//if(abs(z)<1e-4)
 //       hyp2F1=1+a*b*z/c;//gsl_sf_hyperg_2F1(a,b,c,z);
//else 

   num= -sqrt(kappa)*pow(((y2 + kappa)/kappa),-kappa) * gsl_sf_gamma(kappa)*(-kappa*hyp2F1+y2 * ( 2*kappa+1)+kappa);
    denom=y*sqrt(M_PI)*gsl_sf_gamma(3./2.+kappa);

    value=num/denom-u;
    
    return value;
}
double F4(double y,void *params){
    struct f_params *p = (struct f_params*) params;
    double u = p->u;
    double value,denom,num;
    double y2=y*y;

    num=1-(1+y2)*pow((kappa+y2)/(kappa),-kappa);
    denom=1;
    
    value=num/denom-u;
    
    return value;
}
double F5(double y,void *params){
    struct f_params *p = (struct f_params*) params;
    double u = p->u;
    double value,denom,num;
    double y2=y*y;
    double a = 1.;
    double b = -kappa - 1./2.;
    double c = 1./2.;
    double z = -y2/kappa;
    double hyp2F1;
//   if(abs(z)<1e-8){
  //      hyp2F1=1+a*b*z/c + a*(1+a)*b*(1+b)/(2*c*(1+c))*z*z;
   //  }
   // else{
    hyp2F1= pow(1.-z, -a) * tgamma(c) * tgamma(b-a)
    / (tgamma(b)*tgamma(c-a))
    * gsl_sf_hyperg_2F1(a, c-b, a-b+1., 1./(1.-z))
    + pow(1.-z, -b) * tgamma(c) * tgamma(a-b)
    / (tgamma(a) * tgamma(c-b))
    * gsl_sf_hyperg_2F1(b, c-a, b-a+1., 1./(1.-z));
   // }
//fprintf(stderr,"z=%e hyp %e hyp2 %e\n",z,gsl_sf_hyperg_2F1(b, c-a, b-a+1., 1./(1.-z)),gsl_sf_hyperg_2F1(a, c-b, a-b+1., 1./(1.-z)));
//if(abs(z)<1e-5)
//        hyp2F1=1+a*b*z/c;//gsl_sf_hyperg_2F1(a,b,c,z);
//else 

    num=pow((y2+kappa)/kappa,-kappa)*gsl_sf_gamma(kappa+1)*(3*kappa*kappa *(hyp2F1-1)+(1.-4.*kappa*kappa)*y2*y2-3.*kappa*(2.*kappa+1.)*y2);
    denom=3.*pow(kappa,3./2.)*y*sqrt(M_PI)*gsl_sf_gamma(3./2.+kappa);
    value= num/denom-u;
    
    return value;
}
double F6(double y,void *params){
    struct f_params *p = (struct f_params*) params;
    double u = p->u;
    double value,denom,num;
    double y2=y*y;
    double y4 = y2*y2;
    num=(y4-(y4+2.*y2+2.)*kappa)*pow((kappa+y2)/(kappa),-kappa);
    denom=2*kappa;
    value=num/denom+1-u;
    
    return value;
}


double sample_kappa_distr(double Thetae)
{
    
    double S_3, pi_3, pi_4, pi_5, pi_6, y, x1, x2, prob;
    double num, den;
    pi_3 = sqrt(kappa)*sqrt(M_PI)*gsl_sf_gamma(-1./2.+kappa)/(4.*gsl_sf_gamma(kappa));
    pi_4 = kappa /( 2.*kappa - 2.) * sqrt(0.5*Thetae);
    pi_5 = 3 * pow(kappa,3./2.) * sqrt(M_PI) * gsl_sf_gamma(-3./2. + kappa)/(8*gsl_sf_gamma(kappa)) * Thetae;
    pi_6 = kappa * kappa / ( 2. - 3.*kappa + kappa*kappa) * Thetae * sqrt(0.5 * Thetae);
    
    S_3 = pi_3 + pi_4 + pi_5 + pi_6;
    
    pi_3 /= S_3;
    pi_4 /= S_3;
    pi_5 /= S_3;
    pi_6 /= S_3;
    
    do {
        x1 = (double)rand()/(double)RAND_MAX;
        double u = (double)rand()/(double)RAND_MAX;

	if (x1 < pi_3) {
            //printf("try f3\n");
	    y = find_y(u,dF3,F3,Thetae);
        } else if (x1 < pi_3 + pi_4) {
          //  printf("try f4\n");
            y = find_y(u,dF4,F4,Thetae);
        } else if (x1 < pi_3 + pi_4 + pi_5) {
        //    printf("try f5\n");
            y = find_y(u,dF5,F5,Thetae);
        } else {
          //  printf("try f6\n");
            y = find_y(u,dF6,F6,Thetae);
        }
        
        x2 = (double)rand()/(double)RAND_MAX;
        num = sqrt(1. + 0.5 * Thetae * y * y);
        den = (1. + y * sqrt(0.5 * Thetae));
        
        prob = num / den;
        
    } while (x2 >= prob);
    
    return (y);
}


double sample_th_distr(double Thetae) {

  double S_3, pi_3, pi_4, pi_5, pi_6, y, x1, x2, x, prob;
  double num, den;

  pi_3 = sqrt(M_PI) / 4.;
  pi_4 = sqrt(0.5 * Thetae) / 2.;
  pi_5 = 3. * sqrt(M_PI) * Thetae / 8.;
  pi_6 = Thetae * sqrt(0.5 * Thetae);
  S_3 = pi_3 + pi_4 + pi_5 + pi_6;

  pi_3 /= S_3;
  pi_4 /= S_3;
  pi_5 /= S_3;
  pi_6 /= S_3;

  do {
    x1 = (double)rand()/(double)RAND_MAX;

    if (x1 < pi_3) {
      x = gsl_ran_chisq(r, 3);
    } else if (x1 < pi_3 + pi_4) {
      x = gsl_ran_chisq(r, 4);
    } else if (x1 < pi_3 + pi_4 + pi_5) {
      x = gsl_ran_chisq(r, 5);
    } else {
      x = gsl_ran_chisq(r, 6);
    }

    /* this translates between defn of distr in
     Canfield et al. and standard chisq distr */
    y = sqrt(x / 2);

    x2 = (double)rand()/(double)RAND_MAX;
    num = sqrt(1. + 0.5 * Thetae * y * y);
    den = (1. + y * sqrt(0.5 * Thetae));
    //                    printf("th %e
    //                    %e\n",(1+Thetae*y*y)/gamma_max,(1+Thetae*y*y) );

    prob = num / den;

  } while (x2 >= prob);

  return (y);
}

double kappa_fit(double gamma, double Thetae){

	return gamma*sqrt(gamma*gamma-1.)*pow(1+(gamma-1)/(kappa*Thetae),-kappa-1);

}

double pwl_fit(double gamma,double gmin,double gmax){
return (p-1)*pow(gamma,-p)/(4*M_PI*(pow(gmin,1-p)-pow(gmax,1-p)));

}

double th_fit(double gamma, double Thetae){
	return gamma*gamma * sqrt(1-1/(gamma*gamma)) * exp(-gamma/Thetae);

}

double sample_pwl_distr(double gmin,double gmax){

double x1 = (double)rand()/(double)RAND_MAX;

return  pow(pow(gmin,1-p)*(1 - x1 ) + x1 * pow(gmax,1-p),1/(1-p));

}


int main(int argc, const char * argv[]) {
    // insert code here...
    time_t currtime, starttime;
    starttime = time(NULL);
    FILE *fp;
    
    srand(1337);
    
    r = gsl_rng_alloc(gsl_rng_mt19937); /* use Mersenne twister */
    gsl_rng_set(r, 1337);

    double Thetae;
    int distr_size;
    double sample_size_d;
    long int sample_size;
    double gmin=2.5;
    double gmax;
    sscanf(argv[1], "%lf", &Thetae);
    sscanf(argv[2], "%d", &distr_size);
    sscanf(argv[3], "%lf", &sample_size_d);
    sscanf(argv[4], "%lf", &gmax);
    sscanf(argv[5], "%le", &p);

    sample_size =(long int) sample_size_d;

    fprintf(stderr,"%lf %d %ld %lf\n",Thetae,distr_size,sample_size,gmax);
    
    double dg=((log(gmax-1.) - log(gmin)) /distr_size);
    
    fprintf(stderr,"%e\n",dg);
    
    double y;
    double gamma_e;
    double sample=0;
    double *distr;
    double *analytical_distr;
    analytical_distr = (double *)malloc(distr_size * sizeof(double));
    distr = (double *)malloc(distr_size * sizeof(double));

#pragma omp parallel for
    for(int i=0;i<sample_size;i++){
	    if(fmod(i,0.1*sample_size)==0)
        	printf("%d\n",i);

#if KAPPA
	    y=sample_kappa_distr(Thetae);
	    gamma_e = y * y * Thetae + 1.;
#elif PWL
        gamma_e=sample_pwl_distr(gmin+1,gmax+1);
#elif TH
	    y=sample_th_distr(Thetae);
	    gamma_e = y * y * Thetae + 1.;
#endif

    
	    if((int)((log(gamma_e-1)-log(gmin))/dg)<distr_size && (int)((log(gamma_e-1)-log(gmin))/dg)>=0){
#pragma omp atomic
		    distr[(int)((log(gamma_e-1)-log(gmin))/dg)]+=1;
#pragma omp atomic
		    sample++;

        }
    
    }

    currtime = time(NULL);
    fprintf(stderr, "time %g, rate %g ph/s\n",
           (double) (currtime - starttime),
           sample / (currtime - starttime));

    char spec_filename[50];
    sprintf(spec_filename, "distr.txt");
    fp = fopen(spec_filename,"w");

    sample =0;
    double norm=0;
    
    for(int i=0;i<distr_size;i++){
	    sample+=distr[i];
        gamma_e	= exp((i+0.5)*dg + log(gmin))+1.;

#if KAPPA
	    norm += kappa_fit(gamma_e,Thetae) *(gamma_e-1.)*dg;
#elif PWL
            norm += pwl_fit(gamma_e,gmin,gmax) *(gamma_e-1)*dg;
#elif TH
	    norm += th_fit(gamma_e,Thetae) *(gamma_e-1.)*dg;
#endif

	}
    for(int i=0;i<distr_size;i++){
        gamma_e	= exp((i+0.5)*dg + log(gmin))+1.;
        distr[i]/=(double)(sample*(gamma_e-1.)*dg);
#if KAPPA
	    analytical_distr[i]=kappa_fit(gamma_e,Thetae)/norm;
#elif PWL
            analytical_distr[i]=pwl_fit(gamma_e,gmin,gmax)/norm;
#elif TH
	    analytical_distr[i]=th_fit(gamma_e,Thetae)/norm;
#endif
        fprintf(fp,"%e %e %e\n",gamma_e-1, distr[i], analytical_distr[i]);
    }
    
    currtime = time(NULL);
    fprintf(stderr, "time %g, rate %g ph/s\n",
           (double) (currtime - starttime),
           sample / (currtime - starttime));

    return 0;
}
