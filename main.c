//
//  main.c
//  Samplers
//
//  Created by Jordy Davelaar on 16/03/2023.
//  Copyright © 2023 Jordy Davelaar. All rights reserved.
//

#include <stdio.h>
#include "gsl/gsl_sf_hyperg.h"
#include "gsl/gsl_sf_gamma.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include "math.h"
#include <stdlib.h>

#define KAPPA 0
#define PWL 1

double kappa=4.0;
double p=3.0;

//#define  M_PI 3.14159265358979323846

struct f_params {
  double u;
};

double find_y(double u, double (*df)( double),double (*f)( double,void *),double Thetae){
// Get maximum for window
  int status, steps = 0, max_steps = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double solution;
  double low = 1e-5;
  double high =1e5; 

  gsl_function F;
  struct f_params params = {u};
  F.function = f;
  F.params = &params;
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc(T);
  gsl_root_fsolver_set(s, &F, low, high);

  do {
    steps++;
    status = gsl_root_fsolver_iterate(s);
    solution = gsl_root_fsolver_root(s);
    low = gsl_root_fsolver_x_lower(s);
    high = gsl_root_fsolver_x_upper(s);
    status = gsl_root_test_interval(low, high, 1e-10, 0);
  } while (status=GSL_CONTINUE && steps < max_steps);

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

    num=8*y2*y2*pow((kappa+y2)/(kappa),-kappa-1)*gsl_sf_gamma(kappa);
    denom=3*sqrt(M_PI)*pow(kappa,3./2.)*gsl_sf_gamma(kappa-3./2.);
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

    hyp2F1 = pow(1.-z, -a) * tgamma(c) * tgamma(b-a)
    / (tgamma(b)*tgamma(c-a))
    * gsl_sf_hyperg_2F1(a, c-b, a-b+1., 1./(1.-z))
    + pow(1.-z, -b) * tgamma(c) * tgamma(a-b)
    / (tgamma(a) * tgamma(c-b))
    * gsl_sf_hyperg_2F1(b, c-a, b-a+1., 1./(1.-z));

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

    hyp2F1= pow(1.-z, -a) * tgamma(c) * tgamma(b-a)
    / (tgamma(b)*tgamma(c-a))
    * gsl_sf_hyperg_2F1(a, c-b, a-b+1., 1./(1.-z))
    + pow(1.-z, -b) * tgamma(c) * tgamma(a-b)
    / (tgamma(a) * tgamma(c-b))
    * gsl_sf_hyperg_2F1(b, c-a, b-a+1., 1./(1.-z));

    num=pow((y2+kappa)/kappa,-kappa)*gsl_sf_gamma(kappa+1)*(3*kappa*kappa *(hyp2F1-1)+(1.-4.*kappa*kappa)*y2*y2-3.*kappa*(2.*kappa+1.)*y2);
    denom=3.*pow(kappa,3./2.)*y*sqrt(M_PI)*gsl_sf_gamma(3./2.+kappa);
    value=num/denom-u;
    
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


double sample_y_distr(double Thetae)
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
	    y = find_y(u,dF3,F3,Thetae);
        } else if (x1 < pi_3 + pi_4) {
            y = find_y(u,dF4,F4,Thetae);
        } else if (x1 < pi_3 + pi_4 + pi_5) {
            y = find_y(u,dF5,F5,Thetae);
        } else {
            y = find_y(u,dF6,F6,Thetae);
        }
        
        x2 = (double)rand()/(double)RAND_MAX;
        num = sqrt(1. + 0.5 * Thetae * y * y);
        den = (1. + y * sqrt(0.5 * Thetae));
        
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

double sample_gam_distr(double gmin,double gmax){

       double x1 = (double)rand()/(double)RAND_MAX;

       return  pow(pow(gmin,1-p)*(1 - x1 ) + x1 * pow(gmax,1-p),1/(1-p));

}


int main(int argc, const char * argv[]) {
    FILE *fp;
    
    srand(1337);
    
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

    double *distr;
    double *analytical_distr;
    analytical_distr = (double *)malloc(distr_size * sizeof(double));
    distr = (double *)malloc(distr_size * sizeof(double));


#if KAPPA

    for(int i=0;i<distr_size;i++){
        distr[i]=0;
        gamma_e = exp((i+0.5)*dg + log(gmin))+1.;
	analytical_distr[i]=kappa_fit(gamma_e,Thetae);
    }

#pragma omp parallel for
    for(int i=0;i<sample_size;i++){
	if(fmod(i,0.1*sample_size)==0)
        	printf("%d\n",i);
	y=sample_y_distr(Thetae);

	gamma_e = y * y * Thetae + 1.;
    {
	if((int)((log(gamma_e-1)-log(gmin))/dg)<distr_size && (int)((log(gamma_e-1)-log(gmin))/dg)>=0){
	#pragma omp atomic
		distr[(int)((log(gamma_e-1)-log(gmin))/dg)]+=1;
        }
    }
    }

    char spec_filename[50];
    sprintf(spec_filename, "distr_t_%.02lf_thermal.txt",Thetae);
    fp = fopen(spec_filename,"w");
    double sample =0;
    double norm=0;
    for(int i=0;i<distr_size;i++){
	sample+=distr[i];
        gamma_e	= exp((i+0.5)*dg + log(gmin))+1.;
	norm += kappa_fit(gamma_e,Thetae) *(gamma_e-1.)*dg;
	}

    for(int i=0;i<distr_size;i++){
        gamma_e	= exp((i+0.5)*dg + log(gmin))+1.;
        distr[i]/=(double)(sample*(gamma_e-1.)*dg);
	analytical_distr[i]=kappa_fit(gamma_e,Thetae)/norm;
        fprintf(fp,"%e %e %e\n",gamma_e-1, distr[i], analytical_distr[i]);
    }
    
    
#elif PWL

fprintf(stderr,"gmin %e gmax %e pwl %e\n",gmin+1,gmax+1,p);

for(int i=0;i<distr_size;i++){
        distr[i]=0;
        //if(i>0){
       	gamma_e = exp((i+0.5)*dg + log(gmin))+1.;
        analytical_distr[i]=pwl_fit(gamma_e,gmin,gmax);
        }

#pragma omp parallel for
    for(int i=0;i<sample_size;i++){
        if(fmod(i,0.1*sample_size)==0)
                printf("%d\n",i);
        gamma_e=sample_gam_distr(gmin+1,gmax+1);
     	if((int)((log(gamma_e-1)-log(gmin))/dg)<distr_size && (int)((log(gamma_e-1)-log(gmin))/dg)>=0){
       	#pragma omp atomic
                distr[(int)((log(gamma_e-1)-log(gmin))/dg)]+=1;
	}
    }

 char spec_filename[50];
    sprintf(spec_filename, "distr_pwl_%.02lf.txt",p);
    fp = fopen(spec_filename,"w");
    double sample =0;
    double norm=0;
    for(int i=0;i<distr_size;i++){
        sample+=distr[i];
        gamma_e = exp((i+0.5)*dg + log(gmin))+1.;
        norm += analytical_distr[i]*(gamma_e-1)*dg;
        }
    printf("norm %e sample %e\n",norm,sample);
    for(int i=0;i<distr_size;i++){
        gamma_e = exp((i+0.5)*dg + log(gmin))+1.;
        distr[i]/=(double)(sample*(gamma_e-1)*dg);
        analytical_distr[i]/=(double)norm;
        fprintf(fp,"%e %e %e\n",gamma_e-1, distr[i], analytical_distr[i]);
    }
    fclose(fp);
#endif
    return 0;
}
