#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
int np = 30, ne = 512;
//int np = 5,  ne = 10; 
                      // np: number of nodes in each element
					  // ne: number of elements
					  // caution: np 20, ne 256 will take approx 10 minutes on a serial job

void interface_flux(double *qq, double *fstar, double *ib, double speed);
void rhs(double *qq, double *rr, double *dv, double *df, double *ib, double speed);
void lagrange(double xx, double *pts, double *ll);
void lagrange_deriv(double xx, double *dl, double *pts);
void legendre(int n, double *x, int xsize, double *output);
double legendre_scalar(int n, double x);
double dlegendre_scalar(int n, double x);
double ddlegendre_scalar(int n, double x);
double dddlegendre_scalar(int n, double x);
void dlegendre_roots(int polyorder, double *output);
void gausslobatto_quadrature(int polyorder, double *roots, double *weights);
double dot_product(double *v, double *u, int n);
void save_field(double *xx, double *qq, int elem_num, double *roots, int eres);
void initialize(double *qq, double *xx, double xmax, double xmin, char *init_type);

int main(int argc, char **argv){

	double tend = 1E2, speed = 1.;
	// double tend = 1E-1, speed = 1.;
	char *init_type="mixed2";
	double *roots, *weights, *ll, *dl, xmin, xmax, 
		   deltax, jac, xr, xl, cfl, dt, rtime, min_dx;
	int ii, jj, kk, ee, idx, eres;
	long nstep;
	double *dx, *mesh; 
	double *smat, *xx, *qq, *qtemp, *k1, *k2, *k3, *k4, *minv_vec, *mmat, *dv, 
		   *mf, *ib, *df, *fstar;

	// initialize 
	// fortran index structure array[ii,jj,ee] where size(array) = (np, np, ne)
	// c 1d index structure array = [ee*np*np + jj*np + ii]
	roots   = (double*)malloc(np*   sizeof(double));
	weights = (double*)malloc(np*   sizeof(double));
	ll      = (double*)malloc(np*   sizeof(double));
	dl      = (double*)malloc(np*   sizeof(double));
	dx      = (double*)malloc(ne*   sizeof(double));
	mesh    = (double*)malloc((ne+1)*sizeof(double));

	smat	= (double*)malloc(np*np*sizeof(double));    // [jj np, ii np]
	xx		= (double*)malloc(ne*np*sizeof(double));    // [ee ne, ii np]
	qq		= (double*)malloc(ne*np*sizeof(double));    // [ee ne, ii np]	
	qtemp	= (double*)malloc(ne*np*sizeof(double));    // [ee ne, ii np]	
	k1		= (double*)malloc(ne*np*sizeof(double));    // [ee ne, ii np]	
	k2		= (double*)malloc(ne*np*sizeof(double));    // [ee ne, ii np]
	k3		= (double*)malloc(ne*np*sizeof(double));    // [ee ne, ii np]
	k4		= (double*)malloc(ne*np*sizeof(double));    // [ee ne, ii np]
	minv_vec= (double*)malloc(ne*np*sizeof(double));    // [ee ne, ii np]
	mmat	= (double*)malloc(ne*np*np*sizeof(double)); // [ee ne, jj np, ii np]
	dv		= (double*)malloc(ne*np*np*sizeof(double)); // [ee ne, jj np, ii np]
	mf		= (double*)malloc(2*np*sizeof(double));     // [jj 2,  ii np]
	ib		= (double*)malloc(2*np*sizeof(double));     // [jj 2,  ii np]
	fstar	= (double*)malloc(2*ne*sizeof(double));     // [jj 2,  ii ne]
	df		= (double*)malloc(ne*2*np*sizeof(double));  // [ee ne, jj 2, ii np]

	for (ii=0; ii<np; ++ii){
		roots[ii] = 0;
		weights[ii] = 0;
		ll[ii] = 0;
		dl[ii] = 0;
	}
	for (ii=0; ii<ne; ++ii){
		dx[ii] = 0;
		mesh[ii] = 0;
	}
	mesh[ne] = 0;

	
	for (ii=0; ii<np*np; ++ii){
		smat[ii] = 0;
	}
	for (ii=0; ii<ne*np; ++ii){
		xx[ii]	= 0;		 	
		qq[ii]	= 0;		 	
		k1[ii]	= 0;		 	
		k2[ii]	= 0;		 	
		k3[ii]	= 0;		 	
		k4[ii]	= 0;		 	
		minv_vec[ii]	= 0; 
	}
	for (ii=0; ii<ne*np*np; ++ii){
		mmat[ii] = 0;    	
		dv[ii]	 = 0;
	}
	for (ii=0; ii<np*2; ++ii){
		mf[ii] = 0;
		ib[ii] = 0;
	}
	for (ii=0; ii<ne*2; ++ii){
		fstar[ii] = 0;
	}
	for (ii=0; ii<ne*2*np; ++ii){
		df[ii] = 0;
	}

	// mesh setup
	xmin = 0.;
	xmax = 10.;
	deltax = (xmax-xmin)/(double)ne;
	mesh[ne] = xmax;
	for(ee=0;ee<ne;++ee){
		mesh[ee] = xmin+ee*deltax;
	}
	
	// gauss lobatto quadrature point, weight setup
	gausslobatto_quadrature(np, roots, weights);

	// coordinates and element size
	min_dx = xmax - xmin; // initial guess
	for(ee=0;ee<ne;ee++){
		xl = mesh[ee];
		xr = mesh[ee+1];
		dx[ee] = xr-xl; // size of each element
		if(dx[ee] < min_dx){
			min_dx = dx[ee]; // finding minimum dx
		}
		for(ii=0;ii<np;ii++){
			idx = ee*np+ii;
			xx[idx] = xl + 0.5*(1+roots[ii])*dx[ee];
		}
	}

	// mass matrix
	for(ii=0;ii<ne*np*np;ii++){
		mmat[ii] = 0;
	}
	for(ee=0;ee<ne;ee++){
		jac = fabs(dx[ee])/2;
		for(kk=0;kk<np;kk++){
			lagrange(roots[kk], ll, roots);
			for(jj=0;jj<np;jj++){
				for(ii=0;ii<np;ii++){
					idx = ee*np*np+jj*np+ii;
					// mass matrix mmat[ne][np][np] in 1d index representation
					mmat[idx] += jac*weights[kk]*ll[ii]*ll[jj];
				}
			}
		}
	}

	// stiffness matrix
	for(ii=0;ii<np*np;ii++){
		smat[ii] = 0;
	}
	for(kk=0;kk<np;kk++){
		lagrange(roots[kk], ll, roots);
		lagrange_deriv(roots[kk], dl, roots);
		for(jj=0;jj<np;jj++){
			for(ii=0;ii<np;ii++){
				idx = jj*np+ii;
				// stiffness matrix smat[np][np] in 1d index representation
				smat[idx] += weights[kk]*ll[jj]*dl[ii];
			}
		}
	}

	// face integration
	for(ii=0;ii<np*2;ii++){
		mf[ii] = 0;
	}
	lagrange(-1,mf,   roots); // mf[ii] for(ii=0, ii<np,ii++) represents element left face integration
	lagrange( 1,mf+np,roots); // mf[ii] for ii=np, ii<2*np, ii++) reresents element right face integration

	// boundary interpolation
	for(ii=0;ii<np*2;ii++){
		ib[ii] = 0;
	}
	lagrange(-1,ib,   roots); // element left edge interpolation
	lagrange( 1,ib+np,roots); // element right edge interpolation

	
	// divergence operators
	for(ii=0;ii<ne*np*np;ii++){
		dv[ii] = 0;
	}
	for(ii=0;ii<ne*np*2;ii++){
		dv[ii] = 0;
	}
	for(ee=0;ee<ne;ee++){
		for(jj=0;jj<np;jj++){
			// it turn out that mmat is diagonal. i.e., ii != jj, mmat[ee][jj][ii] = 0
			// the inverse of mmat is just the inverse of the diagonal components
			// here, we are extracting the inverse diagonal components only
			minv_vec[ee*np+jj] = 1./mmat[ee*np*np+jj*np+jj];
		}
		for(jj=0;jj<np;jj++){
			for(ii=0;ii<np;ii++){
				dv[ee*np*np+jj*np+ii] = minv_vec[ee*np+ii]*smat[jj*np+ii];
			}
		}
		for(jj=0;jj<2;jj++){
			for(ii=0;ii<np;ii++){
				df[ee*np*2+jj*np+ii]  = minv_vec[ee*np+ii]*mf[jj*np+ii];
			}
		}

	}
	
	// initialize qq field
	initialize(qq, xx, xmax, xmin, init_type);
	cfl = 1./(np*np);
	dt = cfl * min_dx / fabs(speed);
	rtime = 0.;
	nstep = 0;

	printf("Start Time Integration\n");

	// Runge-Kutta 4th order Time integration loop
	
	while(rtime < tend){
		dt = fmin(dt, tend-rtime);

		rhs(qq,	   k1, dv, df, ib, speed);

		for(ii=0;ii<ne*np;ii++)
			qtemp[ii] = qq[ii]+0.5*dt*k1[ii];
		rhs(qtemp, k2, dv, df, ib, speed);

		for(ii=0;ii<ne*np;ii++)
			qtemp[ii] = qq[ii]+0.5*dt*k2[ii];
		rhs(qtemp, k3, dv, df, ib, speed);
		
		for(ii=0;ii<ne*np;ii++)
			qtemp[ii] = qq[ii]+dt*k3[ii];
		rhs(qtemp, k4, dv, df, ib, speed);

		for(ii=0;ii<ne*np;ii++)
			qq[ii] += 1./6.*dt*(k1[ii]+2*k2[ii]+2*k3[ii]+k4[ii]);

		rtime += dt;
		nstep += 1;
		if(nstep%10000 == 0) 
			printf("nstep = %10ld, %5.1f%% complete\n", nstep, rtime/tend*100);
	}

	// timeloop ends here;

	printf("Integration complete\n");

	if(ne > 200){
		eres = 2;
	}
	else if (ne > 60){
		eres = 3;
	}
	else if (ne > 30){
		eres = 6;
	}
	else {
		eres = 10;
	}

	// final report
	printf("-----------------------------------------------\n");
	printf("code type   : c serial\n");
	printf("Final time  : %13.5e\n", rtime);
	printf("CFL         : %13.5e\n", cfl);
	printf("DOF         : %13d\n", ne*np);
	printf("No. of Elem : %13d\n", ne);
	printf("Order       : %13d\n", np);
	printf("eres        : %13d\n", eres);
	printf("time steps  : %13ld\n", nstep);
	printf("-----------------------------------------------\n");

	save_field(xx, qq, ne, roots, eres);

	free(roots);   
	free(weights); 
	free(ll);      
	free(dl);      
	free(dx);      
	free(mesh);    
	free(smat);	
	free(xx);		
	free(qq);		
	free(qtemp);	
	free(k1);		
	free(k2);		
	free(k3);		
	free(k4);		
	free(minv_vec);
	free(mmat);	
	free(dv);		
	free(mf);		
	free(ib);		
	free(fstar);	
	free(df);		

	return 0;
}

////////////////////////////////////////////////////////////////////////////
// parallelization needed functions: interface_flux & rhs

void interface_flux(double *qq, double *fstar, double *ib, double speed){
	int ii, ee;
	double qb[2*ne];

	for(ii=0;ii<2*ne;ii++){
		fstar[ii] = 0;
	}
	// fstar[0:ne-1] stores numerical flux of the left edge on each elements
	// fstar[ne:2*ne-1] stores numerical flux of the right edge on each elements

	for(ii=0;ii<ne;ii++){
		qb[ii]   = dot_product(ib,    qq+ii*np, np); // left edge interpolated value of qq at element ii
		qb[ne+ii] = dot_product(ib+np, qq+ii*np, np); // right edge interpolated value of qq at element ii
	}
	
	ii = 0; // calculating numerical flux (fstar) with periodic boundary condition
	fstar[ii] = -((qb[ne+ne-1]+qb[ii])/2*speed + fabs(speed)*(qb[ne+ne-1]-qb[ii])/2);
	fstar[ne+ne-1] = -fstar[0];
	for(ii=1;ii<ne;ii++){
		fstar[ii] = -((qb[ne+ii-1]+qb[ii])/2*speed + fabs(speed)*(qb[ne+ii-1]-qb[ii])/2);
		fstar[ne+ii-1] = -fstar[ii];
	}
}

void rhs(double *qq, double *rr, double *dv, double *df, double *ib, double speed){
	int ii, jj, ee, idx;
	double fstar[2*ne], term;
	for(ii=0;ii<np*ne;ii++){
		rr[ii] = 0;
	}

	for(ee=0;ee<ne;ee++){
		for(ii=0;ii<np;ii++){
			for(jj=0; jj<np; jj++){
				rr[ee*np+ii] += dv[ee*np*np+jj*np+ii]*speed*qq[ee*np+jj];
			}
		}
	}


	interface_flux(qq, fstar, ib, speed);
	for(ee=0;ee<ne;ee++){
		for(ii=0; ii<np; ii++){
			for(jj=0;jj<2;jj++){
				rr[ee*np+ii] -= df[ee*(np*2)+jj*np+ii]*fstar[ee+jj*ne];
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////
// no modifications needed in the following functions for parallelization

void lagrange(double xx, double *ll, double *pts){
	int ii, jj;
	
	for(ii=0; ii<np; ii++){
		ll[ii] = 1.;
		for(jj=0; jj<np; jj++){
			if(jj != ii){
				ll[ii] *= (xx-pts[jj])/(pts[ii]-pts[jj]);
			}
		}
	}
}

void lagrange_deriv(double xx, double *dl, double *pts){
	int ii, jj, kk;
	double term;

	for(ii=0;ii<np;++ii){
		dl[ii] = 0.;
		for(jj=0;jj<np;++jj){
			if (jj != ii){
				term = 1./(pts[ii]-pts[jj]);
				for(kk=0;kk<np;++kk){
					if ((kk != ii) && (kk != jj)){
						term *= (xx-pts[kk])/(pts[ii]-pts[kk]);
					}
				}
				dl[ii] += term;
			}
		}
	}
}

void legendre(int n, double *x, int xsize, double *output){
	double temp1[xsize], temp2[xsize];
	int ii;
	if(n == 0){
		for(ii=0;ii<xsize;++ii){
			output[ii] = 1.;
		}
	}
	else if (n == 1){
		for(ii=0;ii<xsize;++ii){
			output[ii] = x[ii];
		}
	}
	else{
		legendre(n-1, x, xsize, temp1);
		legendre(n-2, x, xsize, temp2);
		for(ii=0;ii<xsize;++ii){
			output[ii] = ((2*n-1)*x[ii]*temp1[ii]-(n-1)*temp2[ii])/n;
		}
	}
}

double legendre_scalar(int n, double x){
	double output, temp1, temp2;

	if (n == 0){
		output = 1.;
	}
	else if (n == 1){
		output = x;
	}
	else{
		temp1 = legendre_scalar(n-1, x);
		temp2 = legendre_scalar(n-2, x);
		output = ((2*n-1)*x*temp1-(n-1)*temp2)/n;
	}
	return output;
}

double dlegendre_scalar(int n, double x){
	double output, temp1, temp2;
	temp1 = legendre_scalar(n-1, x);
	temp2 = legendre_scalar(n,   x);
	output = n*(temp1-x*temp2)/(1.-(x*x));
	return output;
}

double ddlegendre_scalar(int n, double x){
	double output, temp1, temp2;
	temp1 = dlegendre_scalar(n, x);
	temp2 = legendre_scalar(n,  x);
	output = (2.*x*temp1-n*(n+1)*temp2)/(1.-(x*x));
	return output;
}

double dddlegendre_scalar(int n, double x){
	double output, temp1, temp2;
	temp1 = ddlegendre_scalar(n, x);
	temp2 = dlegendre_scalar(n,  x);
	output = (2.*x*temp1-n*(n+1)*temp2)/(1.-(x*x));
	return output;
}

void dlegendre_roots(int polyorder, double *output){
	double tolerance=1E-16;
	double x, error, dx, y, dy, ddy;
	double temproot[polyorder/2];
	int ii, iters, n;

	n = polyorder;
	for(ii=0;ii<n;++ii){
		output[ii] = 1.;
	}

	if (polyorder < 2){
		printf("Error: Polyorder less than 2 is not supported - dlegendre_roots\n");
		exit(1);
	}
	else{
		for(ii=1; ii<polyorder/2; ++ii){
			x = (1.-(3.*(n-2.)/(8.*(n-1)*(n-1)*(n-1))))*cos((4.*(ii+1.)-3.)/(4.*(n-1.)+1.)*M_PI);
			error = 10*tolerance;
			iters = 0;
			while ((error > tolerance) && (iters < 1000)){
				y   =   dlegendre_scalar(polyorder-1, x);
				dy  =  ddlegendre_scalar(polyorder-1, x);
				ddy = dddlegendre_scalar(polyorder-1, x);
				dx = -2.*y*dy/(2.*(dy*dy)-y*ddy);
				x += dx;
				iters += 1;
				error = fabs(dx);
			}
			output[ii] = x;
		}
		for(ii=0; ii<polyorder/2; ++ii){
			temproot[ii] = output[ii];
		}
		if(polyorder%2 == 0){
			for(ii=0; ii<polyorder/2; ++ii){
				output[ii] = -temproot[ii];
				output[polyorder/2+ii] = temproot[polyorder/2-1-ii];
			}
		}
		else{
			output[polyorder/2] = 0.;
			for(ii=0; ii<polyorder/2; ++ii){
				output[ii] = -temproot[ii];
				output[polyorder/2+1+ii] = temproot[polyorder/2-1-ii];
			}
		}
	}
}

void gausslobatto_quadrature(int polyorder, double *roots, double *weights){
	int n, ii;
	double temp[polyorder];

	if (polyorder == 1){
		roots[0] = 0.;
		weights[0] = 2.;
	}
	else{
		n = polyorder;
		dlegendre_roots(n, roots);
		legendre(n-1, roots, n, temp);
		for(ii=0;ii<n; ++ii){
			weights[ii] = 2./(n*(n-1)*(temp[ii]*temp[ii]));
		}
	}
}

double dot_product(double *v, double *u, int n){
	int ii;
	double result = 0.0;

	for (ii = 0; ii < n; ii++)
		result += v[ii]*u[ii];
	return result;
}

void save_field(double *xx, double *qq, int elem_num, double *roots, int eres){
	double dx, plot_coords[eres], ll[np], xcoord, ycoord;
	int ii, ee;

	dx = 2./(eres-1);
	plot_coords[0] = -1;
	plot_coords[eres-1] = 1;
	for(ii=1;ii<eres-1;++ii){
		plot_coords[ii] = -1+ii*dx;
	}

	FILE *ff;
	ff = fopen("data.txt","w");

	for(ee=0;ee<elem_num;++ee){
		for(ii=0;ii<eres;++ii){
			lagrange(plot_coords[ii], ll, roots);
			xcoord = dot_product(ll, xx+ee*np, np);
			ycoord = dot_product(ll, qq+ee*np, np);
			fprintf(ff, "%21.6f   %21.6f\n", xcoord, ycoord);
		}
	}

	fclose(ff);
	printf("output file saved!\n");
}

void initialize(double *qq, double *xx, double xmax, double xmin, char *init_type){
	int ii;
	double mu, sig;

	if(!strncmp(init_type, "sin", 3)){
		for(ii=0;ii<ne*np;ii++){
			qq[ii] = sin(2*M_PI*xx[ii]/(xmax-xmin));
		}
	}
	else if(!strncmp(init_type, "gaussian", 3)){
		mu = 0.5*(xmin+xmax);
		sig = 0.1*(xmax-xmin);
		for(ii=0;ii<ne*np;ii++){
			qq[ii] = exp(-((xx[ii]-mu)*(xx[ii]-mu))/(2*(sig*sig)));
		}
	}
	else if(!strncmp(init_type, "box", 3)){
		mu = 0.5*(xmin+xmax);
		sig = 0.1*(xmax-xmin);
		for(ii=0;ii<ne*np;ii++){
			qq[ii] = 0;
			if(fabs(mu-xx[ii])<sig){
				qq[ii] = 1;
			}
		}
	}
	else if(!strncmp(init_type, "mixed1", 6)){
		for(ii=0;ii<ne*np;ii++){
			qq[ii] = 0.2*sin(2*M_PI*xx[ii]/(xmax-xmin));
		}
		mu = 0.4*(xmin+xmax);
		sig = 5E-2*(xmax-xmin);
		for(ii=0;ii<ne*np;ii++){
			qq[ii] += 2*exp(-((xx[ii]-mu)*(xx[ii]-mu))/(2*(sig*sig)));
		}
	}
	else if(!strncmp(init_type, "mixed2", 6)){
		mu = 0.3*(xmin+xmax);
		sig = 0.1*(xmax-xmin);
		for(ii=0;ii<ne*np;ii++){
			qq[ii] = 0;
			if(fabs(mu-xx[ii])<sig){
				qq[ii] = 1;
			}
		}
		for(ii=0;ii<ne*np;ii++){
			qq[ii] += 0.2*sin(2*M_PI*xx[ii]/(xmax-xmin));
		}
		mu = 0.6*(xmin+xmax);
		sig = 0.1*(xmax-xmin);
		for(ii=0;ii<ne*np;ii++){
			qq[ii] += exp(-((xx[ii]-mu)*(xx[ii]-mu))/(2*(sig*sig)));
		}
	}
}
