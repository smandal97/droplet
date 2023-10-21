enum{RHO,PPP,UUR};
enum{DEN,TAU,SSR};

#define NVAR  3
#define NG    2
#define GAM   1.666666666666667
#define D_FLR 0.
#define P_FLR 1e-8

struct zone{
	double vars[NVAR];
};

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>


void gridsetup( double * r, double rmin, double rmax, int num );
void initial( struct zone * theZones, double * r, double rmax , int num );
void boundary( struct zone * theZones, double * r, int num );

void flux( double * flux , double * prim );
double getdt( double * prim , double dr );
void prim2cons( double * prim , double * cons , double dV );
void cons2prim( double * cons , double * prim , double dV );
void getUstar( double * prim , double * Ustar , double Sk , double Ss );
void vel( double * Sl , double * Sr , double * Ss , struct zone * pL , struct zone * pR );
void add_source( struct zone * prim , struct zone * cons , double r , double dr , double dVdt );

void riemann( struct zone * pL , struct zone * pR , struct zone * cL , struct zone * cR , double dAdt );





int main( int argc , char * argv[] ){

	//parameters
	int n       = 10000;
	double rmin = 1e-2;
   	double rmax = 1.;
   	double tmin = 0.;
   	double tmax = 0.2;
   	double CFL  = 0.5;


   	//additional variables
   	double t    = tmin;
   	double dr   = (rmax-rmin)/(double)n;
   	double dt   = 1e10;

   	//timers
   	struct timeval start_time, stop_time, elapsed_time;

   	//main arrays
   	struct zone *restrict prim = (struct zone *)malloc( (n+NG)*sizeof(struct zone) );
    struct zone *restrict cons = (struct zone *)malloc( (n+NG)*sizeof(struct zone) );
    double *restrict ri    = (double *)malloc( (n+NG)*sizeof(double) );

    gettimeofday(&start_time,NULL); // Unix timer

   	//initialize and do boundaries
   	gridsetup(ri,rmin,rmax,n);
   	initial(prim,ri,rmax,n);
   	boundary(prim,ri,n);

   	//convergence loop
   	#pragma acc data copy(prim[:n-1],ri[:n-1]), create(cons[:n-1])
   	while(t<=tmax){

   		int i;
   		double dA,dV;
   		//compute conserved variables from primitive ones
   		#pragma acc parallel loop
   		for( i=0 ; i<n+NG ; ++i ){
   			dV = 4.*M_PI*(*(ri+i))*(*(ri+i))*dr + M_PI/3.*pow(dr,3.);
   			prim2cons( (*(prim+i)).vars , (*(cons+i)).vars , dV );
   		}
   		//reduction operation to find dt: DO REDUCTION HERE
   		#pragma acc parallel loop reduction(max:dt)
   		for( i=0 ; i<n+NG ; ++i )
   			//printf("%e\n",fmin(getdt( (*(prim+i)).vars,dr ),1e10));
   			dt = CFL*fmax( getdt( (*(prim+i)).vars,dr ) , dt );
   			//printf("dt = %e\n", dt);
   		//solve Riemann problem for each face
   		#pragma acc parallel loop
   		for( i=0 ; i<n+NG-1 ; ++i ){
   			dA = 4.*M_PI*(*(ri+i)+dr/2.)*(*(ri+i)+dr/2.);
   			riemann( (prim+i) , (prim+i+1) , (cons+i) , (cons+i+1) , dA*dt );
   		}
   		//add source and update primitive variables from conserved ones
   		#pragma acc parallel loop
   		for( i=0 ; i<n+NG ; ++i ){
   			dV = 4.*M_PI*(*(ri+i))*(*(ri+i))*dr + M_PI/3.*pow(dr,3.);
   			add_source( (prim+i) , (cons+i) , (*(ri+i)) , dr , dV*dt );
   			cons2prim( (*(cons+i)).vars , (*(prim+i)).vars , dV );
   		}

   		//reset boundary values

   		//update time and reset dt to high value
   		t += dt;
   		//printf("dt = %.5e\n", dt);
   		dt = 1e10;
   
   	}


   	gettimeofday(&stop_time,NULL);
    timersub(&stop_time, &start_time, &elapsed_time); // Unix time subtract routine


    int i;
    FILE * fp = fopen("output.dat","w");
    for( i=0 ; i<n ; ++i )
    	fprintf(fp,"%e %e %e %e\n", (*(ri+i)), (*(prim+i)).vars[RHO] , (*(prim+i)).vars[PPP] , (*(prim+i)).vars[UUR] );
      fclose(fp);


    //free arrays
    free(prim);
    free(cons);

   	printf("Total time was %f seconds.\n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);
   	return 0.;


}












void gridsetup( double * r, double rmin, double rmax, int num ){
	int i;
	double *r_l;
	for( i=0 ; i<num+NG ; ++i ){
		r_l  = r+i;
		*r_l = rmin + (rmax-rmin)*((double)i+0.5)/(double)num;
	}
}


void initial( struct zone * theZones, double * r, double rmax , int num ){

	int i;
	double rl;
	struct zone * z;
	for( i=0 ; i<num+NG ; ++i ){
		rl = *(r+i);
		z = theZones+i;
		z->vars[RHO] = 1.0;
		z->vars[PPP] = 1.0;
		if(rl>rmax/2.){
			z->vars[RHO] = 0.1;
			z->vars[PPP] = 0.125;
		}
		z->vars[UUR] = 0.0;
	}


}


void boundary( struct zone * theZones, double * r, int num ){

	struct zone * z = theZones+num;
	z->vars[RHO] = 0.1;
	z->vars[PPP] = 0.125;
	z->vars[UUR] = 0.0;
	z = theZones+num+1;
	z->vars[RHO] = 0.1;
	z->vars[PPP] = 0.125;
	z->vars[UUR] = 0.0;

}



void prim2cons( double * prim , double * cons , double dV ){
	double rho = prim[RHO];
	double Pp  = prim[PPP];
	double ur  = prim[UUR];

	double rhoe = Pp/(GAM-1.);

	cons[DEN] = rho*dV;
	cons[TAU] = (.5*rho*ur*ur+rhoe)*dV;
	cons[SSR] = rho*ur*dV;
}



void cons2prim(  double * cons , double * prim , double dV ){
	double rho = cons[DEN]/dV;
   	double Sr  = cons[SSR]/dV;
   	double E   = cons[TAU]/dV;

	double ur = Sr/rho;
   	double rhoe = E - .5*rho*ur*ur;
   	double Pp = (GAM-1.)*rhoe;

   	if( rho < D_FLR ) rho    = D_FLR;
   	if( Pp  < P_FLR*rho ) Pp = P_FLR*rho;

   	prim[RHO] = rho;
   	prim[PPP] = Pp;
   	prim[UUR] = ur;
}



double getdt( double * prim , double dr ){
	double rho = prim[RHO];
	double Pp  = prim[PPP];
	double ur  = prim[UUR];

	double vmax = sqrt(fabs(GAM*Pp/rho)) + fabs(ur);
	return (dr/vmax);
}




void flux( double * flux , double * prim ){
	double rho = prim[RHO];
   	double Pp  = prim[PPP];
   	double ur  = prim[UUR];
   	double rhoe = Pp/(GAM-1.);
 
   	flux[DEN] = rho*ur;
   	flux[SSR] = rho*ur*ur + Pp;
   	flux[TAU] = (.5*rho*ur*ur + rhoe + Pp )*ur;
}



void vel( double * Sl , double * Sr , double * Ss , struct zone * pL , struct zone * pR ){
	double rho1 = pL->vars[RHO];
	double p1   = pL->vars[PPP];
	double u1   = pL->vars[UUR];
	double cs1  = sqrt(fabs(GAM*p1/rho1));

	double rho2 = pR->vars[RHO];
	double p2   = pR->vars[PPP];
	double u2   = pR->vars[UUR];
	double cs2  = sqrt(fabs(GAM*p2/rho2));

	*Ss = ( p2 - p1 + rho1*u1*(-cs1) - rho2*u2*cs2 )/( rho1*(-cs1) - rho2*cs2 );

   	*Sr =  cs1 + u1;
   	*Sl = -cs1 + u1;

   	if( *Sr <  cs2 + u2 ) *Sr =  cs2 + u2;
   	if( *Sl > -cs2 + u2 ) *Sl = -cs2 + u2;
}




void getUstar( double * prim , double * Ustar , double Sk , double Ss ){

	double rho = prim[RHO];
   	double Pp  = prim[PPP];
   	double ur  = prim[UUR];
   	double u2  = ur*ur;

   	double rhoe = Pp/(GAM-1.);

   	double rhostar = rho*(Sk - ur)/(Sk - Ss);
   	double Pstar = Pp*(Ss - ur)/(Sk - Ss);
   	double Us = rhoe*(Sk - ur)/(Sk - Ss);

   	Ustar[DEN] = rhostar;
   	Ustar[TAU] = .5*rhostar*u2 + Us + rhostar*Ss*(Ss - ur) + Pstar;
   	Ustar[SSR] = rhostar*Ss;

}



void add_source( struct zone * prim , struct zone * cons , double r , double dr , double dVdt ){
   double Pp  = prim->vars[PPP];
   double rp  = r+dr/2.;
   double rm  = r-dr/2.;
   double r2  = (rp*rp+rm*rm+rp*rm)/3.;
   cons->vars[SSR] += 2.*Pp*(r/r2)*dVdt;
}










void riemann( struct zone * pL , struct zone * pR , struct zone * cL , struct zone * cR ,  double dAdt ){

	int q;
	double Sl,Sr,Ss;
	double F[NVAR];
	double Fk[NVAR];
	double Uk[NVAR];
	double Ustar[NVAR];
	vel( &Sl , &Sr , &Ss , pL , pR );

	if(Sl>0.){
		flux(F,pL->vars);
	}
	else if(0.>Sr){
		flux(F,pR->vars);
	}
	else{
		if(Ss>0.){
			flux(Fk,pL->vars);
			prim2cons(pL->vars,Uk,1.0);
			getUstar(pL->vars,Ustar,Sl,Ss);
			for( q=0; q<NVAR ; ++q )
				F[q] = Fk[q] + Sl*(Ustar[q]-Uk[q]);
		}
		else{
			flux(Fk,pR->vars);
			prim2cons(pR->vars,Uk,1.0);
			getUstar(pR->vars,Ustar,Sr,Ss);
			for( q=0; q<NVAR ; ++q )
				F[q] = Fk[q] + Sr*(Ustar[q]-Uk[q]);
		}

	}

	for( q=0; q<NVAR ; ++q ){
		cL->vars[q] -= F[q]*dAdt;
		cR->vars[q] += F[q]*dAdt;
	}

}
