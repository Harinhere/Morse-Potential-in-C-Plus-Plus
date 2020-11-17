#include<iostream>
#include<math.h>
#include<fstream>
#include "nr3.h"
using namespace std;
// Global variable declaration
//1D Morse potential parameters
//anion term
double beta=.82,B=0.0928,C=0.009,D=-0.011;
//neutral term and equilibrium separation
double alpha=.8507,A=0.1383,r1e=3.324;
//mass of the particle
double sqm=sqrt(85462.9378),mass=sqm*sqm/2.;
double pi=2.*asin(1.0);



// Trajjectory ccomputation
//Given initial x,p: This function computes
// x(t),p(t),R(x,x(t)),dx(t)/dp,dp(t)/dp



void derivs(const Doub x, VecDoub_I &y, VecDoub_O &dydx)
{
	dydx[0]=y[1];
	dydx[1]=(2.*B*beta*exp(-2.*beta*x)-beta*C*exp(-beta*x))/mass;
}

void rk4(VecDoub_I &y, VecDoub_I &dydx, const Doub x, const Doub h,
	VecDoub_O &yout, void derivs(const Doub, VecDoub_I &, VecDoub_O &))
{
	Int n=y.size();
	VecDoub dym(n),dyt(n),yt(n);
	Doub hh=h*0.5;
	Doub h6=h/6.0;
	Doub xh=x+hh;
	for (Int i=0;i<n;i++) yt[i]=y[i]+hh*dydx[i];
	derivs(xh,yt,dyt);
	for (Int i=0;i<n;i++) yt[i]=y[i]+hh*dyt[i];
	derivs(xh,yt,dym);
	for (Int i=0;i<n;i++) {
		yt[i]=y[i]+h*dym[i];
		dym[i] += dyt[i];
	}
	derivs(x+h,yt,dyt);
	for (Int i=0;i<n;i++)
		yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
}

//derivative function


void traject(const Doub x,const Doub p,const Doub tstep, VecDoub_O &xp )
{
	VecDoub y(2),dydx(2),yout(2);
	y[0]=x,y[1]=p/mass;
	derivs(x,y,dydx);
	rk4(y,dydx,x,tstep,yout,derivs);
	xp[0]=yout[0];
	xp[1]=yout[1]*mass;
}

// Morse potential wave function
double morse(int nv, double x)
{
	double r1,ak,omg,s2pi,s,y,z,r2,pf,ev;
    double t1,t2,t3,t4,t5,t6,t7,t8;
	r1=nv+.50;
	ak=2.*sqm*sqrt(A)/alpha;
	omg=2.*alpha*sqrt(A)/sqm;
	
//eigen energy
	ev=r1*omg*(1.-r1/ak);
//initial eigen function for reference
	s2pi=sqrt(2.*pi);
	s=ak-2.*nv-1.;
	y=ak*exp(-alpha*x);
    z=ak-nv;
	r2=s*log(y)/2.-(z-.5)*log(z)/2.;
	r2=r2-y/2.+z/2.;
    pf=exp(r2)*sqrt(alpha*s/s2pi);
// expansion coefficients
    t1=1.;
    t2=(1./12.)/z;	  
	t3=(1./288.)/pow(z,2);
	t4=(139./51840.)/pow(z,3);
	t5=(571./2488320.)/pow(z,4);
	t6=(163879./209018880.)/pow(z,5);
	t7=(5246819./75246796800.)/pow(z,6);
	t8=(534703531./902961561600.)/pow(z,7);
	pf=pf/sqrt(t1+t2+t3-t4-t5+t6+t7-t8);
//loop for gamma function
	for(int i=1;i<=nv;i++)
	{
		double gn=i;
		pf=pf*sqrt(gn);
	}
//recursion relation for Laguerre
	double fm_1=0.;
	double fm0=1.;
	for(int i=1;i<=nv;i++)
	{
		double cf1,cf2,fm1;
		cf1=2.*i+s-1.-y;
		cf2=s+i-1.;
		fm1=(cf1*fm0-cf2*fm_1)/i;
		fm_1=fm0;
	    fm0=fm1;
	}
	return pf*fm0;		
}
int main()
{
	fstream wave;
	wave.open("wave.txt",ios::out);
//	for(double x=-0.5;x<=0.5;0.1)
//	{
//		x=x+0.01;
//	    double wv=morse(0,x);
//	    cout<<x<<"	"<<wv<<endl;
//		wave<<x<<"	"<<wv<<endl;
		
//	}
		
//	cout<<traject(0.,-1.)<<endl;
// plot a trajectory
	double x0=0.,p0=-1000.,tstep=1.,tmax=1000.;
	VecDoub xp(2);
	int itmax=tmax/tstep;
	for(int i=1;i<=itmax;i++){
		double time=tstep*i;
		traject(x0,p0,tstep,xp);
		x0=xp[0];
		p0=xp[1];
		wave<<time<<"\t"<<x0<<"\t"<<p0<<endl;
		
	}
	
}
    	  
      
	       
     
