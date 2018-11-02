#include <iostream>
#include <math.h>
#include <fstream>
#include <stdlib.h>

using namespace std;


const double muearth=3.986e5;
const double musun=1.327e11;
const double mujupiter=1.266e8;
const double deuropa = 6.71e5;
const double rearth=6.371e3;
const double dse=1.496e8;
const double dsj=7.785e8;
const double pi = 3.14159;

double modul(double a[2]);
double dot(double a[], double b[]);
double anomaly(double e, double th);
double pr1A (double h0);
double pr1B (double h0);
double pr2B (double h0, double dv0, double dv1, double hp);


int main()
{
    double thrust;
    ofstream fout;
    fout.open("pr1A.txt");
    for(int i= 400; i < 10000; i+=10)
    {
        thrust = pr1A(i);	
        fout<< i << ';' << thrust <<endl;
    }
    fout.close();
    fout.open("pr1B.txt");
    for(int i= 400; i < 10000; i+=10)
    {
        thrust = pr1B(i);	
        fout<< i << ';' << thrust <<endl;
    }
    fout.close();
    fout.open("pr2B.txt");
    for(int i= 80; i < 2430; i+=10)
    {
        thrust = pr2B(3939,4,-0.65,i);	
        fout<< i << ';' << thrust <<endl;
    }
    fout.close();
    //pr2B(3939,4,-0.65,450);
    return 0;
    
}

double pr1A (double h0)
{
	long double thrust;
    double a,b,c;
    
    a=musun*((2/dse)-(1/((dse+dsj)/2)));
    b=(muearth)/(rearth+h0);
    c=pow(sqrt(a)-sqrt(musun/dse),2);
    thrust=sqrt(c+(2*b))-sqrt(b);
    cout<<h0<<" "<<thrust<<endl;
    return thrust;
}

double pr1B (double h0)
{
	double thrust;
    double a,b,c;
    
    a= musun*((2/dse)-(1/((dse+dsj)/2)));
    b=(muearth)/(rearth+h0);
    c=pow(sqrt(a)-sqrt(musun/dse),2);
    thrust=sqrt(c+(2*b))-sqrt(b);
    return thrust;
}

double modul(double a[2])
{
	return sqrt(pow(a[0],2)+pow(a[1],2));
}

double dot(double a[], double b[])
{
	double s = 0;
	for (int i = 0; i < 2; ++i)
	{
		s += a[i]*b[i];
	}

	return s;
}

double anomaly(double e, double th)
{
	return acos( (e + cos(th))/(1 + e*cos(th)) );
}

double pr2B (double h0, double dv0, double dv1, double hp)
{
	double uearth = sqrt(musun/dse);

	double u1p = uearth + sqrt( pow(( sqrt(muearth/(rearth+ h0))  + dv0) ,2)  - 2*muearth/(rearth + h0));
	double e1 = (dse*pow(u1p,2))/musun - 1;
	double d1a = dse*(1 + e1)/(1 - e1);
	double a1 = (d1a+dse)/2;
	double u1a = u1p*dse/d1a;
	double u2a = u1a + dv1;
	double e2 = 1 - (d1a*pow(u2a,2))/musun;

	double d2p = pow(d1a*u2a,2)/(musun*(1+e2));
	double a2 = (d1a+d2p)/2;
	double aint = acos((a2*(1-pow(e2,2)) - dse) / (e2*dse));

	double uminus[2] =  {-sin(aint), e2+cos(aint)};
	for (int i=0;i<2;i++)
		uminus[i] = uminus[i]*sqrt(musun/(a2*(1 - pow(e2,2))));

	double uminmod = modul(uminus);

	double uearthvec[2] = {-sin(aint), cos(aint)};
	for (int i=0;i<2;i++)
		uearthvec[i] = uearthvec[i]*uearth;

	double vminus[2];
	for (int i=0;i<2;i++)
		vminus[i] = uminus[i] - uearthvec[i];

	double vminmod = modul(vminus);

	double ehyp = 1 + ((hp+rearth)*pow( vminmod ,2))/muearth;
	double adelta = 2*asin(1/ehyp);

	double rmatrix[2][2] = 
	{
		{cos(adelta), -sin(adelta)},
		{sin(adelta), cos(adelta)}
	};
	double vplus[2] = {
		rmatrix[0][0]*vminus[0] + rmatrix[0][1]*vminus[1],
		rmatrix[1][0]*vminus[0] + rmatrix[1][1]*vminus[1]
	};
	double vplusmod = modul(vplus);

	double uplus[2];
	for(int i=0; i<2; i++)
		uplus[i] = uearthvec[i] + vplus[i];

	double uplusmod = modul(uplus);
	double pos[2] = 
	{
		dse*cos(aint),
		dse*sin(aint)
	};
	double e3vec[2] = 
	{
		((pow(uplusmod,2) - musun/dse)*pos[0] -  dot(pos,uplus)*uplus[0] )/musun,
		((pow(uplusmod,2) - musun/dse)*pos[1] -  dot(pos,uplus)*uplus[1] )/musun
	};

	double e3 = modul(e3vec);

	double d3a = pow(uplusmod*dse,2)/(musun*(1-e3));

	double theta = acos (dot(e3vec,pos)/(modul(e3vec)*modul(pos)));

	double u3a = uplusmod*dse/d3a;
	double d3p = pow(uplusmod*dse,2)/(musun*(1+e3));
	double ujupiter = sqrt(musun/dsj);
	double u3j = uplusmod*dse/dsj;
	double dv2 = sqrt(pow(u3j-ujupiter,2)+ 2*mujupiter/deuropa) - sqrt(mujupiter/deuropa);
	// double intjptr = acos( ( pow(uplusmod*dse,2)/(musun*dsj) - 1))/e3;

	// cout<< "1- " << ((dse+d1a)/(2*dse)) << ' ' <<((dse+d1a)/(2*dse)*e1) << endl;
	// cout<< "2- " << ((d2p+d1a)/(2*dse)) << ' ' <<((d2p+d1a)/(2*dse)*e2) << endl;
	// cout<< "3- " << ((d3p+d3a)/(2*dse)) << ' ' <<((d3p+d3a)/(2*dse)*e3) << endl;
	// cout << "int " << aint << endl;
	// cout << "intjptr(novaref) " << intjptr << endl;
	
	// double t1 = pi*sqrt(pow( (d1a+dse)/2 ,3)/musun);
	// double an2 = anomaly(e2, aint);
	// double t2 = (pi + an2 - e2*sin(an2))*sqrt(pow( (d1a+d2p)/2 ,3)/musun);
	// double an3int = anomaly(e3, theta);
	// double intjptr = acos( ( pow(uplusmod*dse,2)/(musun*dsj) - 1))/e3;
	// double an3j = anomaly(e3, intjptr);
	// double t3 = (an3int - e3*sin(an3int) - an3j + e3*sin(an3j))*sqrt(pow((d3a+d3p)/2 ,3)/musun);
	
	// t1 /= (3600*24*365);
	// t2 /= (3600*24*365);
	// t3 /= (3600*24*365);

	// cout << "Tt = " << t1 + t2 + t3 << endl;
	// cout << "T1 = " << t1 << endl;
	// cout << "T2 = " << t2 << endl;
	// cout << "T3 = " << t3 << endl;
	if (d3a > dsj)
		return dv0+dv1+dv2;
	return 0;
}

