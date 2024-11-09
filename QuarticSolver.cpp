#include "includes.h"

/***********************************************************
THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
************************************************************/
QuarticSolver::QuarticSolver(void)
{

}

QuarticSolver::~QuarticSolver(void)
{
	//delete sInstance;
}

void QuarticSolver::calculateQuadRoots(double a, double b, double c, QuartSolution3& Q)
{
    double sqterm;
    sqterm = b*b - 4*a*c;
    if (sqterm<0)
	{ 
		//imaginary root. return t=1
        Q.u[0]=-1;
        Q.u[1]=-1;
    }
    else
	{
        sqterm=sqrt(sqterm);
        a*=2;
        b*=-1;
        Q.u[0]=(b+sqterm)/(a);
        Q.u[1]=(b-sqterm)/(a);
    }
}
int round(double a)
{
	return int(a + 0.5);
}

bool simpleFuzzyEq(const double & v1, const double & v2, const double &epsilon)
{
	const double diff = v1 - v2;

	return diff < epsilon && diff > -epsilon;
}

int compare_double( const void* a, const void* b )
{
	double* arg1 = (double*) a;
	double* arg2 = (double*) b;
	if( *arg1 < *arg2 ) return -1;
	else if( simpleFuzzyEq(*arg1 , *arg2, EPSILON) ) return 0;
	else return 1;
}  

//Solve a 4th Degree Polynomial (Quartic) equation for 0s.
//taken from a javascript webpage (Explicitly states in source that source may be reused in any way)
//uses the quartic formula! :)
//DOESN"T WORK!!!!!!!!!!!!!! produces a lot of error trajectories
void QuarticSolver::calculateQ(double aq, double bq, double cq, double dq, double eqin, QuartSolution3 & Q)
{
double eq;
double fq;
double gq;
//double hq;


// These are the squares of the floatiables used to calculate the 4 roots-->
double kq;
double lq;
double mq;
double nq;
double mq2;

double compsw;
double kqsw;
double lqsw;
double mqsw;

// Switch used in calculating REAL quartic roots)-->
double sw;

// floatiables for calculating REAL quartic roots)-->

double kans;
double lans;
double theta;


double x1;
double x2;
double x3;
double x4;

double x2a, x2b, x2c, x2d;
//double x1b, x1b2, x2b2, x3b, x3b2, x4b, x4b2;
//more:
double dnm;
double a, b, c, d, f, g, h, k, m, m2,n,n2,r,rc;
double calcy, calcp, calcr, calcq, calcx, calcmod;
double dnmsw;
int i;

// the 'q' suffix  denotes floatiables used in the quartic equation
for (i=0;i<4;i++) 
{
    Q.u[i]=-1.0; //set to complex solutions
}
compsw=0;
kqsw=0;
lqsw=0;
mqsw=0;
dnmsw=0;
sw=0;


dnm=aq;      //note: this assumes aq is non-zero.  Of course it should be (eval 0.25g!)

//Simplifying by dividing all terms by the aq term called 'dnm' meaning denominator
aq=bq/dnm;
bq=cq/dnm;
cq=dq/dnm;
dq=eqin/dnm;
//Which yields an equation of the form X^4 + AX^3 + BX^2 + CX + D = 0

eq= bq-((3*aq*aq)/8);
fq= cq+ (aq*aq*aq/8) -(aq*bq/2);
gq= dq- (3*aq*aq*aq*aq/1024) + (aq*aq*bq/16) - (aq*cq/4);

// SOLVING THE RESULTANT CUBIC EQUATION
// EVALUATING THE 'f'TERM

a=1; b=eq/2; c=((eq*eq)-(4*gq))/16; d= ((fq*fq)/64)*-1;

f = (((3*c)/a) - (((b*b)/(a*a))))/3;
//EVALUATING THE 'g'TERM

g = ((2*((b*b*b)/(a*a*a))-(9*b*c/(a*a)) + ((27*(d/a)))))/27;

//EVALUATING THE 'h'TERM
h = (((g*g)/4) + ((f*f*f)/27));
if (h > 0)
{
	compsw=2;
	m = (-(g/2)+ (sqrt(h)));
	// K is used because math.pow cannot compute negative cube roots?
	k=1;
	if (m < 0) 
	{
		k=-1;
	}
	else 
	{
		k=1;
	}

	//m2 = ( (m*k) * (1.0/3.0));
	//m2 = std::pow( (m*k) ,(1.0/3.0));
	long double mtk = (m*k);
	long double oneOverThree = (1.0/3.0);
	m2 = ::pow( mtk , oneOverThree);

	m2 = m2*k;
	k=1;
	n = (-(g/2)- (sqrt(h)));
	if (n < 0)
	{
		k=-1;
	}
	else
	{
		k=1;
	}

	//n2 = (n*k)*(1.0/3.0);
	long double ntk = n*k;
	n2 = pow(ntk,oneOverThree);

	n2 *= k;
	k=1;
	kq=  ((m2 + n2) - (b/(3*a)));
	kq=sqrt(kq);
// ((S+U)     - (b/(3*a)))
calcmod= sqrt((-1*(m2 + n2)/2 - (b/(3*a)))*(-1*(m2 + n2)/2 - (b/(3*a))) + (((m2 - n2)/2)*sqrt(3.0))*(((m2 - n2)/2)*sqrt(3.0)));
calcy=sqrt((calcmod-(-1*(m2 + n2)/2 - (b/(3*a))))/2);
calcx=(((m2 - n2)/2)*sqrt(3.0))/(2*calcy);
calcp=calcx+calcy;
calcq=calcx-calcy;
calcr=kq;

nq=(aq/4);
x1=kq+calcp+calcq-nq;
x4=kq-calcp-calcq-nq;


Q.u[0]=-x1; //appearently was incorrect by a factor of -1
Q.u[1]=-1; //complex
Q.u[2]=-1; //complex
Q.u[3]=-x4;
}


// FOR H < 0

if (h<=0){
r = sqrt((g*g/4)-h);
k=1;

if (r<0)
  k=-1;
// rc is the cube root of 'r'

long double rtk = r*k;
long double oneOverThree = 1.0/3.0;
rc = pow(rtk,oneOverThree)*k;
k=1;
theta =acos((-g/(2*r)));

kq= (2*(rc*cos(theta/3))-(b/(3*a)));

x2a=rc*-1;
x2b= cos(theta/3.0);
x2c= sqrt(3.0)*(sin(theta/3));
x2d= (b/3.0*a)*-1;

lq=(x2a*(x2b + x2c))-(b/(3*a));

mq=(x2a*(x2b - x2c))-(b/(3*a));

nq=(aq/4.0);
}

if (h<=0){

// psudo-fix 0 bug.. not the best.. but works
if (abs(kq)<1.0/(10000.0))
  kq=0;
if (abs(lq)<1.0/(10000.0))
  lq=0;
if (abs(mq)<-1.0/(10000.0))
  mq=0;
if (kq<0) {return;} else {kq=sqrt(kq);}
if (lq<0) {return;} else {lq=sqrt(lq);}
if (mq<0) {return;} else {mq=sqrt(mq);}

if (kq*lq>0){mq2=((fq*-1)/(8*kq*lq));kans=kq;lans=lq;}
if (kq*mq>0){mq2=((fq*-1)/(8*kq*mq));kans=kq;lans=mq;}
if (lq*mq>0){mq2=((fq*-1)/(8*lq*mq));kans=lq;lans=mq;}




if (compsw==0){
  x1=kans+lans+mq2-nq;
  Q.u[0]=x1;
  x2=kans-lans-mq2-nq;
  Q.u[1]=x2;
  x3=(kans*-1)+lans-mq2-nq;
  Q.u[2]=x3;
  x4=(kans*-1)-lans+mq2-nq;
  Q.u[3]=x4;
}
}
}


//std::vector<double> Targeting::rootsOf(double A, double B, double C, double D, double E)
void QuarticSolver::calculateQ2(double A, double B, double C, double D, double E, QuartSolution3 & Q)
{
	std::vector<double> roots;// = new List<double>();
	if (A == 0.0)
	{
		if (B == 0.0)
		{
			if (C == 0.0)
			{
				if (D == 0.0)
				{
					if (E == 0.0)
						roots.push_back(0);//was Add(0.0);
				}
				else
				{
					roots.push_back(-D / E);
				}
			}
			else
			{
				roots.push_back( (-D + sqrt(D * D - 4 * C * E)) / (2.0 * C) );
				roots.push_back( (-D - sqrt(D * D - 4 * C * E)) / (2.0 * C) );
			}
		}
		else
		{
			C /= B;
			D /= B;
			E /= B;
			double F = (3.0 * D - C * C) / 3.0;
			double G = (2.0 * C * C * C - 9.0 * C * D + 27.0 * E) / 27.0;
			double H = (G * G) / 4.0 + (F * F * F) / 27.0;
			if (H > 0)
			{
				double intermediate = -G / 2.0 + sqrt(H);
				double m = intermediate < 0.0 ? -pow(-intermediate, 1.0 / 3.0) : pow(intermediate, 1.0 / 3.0);
				intermediate -= 2.0 * sqrt(H);
				double n = intermediate < 0.0 ? -pow(-intermediate, 1.0 / 3.0) : pow(intermediate, 1.0 / 3.0);
				roots.push_back(m + n - C / 3.0);
			}
			else
			{
				double intermediate = sqrt(G * G / 4.0 - H);
				double rc = intermediate < 0.0 ? -pow(-intermediate, 1.0 / 3.0) : pow(intermediate, 1.0 / 3.0);
				double theta = acos(-G / (2.0 * intermediate)) / 3.0;

				roots.push_back( 2.0 * rc * cos(theta) - C / 3.0);
				roots.push_back(-rc * (cos(theta) + sqrt(3.0) * sin(theta)) - C / 3.0);
				roots.push_back(-rc * (cos(theta) - sqrt(3.0) * sin(theta)) - C / 3.0 );
			
			}
			if (F + G + H == 0.0)
			{
				double intermediate = E < 0.0 ? pow(-E, 1.0 / 3.0) : -pow(E, 1.0 / 3.0);
				roots.clear();
				roots.push_back( intermediate);
				roots.push_back( intermediate);
				roots.push_back( intermediate);
			}
		}
	}
	else
	{
		B /= A;
		C /= A;
		D /= A;
		E /= A;
		double F = C - (3.0 * B * B) / 8.0;
		double G = D + B * B * B / 8.0 - (B * C) / 2.0;
		double H = E - 3.0 * B * B * B * B / 1024.0 + B * B * C / 16.0 - B * D / 4.0;
		double b = F / 2.0;
		double c = (F * F - 4.0 * H) / 16.0;
		double d = (G * G) / -64.0;
		double f = (3.0 * c - b * b) / 3.0;
		double g = (2.0 * b * b * b - 9.0 * b * c + 27.0 * d) / 27.0;
		double h = (g * g) / 4.0 + (f * f * f) / 27.0;
		double y1;
		double y2r;
		double y2i;
		double y3r;
		double y3i;
		if (h > 0.0)
		{
			double intermediate = -g / 2.0 + sqrt(h);
			double m = intermediate < 0.0 ? -pow(-intermediate, 1.0 / 3.0) : pow(intermediate, 1.0 / 3.0);
			intermediate -= 2.0 * sqrt(h);
			double n = intermediate < 0.0 ? -pow(-intermediate, 1.0 / 3.0) : pow(intermediate, 1.0 / 3.0);
			y1 = m + n - b / 3.0;
			y2r = (m + n) / -2.0 - b / 3.0;
			y2i = ((m - n) / 2.0) * sqrt(3.0);
			y3r = (m + n) / -2.0 - b / 3.0;
			y3i = ((m - n) / 2.0) * sqrt(3.0);
		}
		else
		{
			double intermediate = sqrt((g * g / 4.0 - h));
			double rc = intermediate < 0.0 ? -pow(-intermediate, 1.0 / 3.0) : pow(intermediate, 1.0 / 3.0);
			double theta = acos((-g / (2.0 * intermediate))) / 3.0;
			y1 = 2.0 * rc * cos(theta) - b / 3.0;
			y2r = -rc * (cos(theta) + sqrt(3.0) * sin(theta)) - b / 3.0;
			y2i = 0.0;
			y3r = -rc * (cos(theta) - sqrt(3.0) * sin(theta)) - b / 3.0;
			y3i = 0.0;
		}
		if (f + g + h == 0.0)
		{
			double intermediate = d < 0.0 ? pow(-d, 1.0 / 3.0) : -pow(d, 1.0 / 3.0);
			y1 = intermediate;
			y2r = intermediate;
			y2i = 0.0;
			y3r = intermediate;
			y3i = 0.0;
		}
		double p;
		double q;
		if (h <= 0.0)
		{
			int zeroCheck = 0;
			double cubicRoots[] = { y1, y2r, y3r };
			//Array.Sort(cubicRoots);
			qsort( cubicRoots, 3, sizeof(double), &compare_double );         

			p = sqrt(cubicRoots[1]);
			q = sqrt(cubicRoots[2]);

			if (round(y1) == 0.0)
			{
				p = sqrt(y2r);
				q = sqrt(y3r);
				zeroCheck = 1;
			}
			if (round(y2r) == 0.0)
			{
				p = sqrt(y1);
				q = sqrt(y3r);
				zeroCheck += 2;
			}
			if (round(y3r) == 0.0)
			{
				p = sqrt(y1);
				q = sqrt(y2r);
				zeroCheck += 4;
			}
			switch (zeroCheck)
			{
				case (3):
					p = sqrt(y3r);
					break;
				case (5):
					p = sqrt(y2r);
					break;
				case (6):
					p = sqrt(y1);
					break;
			}
			if (round(y1) < 0.0 || round(y2r) < 0.0 || round(y3r) < 0.0)
			{
				if (E == 0.0)
					roots.push_back(0.0);
			}
			else
			{
				double r;
				if (zeroCheck < 5)
				{
					r = G / (-8.0 * p * q);
				}
				else
				{
					r = 0.0;
				}
				double s = B / 4.0;
				roots.push_back(p + q + r - s);
				roots.push_back(p - q - r - s);
				roots.push_back(-p + q - r - s);
				roots.push_back(-p - q + r - s);
			}
		}
		else
		{
			double r2mod = sqrt(y2r * y2r + y2i * y2i);
			double y2mod = sqrt((r2mod - y2r) / 2.0);
			double x2mod = y2i / (2.0 * y2mod);
			p = x2mod + y2mod;
			double r3mod = sqrt(y3r * y3r + y3i * y3i);
			double y3mod = sqrt((r3mod - y3r) / 2.0);
			double x3mod = y3i / (2.0 * y3mod);
			q = x3mod + y3mod;
			double r = G / (-8.0 * (x2mod * x3mod + y2mod * y3mod));
			double s = B / 4.0;
			roots.push_back(x2mod + x3mod + r - s);
			roots.push_back(-x2mod - x3mod + r - s);
		}
	}


	int rCount = roots.size();
	int vIndex = 0;
	//for (int i = 0; i != roots.Count; i++)
	std::vector<double>::iterator rootsIterator;
	for (rootsIterator = roots.begin(); 
		 rootsIterator != roots.end();vIndex++)
	{
		double root = (*rootsIterator);
		//if (double.IsInfinity(roots[i]) || double.IsNaN(roots[i]))
		if (!_finite(root) || _isnan(root))
		{
			//roots.RemoveAt(i--);
			rootsIterator = roots.erase(rootsIterator);
			//invalidRoots.push_back(i);
		}
		else
		{
			rootsIterator++;
		}
	}

	rCount = roots.size();

	std::sort(roots.begin(), roots.end());
	int i = 0;
	std::vector<double>::iterator rootsIterator2;
	for (rootsIterator2 = roots.begin(); 
		 rootsIterator2 != roots.end();
		 rootsIterator2++)
	{
		double root = (*rootsIterator2);
		Q.u[i++] = root;
	}

}

/*Calculate aiming ideal rotation for firing a projectile at a potentially moving target (assumes pawn physics)
 IN:
 -StartLoc = world location where projectile is starting at
 -EndLoc = world Location we wish to Target (should lie in the targetted actor)
 -ProjSpeed = speed of the projectile being fired
 -Gravity = a vector describing the gravity
 -Target = the actual targetted ACTOR
 -bLeadTarget = Can we track the target?  (the entire point of this void)
 OUT:
 -dest: Location where the projectile will collide with Target
 -returns vector describing direction for projectile to leave at
*/
Vector3 QuarticSolver::-  GetSolutionVectort(Vector3 StartLoc, Vector3 EndLoc, double ProjSpeed, 
										   Vector3 Gravity, Vector3 targetVelocity, 
										   bool FallingPhysics, bool bLeadTarget, Vector3& Dest)
{
  QuartSolution3 Q;
  double best = 0;
  double speed2D = 0;
  double HitTime = 0;
  Vector3 Pr;
  int i;
  Vector3 HitNorm(1,0,0);
  Vector3 HitLoc;
    Vector3 D; //EndLoc-StartLoc
    Vector3 V; //Target.velocity

	//aTarget.FallingPhysics = false;

    D = EndLoc-StartLoc;
    V = targetVelocity;

    best=0;


  if (bLeadTarget && FallingPhysics )//== Phys_Falling)
  {
	calculateQuadRoots(V.Dot( V )- ProjSpeed*ProjSpeed,2*(V.Dot( D ) ),D.Dot(D) ,Q); //use quadratic formula
	for (i=0;i<2;i++)
	{
		if (best <= 0 ||(Q.u[i] > 0 && Q.u[i] < best))
		{
			best=Q.u[i];
		}
	}
	//Pr = normal(D/best + V)*ProjSpeed;
	Pr = ((D/best + V)).Normalize() * ProjSpeed;
	//if (best<=0 || aTarget.Trace(HitLoc,HitNorm,EndLoc+V*best+0.5*Gravity*best *best,EndLoc+vect(1,1,0)*V*best) == none)
	if ( best <= 0 || true )//|| aTarget.Trace(HitLoc,HitNorm,EndLoc+V*best+0.5*Gravity*best *best,EndLoc+vect(1,1,0)*V*best))
	{

	  //will be falling:
		Dest = StartLoc + Pr*best+0.5*Gravity*best*best;
		  gBestFiringSolution = best;

		return Pr.Normalize()*ProjSpeed;
	}
	else if (best>0)  //determine how long actor will be in air
	{

		HitTime = Vector3(HitLoc - (EndLoc+Vector3(1,1,0)*V*best)).Length()
			/Vector3(Vector3(0,0,1)*V*best+0.5*Gravity*best).Length();

	}
	else
		HitTime = 0; //assume most time not in air?
  }


		//ASSUME GROUND TRACKING
		if (bLeadTarget && FallingPhysics)
		{   
			  //trace down from target to get ground normal
			//Target.Trace(HitLoc,HitNorm,EndLoc+normal(Gravity)*5000,EndLoc);
			//D.z=HitLoc.z-StartLoc.Z;  //set destination.z to floor, wipe out velocity.z and re-eval assuming ground
			D.z=EndLoc.z-StartLoc.z;  //set destination.z to floor, wipe out velocity.z and re-eval assuming ground
	
			V.z=0;    //no longer falling - view velcocity in 2D
			if (HitTime>0.5)
			  {  
				  //True if likely in air most of time (in which case keep current V.X and V.y)
				V.z -= HitNorm.z * (V.Dot( HitNorm ));
			}
			else
			  { 
				  //otherwise alter all of velocity vector, but keep current 2D speed
				speed2D = V.Length();
				V = V.Normalize()*speed2D; //assume the same x and y speed if in air most time
				V -= HitNorm * (V.Dot( HitNorm ));   //recalculate players velocity on a slope using hitnormal  (assumes v.x and v.y is "ground speed")
				V=V.Normalize()*speed2D; //assume the same x and y speed if in air most time
			}
		}

    if (bLeadTarget && V != Vector3(0,0,0) )
	{
        calculateQ2(0.25*(Gravity.Dot( Gravity ) ),(-Gravity).Dot(V ) ,(-Gravity).Dot( D  )+
          V.Dot( V  )- ProjSpeed*ProjSpeed,2*(V.Dot(D ) ),D.Dot(D),Q);
        for (i=0;i<4;i++)
		{
			if (best<=0||(Q.u[i]>0 && Q.u[i]<best))
			{
				best=Q.u[i];
			}
		}
    }
    else
	{
		//don't lead. assume stationary target
		calculateQuadRoots(0.25*(Gravity.Dot(Gravity ) ),(-Gravity).Dot( D ) - ProjSpeed*ProjSpeed,D.Dot( D ),Q);
        for (i=0;i<2;i++)
		{
			if (best<=0||(Q.u[i]>0 && Q.u[i]<best))
			{
				best=Q.u[i];
			}
		}

		if (best>0)
		{
			best=sqrt(best);
		}
  }

  if (best<=0)
  { 
	  //projectile is out of range
      //Warning: Out of range adjustments assume gravity is parallel to the z axis and pointed downward!!
        Pr.z =ProjSpeed/sqrt(2.0); //determine z direction of firing
        best = -2*Pr.z/Gravity.z;
		best+=(D.Length()-Pr.z*best)/ProjSpeed; //note p.z = 2D vsize(p)  (this assumes ball travels in a straight line after bounce)
        
		//now recalculate PR to handle velocity prediction (so ball at least partially moves in direction of player)
		Pr = D/best + V - 0.5*Gravity*best;
        //now force maximum height again:
        Pr.z=0;
		Pr = (ProjSpeed/sqrt(2.0))*Pr.Normalize();
        Pr.z = ProjSpeed/sqrt(2.0); //maxmimum
        Dest = StartLoc + Pr*best+0.5*Gravity*best*best;
		  gBestFiringSolution = best;

        return Pr;
		//return Vector3(0,0,0);
  }

  Pr = Vector3(D/best + V - 0.5*Gravity*best).Normalize()*ProjSpeed;

  gBestFiringSolution = best;

    Dest = StartLoc + Pr*best+0.5*Gravity*best*best;
    return Pr;
}
