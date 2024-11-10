#pragma once
/***********************************************************
THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
************************************************************/
////////////////////////////////////////////////////////////
//solutions to quart!
class QuartSolution3
{ 
public:
	double u[4];
	QuartSolution3()
	{
		//memset(u,0,4*sizeof(double));
		u[0] = 0;
		u[1] = 0;
		u[2] = 0;
		u[3] = 0;
	}

};

typedef struct TargetTypeDef3
{
    Vector3 Velocity;
	bool FallingPhysics;
}target3;

class QuarticSolver
{
private:
	//static QuarticSolver* sInstance;
public:

	QuarticSolver(void);
	~QuarticSolver(void);
	//static QuarticSolver& instance();

	void calculateQuadRoots(double a, double b, double c, QuartSolution3& Q);
	void calculateQ(double aq, double bq, double cq, double dq, double eqin, QuartSolution3 & Q);
	void calculateQ2(double A, double B, double C, double D, double E, QuartSolution3 & Q);
	Vector3 GetShootVect(Vector3 StartLoc, Vector3 EndLoc,
		double ProjSpeed, Vector3 Gravity, Vector3 targetVelocity, bool FallingPhysics, bool bLeadTarget, Vector3& Dest);
};

int round(double a);
bool simpleFuzzyEq(const double & v1, const double & v2, const double &epsilon);
int compare_double( const void* a, const void* b );
