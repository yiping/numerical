#include "ode45.h"

using namespace Eigen;
//struct parameters {
//  double m;
//  double len;
//} params;

VectorXF f1(double t, VectorXF y, PARAMS_PTR p)
{
    VectorXF v(4);
    v(0) = t;
    v(1) = y(0);
    v(2) = y(1);
    v(3) = p->m;

    return v;

}

int main()
{
	ode45 integrator;
	VectorXF v = VectorXF::Zero(4);
	v(0)=7.0;
	v(1)=8.5;
	PARAMS params;
	PARAMS * p = &params;
	p->m=1.05;
	integrator.init(0.2, 2.0, v, f1, p);

	Vector3d v3(1, -2, -3);
	Vector3d vv(-1, 4, 1);
	cout<<v3.array().abs()<<endl<<endl;
	cout<<v3.array().abs().max(vv.array())<<endl<<endl;
	cout<<(v3.array()/vv.array()).maxCoeff()<<endl<<endl;
	cout<<max(1,2)<<endl;

	while(1){}

    //integrator.estimateInitStep(f1, params);
}

