#include "ode45.h"

struct parameters {
  double m;
  double len;
} params;

VectorXF f1(double t, VectorXF y, parameters p)
{
    VectorXF v(4);
    v(0) = t;
    v(1) = y(0);
    v(2) = y(1);
    v(3) = p.m;

    return v;

}

int main()
{
	ode45 integrator;
	VectorXF v = VectorXF::Zero(4);
	v(0)=7.0;
	v(1)=8.5;
	integrator.init(0.2, 2.0, v);
    params.m=1.05;

    integrator.estimateInitStep(f1, params);
}

