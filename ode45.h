
#ifndef __ODE45_H__
#define __ODE45_H__

#include <iostream>
#include <stdint.h>
#include <iomanip>

#include "Eigen/Core"

typedef Eigen::Matrix< double, Eigen::Dynamic, 1> VectorXF;




using namespace std;

class ode45
{
public:
	// Dormand Prince (DOPRI) coefficents
	static const double c2, c3, c4, c5, c6, c7,
						a21,
						a31, a32,
						a41, a42, a43,
						a51, a52, a53, a54,
						a61, a62, a63, a64, a65,
						a71, a72, a73, a74, a75, a76,
						b11, b13, b14, b15, b16,
						b21, b23, b24, b25, b26, b27,
						e1, e2, e3, e4, e5, e6, e7;

	ode45(double abs_tol = 1e-6, double rel_tol=1e-6, int maxSteps = 250);

	// Initialize integrator
	void init(double t0, double tf, VectorXF & y0);

	// Compute an initial step size h using y'(t)

	template<typename float_t, typename int_t>
    float_t machine_eps();

	template <typename ODEFUNC, typename USERDATA>
	void estimateInitStep(ODEFUNC func, USERDATA params);
protected:
	double m_abs_tol, m_rel_tol;
	double pow;
	double m_hmin, m_hmax;
	double m_t;
	VectorXF m_y;
	int m_dim;
	int m_maxSteps;
	double m_eps; //machine epsilon


};

template <typename ODEFUNC, typename USERDATA>
void ode45::estimateInitStep(ODEFUNC func, USERDATA params)
{
	VectorXF yprime = func(m_t, m_y, params);
	cout<<yprime<<endl;
}

template<typename float_t, typename int_t>
float_t ode45::machine_eps()
{
    union
    {
        float_t f;
        int_t   i;
    } one, one_plus, little, last_little;

    one.f    = 1.0;
    little.f = 1.0;
    last_little.f = little.f;

    while(true)
    {
        one_plus.f = one.f;
        one_plus.f += little.f;

        if( one.i != one_plus.i )
        {
            last_little.f = little.f;
            little.f /= 2.0;
        }
        else
        {
            return last_little.f;
        }
    }
}

#endif
