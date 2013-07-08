
#ifndef __ODE45_H__
#define __ODE45_H__

#include <Eigen/Core>

typedef Eigen::Matrix< double, Eigen::Dynamic, 1> VectorXF;

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
private:
	double m_abs_tol, m_rel_tol;
	
};

#endif