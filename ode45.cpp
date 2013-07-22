#include "ode45.h"
#include <math.h>

const double ode45::c2 = 1.0/5.0;
const double ode45::c3 = 3.0/10.0;
const double ode45::c4 = 4.0/5.0;
const double ode45::c5 = 8.0/9.0;
const double ode45::c6 = 1.0;
const double ode45::c7 = 1.0;

const double ode45::a21 = 1.0/5.0;

const double ode45::a31 = 3.0/40.0;
const double ode45::a32 = 9.0/40.0;

const double ode45::a41 = 44.0/45.0;
const double ode45::a42 = -56.0/15.0;
const double ode45::a43 = 32.0/9.0;

const double ode45::a51 = 19372.0/6561.0;
const double ode45::a52 = -25360.0/2187.0;
const double ode45::a53 = 64448.0/6561.0;
const double ode45::a54 = -212.0/729.0;

const double ode45::a61 = 9017.0/3168.0;
const double ode45::a62 = -355.0/33.0;
const double ode45::a63 = 46732.0/5247.0;
const double ode45::a64 = 49.0/176.0;
const double ode45::a65 = -5103.0/18656.0;

const double ode45::a71 = 35.0/384.0;
const double ode45::a72 = 0.0;
const double ode45::a73 = 500.0/1113.0;
const double ode45::a74 = 125.0/192.0;
const double ode45::a75 = -2187.0/6784.0;
const double ode45::a76 = 11.0/84.0;

const double ode45::b11 = 35.0/384.0;
const double ode45::b13 = 500.0/1113.0;
const double ode45::b14 = 125.0/192.0;
const double ode45::b15 = -2187.0/6784.0;
const double ode45::b16 = 11.0/84.0;


const double ode45::b21 = 5179.0/57600.0;
const double ode45::b23 = 7571.0/16695.0;
const double ode45::b24 = 393.0/640.0;
const double ode45::b25 = -92097.0/339200.0;
const double ode45::b26 = 187.0/2100.0;
const double ode45::b27 = 1.0/40.0;

const double ode45::e1 = 71.0/57600.0;
const double ode45::e3 = -71.0/16695.0;
const double ode45::e4 = 71.0/1920.0;
const double ode45::e5 = -17253.0/339200.0;
const double ode45::e6 = 22.0/525.0;
const double ode45::e7 = -1.0/40.0;

ode45::ode45(double abs_tol, double rel_tol, int maxSteps): m_abs_tol(abs_tol), m_rel_tol(rel_tol), m_maxSteps(maxSteps)
{
    cout << "ode45 (DOPRI) integrator is created" << endl;
    m_pow = 0.2;
}

void ode45::init(double t0, double tf, VectorXF & y0, ODEFUNC f, PARAMS_PTR p)
{
	m_t = t0;
	m_dim = y0.size();
	m_y.resize(m_dim);
	m_y = y0;

	m_odefunc = f;
	m_params_ptr = p;

    m_eps = machine_eps<double, uint64_t>(t0);
    cout<<"machine epsilon is: "<<setprecision(18)<<m_eps<<endl;

	m_hmin = 16*m_eps;
	if (tf>t0)
		m_hmax = 0.1*(tf - t0);
    else
    {
        cerr<<"tf needs to be greater than t0. "<<endl;
        exit(1);
    }

    m_derivative = m_odefunc(m_t, m_y, m_params_ptr);
    cout<<"y'(0) = "<<endl;
    cout<<setprecision(6)<<m_derivative<<endl;

    // estimate an initial step
    double threshold = m_abs_tol/m_rel_tol;
    VectorXF threshold_vec = threshold*VectorXF::Ones(m_dim);
    VectorXF temp = m_y.cwise().abs().max(threshold_vec);
    double rh = m_derivative.cwise()/temp/0.8*pow(m_rel_tol, m_pow);




}





