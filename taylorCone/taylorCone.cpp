#include "stdafx.h"
#include "taylorCone.h"

#define SCAN_NZ 80
#define SCAN_NR 100

double TaylorCone::curv(double r, double z, double dr, double dz, double ddr, double ddz) {
	if (r > 1e-11) {
		return	(dr * ddz - dz * ddr) / pow(dr * dr + dz * dz, 1.5) + dz / r / sqrt(dr * dr + dz * dz);
	}
	else {
		return ddz / dr / dr * 2.;
	}
}

double TaylorCone::fD2(int order, double h0, double h1, double y0, double y1, double y2, int location) {
	switch (order)	{
	case 1: {
		switch (location) {
		case 0: {
			return (-1. / h0 - 1. / (h0 + h1))*y0 + (1. / h0 + 1 / h1)*y1 + (-1. / h1 + 1 / (h0 + h1))*y2;
			break;
		}
		case 1: {
			return (-1. / h0 + 1. / (h0 + h1))*y0 + (1. / h0 - 1. / h1)*y1 + (1. / h1 - 1. / (h0 + h1))*y2;
			break;
		}
		case 2: {
			return (1. / h0 - 1. / (h0 + h1))*y0 + (-1. / h0 - 1. / h1)*y1 + (1. / h1 + 1. / (h0 + h1))*y2;
			break;
		}
		}
		break;
	}
	case 2: {
		return (2.*y0) / (h0*(h0 + h1)) - (2.*y1) / (h0*h1) + (2. / (h0*h1) - 2. / (h0*(h0 + h1)))*y2;
		break;
	}
	default:
		return 0;
		break;
	}
};
double TaylorCone::fD3(int order, double h0, double h1, double h2, double y0, double y1, double y2, double y3, int location) {
	switch (order) {
	case 1: {
		switch (location) {
		case 0: {	
			return ((-1.*(2.*h0 + h1)) / (h0*(h0 + h1)) - 1. / (h0 + h1 + h2))*y0
				+ ((h0 + h1) / (h0*h1) + (h0 + h1) / (h1*(h1 + h2)))*y1
				+ ((-1.*h0) / (h1*(h0 + h1)) - (1.*h0) / (h1*h2))*y2
				+ (h0 / (h1*h2) + (-1.*h0 - 1.*h1) / (h1*(h1 + h2)) + 1 / (h0 + h1 + h2))*y3;
			break;
		}
		case 3: {	
			return (1 / (h0 + h1) - (1.*h2) / (h0*(h0 + h1)) - 1. / (h0 + h1 + h2))*y0
				+ (1 / h1 + h2 / (h0*h1) - 1. / (h1 + h2))*y1
				+ ((-1.*(h0 + 2.*h1)) / (h1*(h0 + h1)) - 1. / h2 - (1.*h2) / (h1*(h0 + h1)))*y2
				+ (1 / h2 + 1 / (h1 + h2) + 1 / (h0 + h1 + h2))*y3;
			break;
		}		
		}
		break;
	}
	case 2: {	
		switch (location) {
		case 0: {
			return (2. / (h0*(h0 + h1)) + (2.*(2.*h0 + h1)) / (h0*(h0 + h1)*(h0 + h1 + h2)))*y0
				+ (-2. / (h0*h1) - (2.*(2.*h0 + h1)) / (h0*h1*(h1 + h2)))*y1
				+ (2. / (h1*(h0 + h1)) + (2.*(2.*h0 + h1)) / (h1*(h0 + h1)*h2))*y2
				+ ((-2.*(2.*h0 + h1)) / (h1*(h0 + h1)*h2) + (2.*(2.*h0 + h1)) / (h0*h1*(h1 + h2)) - (2.*(2.*h0 + h1)) / (h0*(h0 + h1)*(h0 + h1 + h2)))*y3;
			
			break;
		}
		case 3: {
			return (-4. / (h0*(h0 + h1)) + (2.*(2.*h0 + h1)) / (h0*(h0 + h1)*(h0 + h1 + h2)))*y0
				+ (4. / (h0*h1) + (2.*(h0 - 1.*h1)) / (h0*h1*(h1 + h2)))*y1
				+ (-4. / (h1*(h0 + h1)) - (2.*(h0 + 2.*h1)) / (h1*(h0 + h1)*h2))*y2
				+ ((2.*(h0 + 2.*h1)) / (h1*(h0 + h1)*h2) - (2.*(h0 - 1.*h1)) / (h0*h1*(h1 + h2)) - (2.*(2.*h0 + h1)) / (h0*(h0 + h1)*(h0 + h1 + h2)))*y3;
			break;
		}
		}		
		return 0;
		break;
	}
	default:
		
		return 0;

		break;
	}	
	return 0;
};
double TaylorCone::fD4(int order, double h0, double h1, double h2, double h3, double y0, double y1, double y2, double y3, double y4) {
	switch (order)
	{
	case 1:
		return (h3*(h2 + h3)*(h1 + h2 + h3)*y0
			) / (h0*(h0 + h1)*(h0 + h1 + h2)*(h0 + h1 + h2 + h3)) - (1.*h3*(h2 + h3)*(h0 + h1 + h2 + h3)*y1
				) / (h0*h1*(h1 + h2)*(h1 + h2 + h3)) + (h3*(h1 + h2 + h3)*(h0 + h1 + h2 + h3)*y2
					) / (h1*(h0 + h1)*h2*(h2 + h3)) - (1.*(h2 + h3)*(h1 + h2 + h3)*(h0 + h1 + h2 + h3)*y3
						) / (h2*(h1 + h2)*(h0 + h1 + h2)*h3) + ((pow(h1, 2)*(h2 + 2.*h3) + pow(h2 + h3, 2)*(h2 + 4.*h3) + 2.*h1*(pow(h2, 2) + 4.*h2*h3 + 3.*pow(h3, 2)) + h0 * (pow(h2, 2) + 4.*h2*h3 + 3.*pow(h3, 2) + h1 * (h2 + 2.*h3)))*y4) / (h3*(h2 + h3)*(h1 + h2 + h3)*(h0 + h1 + h2 + h3));
	case 2:
		return (2.*(pow(h2, 2) + 4.*h2*h3 + 3.*pow(h3, 2) + h1 * (h2 + 2.*h3))*y0
			) / (h0*(h0 + h1)*(h0 + h1 + h2)*(h0 + h1 + h2 + h3)) - (2.*(pow(h2, 2) + 4.*h2*h3 + 3.*pow(h3, 2) + h0 * (h2 + 2.*h3) + h1 * (h2 + 2.*h3))*y1
				) / (h0*h1*(h1 + h2)*(h1 + h2 + h3)) + (2.*(pow(h1, 2) + pow(h2, 2) + 4.*h2*h3 + 3.*pow(h3, 2) + 2.*h1*(h2 + 2.*h3) + h0 * (h1 + h2 + 2.*h3))*y2
					) / (h1*(h0 + h1)*h2*(h2 + h3)) - (2.*(pow(h1, 2) + 4.*h1*(h2 + h3) + 3.*pow(h2 + h3, 2) + h0 * (h1 + 2.*(h2 + h3)))*y3
						) / (h2*(h1 + h2)*(h0 + h1 + h2)*h3) + (2.*(pow(h1, 2) + 4.*h1*h2 + 3.*pow(h2, 2) + 6.*h1*h3 + 9.*h2*h3 + 6.*pow(h3, 2) + h0 * (h1 + 2.*h2 + 3.*h3))*y4) / (h3*(h2 + h3)*(h1 + h2 + h3)*(h0 + h1 + h2 + h3));
	default:
		break;
	}
	
};

double TaylorCone::gridDistribution(double t, double q) {
	if (std::abs(t) < 1e-13) { return 0.0; }
	else if (std::abs(t - 1.) < 1e-13) { return 1.0; }
	else {
		//double v = 0.5 + (-0.5 + 1/q)*pow(1. - 1.*t,2) + (-1. + 1.*t)/q;
		//double v1 = 0.5 + (-0.5 + 1 / q)*pow(1. - 1.*(1. - 1.*t), 2) + (-1. + 1.*(1. - 1.*t)) / q;
		return 2. / (1. + pow(-1. + 2. / t, q));
	}

	
}
void TaylorCone::circleDerivativeBegin(const Eigen::MatrixX2d &xy, double &dx, double &ddx, double &dy, double &ddy) {	
	double x0 = xy(0, 0), y0 = xy(0, 1);
	double x1 = xy(1, 0), y1 = xy(1, 1);
	double x2 = xy(2, 0), y2 = xy(2, 1);
	double x3 = xy(3, 0), y3 = xy(3, 1);
	double h0 = sqrt((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0));
	double h1 = sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
	double h2 = sqrt((x3 - x2)*(x3 - x2) + (y3 - y2)*(y3 - y2));

	double curvature = -2.0 / sqrt(x0 * x0 + y0 * y0);
	double slope = std::abs(x0) / std::abs(y0);
	if (xy(xy.rows() - 1, 1) > 0) {	curvature = curvature * -1.;}
	
	dx = fD3(1, h0, h1, h2, x0, x1, x2, x3, 0);
	ddx = fD3(2, h0, h1, h2, x0, x1, x2, x3, 0);	
	dy = dx * slope;
	ddy = (ddx * dy + curvature * pow(dx * dx + dy * dy,1.5)) / dx - dy * (dx * dx + dy * dy) / dx / x0;
};

void TaylorCone::coneDerivativeEnd(const Eigen::MatrixX2d &xy, double c[5], double &dx, double &ddx, double &dy, double &ddy) {
	const int n = xy.rows();
	double xn1 = xy(n - 5, 0), yn1 = xy(n - 5, 1);
	double x0 = xy(n - 4, 0), y0 = xy(n - 4, 1);
	double x1 = xy(n - 3, 0), y1 = xy(n - 3, 1);
	double x2 = xy(n - 2, 0), y2 = xy(n - 2, 1);
	double x3 = xy(n - 1, 0), y3 = xy(n - 1, 1);


	double hn1 = sqrt((x0 - xn1)*(x0 - xn1) + (y0 - yn1)*(y0 - yn1));
	double h0 = sqrt((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0));
	double h1 = sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
	double h2 = sqrt((x3 - x2)*(x3 - x2) + (y3 - y2)*(y3 - y2));

	double r = x3;

	double curvature = (0.75*pow(r, 4.5)*c[1] + 15.75*pow(r, 1.5)*c[3] + pow(r, 6)*(1. + pow(c[0] + (-0.5*pow(r, 4.5)*c[1] - 3.5*pow(r, 1.5)*c[3] - 5.*c[4]) / pow(r, 6), 2))*(c[0] + (-0.5*pow(r, 4.5)*c[1] - 3.5*pow(r, 1.5)*c[3] - 5.*c[4]) / pow(r, 6)) + 30.*c[4]) / (pow(r, 7)*pow(1. + pow(c[0] + (-0.5*pow(r, 4.5)*c[1] - 3.5*pow(r, 1.5)*c[3] - 5.*c[4]) / pow(r, 6), 2), 1.5));
	double slope = c[0] - c[1] / (2.*pow(r, 1.5)) - (7. * c[3]) / (2.*pow(r, 4.5)) - (5. * c[4]) / pow(r, 6);
	
	//dx  = fD3(1, h0, h1, h2, x0, x1, x2, x3, 3);
	dx  = fD4(1, hn1, h0, h1, h2, xn1, x0, x1, x2, x3);
	//ddx = fD3(2, h0, h1, h2, x0, x1, x2, x3, 3);
	ddx = fD4(2, hn1, h0, h1, h2, xn1, x0, x1, x2, x3);
	dy = dx * slope;
	ddy = (ddx * dy + curvature * pow(dx * dx + dy * dy, 1.5)) / dx - dy * (dx * dx + dy * dy) / dx / x3;
}

void TaylorCone::computeCoefabc(double c1, double b0) {
	a[0] = 2.7142175397111330 * c1;
	b[0] = b0;
	c[0] = -0.8604366861256783;

	a[1] = 0.1442586135181731 * a[0] * a[0] - 0.4749808176761397 * b0 * b0 - c[0];
	b[1] = -0.848581976487259 * b0 * c1; // old
	c[1] = c1;

	//a[2] = -2.336057766096800 * c1 - 1.155514902883830 * c1 * c1* c1 + 1.584211046805990 * b0 * b0 * c1;//old	
	c[2] = 0.0;	

	//a[1] = 0.860436686125679 + 1.062749866616300 * c1 * c1	- 0.474980817676140 * b0 * b0;
	a[2] = -2.336057766096800 * c1 - 1.155514902883830 * c1 * c1* c1
		+ 1.584211046805990 * b0 * b0 * c1;
	a[3] = -0.433421293527112 + 1.563930669354330 * c1 * c1 + 1.356140305325190 * pow(c1, 4.0)
		+ 0.478517022152372 * b0 * b0 - 0.132076194633969 * pow(b0, 4.0) + 0.448670228959027 * b0 * b0 * c1 * c1;
	a[4] = -1.723725053118940 * c1 + 4.882389215855330 * pow(c1, 3.0) - 1.096390164067280 * pow(c1, 5.0)
		+ 1.980228762128670 * b0 * b0 * c1 - 0.409457597247434 * pow(b0, 4.0) * c1 - 17.405016393743000 * b0 * b0 * pow(c1, 3.0);


	c[2] = 0;
	c[3] = -0.275783438603136 * c1 + 0.210659453420660 * b0 * b0 * c1 - 0.069483517708871 * c1 * c1 *c1;
	c[4] = -0.045843694202325 + 0.050613544746824 * b0*b0 - 0.013969919726216 * pow(b0, 4.0) - 0.587515210204774 * c1* c1
		+ 0.579439247828955 * b0 * b0 *c1 * c1 - 0.139013862957991 * pow(c1, 4.0);

	//b[0] = b0;
	//b[1] = -0.848581976487259 * b0 * c1;
	b[2] = -1.290655029188520 * b0 * c1 * c1;
	b[3] = 1.887300354735200 * b0 * c1 - 1.441629936818880 * b0 * b0 *b0 * c1 - 5.508686000286550 * b0 * c1 * c1 * c1;
	b[4] = -0.957571779297224 * b0 + 1.057203241210380 * pow(b0, 3.0) - 0.291800238214520 * pow(b0, 5.0) + 16.866940785206700 * b0 * c1 * c1
		- 10.154738911572600 * pow(b0, 3.0)* c1* c1 - 64.833275833138900 * b0 * pow(c1, 4.0);

	
}

Eigen::MatrixX2d TaylorCone::generateCircle(double angle0, double angle1, double radius, int n) {
	Eigen::MatrixX2d xy(n, 2);
	xy.setZero();
	for (int i = 0; i < n; i++) {
		double t = ((double)i) / (n - 1.);

		double theta = angle0 + (angle1 - angle0) * t;// pow(t, 0.8);
		xy(i, 0) = radius * sin(theta) ;
		xy(i, 1) = radius * cos(theta) ;
	}
	if (std::abs(angle0) < 1e-12) { xy(0, 0) = 0.0; }
	if (std::abs(angle1 - M_PI) < 1e-12) { xy(n - 1, 0) = 0.0; }
	return xy;
}

Eigen::MatrixX2d TaylorCone::generateCircle(double r0, double z0, int end, int n) {
	Eigen::MatrixX2d xy(n, 2);
	xy.setZero();
	double radius = sqrt(r0* r0 + z0 * z0);
	double angle0 = acos(z0 / radius);
	double angle1 = 0.;
	
	xy(0, 0) = r0;
	xy(0, 1) = z0;
	if (end != 0) { angle1 = M_PI; };	
	
	for (int i = 1; i < n; i++) {
		double t = ((double)i) / (n - 1.);
		double theta = angle0 + (angle1 - angle0) * t;// pow(t, 0.8);
		xy(i, 0) = radius * sin(theta);
		xy(i, 1) = radius * cos(theta);
	}	
	xy(n - 1, 0) = 0;
	if (end == 0) { xy(n - 1, 1) = radius; }
	else { xy(n - 1, 1) = -radius; }
	return xy;
}

Eigen::MatrixX2d TaylorCone::generateCone(double rc, double rstar, double c[5], int n, double(*foo)(double t, double q)) {
	Eigen::MatrixX2d xy(n, 2);
	for (int i = 0; i < n; i++) {
		double t = ((double)i) / (n - 1.);
		if (foo != nullptr) {			
			t = (*foo)(t, 1.15);
		}		
		double r = rstar * t;
		xy(i, 0) = r;
		xy(i, 1) = c3Cone(r,rc,c);
		
	}
	
	return xy;
};

double TaylorCone::c3Cone(double r, double rc, double c[5]) {
	double f0, f2, f4, f6;
	f0 = (231 * c[4]) / (16.*pow(rc, 5)) + (1045 * c[3]) / (128.*pow(rc, 3.5)) + (195 * c[1]) / (128.*sqrt(rc)) + (5 * c[0] * rc) / 16.;
	f2 = (-495 * c[4]) / (16.*pow(rc, 7)) - (1995 * c[3]) / (128.*pow(rc, 5.5)) - (117 * c[1]) / (128.*pow(rc, 2.5)) + (15 * c[0]) / (16.*rc);
	f4 = (385 * c[4]) / (16.*pow(rc, 9)) + (1463 * c[3]) / (128.*pow(rc, 7.5)) + (65 * c[1]) / (128.*pow(rc, 4.5)) - (5 * c[0]) / (16.*pow(rc, 3));
	f6 = (-105 * c[4]) / (16.*pow(rc, 11)) - (385 * c[3]) / (128.*pow(rc, 9.5)) - (15 * c[1]) / (128.*pow(rc, 6.5)) + c[0] / (16.*pow(rc, 5));
	if (r < rc) {
		return f0 + f2 * r * r + f4 * r * r * r * r + f6 * r * r * r * r * r * r;
	}
	else {
		return c[0] * r + c[1] / sqrt(r) + c[3] / pow(r, 3.5) + c[4] / pow(r, 5.0);
	}
};

double TaylorCone::harmonicGrow(double r, double z, int l, int divide) {
	double R = sqrt(r * r + z * z);
	double cosTh = z / R;
	if (divide == 0) {
		return pow(R, l) * Numeric::legendreP(l, cosTh);
	}
	else {
		return pow(R, l/2.) * Numeric::legendreP(l,2,cosTh);
	}
};

double TaylorCone::dHarmonicGrowdR(double r, double z, int l, int divide) {
	double R = sqrt(r * r + z * z);
	double cosTh = z / R;
	if (divide == 0) {
		return l * pow(R, l - 1.) * Numeric::legendreP(l, cosTh);
	}
	else {
		return l/2. * pow(R, l / 2. - 1.) * Numeric::legendreP(l, 2, cosTh);
	}
};

double TaylorCone::harmonicDecay(double r, double z, int l, int divide) {
	double R = sqrt(r * r + z * z);
	double cosTh = z / R;
	if (divide == 0) {
		return pow(R, -1. -l) * Numeric::legendreP(l, cosTh);
	}
	else {
		return pow(R, -1. - l/2.) * Numeric::legendreP(l, 2, cosTh);
	}
};

double TaylorCone::dHarmonicDecaydR(double r, double z, int l, int divide) {
	double R = sqrt(r * r + z * z);
	double cosTh = z / R;
	if (divide == 0) {
		return (-1. - l) * pow(R, -2. - l) * Numeric::legendreP(l, cosTh);
	}
	else {
		return (-1. - l / 2.) * pow(R, -2. - l / 2.) * Numeric::legendreP(l, 2, cosTh);
	}
};

double TaylorCone::velocityPotentialFarField(double r, double z, const double (&a)[5]) {
	double flip = -1;
	return
		a[0] * harmonicGrow(r, flip * z, 1, 2) +
		a[1] * harmonicDecay(r, flip * z, 0) +
		a[2] * harmonicDecay(r, flip * z, 3, 2) +
		a[3] * harmonicDecay(r, flip * z, 3) +
		a[4] * harmonicDecay(r, flip * z, 5, 2);

};

double TaylorCone::electricPotentialFarField(double r, double z, const double(&b)[5]) {
	double flip = +1;
	return
		b[0] * dHarmonicGrowdR(r, flip * z, 1, 2) +
		b[1] * dHarmonicDecaydR(r, flip * z, 0) +
		b[2] * dHarmonicDecaydR(r, flip * z, 3, 2) +
		b[3] * dHarmonicDecaydR(r, flip * z, 3) +
		b[4] * dHarmonicDecaydR(r, flip * z, 5, 2);

};

void TaylorCone::prepareBem(int type, const Eigen::MatrixX2d &xy, int shift, Bem &bem) {

	bem.settings.indexShift(shift);
	bem.settings.order(2);
	bem.settings.qdOrder(20);
	double dx, ddx, dy, ddy;

	switch (type) {
	case 0: {
		bem.settings.xBC.begin.set(Spline::BC::Odd, 0., 0.);
		bem.settings.yBC.begin.set(Spline::BC::Even, 0., 0.);
		coneDerivativeEnd(xy, c, dx, ddx, dy, ddy);
		bem.settings.xBC.end.set(Spline::BC::Mix, dx, ddx);
		bem.settings.yBC.end.set(Spline::BC::Mix, dy, ddy);
		break;
	}

	case 1: {
		bem.settings.xBC.end.set(Spline::BC::Odd, 0., 0.);
		bem.settings.yBC.end.set(Spline::BC::Even, 0., 0.);
		TaylorCone::circleDerivativeBegin(xy, dx, ddx, dy, ddy);
		bem.settings.xBC.begin.set(Spline::BC::Mix, dx, ddx);
		bem.settings.yBC.begin.set(Spline::BC::Mix, dy, ddy);		
		break;
	}

	case 2: {
		bem.settings.xBC.end.set(Spline::BC::Odd, 0., 0.);
		bem.settings.yBC.end.set(Spline::BC::Even, 0., 0.);
		TaylorCone::circleDerivativeBegin(xy, dx, ddx, dy, ddy);
		bem.settings.xBC.begin.set(Spline::BC::Mix, dx, ddx);
		bem.settings.yBC.begin.set(Spline::BC::Mix, dy, ddy);
		break;
	}

	default:
		break;
	}
	
	bem.initialize(xy);




};

void TaylorCone::setFluidBC(const Bem &bemCone, const Bem &bemPatch, Eigen::VectorXd &fluidBC) const {
	int nCone = bemCone.node().r.rows();
	int nPatch = bemPatch.node().r.rows();
	int nTotal = nCone + nPatch;
	fluidBC.setZero(nTotal);
	//lhs.setZero(nTotal);

	for (int k = 0; k < nTotal; k++) {
		if (k < nCone) {
			double r = bemCone.node().r(k, 0), dr = bemCone.node().r(k, 1);// , ddr = tc.bem0.node().r(k, 2);
			double z = bemCone.node().z(k, 0), dz = bemCone.node().z(k, 1);// , ddz = tc.bem0.node().z(k, 2);
			double nr = -dz / sqrt(dr * dr + dz * dz);
			double nz = dr / sqrt(dr * dr + dz * dz);
			fluidBC(k) = -2. / 3. * (nr * r + nz * z);
		}
		else {
			int kk = k - nCone;
			double r = bemPatch.node().r(kk, 0), dr = bemPatch.node().r(kk, 1);
			double z = bemPatch.node().z(kk, 0), dz = bemPatch.node().z(kk, 1);
			//double nr = -dz / sqrt(dr * dr + dz * dz);
			//double nz = dr / sqrt(dr * dr + dz * dz);
			fluidBC(k) = velocityPotentialFarField(r, z, a);
		}
	}
};

void TaylorCone::setVacuumBC(const Bem &bemCone, const Bem &bemPatch, Eigen::VectorXd &vacuumBC) const {
	int nCone = bemCone.node().r.rows();
	int nPatch = bemPatch.node().r.rows();
	int nTotal = nCone + nPatch;
	vacuumBC.setZero(nTotal);
	for (int k = 0; k < nTotal; k++) {
		if (k < nCone) {			
			vacuumBC(k) = 0.0;
		}
		else {
			int kk = k - nCone;
			double r = bemPatch.node().r(kk, 0), dr = bemPatch.node().r(kk, 1);
			double z = bemPatch.node().z(kk, 0), dz = bemPatch.node().z(kk, 1);			
			vacuumBC(k) = electricPotentialFarField(r, z, b);
		}
	}
};

void TaylorCone::SD2LR(const Eigen::MatrixXd &S, const Eigen::MatrixXd &D, int nSwap, Eigen::MatrixXd &L, Eigen::MatrixXd &R) {

	
	int m = S.rows() - nSwap;

	L.topLeftCorner(nSwap, nSwap) = D.topLeftCorner(nSwap, nSwap);
	L.topRightCorner(nSwap, m) = -S.topRightCorner(nSwap, m);
	L.bottomLeftCorner(m, nSwap) = D.bottomLeftCorner(m, nSwap);
	L.bottomRightCorner(m, m) = -S.bottomRightCorner(m, m);

	R.topLeftCorner(nSwap, nSwap) = S.topLeftCorner(nSwap, nSwap);
	R.topRightCorner(nSwap, m) = -D.topRightCorner(nSwap, m);
	R.bottomLeftCorner(m, nSwap) = S.bottomLeftCorner(m, nSwap);
	R.bottomRightCorner(m, m) = -D.bottomRightCorner(m, m);

}

void TaylorCone::computeResidue(const Bem &bemCone, const Eigen::VectorXd &phi, const Eigen::VectorXd &psin
	, Eigen::VectorXd &residue, Eigen::VectorXd &coord) {
	const int nCone = bemCone.node().r.rows();
	const int o = bemCone.settings.order();
	residue.setZero(nCone);
	coord.setZero(nCone);
	Eigen::VectorXd h; 
	h.setZero(nCone - 1);
	
	for (int k = 0; k < h.size() ; k++) {		
		if (k % o == 0) {
			h(k) = bemCone.e()[k / o].arc() * 1. / o;
		}
		else {
			h(k) = bemCone.e()[k / o].arc() * 1. / o;
		}		
	}	
	for (int k = 0; k < nCone; k++) {
		double r = bemCone.node().r(k, 0), dr = bemCone.node().r(k, 1), ddr = bemCone.node().r(k, 2);
		double z = bemCone.node().z(k, 0), dz = bemCone.node().z(k, 1), ddz = bemCone.node().z(k, 2);
		double nr = -dz / sqrt(dr * dr + dz * dz); // equivalent to -sz
		double nz = dr / sqrt(dr * dr + dz * dz); // equivalent to sr
		double xn = r * nr + z * nz;
		double xs = r * nz + z * (-nr);
		double phis = 0.;

		if (k != 0) {
			if (k == nCone - 1) {								
				double hn3 = h(k-3);
				double hn2 = h(k-2);				
				double hn1 = h(k-1);
				double yn3 = phi(k - 3);
				double yn2 = phi(k - 2);				
				double yn1 = phi(k - 1);				
				double y0 = phi(k);		
				phis = fD3(1, hn3,hn2,hn1, yn3, yn2, yn1, y0,3);
				//phis = fD2(1,  hn2, hn1, yn2, yn1, y0, 2);

			} 
			else if (k == nCone - 2) {
				double hn3 = h(k - 3);
				double hn2 = h(k - 2);
				double hn1 = h(k - 1);
				double h0 = h(k );
				double yn3 = phi(k - 3);
				double yn2 = phi(k - 2);
				double yn1 = phi(k - 1);
				double y0 = phi(k);
				double y1 = phi(k+1);
				phis = fD3(1, hn3, hn2, hn1, yn3, yn2, yn1, y0, 3);
				//phis = fD2(1, hn1, h0, yn1, y0, y1, 2);
			}
			else if (k == nCone - 3) {
				double hn3 = h(k - 3);
				double hn2 = h(k - 2);
				double hn1 = h(k - 1);
				double h0 = h(k);
				double yn3 = phi(k - 3);
				double yn2 = phi(k - 2);
				double yn1 = phi(k - 1);
				double y0 = phi(k);
				double y1 = phi(k + 1);
				phis = fD3(1, hn3, hn2, hn1, yn3, yn2, yn1, y0, 3);
				//phis = fD2(1, hn1, h0, yn1, y0, y1, 2);
			}
			else {
				if (k % o == 0) {					
					double h0 = h(k);
					double h1 = h(k + 1);								
					double h2 = h(k + 2);
					double y0 = phi(k);
					double y1 = phi(k + 1);					
					double y2 = phi(k + 2);					
					double y3 = phi(k + 3);
					//phis = fD2(1, h0 , h1, y0, y1, y2, 0);
					phis = fD3(1, h0, h1, h2, y0, y1, y2, y3, 0);
				}
				else {	
					double h0 = h(k);
					double h1 = h(k + 1);
					double h2 = h(k + 2);
					double y0 = phi(k);
					double y1 = phi(k + 1);
					double y2 = phi(k + 2);
					double y3 = phi(k + 3);
					//phis = fD2(1, h0 , h1, y0, y1, y2, 0);
					phis = fD3(1, h0, h1, h2,y0, y1, y2, y3, 0);
				}
			}

		}
		

		residue(k) += -curv(r, z, dr, dz, ddr, ddz) -0.5 * psin(k) * psin(k);		
		residue(k) += -2. / 9. * xn * xn;
		residue(k) += -1. / 3 * phi(k);
		residue(k) += 0.5 * phis * phis + 2. *xs * phis / 3.;
		/*if (k % o == 1) {
			double rt = bemCone.node().r(k-1, 0), drt = bemCone.node().r(k-1, 1), ddrt = bemCone.node().r(k-1, 2);
			double zt = bemCone.node().z(k-1, 0), dzt = bemCone.node().z(k-1, 1), ddzt = bemCone.node().z(k-1, 2);
			
			residue(k) = -curv(rt, zt, drt, dzt, ddrt, ddzt); };*/
		//residue(k) = 2./ 3. * xs * phis - 1./3. * phi(k) ;
		//residue(k) =  phi(k);
		coord(k) = phis;
		
	}

	
};


void TaylorCone::perturbFluid(const Eigen::MatrixX2d &xyBase, const Bem &bemConeBase, const Bem &bemPatch, int iKnotPerturb, double epsilon, Eigen::VectorXd &output
	, Eigen::MatrixXd &SS, Eigen::MatrixXd &DD) {

	
	Eigen::MatrixXd xy0 = xyBase;
	
	int iNodePerturb = bemConeBase.settings.order() * iKnotPerturb;
	double dr = bemConeBase.node().r(iNodePerturb, 1); 
	double dz = bemConeBase.node().z(iNodePerturb, 1);
	//double perturbr = -dz / sqrt(dr * dr + dz * dz);
	//double perturbz = dr / sqrt(dr * dr + dz * dz);
	double perturbr = 0.0;
	double perturbz = 1.0;
	xy0(iKnotPerturb, 0) +=  perturbr* epsilon;
	xy0(iKnotPerturb, 1) +=  perturbz* epsilon;		
	
	Bem bemConePerturb;
	prepareBem(0, xy0, 0, bemConePerturb);
	int n0 = bemConePerturb.node().r.rows();
	int n1 = bemPatch.node().r.rows();	
	int nTotal = n0 + n1;

	Eigen::MatrixXd S, D,L, R;
	
	if (epsilon < 1e-10) {
		S.setZero(nTotal, nTotal);	D.setZero(nTotal, nTotal);
		R.setZero(nTotal, nTotal);	L.setZero(nTotal, nTotal);
		double distance = 100000000.;
		Bem::assembly(bemConePerturb, bemConePerturb, S, D, distance);
		Bem::assembly(bemConePerturb, bemPatch, S, D, distance);
		Bem::assembly(bemPatch, bemConePerturb, S, D, distance);
		Bem::assembly(bemPatch, bemPatch, S, D, distance);		
		SS = S;		DD = D;		
	}
	else {
		S = SS;		D = DD;
		R.setZero(nTotal, nTotal);	L.setZero(nTotal, nTotal);
		S.topLeftCorner(n0, n0).setZero();
		D.topLeftCorner(n0, n0).setZero();				
		double distance = 100000000.;
		Bem::assembly(bemConePerturb, bemConePerturb, S, D, distance);
	}
	
	Eigen::VectorXd rhs;
	setFluidBC(bemConePerturb, bemPatch, rhs);
	for (int i = 0; i < D.rows(); i++) { D(i, i) = -(D.row(i).sum() - D(i, i)); }
	D.row(n0 - 1) *= 0.;
	D(n0 - 1, n0 - 1) = 1.0;
	D(n0 - 1, n0) = -1.0;
	S.row(n0 - 1) *= 0.;
	TaylorCone::SD2LR(S, D, n0, L, R);
	Eigen::VectorXd phi = L.fullPivLu().solve(R*rhs);
	Eigen::VectorXd coord;
	//TaylorCone::computeResidue(bemConePerturb, phi, output, coord);
	
	if (epsilon < 1e-10) {
		std::ofstream file("./Output/res.txt");
		//std::cout << "phi" << n0 << "\n"; 
		for (int k = 0; k < n0; k++) { file << coord(k) << '\t' << output(k) << '\n'; }
		file.close();
	}
};

void TaylorCone::perturbVacuum(const Eigen::MatrixX2d &xyBase, const Bem &bemConeBase, const Bem &bemFluidPatch, const Bem &bemVacuumPatch
	, int iKnotPerturb, double epsilon, Eigen::VectorXd &output
	, Eigen::MatrixXd &SSF, Eigen::MatrixXd &DDF, Eigen::MatrixXd &SSV, Eigen::MatrixXd &DDV) {

	Eigen::MatrixXd xy0 = xyBase;
	int iNodePerturb = bemConeBase.settings.order() * iKnotPerturb;
	double perturbr = 0.0;
	double perturbz = 1.0;
	xy0(iKnotPerturb, 0) += perturbr * epsilon;
	xy0(iKnotPerturb, 1) += perturbz * epsilon;
	Bem bemConePerturb;
	prepareBem(0, xy0, 0, bemConePerturb);
	int n0 = bemConePerturb.node().r.rows();
	int n1Fluid = bemFluidPatch.node().r.rows();
	int nTotalFluid = n0 + n1Fluid;
	int n1Vacuum = bemVacuumPatch.node().r.rows();
	int nTotalVacuum = n0 + n1Vacuum;

	Eigen::MatrixXd SF, DF, LF, RF;
	Eigen::MatrixXd SV, DV, LV, RV;

	if (epsilon < 1e-10) {
		double distance = 100000000.;
		SF.setZero(nTotalFluid, nTotalFluid);	DF.setZero(nTotalFluid, nTotalFluid);		
		Bem::assembly(bemConePerturb, bemConePerturb, SF, DF, distance);
		Bem::assembly(bemConePerturb, bemFluidPatch, SF, DF, distance);
		Bem::assembly(bemFluidPatch, bemConePerturb, SF, DF, distance);
		Bem::assembly(bemFluidPatch, bemFluidPatch, SF, DF, distance);
		SSF = SF;		DDF = DF;

		SV.setZero(nTotalVacuum, nTotalVacuum);	DV.setZero(nTotalVacuum, nTotalVacuum);		
		Bem::assembly(bemConePerturb, bemConePerturb, SV, DV, distance);
		Bem::assembly(bemConePerturb, bemVacuumPatch, SV, DV, distance);
		Bem::assembly(bemVacuumPatch, bemConePerturb, SV, DV, distance);
		Bem::assembly(bemVacuumPatch, bemVacuumPatch, SV, DV, distance);
		SSV = SV;		DDV = DV;

	}
	else {
		double distance = 100000000.;
		SF = SSF;		DF = DDF;		
		SF.topLeftCorner(n0, n0).setZero();
		DF.topLeftCorner(n0, n0).setZero();		
		Bem::assembly(bemConePerturb, bemConePerturb, SF, DF, distance);
		SV = SSV;		DV = DDV;
		SV.topLeftCorner(n0, n0).setZero();
		DV.topLeftCorner(n0, n0).setZero();
		Bem::assembly(bemConePerturb, bemConePerturb, SV, DV, distance);
	}

	RF.setZero(nTotalFluid, nTotalFluid);	LF.setZero(nTotalFluid, nTotalFluid);
	RV.setZero(nTotalVacuum, nTotalVacuum);	LV.setZero(nTotalVacuum, nTotalVacuum);
	

	
	Eigen::VectorXd rhsF;
	setFluidBC(bemConePerturb, bemFluidPatch, rhsF);	
	for (int i = 0; i < DF.rows(); i++) { DF(i, i) = -(DF.row(i).sum() - DF(i, i)); }
	DF.row(n0 - 1) *= 0.;
	DF(n0 - 1, n0 - 1) = 1.0;
	DF(n0 - 1, n0) = -1.0;
	SF.row(n0 - 1) *= 0.;
	TaylorCone::SD2LR(SF, DF, n0, LF, RF);
	
	Eigen::VectorXd lhsV;
	setVacuumBC(bemConePerturb, bemVacuumPatch, lhsV);
	DV = DV * -1;
	for (int i = 0; i < DV.rows(); i++) { DV(i, i) = -(DV.row(i).sum() - DV(i, i)); }
	DV.row(n0 - 1) *= 0.;
	DV(n0 - 1, n0 - 1) = 1.0;
	DV(n0 - 1, n0) = -1.0;
	SV.row(n0 - 1) *= 0.;
	TaylorCone::SD2LR(SV, DV, n0, LV, RV);


	Eigen::VectorXd lhsF = LF.partialPivLu().solve(RF*rhsF);
	Eigen::VectorXd rhsV = RV.partialPivLu().solve(LV*lhsV);
	Eigen::VectorXd phi = rhsF; 
	phi.head(n0) = lhsF.head(n0);
	Eigen::VectorXd phin = lhsF; 
	phin.head(n0) = rhsF.head(n0);
	Eigen::VectorXd psi = rhsV; 
	psi.head(n0) = lhsV.head(n0);
	Eigen::VectorXd psin = lhsV; 
	psin.head(n0) = rhsV.head(n0);

	
	Eigen::VectorXd coord;
	TaylorCone::computeResidue(bemConePerturb, phi, psin, output, coord);
	
	if (epsilon < 1e-10) {
		std::ofstream file("./Output/phi.txt");
		
		Eigen::IOFormat fmt(Eigen::FullPrecision, 0, "\t", "\n", "", "", "", "");
		file << phi.format(fmt);
		file.close();
		file.open("./Output/psin.txt");		
		file << psin.format(fmt);;
		file.close();
		file.open("./Output/r.txt");
		file << bemConePerturb.node().r.format(fmt);
		file.close();
		file.open("./Output/z.txt");
		file << bemConePerturb.node().z.format(fmt);
		file.close();
		

		//std::cout << "psin\t" << rhsV(0)/2. << "\n";
		//for (int k = 0; k < n0; k++) { file << coord(k) << '\t' << rhsV(k) << '\n'; }
		//file.close();

		scan(bemConePerturb,bemFluidPatch,phin,phi,coord,-6,0.0, SCAN_NZ, SCAN_NR);

	}
};

void TaylorCone::scan(const Bem &bemCone, const Bem &bemPatch, 
	const Eigen::VectorXd &q, const Eigen::VectorXd &p, const Eigen::VectorXd &ps,
	double zLowerBound, double rUpperBound, int zGrid, int rGrid) {

	std::ofstream file("./Output/field.txt");	
	Eigen::IOFormat fmt(Eigen::FullPrecision, 0, "\t", "\n", "", "", "", "");
	Eigen::MatrixXd output;
	output.setZero(zGrid * rGrid, 5);

#pragma omp parallel for 
	for (int j = 0; j < zGrid * rGrid; j++) {
		int ir = j / zGrid;
		int iz = j % zGrid;
		double rp = bemCone.node().r(ir, 0);
		double zpUpper = bemCone.node().z(ir, 0);

		double drUpper = bemCone.node().r(ir, 1);
		double dzUpper = bemCone.node().z(ir, 1);
		double nrUpper = -dzUpper / sqrt(drUpper * drUpper + dzUpper * dzUpper);
		double nzUpper = drUpper / sqrt(drUpper * drUpper + dzUpper * dzUpper);
		double srUpper = nzUpper;
		double szUpper = -nrUpper;

		double zp = 0;
		double fieldValue, fieldValueDr, fieldValueDz;

		if (iz == 0) {
			fieldValue = p(ir);			
			zp = zpUpper;
			fieldValueDr = ps(ir) * srUpper + q(ir) * nrUpper;
			fieldValueDz = ps(ir) * szUpper + q(ir) * nzUpper;
		}
		else {
			double t = (double)iz / (double)(zGrid - 1.);
			t = gridDistribution(t, 1.1);
			zp = zpUpper + (zLowerBound)* t;
			fieldValue = Bem::assembly(rp, zp, bemCone, q, p);
			fieldValue += Bem::assembly(rp, zp, bemPatch, q, p);

			if (rp < 1e-13) { fieldValueDr = 0; }
			else {
				fieldValueDr = Bem::assemblyDr(rp, zp, bemCone, q, p);
				fieldValueDr += Bem::assemblyDr(rp, zp, bemPatch, q, p);
			}

			fieldValueDz = Bem::assemblyDz(rp, zp, bemCone, q, p);
			fieldValueDz += Bem::assemblyDz(rp, zp, bemPatch, q, p);

		}

		output(j, 0) = rp;
		output(j, 1) = zp;
		output(j, 2) = fieldValue;
		output(j, 3) = fieldValueDr;
		output(j, 4) = fieldValueDz;



	}

	file << output.format(fmt);
	//for (int ir = 0; ir < rGrid; ir++) {
	//	double rp = bemCone.node().r(ir, 0);
	//	double zpUpper = bemCone.node().z(ir, 0);		
	//	double fieldValue;
	//	fieldValue = p(ir);
	//	file << rp << '\t' << zpUpper  << '\t' << fieldValue << '\n';

	//	for (int iz = 1; iz < zGrid; iz++) {
	//		double t = (double)iz / (double)(zGrid - 1.);
	//		t = gridDistribution(t, 1.1);
	//		//double zp = zpUpper+ (zLowerBound - zpUpper) * t ;			
	//		double zp = zpUpper + (zLowerBound) * t;
	//		fieldValue = Bem::assembly(rp, zp, bemCone,q, p);
	//		fieldValue += Bem::assembly(rp, zp, bemPatch, q, p);
	//		file << rp << '\t' << zp << '\t' << fieldValue << '\n';
	//	}

	//}
	file.close();

	


};