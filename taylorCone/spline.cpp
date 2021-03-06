#include "stdafx.h"
#include "spline.h"
#include <Eigen/Sparse>
#include "numeric.h"

Spline& Spline::operator = (const Spline &sp) {
	if (this == &sp)
		return *this;
	else {
		_node.resize(sp.node().rows(), 2);
		_node = sp.node();
		_h.resize(sp.h().size());
		_h = sp.h();
		_x.resize(sp.x().rows(), 6);
		_x = sp.x();
		_y.resize(sp.y().rows(), 6);
		_y = sp.y();
		return *this;
	}
};

void Spline::node(const Eigen::MatrixX2d &xy) {	
	_node.resize(xy.rows(), xy.cols());
	_node = xy;
	h(_node);
	setComponent(0, _x);
	setComponent(1, _y);
}

void Spline::h(const Eigen::MatrixX2d &node) {
	_h.setZero(node.rows() - 1);
	for (int i = 0; i < _h.size(); i++) {
		_h(i) = sqrt(pow(node(i + 1, 0) - node(i, 0), 2.0) + pow(node(i + 1, 1) - node(i, 1), 2.0));
	}
};

void Spline::setComponent(int i, Coef &x) {
	x.resize(_node.rows(), 6);
	x.setZero();
	x.col(0) = _node.col(i);
};

void Spline::x(BC bc0, BC bc1, double a0, double b0, double a1, double b1) {
	computeCoef(_x, bc0, bc1, a0, b0, a1, b1);
};

void Spline::y(BC bc0, BC bc1, double a0, double b0, double a1, double b1) {
	computeCoef(_y, bc0, bc1, a0, b0, a1, b1);
};

void Spline::computeCoef(Coef &x, BC bc0, BC bc1, double a0, double b0, double a1, double b1) {
	int N = _node.rows() -1;	
	Eigen::VectorXd rhs(2 * N + 2);
	const Eigen::VectorXd X = x.col(0);
	rhs.setZero();
	Eigen::SparseMatrix<double> A(2 * N + 2, 2 * N + 2);
	A.reserve(Eigen::VectorXi::Constant(A.cols(), 6));

	switch (bc0) {
	case BC::Even: {//printf("Even at 0\t");
		rhs(0) = 0;
		rhs(1) = -10.* X(0) + 10. * X(1);	
		double h0 = _h(0);
		A.insert(0, 0) = 1.0;
		A.insert(1, 0) = 6. * h0;
		A.insert(1, 1) = 3. * h0 * h0;
		A.insert(1, 2) = 4. * h0;
		A.insert(1, 3) = - h0 * h0;
		break;
	}	
	case BC::Odd: {//printf("Odd at 0\t");
		rhs(0) = 15.* X(0) - 15. * X(1);
		rhs(1) = 0;
		double h0 = _h(0);
		A.insert(0, 0) = - 8. * h0;
		A.insert(0, 1) = - 3. * h0 * h0;
		A.insert(0, 2) = - 7. * h0;
		A.insert(0, 3) = 2. * h0 * h0;
		A.insert(1, 1) = 1.0;
		break;
	}	
	case BC::Mix: {
		rhs(0) = a0;
		rhs(1) = b0/2;
		A.insert(0, 0) = 1.0;
		A.insert(1, 1) = 1.0;
		break;		
	}
	default:
		break;
	}
	switch (bc1) {
	case BC::Even: {		//printf("even at 1\n");
		rhs(rhs.size() - 2) = 0;
		rhs(rhs.size() - 1) = 10.* X(X.size() - 2) - 10. * X(X.size() - 1);
		double hm1 = _h(_h.size()-1);
		A.insert(A.rows() - 2, A.rows() - 2) = 1.0;
		
		A.insert(A.rows() - 1, A.rows() - 4) = -4. * hm1;
		A.insert(A.rows() - 1, A.rows() - 3) = -1. * hm1 * hm1;
		A.insert(A.rows() - 1, A.rows() - 2) = -6. * hm1;
		A.insert(A.rows() - 1, A.rows() - 1) = +3. * hm1 * hm1;
		break;
	}
	case BC::Odd: {//printf("odd at 1\n");
		rhs(rhs.size() - 2) = 15.* X(X.size() - 2) - 15. * X(X.size() - 1);
		rhs(rhs.size() - 1) = 0;
		double hm1 = _h(_h.size() - 1);
		A.insert(A.rows() - 2, A.rows() - 4) = -7. * hm1;
		A.insert(A.rows() - 2, A.rows() - 3) = -2. * hm1 * hm1;
		A.insert(A.rows() - 2, A.rows() - 2) = -8. * hm1;
		A.insert(A.rows() - 2, A.rows() - 1) = +3. * hm1 * hm1;
		A.insert(A.rows() - 1, A.rows() - 1) = 1.0;
		break;
	}
	case BC::Mix: {
		rhs(rhs.size() - 2) = a1;
		rhs(rhs.size() - 1) = b1 / 2.;
		A.insert(A.rows() - 2, A.rows() - 2) = 1.0;
		A.insert(A.rows() - 1, A.rows() - 1) = 1.0;
		break;
	}
	default:
		break;
	}
	
	for (int j = 1; j <= N-1; j++) {
		double lb = _h(j) / _h(j - 1);
		double lb2 = lb * lb;
		double lb3 = lb * lb * lb;
		double lb4 = lb * lb * lb * lb;
		double hj = _h(j);
		double hj2 = _h(j) * _h(j);

		rhs(2 * j) = -10. * ( lb3 * X(j - 1) - (1. + lb3) * X(j) + X(j + 1) );
		rhs(2 * j + 1) = 15. * (-lb4 * X(j - 1) + (lb4 - 1.)*X(j) + X(j + 1));

		A.insert(2 * j, 2 * j - 2) =  4.0 * hj * lb2;
		A.insert(2 * j, 2 * j - 1) =  1.0 * hj2 * lb;
		A.insert(2 * j, 2 * j    ) =  6.0 * hj * (lb2 - 1.);
		A.insert(2 * j, 2 * j + 1) = -3.0 * hj2 * (1. + lb);
		A.insert(2 * j, 2 * j + 2) = -4.0 * hj;
		A.insert(2 * j, 2 * j + 3) =  1.0 * hj2;

		A.insert(2 * j + 1, 2 * j - 2) =  7.0 * hj * lb3;
		A.insert(2 * j + 1, 2 * j - 1) =  2.0 * hj2 * lb2;
		A.insert(2 * j + 1, 2 * j    ) =  8.0 * hj * (1. + lb3 );
		A.insert(2 * j + 1, 2 * j + 1) =  3.0 * hj2 * (1. - lb2);
		A.insert(2 * j + 1, 2 * j + 2) =  7.0 * hj;
		A.insert(2 * j + 1, 2 * j + 3) = -2.0 * hj2;

	}
	
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver(A);
	solver.compute(A);	
	Eigen::VectorXd lhs = solver.solve(rhs);

	//std::cout << lhs <<std::endl;

	Eigen::Map<Eigen::VectorXd, 0, Eigen::InnerStride<2> > c1(lhs.data()  , lhs.size() / 2);
	Eigen::Map<Eigen::VectorXd, 0, Eigen::InnerStride<2> > c2(lhs.data()+1, lhs.size() / 2);

	x.col(1) = c1;
	x.col(2) = c2;

	for (int j = 0; j <= N - 1; j++) {
		double hj = _h(j);
		double c1j = c1(j), c2j = c2(j), c1jp1 = c1(j+1), c2jp1 = c2(j+1);
		
		x(j, 1) = c1j * hj;

		x(j, 2) = c2j * hj * hj;

		x(j, 3) = -6. * hj * c1j - 3. * hj * hj * c2j - 4. * hj * c1jp1 + 1. * hj * hj * c2jp1 
				+ 10.* (X(j + 1) - X(j));
		
		x(j, 4) = +8. * hj * c1j + 3. * hj * hj * c2j + 7. * hj * c1jp1 - 2. * hj * hj * c2jp1 
				- 15.* (X(j + 1) - X(j));
		
		x(j, 5) = -3. * hj * c1j - 1. * hj * hj * c2j - 3. * hj * c1jp1 + 1. * hj * hj * c2jp1 
				+ 6. * (X(j + 1) - X(j));
	}	
	
	//Eigen::IOFormat fmt(3, 10, ",", "\n", "[", "]");	
	//std::cout << Eigen::MatrixXd(A).format(fmt) <<std::endl;
	
};

const Eigen::Vector3d Spline::d(const Coef &x, int i, double t) const {
	double t2 = t * t;
	double t3 = t * t * t;
	double t4 = t * t * t * t;
	double t5 = t * t * t * t * t;
	double x0 = x(i, 0);
	double x1 = x(i, 1);
	double x2 = x(i, 2);
	double x3 = x(i, 3);
	double x4 = x(i, 4);
	double x5 = x(i, 5);
	Eigen::Vector3d tmp(0,0,0);
	tmp(0) = x0 + t * x1 + t2 * x2 + t3 * x3 + t4 * x4 + t5 * x5;
	tmp(1) = x1 + 2. * t* x2 + 3. * t2 * x3 + 4. * t3 * x4 + 5. * t4 * x5;
	tmp(2) = 2. * x2 + 6. * t * x3 + 12. * t2 * x4 + 20. * t3 * x5;
	return tmp;
};

double Spline::localArc(int i, double t, int nqd) const {	
	if (i < _node.rows() - 1 && i >= 0) {		
		//Eigen::Map<const Eigen::Matrix2Xd> v(Numeric::qd[nqd],  2, nqd);				
		const double*qdx = Numeric::qd_GL_x[nqd];
		const double*qdw = Numeric::qd_GL_w[nqd];
		double arc = 0;
		for (int k = 0; k < nqd; k++) {
			double ab = t * qdx[k];
			arc += qdw[k] * sqrt(pow((d(_x, i, ab))(1),2.0) + pow((d(_y, i, ab))(1), 2.0));
		}		
		return t * arc;		
	}
	else {
		printf("Spline index out of range! ");
		return 0.;
	}
};

double Spline::arc2t(int i, double arc, double eps, int nqd) const {
	double x0 = 0.5;
	double epsilon = eps;
	double f0 = localArc(i, x0, nqd) - arc;	

	int counter = 0;
	while (std::abs(f0) > epsilon ) {		
		double df0 = sqrt(pow((d(_x, i, x0))(1), 2.0) + pow((d(_y, i, x0))(1), 2.0));
		x0 = x0 - f0 / df0;
		f0 = localArc(i, x0, nqd) - arc;
		counter++;
		if (counter > 10) { printf("arc2t not converging .. try lowering error tolerance"); }
	}
	return x0;
};