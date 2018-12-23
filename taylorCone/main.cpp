#include "stdafx.h"
#include "numeric.h"
#include "bem.h"
#include "taylorCone.h"
#include <omp.h>
#include <Eigen/Dense>
#include <string>
//#define COMPUTE
//#define TEST1
//#define TEST2
double inputC1 = -0.87318;
double inputB0 = 1.7576;
//double inputC1 =0.0;
//double inputB0 = 1.15;
int nConeKnots = 200 + 1;
double rTruncate = 50.;
void parse(int &argc, char ** &argv) {
	switch (argc)	{	
	case 3: {
		inputC1 = atof(argv[1]);
		inputB0 = atof(argv[2]);
		break;
	}
	case 4: {
		inputC1 = atof(argv[1]);
		inputB0 = atof(argv[2]);
		nConeKnots = atoi(argv[3]);
		break;
	}
	case 5: {
		inputC1 = atof(argv[1]);
		inputB0 = atof(argv[2]);
		nConeKnots = atoi(argv[3]);
		rTruncate = atof(argv[4]);
		break;
	}
	default:
		break;
	}

}

Eigen::IOFormat fmt(Eigen::FullPrecision, 0, "\t", "\n", "", "", "", "");
void printSpline(const Spline &sp, const std::string &name) {
	std::ofstream file(name);	
	file << sp.x().format(fmt) << '\n'	<< sp.y().format(fmt) << '\n' << sp.h().format(fmt) << '\n';
	file.close();
}

int main(int argc, char** argv) {
	parse(argc, argv);
	TaylorCone tc(inputC1, inputB0);
	printf("--------------------a b c--------------------\n");
	printf("a: %+18.15f\t%+18.15f\t%+18.15f\t%+18.15f\t%+18.15f\t\n", tc.a[0], tc.a[1], tc.a[2], tc.a[3], tc.a[4]);
	printf("b: %+18.15f\t%+18.15f\t%+18.15f\t%+18.15f\t%+18.15f\t\n", tc.b[0], tc.b[1], tc.b[2], tc.b[3], tc.b[4]);
	printf("c: %+18.15f\t%+18.15f\t%+18.15f\t%+18.15f\t%+18.15f\t\n", tc.c[0], tc.c[1], tc.c[2], tc.c[3], tc.c[4]);
	printf("--------------------a b c--------------------\n");	
	
	int nElectricKnots = (int)floor((nConeKnots - 1) * 0.9) + 1;
	int nVelocityKnots = (int)floor((nConeKnots - 1) * 2.3) + 1;
	tc._xy0 = TaylorCone::generateCone(4, rTruncate,tc.c, nConeKnots,TaylorCone::gridDistribution);
	tc._xy1 = TaylorCone::generateCircle(tc._xy0(tc._xy0.rows() - 1, 0), tc._xy0(tc._xy0.rows() - 1, 1), 1, (int)(nElectricKnots));
	tc._xy2 = TaylorCone::generateCircle(tc._xy0(tc._xy0.rows() - 1, 0), tc._xy0(tc._xy0.rows() - 1, 1), 0, (int)(nVelocityKnots));

	tc.prepareBem(0, tc._xy0, 0, tc.bem0);
	int n0 = tc.bem0.node().r.rows();
	printSpline(tc.bem0.sp(), "./Output/sp0.txt");


	tc.prepareBem(1, tc._xy1, n0, tc.bem1);
	int n1 = tc.bem1.node().r.rows();
	printSpline(tc.bem1.sp(), "./Output/sp1.txt");

	tc.prepareBem(2, tc._xy2, n0, tc.bem2);
	int n2 = tc.bem2.node().r.rows();
	printSpline(tc.bem2.sp(), "./Output/sp2.txt");	
	
	//int nTotal = n0 + n1;
	
	Eigen::VectorXd res, resPerturb,tmp;
	double epsilon = 0.0001;
	Eigen::MatrixXd J(n0, nConeKnots - 1);
	//Eigen::MatrixXd J(nConeKnots, nConeKnots - 1);

	
	Eigen::MatrixXd SF0, DF0,SV0, DV0;
	tc.perturbVacuum(tc._xy0, tc.bem0, tc.bem1, tc.bem2, 0, 0., res, SF0, DF0, SV0, DV0);
	Eigen::VectorXd abc(nConeKnots) ,abcP(nConeKnots);
	abc = Eigen::Map<Eigen::VectorXd, 0, Eigen::InnerStride<2> >(res.data(), nConeKnots);
	
	
	
	printf("error = %2.2e\n", res.norm());
	for (int j = 0; j < J.cols(); j++) {		
		printf("\r%03d", j);
		fflush(stdout);
		tc.perturbVacuum(tc._xy0, tc.bem0, tc.bem1, tc.bem2, j, epsilon, resPerturb, SF0, DF0, SV0, DV0);
		Eigen::VectorXd fk = (resPerturb - res) / epsilon;
		abcP = Eigen::Map<Eigen::VectorXd, 0, Eigen::InnerStride<2> >(fk.data(), nConeKnots);
		//J.col(j) = abcP;		
		J.col(j) = fk;		
	}
	printf("\n-----------------\n");
		
	
	int counter = 1;

	while (counter < 33)	{
		tmp = J.colPivHouseholderQr().solve(-res);
		//tmp = J.fullPivLu().solve(-abc);
		for (int k = 0; k < tmp.size(); k++) {
			int iNodePerturb = tc.bem0.settings.order() * k;
			//double dr = tc.bem0.node().r(iNodePerturb, 1);
			//double dz = tc.bem0.node().z(iNodePerturb, 1);
			double perturbr = 0.0;
			double perturbz = 1.0;
			//double perturbr = -dz / sqrt(dr * dr + dz * dz);
			//double perturbz = dr / sqrt(dr * dr + dz * dz);
			tc._xy0(k, 0) += perturbr * tmp(k) * 0.99;
			tc._xy0(k, 1) += perturbz * tmp(k) * 0.99;
		}
		tc.prepareBem(0, tc._xy0, 0, tc.bem0);	
		printSpline(tc.bem0.sp(), "./Output/sp0.txt");
	/*	std::ofstream file("./Output/answer0.txt");
		for (int k = 0; k < tmp.size(); k++) {
			file << tc._xy0(k, 0) << '\t' << tc._xy0(k, 1) << '\n';
		}
		file.close();*/
		res.setZero();

		//tc.perturbFluid(tc._xy0, tc.bem0, tc.bem1, 0, 0., res,FS0,FD0);
		tc.perturbVacuum(tc._xy0, tc.bem0, tc.bem1, tc.bem2, 0, 0., res, SF0, DF0, SV0, DV0);
		abc = Eigen::Map<Eigen::VectorXd, 0, Eigen::InnerStride<2> >(res.data(), nConeKnots);


		printf("error\t %2.2e\n", res.cwiseAbs().maxCoeff());
		
		if (counter % 4 == 0) {
			J.setZero();
			//tc.perturbFluid(tc._xy0, tc.bem0, tc.bem1, 0, 0., res);
			for (int j = 0; j < J.cols(); j++) {
				printf("\r%03d", j);
				fflush(stdout);
				tc.perturbVacuum(tc._xy0, tc.bem0, tc.bem1, tc.bem2, j, epsilon, resPerturb, SF0, DF0, SV0, DV0);
				Eigen::VectorXd fk = (resPerturb - res) / epsilon;
				abcP = Eigen::Map<Eigen::VectorXd, 0, Eigen::InnerStride<2> >(fk.data(), nConeKnots);
				//J.col(j) = abcP;
				J.col(j) = fk;
			}
		}
		
		counter = counter + 1;
		//std::cout << res << "\n";
		std::cout << 2 * tc.bem0.sp().y()(0, 2)/ pow(tc.bem0.sp().x()(0, 1),2.0) <<"\n";
	}
	
















	// ------------------------------------	
	
	//Eigen::MatrixXd S, D, L, R;		
	//S.setZero(nTotal, nTotal);	D.setZero(nTotal, nTotal);
	//R.setZero(nTotal, nTotal);	L.setZero(nTotal, nTotal);	
	//Bem::assembly(tc.bem0, tc.bem0, S, D);	Bem::assembly(tc.bem0, tc.bem1, S, D);
	//Bem::assembly(tc.bem1, tc.bem0, S, D);	Bem::assembly(tc.bem1, tc.bem1, S, D);	
	//Eigen::VectorXd rhs, lhs;
	//tc.setFluidBC(tc.bem0, tc.bem1, rhs);
	//for (int i = 0; i < D.rows(); i++) {	D(i, i) = -(D.row(i).sum() - D(i, i));}
	//D.row(n0 - 1) *= 0.;
	//D(n0 - 1, n0 - 1) = 1.0;
	//D(n0 - 1, n0) = -1.0;
	//S.row(n0 - 1) *= 0.;	
	//TaylorCone::SD2LR(S, D, n0, L, R);
	//Eigen::VectorXd phi = L.fullPivLu().solve(R*rhs);	
	//Eigen::VectorXd residue, coord;
	//TaylorCone::computeResidue(tc.bem0, phi, residue, coord);

	



	//std::ofstream file("./Output/answer0.txt");
	//for (int k = 0; k < n0; k++) {	file << coord(k) <<'\t' << residue(k) << '\n';	}
	//file.close();

	//file.open("./Output/xy0.txt");
	//file << tc._xy0 << '\n';
	//file.close();
	//file.open("./Output/xy1.txt");
	//file << tc._xy1 << '\n';
	//file.close();
	//file.open("./Output/xy2.txt");
	//file << tc._xy2 << '\n';
	//file.close();


#ifdef TEST1
	int elementOrder = 2;
	for (int global = 7; global < 21; global++) {
		int n = (int)pow(sqrt(2.0), global) + 1;		
		Eigen::MatrixX2d xy0 = circle(0.0, M_PI/2, n);
		Bem bem0;
		bem0.settings.indexShift(0);
		bem0.settings.order(elementOrder);
		bem0.settings.qdOrder(20);
		bem0.settings.xBC.end.set(Spline::BC::Even,0,0);		
		bem0.settings.yBC.end.set(Spline::BC::Odd,0,0);
		bem0.initialize(xy0);
		printSpline(bem0.sp(),"./Output/sp0.txt");
		int n0 = bem0.node().r.rows();
		

		Eigen::MatrixX2d xy1 = line(xy0(xy0.rows() - 1, 0 ), xy0(xy0.rows() - 1, 1), 0.0, 0.0,  n );
		Bem bem1;
		bem1.settings.indexShift(bem0.node().r.rows() );
		bem1.settings.order(elementOrder);
		bem1.settings.qdOrder(20);
		bem1.settings.xBC.end.set(Spline::BC::Odd, 0, 0);		
		bem1.settings.yBC.end.set(Spline::BC::Even, 0, 0);
		bem1.initialize(xy1);
		printSpline(bem1.sp(), "./Output/sp1.txt");
		int n1 = bem1.node().r.rows();

		Eigen::MatrixXd S, D, B, L,R;
		int nTotal = n0  + n1;
		S.setZero(nTotal, nTotal);
		B.setZero(nTotal, nTotal);
		D.setZero(nTotal, nTotal);
		R.setZero(nTotal, nTotal);
		L.setZero(nTotal, nTotal);
		
		B.setIdentity();
		B *= 0.5;
		B(n0 - 1, n0 - 1) = 0.25;
		B(n0, n0) = 0.25;
		
		Bem::assembly(bem0, bem0, S, D);								
		Bem::assembly(bem0, bem1, S, D);						
		Bem::assembly(bem1, bem0, S, D);		
	    Bem::assembly(bem1, bem1, S, D);		
		
		D = D + B;				
				
		Eigen::VectorXd rhs, lhs, errorCurvature;
		rhs.setZero(nTotal);
		lhs.setZero(nTotal);
		errorCurvature.setZero(n0);
		
		for (int k = 0; k < nTotal; k++) {
			if (k < n0) {
				double r = bem0.node().r(k, 0), dr = bem0.node().r(k, 1), ddr = bem0.node().r(k, 2);
				double z = bem0.node().z(k, 0), dz = bem0.node().z(k, 1), ddz = bem0.node().z(k, 2);				
				double nr = -dz / sqrt(dr * dr + dz * dz);
				double nz = dr / sqrt(dr * dr + dz * dz);				
				rhs(k) = phin(nr,nz,r,z);
				lhs(k) = phi(r, z);
				errorCurvature(k) = curvAnalytic(acos(z / sqrt(r *r + z * z)))-  curv(r, z, dr, dz, ddr, ddz);					
			}
			else {
				int kk = k - n0;
				double r = bem1.node().r(kk, 0), dr = bem1.node().r(kk, 1);
				double z = bem1.node().z(kk, 0), dz = bem1.node().z(kk, 1);
				double nr = -dz / sqrt(dr * dr + dz * dz);
				double nz = dr / sqrt(dr * dr + dz * dz);
				rhs(k) = phi(r, z);
				lhs(k) = phin(nr, nz, r, z);
			}
		}		
		D.row(n0 - 1) *= 0.;		
		D(n0 -1 , n0 - 1) = 1.0;
		D(n0 -1,  n0 ) = -1.0;
		S.row(n0 - 1) *= 0.;		

		L.topLeftCorner(n0, n0) = D.topLeftCorner(n0, n0);
		L.topRightCorner(n0, n1) = -S.topRightCorner(n0, n1);
		L.bottomLeftCorner(n1, n0) = D.bottomLeftCorner(n1, n0);
		L.bottomRightCorner(n1, n1) = -S.bottomRightCorner(n1, n1);

		R.topLeftCorner(n0, n0) = S.topLeftCorner(n0, n0);
		R.topRightCorner(n0, n1) = -D.topRightCorner(n0, n1);
		R.bottomLeftCorner(n1, n0) = S.bottomLeftCorner(n1, n0);
		R.bottomRightCorner(n1, n1) = -D.bottomRightCorner(n1, n1);

		Eigen::VectorXd answer =  R.fullPivLu().solve(L*lhs);
		Eigen::VectorXd errorDirchlet(n0), errorNeumann(n1);
		for (int kk = 0; kk < n0; kk++) { errorDirchlet(kk) = answer(kk) - rhs(kk); }		
		for (int kk = 0; kk < n1; kk++) { errorNeumann(kk) = answer(n0 + kk) - rhs(n0 + kk); }
		printf("{%04d, %16.16f, %16.16f, %16.16f},\n", 
			nTotal,
			errorDirchlet.cwiseAbs().maxCoeff(), 
			errorNeumann.cwiseAbs().maxCoeff(), 
			errorCurvature.cwiseAbs().maxCoeff()
		);
}
#endif // TEST2
#ifdef TEST1

		Eigen::MatrixX2d xy0 = circle(M_PI, n);
		Bem bem0;
		bem0.settings.indexShift(0);
		bem0.settings.order(2);
		bem0.settings.qdOrder(20);		
		bem0.initialize(xy0);
		printSpline(bem0.sp(), "./Output/sp0.txt");

		Eigen::MatrixXd S, D, B;
		int nTotal = bem0.node().r.rows();
		std::cout << nTotal;
		S.setZero(nTotal, nTotal);
		B.setZero(nTotal, nTotal);
		D.setZero(nTotal, nTotal);
		Bem::assembly(bem0, bem0, S, D);

		for (int i = 0; i < D.rows(); i++) {
			D(i, i) = -(D.row(i).sum() - D(i, i));
		}
			   
		B.resize(D.rows(), D.cols());
		B.setIdentity();
		B *= 0.5;

		const Eigen::VectorXd &z = bem0.node().z.col(0);
		const Eigen::VectorXd &r = bem0.node().r.col(0);
		const Eigen::VectorXd &dr = bem0.node().r.col(1);
		const Eigen::VectorXd &dz = bem0.node().z.col(1);
		const Eigen::VectorXd &cosTh = z.array() *(r.array().square() + z.array().square()).rsqrt();
		const Eigen::VectorXd &R2 = (r.array().square() + z.array().square());

		Eigen::VectorXd nr = -dz.array() *(dr.array().square() + dz.array().square()).rsqrt();
		Eigen::VectorXd nz = dr.array() *(dr.array().square() + dz.array().square()).rsqrt();
		Eigen::VectorXd phi = R2.array() * (0.5 * (3. * cosTh.array().square() - 1.0));

		const Eigen::VectorXd lhs = D * phi;
		
		printf("{%04d, %16.16f},\n", 2 * (n - 1) + 1, (S.fullPivLu().solve(lhs) -
			Eigen::VectorXd(nr.array() * (-r.array()) + nz.array() * (2. * z.array()))
			).cwiseAbs().maxCoeff());

#endif // TEST1
			   		 	  	  
#ifdef COMPUTE
		Eigen::MatrixXd S, D, B;	
		Bem::assembly(bem0, bem0,S, D);
		//std::cout << (D).row(D.rows() - 2).sum()+0.5 << std::endl;
		//std::cout << (D).rowwise().sum().array() + 0.5  << std::endl;

		for (int i = 0; i < D.rows(); i++){

			D(i, i) = -(D.row(i).sum() - D(i, i));
		}

		
		
		B.resize(D.rows(), D.cols());
		B.setIdentity();
		B *= 0.5;
			
		const Eigen::VectorXd &z = bem.node().z.col(0);
		const Eigen::VectorXd &r = bem.node().r.col(0);
		const Eigen::VectorXd &dr = bem.node().r.col(1);
		const Eigen::VectorXd &dz = bem.node().z.col(1);
		const Eigen::VectorXd &cosTh = z.array() *(r.array().square() + z.array().square()).rsqrt();
		const Eigen::VectorXd &R2 = (r.array().square() + z.array().square());

		
		Eigen::VectorXd nr = -dz.array() *(dr.array().square() + dz.array().square()).rsqrt();
		Eigen::VectorXd nz = dr.array() *(dr.array().square() + dz.array().square()).rsqrt();
		Eigen::VectorXd phi = R2.array() * (0.5 * (3. * cosTh.array().square() - 1.0));

		const Eigen::VectorXd lhs = D * phi;

		//std::cout << lhs;
		printf("{%04d, %16.16f},\n", 2 * (n -1) + 1, (S.fullPivLu().solve(lhs) -		
			  Eigen::VectorXd (nr.array() * (-r.array())  + nz.array() * (2. * z.array()))
			).cwiseAbs().maxCoeff());
		
		//const Eigen::VectorXd lhs = (S) * (z);
		//std::cout << ((B - D).fullPivLu().solve(lhs)) ;
		//printf("{%04d, %16.16f},\n", n - 1, ((B - D).fullPivLu().solve(lhs) -0.5 * z).cwiseAbs().maxCoeff());

#endif // COMPUTE
		
		

	return 0;
}


