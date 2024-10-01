#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#include "bicgstab.h"
#include "ilduc.h"

using namespace std;
const double PI = 3.14159265358979323846;

const double L = 100.0;
const double T = 1.0;
const double A = 1.0; //площадь поперечного сечения
const double q = 1.0;

const double Swr = 0.2, Sor = 0.15;
const double mu_o = 50.0, mu_w = 0.5;
const double No = 2.0, Nw = 2.0;
const double kro_ = 1.0, krw_ = 0.6;

const double phi = 0.1;

const double dt = 0.01;
const double h = 1.0; // 0.1;
const double h2 = h * h;
const double tol = 1e-8;

int OutIter = (int)(0.1 / dt);
//num of nodes (1D - task):
const int N = (int)(L / h) + 1;
const int Nt = (int)(T / dt); //time steps

double Sae(double Sa, double Sar) { //Sar - остаточная нас. фазы a
	return (Sa - Sar) / (1.0 - Sor - Swr);
}

//fw(i) = func(Sw(i))
double fw(double Sw) {
	double alpha = mu_w * kro_ / mu_o / krw_;
	double Soe = Sae(1.0 - Sw, Sor);
	double Swe = Sae(Sw, Swr);
	double c = 1.0 + alpha * std::pow(Soe, No) / std::pow(Swe, Nw);
	return 1.0 / c;
}

double dfwds (double Sw) {
	double s1 = (1.0 - Sw - Sor);
	double s2 = (1.0 - Swr - Sor);
	double swwr = Sw - Swr;
	double alpha = mu_w * kro_ / mu_o / krw_;
	return 2.0 * alpha * s1 * s2 * swwr / std::pow(swwr*swwr + alpha*s1*s1, 2.0);
}
//Jacobi Matrix struct:
//func fi - is our difference finite equasion fr node i
/*
df0/ds0 df0/ds1    0    0......0.....................0....
   0    df1/ds1 df1/ds2 0..... 0.....................0....
   .......................................................
   .......................................................
   0   ...........df(s(i-1))/ds(i-1)  df(s(i))/ds(i) 0 ...

*/
//Difference of equasion: 
//df/ds(i-1)
double dfds1(double Sw) {
	return -q * dfwds(Sw) / A / h;
}
//df/ds(i)
double dfds2(double Sw) {
	return phi/dt + q * dfwds(Sw) / A / h;
}

void generateVtkFile(std::string folder_path, std::vector<double> x, double h, int Nx, int Ny) {
	std::ofstream vtkFile(folder_path);

	vtkFile << "# vtk DataFile Version 3.0" << std::endl;
	vtkFile << "vtk output file" << std::endl;
	vtkFile << "ASCII" << std::endl;
	vtkFile << "DATASET STRUCTURED_GRID" << std::endl;
	vtkFile << "DIMENSIONS " << Nx << " " << Ny << " 1" << std::endl;
	vtkFile << "POINTS " << Nx * Ny << " double" << std::endl;


	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			vtkFile << (double)(i * h) << " " << (double)(j * h) << " 0" << std::endl;
		}
	}

	vtkFile << "POINT_DATA " << Nx * Ny << std::endl;
	vtkFile << "SCALARS Sw double" << std::endl;
	vtkFile << "LOOKUP_TABLE default" << std::endl;

	// Write solution values
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {

			vtkFile << x[i] << std::endl;
		}
	}

	vtkFile << "SCALARS So double" << std::endl;
	vtkFile << "LOOKUP_TABLE default" << std::endl;

	// Write solution values
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {

			vtkFile << 1.0 - x[i] << std::endl;
		}
	}


	vtkFile.close();
}

void generateCsvFile(std::string folder_path, std::vector<double> x, double h, int Nx, int Ny, double time) {
	std::string file_Sw = folder_path + "Sw.csv";
	std::string file_an_Sw = folder_path + "analytSw.csv";

	std::ofstream out_Sw(file_Sw);
	std::ofstream out_an_Sw(file_an_Sw);

	double Swf = 0.252; //
	double s = 1.0;
	int Np_s = 100.0; // num of intervals for ds
	double ds = (1.0 - Swf) / (double)Np_s;
	double xc = 0.0;
	for (int i = 0; i < 2*Np_s; i++) {
		if(i < Np_s){
			xc = q * time * dfwds(s) / phi / A;
			out_an_Sw << xc <<";"<< s << std::endl;
			s -= ds;
		}
		else {
			out_an_Sw << xc + h << ";" << Swr << std::endl;
			xc += h;
		}
	}

	for (int i = 0; i < Nx; ++i) {
		out_Sw << (double)(i * h) << ";" << x[i] << ";" << std::endl;
	}

}

int main() {

	std::vector<double> x; //solution at time n = 0
	std::vector<unsigned int> ia, ja; //.mtx matrix format
	std::vector<double> a, b; //Jacobi matrix and right half

	b.resize(N, 0.0);
	//x.resize(N, Swr); //solution at time n = 0
	x.resize(N, Swr); //solution at time n = 0
	x[0] = 1.0;

	BICGSTAB<ILDUC> Solver;
	//Solver.GetParameters().Save("params_default.txt");
	Solver.GetParameters().Load("params.txt");
	std::cout << "Loaded parameters: " << std::endl;
	Solver.GetParameters().Print();

	for(int tt = 0; tt < Nt; tt++){

		double error = 2.0 * tol;
		std::vector<double> xk = x;
		xk[0] = 1.0;
		while (error > tol) {
			a.clear();
			ia.clear();
			ja.clear();
			ia.resize(1, 0);
			std::fill(b.begin(), b.end(), 0);

			for (int i = 0; i < N; i++) {

				int bpos = ia.size() - 1;

					if (i == 0) {

						//Condition Sw = 1.0 on left boundary:
						///* 
							ja.push_back(i);
							a.push_back(1.0);

							b[bpos] = -(xk[i] - 1.0);
						//*/

						//Condition Swi - Sw(i + 1) = 0 on left boundary:
						/*
						ja.push_back(i);
						a.push_back(1.0);

						ja.push_back(i + 1);
						a.push_back(-1.0);

						b[bpos] = -(xk[i] - xk[i + 1]);
						*/
					}
					else {

						ja.push_back(i);
						a.push_back(dfds2(xk[i]));

						ja.push_back(i - 1);
						a.push_back(dfds1(xk[i - 1])); //0

						//std::cout << "dfds1 = " << dfds1(xk[i - 1]) << std::endl;
						//std::cout << "dfds2 = " << dfds2(xk[i]) << std::endl;

						b[bpos] += -(phi * (xk[i] - x[i]) / dt + q * (fw(xk[i]) - fw(xk[i - 1])) / A / h);
					}
					ia.push_back((int)ja.size()); //close row
			}

			CSRMatrix A(ia, ja, a);
			A.Save("matrix.txt");
			/*
			if (A.Symmetric(1.0e-5)) std::cout << "Symmetric matrix \n";
			else std::cout << " Not a symmetric matrix \n";
			*/

			std::vector<double> xk_new;
			if (Solver.Setup(A) && Solver.Solve(b, xk_new))
			{
				double r = Resid(A, b, xk_new);
				std::cout << "Final residual " << r << std::endl;
				if (r > 1000000.0) {
					std::cout << "ERROR: So big residual. Solution failed!" << std::endl;
					exit(-3);
				}

			}
			else {
				std::cout << "ERROR: Solution failed!" << std::endl;
				exit(-1);
			}

			//Calculate error:
			double err = 0.0;
			for (int i = 0; i < N; i++) {
				double dx = std::abs(xk_new[i]);
				if (dx > err) err = dx;
				xk_new[i] += xk[i]; //x(k + 1) = dx(k) + xk
			}
			error = err;
			xk = xk_new;
			std::cout << "error = " << error << std::endl;
		} // while()
		std::cout << " New Iter = " << tt << std::endl;
		x = xk; //x(n + 1) = x(k + 1);

		if(tt%OutIter == 0){
			std::string filename = "output" + std::to_string(tt) + ".vtk";
			generateVtkFile(filename, x, h, N, 20);
			filename = "output" + std::to_string(tt);
			generateCsvFile(filename, x, h, N, 20, (double)tt*dt);
		}
	}
	return 0;
}
