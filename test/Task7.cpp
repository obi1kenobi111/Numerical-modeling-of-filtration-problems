#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#include "bicgstab.h"
#include "ilduc.h"

using namespace std;
const double PI = 3.14159265358979323846;

const double L = 10000.0;
const double T = 10.0;
const double A = 1.0; //площадь поперечного сечения

//const double phi_m = 0.8, phi_f = 0.9;
const double phi_m = 0.2, phi_f = 0.01;
const double km = 0.01, kf = 50;
const double ct = 1e-6;

const double mu = 0.5; // ?
const double lambda = 1e-7; // 0.000001; // ?

const double dt = 1.0;
const double h = 100.0;
const double h2 = h * h;

const double p0 = 1000.0; // initial cond
const double p1 = 500.0; // left boundary cond
//num of nodes (1D - task):
const int N = (int)(L / h) + 1;
const int Nt = (int)(T / dt); //time steps

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
	vtkFile << "SCALARS Pm double" << std::endl;
	vtkFile << "LOOKUP_TABLE default" << std::endl;

	// Write solution values
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {

			vtkFile << x[i] << std::endl;
		}
	}

	vtkFile << "SCALARS Pf double" << std::endl;
	vtkFile << "LOOKUP_TABLE default" << std::endl;

	// Write solution values
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {

			vtkFile << x[i + Nx] << std::endl;
		}
	}

	vtkFile.close();
}

void generateCsvFile(std::string folder_path, std::vector<double> x, double h, int Nx, int Ny){
	std::string file_pf = folder_path + "pf.csv";
	std::string file_pm = folder_path + "pm.csv";

	std::ofstream out_pm(file_pm);
	std::ofstream out_pf(file_pf);

	for (int i = 0; i < Nx; ++i) {
			out_pm << (double)(i*h) << ";" << x[i] << std::endl;
	}

	// Write solution values
	for (int i = 0; i < Nx; ++i) {
			out_pf << (double)(i * h) << ";" << x[i + Nx] << std::endl;
	}

}

int main() {

	std::vector<double> x; //solution at time n = 0
	std::vector<unsigned int> ia, ja; //.mtx matrix format
	std::vector<double> a, b; //matrix and right half

	int Spm = 0;
	int Spf = N;
#define Ipm(i) (Spm + (i)) 
#define Ipf(i) (Spf + (i))
//Доступ к решению на n-ом шаге:
#define pm(i) (x[Ipm(i)])
#define pf(i) (x[Ipf(i)])

	int Np = 2 * N;
	b.resize(Np, 0.0);
	x.resize(Np, p0); //solution at time n = 0

	BICGSTAB<ILDUC> Solver;
	//Solver.GetParameters().Save("params_default.txt");
	Solver.GetParameters().Load("params.txt");
	std::cout << "Loaded parameters: " << std::endl;
	Solver.GetParameters().Print();

	for(int tt = 0; tt < Nt; tt++){

		a.clear();
		ia.clear();
		ja.clear();
		ia.resize(1, 0);
		std::fill(b.begin(), b.end(), 0);

		//Eq. №1 - dpm/dt:
		for (int i = 0; i < N; i++) {

				int bpos = ia.size() - 1;

				if (i == 0 || i == N - 1) {

					ja.push_back(Ipm(i));
					a.push_back(1.0);

					if(i == N - 1){
						ja.push_back(Ipm(i - 1));
						a.push_back(-1.0);
					}
					else if(i == 0) b[bpos] = p1;
				}
				else {
					//dpm/dt:
					/*
					ja.push_back(Ipm(i));
					a.push_back(phi_m * ct);
					*/
					double k1 = phi_m * ct / dt; 
					b[bpos] += pm(i);

					//We will divide all eq. on term k1 =>
					
					//-cm * d2(pm)/dx2 + dpm/dt:
					double cm = km / mu;
					ja.push_back(Ipm(i));
					a.push_back(2.0 * cm / h2 / k1 + 1.0 + 1.0 * km * lambda / mu / k1); //c1 - c1*p(n + 1))ij  - dp/dt

					ja.push_back(Ipm(i + 1));
					a.push_back(-1.0 * cm / h2 / k1);

					ja.push_back(Ipm(i - 1));
					a.push_back(-1.0 * cm / h2 / k1);

					// - (km/nu)*(pf - pm):
					ja.push_back(Ipf(i));
					a.push_back(-1.0 * km * lambda / mu / k1);

					//ja.push_back(Ipm(i));
					//a.push_back(1.0 * km * lambda / mu / k1);

				}

				ia.push_back((int)ja.size()); //close row
		}

		//Eq. №2 - dpf/dt:
		for (int i = 0; i < N; i++) {

			int bpos = ia.size() - 1;

			if (i == 0 || i == N - 1) {

				ja.push_back(Ipf(i));
				a.push_back(1.0);

				if (i == N - 1) {
					ja.push_back(Ipf(i - 1));
					a.push_back(-1.0);
				}
				else if (i == 0) b[bpos] = p1;
			}
			else {

				//dpmf/dt:
				/*
				ja.push_back(Ipf(i));
				a.push_back(phi_f * ct);
				*/

				double k2 = phi_f * ct / dt;
				b[bpos] += pf(i);

				//We will divide all eq. on term k2 =>
				// 
				//-cf * d2(pf)/dx2 + dpmf/ft
				double cf = kf / mu;
				ja.push_back(Ipf(i));
				a.push_back(2.0 * cf / h2 / k2 + 1.0 + 1.0 * kf * lambda / mu / k2);

				ja.push_back(Ipf(i + 1));
				a.push_back(-1.0 * cf / h2 / k2);

				ja.push_back(Ipf(i - 1));
				a.push_back(-1.0 * cf / h2 / k2);

				//Added to main term for stability:
				// - (km/nu)*(pf - pm): 
				//ja.push_back(Ipf(i));
				//a.push_back(1.0 * kf * lambda / mu / k2);

				ja.push_back(Ipm(i));
				a.push_back(-1.0 * kf * lambda / mu / k2);

			}

			ia.push_back((int)ja.size()); //close row
		}

		CSRMatrix A(ia, ja, a);
		A.Save("matrix.txt");
		/*
		if (A.Symmetric(1.0e-5)) std::cout << "Symmetric matrix \n";
		else std::cout << " Not a symmetric matrix \n";
		*/

		std::vector<double> x_new;
		if (Solver.Setup(A) && Solver.Solve(b, x_new))
		{
			double r = Resid(A, b, x_new);
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

		std::string filename = "output" + std::to_string(tt) + ".vtk";
		generateVtkFile(filename, x_new, h, N, 20);
		filename = "output" + std::to_string(tt);
		generateCsvFile(filename, x, h, N, 10);

		x = x_new;
	}
	return 0;
}
