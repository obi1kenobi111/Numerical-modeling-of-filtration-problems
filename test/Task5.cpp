#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#include "bicgstab.h"
#include "ilduc.h"

using namespace std;
const double PI = 3.14159265358979323846;

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
            vtkFile << (double)(i*h) << " " << (double)(j * h) << " 0" << std::endl;
        }
    }

    vtkFile << "POINT_DATA " << Nx * Ny << std::endl;
    vtkFile << "SCALARS PRESSURE double" << std::endl;
    vtkFile << "LOOKUP_TABLE default" << std::endl;

    // Write solution values
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {

            vtkFile << x[i * Ny + j] << std::endl;
        }
    }

    vtkFile.close();
}

int main() {
	const double L = 10000.0; // 1.0;
	const double h = 1000.0; // 0.1;
	const double h2 = h * h;
	//num of nodes:
	int Nx = 11;// (int)(L / h) + 1;
	int Ny = 11;// Nx;
	std::vector<double> x; //solution at time n = 0
	std::vector<unsigned int> ia, ja; //.mtx matrix format
	std::vector<double> a, b; //matrix and right half
	int Np = Nx * Ny;

#define Ip(i, j) ((i)*Ny + (j)) 
#define p(i,j) (x[Ip((i),(j))])

	double rw = 0.25;
	double re = 0.2 * h;
	double pbh = 600.0; //[psi] 800
	double p0 = 1000.0; // initial pressure [psi]
	double p_left = 2000.0; //[psi]
	double qc = 2.0*PI/std::log(re/rw);

	double dt = 1.0;
	int Nt = 20; //time steps

	double q = 0.05;//1e-5;
	double bw = 1.0;
	double phi = 0.2;
	double ct = 1e-6;
	double k = 50.0 / 300.0;
	double mu_w = 1.0;

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

		for (int i = 0; i < Nx; i++) {
			for (int j = 0; j < Ny; j++) {

				int bpos = ia.size() - 1;

				if (i == 0 || j == 0 || i == Nx - 1 || j == Ny - 1) {

					ja.push_back(Ip(i, j));
					a.push_back(1.0);

					if (i == 0) {
						ja.push_back(Ip(i + 1, j));
						a.push_back(-1.0);
					}
					else if(j == 0){
						ja.push_back(Ip(i, j + 1));
						a.push_back(-1.0);
					}
					else if (i == Nx - 1) { //Direchlet condition
						b[bpos] = p_left;
					}
					else if (j == Ny - 1) {
						ja.push_back(Ip(i, j - 1));
						a.push_back(-1.0);
					}
				}
				else {
					//dp/dt:
					/*
					ja.push_back(Ip(i, j));
					a.push_back(c1);
					b[bpos] += c1 * p(i, j);
					*/

					double c1 = phi * ct / bw;
					b[bpos] += c1* p(i, j);

					//d2p/dx2 + d2p/dy2  + c1*p(n + 1):
					double c2 = dt * k / mu_w / bw;
					ja.push_back(Ip(i, j));
					a.push_back(4.0 * c2 / h2 + c1); //c1 - c1*p(n + 1))ij  - dp/dt

					ja.push_back(Ip(i + 1, j));
					a.push_back(-1.0 * c2 / h2);

					ja.push_back(Ip(i - 1, j));
					a.push_back(-1.0 * c2 / h2);

					ja.push_back(Ip(i, j + 1));
					a.push_back(-1.0 * c2 / h2);

					ja.push_back(Ip(i, j - 1));
					a.push_back(-1.0 * c2 / h2);

					//нагнетательная скважина in point (5000ft, 5000ft) 
					///*
					if (i == 5 && j == 5) {
						b[bpos] += q * dt;
					}

					//производящая скважина at point (9000ft, 9000ft)
					if (i == 9 && j == 9) {
						ja.push_back(Ip(i, j));
						a.push_back(qc * dt);

						b[bpos] += qc * pbh * dt;
					}
					//*/

					/*
					if (i == 2 && j == 2) {
						ja.push_back(Ip(i, j));
						a.push_back(qc * dt);

						b[bpos] += qc * (pbh) * dt;
					}

					//well injector at point (0.75m, 0.75m)
					if (i == 6 && j == 6) {
						ja.push_back(Ip(i, j));
						a.push_back(qc * dt);

						b[bpos] += qc * (-pbh)*dt;
					}
					*/
				}

				ia.push_back((int)ja.size()); //close row
			}
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
		generateVtkFile(filename, x_new, h, Nx, Ny);

		x = x_new;
	}
	return 0;
}
