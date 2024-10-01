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
	const double L = 1.0;
	const double h = 0.125;
	const double h2 = h * h;
	//num of nodes:
	int Nx = 9;// (int)(L / h) + 1;
	int Ny = 9;// Nx;
	std::vector<double> x; //solution at time n = 0
	std::vector<unsigned int> ia, ja; //.mtx matrix format
	std::vector<double> a, b; //matrix and right half
	int Np = Nx * Ny;
	b.resize(Np, 0.0);
	x.resize(Np, 0.0); //solution at time n = 0

#define Ip(i, j) ((i)*Ny + (j)) 
#define p(i,j) (x[Ip((i),(j))])
	double rw = 0.001;
	double re = 0.2 * h;
	double pb = 1.0;
	double qc = 2.0*PI/std::log(re/rw);
	// N x N regular grid with dx = dy = h
		//   *------------*
		//   |            |
		//   |  p(i,j)    |
		//   |            |
		//   *------------*
		// p_xx: - p(i-1,j) + 2 p(i,j) - p(i+1,j)
		// p_xx + p_yy
		//            -1 p(i,j+1)
		//-1 p(i-1,j) +4 p(i,j)   - 1 p(i+1,j)
		//            -1 p(i,j-1)
	const double coeffs[9] =
	{
		 0.0,-1.0, 0.0,
		-1.0, 0.0,-1.0,
		 0.0,-1.0, 0.0,
	};
#define Is(ii,jj) (((ii)+1)*3 + ((jj)+1))

	a.clear();
	ia.clear();
	ja.clear();
	ia.resize(1, 0);
	std::fill(b.begin(), b.end(), 0);


	BICGSTAB<ILDUC> Solver;
	//Solver.GetParameters().Save("params_default.txt");
	Solver.GetParameters().Load("params.txt");
	std::cout << "Loaded parameters: " << std::endl;
	Solver.GetParameters().Print();

	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {

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
				else if (i == Nx - 1) {
					ja.push_back(Ip(i - 1, j));
					a.push_back(-1.0);
				}
				else if (j == Ny - 1) {
					ja.push_back(Ip(i, j - 1));
					a.push_back(-1.0);
				}
			}
			else {
				ja.push_back(Ip(i, j));
				a.push_back(4.0);

				ja.push_back(Ip(i + 1, j));
				a.push_back(-1.0);

				ja.push_back(Ip(i - 1, j));
				a.push_back(-1.0);

				ja.push_back(Ip(i, j + 1));
				a.push_back(-1.0);

				ja.push_back(Ip(i, j - 1));
				a.push_back(-1.0);

				int bpos = ia.size() - 1;
				//well producer in point (0.25m, 0.25m) 
				if (i == 2 && j == 2) {
					ja.push_back(Ip(i, j));
					a.push_back(qc * h2);

					b[bpos] += qc * (pb)*h2;
				}
				
				//well injector at point (0.75m, 0.75m)
				if (i == 6 && j == 6) {
					ja.push_back(Ip(i, j));
					a.push_back(qc * h2);

					b[bpos] += qc * (-pb) * h2;
				}

			}

			ia.push_back((int)ja.size()); //close row
		}
	}

	//std::cout << "ia.size() = " << ia.size() << std::endl;
	//std::cout << "ja.size() = " << ja.size() << std::endl;
	CSRMatrix A(ia, ja, a);

	if (A.Symmetric(1.0e-5)) std::cout << "Symmetric matrix \n";
	else std::cout << " Not a symmetric matrix \n";

	std::vector<double> x_new;
	if (Solver.Setup(A) && Solver.Solve(b, x_new))
	{
		double r = Resid(A, b, x_new);
		std::cout << "Final residual " << r << std::endl;
		if (r > 1000000.0) {
			std::cout << "ERROR: So big residual. Solution failed!" << std::endl;
			exit(-3);
		}
		//SaveVector(std::string("solution"), x_new);
		//SaveMyVector(std::string("result.txt"), 5, x_new);

	}
	else {
		std::cout << "ERROR: Solution failed!" << std::endl;
		exit(-1);
	}

	generateVtkFile("output.vtk", x_new, h, Nx, Ny);
	return 0;
}
