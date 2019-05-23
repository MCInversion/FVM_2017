// FVM-3D-C++.cpp : Defines the entry point for the console application.
//

//#include "stdafx.h"
#include "pch.h"
#include <stdio.h>

#define _USE_MATH_DEFINES
#include <math.h>
#include <mpi.h>


#define RAD M_PI / 180

#define n1 60       /* number of divisions in R direction*/
#define n2 40       /* number of divisions in B direction*/
#define n3 20       /* number of divisions in L direction*/

#define NPROC 3

#define Lu 60.0     /* upper boundary for L*/
#define Ld 10.0     /* lower boundary for L*/
#define Bu 70.0 	/* upper boundary for B*/
#define Bd 20.0 	/* lower boundary for B*/

#define Ru 2.0 		/* upper boundary for R*/
#define Rd 1.0		/* lower boundary for R*/

#define GM 1.0
/* #define GM 3.986005e+14 */

/* parameters for SOR */
#define omega 1.9
#define tol 1.0e-20     /* tolerance */
#define max_it 10000	 /* number of iterations */

int main(int argc, char **argv) {
	// MPI Vars:
	int nprocs, myrank;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	long i, j, k, it;
	double deltaR, deltaL, deltaB, deltaL1, deltaB1;
	double res_local, pom_local, pom_global, pom1_local, z, sigma, res3_local, res_global, res3_global;

	typedef double pole1 [n1 / NPROC + 3][n2 + 2][n3 + 2];
	typedef double pole2 [n1 / NPROC + 3][n2 + 2][n3 + 2];
	typedef double pole0 [n2 + 2][n3 + 2];
	static pole1 B, L, R;
	static pole2 T_B, T_L, T_R, s, u, deltag;
	static pole2 an, as, aw, ae, au, ad, ap, b, res2;

	deltaL = (Lu - Ld) / n3;
	deltaB = (Bu - Bd) / n2;
	deltaR = (Ru - Rd) / n1;
	deltaL1 = deltaL * RAD * Rd;
	deltaB1 = deltaB * RAD * Rd;
	// printf("%.2lf %.2lf %.2lf\n\n", deltaL1, deltaB1, deltaR);

	// local/last packet sizes
	int istart, iend, nlocal, nlast;
	nlocal = n1 / nprocs;
	nlast = n1 - (nprocs - 1) * nlocal;

	// packet indexing
	istart = nlocal * myrank;
	if (myrank == nprocs - 1) {
		iend = n1 - 1;
	}
	else {
		iend = istart + nlocal - 1;
	}
		

	if (myrank == 0) {
		printf("===============================================\n");
		printf("--------- MPI params --------------------------\n");
		printf("nprocs = %d \n", nprocs);
		printf("p%d: nlocal = %d, nlast = %d\n", myrank, nlocal, nlast);
		printf("-----------------------------------------------\n");
	}
	printf("p%d: istart = %d, iend = %d\n\n", myrank, istart, iend);

	/*----------------------------------------------------------------------------*/
	/*vytvorenie 3D siete*/

	for (i = 1; i <= nlocal + 1; i++) {
		for (j = 1; j <= n2 + 1; j++) {
			for (k = 1; k <= n3 + 1; k++) {
				L[i][j][k] = Ld + (k - 1) * deltaL;
				B[i][j][k] = Bd + (j - 1) * deltaB;
				R[i][j][k] = Rd + (istart + i - 1) * deltaR;
			}
		}
	}

	/*----------------------------------------------------------------------------*/
	/*vypocet tazisk v BLR*/

	for (i = 1; i <= nlocal; i++) {
		for (j = 1; j <= n2; j++) {
			for (k = 1; k <= n3; k++) {
				T_L[i][j][k] = 0.5 * (L[i][j][k] + L[i][j][k + 1]);
				T_B[i][j][k] = 0.5 * (B[i][j][k] + B[i][j + 1][k]);
				T_R[i][j][k] = R[i][j][k] + 0.5 * deltaR;
			}
		}
	}

	/*----------------------------------------------------------------------------*/
	/*pridanie okrajovych tazisk*/

	for (i = 1; i <= nlocal; i++) {
		for (j = 1; j <= n2; j++) {
			T_B[i][j][0] = T_B[i][j][1];
			T_B[i][j][n3 + 1] = T_B[i][j][n3];

			T_L[i][j][0] = T_L[i][j][1] - deltaL;
			T_L[i][j][n3 + 1] = T_L[i][j][n3] + deltaL;

			T_R[i][j][0] = T_R[i][j][1];
			T_R[i][j][n3 + 1] = T_R[i][j][n3];
		}
	}

	for (i = 1; i <= nlocal; i++) {
		for (k = 1; k <= n3; k++) {
			T_B[i][0][k] = T_B[i][1][k] - deltaB;
			T_B[i][n2 + 1][k] = T_B[i][n2][k] + deltaB;

			T_L[i][0][k] = T_L[i][1][k];
			T_L[i][n2 + 1][k] = T_L[i][n2][k];

			T_R[i][0][k] = T_R[i][1][k];
			T_R[i][n2 + 1][k] = T_R[i][n2][k];
		}
	}

	if (myrank == 0) {
		for (j = 1; j <= n2; j++) {
			for (k = 1; k <= n3; k++) {
				T_B[0][j][k] = T_B[1][j][k];
				T_L[0][j][k] = T_L[1][j][k];
				T_R[0][j][k] = T_R[1][j][k] - deltaR;
			}
		}
	}	
	if (myrank == nprocs - 1) {
		for (j = 1; j <= n2; j++) {
			for (k = 1; k <= n3; k++) {
				T_B[nlast + 1][j][k] = T_B[nlast][j][k];
				T_L[nlast + 1][j][k] = T_L[nlast][j][k];
				T_R[nlast + 1][j][k] = T_R[nlast][j][k] + deltaR;
			}
		}
	}

	/*----------------------------------------------------------------------------*/
	/*vynulovanie poli*/

	for (i = 0; i <= nlocal + 1; i++) {
		for (j = 0; j <= n2 + 1; j++) {
			for (k = 0; k <= n3 + 1; k++) {
				u[i][j][k] = 0.0;
				aw[i][j][k] = 0.0;
				ae[i][j][k] = 0.0;
				as[i][j][k] = 0.0;
				an[i][j][k] = 0.0;
				au[i][j][k] = 0.0;
				ad[i][j][k] = 0.0;
				ap[i][j][k] = 0.0;
				b[i][j][k] = 0.0;
				s[i][j][k] = 0.0;
				res2[i][j][k] = 0.0;
			}
		}
	}

	/*----------------------------------------------------------------------------*/
	/*vypocet koeficientov*/

	for (i = 1; i <= nlocal; i++) {
		for (j = 1; j <= n2; j++) {
			for (k = 1; k <= n3; k++) {
				//x = r*Cos(B)*Cos(L)
				//y = r*Cos(B)*Sin(L)
				//z = r*Sin(B)

				//koeficienty pri susedovi sa pocitaju ako plocha steny lomeno vzdialenost k susedovi

				//Obsah W a E steny = deltaB/2*(Rspodne^2-Rvrchne^2), ... Obsah W a E steny sa meni len s R. RAD preraba stupne na radiany.
				// vzdialenosti medzi bodmi v smere W a E = deltaL*R*cos(B), ... vzdialenosti sa menia s R a B lebo cim sme blizsie k polom tak sa vzdianenosti zmensuju (lebo su kratsie rovnobezky)
				aw[i][j][k] = (deltaB / 2. * RAD * (R[i + 1][j][k] * R[i + 1][j][k] - R[i][j][k] * R[i][j][k])) / (deltaL * RAD * T_R[i][j][k] * cos(T_B[i][j][k] * RAD));
				ae[i][j][k] = aw[i][j][k];

				//Obsah N a S steny = deltaL/2*(Rspodne^2-Rvrchne^2)*cos(B), ... Obsah N a S steny sa meni s R ale aj s B lebo cim sme blizsie k polom tak sa obsahy stien zmensuju (lebo su kratsie rovnobezky)
				// vzdialenosti medzi bodmi v smere N a S = deltaB*R, ... vzdialenosti sa menia s R lebo poludniky su vzdy rovnako dlhe
				an[i][j][k] = (deltaL / 2. * RAD * (R[i + 1][j][k] * R[i + 1][j][k] - R[i][j][k] * R[i][j][k]) * cos(B[i][j + 1][k] * RAD)) / (deltaB * RAD * T_R[i][j][k]);
				as[i][j][k] = (deltaL / 2. * RAD * (R[i + 1][j][k] * R[i + 1][j][k] - R[i][j][k] * R[i][j][k]) * cos(B[i][j][k] * RAD)) / (deltaB * RAD * T_R[i][j][k]);

				//Obsah D steny = deltaL * (sin(Be) - sin(Bw)) * Rspodne^2
				//Obsah U steny = deltaL * (sin(Be) - sin(Bw)) * Rvrchne^2
				//vzdialenost medzi bodmi v smere U a D = deltaR
				ad[i][j][k] = deltaL * RAD * (sin(B[i + 1][j + 1][k] * RAD) - sin(B[i + 1][j][k] * RAD)) * R[i][j][k] * R[i][j][k] / deltaR;
				au[i][j][k] = deltaL * RAD * (sin(B[i][j + 1][k] * RAD) - sin(B[i][j][k] * RAD)) * R[i + 1][j][k] * R[i + 1][j][k] / deltaR;

				//vsetky aw, ae, an, as, ad, au by mali byt s minuskou ale v SOR sa to potom kompenzuje.
				ap[i][j][k] = -(aw[i][j][k] + ae[i][j][k] + as[i][j][k] + an[i][j][k] + au[i][j][k] + ad[i][j][k]);
			}
		}
	}

	/*----------------------------------------------------------------------------*/
	/*nacitanie okrajovych podmienok - Neumann*/

	// presne riesenie je GM/R, derivacia v smere R je -GM/R^2
	// na spodnej hranici je dana neumannova okrajova podmienka vsade inde dirichletova okrajova podmienka

	if (myrank == 0) {
		for (j = 1; j <= n2; j++) {
			for (k = 1; k <= n3; k++) {
				//deltag je derivacia v radialnom smere (t.j. v smere normaly k D stene)
				deltag[1][j][k] = GM / (R[1][j][k] * R[1][j][k]);
				//kedze derivaciu v smere normaly k dolnej stene pozname, mozme ju prenasobit obsahom D steny a prehodit na pravu stranu
				b[1][j][k] = b[1][j][k] - deltaL * RAD * (sin(B[1][j + 1][k] * RAD) - sin(B[1][j][k] * RAD)) * R[1][j][k] * R[1][j][k] * deltag[1][j][k];
				//kedze sme ju prehodili na pravu stranu, treba ju nulovat medzi neznamymi
				ad[1][j][k] = 0;
				ap[1][j][k] = -(aw[1][j][k] + ae[1][j][k] + as[1][j][k] + an[1][j][k] + au[1][j][k]);
			}
		}
	}

	/*----------------------------------------------------------------------------*/
	/*nacitanie okrajovych podmienok - Dirichlet*/

	if (myrank == nprocs - 1) {
		for (j = 1; j <= n2; j++) {
			for (k = 1; k <= n3; k++) {
				//pom je presne riesenie v bode [n1+1][j][k]
				pom_local = GM / T_R[nlast + 1][j][k];
				//kedze mame presne riesenie jedneho suseda, uz to nie je neznama a mozme ju prehodit na druhu stranu
				b[nlast][j][k] = b[nlast][j][k] - pom_local * au[nlast][j][k];
				//kedze sme ju prehodili na pravu stranu, treba ju nulovat medzi neznamymi
				au[nlast][j][k] = 0;
			}
		}
	}

	for (i = 1; i <= nlocal; i++) {
		for (k = 1; k <= n3; k++) {
			//pom je presne riesenie v bode [i][0][k]
			pom_local = GM / T_R[i][0][k];
			//kedze mame presne riesenie jedneho suseda, uz to nie je neznama a mozme ju prehodit na druhu stranu
			b[i][1][k] = b[i][1][k] - pom_local * as[i][1][k];
			//kedze sme ju prehodili na pravu stranu, treba ju nulovat medzi neznamymi
			as[i][1][k] = 0;

			//pom je presne riesenie v bode [i][n2+1][k]
			pom_local = GM / T_R[i][n2 + 1][k];
			//kedze mame presne riesenie jedneho suseda, uz to nie je neznama a mozme ju prehodit na druhu stranu
			b[i][n2][k] = b[i][n2][k] - pom_local * an[i][n2][k];
			//kedze sme ju prehodili na pravu stranu, treba ju nulovat medzi neznamymi
			an[i][n2][k] = 0;
		}
	}

	for (i = 1; i <= nlocal; i++) {
		for (j = 1; j <= n2; j++) {
			//pom je presne riesenie v bode [i][j][0]
			pom_local = GM / T_R[i][j][0];
			//kedze mame presne riesenie jedneho suseda, uz to nie je neznama a mozme ju prehodit na druhu stranu
			b[i][j][1] = b[i][j][1] - pom_local * aw[i][j][1];
			//kedze sme ju prehodili na pravu stranu, treba ju nulovat medzi neznamymi
			aw[i][j][1] = 0;

			//pom je presne riesenie v bode [i][j][n3+1]
			pom_local = GM / T_R[i][j][n3 + 1];
			//kedze mame presne riesenie jedneho suseda, uz to nie je neznama a mozme ju prehodit na druhu stranu
			b[i][j][n3] = b[i][j][n3] - pom_local * ae[i][j][n3];
			//kedze sme ju prehodili na pravu stranu, treba ju nulovat medzi neznamymi
			ae[i][j][n3] = 0;
		}
	}


	/*----------------------------------------------------------------------------*/
	/*riesenie sustavy rovnic*/

	pom_local = 0.0;
	it = 0;	
	MPI_Status status;
	// int nreq = 2 * (nprocs - 1);
	// MPI_Status* statuses = new MPI_Status [nreq];
	// MPI_Request* requests = new MPI_Request [nreq];

	if (myrank == 0) printf("p%d:    it:       res: \n", myrank);

	do {
		it++;
		// even blocks ================================

		// send boundary coeffs to neighboring processes
		if (myrank > 0) {
			MPI_Send(&u[1][0][0], (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank - 1, 0, MPI_COMM_WORLD);
		}
		if (myrank < nprocs - 1) {
			MPI_Send(&u[nlocal][0][0], (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank + 1, 1, MPI_COMM_WORLD);
		}

		if (myrank < nprocs - 1) {
			MPI_Recv(&u[nlocal + 1][0][0], (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank + 1, 0, MPI_COMM_WORLD, &status);
		}
		if (myrank > 0) {
			MPI_Recv(&u[0][0][0], (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank - 1, 1, MPI_COMM_WORLD, &status);
		}
		/*
		if (myrank > 0) {
			MPI_Isend(&u[1][0][0], (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank - 1, 0, MPI_COMM_WORLD, &requests[myrank - 1]);
		}
		if (myrank < nprocs - 1) {
			MPI_Isend(&u[nlocal][0][0], (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank + 1, 1, MPI_COMM_WORLD, &requests[nprocs + myrank]);
		}

		if (myrank < nprocs - 1) {
			MPI_Irecv(&u[nlocal + 1][0][0], (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank + 1, 0, MPI_COMM_WORLD, &requests[myrank + 1]);
		}
		if (myrank > 0) {
			MPI_Irecv(&u[0][0][0], (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank - 1, 1, MPI_COMM_WORLD, &requests[nprocs + myrank - 2]);
		}
		MPI_Waitall(nreq, requests, statuses);*/

		for (i = 1; i <= nlocal; i++) {
			for (j = 1; j <= n2; j++) {
				for (k = 1; k <= n3; k++) {
					if ((istart + i + j + k) % 2 == 0) {
						z = (b[i][j][k] - u[i][j][k + 1] * ae[i][j][k]
							- u[i][j][k - 1] * aw[i][j][k]
							- u[i][j + 1][k] * an[i][j][k]
							- u[i][j - 1][k] * as[i][j][k]
							- u[i + 1][j][k] * au[i][j][k]
							- u[i - 1][j][k] * ad[i][j][k]) / ap[i][j][k];
						u[i][j][k] = u[i][j][k] + omega * (z - u[i][j][k]);
					}
				}
			}
		}

		// odd blocks

		// send boundary coeffs to neighboring processes
		if (myrank > 0) {
			MPI_Send(&u[1][0][0], (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank - 1, 0, MPI_COMM_WORLD);
		}
		if (myrank < nprocs - 1) {
			MPI_Send(&u[nlocal][0][0], (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank + 1, 1, MPI_COMM_WORLD);
		}

		if (myrank < nprocs - 1) {
			MPI_Recv(&u[nlocal + 1][0][0], (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank + 1, 0, MPI_COMM_WORLD, &status);
		}
		if (myrank > 0) {
			MPI_Recv(&u[0][0][0], (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank - 1, 1, MPI_COMM_WORLD, &status);
		}
		/*
		if (myrank > 0) {
			MPI_Isend(&u[1][0][0], (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank - 1, 0, MPI_COMM_WORLD, &requests[myrank - 1]);
		}
		if (myrank < nprocs - 1) {
			MPI_Isend(&u[nlocal][0][0], (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank + 1, 1, MPI_COMM_WORLD, &requests[nprocs + myrank]);
		}

		if (myrank < nprocs - 1) {
			MPI_Irecv(&u[nlocal + 1][0][0], (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank + 1, 0, MPI_COMM_WORLD, &requests[myrank + 1]);
		}
		if (myrank > 0) {
			MPI_Irecv(&u[0][0][0], (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank - 1, 1, MPI_COMM_WORLD, &requests[nprocs + myrank - 2]);
		}
		MPI_Waitall(nreq, requests, statuses);*/

		for (i = 1; i <= nlocal; i++) {
			for (j = 1; j <= n2; j++) {
				for (k = 1; k <= n3; k++) {
					if ((istart + i + j + k) % 2 == 1) {
						z = (b[i][j][k] - u[i][j][k + 1] * ae[i][j][k]
							- u[i][j][k - 1] * aw[i][j][k]
							- u[i][j + 1][k] * an[i][j][k]
							- u[i][j - 1][k] * as[i][j][k]
							- u[i + 1][j][k] * au[i][j][k]
							- u[i - 1][j][k] * ad[i][j][k]) / ap[i][j][k];
						u[i][j][k] = u[i][j][k] + omega * (z - u[i][j][k]);
					}
				}
			}
		}

		// send boundary coeffs to neighboring processes
		if (myrank > 0) {
			MPI_Send(&u[1][0][0], (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank - 1, 0, MPI_COMM_WORLD);
		}
		if (myrank < nprocs - 1) {
			MPI_Send(&u[nlocal][0][0], (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank + 1, 1, MPI_COMM_WORLD);
		}

		if (myrank < nprocs - 1) {
			MPI_Recv(&u[nlocal + 1][0][0], (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank + 1, 0, MPI_COMM_WORLD, &status);
		}
		if (myrank > 0) {
			MPI_Recv(&u[0][0][0], (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank - 1, 1, MPI_COMM_WORLD, &status);
		}
		/*
		if (myrank > 0) {
			MPI_Isend(&u[1][0][0], (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank - 1, 0, MPI_COMM_WORLD, &requests[myrank - 1]);
		}
		if (myrank < nprocs - 1) {
			MPI_Isend(&u[nlocal][0][0], (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank + 1, 1, MPI_COMM_WORLD, &requests[nprocs + myrank]);
		}

		if (myrank < nprocs - 1) {
			MPI_Irecv(&u[nlocal + 1][0][0], (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank + 1, 0, MPI_COMM_WORLD, &requests[myrank + 1]);
		}
		if (myrank > 0) {
			MPI_Irecv(&u[0][0][0], (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank - 1, 1, MPI_COMM_WORLD, &requests[nprocs + myrank - 2]);
		}
		MPI_Waitall(nreq, requests, statuses);*/

		res_local = res3_local = 0.0;

		for (i = 1; i <= nlocal; i++) {
			for (j = 1; j <= n2; j++) {
				for (k = 1; k <= n3; k++) {
					pom_local = (u[i][j][k] * ap[i][j][k]
						+ u[i][j][k + 1] * ae[i][j][k]
						+ u[i][j][k - 1] * aw[i][j][k]
						+ u[i][j + 1][k] * an[i][j][k]
						+ u[i][j - 1][k] * as[i][j][k]
						+ u[i + 1][j][k] * au[i][j][k]
						+ u[i - 1][j][k] * ad[i][j][k] - b[i][j][k]);

					res_local += pom_local * pom_local;

					pom1_local = GM / T_R[i][j][k] - u[i][j][k];
					res3_local += pom1_local * pom1_local;
				}
			}
		}

		MPI_Allreduce(&res3_local, &res3_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&res_local, &res_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		res3_global = sqrt(res3_global / (n1 * n2 * n3));
		if (myrank == 0) printf("p%d:\t%d\t%.12lf\n", myrank, it, res_global);

	} while ((res_global > tol) && (it < max_it));


	/*----------------------------------------------------------------------------*/
	/*porovnanie s presnym riesenim a zapis vysledkov*/

	sigma = 0.0;
	pom_local = 0.0;

	for (i = 1; i <= nlocal; i++) {
		for (j = 1; j <= n2; j++) {
			for (k = 1; k <= n3; k++) {
				res2[i][j][k] = (GM / T_R[i][j][k]) - u[i][j][k];
				pom_local +=
					res2[i][j][k] * res2[i][j][k] * (deltaL * (pow(T_R[i][j][k] - deltaR / 2., 3.)
						- pow(T_R[i][j][k] + deltaR / 2., 3.)) * (sin(T_B[i][j][k] + deltaB / 2.)
							- sin(T_B[i][j][k] - deltaB / 2.))) / 3.;
				//	fprintf(fw,"%d %d %d\t%.7lf\t%.9lf\n",i,j,k,u[i][j][k],res2[i][j][k]);
			}
		}
	}

	MPI_Allreduce(&pom_local, &pom_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	if (myrank == 0) {
		sigma = sqrt(pom_global);
		printf("sigma = %.20lf\n", sigma);
		FILE *fw; // *fr,
		fopen_s(&fw, "stat.txt", "w");

		i = 1;
		for (j = 1; j <= n2; j++) {
			for (k = 1; k <= n3; k++) {
				res2[i][j][1] = (GM / T_R[i][j][k]) - u[i][j][k];
				fprintf(fw, "%d %d %d\t%.7lf\t%.9lf\t%.7lf\n", i, j, k, u[i][j][j], res2[i][j][j], (GM / T_R[i][j][k]));
			}
		}

		fclose(fw);
	}

	MPI_Finalize();

	// getchar();

	return 0;
}

