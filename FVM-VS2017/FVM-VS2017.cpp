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
	const int NPROC = 4;
	double deltaR, deltaL, deltaB, deltaL1, deltaB1;
	double res_local, pom_local, pom_global, pom1_local, z, sigma, res3_local, res_global, res3_global;

	// typedef double pole1 [n1 / NPROC + 3][n2 + 2][n3 + 2];
	// typedef double pole2 [n1 / NPROC + 3][n2 + 2][n3 + 2];
	// static pole1 B, L, R;
	// static pole2 T_B, T_L, T_R, s, u, deltag;
	// static pole2 an, as, aw, ae, au, ad, ap, b, res2;

	double ***B = new double **[n1 / NPROC + 3];
	double ***L = new double **[n1 / NPROC + 3];
	double ***R = new double **[n1 / NPROC + 3];

	double ***T_B = new double **[n1 / NPROC + 3];
	double ***T_L = new double **[n1 / NPROC + 3];
	double ***T_R = new double **[n1 / NPROC + 3];

	double ***s = new double **[n1 / NPROC + 3];
	double ***u = new double **[n1 / NPROC + 3];
	double ***deltag = new double **[n1 / NPROC + 3];

	double ***an = new double **[n1 / NPROC + 3];
	double ***as = new double **[n1 / NPROC + 3];
	double ***aw = new double **[n1 / NPROC + 3];
	double ***ae = new double **[n1 / NPROC + 3];
	double ***au = new double **[n1 / NPROC + 3];
	double ***ad = new double **[n1 / NPROC + 3];
	double ***ap = new double **[n1 / NPROC + 3];

	double ***b = new double **[n1 / NPROC + 3];
	double ***res2 = new double **[n1 / NPROC + 3];

	deltaL = (Lu - Ld) / n3;
	deltaB = (Bu - Bd) / n2;
	deltaR = (Ru - Rd) / n1;
	deltaL1 = deltaL * RAD * Rd;
	deltaB1 = deltaB * RAD * Rd;
	printf("%.2lf %.2lf %.2lf\n\n", deltaL1, deltaB1, deltaR);

	// local/last packet sizes
	int istart, iend, nlocal = 0, nlast = 0;
	nlocal = (n1 / nprocs) + 1;
	nlast = n1 - (nprocs - 1) * nlocal;

	// packet indexing
	istart = nlocal * myrank;
	if (myrank == nprocs - 1)
		iend = n1 - 1;
	else
		iend = istart + nlocal - 1;

	if (myrank == 0) {
		printf("===============================================\n");
		printf("--------- MPI params --------------------------\n");
		printf("nprocs = %d \n", nprocs);
		printf("p%d: nlocal = %d, nlast = %d\n", myrank, nlocal, nlast);
	}
	printf("p%d: istart = %d, iend = %d\n\n", myrank, istart, iend);

	/*----------------------------------------------------------------------------*/
	/*vytvorenie 3D siete*/

	for (i = 0; i <= nlocal + 1; i++) {
		L[i] = new double *[n2 + 2];
		B[i] = new double *[n2 + 2];
		R[i] = new double *[n2 + 2];
		for (j = 1; j <= n2 + 1; j++) {
			L[i][j] = new double[n3 + 2];
			B[i][j] = new double[n3 + 2];
			R[i][j] = new double[n3 + 2];
			for (k = 1; k <= n3 + 1; k++) {
				L[i][j][k] = Ld + (k - 1) * deltaL;
				B[i][j][k] = Bd + (j - 1) * deltaB;
				R[i][j][k] = Rd + (istart + i - 1) * deltaR;
			}
		}
	}

	// process layer boundaries (upper and lower)
	double **L_bdL = L[0];
	double **L_bdU = L[nlocal + 1];

	double **B_bdL = B[0];
	double **B_bdU = B[nlocal + 1];

	double **R_bdL = R[0];
	double **R_bdU = R[nlocal + 1];

	/*----------------------------------------------------------------------------*/
	/*vypocet tazisk v BLR*/

	for (i = 1; i <= nlocal; i++) {
		T_L[i] = new double *[n2 + 2];
		T_B[i] = new double *[n2 + 2];
		T_R[i] = new double *[n2 + 2];
		for (j = 1; j <= n2; j++) {
			T_L[i][j] = new double[n3 + 2];
			T_B[i][j] = new double[n3 + 2];
			T_R[i][j] = new double[n3 + 2];
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
		T_L[0] = new double *[n2 + 2];
		T_B[0] = new double *[n2 + 2];
		T_R[0] = new double *[n2 + 2];
	}
	for (j = 1; j <= n2; j++) {
		if (myrank == 0) {
			T_L[0][j] = new double[n3 + 2];
			T_B[0][j] = new double[n3 + 2];
			T_R[0][j] = new double[n3 + 2];
		}
		for (k = 1; k <= n3; k++) {
			if (myrank == 0) {
				T_B[0][j][k] = T_B[1][j][k];
				T_L[0][j][k] = T_L[1][j][k];
				T_R[0][j][k] = T_R[1][j][k] - deltaR;
			}

			if (myrank == nprocs - 1) {
				T_B[nlocal + 1][j][k] = T_B[nlocal][j][k];
				T_L[nlocal + 1][j][k] = T_L[nlocal][j][k];
				T_R[nlocal + 1][j][k] = T_R[nlocal][j][k] + deltaR;
			}
		}
	}

	// process layer (centers of mass) boundaries (upper and lower)
	double **TL_bdU = T_L[nlocal + 1];
	double **TL_bdL = T_L[0];

	double **TB_bdU = T_B[nlocal + 1];
	double **TB_bdL = T_B[0];

	double **TR_bdU = T_R[nlocal + 1];
	double **TR_bdL = T_R[0];

	/*----------------------------------------------------------------------------*/
	/*vynulovanie poli*/

	for (i = 0; i <= nlocal + 1; i++) {
		aw[i] = new double *[n2 + 2]; ae[i] = new double *[n2 + 2];
		an[i] = new double *[n2 + 2]; as[i] = new double *[n2 + 2];
		au[i] = new double *[n2 + 2]; ad[i] = new double *[n2 + 2];

		ap[i] = new double *[n2 + 2]; u[i] = new double *[n2 + 2];
		b[i] = new double *[n2 + 2]; s[i] = new double *[n2 + 2]; res2[i] = new double *[n2 + 2];
		for (j = 0; j <= n2 + 1; j++) {
			aw[i][j] = new double[n3 + 2]; ae[i][j] = new double[n3 + 2];
			an[i][j] = new double[n3 + 2]; as[i][j] = new double[n3 + 2];
			au[i][j] = new double[n3 + 2]; ad[i][j] = new double[n3 + 2];

			ap[i][j] = new double [n3 + 2]; u[i][j] = new double [n3 + 2];
			b[i][j] = new double [n3 + 2]; s[i][j] = new double [n3 + 2]; res2[i][j] = new double [n3 + 2];
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

	// send boundary layers to neighboring processes
	if (myrank > 0) {
		MPI_Send(L_bdL, (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank - 1, 0, MPI_COMM_WORLD);
		MPI_Send(B_bdL, (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank - 1, 1, MPI_COMM_WORLD);
		MPI_Send(R_bdL, (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank - 1, 2, MPI_COMM_WORLD);

		MPI_Send(TL_bdL, (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank - 1, 3, MPI_COMM_WORLD);
		MPI_Send(TB_bdL, (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank - 1, 4, MPI_COMM_WORLD);
		MPI_Send(TR_bdL, (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank - 1, 5, MPI_COMM_WORLD);
	}
	if (myrank < nprocs - 1) {
		MPI_Send(L_bdU, (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank + 1, 6, MPI_COMM_WORLD);
		MPI_Send(B_bdU, (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank + 1, 7, MPI_COMM_WORLD);
		MPI_Send(R_bdU, (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank + 1, 8, MPI_COMM_WORLD);

		MPI_Send(TL_bdU, (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank + 1, 9, MPI_COMM_WORLD);
		MPI_Send(TB_bdU, (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank + 1, 10, MPI_COMM_WORLD);
		MPI_Send(TR_bdU, (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank + 1, 11, MPI_COMM_WORLD);
	}

	// recv buffers
	double **recv_L_bdL;
	double **recv_L_bdU;
	double **recv_B_bdL;
	double **recv_B_bdU;
	double **recv_R_bdL;
	double **recv_R_bdU;

	double **recv_TL_bdU;
	double **recv_TL_bdL;
	double **recv_TB_bdU;
	double **recv_TB_bdL;
	double **recv_TR_bdU;
	double **recv_TR_bdL;

	MPI_Status status;

	if (myrank < nprocs - 1) {
		MPI_Recv(recv_L_bdL, (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank + 1, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(recv_B_bdL, (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank + 1, 1, MPI_COMM_WORLD, &status);
		MPI_Recv(recv_R_bdL, (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank + 1, 2, MPI_COMM_WORLD, &status);

		MPI_Recv(recv_TL_bdL, (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank + 1, 3, MPI_COMM_WORLD, &status);
		MPI_Recv(recv_TB_bdL, (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank + 1, 4, MPI_COMM_WORLD, &status);
		MPI_Recv(recv_TR_bdL, (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank + 1, 5, MPI_COMM_WORLD, &status);
	}
	if (myrank > 0) {
		MPI_Recv(recv_L_bdU, (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank - 1, 6, MPI_COMM_WORLD, &status);
		MPI_Recv(recv_B_bdU, (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank - 1, 7, MPI_COMM_WORLD, &status);
		MPI_Recv(recv_R_bdU, (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank - 1, 8, MPI_COMM_WORLD, &status);

		MPI_Recv(recv_TL_bdU, (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank - 1, 9, MPI_COMM_WORLD, &status);
		MPI_Recv(recv_TB_bdU, (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank - 1, 10, MPI_COMM_WORLD, &status);
		MPI_Recv(recv_TR_bdU, (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank - 1, 11, MPI_COMM_WORLD, &status);
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

	/*----------------------------------------------------------------------------*/
	/*nacitanie okrajovych podmienok - Dirichlet*/

	for (j = 1; j <= n2; j++) {
		for (k = 1; k <= n3; k++) {
			//pom je presne riesenie v bode [n1+1][j][k]
			pom_local = GM / T_R[n1 + 1][j][k];
			//kedze mame presne riesenie jedneho suseda, uz to nie je neznama a mozme ju prehodit na druhu stranu
			b[n1][j][k] = b[n1][j][k] - pom_local * au[n1][j][k];
			//kedze sme ju prehodili na pravu stranu, treba ju nulovat medzi neznamymi
			au[n1][j][k] = 0;
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

	// process layer (coefficients au, ad) boundaries (upper and lower)
	double **au_bd = au[nlocal];
	double **ad_bd = ad[0];

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

	// send boundary coeffs to neighboring processes
	if (myrank > 0) {
		MPI_Send(ad_bd, (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank - 1, 0, MPI_COMM_WORLD);
	}
	if (myrank < nprocs - 1) {
		MPI_Send(au_bd, (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank + 1, 6, MPI_COMM_WORLD);
	}

	// recv buffers
	double **recv_ad_bd;
	double **recv_au_bd;

	if (myrank < nprocs - 1) {
		MPI_Recv(recv_ad_bd, (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank + 1, 0, MPI_COMM_WORLD, &status);
	}
	if (myrank > 0) {
		MPI_Recv(recv_au_bd, (n2 + 2) * (n3 + 2), MPI_DOUBLE, myrank - 1, 6, MPI_COMM_WORLD, &status);
	}


	/*----------------------------------------------------------------------------*/
	/*riesenie sustavy rovnic*/

	pom_local = 0.0;
	it = 0;	

	do {
		it = it + 1;
		for (i = 1; i <= nlocal; i++) {
			for (j = 1; j <= n2; j++) {
				for (k = 1; k <= n3; k++) {
					if ((i + j + k) % 2 == 0) {
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

		for (i = 1; i <= nlocal; i++) {
			for (j = 1; j <= n2; j++) {
				for (k = 1; k <= n3; k++) {
					if ((i + j + k) % 2 == 1) {
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

					res_local = res_local + pom_local * pom_local;

					pom1_local = GM / T_R[i][j][k] - u[i][j][k];
					res3_local = res3_local + pom1_local * pom1_local;
				}
			}
		}

		MPI_Reduce(&res3_local, &res3_global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&res_local, &res_global, 1, MPI_DOUBLE, MPI_SUM, 1, MPI_COMM_WORLD);

		res3_global = sqrt(res3_global / (n1 * n2 * n3));
		printf("\t\t%d\t%.12lf\n", it, res_global);

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

	MPI_Reduce(&pom_local, &pom_global, 1, MPI_DOUBLE, MPI_SUM, 2, MPI_COMM_WORLD);

	sigma = sqrt(pom_global);
	printf("sigma = %.20lf\n", sigma);

	double ***u_global = new double**[n1 + 2];
	double ***res2_global = new double**[n1 + 2];
	double ***TR_global = new double**[n1 + 2];

	for (i = 0; i <= n1; i++) {
		u_global[i] = new double*[n2 + 2];
		res2_global[i] = new double*[n2 + 2];
		TR_global[i] = new double*[n2 + 2];
		for (j = 0; j <= n2; j++) {
			u_global[i][j] = new double[n3 + 2];
			res2_global[i][j] = new double[n3 + 2];
			TR_global[i][j] = new double[n3 + 2];
		}
	}

	MPI_Allgather(u, nlocal * n2 * n3, MPI_DOUBLE, u_global + myrank * (n1 / NPROC + 3), nlocal * n2 * n3, MPI_DOUBLE, MPI_COMM_WORLD);
	MPI_Allgather(res2, nlocal * n2 * n3, MPI_DOUBLE, res2_global + myrank * (n1 / NPROC + 3), nlocal * n2 * n3, MPI_DOUBLE, MPI_COMM_WORLD);
	MPI_Allgather(T_R, nlocal * n2 * n3, MPI_DOUBLE, TR_global + myrank * (n1 / NPROC + 3), nlocal * n2 * n3, MPI_DOUBLE, MPI_COMM_WORLD);

	if (myrank == 0) {
		FILE *fw; // *fr,
		fopen_s(&fw, "stat.txt", "w");

		i = 1;
		for (j = 1; j <= n2; j++) {
			for (k = 1; k <= n3; k++) {
				res2_global[i][j][1] = (GM / TR_global[i][j][k]) - u_global[i][j][k];
				fprintf(fw, "%d %d %d\t%.7lf\t%.9lf\t%.7lf\n", i, j, k, u_global[i][j][j], res2_global[i][j][j], (GM / TR_global[i][j][k]));
			}
		}

		i = n1 / 2;
		for (j = 1; j <= n2; j++) {
			for (k = 1; k <= n3; k++) {
				res2_global[i][j][1] = (GM / TR_global[i][j][k]) - u_global[i][j][k];
				fprintf(fw, "%d %d %d\t%.7lf\t%.9lf\t%.7lf\n", i, j, k, u_global[i][j][j], res2_global[i][j][j], (GM / TR_global[i][j][k]));
			}
		}

		i = n1;
		for (j = 1; j <= n2; j++) {
			for (k = 1; k <= n3; k++) {
				res2_global[i][j][1] = (GM / TR_global[i][j][k]) - u_global[i][j][k];
				fprintf(fw, "%d %d %d\t%.7lf\t%.9lf\t%.7lf\n", i, j, k, u[i][j][j], res2_global[i][j][j], (GM / TR_global[i][j][k]));
			}
		}

		fclose(fw);
	}

	/*
	double **au_bd = au[nlocal];
	double **ad_bd = ad[0];
	*/

	for (i = 0; i <= nlocal; i++) {
		for (j = 0; j <= n2; j++) {
			delete[] B[i][j]; delete[] L[i][j]; delete[] R[i][j];
			delete[] T_B[i][j]; delete[] T_L[i][j]; delete[] T_R[i][j];
			delete[] s[i][j]; delete[] u[i][j]; delete[] deltag[i][j]; delete[] b[i][j]; delete[] res2[i][j];
			delete[] an[i][j]; delete[] as[i][j]; delete[] aw[i][j]; delete[] an[i][j];
			delete[] au[i][j]; delete[] ad[i][j]; delete[] ap[i][j];
			if (i == 0) {
				delete[] recv_L_bdL[j]; delete L_bdL[j]; delete[] recv_L_bdU[j]; delete L_bdU[j];
				delete[] recv_B_bdL[j]; delete B_bdL[j]; delete[] recv_B_bdU[j]; delete B_bdU[j];
				delete[] recv_R_bdL[j]; delete R_bdL[j]; delete[] recv_R_bdU[j]; delete R_bdU[j];
				delete[] recv_TL_bdL[j]; delete TL_bdL[j]; delete[] recv_TL_bdU[j]; delete TL_bdU[j];
				delete[] recv_TB_bdL[j]; delete TB_bdL[j]; delete[] recv_TB_bdU[j]; delete TB_bdU[j];
				delete[] recv_TR_bdL[j]; delete TR_bdL[j]; delete[] recv_TR_bdU[j]; delete TR_bdU[j];
				delete[] recv_ad_bd[j]; delete ad_bd[j]; delete recv_au_bd[j]; delete au_bd[j];
			}
		}
		delete[] B[i]; delete[] L[i]; delete[] R[i];
		delete[] T_B[i]; delete[] T_L[i]; delete[] T_R[i];
		delete[] s[i]; delete[] u[i]; delete[] deltag[i]; delete[] b[i]; delete[] res2[i];
		delete[] an[i]; delete[] as[i]; delete[] aw[i]; delete[] an[i];
		delete[] au[i]; delete[] ad[i]; delete[] ap[i];
	}
	delete[] B; delete[] L; delete[] R;
	delete[] T_B; delete[] T_L; delete[] T_R;
	delete[] s; delete[] u; delete[] deltag; delete[] b; delete[] res2;
	delete[] an; delete[] as; delete[] aw; delete[] an;
	delete[] au; delete[] ad; delete[] ap;

	delete[] recv_L_bdL; delete L_bdL; delete[] recv_L_bdU; delete L_bdU;
	delete[] recv_B_bdL; delete B_bdL; delete[] recv_B_bdU; delete B_bdU;
	delete[] recv_R_bdL; delete R_bdL; delete[] recv_R_bdU; delete R_bdU;
	delete[] recv_TL_bdL; delete TL_bdL; delete[] recv_TL_bdU; delete TL_bdU;
	delete[] recv_TB_bdL; delete TB_bdL; delete[] recv_TB_bdU; delete TB_bdU;
	delete[] recv_TR_bdL; delete TR_bdL; delete[] recv_TR_bdU; delete TR_bdU;
	delete[] recv_ad_bd; delete ad_bd; delete recv_au_bd; delete au_bd;


	MPI_Finalize();

	getchar();

	return 0;
}

