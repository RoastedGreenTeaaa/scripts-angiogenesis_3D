#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <windows.h>

double random_Deci(double a, double b){
    return a + (double)((rand()%((int)(1e+6*(b - a))))/1e+6);
}
int main(int argc, const char * argv[]) {
    long cpu_time;
    double sec;
    //srand( (int)time(NULL) );
    srand(1);
    
    int N = 2500, T = 25000, n_cell = 0, delta = 10, videoDel = 25, cell_plus = 1;
    double initial_cell = 15.0;
    int D = 10, L = 8;
    double f = 0.7;
    double a_r = 1.0/pow(1.0 - f, 2.0/3.0), b_r = pow(1.0 - f, 1.0/3.0), c_r = b_r;
    double a_a = 1.25*a_r, b_a = 1.25*a_r, c_a = b_a;
    double fr = 0.02, fa = 0.001;
    double gamma1 = 0.1, d = 0.002;
    double fp = 0.001, gamma2 = 0.01;
    char directory[64];
    strcpy(directory, "sim1");

    double rx[N], ry[N], rz[N], vx[N], vy[N], vz[N];
    double rx_new[N], ry_new[N], rz_new[N], vx_new[N], vy_new[N], vz_new[N];
    double ex1[N], ex2[N], ex3[N], ey1[N], ey2[N], ey3[N], ez1[N], ez2[N], ez3[N], ex_norm, ey_norm, ez_norm;
    double ex1_, ex2_, ex3_, ey1_, ey2_, ey3_, ez1_, ez2_, ez3_;
    int n_count[N];
    double psi[N], nx[N], ny[N], nz[N], n_norm;
    double RM[N][3][3], RM_new[N][3][3], RM_new2[N][3][3], RM_Save[N][3][3], RM_Save2[N][3][3];
    double Ux_a, Uy_a, Uz_a, x_a, y_a, z_a, A_ellipsoid;
    double Ux_r, Uy_r, Uz_r, x_r, y_r, z_r, R_ellipsoid;
    double Fx[N], Fy[N], Fz[N], c1, c2;
    double x_ij, y_ij, z_ij, norm_ij;
    double RMj_inv[N][9], X, Y, Z;
    double RM_ex[N][9], dt = (2*M_PI)/L, nx_l, ny_l, nz_l, nx_l2, ny_l2, nz_l2;
    double samp_norm, sr1, sr2, sr3;
    double c3[N];
    double v_norm, inner;
    double chi, nx_v, ny_v, nz_v, RM_v[3][3];
    double theta_l, phi_lk;
    double dc, det;
    double R1 = 2.0*a_a + 0.1;
    char fileName1[64];
    int fileCount = 0;

    for (int i=0; i<N; i++) {
        rx[i] = random_Deci(-1.0*initial_cell, initial_cell);
        ry[i] = random_Deci(-1.0*initial_cell, initial_cell);
        rz[i] = random_Deci(-1.0*initial_cell, initial_cell);
        vx[i] = 0.0;
        vy[i] = 0.0;
        vz[i] = 0.0;
        psi[i] = random_Deci(-M_PI, M_PI);        
        nx[i] = random_Deci(-100.0, 100.0);
        ny[i] = random_Deci(-100.0, 100.0);
        nz[i] = random_Deci(-100.0, 100.0);
        n_norm = sqrt(nx[i]*nx[i] + ny[i]*ny[i] + nz[i]*nz[i]);
        nx[i] = nx[i]/(n_norm + 1e-10);
        ny[i] = ny[i]/(n_norm + 1e-10);
        nz[i] = nz[i]/(n_norm + 1e-10);
        ex1[i] = 1;
        ex2[i] = 0;
        ex3[i] = 0;
        ey1[i] = 0;
        ey2[i] = 1;
        ey3[i] = 0;
        ez1[i] = 0;
        ez2[i] = 0;
        ez3[i] = 1;
        
        RM_Save[i][0][0] = 1 + (1 - cos(psi[i]))*(nx[i]*nx[i] - 1);
        RM_Save[i][0][1] = - sin(psi[i])*nz[i] + (1 - cos(psi[i]))*nx[i]*ny[i];
        RM_Save[i][0][2] = sin(psi[i])*ny[i] + (1 - cos(psi[i]))*nx[i]*nz[i];                    
        RM_Save[i][1][0] = sin(psi[i])*nz[i] + (1 - cos(psi[i]))*nx[i]*ny[i];
        RM_Save[i][1][1] = 1 + (1 - cos(psi[i]))*(ny[i]*ny[i] - 1);
        RM_Save[i][1][2] = - sin(psi[i])*nx[i] + (1 - cos(psi[i]))*ny[i]*nz[i];                    
        RM_Save[i][2][0] = - sin(psi[i])*ny[i] + (1 - cos(psi[i]))*nx[i]*nz[i];
        RM_Save[i][2][1] = sin(psi[i])*nx[i] + (1 - cos(psi[i]))*ny[i]*nz[i];
        RM_Save[i][2][2] = 1 + (1 - cos(psi[i]))*(nz[i]*nz[i] - 1);
        
        ex1_ = RM_Save[i][0][0]*ex1[i] + RM_Save[i][0][1]*ex2[i] + RM_Save[i][0][2]*ex3[i];
        ex2_ = RM_Save[i][1][0]*ex1[i] + RM_Save[i][1][1]*ex2[i] + RM_Save[i][1][2]*ex3[i];
        ex3_ = RM_Save[i][2][0]*ex1[i] + RM_Save[i][2][1]*ex2[i] + RM_Save[i][2][2]*ex3[i];
        ey1_ = RM_Save[i][0][0]*ey1[i] + RM_Save[i][0][1]*ey2[i] + RM_Save[i][0][2]*ey3[i];
        ey2_ = RM_Save[i][1][0]*ey1[i] + RM_Save[i][1][1]*ey2[i] + RM_Save[i][1][2]*ey3[i];
        ey3_ = RM_Save[i][2][0]*ey1[i] + RM_Save[i][2][1]*ey2[i] + RM_Save[i][2][2]*ey3[i];
        ez1_ = RM_Save[i][0][0]*ez1[i] + RM_Save[i][0][1]*ez2[i] + RM_Save[i][0][2]*ez3[i];
        ez2_ = RM_Save[i][1][0]*ez1[i] + RM_Save[i][1][1]*ez2[i] + RM_Save[i][1][2]*ez3[i];
        ez3_ = RM_Save[i][2][0]*ez1[i] + RM_Save[i][2][1]*ez2[i] + RM_Save[i][2][2]*ez3[i];
        ex1[i] = ex1_;
        ex2[i] = ex2_;
        ex3[i] = ex3_;
        ey1[i] = ey1_;
        ey2[i] = ey2_;
        ey3[i] = ey3_;
        ez1[i] = ez1_;
        ez2[i] = ez2_;
        ez3[i] = ez3_;
        ex_norm = sqrt(ex1[i]*ex1[i] + ex2[i]*ex2[i] + ex3[i]*ex3[i]);
        ey_norm = sqrt(ey1[i]*ey1[i] + ey2[i]*ey2[i] + ey3[i]*ey3[i]);
        ez_norm = sqrt(ez1[i]*ez1[i] + ez2[i]*ez2[i] + ez3[i]*ez3[i]);
        ex1[i] = ex1[i]/ex_norm;
        ex2[i] = ex2[i]/ex_norm;
        ex3[i] = ex3[i]/ex_norm;
        ey1[i] = ey1[i]/ey_norm;
        ey2[i] = ey2[i]/ey_norm;
        ey3[i] = ey3[i]/ey_norm;
        ez1[i] = ez1[i]/ez_norm;
        ez2[i] = ez2[i]/ez_norm;
        ez3[i] = ez3[i]/ez_norm;
    }
    
    if(CreateDirectory(directory, NULL)==0) {
    }
    sprintf(fileName1, "%s\\%08d.txt", directory, fileCount);
    FILE* file = fopen(fileName1, "w");
    
    for (int i=0; i<N; i++) {
        fprintf(file, "%f   ", rx[i]);
        fprintf(file, "%f   ", ry[i]);
        fprintf(file, "%f   ", rz[i]);
        fprintf(file, "%f   ", ex1[i]);
        fprintf(file, "%f   ", ex2[i]);
        fprintf(file, "%f   ", ex3[i]);
        fprintf(file, "%f   ", ey1[i]);
        fprintf(file, "%f   ", ey2[i]);
        fprintf(file, "%f   ", ey3[i]);
        fprintf(file, "%f   ", ez1[i]);
        fprintf(file, "%f   ", ez2[i]);
        fprintf(file, "%f   ", ez3[i]);
        fprintf(file, "%f   ", a_r);
        fprintf(file, "%f   ", b_r);
        fprintf(file, "%f   \n", c_r);
    }
    fclose(file);
    fileCount = fileCount + videoDel;
    
    for (int t=0; t<T; t++) {                
        if (((t+1)%delta==0)&&(n_cell<N)){
            for (int h=0; h<cell_plus; h++) {
                rx[n_cell+h] = random_Deci(-initial_cell, initial_cell);
                ry[n_cell+h] = random_Deci(-initial_cell, initial_cell);
                rz[n_cell+h] = random_Deci(-initial_cell, initial_cell);
                vx[n_cell+h] = 0.0;
                vy[n_cell+h] = 0.0;
                vz[n_cell+h] = 0.0;
                psi[n_cell+h] = random_Deci(-0.75*M_PI, 0.75*M_PI);                
                nx[n_cell+h] = random_Deci(-100.0, 100.0);
                ny[n_cell+h] = random_Deci(-100.0, 100.0);
                nz[n_cell+h] = random_Deci(-100.0, 100.0);
                n_norm = sqrt(nx[n_cell+h]*nx[n_cell+h] + ny[n_cell+h]*ny[n_cell+h] + nz[n_cell+h]*nz[n_cell+h]);
                nx[n_cell+h] = nx[n_cell+h]/(n_norm + 1e-10);
                ny[n_cell+h] = ny[n_cell+h]/(n_norm + 1e-10);
                nz[n_cell+h] = nz[n_cell+h]/(n_norm + 1e-10);
                ex1[n_cell+h] = 1;
                ex2[n_cell+h] = 0;
                ex3[n_cell+h] = 0;
                ey1[n_cell+h] = 0;
                ey2[n_cell+h] = 1;
                ey3[n_cell+h] = 0;
                ez1[n_cell+h] = 0;
                ez2[n_cell+h] = 0;
                ez3[n_cell+h] = 1;
                RM_Save[n_cell+h][0][0] = 1 + (1 - cos(psi[n_cell+h]))*(nx[n_cell+h]*nx[n_cell+h] - 1);
                RM_Save[n_cell+h][0][1] = - sin(psi[n_cell+h])*nz[n_cell+h] + (1 - cos(psi[n_cell+h]))*nx[n_cell+h]*ny[n_cell+h];
                RM_Save[n_cell+h][0][2] = sin(psi[n_cell+h])*ny[n_cell+h] + (1 - cos(psi[n_cell+h]))*nx[n_cell+h]*nz[n_cell+h];                            
                RM_Save[n_cell+h][1][0] = sin(psi[n_cell+h])*nz[n_cell+h] + (1 - cos(psi[n_cell+h]))*nx[n_cell+h]*ny[n_cell+h];
                RM_Save[n_cell+h][1][1] = 1 + (1 - cos(psi[n_cell+h]))*(ny[n_cell+h]*ny[n_cell+h] - 1);
                RM_Save[n_cell+h][1][2] = - sin(psi[n_cell+h])*nx[n_cell+h] + (1 - cos(psi[n_cell+h]))*ny[n_cell+h]*nz[n_cell+h];                            
                RM_Save[n_cell+h][2][0] = - sin(psi[n_cell+h])*ny[n_cell+h] + (1 - cos(psi[n_cell+h]))*nx[n_cell+h]*nz[n_cell+h];
                RM_Save[n_cell+h][2][1] = sin(psi[n_cell+h])*nx[n_cell+h] + (1 - cos(psi[n_cell+h]))*ny[n_cell+h]*nz[n_cell+h];
                RM_Save[n_cell+h][2][2] = 1 + (1 - cos(psi[n_cell+h]))*(nz[n_cell+h]*nz[n_cell+h] - 1);
                ex1_ = RM_Save[n_cell+h][0][0]*ex1[n_cell+h];
                ex2_ = RM_Save[n_cell+h][1][0]*ex1[n_cell+h];
                ex3_ = RM_Save[n_cell+h][2][0]*ex1[n_cell+h];
                ey1_ = RM_Save[n_cell+h][0][1]*ey2[n_cell+h];
                ey2_ = RM_Save[n_cell+h][1][1]*ey2[n_cell+h];
                ey3_ = RM_Save[n_cell+h][2][1]*ey2[n_cell+h];
                ez1_ = RM_Save[n_cell+h][0][2]*ez3[n_cell+h];
                ez2_ = RM_Save[n_cell+h][1][2]*ez3[n_cell+h];
                ez3_ = RM_Save[n_cell+h][2][2]*ez3[n_cell+h];
                ex1[n_cell+h] = ex1_;
                ex2[n_cell+h] = ex2_;
                ex3[n_cell+h] = ex3_;
                ey1[n_cell+h] = ey1_;
                ey2[n_cell+h] = ey2_;
                ey3[n_cell+h] = ey3_;
                ez1[n_cell+h] = ez1_;
                ez2[n_cell+h] = ez2_;
                ez3[n_cell+h] = ez3_;
                ex_norm = sqrt(ex1[n_cell+h]*ex1[n_cell+h] + ex2[n_cell+h]*ex2[n_cell+h] + ex3[n_cell+h]*ex3[n_cell+h]);
                ey_norm = sqrt(ey1[n_cell+h]*ey1[n_cell+h] + ey2[n_cell+h]*ey2[n_cell+h] + ey3[n_cell+h]*ey3[n_cell+h]);
                ez_norm = sqrt(ez1[n_cell+h]*ez1[n_cell+h] + ez2[n_cell+h]*ez2[n_cell+h] + ez3[n_cell+h]*ez3[n_cell+h]);
                ex1[n_cell+h] = ex1[n_cell+h]/ex_norm;
                ex2[n_cell+h] = ex2[n_cell+h]/ex_norm;
                ex3[n_cell+h] = ex3[n_cell+h]/ex_norm;
                ey1[n_cell+h] = ey1[n_cell+h]/ey_norm;
                ey2[n_cell+h] = ey2[n_cell+h]/ey_norm;
                ey3[n_cell+h] = ey3[n_cell+h]/ey_norm;
                ez1[n_cell+h] = ez1[n_cell+h]/ez_norm;
                ez2[n_cell+h] = ez2[n_cell+h]/ez_norm;
                ez3[n_cell+h] = ez3[n_cell+h]/ez_norm;
            }
            n_cell = n_cell + cell_plus;
        }
        for (int i=0; i<n_cell; i++) {
            psi[i] = 0;
            Fx[i] = 0.0;
            Fy[i] = 0.0;
            Fz[i] = 0.0;
            nx[i] = 0.0;
            ny[i] = 0.0;
            nz[i] = 0.0;
            n_count[i] = 0;
            RM_new[i][0][0] = 1.0;
            RM_new[i][0][1] = 0.0;
            RM_new[i][0][2] = 0.0;
            RM_new[i][1][0] = 0.0;
            RM_new[i][1][1] = 1.0;
            RM_new[i][1][2] = 0.0;
            RM_new[i][2][0] = 0.0;
            RM_new[i][2][1] = 0.0;
            RM_new[i][2][2] = 1.0;
            RM_ex[i][0] = 1.0 + (1.0 - cos(dt))*(ex1[i]*ex1[i] - 1.0);
            RM_ex[i][1] = - sin(dt)*ex3[i] + (1.0 - cos(dt))*ex1[i]*ex2[i];
            RM_ex[i][2] = sin(dt)*ex2[i] + (1.0 - cos(dt))*ex1[i]*ex3[i];
            RM_ex[i][3] = sin(dt)*ex3[i] + (1.0 - cos(dt))*ex1[i]*ex2[i];
            RM_ex[i][4] = 1.0 + (1.0 - cos(dt))*(ex2[i]*ex2[i] - 1.0);
            RM_ex[i][5] = - sin(dt)*ex1[i] + (1.0 - cos(dt))*ex2[i]*ex3[i];
            RM_ex[i][6] = - sin(dt)*ex2[i] + (1.0 - cos(dt))*ex1[i]*ex3[i];
            RM_ex[i][7] = sin(dt)*ex1[i] + (1.0 - cos(dt))*ex2[i]*ex3[i];
            RM_ex[i][8] = 1.0 + (1.0 - cos(dt))*(ex3[i]*ex3[i] - 1.0);
            RMj_inv[i][0] = RM_Save[i][1][1]*RM_Save[i][2][2] - RM_Save[i][1][2]*RM_Save[i][2][1];
            RMj_inv[i][1] = RM_Save[i][0][2]*RM_Save[i][2][1] - RM_Save[i][0][1]*RM_Save[i][2][2];
            RMj_inv[i][2] = RM_Save[i][0][1]*RM_Save[i][1][2] - RM_Save[i][0][2]*RM_Save[i][1][1];
            RMj_inv[i][3] = RM_Save[i][1][2]*RM_Save[i][2][0] - RM_Save[i][1][0]*RM_Save[i][2][2];
            RMj_inv[i][4] = RM_Save[i][0][0]*RM_Save[i][2][2] - RM_Save[i][0][2]*RM_Save[i][2][0];
            RMj_inv[i][5] = RM_Save[i][0][2]*RM_Save[i][1][0] - RM_Save[i][0][0]*RM_Save[i][1][2];
            RMj_inv[i][6] = RM_Save[i][1][0]*RM_Save[i][2][1] - RM_Save[i][1][1]*RM_Save[i][2][0];
            RMj_inv[i][7] = RM_Save[i][0][1]*RM_Save[i][2][0] - RM_Save[i][0][0]*RM_Save[i][2][1];
            RMj_inv[i][8] = RM_Save[i][0][0]*RM_Save[i][1][1] - RM_Save[i][0][1]*RM_Save[i][1][0];
            c3[i] = 0.0;
        }
        for (int i=0; i<n_cell; i++) {
            for (int j=i+1; j<n_cell; j++) {
                x_ij = rx[j] - rx[i];
                y_ij = ry[j] - ry[i];
                z_ij = rz[j] - rz[i];
                norm_ij = sqrt(x_ij*x_ij + y_ij*y_ij + z_ij*z_ij);
                if (norm_ij >= R1) {
                    continue;
                }
                x_ij = x_ij/(norm_ij + 1e-10);
                y_ij = y_ij/(norm_ij + 1e-10);
                z_ij = z_ij/(norm_ij + 1e-10);
                c1 = 0;
                c2 = 0;
                for (int l=0; l<L; l++) {
                    theta_l = (2.0*M_PI*l)/L;
                    for (int k=0; k<=D; k++) {
                        phi_lk = (M_PI*k)/D;
                        Ux_a = a_a*cos(phi_lk);
                        Uy_a = b_a*sin(phi_lk)*cos(theta_l);
                        Uz_a = c_a*sin(phi_lk)*sin(theta_l);
                        x_a = rx[i] + RM_Save[i][0][0]*Ux_a + RM_Save[i][0][1]*Uy_a + RM_Save[i][0][2]*Uz_a;
                        y_a = ry[i] + RM_Save[i][1][0]*Ux_a + RM_Save[i][1][1]*Uy_a + RM_Save[i][1][2]*Uz_a;
                        z_a = rz[i] + RM_Save[i][2][0]*Ux_a + RM_Save[i][2][1]*Uy_a + RM_Save[i][2][2]*Uz_a;
                        sr1 = rx[j] - x_a;
                        sr2 = ry[j] - y_a;
                        sr3 = rz[j] - z_a;
                        samp_norm = sqrt(sr1*sr1 + sr2*sr2 + sr3*sr3);
                        if (samp_norm > (a_a + 0.25*a_a)) {
                            continue;
                        }
                        X = RMj_inv[j][0]*(x_a-rx[j]) + RMj_inv[j][1]*(y_a-ry[j]) + RMj_inv[j][2]*(z_a-rz[j]);
                        Y = RMj_inv[j][3]*(x_a-rx[j]) + RMj_inv[j][4]*(y_a-ry[j]) + RMj_inv[j][5]*(z_a-rz[j]);
                        Z = RMj_inv[j][6]*(x_a-rx[j]) + RMj_inv[j][7]*(y_a-ry[j]) + RMj_inv[j][8]*(z_a-rz[j]);                            
                        A_ellipsoid = (X/a_a)*(X/a_a) + (Y/b_a)*(Y/b_a) + (Z/c_a)*(Z/c_a);
                        if (A_ellipsoid < 1.0){
                            c1 = fa;
                            c3[i] = d;
                            c3[j] = d;
                            l = L;
                            break;
                        }
                    }
                }
                for (int l=0; l<L; l++) {
                    theta_l = (2.0*M_PI*l)/L;
                    if (l==0) {
                        nx_l = ez1[i];
                        ny_l = ez2[i];
                        nz_l = ez3[i];
                    }else{
                        nx_l2 = RM_ex[i][0]*nx_l + RM_ex[i][1]*ny_l + RM_ex[i][2]*nz_l;
                        ny_l2 = RM_ex[i][3]*nx_l + RM_ex[i][4]*ny_l + RM_ex[i][5]*nz_l;
                        nz_l2 = RM_ex[i][6]*nx_l + RM_ex[i][7]*ny_l + RM_ex[i][8]*nz_l;
                        nx_l = nx_l2;
                        ny_l = ny_l2;
                        nz_l = nz_l2;
                    }                        
                    for (int k=0; k<=D; k++) {
                        phi_lk = (M_PI*k)/D;
                        Ux_r = a_r*cos(phi_lk);
                        Uy_r = b_r*sin(phi_lk)*cos(theta_l);
                        Uz_r = c_r*sin(phi_lk)*sin(theta_l);
                        x_r = rx[i] + RM_Save[i][0][0]*Ux_r + RM_Save[i][0][1]*Uy_r + RM_Save[i][0][2]*Uz_r;
                        y_r = ry[i] + RM_Save[i][1][0]*Ux_r + RM_Save[i][1][1]*Uy_r + RM_Save[i][1][2]*Uz_r;
                        z_r = rz[i] + RM_Save[i][2][0]*Ux_r + RM_Save[i][2][1]*Uy_r + RM_Save[i][2][2]*Uz_r;
                        sr1 = rx[j] - x_r;
                        sr2 = ry[j] - y_r;
                        sr3 = rz[j] - z_r;
                        samp_norm = sqrt(sr1*sr1 + sr2*sr2 + sr3*sr3);
                        if (samp_norm > (a_r + 0.25*a_r)) {
                            continue;
                        }
                        X = RMj_inv[j][0]*(x_r-rx[j]) + RMj_inv[j][1]*(y_r-ry[j]) + RMj_inv[j][2]*(z_r-rz[j]);
                        Y = RMj_inv[j][3]*(x_r-rx[j]) + RMj_inv[j][4]*(y_r-ry[j]) + RMj_inv[j][5]*(z_r-rz[j]);
                        Z = RMj_inv[j][6]*(x_r-rx[j]) + RMj_inv[j][7]*(y_r-ry[j]) + RMj_inv[j][8]*(z_r-rz[j]);
                        R_ellipsoid = (X/a_r)*(X/a_r) + (Y/b_r)*(Y/b_r) + (Z/c_r)*(Z/c_r);
                        if (R_ellipsoid < 1.0){
                            c2 = fr;
                            psi[i] = - fp*sin(2.0*phi_lk);
                            nx[i] = nx_l;
                            ny[i] = ny_l;
                            nz[i] = nz_l;
                            RM[i][0][0] = 1.0 + (1.0 - cos(psi[i]))*(nx[i]*nx[i] - 1.0);
                            RM[i][0][1] = - sin(psi[i])*nz[i] + (1.0 - cos(psi[i]))*nx[i]*ny[i];
                            RM[i][0][2] = sin(psi[i])*ny[i] + (1.0 - cos(psi[i]))*nx[i]*nz[i];                                            
                            RM[i][1][0] = sin(psi[i])*nz[i] + (1.0 - cos(psi[i]))*nx[i]*ny[i];
                            RM[i][1][1] = 1.0 + (1.0 - cos(psi[i]))*(ny[i]*ny[i] - 1.0);
                            RM[i][1][2] = - sin(psi[i])*nx[i] + (1.0 - cos(psi[i]))*ny[i]*nz[i];                                            
                            RM[i][2][0] = - sin(psi[i])*ny[i] + (1.0 - cos(psi[i]))*nx[i]*nz[i];
                            RM[i][2][1] = sin(psi[i])*nx[i] + (1.0 - cos(psi[i]))*ny[i]*nz[i];
                            RM[i][2][2] = 1.0 + (1.0 - cos(psi[i]))*(nz[i]*nz[i] - 1.0);
                            for (int p=0; p<3; p++) {
                                for (int q=0; q<3; q++) {
                                    dc = 0.0;
                                    for(int o=0; o<3; o++){
                                        dc = dc + RM_new[i][p][o] * RM[i][o][q];
                                    }
                                    RM_new2[i][p][q]= dc;
                                }
                            }                                
                            for (int p=0; p<3; p++) {
                                for (int q=0; q<3; q++) {
                                    RM_new[i][p][q] = RM_new2[i][p][q];
                                }
                            }                                
                            n_count[i] = n_count[i] + 1;
                        }
                    }
                }
                Fx[i] = Fx[i] + (c1 - c2)*x_ij;
                Fy[i] = Fy[i] + (c1 - c2)*y_ij;
                Fz[i] = Fz[i] + (c1 - c2)*z_ij;
                Fx[j] = Fx[j] + (c2 - c1)*x_ij;
                Fy[j] = Fy[j] + (c2 - c1)*y_ij;
                Fz[j] = Fz[j] + (c2 - c1)*z_ij;
                for (int l=0; l<L; l++) {
                    theta_l = (2.0*M_PI*l)/L;
                    if (l==0) {
                        nx_l = ez1[j];
                        ny_l = ez2[j];
                        nz_l = ez3[j];
                    }else{
                        nx_l2 = RM_ex[j][0]*nx_l + RM_ex[j][1]*ny_l + RM_ex[j][2]*nz_l;
                        ny_l2 = RM_ex[j][3]*nx_l + RM_ex[j][4]*ny_l + RM_ex[j][5]*nz_l;
                        nz_l2 = RM_ex[j][6]*nx_l + RM_ex[j][7]*ny_l + RM_ex[j][8]*nz_l;
                        nx_l = nx_l2;
                        ny_l = ny_l2;
                        nz_l = nz_l2;
                    }                        
                    for (int k=0; k<=D; k++) {
                        phi_lk = (M_PI*k)/D;
                        Ux_r = a_r*cos(phi_lk);
                        Uy_r = b_r*sin(phi_lk)*cos(theta_l);
                        Uz_r = c_r*sin(phi_lk)*sin(theta_l);
                        x_r = rx[j] + RM_Save[j][0][0]*Ux_r + RM_Save[j][0][1]*Uy_r + RM_Save[j][0][2]*Uz_r;
                        y_r = ry[j] + RM_Save[j][1][0]*Ux_r + RM_Save[j][1][1]*Uy_r + RM_Save[j][1][2]*Uz_r;
                        z_r = rz[j] + RM_Save[j][2][0]*Ux_r + RM_Save[j][2][1]*Uy_r + RM_Save[j][2][2]*Uz_r;
                        sr1 = rx[i] - x_r;
                        sr2 = ry[i] - y_r;
                        sr3 = rz[i] - z_r;
                        samp_norm = sqrt(sr1*sr1 + sr2*sr2 + sr3*sr3);
                        if (samp_norm > (a_r + 0.25*a_r)) {
                            continue;
                        }
                        X = RMj_inv[i][0]*(x_r-rx[i]) + RMj_inv[i][1]*(y_r-ry[i]) + RMj_inv[i][2]*(z_r-rz[i]);
                        Y = RMj_inv[i][3]*(x_r-rx[i]) + RMj_inv[i][4]*(y_r-ry[i]) + RMj_inv[i][5]*(z_r-rz[i]);
                        Z = RMj_inv[i][6]*(x_r-rx[i]) + RMj_inv[i][7]*(y_r-ry[i]) + RMj_inv[i][8]*(z_r-rz[i]);
                        R_ellipsoid = (X/a_r)*(X/a_r) + (Y/b_r)*(Y/b_r) + (Z/c_r)*(Z/c_r);
                        if (R_ellipsoid < 1.0){
                            psi[j] = - fp*sin(2.0*phi_lk);
                            nx[j] = nx_l;
                            ny[j] = ny_l;
                            nz[j] = nz_l;
                            RM[j][0][0] = 1.0 + (1.0 - cos(psi[j]))*(nx[j]*nx[j] - 1.0);
                            RM[j][0][1] = - sin(psi[j])*nz[j] + (1.0 - cos(psi[j]))*nx[j]*ny[j];
                            RM[j][0][2] = sin(psi[j])*ny[j] + (1.0 - cos(psi[j]))*nx[j]*nz[j];                                            
                            RM[j][1][0] = sin(psi[j])*nz[j] + (1.0 - cos(psi[j]))*nx[j]*ny[j];
                            RM[j][1][1] = 1.0 + (1.0 - cos(psi[j]))*(ny[j]*ny[j] - 1.0);
                            RM[j][1][2] = - sin(psi[j])*nx[j] + (1.0 - cos(psi[j]))*ny[j]*nz[j];                                            
                            RM[j][2][0] = - sin(psi[j])*ny[j] + (1.0 - cos(psi[j]))*nx[j]*nz[j];
                            RM[j][2][1] = sin(psi[j])*nx[j] + (1.0 - cos(psi[j]))*ny[j]*nz[j];
                            RM[j][2][2] = 1.0 + (1.0 - cos(psi[j]))*(nz[j]*nz[j] - 1.0);
                            for (int p=0; p<3; p++) {
                                for (int q=0; q<3; q++) {
                                    dc = 0.0;
                                    for(int o=0; o<3; o++){
                                        dc = dc + RM_new[j][p][o] * RM[j][o][q];
                                    }
                                    RM_new2[j][p][q]= dc;
                                }
                            }                                
                            for (int p=0; p<3; p++) {
                                for (int q=0; q<3; q++) {
                                    RM_new[j][p][q] = RM_new2[j][p][q];
                                }
                            }
                            n_count[j] = n_count[j] + 1;
                        }
                    }
                }
                
            }
        }
for (int i=0; i<n_cell; i++) {
            rx_new[i] = rx[i] + vx[i];
            ry_new[i] = ry[i] + vy[i];
            rz_new[i] = rz[i] + vz[i];
            if((vx[i]>-1e-8)&&(vx[i]<1e-8)&&(vy[i]>-1e-8)&&(vy[i]<1e-8)&&(vz[i]>-1e-8)&&(vz[i]<1e-8)) {
                vx_new[i] = ((1.0 - gamma1)*vx[i]) + Fx[i];
                vy_new[i] = ((1.0 - gamma1)*vy[i]) + Fy[i];
                vz_new[i] = ((1.0 - gamma1)*vz[i]) + Fz[i];
            }else{
                v_norm = sqrt(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
                vx_new[i] = ((1.0 - gamma1)*vx[i]) + Fx[i] + ((c3[i]*vx[i])/(v_norm + 1e-10));
                vy_new[i] = ((1.0 - gamma1)*vy[i]) + Fy[i] + ((c3[i]*vy[i])/(v_norm + 1e-10));
                vz_new[i] = ((1.0 - gamma1)*vz[i]) + Fz[i] + ((c3[i]*vz[i])/(v_norm + 1e-10));
            }
            if (n_count[i] == 0) {
                RM_new[i][0][0] = 1.0;
                RM_new[i][0][1] = 0.0;
                RM_new[i][0][2] = 0.0;
                RM_new[i][1][0] = 0.0;
                RM_new[i][1][1] = 1.0;
                RM_new[i][1][2] = 0.0;
                RM_new[i][2][0] = 0.0;
                RM_new[i][2][1] = 0.0;
                RM_new[i][2][2] = 1.0;
            }
            if ((-vx_new[i]<1e-5)&&(vx_new[i]<1e-5)&&(-vy_new[i]<1e-5)&&(vy_new[i]<1e-5)&&(-vz_new[i]<1e-5)&&(vz_new[i]<1e-5)) {
                RM_v[0][0] = 1.0;
                RM_v[0][1] = 0.0;
                RM_v[0][2] = 0.0;
                RM_v[1][0] = 0.0;
                RM_v[1][1] = 1.0;
                RM_v[1][2] = 0.0;
                RM_v[2][0] = 0.0;
                RM_v[2][1] = 0.0;
                RM_v[2][2] = 1.0;
            }else{
                inner = ex1[i]*vx_new[i] + ex2[i]*vy_new[i] + ex3[i]*vz_new[i];
                v_norm = sqrt(vx_new[i]*vx_new[i] + vy_new[i]*vy_new[i] + vz_new[i]*vz_new[i]);
                inner = inner/(v_norm + 1e-10);
                chi = gamma2*sin(2.0*acos(inner));                
                nx_v = ex2[i]*vz_new[i] - ex3[i]*vy_new[i];
                ny_v = ex3[i]*vx_new[i] - ex1[i]*vz_new[i];
                nz_v = ex1[i]*vy_new[i] - ex2[i]*vx_new[i];
                v_norm = sqrt(nx_v*nx_v + ny_v*ny_v + nz_v*nz_v);
                nx_v = nx_v/(v_norm + 1e-10);
                ny_v = ny_v/(v_norm + 1e-10);
                nz_v = nz_v/(v_norm + 1e-10);                
                RM_v[0][0] = 1.0 + (1.0 - cos(chi))*(nx_v*nx_v - 1.0);
                RM_v[0][1] = - sin(chi)*nz_v + (1.0 - cos(chi))*nx_v*ny_v;
                RM_v[0][2] = sin(chi)*ny_v + (1.0 - cos(chi))*nx_v*nz_v;
                RM_v[1][0] = sin(chi)*nz_v + (1.0 - cos(chi))*nx_v*ny_v;
                RM_v[1][1] = 1.0 + (1.0 - cos(chi))*(ny_v*ny_v - 1.0);
                RM_v[1][2] = - sin(chi)*nx_v + (1.0 - cos(chi))*ny_v*nz_v;                           
                RM_v[2][0] = - sin(chi)*ny_v + (1.0 - cos(chi))*nx_v*nz_v;
                RM_v[2][1] = sin(chi)*nx_v + (1.0 - cos(chi))*ny_v*nz_v;
                RM_v[2][2] = 1.0 + (1.0 - cos(chi))*(nz_v*nz_v - 1.0);
                for (int p=0; p<3; p++) {
                    for (int q=0; q<3; q++) {
                        dc = 0.0;
                        for(int o=0; o<3; o++){
                            dc = dc + RM_v[p][o] * RM_new[i][o][q];
                        }
                        RM_new2[i][p][q]= dc;
                    }
                }                
                for (int p=0; p<3; p++) {
                    for (int q=0; q<3; q++) {
                        RM_new[i][p][q] = RM_new2[i][p][q];
                    }
                }                
            }
        }
        for (int i=0; i<n_cell; i++) {
            ez1_ = RM_new[i][0][0]*ez1[i] + RM_new[i][0][1]*ez2[i] + RM_new[i][0][2]*ez3[i];
            ez2_ = RM_new[i][1][0]*ez1[i] + RM_new[i][1][1]*ez2[i] + RM_new[i][1][2]*ez3[i];
            ez3_ = RM_new[i][2][0]*ez1[i] + RM_new[i][2][1]*ez2[i] + RM_new[i][2][2]*ez3[i];
            ey1_ = RM_new[i][0][0]*ey1[i] + RM_new[i][0][1]*ey2[i] + RM_new[i][0][2]*ey3[i];
            ey2_ = RM_new[i][1][0]*ey1[i] + RM_new[i][1][1]*ey2[i] + RM_new[i][1][2]*ey3[i];
            ey3_ = RM_new[i][2][0]*ey1[i] + RM_new[i][2][1]*ey2[i] + RM_new[i][2][2]*ey3[i];
            ex1_ = RM_new[i][0][0]*ex1[i] + RM_new[i][0][1]*ex2[i] + RM_new[i][0][2]*ex3[i];
            ex2_ = RM_new[i][1][0]*ex1[i] + RM_new[i][1][1]*ex2[i] + RM_new[i][1][2]*ex3[i];
            ex3_ = RM_new[i][2][0]*ex1[i] + RM_new[i][2][1]*ex2[i] + RM_new[i][2][2]*ex3[i];
            ex1[i] = ex1_;
            ex2[i] = ex2_;
            ex3[i] = ex3_;
            ey1[i] = ey1_;
            ey2[i] = ey2_;
            ey3[i] = ey3_;
            ez1[i] = ez1_;
            ez2[i] = ez2_;
            ez3[i] = ez3_;
            ez_norm = sqrt(ez1[i]*ez1[i] + ez2[i]*ez2[i] + ez3[i]*ez3[i]);
            ey_norm = sqrt(ey1[i]*ey1[i] + ey2[i]*ey2[i] + ey3[i]*ey3[i]);
            ex_norm = sqrt(ex1[i]*ex1[i] + ex2[i]*ex2[i] + ex3[i]*ex3[i]);
            ez1[i] = ez1[i]/(ez_norm + 1e-10);
            ez2[i] = ez2[i]/(ez_norm + 1e-10);
            ez3[i] = ez3[i]/(ez_norm + 1e-10);
            ey1[i] = ey1[i]/(ey_norm + 1e-10);
            ey2[i] = ey2[i]/(ey_norm + 1e-10);
            ey3[i] = ey3[i]/(ey_norm + 1e-10);
            ex1[i] = ex1[i]/(ex_norm + 1e-10);
            ex2[i] = ex2[i]/(ex_norm + 1e-10);
            ex3[i] = ex3[i]/(ex_norm + 1e-10);            
        }
        for (int i=0; i<n_cell; i++) {
            for (int p=0; p<3; p++) {
                for (int q=0; q<3; q++) {
                    dc = 0.0;
                    for(int o=0; o<3; o++){
                        dc = dc + RM_new[i][p][o] * RM_Save[i][o][q];
                    }
                    RM_Save2[i][p][q]= dc;
                }
            }            
            for (int p=0; p<3; p++) {
                for (int q=0; q<3; q++) {
                    RM_Save[i][p][q] = RM_Save2[i][p][q];
                }
            }
        }
        if ((t+1)%videoDel==0) {
            sprintf(fileName1, "%s\\%08d.txt", directory, fileCount);
            FILE* file = fopen(fileName1, "w");
            for (int i=0; i<n_cell; i++) {
                fprintf(file, "%f   ", rx_new[i]);
                fprintf(file, "%f   ", ry_new[i]);
                fprintf(file, "%f   ", rz_new[i]);
                fprintf(file, "%f   ", ex1[i]);
                fprintf(file, "%f   ", ex2[i]);
                fprintf(file, "%f   ", ex3[i]);
                fprintf(file, "%f   ", ey1[i]);
                fprintf(file, "%f   ", ey2[i]);
                fprintf(file, "%f   ", ey3[i]);
                fprintf(file, "%f   ", ez1[i]);
                fprintf(file, "%f   ", ez2[i]);
                fprintf(file, "%f   ", ez3[i]);
                fprintf(file, "%f   ", a_r);
                fprintf(file, "%f   ", b_r);
                fprintf(file, "%f   \n", c_r);
            }
        }        
        if ((t+1)%videoDel==0) {
            fclose(file);
            fileCount = fileCount + videoDel;
        }
        for (int i=0; i<n_cell; i++) {
            rx[i] = rx_new[i];
            ry[i] = ry_new[i];
            rz[i] = rz_new[i];
            vx[i] = vx_new[i];
            vy[i] = vy_new[i];
            vz[i] = vz_new[i];
        }        
    }
    cpu_time = clock();
    sec = (double)cpu_time / CLOCKS_PER_SEC;
    printf("%f\n", sec);    
    return 0;
}
