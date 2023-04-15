#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "Conjugate_Gradient.h"

int main()
{
    char junk[20];
    int i, j, k, p, q, temp;
    FILE *inputPtr;
    if ((inputPtr = fopen("FEMInput.txt", "r")) == NULL)
        printf("Error: no such file!");

    int nodenum;
    fscanf(inputPtr, "%s%d", &junk, &nodenum);
    fscanf(inputPtr, "%s%s%s%s", &junk, &junk, &junk, &junk);

    double node[nodenum][3];
    for (i = 0; i < nodenum; i++)
    {
        fscanf(inputPtr, "%d%lf%lf%lf", &junk, &node[i][0], &node[i][1], &node[i][2]);
    }

    int elementnum;

    fscanf(inputPtr, "%s%d%s%s%s%s%s%s%s%s%s%s%s%s", &junk, &elementnum, &junk, &junk, &junk, &junk, &junk, &junk, &junk, &junk, &junk, &junk, &junk, &junk);
    double E[elementnum];
    double v[elementnum];
    double rho[elementnum];
    int NodeX[elementnum][8];

    for (i = 0; i < elementnum; i++)
    {
        fscanf(inputPtr, "%d%lf%lf%lf", &junk, &E[i], &v[i], &rho[i]);
        for(j = 0; j < 8; j++)
        {
            fscanf(inputPtr, "%d", &NodeX[i][j]);
        }
        rho[i] = rho[i];
    }

    Matrix pos,vel,acc;
    pos = make_mat(nodenum*3,1);
    vel = make_mat(nodenum*3,1);
    acc = make_mat(nodenum*3,1);
    fscanf(inputPtr,"%s %s %s %s %s", &junk,&junk,&junk,&junk,&junk);
    for (i = 0; i < nodenum; i++)
    {
        fscanf(inputPtr, "%d%lf%lf%lf", &junk, &pos.mat[i*3][0], &pos.mat[i*3+1][0], &pos.mat[i*3+1][0]);
    }
    fscanf(inputPtr,"%s %s %s %s %s", &junk,&junk,&junk,&junk,&junk);
    for (i = 0; i < nodenum; i++)
    {
        fscanf(inputPtr, "%d%lf%lf%lf", &junk, &vel.mat[i*3][0], &vel.mat[i*3+1][0], &vel.mat[i*3+2][0]);
    }
    fscanf(inputPtr,"%s %s %s %s %s", &junk,&junk,&junk,&junk,&junk);
    for (i = 0; i < nodenum; i++)
    {
        fscanf(inputPtr, "%d%lf%lf%lf", &junk, &acc.mat[i*3][0], &acc.mat[i*3+1][0], &acc.mat[i*3+2][0]);
    }

    Matrix D;
    D = make_mat(6,6);
    Matrix KLoc;
    KLoc = make_mat(24,24);
    double xi[3], W[3];
    xi[0] = 0.7745966692;
    xi[1] = -0.7745966692;
    xi[2] = 0;
    W[0] = 0.5555555556;
    W[1] = 0.5555555556;
    W[2] = 0.8888888889;
    Matrix dphi;
    dphi = make_mat(3,8);
    Matrix position;
    position = make_mat(8,3);
    Matrix J;
    J = make_mat(3,3);
    Matrix Ji;
    Ji = make_mat(3,3);
    Matrix B, Bnew;
    B = make_mat(3,8);
    Bnew = make_mat(6,24);
    Matrix MLoc, MT;
    MT = make_mat(24,24);
    MLoc = make_mat(24, 24);
    Matrix phi;
    phi = make_mat(1,8);
    int m,n;
    FILE *Mlocout = fopen("M_Local.txt", "w");

    int e, ngp = 3;
    double r, s, t;

    FILE *Klocout = fopen("K_Local.txt", "w");

    for (e = 0; e < elementnum; e++)
    {
        for (i = 0; i < ngp; i++)
        {
            for (j = 0; j < ngp; j++)
            {
                for (k = 0; k < ngp; k++)
                {
                    r = xi[i];
                    s = xi[j];
                    t = xi[k];

                    phi.mat[0][0] = 0.125 * (1-r) * (1-s) * (1-t);
                    phi.mat[0][1] = 0.125 * (1+r) * (1-s) * (1-t);
                    phi.mat[0][2] = 0.125 * (1+r) * (1+s) * (1-t);
                    phi.mat[0][3] = 0.125 * (1-r) * (1+s) * (1-t);
                    phi.mat[0][4] = 0.125 * (1-r) * (1-s) * (1+t);
                    phi.mat[0][5] = 0.125 * (1+r) * (1-s) * (1+t);
                    phi.mat[0][6] = 0.125 * (1+r) * (1+s) * (1+t);
                    phi.mat[0][7] = 0.125 * (1-r) * (1+s) * (1+t);

                    dphi.mat[0][0] = -0.125 * (1-s) * (1-t);
                    dphi.mat[0][1] = 0.125 * (1-s) * (1-t);
                    dphi.mat[0][2] = 0.125 * (1+s) * (1-t);
                    dphi.mat[0][3] = -0.125 * (1+s) * (1-t);
                    dphi.mat[0][4] = -0.125 * (1-s) * (1+t);
                    dphi.mat[0][5] = 0.125 * (1-s) * (1+t);
                    dphi.mat[0][6] = 0.125 * (1+s) * (1+t);
                    dphi.mat[0][7] = -0.125 * (1+s) * (1+t);

                    dphi.mat[1][0] = -0.125 * (1-r) * (1-t);
                    dphi.mat[1][1] = -0.125 * (1+r) * (1-t);
                    dphi.mat[1][2] = 0.125 * (1+r) * (1-t);
                    dphi.mat[1][3] = 0.125 * (1-r) * (1-t);
                    dphi.mat[1][4] = -0.125 * (1-r) * (1+t);
                    dphi.mat[1][5] = -0.125 * (1+r) * (1+t);
                    dphi.mat[1][6] = 0.125 * (1+r) * (1+t);
                    dphi.mat[1][7] = 0.125 * (1-r) * (1+t);

                    dphi.mat[2][0] = -0.125 * (1-r) * (1-s);
                    dphi.mat[2][1] = -0.125 * (1+r) * (1-s);
                    dphi.mat[2][2] = -0.125 * (1+r) * (1+s);
                    dphi.mat[2][3] = -0.125 * (1-r) * (1+s);
                    dphi.mat[2][4] = 0.125 * (1-r) * (1-s);
                    dphi.mat[2][5] = 0.125 * (1+r) * (1-s);
                    dphi.mat[2][6] = 0.125 * (1+r) * (1+s);
                    dphi.mat[2][7] = 0.125 * (1-r) * (1+s);

                    for (p = 0; p < 8; p++)
                    {
                        position.mat[p][0] = node[NodeX[e][p]-1][0];
                        position.mat[p][1] = node[NodeX[e][p]-1][1];
                        position.mat[p][2] = node[NodeX[e][p]-1][2];
                    }

                    J = Matrix_Mult(dphi, position);

                    Ji.mat[0][0] = pow(det(J),-1) * (J.mat[1][1]*J.mat[2][2] - J.mat[1][2]*J.mat[2][1]); //mali
                    Ji.mat[0][1] = pow(det(J),-1) * (J.mat[0][2]*J.mat[2][1] - J.mat[0][1]*J.mat[2][2]);
                    Ji.mat[0][2] = pow(det(J),-1) * (J.mat[0][1]*J.mat[1][2] - J.mat[0][2]*J.mat[1][1]);
                    Ji.mat[1][0] = pow(det(J),-1) * (J.mat[1][2]*J.mat[2][0] - J.mat[1][0]*J.mat[2][2]);
                    Ji.mat[1][1] = pow(det(J),-1) * (J.mat[0][0]*J.mat[2][2] - J.mat[0][2]*J.mat[2][0]);
                    Ji.mat[1][2] = pow(det(J),-1) * (J.mat[0][2]*J.mat[1][0] - J.mat[0][0]*J.mat[1][2]);
                    Ji.mat[2][0] = pow(det(J),-1) * (J.mat[1][0]*J.mat[2][1] - J.mat[1][1]*J.mat[2][0]);
                    Ji.mat[2][1] = pow(det(J),-1) * (J.mat[0][1]*J.mat[2][0] - J.mat[0][0]*J.mat[2][1]);
                    Ji.mat[2][2] = pow(det(J),-1) * (J.mat[0][0]*J.mat[1][1] - J.mat[0][1]*J.mat[1][0]);

                    B = Matrix_Mult(Ji,dphi);

                    for (p = 0; p < 8; p++)
                    {
                        Bnew.mat[0][p*3] = B.mat[0][p];
                        Bnew.mat[1][p*3+1] = B.mat[1][p];
                        Bnew.mat[2][p*3+2] = B.mat[2][p];
                        Bnew.mat[3][p*3] = B.mat[1][p];
                        Bnew.mat[3][p*3+1] = B.mat[0][p];
                        Bnew.mat[4][p*3+1] = B.mat[2][p];
                        Bnew.mat[4][p*3+2] = B.mat[1][p];
                        Bnew.mat[5][p*3] = B.mat[2][p];
                        Bnew.mat[5][p*3+2] = B.mat[0][p];
                    }

                    D.mat[0][0] = E[e]*pow(1+v[e],-1)*pow(1-2*v[e],-1) * (1-v[e]);
                    D.mat[0][1] = E[e]*pow(1+v[e],-1)*pow(1-2*v[e],-1) * (v[e]);
                    D.mat[0][2] = E[e]*pow(1+v[e],-1)*pow(1-2*v[e],-1) * (v[e]);
                    D.mat[1][0] = E[e]*pow(1+v[e],-1)*pow(1-2*v[e],-1) * (v[e]);
                    D.mat[1][1] = E[e]*pow(1+v[e],-1)*pow(1-2*v[e],-1) * (1-v[e]);
                    D.mat[1][2] = E[e]*pow(1+v[e],-1)*pow(1-2*v[e],-1) * (v[e]);
                    D.mat[2][0] = E[e]*pow(1+v[e],-1)*pow(1-2*v[e],-1) * (v[e]);
                    D.mat[2][1] = E[e]*pow(1+v[e],-1)*pow(1-2*v[e],-1) * (v[e]);
                    D.mat[2][2] = E[e]*pow(1+v[e],-1)*pow(1-2*v[e],-1) * (1-v[e]);
                    D.mat[3][3] = E[e]*pow(1+v[e],-1)*pow(1-2*v[e],-1) * (0.5-v[e]);
                    D.mat[4][4] = E[e]*pow(1+v[e],-1)*pow(1-2*v[e],-1) * (0.5-v[e]);
                    D.mat[5][5] = E[e]*pow(1+v[e],-1)*pow(1-2*v[e],-1) * (0.5-v[e]);

                    KLoc = Matrix_Add(KLoc,Matrix_Scalar(Matrix_Mult(Matrix_Mult(Matrix_Trans(Bnew),D),Bnew),W[i]*W[j]*W[k]*det(J)));

                    for (m = 0; m < 8; m++)
                    {
                        for (n = 0; n < 8; n++)
                        {
                            MT.mat[m*3][n*3] = W[i]*W[j]*W[k]*phi.mat[0][m]*phi.mat[0][n]*rho[e]*det(J);
                            MT.mat[m*3+1][n*3+1] = W[i]*W[j]*W[k]*phi.mat[0][m]*phi.mat[0][n]*rho[e]*det(J);
                            MT.mat[m*3+2][n*3+2] = W[i]*W[j]*W[k]*phi.mat[0][m]*phi.mat[0][n]*rho[e]*det(J);
                        }
                    }
                    MLoc = Matrix_Add(MLoc,MT);
                    zero_mat(MT);
                }
            }
        }
        fprintf(Klocout, "Element %d\n", e+1);
        for (p = 0; p < 24; p++)
        {
            for (q = 0; q < 24; q++)
            {
                fprintf(Klocout, "%lf\t", KLoc.mat[p][q]);
            }
            fprintf(Klocout, "\n");
        }
        fprintf(Klocout, "\n");
        fprintf(Mlocout, "Element %d\n", e+1);
        for (p = 0; p < 24; p++)
        {
            for (q = 0; q < 24; q++)
            {
                fprintf(Mlocout, "%lf\t", MLoc.mat[p][q]);
            }
            fprintf(Mlocout, "\n");
        }
        fprintf(Mlocout, "\n");
        zero_mat(phi);
        zero_mat(dphi);
        zero_mat(position);
        zero_mat(J);
        zero_mat(Ji);
        zero_mat(Bnew);
        zero_mat(D);
        zero_mat(KLoc);
        zero_mat(MLoc);
    }

    fclose(Klocout);

    Matrix KGlobal, L;
    KGlobal = make_mat(nodenum*3, nodenum*3);
    Matrix MGlobal;
    MGlobal = make_mat(nodenum*3, nodenum*3);
    L = make_mat(24, nodenum*3);
    FILE *Klocin = fopen("K_Local.txt", "r");
    FILE *Mlocin = fopen("M_Local.txt", "r");

    Matrix KGadd, MGadd;
    KGadd = make_mat(nodenum*3, nodenum*3);
    MGadd = make_mat(nodenum*3, nodenum*3);

    for(i = 0; i < elementnum; i++)
    {
        for(j = 0; j < 8; j++)
        {
            L.mat[j*3][NodeX[i][j]*3-3] = 1;
            L.mat[j*3+1][NodeX[i][j]*3-2] = 1;
            L.mat[j*3+2][NodeX[i][j]*3-1] = 1;
        }
        fscanf(Klocin, "%s %s", &junk, &junk);
        for (j = 0; j < 24; j++)
        {
            for (k = 0; k < 24; k++)
            {
                fscanf(Klocin, "%lf", &KLoc.mat[j][k]);
            }
        }
        KGadd = Matrix_Mult(Matrix_Mult(Matrix_Trans(L),KLoc),L);
        KGlobal = Matrix_Add(KGlobal, KGadd);
        fscanf(Mlocin, "%s %s", &junk, &junk);
        for (j = 0; j < 24; j++)
        {
            for (k = 0; k < 24; k++)
            {
                fscanf(Mlocin, "%lf", &MLoc.mat[j][k]);
            }
        }
        MGadd = Matrix_Mult(Matrix_Mult(Matrix_Trans(L),MLoc),L);
        MGlobal = Matrix_Add(MGlobal, MGadd);
        //print_mat(KGadd);
        zero_mat(L);
        zero_mat(MGadd);
        zero_mat(KGadd);
    }

    FILE *KM = fopen("KM.csv", "w");
    fprintf(KM, "K\n");
    for (i = 0; i < nodenum*3; i++)
    {
        for (j = 0; j < nodenum*3; j++)
        {
            fprintf(KM, "%lf,", KGlobal.mat[i][j]);
        }
        fprintf(KM, "\n");
    }
    fprintf(KM, "M\n");
    for (i = 0; i < nodenum*3; i++)
    {
        for (j = 0; j < nodenum*3; j++)
        {
            fprintf(KM, "%lf,", MGlobal.mat[i][j]);
        }
        fprintf(KM, "\n");
    }
    fprintf(KM, "pos\n");
    for (i = 0; i < nodenum*3; i++)
    {
        fprintf(KM, "%lf\n", pos.mat[i][0]);
    }

    Matrix ga;
    ga = make_mat(nodenum*3,1);

    Matrix f;
    f = make_mat(nodenum*3,1);

    double dt = 0.1;
    Matrix TK;
    TK = make_mat(nodenum*3, nodenum*3);
    TK = Matrix_Add(KGlobal,Matrix_Scalar(MGlobal, 4*pow(dt,-2)));
    Matrix TF;
    TF = make_mat(nodenum*3, 1);
    Matrix posnew, velnew, accnew, posnewa;
    posnew = make_mat(nodenum*3,1);
    posnewa = make_mat(270,1);
    velnew = make_mat(nodenum*3,1);
    accnew = make_mat(nodenum*3,1);
    FILE *out = fopen("output.csv", "w");
    fprintf(out, "t,%lf\n,u,v,a\n",0);
    for(j = 0; j < nodenum*3; j++)
    {
        fprintf(out, "%d,%lf,", j+1, pos.mat[j][0]);
        fprintf(out, "%lf,", vel.mat[j][0]);
        fprintf(out, "%lf\n", acc.mat[j][0]);
    }
    fprintf(out,"\n");

    FILE *resfile = fopen("highrise.res", "w");
    fprintf(resfile, "GiD Post Results File 1.0\n\n");
    fprintf(resfile, "Result \"Displacement\" \"GNS-analysis\"      1 Vector OnNodes\nComponentNames \"x_disp\",\"y_disp\",\"z_disp\"\nValues\n");
    for(j = 0; j < nodenum; j++)
    {
        fprintf(resfile, "%d\t%lf\t%lf\t%lf\n", j+1, pos.mat[j*3][0],pos.mat[j*3+1][0],pos.mat[j*3+2][0]);
    }
    fprintf(resfile, "End Values\n\n");
    Matrix TKA,TFA;
    TKA = make_mat(270, 270);
    TFA = make_mat(270, 1);
    FILE* Kmatout = fopen("Kmatout.csv", "w");
    FILE* Fmatout = fopen("Fmatout.csv", "w");
    FILE* top = fopen("topdef.csv", "w");
    for (i = 0; i < 100; i++)
    {
        printf("%d\n", i);
        f = Matrix_Mult(MGlobal,ga);
        TF = Matrix_Mult(Matrix_Scalar(MGlobal,4*pow(dt,-2)),pos);
        TF = Matrix_Add(TF,Matrix_Mult(Matrix_Scalar(MGlobal,4*pow(dt,-1)),vel));
        TF = Matrix_Add(TF,Matrix_Mult(MGlobal,acc));
        TF = Matrix_Add(TF,f);
        for (p = 0; p < 270; p++)
        {
            for (q = 0; q < 270; q++)
            {
                TKA.mat[p][q] = TK.mat[p][q];
                fprintf(Kmatout, "%lf,", TKA.mat[p][q]);
            }
            TFA.mat[p][0] = TF.mat[p][0];
            fprintf(Fmatout, "%lf\n", TFA.mat[p][0]);
            fprintf(Kmatout, "\n");
        }
        fprintf(Kmatout, "\n");
        fprintf(Fmatout, "\n");
        posnewa = Gauss_Jordan(TKA,TFA);
        for(p=0;p<270;p++)
        {
            posnew.mat[p][0] = posnewa.mat[p][0];
        }
        velnew = Matrix_Scalar(vel,-1);
        velnew = Matrix_Add(velnew,Matrix_Scalar(pos,-2*pow(dt,-1)));
        velnew = Matrix_Add(velnew,Matrix_Scalar(posnew,2*pow(dt,-1)));
        accnew = Matrix_Add(accnew,Matrix_Scalar(posnew,4*pow(dt,-2)));
        accnew = Matrix_Add(accnew,Matrix_Scalar(pos,-4*pow(dt,-2)));
        accnew = Matrix_Add(accnew,Matrix_Scalar(vel,-4*pow(dt,-1)));
        accnew = Matrix_Add(accnew,Matrix_Scalar(acc,-1));
        fprintf(out, "t,%lf\n,u,v,a\n",(i+1)*dt);
        fprintf(top, "%lf,%lf\n", dt*i, pos.mat[0][0]);
        for(j = 0; j < nodenum*3; j++)
        {
            fprintf(out, "%d,%lf,", j+1, posnew.mat[j][0]);
            pos.mat[j][0] = posnew.mat[j][0];
            fprintf(out, "%lf,", velnew.mat[j][0]);
            vel.mat[j][0] = velnew.mat[j][0];
            fprintf(out, "%lf\n", accnew.mat[j][0]);
            acc.mat[j][0] = accnew.mat[j][0];
        }
        fprintf(out,"\n");

        fprintf(resfile, "Result \"Displacement\" \"GNS-analysis\"      %d Vector OnNodes\nComponentNames \"x_disp\",\"y_disp\",\"z_disp\"\nValues\n", i+2);
        for(j = 0; j < nodenum; j++)
        {
            fprintf(resfile, "%d\t%lf\t%lf\t%lf\n", j+1, posnew.mat[j*3][0],posnew.mat[j*3+1][0],posnew.mat[j*3+2][0]);
        }
        fprintf(resfile, "End Values\n\n");
        zero_mat(posnew);
        zero_mat(velnew);
        zero_mat(accnew);
        zero_mat(posnewa);
    }
    printf("Enter");
    getch();
}
