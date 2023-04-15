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

    fscanf(inputPtr, "%s%d%s%s%s%s%s%s%s%s%s%s%s", &junk, &elementnum, &junk, &junk, &junk, &junk, &junk, &junk, &junk, &junk, &junk, &junk, &junk);
    double E[elementnum];
    double v[elementnum];
    int NodeX[elementnum][8];

    for (i = 0; i < elementnum; i++)
    {
        fscanf(inputPtr, "%d%lf%lf", &junk, &E[i], &v[i]);
        for(j = 0; j < 8; j++)
        {
            fscanf(inputPtr, "%d", &NodeX[i][j]);
        }
    }

    int bcnum, g, h, a, b;
    fscanf(inputPtr, "%s%d%s%s%s", &junk, &bcnum, &junk, &junk, &junk);
    h = bcnum;
    g = nodenum * 3 - h;
    int c[h], d[g];
    Matrix uh;
    uh = make_mat(h,1);
    for (i = 0; i < h; i++)
    {
        fscanf(inputPtr, "%d %d %lf", &a, &b, &uh.mat[i][0]);
        if (b == 1)
            c[i] = a * 3 - 3;
        if (b == 2)
            c[i] = a * 3 - 2;
        if (b == 3)
            c[i] = a * 3 - 1;
    }
    j = 0;
    k = 0;
    for (i = 0; i < nodenum * 3; i++)
    {
        if (c[j] == i)
            j += 1;
        else
        {
            d[k] = i;
            k += 1;
        }
    }

    Matrix fg;
    fg = make_mat(g,1);

    int loadnum;
    fscanf(inputPtr, "%s %d", &junk, &loadnum);
    fscanf(inputPtr, "%s %s %s", &junk, &junk, &junk);
    Matrix load;
    load = make_mat(nodenum*3,1);
    double temload;
    for (i = 0; i < loadnum; i++)
    {
        fscanf(inputPtr, "%d %d %lf", &a, &b, &temload);
        if (b == 1)
            load.mat[a*3-3][0] = temload;
        if (b == 2)
            load.mat[a*3-2][0] = temload;
        if (b == 3)
            load.mat[a*3-1][0] = temload;
    }

    for (i = 0; i < g; i++)
    {
        a = d[i];
        fg.mat[i][0] = load.mat[a][0];
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
        zero_mat(dphi);
        zero_mat(position);
        zero_mat(J);
        zero_mat(Ji);
        zero_mat(Bnew);
        zero_mat(D);
        zero_mat(KLoc);
    }

    fclose(Klocout);

    Matrix KGlobal, L;
    KGlobal = make_mat(nodenum*3, nodenum*3);
    L = make_mat(24, nodenum*3);
    FILE *Klocin = fopen("K_Local.txt", "r");
    double place;

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
        KGlobal = Matrix_Add(KGlobal, Matrix_Mult(Matrix_Mult(Matrix_Trans(L),KLoc),L));
    }

    FILE *KGOut = fopen("KGlobal.txt", "w");

    Matrix Kgg, Kgh;
    Kgg = make_mat(g,g);
    Kgh = make_mat(g,h);

    for (i = 0; i < g; i++)
    {
        for (j = 0; j < g; j++)
        {
            a = d[i];
            b = d[j];
            Kgg.mat[i][j] = KGlobal.mat[a][b];
        }
    }

    for (p = 0; p < g; p++)
    {
        for (q = 0; q < g; q++)
        {
            fprintf(KGOut, "%lf\t", Kgg.mat[p][q]);
        }
        fprintf(KGOut, "\n");
    }

    for (i = 0; i < g; i++)
    {
        for (j = 0; j < h; j++)
        {
            a = d[i];
            b = c[j];
            Kgh.mat[i][j] = KGlobal.mat[a][b];
        }
    }

    Matrix ug = make_mat(g,1);
    Matrix M = make_mat(g,1);
    M = Matrix_Sub(fg,Matrix_Mult(Kgh,uh));
    ug = Conjugate_Gradient(Kgg,Matrix_Sub(fg,Matrix_Mult(Kgh,uh)));
    FILE *output = fopen("output.csv", "w");
    double u[nodenum*3];
    for (i = 0; i < g; i++)
    {
        a = d[i];
        u[a] = -ug.mat[i][0];
    }
    for (i = 0; i < h; i++)
    {
        a = c[i];
        u[a] = 0;
    }
    for (i=0; i < (nodenum); i++)
    {
        fprintf(output, "%d, %lf, %lf, %lf\n", i+1, -u[i*3], -u[i*3+1], -u[i*3+2]);
    }
}
