#include<stdio.h>

Matrix Matrix_Trans(Matrix A)
{
    Matrix C;
    int i,j;
    C.row=A.col;
    C.col=A.row;
    C.mat=(double **)malloc(sizeof(double *)*C.row);
    for(i=0;i<C.row;i++)
        C.mat[i]=(double *)malloc(sizeof(double)*C.col); //declaration of size of rows and columns

    for(i=0;i<C.row;i++)
    {
        for(j=0;j<C.col; j++)
            C.mat[i][j]=A.mat[j][i];
    }

    return C;
}
