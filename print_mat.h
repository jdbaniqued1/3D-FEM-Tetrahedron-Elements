#include<stdio.h>

void print_mat(Matrix A)
{
    int i,j;
    for(i=0;i<A.row;i++)
    {
        for(j=0;j<A.col;j++)
        printf("%lf\t",A.mat[i][j]);
        printf("\n");
    }
}
