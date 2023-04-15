#include<stdio.h>
#include<math.h>
//expansion by row 1, made by Maranan Christian Paul, woot!
//Notes: xrow is the row no. of the element(xrow>=1) concerned, same goes for xcol.
double cofact(Matrix A,int xrow,int xcol);

double det(Matrix A)
{
    double sum=0;
    int i;
    if(A.row!=A.col)
    {
        printf("Matrix is not a square matrix.\n");
        printf("Terminating program.\n");
        return sum;
    }

    else if(A.row==1)
    {
        return A.mat[0][0];
    }


    for(i=0;i<A.col;i++)
    {
        sum = sum + ((A.mat[0][i]) * (cofact(A,1,i+1)));
    }
    return sum;
}

double cofact(Matrix A,int xrow,int xcol)
{
    Matrix B;
    int i,j,k=0,l=0;
    double cof;
    if(A.row!=A.col)
    {
        printf("Matrix is not a square matrix.\n");
        printf("Terminating program.\n");
        return cof;
    }

    B.col=A.col - 1;
    B.row=A.row - 1;
    B.mat=(double **)malloc(sizeof(double *)*B.row);
    for(i=0;i<B.row;i++)
        B.mat[i]=(double *)malloc(sizeof(double)*B.col);

    for(i=0;i<A.row;i++)
    {
        if(i == xrow - 1)
            k--;
        else
        {
            l=0;
            for(j=0;j<A.col;j++)
            {
                if(j == xcol - 1)
                    l--;
                else
                {
                    B.mat[k][l]=A.mat[i][j];
                }

                l++;
            }
        }
        k++;
    }
        return pow(-1,xrow+xcol)*det(B);

}

