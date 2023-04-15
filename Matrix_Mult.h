Matrix Matrix_Mult(Matrix A, Matrix B)
{
    Matrix C;
    int i,j, k;
    if (A.col != B.row)
	{
	printf("The number of columns of matrix A and rows of matrix B are not equal, operation cannot proceed");
     	printf("Terminating program.\n");
        return C;
	}

    C.row=A.row;
    C.col=B.col;
    C.mat=(double **)malloc(sizeof(double *)*C.row);
    for(i=0;i<C.row;i++)
        C.mat[i]=(double *)calloc(C.col,sizeof(double)); //calloc is used for initializing the elements to zero.

    for(i=0;i<C.row;i++)
    {
        for(j=0;j<C.col; j++)
        {
            for (k=0; k <A.col; k++)
                C.mat[i][j] += A.mat[i][k] * B.mat[k][j];
        }
    }

    return C;
}
