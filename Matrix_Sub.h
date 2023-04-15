Matrix Matrix_Sub(Matrix A, Matrix B)
{
    Matrix C;
    int i,j;


    if(A.row!=B.row)
    {
        printf("Row sizes are not equal\n");
        printf("Terminating program.\n");
        return C;
    }
    else if(A.col != B.col)
    {
        printf("Column sizes are not equal\n");
        printf("Terminating program.\n");
        return C;
    }

    C.row=A.row;
    C.col=A.col;
    C.mat=(double **)malloc(sizeof(double *)*C.row);
    for(i=0;i<C.row;i++)
        C.mat[i]=(double *)malloc(sizeof(double)*C.col);

    for(i=0;i<C.row;i++)
    {
        for(j=0;j<C.col; j++)
            C.mat[i][j]=A.mat[i][j]-B.mat[i][j];
    }
    return C;
}
