Matrix Matrix_Scalar(Matrix A, double num)
{
    Matrix C;
    int i,j;
    C.row=A.row;
    C.col=A.col;
    C.mat=(double **)malloc(sizeof(double *)*C.row);
    for(i=0;i<C.row;i++)
        C.mat[i]=(double *)malloc(sizeof(double)*C.col);

    for(i=0;i<C.row;i++)
    {
        for(j=0;j<C.col; j++)
            C.mat[i][j]=A.mat[i][j] * num;
    }

    return C;
}
