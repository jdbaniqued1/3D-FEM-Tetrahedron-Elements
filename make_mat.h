Matrix make_mat(int row, int col)
{
    Matrix A;
    int i,j;
    A.row = row;
    A.col = col;

    A.mat=(double **)malloc(sizeof(double *)*A.row);
    for(i=0;i<A.row;i++)
    {
        A.mat[i]=(double*)malloc(sizeof(double)*A.col);
    }

    for (i=0;i<A.row;i++)
        for(j=0;j<A.col;j++)
        {
           A.mat[i][j] = 0;
        }
    return A;
}
