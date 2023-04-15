void zero_mat(Matrix A)
{
    int i,j;

    for (i=0;i<A.row;i++)
        for(j=0;j<A.col;j++)
        {
           A.mat[i][j] = 0;
        }
}
