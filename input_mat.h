Matrix input_mat()
{
    Matrix A;
    int i,j;
    printf("Enter number of rows:");
    scanf("%d",&A.row);
    printf("Enter number of columns:");
    scanf("%d",&A.col);


    A.mat=(double **)malloc(sizeof(double *)*A.row);
    for(i=0;i<A.row;i++)
    {
        A.mat[i]=(double*)malloc(sizeof(double)*A.col);
    }

    for (i=0;i<A.row;i++)
        for(j=0;j<A.col;j++)
        {
           printf("Matrix [%d][%d]:",i,j);
           scanf("%lf",&A.mat[i][j]);
        }
    return A;
}
