Matrix Aug(Matrix A,Matrix B)
{
    Matrix C;
    int i,j,k;
    C.row=A.row;
    C.col=A.col+B.col;


    C.mat=(double **) malloc(sizeof(double*)*C.row);
    for (i=0;i<C.row;i++)
        C.mat[i]=(double *) malloc(sizeof (double)*C.col);

//Matrix Augmentation
   for(i=0;i<C.row;i++)
    {
        k=0;
        for(j=0;j<(C.col);j++)
        {
           if(j>=A.col)
            {
                C.mat[i][j]=B.mat[i][k];
                k++;
            }
        else

            C.mat[i][j]=A.mat[i][j];
        }
    }

    return C;//augmented matrix
}
