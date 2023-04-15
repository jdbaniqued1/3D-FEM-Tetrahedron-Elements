
Matrix Gauss_Jordan(Matrix A,Matrix B)
{

   Matrix D;
   Matrix X;
   double *ptr,store;
    int i,j,a,b,h;
    int k;

    if(A.row!=B.row)
    {
        printf("Number of rows of A and B must be equal!\n");
   }

    X.row=B.row;
    X.col=B.col;

    X.mat=(double **) malloc(sizeof(double*)*X.row);
    for (i=0;i<X.row;i++)
    X.mat[i]=(double*)malloc(sizeof(double)*X.col);

    D=Aug(A,B);

//Reduced Echelon Form
    for(i=0;i<D.row;i++)
    {
        store=D.mat[i][i];
        for(j=0;j<D.col;j++)
        {
            if(D.mat[i][i]==0)
            {
                ptr=D.mat[i];
                D.mat[i]=D.mat[i+1];
                D.mat[i+1]=ptr;

            }

            D.mat[i][j]=D.mat[i][j]/store;
        }

       for(a=0;a<D.row;a++)
        {
            if(a==i)
                continue;
            store=D.mat[a][i];
            for(b=0;b<D.col;b++)
            {
                D.mat[a][b]=D.mat[a][b]-store*D.mat[i][b];
            }
       }
    }

//MATRIX SEPARATION
for(i=0;i<X.row;i++)
    {
    k=A.col;
        for(j=0;j<X.col;j++)
        {
        X.mat[i][j]=D.mat[i][k];
        k++;
        }
    }
        return X;    //answer
}
