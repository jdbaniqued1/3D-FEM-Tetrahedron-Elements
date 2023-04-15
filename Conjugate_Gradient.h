Matrix Conjugate_Gradient(Matrix A, Matrix b)
{
    int m, n = A.row, temp, i, j, mmax;
    mmax = 10000;
    double e, a, B;
    e = 0.0001;

    Matrix x;
    x.row = n;
    x.col = 1;

    x.mat=(double **) malloc(sizeof(double*)*x.row);
    for (i=0;i<x.row;i++)
        x.mat[i]=(double *) malloc(sizeof (double)*x.col);

    for (i=0; i<n; i++)
    {
        x.mat[i][0] = 0.3;
    }

    //Step 1
    m = 1;

    //Step 2
    Matrix r;
    r.row = n;
    r.col = 1;

    r.mat=(double **) malloc(sizeof(double*)*r.row);
    for (i=0;i<r.row;i++)
        r.mat[i]=(double *) malloc(sizeof (double)*r.col);

    r = Matrix_Mult(A,x);
    r = Matrix_Sub(b,r);

    //Step 3
    Matrix d;
    d.row = n;
    d.col = 1;

    d.mat=(double **) malloc(sizeof(double*)*d.row);
    for (i=0;i<d.row;i++)
        d.mat[i]=(double *) malloc(sizeof (double)*d.col);

    Matrix D;
    D.row = n;
    D.col = 1;

    D.mat=(double **) malloc(sizeof(double*)*D.row);
    for (i=0;i<D.row;i++)
        D.mat[i]=(double *) malloc(sizeof (double)*D.col);

    for (i=0; i<n; i++)
    {
        d.mat[i][0] = r.mat[i][0];
    }

    //Step 4
    double c0, cold, cnew;
    Matrix rT;
    rT.row = 1;
    rT.col = n;

    rT.mat=(double **) malloc(sizeof(double*)*rT.row);
    for (i=0;i<rT.row;i++)
        rT.mat[i]=(double *) malloc(sizeof (double)*rT.col);

    for(i=0; i<n; i++)
    {
        rT.mat[0][i] = r.mat[i][0];
    }

    Matrix Scalar;
    Scalar.row = 1;
    Scalar.col = 1;

    Scalar.mat=(double **) malloc(sizeof(double*)*Scalar.row);
    for (i=0;i<Scalar.row;i++)
        Scalar.mat[i]=(double *) malloc(sizeof (double)*Scalar.col);

    Scalar = Matrix_Mult(rT,r);
    cnew = Scalar.mat[0][0];

    //Step 5
    c0 = cnew;

    Matrix q;
    q.row = n;
    q.col = 1;

    q.mat=(double **) malloc(sizeof(double*)*q.row);
    for (i=0;i<q.row;i++)
        q.mat[i]=(double *) malloc(sizeof (double)*q.col);

    Matrix dT;
    dT.row = 1;
    dT.col = n;

    dT.mat=(double **) malloc(sizeof(double*)*dT.row);
    for (i=0;i<dT.row;i++)
        dT.mat[i]=(double *) malloc(sizeof (double)*dT.col);

    for(i=0; i<n; i++)
    {
        dT.mat[0][i] = d.mat[i][0];
    }

    //Step 6
    while (m < mmax && cnew > (pow(e,2)*c0))
    {
        //Step 7
        q = Matrix_Mult(A,d);

        //Step 8
        Scalar = Matrix_Mult(dT,q);
        a = cnew/Scalar.mat[0][0];

        //Step 9
        D = Matrix_Scalar(d,a);
        x = Matrix_Add(x,D);

        //Step 10
        q = Matrix_Scalar(q,a);
        r = Matrix_Sub(r,q);

        //Step 11
        cold = cnew;

        //Step 12
        for(i=0; i<n; i++)
        {
            rT.mat[0][i] = r.mat[i][0];
        }
        Scalar = Matrix_Mult(rT,r);
        cnew = Scalar.mat[0][0];

        //Step 13
        B = cnew/cold;

        //Step 14
        d = Matrix_Scalar(d,B);
        d = Matrix_Add(r,d);
        for(i=0; i<n; i++)
        {
            dT.mat[0][i] = d.mat[i][0];
        }

        //Step 15
        m = m+1;
    }

    return x;
}
