#include <stdio.h>
#include <stdlib.h>

//����ת��
void mat_transpose(double *A, double *AT, int nrow, int ncol)
{
    for(int i = 0; i < nrow; i++)
        for(int j = 0; j < ncol; j++)
            AT[j * nrow + i] = A[i * ncol + j];
}

//����˷� A(n*m) �� B(m*p) = C(n*p)
void mat_mult(double *A, double *B, double *C, int n, int m, int p)
{
    for(int i = 0; i < n; i++){
        for(int j = 0; j < p; j++){
            C[i * p + j] = 0.0;
            for(int k = 0; k < m; k++){
                C[i * p + j] += A[i * m + k] * B[k * p + j];
            }
        }
    }
}

//��˹Լ���������棬n=6
int mat_inv(double *A, double *invA, int n)
{
    int i,j,k;
    double aug[6][12];
    for(i=0;i<n;i++){
        for(j=0;j<n;j++) aug[i][j] = A[i*n+j];
        for(j=n;j<2*n;j++) aug[i][j] = (i == j-n) ? 1.0 : 0.0;
    }
    for(k=0;k<n;k++){
        double pivot = aug[k][k];
        if(fabs(pivot) < 1e-10) return 0;
        for(j=k;j<2*n;j++) aug[k][j] /= pivot;
        for(i=0;i<n;i++){
            if(i != k && fabs(aug[i][k])>1e-10){
                double fac = aug[i][k];
                for(j=k;j<2*n;j++) aug[i][j] -= fac * aug[k][j];
            }
        }
    }
    for(i=0;i<n;i++)
        for(j=0;j<n;j++) invA[i*n+j] = aug[i][j+n];
    return 1;
}

int main()
{
    const int N = 400;  //��������
    const int p = 6;    //������������0,Price,Income,Advertising,Good,Medium
    double *X  = (double*)malloc(N * p * sizeof(double));
    double *y  = (double*)malloc(N * sizeof(double));

    // ========== Carseats 400������ ==========
    int row = 0;
    //��ʽ��X[row] = [1, Price, Income, Advertising, ShelveLoc_Good, ShelveLoc_Medium], y[row]=Sales
    //1
    X[row*p+0]=1; X[row*p+1]=120; X[row*p+2]=73; X[row*p+3]=11; X[row*p+4]=0; X[row*p+5]=0; y[row]=9.50; row++;
    //2
    X[row*p+0]=1; X[row*p+1]=83; X[row*p+2]=48; X[row*p+3]=16; X[row*p+4]=1; X[row*p+5]=0; y[row]=11.22; row++;
    //3
    X[row*p+0]=1; X[row*p+1]=80; X[row*p+2]=35; X[row*p+3]=10; X[row*p+4]=0; X[row*p+5]=1; y[row]=10.06; row++;
    //4
    X[row*p+0]=1; X[row*p+1]=97; X[row*p+2]=100;X[row*p+3]=4; X[row*p+4]=0; X[row*p+5]=1; y[row]=7.40; row++;
    //5
    X[row*p+0]=1; X[row*p+1]=128;X[row*p+2]=64; X[row*p+3]=3; X[row*p+4]=0; X[row*p+5]=0; y[row]=4.15; row++;
    //6
    X[row*p+0]=1; X[row*p+1]=72; X[row*p+2]=113;X[row*p+3]=13;X[row*p+4]=0; X[row*p+5]=0; y[row]=10.81; row++;
    //7
    X[row*p+0]=1; X[row*p+1]=108;X[row*p+2]=105;X[row*p+3]=0; X[row*p+4]=0; X[row*p+5]=1; y[row]=6.63; row++;
    //8
    X[row*p+0]=1; X[row*p+1]=120;X[row*p+2]=81; X[row*p+3]=15;X[row*p+4]=1; X[row*p+5]=0; y[row]=11.85; row++;
    //9
    X[row*p+0]=1; X[row*p+1]=124;X[row*p+2]=110;X[row*p+3]=0; X[row*p+4]=0; X[row*p+5]=1; y[row]=6.54; row++;
    //10
    X[row*p+0]=1; X[row*p+1]=124;X[row*p+2]=113;X[row*p+3]=0; X[row*p+4]=0; X[row*p+5]=1; y[row]=4.69; row++;
 

    double *XT    = (double*)malloc(p * N * sizeof(double));
    double *XTX   = (double*)malloc(p * p * sizeof(double));
    double *XTy   = (double*)malloc(p * sizeof(double));
    double *invXTX= (double*)malloc(p * p * sizeof(double));
    double *beta  = (double*)malloc(p * sizeof(double));

    mat_transpose(X, XT, N, p);
    mat_mult(XT, X, XTX, p, N, p);
    mat_mult(XT, y, XTy, p, N, 1);

    if(!mat_inv(XTX, invXTX, p)){
        printf("���󲻿��棡\n");
        return -1;
    }
    mat_mult(invXTX, XTy, beta, p, p, 1);

    printf("==== Carseats��Ԫ���Իع��� ====\n");
    printf("�ؾ� ��0        = %.4f\n", beta[0]);
    printf("Priceϵ�� ��1   = %.4f\n", beta[1]);
    printf("Incomeϵ�� ��2  = %.4f\n", beta[2]);
    printf("Advertising ��3 = %.4f\n", beta[3]);
    printf("ShelveLoc_Good ��4= %.4f\n", beta[4]);
    printf("ShelveLoc_Medium ��5= %.4f\n", beta[5]);
    printf("\n�ع鷽�̣�\nSales = %.4f + %.4f*Price + %.4f*Income + %.4f*Advertising + %.4f*ShelveLoc_Good + %.4f*ShelveLoc_Medium\n",
           beta[0], beta[1], beta[2], beta[3], beta[4], beta[5]);

    free(X);free(y);free(XT);free(XTX);free(XTy);free(invXTX);free(beta);
    return 0;
}

