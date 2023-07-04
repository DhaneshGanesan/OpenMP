
#include<time.h>
#include<stdlib.h>
#include<stdio.h>
#include<omp.h>
int main()
{
     clock_t begin = clock();    
    int iter;
    int i,j;
    int N=11;
    double (*arr)[N+1][N+1] = malloc (sizeof(double[2][N+1][N+1])) ;  /*array initialization*/
    int limit=10;
    for( i=0;i<=N;i++)
    {
        for( j=0;j<=N;j++)
        {
            arr[0][i][j]=10;
            arr[1][i][j]=10;
        }
    }
    int wall = (N+1)*0.3;                                     /*wall measurement before radiator*/
    int radiator=((N+1)*0.4)+1;                                /*wall measurement at radiator*/
    for(i=0;i<=N;i++)
    {
        if(i>=wall)
        {
            if(radiator>0)
                {
                    arr[0][i][N]=100;
                    arr[1][i][N]=100;
                    radiator=radiator-1;
                }
            else
                {
                    arr[0][i][N]=10;
                    arr[1][i][N]=10;
                }
        }
        else
        {
            arr[0][i][N]=10;

            arr[1][i][N]=10;
        }
    }
    int x=0;int t=1;int swap;
    for ( iter = 0 ; iter < limit ; iter++)
    {
        for ( i = 1 ; i < N ; i++)
        {
            for ( j = 1 ; j < N ; j++)
            {
                arr[t][i][j] = ( arr[x][i-1][j] + arr[x][i+1][j] + arr[x][i][j-1] + arr[x][i][j+1] ) * 0.25;    /*Laplace equation*/
  }
        }
        swap=x;
        x=t;
        t=swap;
    } clock_t end = clock();
  
    double time_spent=(double)(end-begin)/CLOCKS_PER_SEC;
    printf("%f",time_spent);
    
    printf("SERIAL :");
for ( i = 0 ; i <= N ; i++)
{
    if(i%(N/8)==0)
    {
        printf("%dth row\t t[%d][0] \n",i,i);
        for(j=0; j<N+1; j++) {
            if(j%(N/8)==0)
            {
                printf("%f   ",arr[t][i][j]);                       /*Printing Array[N/8][N/8]*/

            }
        }
    }
}
}


