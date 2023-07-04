
#include<time.h>
#include<stdlib.h>
#include<stdio.h>
#include<omp.h>
int main()
{
    double t1=omp_get_wtime();
    int iter;
    int i,j;
    int N=11;
    double (*arr)[N+1][N+1] = malloc (sizeof(double[2][N+1][N+1])) ;      /*array initialization*/
    int limit=10;
    omp_set_num_threads(1); /*thread count initialization*/
    #pragma omp for 
     for( i=0;i<=N;i++)
    {
        for( j=0;j<=N;j++)
        {
            arr[0][i][j]=10;
            arr[1][i][j]=10;
        }
    }
    int wall = (N+1)*0.3;
    int radiator=((N+1)*0.4)+1;
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
        #pragma omp parallel for private(i,j)
        for ( i = 1 ; i < N ; i++)
        {
            for ( j = 1 ; j < N ; j++)
            {
                arr[t][i][j] = ( arr[x][i-1][j] + arr[x][i+1][j] + arr[x][i][j-1] + arr[x][i][j+1] ) * 0.25;
  }
        }
        swap=x;
        x=t;
        t=swap;
    }
        double t2=omp_get_wtime();
  double time_spend=(double)t2-t1;              /*Time calculation*/
    printf("%f",time_spent);
    
    printf("Parallel :");
for ( i = 0 ; i <= N ; i++)
{
    if(i%(N/8)==0)
    {
        printf("%dth row\t t[%d][0] \n",i,i);
        for(j=0; j<N+1; j++) {
            if(j%(N/8)==0)
            {    
                printf("%f   ",arr[t][i][j]);

            }
        }
    }
}
}


