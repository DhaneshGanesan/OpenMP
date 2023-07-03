#include<stdio.h>
#include<mpi.h>
#include<stdlib.h>
#include <math.h>

#define MAX_SEARCH_ITERATIONS 1000

double get_final_temperature(double, double, double);

// Generate N between start and end
void linspace(double start, double end, int n, double* res) {
    double stepsize = (end - start)/(n-1);
    res[0] = start;
    for (int i = 1; i < n-1; i++) {
        res[i] = res[i-1] + stepsize;
    }
    res[n-1] = end;
}

/* bisection */
double bisection_search(double a, double b, double tol, double target, double fx, double fy) {
    int iter=0;
    while(iter<MAX_SEARCH_ITERATIONS) {
        // double mid = ...; Calculate next point to test
        double mid = (a + b)/2.0;

        double temp = get_final_temperature(mid, fx, fy);

        if (a>b) {
            printf("Invalid inputs\n");
            return -1;
        }

        //if(temperature is within tol of target) {
        //    printf("Reached target %f with temp %f and tol %f in %d iterations\n", target, temp, tol, iter);
        //    return mid;
        //} else if(temp is greater than target) {
            //do something (change either a or b)
        //} else {
            //do something else (change either a or b)
        //}

        if (fabs(temp - target) < tol) {
            return mid;
        } else if (temp > target) {
            b = mid;
        } else {
            a = mid;
        }

        if(abs(b-a) < tol) {
            printf("Interval too small %d iterations\n", iter);
            //what should be done here?
            // In this case both points have been squeezed sufficiently closed together, return mid
            return mid;
        }

        printf("Reached a temperature of %f with radiator set at %f\n", temp, mid);
        iter++;
    }

    printf("Search failed\n");
    return -1;
}

double nsection_search(int n, double a, double b, double tol, double target, double fx, double fy) {
    int iter=0;
    //Allocate arrays 
    double *temps = (double*)malloc(sizeof(double)*n);
    double *mids = (double*)malloc(sizeof(double)*n);

    if(a>b) {
            printf("Invalid inputs\n");
            return -1;
    }

    while(iter<MAX_SEARCH_ITERATIONS) {
        printf("%d\n", iter);
        // Populate array mids with n equally spaced points 
        linspace(a, b, n, mids);

        // Initialize the temps array 
        for(int i=0; i<n; i++) {
            temps[i] = get_final_temperature(mids[i], fx, fy);
        }


        for(int i=0; i<n; i++) {
            //if(temps[i] is close enough to target) {
            //    printf("Reached target %f with temp %f and tol %f in %d iterations\n", target, mids[i], tol, iter);
            //    return mids[i];
            //} else if(temps[i] < target) {
            //    printf("Temperature of %f is less than target. Setting a to %f\n", temps[i], mids[i]);
            //    do something (change either a or b)
            //} else {
            //    printf("Temperature of %f is greater than target. Setting b to %f\n", temps[i], mids[i]);
            //    do something else; (change either a or b)
            //}

            if (fabs(temps[i] - target) < tol) {
                printf("Reached target %f with temp %f and tol %f in %d iterations\n", target, mids[i], tol, iter);
                return mids[i];
            } else if (temps[i] < target) {
                printf("Temperature of %f is less than target. Setting a to %f\n", temps[i], mids[i]);
                a = mids[i];
            } else if (temps[i] > target) {
                printf("Temperature of %f is greater than target. Setting b to %f\n", temps[i], mids[i]);
                b = mids[i];
                break;
            }
        }

        
        if (fabs(b -a ) < tol) {
            printf("Interval too small %d iterations\n", iter);
          
            return target;
        }

        iter++;
    }
    printf("Search failed\n");
    return -1;
}

double nsection_search_mpi(int rank, int n, double a, double b, double tol, double target, double fx, double fy) {
    int iter=0;

    if (a>b) {
        printf("Invalid inputs\n");
        return -1;
    }

    // Allocate arrays of size n to store points to test 
    double *temps = (double*)malloc(sizeof(double)*n);
    double *mids = (double*)malloc(sizeof(double)*n);
    
    // Perhaps some local variables initialised here
    double myRadTemp, myTemp, result, stepsize;
    MPI_Request request;
    int found = 0;
    int i;
    while(iter<MAX_SEARCH_ITERATIONS && !found) {
        // The interval (a, b) should be divided into (n+1) parts here, finding n points that do this.
        // Our function of interest will be evaluated at those n points, one on each rank.
        // Calculate the point to be evaluated

        stepsize = (b - a)/(n+1);
        myRadTemp = a + stepsize*rank;

        myTemp = get_final_temperature(myRadTemp, fx, fy);

        // Logic to decide whether one of the evaluated points is close enough to the target.
        if (fabs(myTemp - target) < tol) {
            printf("Process: %d: Reached target %f with temp %f and tol %f in %d iterations\n", rank, target, myRadTemp, tol, iter);
            result = myRadTemp;
            found = 1;
        // If not, change a and b.
        } else if (myTemp < target) {
            printf("Process %d: Temperature of %f is less than target. Setting a to %f\n", rank, myTemp, myRadTemp);
            a = myRadTemp;
        } else {
            printf("Process %d: Temperature of %f is greater than target. Setting b to %f\n", rank, myTemp, myRadTemp);
            b = myRadTemp;
        }

        // Perhaps the interval (a, b) has become smaller than tol.
        if (fabs(b - a ) < tol) {
            printf("Interval too small %d iterations\n", iter);
            // What should be done here?
            result = myRadTemp;
            found = 1;
        }

        //Tell everyone whether we've found a solution, so the loop might exit on all ranks
        if (found == 1) {
            for (int i = 0; i < n; i++) {
                if (i != rank) {
                    MPI_Isend(&found, 1, MPI_INT, i, 10*i, MPI_COMM_WORLD, &request);
                    MPI_Isend(&result, 1, MPI_DOUBLE, i, 20*i, MPI_COMM_WORLD, &request);
                }
            }
        } else if (!found){
            MPI_Irecv(&found, 1, MPI_INT, MPI_ANY_SOURCE, 10*rank, MPI_COMM_WORLD, &request);
            MPI_Irecv(&result, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 20*rank, MPI_COMM_WORLD, &request);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        iter++;
    }

    if (!found) {
        printf("Search failed\n");
        result = -1;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    return result;
}


int main(int argc, char *argv[])
{
    int rank, size;
    //Do some MPI initialisation here
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Start timing
    double a = atoi(argv[1]);
    double b = atoi(argv[2]);
    double target = atoi(argv[3]);

    double startTime, endTime;

    startTime = MPI_Wtime();

    // double n = size;
    int n = size;
    double px=5.0, py=5.0, lx=10.0, ly=10.0;
    // double tol = 0.000000000000001;
    double tol = 1e-5;
    double fx=px/lx, fy=py/ly;

    //Only on one rank, print a "Starting inverse solver..." message with problem parameters
    if (rank == 0) {
        printf("Starting inverse solver with a=%f, b=%f, tol=%f\n", a, b, tol);
    }

    // Start by implementing a bisection search function - test that it works
    // Next, move to a serial nsection search function - test that it produces the same result

    double result = nsection_search_mpi(rank, n, a, b, tol, target, fx, fy);

    endTime = MPI_Wtime();

    // On a single rank, print the final result and the time taken
    if (rank == 0) {
        printf("Time taken: %f (s), radiator temperature: %f\n", endTime - startTime, result);
    }

    MPI_Finalize();

}
