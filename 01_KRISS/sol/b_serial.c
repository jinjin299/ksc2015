#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>


double genvv(double x){
	return (x*x+pow(x,4)+pow(x,6)+exp(-x*x)+cos(x)+sin(x)+tan(x));
}

int main(int argc, char **argv){

	int i,n1,n2,j,jsta,jend;
	int iter,niter;
	double xi,xf,dx;
	double tmr;
	double *ar, *br;
	/* Do not change */
	n1 = 0;
	n2 = 100000000;
	niter = 3;
	/* Do not change */
    int wsize, wrank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
	ar = (double*) malloc(sizeof(double)*n2);
	br = (double*) malloc(sizeof(double)*n2);
	jsta = n1; 
	jend = n2;
	jsta = n1+1;
	jend = n2-1;
    int prank = wrank - 1;
    int nrank = wrank + 1;
    int sdex, edex;
    int nr1, nr2;

    sdex = jsta + wrank * (jend - jsta)/wsize;
    edex = jsta + (wrank + 1) * (jend - jsta)/wsize;
    if (wrank == wsize -1)
        edex = jend;

    nr1 = sdex;
    nr2 = edex;
    if (wrank == 0)
        nr1 = n1;
    if (wrank == wsize -1)
        nr2 = n2;
	xi = 0.L;
	xf = 1.L;
	dx = (xf-xi)/(double)(n2-n1-1);
    //Init
	for(i=nr1;i<nr2;i++){
		br[i] = xi+(double)(i-n1)*dx;
		ar[i] = 0.0;
	}
    MPI_Request rq1, rq2;
    // Operation 
	for(iter=0;iter<niter;iter++){
        // Communication b/w gap
        if (wrank == 0)
        {
            MPI_Isend(br+edex-1, 1, MPI_DOUBLE, nrank, 0, MPI_COMM_WORLD,rq1);
            MPI_Recv(br+edex, 1, MPI_DOUBLE, nrank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //MPI_Send(br+edex-1, 1, MPI_DOUBLE, nrank, 0, MPI_COMM_WORLD);
            //MPI_Recv(br+edex, 1, MPI_DOUBLE, nrank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else if (wrank == wsize -1)
        {
            MPI_Recv(br+sdex-1, 1, MPI_DOUBLE, prank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Isend(br+sdex, 1, MPI_DOUBLE, prank, 0, MPI_COMM_WORLD, rq2);
            //MPI_Recv(br+sdex-1, 1, MPI_DOUBLE, prank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //MPI_Send(br+sdex, 1, MPI_DOUBLE, prank, 0, MPI_COMM_WORLD);
        }
        else
        {
            MPI_Recv(br+sdex-1, 1, MPI_DOUBLE, prank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Isend(br+sdex, 1, MPI_DOUBLE, prank, 0, MPI_COMM_WORLD, rq2);

            MPI_Isend(br+edex-1, 1, MPI_DOUBLE, nrank, 0, MPI_COMM_WORLD,rq1);
            MPI_Recv(br+edex, 1, MPI_DOUBLE, nrank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //MPI_Recv(br+sdex-1, 1, MPI_DOUBLE, prank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //MPI_Send(br+sdex, 1, MPI_DOUBLE, prank, 0, MPI_COMM_WORLD);
            //MPI_Send(br+edex-1, 1, MPI_DOUBLE, nrank, 0, MPI_COMM_WORLD);
            //MPI_Recv(br+edex, 1, MPI_DOUBLE, nrank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if(wrank==0) printf("Send done");


		for(j=sdex;j<edex;j++){
			/* Do not change */
			ar[j] = (br[j-1]+br[j+1])/4.L+br[j]/2.L+1.L/genvv(br[j]);
			/* Do not change */
		}
		for(i=nr1;i<nr2;i++){
			/* Do not change */
			br[i] = ar[i];
			/* Do not change */
		}
	}
    double ltmr;
	ltmr = 0.L;
	for(j=sdex;j<edex;j++){
		ltmr += ar[j];
	}
    MPI_Reduce(&ltmr, &tmr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if (wrank==0) printf("tmr = %16.7f\n",tmr);
	
	free(ar);
	free(br);
    MPI_Finalize();
	return 0;

}
