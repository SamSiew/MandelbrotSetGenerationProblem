//////////////////////////////////////////////////////////////////////////////////////
// mandelbrot.c program: Mandelbort Set Fractal (Color Serial Code Implementation).
// --------------------------------
//  1. Draws Mandelbrot set for Fc(z) = z*z +c
//  using Mandelbrot algorithm ( boolean escape time )
//	This code is modified from the original version as available at:
//	http://rosettacode.org/wiki/Mandelbrot_set#PPM_non_interactive
// -------------------------------
// 2. Technique of creating ppm file is  based on the code of Claudio Rocchini
// http://en.wikipedia.org/wiki/Image:Color_complex_plot.jpg
// create 24 bit color graphic file ,  portable pixmap file = PPM
// see http://en.wikipedia.org/wiki/Portable_pixmap
// to see the file use external application ( graphic viewer)
//////////////////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

double* doXwork(int iX,double Cy,double CxMin, double PixelWidth, int IterationMax, double ER2);

// Main program
int main(int argc, char *argv[])
{
    /* screen ( integer) coordinate */
    int size,rank,iX,iY,iXworkload,iYworkload,counter=0,totalCounter = 0;
    double totaltime = 0;
    const int iXmax = 8000; // default
    const int iYmax = 8000; // default

    /* world ( double) coordinate = parameter plane*/
    double Cx, Cy;
    const double CxMin = -2.5;
    const double CxMax = 1.5;
    const double CyMin = -2.0;
    const double CyMax = 2.0;

    /* */
    double PixelWidth = (CxMax - CxMin)/iXmax;
    double PixelHeight = (CyMax - CyMin)/iYmax;

    /* color component ( R or G or B) is coded from 0 to 255 */
    /* it is 24 bit color RGB file */
    const int MaxColorComponentValue = 255;
    FILE * fp;
    char *filename = "Mandelbrot.ppm";
    char *comment = "# ";	/* comment should start with # */

    // RGB color array
    static unsigned char color[3];

    /* Z = Zx + Zy*i;	Z0 = 0 */
    double Zx, Zy;
    double Zx2, Zy2; /* Zx2 = Zx*Zx;  Zy2 = Zy*Zy  */
    /*  */
    int Iteration;
    const int IterationMax = 2000; // default

    /* bail-out value , radius of circle ;  */
    const double EscapeRadius = 400;
    double ER2 = EscapeRadius * EscapeRadius;

    /* Clock information */
    clock_t start, end;
    double cpu_time_used;

    MPI_Status status;

    // Creation of parallel processes
    MPI_Init(&argc, &argv);

    // find out process ID,
    // and how many processes were started
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    start = MPI_Wtime();
    int rankworkload[size];
    if(rank == 0)
    {
        /*create new file,give it a name and open it in binary mode  */
        fp = fopen(filename, "wb"); /* b -  binary mode */

        /*write ASCII header to the file (PPM file format)*/
        fprintf(fp,"P6\n %s\n %d\n %d\n %d\n", comment, iXmax, iYmax, MaxColorComponentValue);

        printf("File: %s successfully opened for writing.\n", filename);
        printf("Computing Mandelbrot Set. Please wait...\n");

        int i, offset;
        //master work
        if(size > 1){
            offset = iYmax % (size / 2);
            iYworkload = iYmax / (size / 2);
            iXworkload = iXmax / 2;
            //allocate workload and allows fifo workload
            for (i = (size-1);  i >= 0; i--) {
                rankworkload[i] = iYworkload;
                if(offset > 0){
                    rankworkload[i] += 1;
                    offset--;
                }
            }

        } else{
            rankworkload[0] = iYmax;
            iXworkload = iXmax;
        }


    }
    MPI_Bcast(&rankworkload[0], size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&iXworkload, 1, MPI_INT, 0, MPI_COMM_WORLD);
    /* compute and write image data bytes to the file */

    int iYstartfrom = ((rank/2) * rankworkload[rank]);
    int iYnumofjob = rankworkload[rank] * ((rank/2) + 1);

    int iXstartfrom = 0;
    int iXnumofjob = iXworkload;

    if(rank % 2 != 0){
        iXstartfrom = iXworkload;
        iXnumofjob = iXmax;
    }

    int isNotflag = 1;

    if(rank % 2 == 0)
    {
        int i,j;
        double *retcolor , recvcolor[iXnumofjob - iXstartfrom][3];
        for (iY = iYstartfrom; iY < iYnumofjob; iY++) {
            Cy = CyMin + (iY * PixelHeight);

            if (fabs(Cy) < (PixelHeight / 2)) {
                Cy = 0.0; /* Main antenna */
            }

            if (rank < size - 1) {
                MPI_Send(&isNotflag, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
                MPI_Send(&Cy, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
            }

            for (iX = iXstartfrom; iX < iXnumofjob; iX++) {
                retcolor = doXwork(iX, Cy, CxMin, PixelWidth, IterationMax, ER2);
                /* write color to the file */
                if (rank != 0) {
                    MPI_Send(&retcolor[0], 3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
                    counter++;
                } else {
                    for (j = 0; j < 3; j++) {
                        color[j] = retcolor[j];
                    }
                    fwrite(color, 1, 3, fp);
                }
            }
            if(size > 1){
                MPI_Recv(&recvcolor,(iXnumofjob-iXstartfrom)*3,MPI_DOUBLE,rank + 1,0,MPI_COMM_WORLD,&status);
                for (i = 0; i < (iXnumofjob-iXstartfrom) ; i++) {
                    if(rank != 0){
                        MPI_Send(&recvcolor[i][0], 3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
                        counter++;
                    } else{
                        for (j = 0; j < 3; j++) {
                            color[j] = recvcolor[i][j];
                        }
                        fwrite(color, 1, 3, fp);
                    }
                }
            }
        }
        if(size > 1){
            isNotflag = 0;
            MPI_Send(&isNotflag, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
        }
    }
    else {
        int i;
        double *retcolor, totalcolor[iXnumofjob-iXstartfrom][3];
        while (isNotflag) {
            MPI_Recv(&isNotflag, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, &status);
            if(isNotflag == 0){
                break;
            }
            MPI_Recv(&Cy, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
            for (iX = iXstartfrom; iX < iXnumofjob; iX++) {
                retcolor = doXwork(iX, Cy, CxMin, PixelWidth, IterationMax, ER2);
                /* write color to the file */
                for (i = 0; i < 3; ++i) {
                    totalcolor[iX-iXstartfrom][i] = retcolor[i];
                }
            }
            MPI_Send(&totalcolor[0],(iXnumofjob-iXstartfrom)*3,MPI_DOUBLE,rank - 1,0, MPI_COMM_WORLD);
        }
    }
    cpu_time_used = MPI_Wtime() - start;
    MPI_Reduce(&cpu_time_used, &totaltime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&counter, &totalCounter, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if(rank == 0){
        int j;
        while(totalCounter != 0){
			double recvcolor[3];
            MPI_Recv(&recvcolor, 3, MPI_DOUBLE, MPI_ANY_SOURCE, 0,MPI_COMM_WORLD,&status);
            for (j = 0; j < 3; j++) {
                color[j] = recvcolor[j];
            } 
            fwrite(color, 1, 3, fp);
            totalCounter -= 1;
        }

        fclose(fp);

        printf("Completed Computing Mandelbrot Set.\n");
        printf("File: %s successfully closed.\n", filename);
        printf("Mandelbrot computational process time: %lf\n", totaltime);
    }

    MPI_Finalize();

    return 0;
}

double* doXwork(int iX, double Cy, double CxMin, double PixelWidth, int IterationMax, double ER2) {
    int Iteration;
    double Cx, Zx, Zy, Zx2, Zy2;
    static double retcolor[3];
    Cx = CxMin + (iX * PixelWidth);
    /* initial value of orbit = critical point Z= 0 */
    Zx = 0.0;
    Zy = 0.0;
    Zx2 = Zx * Zx;
    Zy2 = Zy * Zy;

    /* */
    for (Iteration = 0; Iteration < IterationMax && ((Zx2 + Zy2) < ER2); Iteration++) {
        Zy = (2 * Zx * Zy) + Cy;
        Zx = Zx2 - Zy2 + Cx;
        Zx2 = Zx * Zx;
        Zy2 = Zy * Zy;
    };

    /* compute  pixel color (24 bit = 3 bytes) */
    if (Iteration == IterationMax) {
        // Point within the set. Mark it as black
        retcolor[0] = 0;
        retcolor[1] = 0;
        retcolor[2] = 0;
    } else {
        // Point outside the set. Mark it as white
        double c = 3 * log((double) Iteration) / log((double) (IterationMax) - 1.0);
        if (c < 1) {
            retcolor[0] = 0;
            retcolor[1] = 0;
            retcolor[2] = 255 * c;
        } else if (c < 2) {
            retcolor[0] = 0;
            retcolor[1] = 255 * (c - 1);
            retcolor[2] = 255;
        } else {
            retcolor[0] = 255 * (c - 2);
            retcolor[1] = 255;
            retcolor[2] = 255;
        }
    }
    return retcolor;
}
