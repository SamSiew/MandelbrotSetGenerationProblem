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
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

// Main program
int main()
{
    /* screen ( integer) coordinate */
    int size,rank,iX,iY,workload,offset;
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
    MPI_Init(NULL,NULL);

    // find out process ID,
    // and how many processes were started
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    start = MPI_Wtime();

    if(rank == 0)
    {
        /*create new file,give it a name and open it in binary mode  */
        fp = fopen(filename, "wb"); /* b -  binary mode */

        /*write ASCII header to the file (PPM file format)*/
        fprintf(fp,"P6\n %s\n %d\n %d\n %d\n", comment, iXmax, iYmax, MaxColorComponentValue);

        printf("File: %s successfully opened for writing.\n", filename);
        printf("Computing Mandelbrot Set. Please wait...\n");
    }
    workload = iYmax / size;
    offset = iYmax - (workload * size);
    unsigned char* retcolor;
    /* compute and write image data bytes to the file */
    int startfrom = (rank * workload), numofjob = startfrom + workload, array_size=0, j ;
    if(rank == size - 1){
        numofjob = numofjob + offset;
    }
    if(rank < size - 1){
        retcolor = malloc(iXmax * workload * sizeof(color)* sizeof(unsigned char));
    }

    if(rank != 0 && rank == size -1){
        retcolor = malloc(iXmax * (workload+offset) * sizeof(color)*sizeof(unsigned char));
    }

    for (iY = startfrom; iY < numofjob; iY++) {
        Cy = CyMin + (iY * PixelHeight);
        if (fabs(Cy) < (PixelHeight / 2)) {
            Cy = 0.0; /* Main antenna */
        }

        for (iX = 0; iX < iXmax; iX++) {
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
                color[0] = 0;
                color[1] = 0;
                color[2] = 0;
            } else {
                // Point outside the set. Mark it as white
                double c = 3 * log((double) Iteration) / log((double) (IterationMax) - 1.0);
                if (c < 1) {
                    color[0] = 0;
                    color[1] = 0;
                    color[2] = 255 * c;
                } else if (c < 2) {
                    color[0] = 0;
                    color[1] = 255 * (c - 1);
                    color[2] = 255;
                } else {
                    color[0] = 255 * (c - 2);
                    color[1] = 255;
                    color[2] = 255;
                }
            }
            if(rank != 0 ){
                for (j = 0; j < 3; j++){
                    retcolor[array_size + j] = color[j];
                }
                array_size += 3;
            } else{
                fwrite(color, 1, 3, fp);
            }
        }
    }
    if(rank != 0 && rank < size - 1) {
        MPI_Send(retcolor, iXmax * workload * sizeof(color), MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD);
    }
    else {
        if(rank != 0 && rank == size - 1){
            MPI_Send(retcolor, iXmax * (workload + offset) * sizeof(color), MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD);
        }
    }
    if(rank == 0){
        int i,k,element_offset = 0;
        unsigned char* recvcolor = malloc(iXmax * (iYmax-workload) * sizeof(color)*sizeof(unsigned char));
        for(i = 1; i < size; i++){
            if (i == size - 1){
                MPI_Recv(recvcolor + element_offset, iXmax * (workload + offset) * sizeof(color), MPI_UNSIGNED_CHAR, i, 0, MPI_COMM_WORLD, &status);
            }
            else {
                MPI_Recv(recvcolor + element_offset, iXmax * workload * sizeof(color), MPI_UNSIGNED_CHAR, i, 0, MPI_COMM_WORLD, &status);
                element_offset += iXmax * workload * sizeof(color);
            }
        }
        if(size > 1){
            for (k = 0; k < iXmax * (iYmax-workload) * sizeof(color); k+=3){
                color[0] = recvcolor[k];
                color[1] = recvcolor[k + 1];
                color[2] = recvcolor[k + 2];

                fwrite(&color, 1, 3, fp);
            }
        }

        fclose(fp);
    }

    end = MPI_Wtime();

    cpu_time_used = ((double)(end - start));

    MPI_Reduce(&cpu_time_used, &totaltime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if(rank == 0){
        printf("Completed Computing Mandelbrot Set.\n");
        printf("File: %s successfully closed.\n", filename);
        printf("Mandelbrot computational process time: %lf\n", totaltime);
    }

    MPI_Finalize();

    return 0;
}
