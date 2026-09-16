/************************************************************
read3B4XRT.c
Example program to read TMPA-RT (3B4xRT) Level-3 3 hour products.

First, this program calls either read3B40RT, read3B41RT,
or read3B42RT.  Each of these functions does the following:

1. Reads the header and checks that the file is of
   the expected type (3B40RT, 3B41RT, or 3B42RT).

2. Reads and unscales the precipitation and error data.

3. Reads the other data fields with no unscaling, except for
   the uncalibrated precipitation.

4. Checks that there is no unread data in the file.

Second, this program dumps the header and prints out
values for the 4 points around the intersection
of the dateline and the equator.

CHANGE LOG
  DATE            PROGRAMMER                  CHANGE DESCRIPTION
  ----            -----------                 ------------------
  2/20/02         J. Stout/GMU                Original version
  2/3/09          G.J. Huffman/SSAI           Upgrade for next
                                              version Ð more fields
************************************************************/


#include <stdio.h>
#include <stdlib.h>
#define LONMAX 1440
#define LATMAX90 720
#define LATMAX60 480
#define I2MISSING -31999
#define FMISSING -99999.
#define SCALE 100

int read3B40RT (
    char *filename,
    char header[],
    float precip[][LONMAX],
    float error[][LONMAX],
    char npixel[][LONMAX],
    char nambig[][LONMAX],
    char nrain[][LONMAX],
    char source[][LONMAX]);
int read3B41RT (
    char *filename,
    char header[],
    float precip[][LONMAX],
    float error[][LONMAX],
    char npixel[][LONMAX]);
int read3B42RT (
    char *filename,
    char header[],
    float precip[][LONMAX],
    float error[][LONMAX],
    char source[][LONMAX],
    float p_uncal[][LONMAX]);
int check3B4Xheader (char alg[], int nchar, char header[]);
int check_at_EOF (FILE *fp);
void dump3B4XRTheader (int nchar, char header[]);
void print3B4XRTdata (
    char algtype[],
    int latmax,
    float p[][LONMAX],
    float e[][LONMAX],
    char npixels[][LONMAX],
    char nambig[][LONMAX],
    char nrain[][LONMAX],
    char source[][LONMAX],
    float p_uncal[][LONMAX]);
int max1Barray (int latmax, char a[][LONMAX]);
int read3B4XRT2B (
    int latmax,
    FILE *fp,
    float buffer[][LONMAX]);
int read3B4XRT1B (
    int latmax,
    FILE *fp,
    char buffer[][LONMAX]);

void main (int argc, char *argv[])
{
    int nhchar;
    char header[LONMAX*2];
    float precip90[LATMAX90][LONMAX];
    float error90[LATMAX90][LONMAX];
    char npixels90[LATMAX90][LONMAX];
    char n_ambig_pix90[LATMAX90][LONMAX];
    char n_rain_pix90[LATMAX90][LONMAX];
    char source90[LATMAX90][LONMAX];
    float precip60[LATMAX60][LONMAX];
    float error60[LATMAX60][LONMAX];
    char npixels60[LATMAX60][LONMAX];
    char source60[LATMAX60][LONMAX];
    float p_uncal60[LATMAX60][LONMAX];
    int status;
    nhchar = LONMAX*2;

    if (argc < 2)
    {
        printf("USAGE: read3B4XRT <type> <file>\n");
        printf("       <type> = algorithm type (3B40RT,3B41RT,3B42RT)\n");
        printf("       <file> = file name (3B42RT.2001121809.bin)\n");
        exit (1);
    }

    if (strcmp(argv[1], "3B40RT") == 0)
    {
        status = read3B40RT (argv[2], header, precip90, error90,
          npixels90, n_ambig_pix90, n_rain_pix90, source90);
        if (status == 0)
        {
            dump3B4XRTheader (nhchar, header);
            print3B4XRTdata ("3B40RT", LATMAX90, precip90, error90,
              npixels90, n_ambig_pix90, n_rain_pix90, source90,
              p_uncal60);
        }
    }
    else if (strcmp(argv[1], "3B41RT") == 0)
    {
        status = read3B41RT (argv[2], header, precip60, error60,
          npixels60);
        if (status == 0)
        {
            dump3B4XRTheader (nhchar, header);
            print3B4XRTdata ("3B41RT", LATMAX60, precip60, error60,
              npixels60, n_ambig_pix90, n_rain_pix90, source60,
              p_uncal60);
        }
    }
    else if (strcmp(argv[1], "3B42RT") == 0)
    {
        status = read3B42RT (argv[2], header, precip60, error60,
          source60, p_uncal60);
        if (status == 0)
        {
            dump3B4XRTheader (nhchar, header);
            print3B4XRTdata ("3B42RT", LATMAX60, precip60, error60,
              npixels60, n_ambig_pix90, n_rain_pix90, source60,
              p_uncal60);
        }
    }
    else
    {
        printf("ERROR in read_rt_file: bad algor. type: %s\n",argv[1]);
        printf("     it must be one of 3B40RT,3B41RT,3B42RT\n");
        exit (1);
    }
}


/************************************************************
read3B40RT

DESCRIPTION
    Reads a 3B40RT file from the TRMM realtime system

INPUT
    filename    Name of the 3B40RT file

OUTPUT
    header      The header
    precip      Precipitation
    error       Error of precipitation
    npixel      Number of pixels
    nambig      Number of ambiguous pixels
    nrain       Number of raining pixels
    source      Sensor type that provided estimate

RETURN VALUES
    0           Success
    -1          Error - could not open file
    -2          Error - file is wrong algorithm
    -3          Error - file smaller than expected
    -4          Error - file larger than expected
    -5          Error - could not close file

************************************************************/
int read3B40RT (
    char *filename,
    char header[],
    float precip[][LONMAX],
    float error[][LONMAX],
    char npixel[][LONMAX],
    char nambig[][LONMAX],
    char nrain[][LONMAX],
    char source[][LONMAX])
{
    FILE *fp;
    int status;
    int nread;

    fp = fopen (filename, "r");
    if (fp == NULL)
    {
        printf ("ERROR: could not open %s\n",filename);
        return -1;
    }

    nread = fread (header, 1, LONMAX*2, fp);
    status = check3B4Xheader ("3B40RT", LONMAX*2, header);
    if (status < 0)
    {
        printf ("ERROR: %s is not 3B40RT\n",filename);
        return -2;
    }
    nread = read3B4XRT2B (LATMAX90, fp, precip);
    nread = read3B4XRT2B (LATMAX90, fp, error);
    nread = read3B4XRT1B (LATMAX90, fp, npixel);
    nread = read3B4XRT1B (LATMAX90, fp, nambig);
    nread = read3B4XRT1B (LATMAX90, fp, nrain);
    nread = read3B4XRT1B (LATMAX90, fp, source);
    if (nread < 1)
    {
        printf ("ERROR: %s is smaller than expected\n",filename);
        return -3;
    }
    status = check_at_EOF (fp);
    if (status < 0)
    {
        printf ("ERROR: %s is larger than expected\n",filename);
        return -4;
    }
    status = fclose (fp);
    if (status < 0)
    {
        printf ("ERROR closing %s\n",filename);
        return -5;
    }
    return 0;
}

/************************************************************
read3B41RT

DESCRIPTION
    Reads a 3B41RT file from the TRMM realtime system

INPUT
    filename    Name of the 3B41RT file

OUTPUT
    header      The header
    precip      Precipitation
    error       Error of precipitation
    npixel      Number of pixels

RETURN VALUES
    0           Success
    -1          Error - could not open file
    -2          Error - file is wrong algorithm
    -3          Error - file smaller than expected
    -4          Error - file larger than expected
    -5          Error - could not close file

************************************************************/
int read3B41RT (
    char *filename,
    char header[],
    float precip[][LONMAX],
    float error[][LONMAX],
    char npixel[][LONMAX])
{
    FILE *fp;
    int status;
    int nread;

    fp = fopen (filename, "r");
    if (fp == NULL)
    {
        printf ("ERROR: could not open %s\n",filename);
        return -1;
    }

    nread = fread (header, 1, LONMAX*2, fp);
    status = check3B4Xheader ("3B41RT", LONMAX*2, header);
    if (status < 0)
    {
        printf ("ERROR: %s is not 3B40RT\n",filename);
        return -2;
    }
    nread = read3B4XRT2B (LATMAX60, fp, precip);
    nread = read3B4XRT2B (LATMAX60, fp, error);
    nread = read3B4XRT1B (LATMAX60, fp, npixel);
    if (nread < 1)
    {
        printf ("ERROR: %s is smaller than expected\n",filename);
        return -3;
    }
    status = check_at_EOF (fp);
    if (status < 0)
    {
        printf ("ERROR: %s is larger than expected\n",filename);
        return -4;
    }
    status = fclose (fp);
    if (status < 0)
    {
        printf ("ERROR closing %s\n",filename);
        return -5;
    }
    return 0;
}

/************************************************************
read3B42RT

DESCRIPTION
    Reads a 3B42RT file from the TRMM realtime system

INPUT
    filename    Name of the 3B42RT file

OUTPUT
    header      The header
    precip      Precipitation
    error       Error of precipitation
    source      Data type used for grid box
    p_uncal     Precipitation before calibration to 3B42 V.6

RETURN VALUES
    0           Success
    -1          Error - could not open file
    -2          Error - file is wrong algorithm
    -3          Error - file smaller than expected
    -4          Error - file larger than expected
    -5          Error - could not close file

************************************************************/
int read3B42RT (
    char *filename,
    char header[],
    float precip[][LONMAX],
    float error[][LONMAX],
    char source[][LONMAX],
    float p_uncal[][LONMAX])
{
    FILE *fp;
    int status;
    int nread;

    fp = fopen (filename, "r");
    if (fp == NULL)
    {
        printf ("ERROR: could not open %s\n",filename);
        return -1;
    }

    nread = fread (header, 1, LONMAX*2, fp);
    status = check3B4Xheader ("3B42RT", LONMAX*2, header);
    if (status < 0)
    {
        printf ("ERROR: %s is not 3B40RT\n",filename);
        return -2;
    }
    nread = read3B4XRT2B (LATMAX60, fp, precip);
    nread = read3B4XRT2B (LATMAX60, fp, error);
    nread = read3B4XRT1B (LATMAX60, fp, source);
    nread = read3B4XRT2B (LATMAX60, fp, p_uncal);
    if (nread < 1)
    {
        printf ("ERROR: %s is smaller than expected\n",filename);
        return -3;
    }
    status = check_at_EOF (fp);
    if (status < 0)
    {
        printf ("ERROR: %s is larger than expected\n",filename);
        return -4;
    }
    status = fclose (fp);
    if (status < 0)
    {
        printf ("ERROR closing %s\n",filename);
        return -5;
    }
    return 0;
}

void print3B4XRTdata (
    char algtype[],
    int latmax,
    float precip[][LONMAX],
    float error[][LONMAX],
    char npixel[][LONMAX],
    char nambig[][LONMAX],
    char nrain[][LONMAX],
    char source[][LONMAX],
    float p_uncal[][LONMAX])
{
    int i, j, i1, j1;

    i1 = LONMAX / 2;
    j1 = latmax / 2;
    printf("\nVALUES FOR THE 4 POINTS AROUND THE INTERSECTION\n");
    printf("OF THE DATELINE AND THE EQUATOR:\n");
    if (strcmp(algtype, "3B40RT") == 0)
    {
        printf("\nilon, ilat, precip, error, npixel, nambig, nrain, source\n");
        for (j = j1-1; j < j1+1; j++)
        for (i = i1-1; i < i1+1; i++)
            printf("%i %i %f %f %i %i %i %i\n", j, i,
            precip[j][i], error[j][i], npixel[j][i],
            nambig[j][i], nrain[j][i], source[j][i]);
    }
    else if (strcmp(algtype, "3B41RT") == 0)
    {
        printf("\nilon, ilat, precip, error, npixel\n");
        for (j = j1-1; j < j1+1; j++)
        for (i = i1-1; i < i1+1; i++)
            printf("%i %i %f %f %i\n", j, i,
            precip[j][i], error[j][i], npixel[j][i]);
    }
    else if (strcmp(algtype, "3B42RT") == 0)
    {
        printf("\nilon, ilat, precip, error, source, p_uncal\n");
        for (j = j1-1; j < j1+1; j++)
        for (i = i1-1; i < i1+1; i++)
            printf("%i %i %f %f %i %f\n", j, i,
            precip[j][i], error[j][i], npixel[j][i], p_uncal[j][i]);
    }
}


int max1Barray (int latmax, char a[][LONMAX])
{
    int i, j;
    int max;
    max = -9999;
    for (i=0; i<latmax; i++) for (j=0; j<LONMAX; j++)
        if (a[i][j] > max) max = a[i][j];
    return max;
}

int check3B4Xheader (char alg[], int nchar, char header[])
/* returns 0 if alg found in the right place in header
   returns -1 if no equal sign found in header
   returns -2 if characters after first equal sign do not match alg */
{
    int i, j, status;

    /* find first equal sign*/
    for (i=0; i<nchar; i++)
    {
        if (header[i] == '=') break;
    }
    if (i == nchar-1)
    {
        printf("ERROR: no equal sign in header\n");
        return -1;
    }

    /* advance beyond any blanks */
    for (i=i; i<nchar; i++)
    {
        printf ("i,header[i+1] = %i %c\n",i,header[i+1]);
        if (header[i+1] != ' ') break;
    }
    if (i == nchar-1)
    {
        printf("ERROR: nothing but blanks after first equal sign in header\n");
        return -1;
    }

    /* compare 6 characters after equal sign to expected alg */
    status = 0;
    for (j=0; j<6; j++)
    {
        if (header[i+j+1] != alg[j])
        {
            printf("ERROR: did not find %s in header\n",alg);
            status = -2;
        }
    }
    return status;
}

int check_at_EOF (FILE *fp)
{
    int c;
    /* read one more byte, should be EOF */
    c = getc(fp);
    if (c != -1)
    {
        printf("ERROR: unread data in file.  First unread char = %i\n",c);
        return -2;
    }
    return 0;
}

void dump3B4XRTheader (int nchar, char header[])
{
    int i;
    printf("\nHeader:\n");
    for (i=0; i<nchar; i++)
    {
        putchar(header[i]);
    }
    putchar('\n');
}

int read3B4XRT2B (
    int latmax,
    FILE *fp,
    float out[][LONMAX])
{
    int row, col;
    int nread;
    short int buffer[LONMAX];
    for (row=0; row<latmax; row++)
    {
        nread = fread (buffer, 2, LONMAX, fp);
        for (col=0; col<LONMAX; col++)
        {
            if (buffer[col] == I2MISSING)
                out[row][col] = FMISSING;
            else
                out[row][col] = buffer[col] / SCALE;
        }
    }
    return nread;
}

int read3B4XRT1B (
    int latmax,
    FILE *fp,
    char buffer[][LONMAX])
{
    int row;
    int nread;
    for (row=0; row<latmax; row++)
    {
        nread = fread (buffer[row], 1, LONMAX, fp);
    }
    return nread;
}

