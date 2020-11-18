
/* Ensure C understands `bool` */
/* #include<stdbool.h> */
#define bool _Bool

float*
MLEM_TOF_Reco(int niterations,
              bool TOF,
              float TOF_resolution,
              float FOV_XY,
              float FOV_Z,
              int NXY,
              int NZ,
              int ncoinc,
              float * LOR_X1,
              float * LOR_Y1,
              float * LOR_Z1,
              float * LOR_T1,
              float * LOR_X2,
              float * LOR_Y2,
              float * LOR_Z2,
              float * LOR_T2,
              float * SENS,
              const char * outfile_prefix,
              int out_niter);

float ToFFunction(float dist, float deltaT, float TOF_resolution);
