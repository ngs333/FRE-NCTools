/***********************************************************************
 *                   GNU Lesser General Public License
 *
 * This file is part of the GFDL FRE NetCDF tools package (FRE-NCTools).
 *
 * FRE-NCtools is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at
 * your option) any later version.
 *
 * FRE-NCtools is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with FRE-NCTools.  If not, see
 * <http://www.gnu.org/licenses/>.
 **********************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "globals.h"
#include "mosaic_util.h"
#include "create_xgrid.h"
#include "create_xgrid_util.h"
#include "create_xgrid_acc.h"
#include "constant.h"

#define AREA_RATIO_THRESH (1.e-6)
#define MASK_THRESH       (0.5)
#define EPSLN8            (1.e-8)
#define EPSLN30           (1.0e-30)
#define EPSLN10           (1.0e-10)
#define MAX_V 8

/*******************************************************************************
  prepare_create_xgrid_2dx2d_order2_acc
*******************************************************************************/
int prepare_create_xgrid_2dx2d_order2_acc(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
                                          const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                                          Minmaxavg_lists *out_minmaxavg_lists, const double *mask_in,
                                          int *counts_per_ij1, int *ij2_start, int *ij2_end)
{

#define MAX_V 8
  int nx1, nx2, ny1, ny2, nx1p, nx2p;
  size_t approx_nxgrid;

  nx1 = *nlon_in;
  ny1 = *nlat_in;
  nx2 = *nlon_out;
  ny2 = *nlat_out;
  nx1p = nx1 + 1;
  nx2p = nx2 + 1;

  approx_nxgrid = 0;

#pragma acc data present(lon_out[0:(nx2+1)*(ny2+1)], lat_out[0:(nx2+1)*(ny2+1)])
#pragma acc data present(lon_in[0:(nx1+1)*(ny1+1)], lat_in[0:(nx1+1)*(ny1+1)], mask_in[0:nx1*ny1])
#pragma acc data present(out_minmaxavg_lists->lon_list[0:MAX_V*nx2*ny2], out_minmaxavg_lists->lat_list[0:MAX_V*nx2*ny2])
#pragma acc data present(out_minmaxavg_lists->n_list[0:nx2*ny2], out_minmaxavg_lists->lon_avg[0:nx2*ny2])
#pragma acc data present(out_minmaxavg_lists->lat_min_list[0:nx2*ny2], out_minmaxavg_lists->lat_max_list[0:nx2*ny2])
#pragma acc data present(out_minmaxavg_lists->lon_min_list[0:nx2*ny2], out_minmaxavg_lists->lon_max_list[0:nx2*ny2])
#pragma acc data present(counts_per_ij1[0:nx1*ny1], ij2_start[0:nx1*ny1], ij2_end[0:nx1*ny1])
#pragma acc data copy(approx_nxgrid)
#pragma acc parallel
{
#pragma acc loop independent reduction(+:approx_nxgrid)
  for( int ij1=0 ; ij1 < nx1*ny1 ; ij1++) {

    int i1, j1;
    int icount=0 ;
    int ij2_max=0 , ij2_min=nx2*ny2+1;

    i1 = ij1%nx1;
    j1 = ij1/nx1;

    counts_per_ij1[ij1]=0;

    if( mask_in[ij1] > MASK_THRESH ) {

      int n0, n1, n2, n3, n1_in;
      double lat_in_min, lat_in_max, lon_in_min, lon_in_max, lon_in_avg;
      double x1_in[MV], y1_in[MV];

      n0 = j1*nx1p+i1;       n1 = j1*nx1p+i1+1;
      n2 = (j1+1)*nx1p+i1+1; n3 = (j1+1)*nx1p+i1;
      x1_in[0] = lon_in[n0]; y1_in[0] = lat_in[n0];
      x1_in[1] = lon_in[n1]; y1_in[1] = lat_in[n1];
      x1_in[2] = lon_in[n2]; y1_in[2] = lat_in[n2];
      x1_in[3] = lon_in[n3]; y1_in[3] = lat_in[n3];
      lat_in_min = minval_double(4, y1_in);
      lat_in_max = maxval_double(4, y1_in);
      n1_in = fix_lon(x1_in, y1_in, 4, M_PI);
      lon_in_min = minval_double(n1_in, x1_in);
      lon_in_max = maxval_double(n1_in, x1_in);
      lon_in_avg = avgval_double(n1_in, x1_in);

#pragma acc loop independent reduction(+:approx_nxgrid) reduction(+:icount) reduction(min:ij2_min) reduction(max:ij2_max)
      for(int ij2=0; ij2<nx2*ny2; ij2++) {

        int i2, j2, l;
        double dx, lon_out_min, lon_out_max;
        double x2_in[MAX_V], y2_in[MAX_V],  x_out[MV], y_out[MV];;

        i2 = ij2%nx2;
        j2 = ij2/nx2;

        if(out_minmaxavg_lists->lat_min_list[ij2] >= lat_in_max || out_minmaxavg_lists->lat_max_list[ij2] <= lat_in_min ) continue;

        /* adjust x2_in according to lon_in_avg*/
        lon_out_min = out_minmaxavg_lists->lon_min_list[ij2];
        lon_out_max = out_minmaxavg_lists->lon_max_list[ij2];
        dx = out_minmaxavg_lists->lon_avg[ij2] - lon_in_avg;

        if(dx < -M_PI ) {
          lon_out_min += TPI;
          lon_out_max += TPI;
        }

        else if (dx >  M_PI) {
          lon_out_min -= TPI;
          lon_out_max -= TPI;
        }

        /* x2_in should in the same range as x1_in after lon_fix, so no need to consider cyclic condition */
        if(lon_out_min >= lon_in_max || lon_out_max <= lon_in_min ) continue;


        //Note, the check for AREA_RATIO_THRESH has been removed
        //Thus, the computed value of approx_nxgrid will be equal to or greater than nxgrid
        approx_nxgrid++;
        icount++;
        ij2_min = min(ij2_min, ij2);
        ij2_max = max(ij2_max, ij2);

      } //ij2

      counts_per_ij1[ij1] = icount;
      ij2_start[ij1] = ij2_min ;
      ij2_end[ij1] = ij2_max;

    } // mask
  } //ij1
} //kernel

  return approx_nxgrid;

}


/*******************************************************************************
  create_xgrid_2dx2d_order2 OPENACC version
*******************************************************************************/
int create_xgrid_2dx2d_order2_acc(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
                                  const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                                  Minmaxavg_lists *out_minmaxavg_lists, const double *mask_in, const int approx_nxgrid,
                                  const int *counts_per_ij1, const int *ij2_start, const int *ij2_end,
                                  int **i_in, int **j_in, int **i_out, int **j_out,
                                  double **xgrid_area, double **xgrid_clon, double **xgrid_clat)
{

#define MAX_V 8
  int nx1, nx2, ny1, ny2, nx1p, nx2p;
  double *area_in, *area_out;

  int ixgrid2;
  int *i_in2, *j_in2, *i_out2, *j_out2 ;
  double *xgrid_area2, *xgrid_clon2, *xgrid_clat2;

  size_t nxgrid;

  nx1 = *nlon_in;
  ny1 = *nlat_in;
  nx2 = *nlon_out;
  ny2 = *nlat_out;
  nx1p = nx1 + 1;
  nx2p = nx2 + 1;

  //Temporarily holds information about exchange grid cells.
  //Not all elements of these arrays will be filled in because
  //approx_nxgrid >= nxgrid
  i_in2 = (int *)malloc(approx_nxgrid*sizeof(int));
  j_in2 = (int *)malloc(approx_nxgrid*sizeof(int));
  i_out2 = (int *)malloc(approx_nxgrid*sizeof(int));
  j_out2 = (int *)malloc(approx_nxgrid*sizeof(int));
  xgrid_area2 = (double *)malloc(approx_nxgrid*sizeof(double));
  xgrid_clon2 = (double *)malloc(approx_nxgrid*sizeof(double));
  xgrid_clat2 = (double *)malloc(approx_nxgrid*sizeof(double));
  for(int i=0 ; i<approx_nxgrid ; i++) i_in2[i]=-99;

  area_in  = (double *)malloc(nx1*ny1*sizeof(double));
  area_out = (double *)malloc(nx2*ny2*sizeof(double));
  get_grid_area(nlon_in, nlat_in, lon_in, lat_in, area_in);
  get_grid_area(nlon_out, nlat_out, lon_out, lat_out, area_out);

  nxgrid = 0;

#pragma acc data present(lon_out[0:(nx2+1)*(ny2+1)], lat_out[0:(nx2+1)*(ny2+1)])
#pragma acc data present(lon_in[0:(nx1+1)*(ny1+1)], lat_in[0:(nx1+1)*(ny1+1)], mask_in[0:nx1*ny1])
#pragma acc data present(out_minmaxavg_lists->lon_list[0:MAX_V*nx2*ny2], out_minmaxavg_lists->lat_list[0:MAX_V*nx2*ny2])
#pragma acc data present(out_minmaxavg_lists->n_list[0:nx2*ny2], out_minmaxavg_lists->lon_avg[0:nx2*ny2])
#pragma acc data present(out_minmaxavg_lists->lat_min_list[0:nx2*ny2], out_minmaxavg_lists->lat_max_list[0:nx2*ny2])
#pragma acc data present(out_minmaxavg_lists->lon_min_list[0:nx2*ny2], out_minmaxavg_lists->lon_max_list[0:nx2*ny2])
#pragma acc data present(counts_per_ij1[0:nx1*ny1], ij2_start[0:nx1*ny1], ij2_end[0:nx1*ny1])
#pragma acc data copyout(xgrid_area2[0:approx_nxgrid], xgrid_clon2[0:approx_nxgrid], xgrid_clat2[0:approx_nxgrid], \
                         j_in2[0:approx_nxgrid], j_out2[0:approx_nxgrid], i_out2[0:approx_nxgrid])
#pragma acc data copyin(area_in[0:nx1*ny1], area_out[0:nx2*ny2])
#pragma acc data copy(nxgrid, i_in2[0:approx_nxgrid])
#pragma acc parallel
{
#pragma acc loop independent reduction(+:nxgrid)
  for( int ij1=0 ; ij1<nx1*ny1 ; ij1++) {
    if( mask_in[ij1] > MASK_THRESH ) {

      int n0, n1, n2, n3, n1_in;
      int i1, j1;
      int ij1_start, ixgrid;
      double lat_in_min, lat_in_max, lon_in_min, lon_in_max, lon_in_avg;
      double x1_in[MV], y1_in[MV];

      i1 = ij1%nx1;
      j1 = ij1/nx1;

      n0 = j1*nx1p+i1;       n1 = j1*nx1p+i1+1;
      n2 = (j1+1)*nx1p+i1+1; n3 = (j1+1)*nx1p+i1;
      x1_in[0] = lon_in[n0]; y1_in[0] = lat_in[n0];
      x1_in[1] = lon_in[n1]; y1_in[1] = lat_in[n1];
      x1_in[2] = lon_in[n2]; y1_in[2] = lat_in[n2];
      x1_in[3] = lon_in[n3]; y1_in[3] = lat_in[n3];
      lat_in_min = minval_double(4, y1_in);
      lat_in_max = maxval_double(4, y1_in);

      n1_in = fix_lon(x1_in, y1_in, 4, M_PI);
      lon_in_min = minval_double(n1_in, x1_in);
      lon_in_max = maxval_double(n1_in, x1_in);
      lon_in_avg = avgval_double(n1_in, x1_in);

      ixgrid=0;
      // ij1_start, the total number of exchange grid cells computed for input cell ij1
      // is an approximation.
      ij1_start=0;
      if(ij1>0) {
#pragma acc loop seq
        for(int i=0 ; i<ij1 ; i++) ij1_start+=counts_per_ij1[i];
      }

#pragma acc loop seq reduction(+:nxgrid)
      for(int ij2=ij2_start[ij1]; ij2<=ij2_end[ij1]; ij2++) {

        int n_out, i2, j2, n2_in, l;
        double xarea, dx, lon_out_min, lon_out_max;
        double x2_in[MAX_V], y2_in[MAX_V],  x_out[MV], y_out[MV];;

        if(out_minmaxavg_lists->lat_min_list[ij2] >= lat_in_max || out_minmaxavg_lists->lat_max_list[ij2] <= lat_in_min ) continue;

        i2 = ij2%nx2;
        j2 = ij2/nx2;

        /* adjust x2_in according to lon_in_avg*/
        n2_in = out_minmaxavg_lists->n_list[ij2];
#pragma acc loop seq
        for(l=0; l<n2_in; l++) {
          x2_in[l] = out_minmaxavg_lists->lon_list[ij2*MAX_V+l];
          y2_in[l] = out_minmaxavg_lists->lat_list[ij2*MAX_V+l];
        }
        lon_out_min = out_minmaxavg_lists->lon_min_list[ij2];
        lon_out_max = out_minmaxavg_lists->lon_max_list[ij2];
        dx = out_minmaxavg_lists->lon_avg[ij2] - lon_in_avg;
        if(dx < -M_PI ) {
          lon_out_min += TPI;
          lon_out_max += TPI;
#pragma acc loop seq
          for (l=0; l<n2_in; l++) x2_in[l] += TPI;
        }
        else if (dx >  M_PI) {
          lon_out_min -= TPI;
          lon_out_max -= TPI;
#pragma acc loop seq
          for (l=0; l<n2_in; l++) x2_in[l] -= TPI;
        }

        n_out = 1;
        if (  (n_out = clip_2dx2d( x1_in, y1_in, n1_in, x2_in, y2_in, n2_in, x_out, y_out )) > 0) {
          double min_area;
          xarea = poly_area (x_out, y_out, n_out ) * mask_in[ij1];
          min_area = min(area_in[ij1], area_out[ij2]);
          if( xarea/min_area > AREA_RATIO_THRESH ) {
            xgrid_area2[ij1_start+ixgrid] = xarea;
            xgrid_clon2[ij1_start+ixgrid] = poly_ctrlon(x_out, y_out, n_out, lon_in_avg);
            xgrid_clat2[ij1_start+ixgrid] = poly_ctrlat (x_out, y_out, n_out );
            i_in2[ij1_start+ixgrid] = i1;
            j_in2[ij1_start+ixgrid] = j1;
            i_out2[ij1_start+ixgrid] = i2;
            j_out2[ij1_start+ixgrid] = j2;
            ixgrid++;
            nxgrid++;
          } //if
        } //if
      } //ij2
    } //mask
  } //ij1
} //kernel

 free(area_in);
 free(area_out);

 // record data
 *i_in=(int *)malloc(nxgrid*sizeof(int));
 *j_in=(int *)malloc(nxgrid*sizeof(int));
 *i_out=(int *)malloc(nxgrid*sizeof(int));
 *j_out=(int *)malloc(nxgrid*sizeof(int));
 *xgrid_area=(double *)malloc(nxgrid*sizeof(double));
 *xgrid_clon=(double *)malloc(nxgrid*sizeof(double));
 *xgrid_clat=(double *)malloc(nxgrid*sizeof(double));

 ixgrid2=0;
 for(int i=0 ; i<approx_nxgrid ; i++) {
   if( i_in2[i]!=-99) {
     (*i_in)[ixgrid2] = i_in2[i];
     (*j_in)[ixgrid2] = j_in2[i];
     (*i_out)[ixgrid2] = i_out2[i];
     (*j_out)[ixgrid2] = j_out2[i];
     (*xgrid_area)[ixgrid2] = xgrid_area2[i];
     (*xgrid_clon)[ixgrid2] = xgrid_clon2[i];
     (*xgrid_clat)[ixgrid2] = xgrid_clat2[i];
     ixgrid2++;
   }
 }

  return nxgrid;

};/* get_xgrid_2Dx2D_order2 */