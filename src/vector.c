// Copyright (c) 2015, UChicago Argonne, LLC. All rights reserved.

// Copyright 2015. UChicago Argonne, LLC. This software was produced 
// under U.S. Government contract DE-AC02-06CH11357 for Argonne National 
// Laboratory (ANL), which is operated by UChicago Argonne, LLC for the 
// U.S. Department of Energy. The U.S. Government has rights to use, 
// reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR 
// UChicago Argonne, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR 
// ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is 
// modified to produce derivative works, such modified software should 
// be clearly marked, so as not to confuse it with the version available 
// from ANL.

// Additionally, redistribution and use in source and binary forms, with 
// or without modification, are permitted provided that the following 
// conditions are met:

//     * Redistributions of source code must retain the above copyright 
//       notice, this list of conditions and the following disclaimer. 

//     * Redistributions in binary form must reproduce the above copyright 
//       notice, this list of conditions and the following disclaimer in 
//       the documentation and/or other materials provided with the 
//       distribution. 

//     * Neither the name of UChicago Argonne, LLC, Argonne National 
//       Laboratory, ANL, the U.S. Government, nor the names of its 
//       contributors may be used to endorse or promote products derived 
//       from this software without specific prior written permission. 

// THIS SOFTWARE IS PROVIDED BY UChicago Argonne, LLC AND CONTRIBUTORS 
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL UChicago 
// Argonne, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
// BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN 
// ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
// POSSIBILITY OF SUCH DAMAGE.

#include "utils.h"


void 
vector(
    const float *data1, int dy1, int dt1, int dx1,
    const float *data2, int dy2, int dt2, int dx2,
    const float *center1, const float *center2,
    const float *theta1, const float *theta2,
    float *recon1, float *recon2, float *recon3, 
    int ngridx, int ngridy, int num_iter)
{
    printf ("Hey!\n");
    printf ("dx=%d\n", dx1);

    float *gridx = (float *)malloc((ngridx+1)*sizeof(float));
    float *gridy = (float *)malloc((ngridy+1)*sizeof(float));
    float *coordx = (float *)malloc((ngridy+1)*sizeof(float));
    float *coordy = (float *)malloc((ngridx+1)*sizeof(float));
    float *ax = (float *)malloc((ngridx+ngridy)*sizeof(float));
    float *ay = (float *)malloc((ngridx+ngridy)*sizeof(float));
    float *bx = (float *)malloc((ngridx+ngridy)*sizeof(float));
    float *by = (float *)malloc((ngridx+ngridy)*sizeof(float));
    float *coorx = (float *)malloc((ngridx+ngridy)*sizeof(float));
    float *coory = (float *)malloc((ngridx+ngridy)*sizeof(float));
    float *dist = (float *)malloc((ngridx+ngridy)*sizeof(float));
    int *indi = (int *)malloc((ngridx+ngridy)*sizeof(int));

    assert(coordx != NULL && coordy != NULL &&
        ax != NULL && ay != NULL && by != NULL && bx != NULL &&
        coorx != NULL && coory != NULL && dist != NULL && indi != NULL);

    int s, p, d, i, n;
    int quadrant;
    float theta_p, sin_p, cos_p;
    float srcx, srcy, detx, dety, dv, vx, vy;
    float mov, xi, yi;
    int asize, bsize, csize;
    float *simdata;
    float upd;
    int ind_data, ind_recon;
    float *sum_dist;
    float sum_dist2;
    float *update;

    for (i=0; i<num_iter; i++) 
    {
        simdata = (float *)calloc((dt1*dy1*dx1), sizeof(float));

        // For each slice
        for (s=0; s<dy1; s++)
        {
            preprocessing(ngridx, ngridy, dx1, center1[s], 
                &mov, gridx, gridy); // Outputs: mov, gridx, gridy

            sum_dist = (float *)calloc((ngridx*ngridy), sizeof(float));
            update = (float *)calloc((ngridx*ngridy), sizeof(float));
            
            // For each projection angle 
            for (p=0; p<dt1; p++) 
            {
                // Calculate the sin and cos values 
                // of the projection angle and find
                // at which quadrant on the cartesian grid.
                theta_p = fmod(theta1[p], 2*M_PI);
                quadrant = calc_quadrant(theta_p);
                sin_p = sinf(theta_p);
                cos_p = cosf(theta_p);

                // For each detector pixel 
                for (d=0; d<dx1; d++) 
                {
                    // Calculate coordinates
                    xi = -ngridx-ngridy;
                    yi = (1-dx1)/2.0+d+mov;

                    srcx = xi*cos_p-yi*sin_p;
                    srcy = xi*sin_p+yi*cos_p;
                    detx = -xi*cos_p-yi*sin_p;
                    dety = -xi*sin_p+yi*cos_p;

                    dv = sqrt(pow(srcx-detx, 2)+pow(srcy-dety, 2));
                    vx = (srcx-detx)/dv;
                    vy = (srcy-dety)/dv;

                    calc_coords(
                        ngridx, ngridy, xi, yi, sin_p, cos_p, gridx, gridy, 
                        coordx, coordy);

                    // Merge the (coordx, gridy) and (gridx, coordy)
                    trim_coords(
                        ngridx, ngridy, coordx, coordy, gridx, gridy, 
                        &asize, ax, ay, &bsize, bx, by);

                    // Sort the array of intersection points (ax, ay) and
                    // (bx, by). The new sorted intersection points are 
                    // stored in (coorx, coory). Total number of points 
                    // are csize.
                    sort_intersections(
                        quadrant, asize, ax, ay, bsize, bx, by, 
                        &csize, coorx, coory);

                    // Calculate the distances (dist) between the 
                    // intersection points (coorx, coory). Find the 
                    // indices of the pixels on the reconstruction grid.
                    calc_dist(
                        ngridx, ngridy, csize, coorx, coory, 
                        indi, dist);

                    // Calculate simdata 
                    calc_simdata2(s, p, d, ngridx, ngridy, dt1, dx1,
                        csize, indi, dist, vx, vy, recon1, recon2,
                        simdata); // Output: simulated data

                    // // Calculate simdata 
                    // calc_simdata(s, p, d, ngridx, ngridy, dt1, dx1,
                    //     csize, indi, dist, recon1,
                    //     simdata); // Output: simdata

                    // Calculate dist*dist
                    sum_dist2 = 0.0;
                    for (n=0; n<csize-1; n++) 
                    {
                        sum_dist2 += dist[n]*dist[n];
                        sum_dist[indi[n]] += dist[n];
                    }

                    // Update
                    if (sum_dist2 != 0.0) 
                    {
                        ind_data = d+p*dx1+s*dt1*dx1;
                        upd = (data1[ind_data]-simdata[ind_data])/sum_dist2;
                        for (n=0; n<csize-1; n++) 
                        {
                            update[indi[n]] += upd*dist[n];
                        }
                    }
                }
            }

            for (n = 0; n < ngridx*ngridy; n++) {
                if (sum_dist[n] != 0.0) {
                    ind_recon = s*ngridx*ngridy;
                    recon1[n+ind_recon] += update[n]/sum_dist[n];
                    recon2[n+ind_recon] += update[n]/sum_dist[n];
                }
            }

            free(sum_dist);
            free(update);
        }

        free(simdata);



        simdata = (float *)calloc((dt2*dy2*dx2), sizeof(float));

        // For each slice
        for (s=0; s<dy2; s++)
        {
            preprocessing(ngridx, ngridy, dx2, center1[s], 
                &mov, gridx, gridy); // Outputs: mov, gridx, gridy

            sum_dist = (float *)calloc((ngridx*ngridy), sizeof(float));
            update = (float *)calloc((ngridx*ngridy), sizeof(float));
            
            // For each projection angle 
            for (p=0; p<dt2; p++) 
            {
                // Calculate the sin and cos values 
                // of the projection angle and find
                // at which quadrant on the cartesian grid.
                theta_p = fmod(theta2[p], 2*M_PI);
                quadrant = calc_quadrant(theta_p);
                sin_p = sinf(theta_p);
                cos_p = cosf(theta_p);

                // For each detector pixel 
                for (d=0; d<dx2; d++) 
                {
                    // Calculate coordinates
                    xi = -ngridx-ngridy;
                    yi = (1-dx2)/2.0+d+mov;

                    srcx = xi*cos_p-yi*sin_p;
                    srcy = xi*sin_p+yi*cos_p;
                    detx = -xi*cos_p-yi*sin_p;
                    dety = -xi*sin_p+yi*cos_p;

                    dv = sqrt(pow(srcx-detx, 2)+pow(srcy-dety, 2));
                    vx = (srcx-detx)/dv;
                    vy = (srcy-dety)/dv;

                    calc_coords(
                        ngridx, ngridy, xi, yi, sin_p, cos_p, gridx, gridy, 
                        coordx, coordy);

                    // Merge the (coordx, gridy) and (gridx, coordy)
                    trim_coords(
                        ngridx, ngridy, coordx, coordy, gridx, gridy, 
                        &asize, ax, ay, &bsize, bx, by);

                    // Sort the array of intersection points (ax, ay) and
                    // (bx, by). The new sorted intersection points are 
                    // stored in (coorx, coory). Total number of points 
                    // are csize.
                    sort_intersections(
                        quadrant, asize, ax, ay, bsize, bx, by, 
                        &csize, coorx, coory);

                    // Calculate the distances (dist) between the 
                    // intersection points (coorx, coory). Find the 
                    // indices of the pixels on the reconstruction grid.
                    calc_dist(
                        ngridx, ngridy, csize, coorx, coory, 
                        indi, dist);

                    // Calculate simdata 
                    calc_simdata2(s, p, d, ngridx, ngridy, dt2, dx2,
                        csize, indi, dist, vx, vy, recon1, recon3,
                        simdata); // Output: simulated data

                    // // Calculate simdata 
                    // calc_simdata(s, p, d, ngridx, ngridy, dt1, dx1,
                    //     csize, indi, dist, recon1,
                    //     simdata); // Output: simdata

                    // Calculate dist*dist
                    sum_dist2 = 0.0;
                    for (n=0; n<csize-1; n++) 
                    {
                        sum_dist2 += dist[n]*dist[n];
                        sum_dist[indi[n]] += dist[n];
                    }

                    // Update
                    if (sum_dist2 != 0.0) 
                    {
                        ind_data = d+p*dx2+s*dt2*dx2;
                        upd = (data2[ind_data]-simdata[ind_data])/sum_dist2;
                        for (n=0; n<csize-1; n++) 
                        {
                            update[indi[n]] += upd*dist[n];
                        }
                    }
                }
            }

            for (n = 0; n < ngridx*ngridy; n++) {
                if (sum_dist[n] != 0.0) {
                    ind_recon = s*ngridx*ngridy;
                    recon1[n+ind_recon] += update[n]/sum_dist[n];
                    recon3[n+ind_recon] += update[n]/sum_dist[n];
                }
            }

            free(sum_dist);
            free(update);
        }

        free(simdata);

    }

    free(gridx);
    free(gridy);
    free(coordx);
    free(coordy);
    free(ax);
    free(ay);
    free(bx);
    free(by);
    free(coorx);
    free(coory);
    free(dist);
    free(indi);
}
