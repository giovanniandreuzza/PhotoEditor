/*
 Created by Sebastiano Vascon on 23/03/20.
*/

#include <stdio.h>
#include "ip_lib.h"
#include "bmp.h"

void ip_mat_show(ip_mat * t){
    unsigned int i,l,j;
    printf("Matrix of size %d x %d x %d (hxwxk)\n",t->h,t->w,t->k);
    for (l = 0; l < t->k; l++) {
        printf("Slice %d\n", l);
        for(i=0;i<t->h;i++) {
            for (j = 0; j < t->w; j++) {
                printf("%f ", get_val(t,i,j,l));
            }
            printf("\n");
        }
        printf("\n");
    }
}

void ip_mat_show_stats(ip_mat * t){
    unsigned int k;

    compute_stats(t);

    for(k=0;k<t->k;k++){
        printf("Channel %d:\n", k);
        printf("\t Min: %f\n", t->stat[k].min);
        printf("\t Max: %f\n", t->stat[k].max);
        printf("\t Mean: %f\n", t->stat[k].mean);
    }
}

ip_mat * bitmap_to_ip_mat(Bitmap * img){
    unsigned int i=0,j=0;

    unsigned char R,G,B;

    unsigned int h = img->h;
    unsigned int w = img->w;

    ip_mat * out = ip_mat_create(h, w,3,0);

    for (i = 0; i < h; i++)              /* rows */
    {
        for (j = 0; j < w; j++)          /* columns */
        {
            bm_get_pixel(img, j,i,&R, &G, &B);
            set_val(out,i,j,0,(float) R);
            set_val(out,i,j,1,(float) G);
            set_val(out,i,j,2,(float) B);
        }
    }

    return out;
}

Bitmap * ip_mat_to_bitmap(ip_mat * t){

    Bitmap *b = bm_create(t->w,t->h);

    unsigned int i, j;
    for (i = 0; i < t->h; i++)              /* rows */
    {
        for (j = 0; j < t->w; j++)          /* columns */
        {
            bm_set_pixel(b, j,i, (unsigned char) get_val(t,i,j,0),
                    (unsigned char) get_val(t,i,j,1),
                    (unsigned char) get_val(t,i,j,2));
        }
    }
    return b;
}

float get_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k){
    if(i<a->h && j<a->w &&k<a->k){  /* j>=0 and k>=0 and i>=0 is non sense*/
        return a->data[i][j][k];
    }else{
        printf("Errore get_val!!!");
        exit(1);
    }
}

void set_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k, float v){
    if(i<a->h && j<a->w &&k<a->k){
        a->data[i][j][k]=v;
    }else{
        printf("Errore set_val!!!");
        exit(1);
    }
}

float get_normal_random(){
    float y1 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float y2 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    return cos(2*PI*y2)*sqrt(-2.*log(y1));
}

ip_mat * ip_mat_create(unsigned int h, unsigned int w,unsigned int k, float v) {
    int i, j, m;
    ip_mat * ip_mat_new = (ip_mat *) malloc (sizeof(ip_mat *));
    
    stats * stat = (stats *) malloc (sizeof(stats *));
    stat -> min = v;
    stat -> max = v;
    stat -> mean = v;
    
    float *** data = (float ***) malloc(h * sizeof(float **));
    for (i = 0; i < h; i++) {
        data[i] = (float **) malloc (w * sizeof(float *));
        for (j = 0; j < w; j++) {
            data[i][j] = (float *) malloc (k * sizeof(float));
        }
    }

    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j++) {
            for (m = 0; m < k; m++) {
                data[i][j][m] = v;
            }
        }
    }
    
    ip_mat_new -> w = w;
    ip_mat_new -> h = h;
    ip_mat_new -> k = k;
    ip_mat_new -> stat = stat;
    ip_mat_new -> data = data;
    return ip_mat_new;
}

void compute_stats(ip_mat * t) {
    float min = 255, max = 0, mean = 0, cont = 0;
    int i, j, m;
    
    for (i = 0; i < t -> h; i++) {
        for (j = 0; j < t -> w; j++) {
            for (m = 0; m < t -> k; m++) {
                float value = (t -> data)[i][j][m];
                if(value > max) {
                    max = value;
                }
                if(value < min) {
                    min = value;
                }
                mean += value;
                cont++;
            }
        }
    }
    
    t -> stat -> min = min;
    t -> stat -> max = max;
    t -> stat -> mean = mean / cont;
}

ip_mat * ip_mat_concat(ip_mat * a, ip_mat * b, int dimensione) {
    int i, j, m;
    ip_mat * ip_mat_new;
    if (dimensione == 0) {
        if (a -> w == b -> w && a -> k == b -> k) {
            ip_mat_new = ip_mat_create(a -> h + b -> h, a -> w, a -> k, 0);
            for (i = 0; i < a -> h + b -> h; i++) {
                for (j = 0; j < a -> w; j++) {
                    for (m = 0; m < a -> k; m++) {
                        if (i < a -> h) {
                            ip_mat_new -> data[i][j][m] = (a -> data)[i][j][m];
                        } else {
                            ip_mat_new -> data[i][j][m] = (b -> data)[i - a -> h][j][m];
                        }
                    }
                }
            }
        } else {
            ip_mat_new = NULL;
        }
    } else if (dimensione == 1) {
        if (a -> h == b -> h && a -> k == b -> k) {
            ip_mat_new = ip_mat_create(a -> h, a -> w + b -> w, a -> k, 0);
            for (i = 0; i < a -> h; i++) {
                for (j = 0; j < a -> w + b -> w; j++) {
                    for (m = 0; m < a -> k; m++) {
                        if (j < a -> w) {
                            ip_mat_new -> data[i][j][m] = (a -> data)[i][j][m];
                        } else {
                            ip_mat_new -> data[i][j][m] = (b -> data)[i][j - a -> w][m];
                        }
                    }
                }
            }
        } else {
            ip_mat_new = NULL;
        }
    } else {
        if (a -> h == b -> h && a -> w == b -> w) {
            ip_mat_new = ip_mat_create(a -> h, a -> w, a -> k + b -> k, 0);
            for (i = 0; i < a -> h; i++) {
                for (j = 0; j < a -> w; j++) {
                    for (m = 0; m < a -> k + b -> k; m++) {
                        if (m < a -> k) {
                            ip_mat_new -> data[i][j][m] = (a -> data)[i][j][m];
                        } else {
                            ip_mat_new -> data[i][j][m] = (b -> data)[i][j][m - a -> k];
                        }
                    }
                }
            }
        } else {
            ip_mat_new = NULL;
        }
    }
    
    return ip_mat_new;
}
