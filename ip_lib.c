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

float get_normal_random(float media, float std){

    float y1 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float y2 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float num = cos(2*PI*y2)*sqrt(-2.*log(y1));

    return media + num*std;
}

/* Libera la memoria (data, stat e la struttura) */
void ip_mat_free(ip_mat *a){/*da verificare in base alla create*/
    unsigned i, j;

    if(a != NULL){
        free(a->stat);

        for(i=0; i<a->h; i++){
            for(j=0; j<a->w; j++){
                free(a->data[i][j]);
            }
            free(a->data[i]);
        }
        free(a->data);
        free(a);
    }    
}

/* Inizializza una ip_mat con dimensioni w h e k.
 * Ogni elemento Ã¨ generato da una gaussiana con media mean e varianza var */
void ip_mat_init_random(ip_mat * t, float mean, float var){
    unsigned i, j, k;

    for(i=0; i<t->h; i++)
        for(j=0; j<t->w; j++)
            for(k=0; k<t->k; k++){
                set_val(t, i, j, k, get_normal_random(mean, var));
            }
}

/* Esegue la somma di due ip_mat (tutte le dimensioni devono essere identiche)
 * e la restituisce in output. */
ip_mat * ip_mat_sum(ip_mat * a, ip_mat * b){
    unsigned i, j, k;
    ip_mat *out; 

    if((a->h != b->h) && (a->w != b->w) && (a->k != b->k)){
        printf("The dimensions are different");
        exit(1);
    }
    else{
    	out = ip_mat_create(a->h, a->w, a->k, 0.0);
    	
        for(i=0; i<out->h; i++)
            for(j=0; j<out->w; j++)
                for(k=0; k<out->k; k++)
                    set_val(out, i, j, k, (get_val(a, i, j, k) + get_val(b, i, j, k)));

        compute_stats(out);
        
        return out;
    }
}

/* Calcola la media di due ip_mat a e b e la restituisce in output.*/
ip_mat * ip_mat_mean(ip_mat * a, ip_mat * b){
    unsigned i, j, k;
    ip_mat *out;

    out = ip_mat_create(a->h, a->w, a->k, 0.0);

    if((a->h != b->h) && (a->w != b->w) && (a->k != b->k)){
        printf("The dimensions are different");
        exit(1);
    }
    else{
        for(i=0; i<out->h; i++)
            for(j=0; j<out->w; j++)
                for(k=0; k<out->k; k++)
                    set_val(out, i, j, k, (get_val(a, i, j, k) + get_val(b, i, j, k)/2));

        compute_stats(out);
        return out;
    }

    compute_stats(out);
    return out;
}

/* Operazione di corruzione con rumore gaussiano:
 * Aggiunge del rumore gaussiano all'immagine, il rumore viene enfatizzato
 * per mezzo della variabile amount.
 * out = a + gauss_noise*amount
 * */
ip_mat * ip_mat_corrupt(ip_mat * a, float amount){
    unsigned i, j, k;
    ip_mat *out;

    out = ip_mat_create(a->h, a->w, a->k, 0.0);

    for(i=0; i<a->h; i++)
        for(j=0; j<a->w; j++)
        {
            for(k=0; k<a->k; k++)
                set_val(out, i, j, k, (get_val(a, i, k, k) + get_normal_random(0.0, (amount / 2)));
        }
            
    /*bisogna aggiungere controllo che i valori siano tra 0 e 255 ?*/

    compute_stats(out);
    return out;
}

/* Crea un filtro per rilevare i bordi */
ip_mat * create_edge_filter(){
    unsigned i, j, k;
    ip_mat *out;

    out = ip_mat_create(3, 3, 1, -1);

    set_val(out, 1, 1, 0, 8);

    compute_stats(out);
    return out;
}

/* Effettua una riscalatura dei dati tale che i valori siano in [0,new_max].
 * Utilizzate il metodo compute_stat per ricavarvi il min, max per ogni canale.
 *
 * I valori sono scalati tramite la formula valore-min/(max - min)
 *
 * Si considera ogni indice della terza dimensione indipendente, quindi l'operazione
 * di scalatura va ripetuta per ogni "fetta" della matrice 3D.
 * Successivamente moltiplichiamo per new_max gli elementi della matrice in modo da ottenere un range
 * di valori in [0,new_max].
 * */
void rescale(ip_mat * t, float new_max){
    unsigned i, j, k;

    compute_stats(t);
    float min = t->stat->min;
    float max = t->stat->max;

    for(i=0; i<t->h; i++)
        for(j=0; j<t->w; j++)
            for(k=0; k<t->k; k++)
            {
                float val = get_val(t, i, j, k);
                float scal = ((val - min) / (max - min)) * new_max;
                
                set_val(t, i, j, k, scal);
            }
}
