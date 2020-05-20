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

    compute_stats(out);

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
    if(i<a->h && j<a->w &&k<a->k){
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

/* Crea una copia di una ip_mat e lo restituisce in output */
ip_mat * ip_mat_copy(ip_mat * in) /*che sia da controllare che *in != NULL ?*/
{
    unsigned i, j, k;
    ip_mat *out;

    out = ip_mat_create(a->h, a->w, a->k, 0.0);

    for(i = 0; i < in->h; i++)
        for(j = 0; j < in->w; j++)
            for(k = 0; k < in->k; k++)
                set_val(out, i, j, k, get_val(in, i, j, k));

    compute_stats(out);

    return out;
}

/* Esegue la sottrazione di due ip_mat (tutte le dimensioni devono essere identiche)
 * e la restituisce in output.
 * */
ip_mat * ip_mat_sub(ip_mat * a, ip_mat * b)
{
    unsigned i, j, k;
    ip_mat *out;

    if((a->h != b->h) && (a->w != b->w) && (a->k != b->k)) /*controllo che le dimensioni siano identiche*/
    {
        printf("The dimensions are different")
        exit(1); /*inserire codice errore*/
    }
    else
    {
        out = ip_mat_create(a->h, a->w, a->k, 0.0);

        for(i = 0; i < out->h; i++)
            for(j = 0; j < out->w; j++)
                for(k = 0; k < out->k; k++)
                    set_val(out, i, j, k, (get_val(a, i, k, k) - get_val(b, i, k, k)));

        compute_stats(out);

        return out;
    }
}

/* Converte un'immagine RGB ad una immagine a scala di grigio.
 * Quest'operazione viene fatta calcolando la media per ogni pixel sui 3 canali
 * e creando una nuova immagine avente per valore di un pixel su ogni canale la media appena calcolata.
 * Avremo quindi che tutti i canali saranno uguali.
 * */
/*la nuova ip_mat ha 3 canali con valori tutti uguali*/
ip_mat * ip_mat_to_gray_scale(ip_mat * in)
{
    unsigned i, j, k, n;
    float sum, average;
    ip_mat *out;

    out = ip_mat_create(a->h, a->w, a->k, 0.0); /*che sia da fare il controllo che *in non sia NULL?*/

    for(i = 0; i < in->h; i++)
        for(j = 0; j < in->w; j++)
        {
            /*effetuo la somma dei tre livelli*/
            sum = 0.0;
            for(k = 0; k < in->k; k++)
                sum += get_val(in, i, j, k);

            for(n = 0; n < 3; n++)
                set_val(out, i, j, n, (sum / 3));
        }
        
    compute_stats(out);
    
    return out;
}

/* Crea un filtro per aggiungere profondità  */
ip_mat * create_emboss_filter()
{
    unsigned i, j, k;
    ip_mat *out;
    float emboss_val[][] = {
        {-2.0, -1.0, 0.0}, 
        {-1.0, 1.0, 1.0}, 
        {0.0, 1.0, 2.0}
    };

    out = ip_mat_create(3, 3, 3, 0.0); /*ip_mat 3X3X3 inizializzata a 0.0*/

    for(i = 0; i < 3; i++)
        for(j = 0; i < 3; j++)
            for(k = 0; k < 3; k++)
                set_val(out, i, j, k, emboss_val[i][j]);

    return out;
}

/* Nell'operazione di clamping i valori <low si convertono in low e i valori >high in high.*/
void clamp(ip_mat * t, float low, float high)
{
    unsigned i, j, k;

    for(i = 0; i < t->w; i++)
        for(j = 0; j < t->h; j++)
            for(k = 0; k < t->k; k++)
            {
                if(get_val(t, i, j, k) < low)
                    set_val(t, i, j, k, low);
                if(get_val(t, i, j, k) > high)
                    set_val(t, i, j, k, high);
            }
}

/* Effettua la convoluzione di un ip_mat "a" con un ip_mat "f".
 * La funzione restituisce un ip_mat delle stesse dimensioni di "a".
 * */
ip_mat * ip_mat_convolve(ip_mat * a, ip_mat * f)
{
    unsigned i, j, k;
    ip_mat *out;

    out = ip_mat_create(a->h, a->w, a->k, 0.0);


}