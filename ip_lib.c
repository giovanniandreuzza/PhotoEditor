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

/* Restituisce una sotto-matrice, ovvero la porzione individuata da:
 * t->data[row_start...row_end][col_start...col_end][0...k]
 * La terza dimensione la riportiamo per intero, stiamo in sostanza prendendo un sottoinsieme
 * delle righe e delle colonne.
 * */						
ip_mat * ip_mat_subset(ip_mat * t, unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end){
	/*(row_end-row_start)+1 sarà la dimensione di "righe", stessa cosa per le colonne */
	unsigned int i,j,l;
	ip_mat * sub = ip_mat_create((row_end-row_start+1), (col_end-col_start+1), t->k, 0.0);
	for(i = row_start; i <= (row_end-row_start); i++){
		for(j = col_start; j <= (col_end-col_start); j++){
			for(l = 0; l < t->k; l++){
				sub->data[i-row_start][j-col_start][l] = t->data[i][j][l];
			}
		}
	}
	return sub;
}

/* Moltiplica un ip_mat per uno scalare c. Si moltiplica c per tutti gli elementi di "a"
 * e si salva il risultato in un nuovo tensore in output. */
ip_mat * ip_mat_mul_scalar(ip_mat *a, float c){
	unsigned int i,l,j;
	ip_mat * tensore = ip_mat_create(a->h, a->w, a->k, 0.0);
	for(i = 0; i< a->h; i++){
		for(j = 0; j < a->w; j++){
			for(l = 0; l < a->k; l++){
				tensore->data[i][j][l] = (a->data[i][j][l])*c;
			}
		}
	}
	return tensore;
}

/* Effettua la fusione (combinazione convessa) di due immagini */ 
ip_mat * ip_mat_blend(ip_mat * a, ip_mat * b, float alpha){
	unsigned int i,j,l;
	ip_mat * fusion = ip_mat_create(a->h, a->w, a->k, 0);
	for(i = 0; i < a->h; i++){
		for(j = 0; j < a->w; j++){
			for(l = 0; l < a->k; l++){
				fusion->data[i][j][l] = (alpha*(a->data[i][j][l])+((1-alpha)*(b->data[i][j][l])));
			}
		}
	}
	return fusion;
}

/* Aggiunge un padding all'immagine. Il padding verticale è pad_h mentre quello
 * orizzontale è pad_w.
 * L'output sarà un'immagine di dimensioni:
 *      out.h = a.h + 2*pad_h;
 *      out.w = a.w + 2*pad_w;
 *      out.k = a.k
 * con valori nulli sui bordi corrispondenti al padding e l'immagine "a" riportata
 * nel centro
 * */ 
ip_mat * ip_mat_padding(ip_mat * a, int pad_h, int pad_w){
	unsigned int i,j,l;
	int supp_i = 0, supp_j = 0;
	ip_mat * out = ip_mat_create(a->h + 2*pad_h, a->w + 2*pad_w, a->k, 0);
	for(i = 0; i < out->h; i++){
		for(j = 0; j < out->w; j++){
			for(l = 0; l < out->k; l++){
				if(supp_i < pad_h || i > (out->h-pad_h) || supp_j < pad_w || j > (out->w-pad_w)){	/*setto la cella della matrice a zero solo se l'indice è minore della*/
					out->data[i][j][l] = 0;												/*DIM del padding o maggiore della dim della matrice originaria*/
				}else{
					out->data[i][j][l] = a->data[i-pad_h][j-pad_w][l];
				}
			}
			supp_j++;
		}
		supp_i++;
	}
	return out;
}

/* Crea un filtro medio per la rimozione del rumore */
ip_mat * create_average_filter(int w, int h, int k){
	int i,j,l;
	ip_mat * average = ip_mat_create(h, w, k,0);
	for(i = 0; i<h; i++){
		for(j = 0; j<w; j++){
			for(l = 0; l<k; l++){
				average->data[i][j][l] = (1/(h*w));
			}		
		}
	}
	return average;
}

/* Samuele*/


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
    ip_mat *out = NULL;

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

    out = ip_mat_create(3, 3, 3, 0.0);

    for(i = 0; i < 3; i++)
        for(j = 0; j < 3; j++)
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
    unsigned i, j, k, ii, jj;
    unsigned pad_h = ((f->h) - 1)/2;
    unsigned pad_w = ((f->h) - 1)/2;

    ip_mat *out;
    ip_mat *pad;

    out = ip_mat_create(a->h, a->w, a->k, 0.0);/**che sia da fare il controllo che *a != NULL?*/
    pad = ip_mat_padding(a, pad_h, pad_w)
    
    for(k = 0; k < out->k; k++)
        for(i = 0; i < out->h; i++)
            for(j = 0; j < out->w; j++)        
            {               
                ip_mat *supp;
                float val = 0.0;

                supp = ip_mat_subset(pad, i, i + f->h, j, j + f->w);

                for(ii = 0; ii < supp->h; ii++)
                    for(jj = 0; jj < supp->w; jj++)
                        val += get_val(a, i + ii, j + jj; k) * get_val(sub, ii, jj, k);
 
                set_val(out, i, j, k, val);
                ip_mat_free(supp);
            }

    ip_mat_free(pad);
    /*clamp(out, 0, 255); da verificare se va inserito*/
    compute_stats(out);

    return out;
}
