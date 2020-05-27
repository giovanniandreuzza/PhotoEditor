/*
 Created by Sebastiano Vascon on 23/03/20.
*/

#include <stdio.h>
#include <math.h>
#include "ip_lib.h"

#include "bmp.h"

void ip_mat_show(ip_mat * t) {
    unsigned int i, l, j;
    printf("Matrix of size %d x %d x %d (hxwxk)\n", t -> h, t -> w, t -> k);
    for (l = 0; l < t -> k; l++) {
        printf("Slice %d\n", l);
        for (i = 0; i < t -> h; i++) {
            for (j = 0; j < t -> w; j++) {
                printf("%f ", get_val(t, i, j, l));
            }
            printf("\n");
        }
        printf("\n");
    }
}

void ip_mat_show_stats(ip_mat * t) {
    unsigned int k;

    compute_stats(t);

    for (k = 0; k < t -> k; k++) {
        printf("Channel %d:\n", k);
        printf("\t Min: %f\n", t -> stat[k].min);
        printf("\t Max: %f\n", t -> stat[k].max);
        printf("\t Mean: %f\n", t -> stat[k].mean);
    }
}

ip_mat * bitmap_to_ip_mat(Bitmap * img) {
    unsigned int i = 0, j = 0;

    unsigned char R, G, B;

    unsigned int h = img -> h;
    unsigned int w = img -> w;

    ip_mat * out = ip_mat_create(h, w, 3, 0);

    for (i = 0; i < h; i++) /* rows */ {
        for (j = 0; j < w; j++) /* columns */ {
            bm_get_pixel(img, j, i, & R, & G, & B);
            set_val(out, i, j, 0, (float) R);
            set_val(out, i, j, 1, (float) G);
            set_val(out, i, j, 2, (float) B);
        }
    }

    return out;
}

Bitmap * ip_mat_to_bitmap(ip_mat * t) {

    Bitmap * b = bm_create(t -> w, t -> h);

    unsigned int i, j;
    for (i = 0; i < t -> h; i++) /* rows */ {
        for (j = 0; j < t -> w; j++) /* columns */ {
            bm_set_pixel(b, j, i, (unsigned char) get_val(t, i, j, 0),
                (unsigned char) get_val(t, i, j, 1),
                (unsigned char) get_val(t, i, j, 2));
        }
    }
    return b;
}

float get_val(ip_mat * a, unsigned int i, unsigned int j, unsigned int k) {
    if (i < a -> h && j < a -> w && k < a -> k) {
        /* j>=0 and k>=0 and i>=0 is non sense*/
        return a -> data[i][j][k];
    } else {
        printf("Errore get_val!!!");
        exit(1);
    }
}

void set_val(ip_mat * a, unsigned int i, unsigned int j, unsigned int k, float v) {
    if (i < a -> h && j < a -> w && k < a -> k) {
        a -> data[i][j][k] = v;
    } else {
        printf("Errore set_val!!!");
        exit(1);
    }
}

float get_normal_random(float media, float std) {
    float y1 = ((float)(rand()) + 1.) / ((float)(RAND_MAX) + 1.);
    float y2 = ((float)(rand()) + 1.) / ((float)(RAND_MAX) + 1.);
    float num = cosf(2 * PI * y2) * sqrt(-2. * log(y1));

    return media + num * std;
}

ip_mat * ip_mat_create(unsigned int h, unsigned int w, unsigned int k, float v) {
    unsigned int i, j, m;
    float *** data = (float ***) malloc(h * sizeof(float ** ));
    
    ip_mat * ip_mat_new = (ip_mat *) malloc(sizeof(ip_mat));

    stats * stat = (stats *) malloc(sizeof(stats));
    stat -> min = v;
    stat -> max = v;
    stat -> mean = v;
  
    for (i = 0; i < h; i++) {
        data[i] = (float ** ) malloc(w * sizeof(float * ));
        for (j = 0; j < w; j++) {
            data[i][j] = (float * ) malloc(k * sizeof(float));
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

void ip_mat_free(ip_mat * a) {
    unsigned int i, j;

    if (a != NULL) {
        free(a -> stat);

        for (i = 0; i < a -> h; i++) {
            for (j = 0; j < a -> w; j++) {
                free(a -> data[i][j]);
            }
            free(a -> data[i]);
        }
        
        free(a -> data);
        free(a);
    }
}

void compute_stats(ip_mat * t) {
    float min = 255, max = 0, mean = 0, cont = 0;
    unsigned int i, j, m;

    for (i = 0; i < t -> h; i++) {
        for (j = 0; j < t -> w; j++) {
            for (m = 0; m < t -> k; m++) {
                float value = (t -> data)[i][j][m];
                if (value > max) {
                    max = value;
                }
                if (value < min) {
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

void ip_mat_init_random(ip_mat * t, float mean, float var) {
    unsigned i, j, k;

    for (i = 0; i < t -> h; i++) {
        for (j = 0; j < t -> w; j++) {
            for (k = 0; k < t -> k; k++) {
                set_val(t, i, j, k, get_normal_random(mean, var));
            }
        }
    }
}

ip_mat * ip_mat_copy(ip_mat * in ) {

    unsigned i, j, k;
    ip_mat * out;

    out = ip_mat_create(in -> h, in -> w, in -> k, 0.0);

    for (i = 0; i < in -> h; i++) {
        for (j = 0; j < in -> w; j++) {
            for (k = 0; k < in -> k; k++) {
                set_val(out, i, j, k, get_val( in , i, j, k));
	    	}
		}
    }

    compute_stats(out);

    return out;
}

ip_mat * ip_mat_subset(ip_mat * t, unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end) {
    /*(row_end-row_start)+1 sarà la dimensione di "righe", stessa cosa per le colonne */
    unsigned int i, j, l;
    
    ip_mat * sub = ip_mat_create((row_end - row_start + 1), (col_end - col_start + 1), t -> k, 0.0);
    
    for (i = row_start; i <= (row_end - row_start); i++) {
        for (j = col_start; j <= (col_end - col_start); j++) {
            for (l = 0; l < t -> k; l++) {
                set_val(sub, i - row_start, j - col_start, l, (get_val(t, i, j, l)));
            }
        }
    }
    
    return sub;
}

ip_mat * ip_mat_concat(ip_mat * a, ip_mat * b, int dimensione) {
    unsigned int i, j, m;
    ip_mat * ip_mat_new;
    if (dimensione == 0) {
        if (a -> w == b -> w && a -> k == b -> k) {
            ip_mat_new = ip_mat_create(a -> h + b -> h, a -> w, a -> k, 0);
            for (i = 0; i < a -> h + b -> h; i++) {
                for (j = 0; j < a -> w; j++) {
                    for (m = 0; m < a -> k; m++) {
                        if (i < a -> h) {
                            /*ip_mat_new -> data[i][j][m] = (a -> data)[i][j][m];*/
                            set_val(ip_mat_new, i, j, m, (a -> data)[i][j][m]);
                        } else {
                            /*ip_mat_new -> data[i][j][m] = (b -> data)[i - a -> h][j][m];*/
                            set_val(ip_mat_new, i, j, m, (b -> data)[i - a -> h][j][m]);
                        }
                    }
                }
            }
        } else {
            printf("Errore...\n");
            exit(1);
        }
    } else if (dimensione == 1) {
        if (a -> h == b -> h && a -> k == b -> k) {
            ip_mat_new = ip_mat_create(a -> h, a -> w + b -> w, a -> k, 0);
            for (i = 0; i < a -> h; i++) {
                for (j = 0; j < a -> w + b -> w; j++) {
                    for (m = 0; m < a -> k; m++) {
                        if (j < a -> w) {
                            /*ip_mat_new -> data[i][j][m] = (a -> data)[i][j][m];*/
                            set_val(ip_mat_new, i, j, m, (a -> data)[i][j][m]);
                        } else {
                            /*ip_mat_new -> data[i][j][m] = (b -> data)[i][j - a -> w][m];*/
                            set_val(ip_mat_new, i, j, m, (b -> data)[i][j - a -> w][m]);
                        }
                    }
                }
            }
        } else {
            printf("Errore...\n");
            exit(1);
        }
    } else {
        if (a -> h == b -> h && a -> w == b -> w) {
            ip_mat_new = ip_mat_create(a -> h, a -> w, a -> k + b -> k, 0);
            for (i = 0; i < a -> h; i++) {
                for (j = 0; j < a -> w; j++) {
                    for (m = 0; m < a -> k + b -> k; m++) {
                        if (m < a -> k) {
                            /*ip_mat_new -> data[i][j][m] = (a -> data)[i][j][m];*/
                            set_val(ip_mat_new, i, j, m, (a -> data)[i][j][m]);
                        } else {
                            /*ip_mat_new -> data[i][j][m] = (b -> data)[i][j][m - a -> k];*/
                            set_val(ip_mat_new, i, j, m, (b -> data)[i][j][m - a -> k]);
                        }
                    }
                }
            }
        } else {
            printf("Errore...\n");
            exit(1);
        }
    }

    return ip_mat_new;
}

ip_mat * ip_mat_sum(ip_mat * a, ip_mat * b) {
    unsigned i, j, k;
    ip_mat * out;

    if ((a -> h != b -> h) && (a -> w != b -> w) && (a -> k != b -> k)) {
        printf("Immagini di dimensioni diverse\n");
        exit(1);
    } else {
        out = ip_mat_create(a -> h, a -> w, a -> k, 0.0);

        for (i = 0; i < out -> h; i++)
            for (j = 0; j < out -> w; j++)
                for (k = 0; k < out -> k; k++)
                    set_val(out, i, j, k, (get_val(a, i, j, k) + get_val(b, i, j, k)));

        compute_stats(out);

        return out;
    }
}

ip_mat * ip_mat_sub(ip_mat * a, ip_mat * b) {

    unsigned i, j, k;
    ip_mat * out = NULL;

    if ((a -> h != b -> h) && (a -> w != b -> w) && (a -> k != b -> k)) /*controllo che le dimensioni siano identiche*/ {
        printf("Immagini di dimensioni diverse\n");
        exit(1); /*inserire codice errore*/
    } else {
        out = ip_mat_create(a -> h, a -> w, a -> k, 0.0);

        for (i = 0; i < out -> h; i++) {
            for (j = 0; j < out -> w; j++) {
                for (k = 0; k < out -> k; k++) {
                    set_val(out, i, j, k, (get_val(a, i, j, k) - get_val(b, i, j, k)));
                }
            }
        }

        compute_stats(out);

        return out;
    }
}

ip_mat * ip_mat_mul_scalar(ip_mat * a, float c) {
    unsigned int i, l, j;
    
    ip_mat * tensore = ip_mat_create(a -> h, a -> w, a -> k, 0.0);
    
    for (i = 0; i < a -> h; i++) {
        for (j = 0; j < a -> w; j++) {
            for (l = 0; l < a -> k; l++) {
                set_val(tensore, i, j, l, (get_val(a, i, j, l))*c);
            }
        }
    }
    
    return tensore;
}

ip_mat * ip_mat_add_scalar(ip_mat * a, float c) {
    unsigned int i, j, m;
    ip_mat * ip_mat_new = ip_mat_create(a -> h, a -> w, a -> k, 0.0);
    for (i = 0; i < a -> h; i++) {
        for (j = 0; j < a -> w; j++) {
            for (m = 0; m < a -> k; m++) {
                set_val(ip_mat_new, i, j, m, (a -> data)[i][j][m] + c);
            }
        }
    }
    clamp(ip_mat_new, 0.0, 255.0);
    return ip_mat_new;
}

ip_mat * ip_mat_mean(ip_mat * a, ip_mat * b) {
    unsigned i, j, k;
    ip_mat * out;

    out = ip_mat_create(a -> h, a -> w, a -> k, 0.0);

    if ((a -> h != b -> h) && (a -> w != b -> w) && (a -> k != b -> k)) {
        printf("Immagini di dimensioni diverse\n");
        exit(1);
    } else {
        for (i = 0; i < out -> h; i++) {
            for (j = 0; j < out -> w; j++) {
                for (k = 0; k < out -> k; k++) {
                    set_val(out, i, j, k, (get_val(a, i, j, k) + get_val(b, i, j, k) / 2));
				}
			}
		}
        compute_stats(out);
        return out;
    }

    compute_stats(out);
    return out;
}

ip_mat * ip_mat_to_gray_scale(ip_mat * in ) {
    unsigned i, j, k, n;
    float sum;
    ip_mat * out;

    out = ip_mat_create(in -> h, in -> w, in -> k, 0.0); 
    
    for (i = 0; i < in -> h; i++) {
        for (j = 0; j < in -> w; j++) {
            sum = 0.0;
            for (k = 0; k < in -> k; k++) {
                sum += get_val( in , i, j, k);
            }

            for (n = 0; n < 3; n++) {
                set_val(out, i, j, n, (sum / 3));
            }
        }
    }

    compute_stats(out);

    return out;
}

ip_mat * ip_mat_blend(ip_mat * a, ip_mat * b, float alpha) {
	if ((a -> h == b -> h) && (a -> w == b -> w) && (a -> k == b -> k)) {
        unsigned int i, j, l;
        
        ip_mat * fusion = ip_mat_create(a -> h, a -> w, a -> k, 0.0);
        
        for (i = 0; i < a -> h; i++) {
            for (j = 0; j < a -> w; j++) {
                for (l = 0; l < a -> k; l++) {
                    set_val(fusion, i, j, l, ((alpha * get_val(a, i, j, l)) + ((1 - alpha) * get_val(b, i, j, l))));
                }
            }
        }
        return fusion;
    } else {
		printf("Immagini di dimensioni diverse\n");
		exit(1);
	}
}

ip_mat * ip_mat_brighten(ip_mat * a, float bright) {
	if(a){
    	return ip_mat_add_scalar(a, bright);
	}
	else{
		printf("Errore: puntatore == NULL\n");
		exit(1);
	}
}

ip_mat * ip_mat_corrupt(ip_mat * a, float amount) {
    unsigned i, j, k;
    ip_mat * out;

    out = ip_mat_create(a -> h, a -> w, a -> k, 0.0);

    for (i = 0; i < a -> h; i++) {
        for (j = 0; j < a -> w; j++) {
            for (k = 0; k < a -> k; k++) {
                set_val(out, i, j, k, (get_val(a, i, j, k) + get_normal_random(0.0, amount)));
            }            
        }
    }
    
    compute_stats(out);
    return out;
}

ip_mat * ip_mat_convolve(ip_mat * a, ip_mat * f) {
    unsigned i, j, k, ii, jj;
    unsigned pad_h = ((f -> h) - 1) / 2;
    unsigned pad_w = ((f -> h) - 1) / 2;

    ip_mat * out;
    ip_mat * pad;

    out = ip_mat_create(a -> h, a -> w, a -> k, 0.0);
    pad = ip_mat_padding(a, pad_h, pad_w);

    for (i = 0; i < ((pad -> h) - (f -> h) + 1); i++) {
        for (j = 0; j < ((pad -> w) - (f -> w) + 1); j++) {
            for (k = 0; k < pad -> k; k++) {
            
            	float val = 0.0;

                for (ii = 0; ii < (f -> h) ; ii++) {
                    for (jj = 0; jj < (f -> w); jj++) {
                    
                        val = val+(pad ->data[i+ii][j+jj][k] * f ->data[ii][jj][k]);
                    }
                }
					
                set_val(out, i, j, k, val);
            }
        }
    }

    ip_mat_free(pad);
    compute_stats(out);

    return out;
}

ip_mat * ip_mat_padding(ip_mat * a, unsigned int pad_h, unsigned int pad_w) {

    unsigned int i, j, l;
    unsigned int supp_i = 0, supp_j = 0;
    
    ip_mat * out = ip_mat_create(a -> h + 2 * pad_h, a -> w + 2 * pad_w, a -> k, 0.0);
    
    for (i = 0; i < out -> h; i++) {
    	supp_i = i;
        for (j = 0; j < out -> w; j++) {
        	supp_j = j;
            for (l = 0; l < out -> k; l++) {
                if (supp_i < pad_h || i >= (a -> h + pad_h) || supp_j < pad_w || j >= (a -> w + pad_w)) {
                    set_val(out, i, j, l, 0.0); 
                } else {
                    set_val(out, i, j, l, get_val(a, (i - pad_h), (j - pad_w), l));
                }
            }
        }
    }
    return out;
}

ip_mat * create_sharpen_filter() {

    unsigned i, j, k;
    ip_mat * out;
    
    float sharpen_val[][3] = {
        {0.0, -1.0, 0.0},
        {-1.0, 5.0, -1.0},
        {0.0, -1.0, 0.0}
    };

    out = ip_mat_create(3, 3, 3, 0.0);

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            for (k = 0; k < 3; k++) {
                set_val(out, i, j, k, sharpen_val[i][j]);
	    	}
		}
    }

    return out;
}

ip_mat * create_edge_filter() {
	int i;
    ip_mat * out;

    out = ip_mat_create(3, 3, 3, -1);
	
	for (i = 0; i < 3; i++){
		set_val(out, 1, 1, i, 8);
	}

    compute_stats(out);
    return out;
}

ip_mat * create_emboss_filter() {
    unsigned i, j, k;
    ip_mat * out;
    float emboss_val[][3] = {
        {-2.0, -1.0, 0.0},
        {-1.0, 1.0, 1.0},
        {0.0, 1.0, 2.0}
    };

    out = ip_mat_create(3, 3, 3, 0.0);

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            for (k = 0; k < 3; k++) {
                set_val(out, i, j, k, emboss_val[i][j]);
            }
        }
    }

    return out;
}

ip_mat * create_average_filter(unsigned int w, unsigned int h, unsigned int k) {
    unsigned int i, j, l;
    float avg = (1.0/(h*w));
    
    ip_mat * average = ip_mat_create(h, w, k, 0.0);
    
    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j++) {
            for (l = 0; l < k; l++) {
                set_val(average, i, j, l, avg);
                
            }
        }
    }
    return average;
}

ip_mat * create_gaussian_filter(unsigned int h, unsigned int w, unsigned int k, float sigma) {
    unsigned int i, j, m, cx = (int) w / 2, cy = (int) h / 2;
    float somma = 0.0;
    
    ip_mat * ip_mat_new = ip_mat_create(h, w, k, 0.0);
    
    for (m = 0; m < k; m++) {
        somma = 0.0;
        for (j = 0; j < w; j++) {
            for (i = 0; i < h; i++) {
                int x = j - cx;
                int y = i - cy;
                set_val(ip_mat_new, i, j, m, (1.0 / (2 * PI * (sigma * sigma))) * (exp(-((x * x) + (y * y)) / (2 * sigma * sigma))));
                somma = somma +  ip_mat_new->data[i][j][m];
            }
        }
    }
    
    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j++) {
            for (m = 0; m < k; m++) {
                set_val(ip_mat_new, i, j, m, (get_val(ip_mat_new, i, j, m)) / somma);
            }
        }
    }
    compute_stats(ip_mat_new);
    return ip_mat_new;
}

void rescale(ip_mat * t, float new_max) {
    unsigned i, j, k;
    
    float min = t -> stat -> min;
    float max = t -> stat -> max;
    
    compute_stats(t);
    
    for (i = 0; i < t -> h; i++) {
        for (j = 0; j < t -> w; j++) {
            for (k = 0; k < t -> k; k++) {
                float val = get_val(t, i, j, k);
                float scal = ((val - min) / (max - min)) * new_max;

                set_val(t, i, j, k, scal);
            }
        }
    }
}

void clamp(ip_mat * t, float low, float high) {
    unsigned i, j, k;
    for (i = 0; i < t -> h; i++) {
        for (j = 0; j < t -> w; j++) {
            for (k = 0; k < t -> k; k++) {
                if (get_val(t, i, j, k) < low) {
                    set_val(t, i, j, k, low);
                }
                if (get_val(t, i, j, k) > high) {
                    set_val(t, i, j, k, high);
                }
            }
        }
    }
}
