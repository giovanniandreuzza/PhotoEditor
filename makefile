main_iplib: main_iplib.o bmp.o ip_lib.o
	gcc main_iplib.o bmp.o ip_lib.o -omain_iplib -lm -Wall --ansi --pedantic -g3 -O3 -std=gnu89 -Wextra -fsanitize=address -fsanitize=undefined

main_iplib.o: main_iplib.c ip_lib.h bmp.h
	gcc main_iplib.c -omain_iplib.o -c -lm -Wall --ansi --pedantic -g3 -O3 -std=gnu89 -Wextra -fsanitize=address -fsanitize=undefined

ip_lib.o: ip_lib.c bmp.h
	gcc ip_lib.c -oip_lib.o -c -lm -Wall --ansi --pedantic -g3 -O3 -std=gnu89 -Wextra -fsanitize=address -fsanitize=undefined
