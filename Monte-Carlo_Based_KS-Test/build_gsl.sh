input=$1

prefix="${input%.c}"
gcc -Wall -I/usr/local/include -c $input -lpthread
gcc -L/usr/local/lib $prefix".o" -lgsl -lgslcblas -lm -lpthread -O3 -o $prefix
