input=$1

prefix="${input%.c}"
gcc -Wall -I/usr/local/include -c $input
gcc -L/usr/local/lib $prefix".o" -lgsl -lgslcblas -lm -O3 -o $prefix
