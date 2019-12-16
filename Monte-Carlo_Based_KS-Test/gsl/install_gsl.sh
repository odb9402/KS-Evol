wget ftp://ftp.gnu.org/gnu/gsl/gsl-2.4.tar.gz
tar -zxvf gsl-2.4.tar.gz

cd gsl-2.4

mkdir gsl
./configure
make
make check
make install

