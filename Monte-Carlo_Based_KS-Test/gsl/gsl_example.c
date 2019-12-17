#include <stdio.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

int main (void)
{
    gsl_rng* gsl_rng_bin = gsl_rng_alloc(gsl_rng_taus);
    for(int i=0; i<10; i++){
        printf("%d\n",gsl_ran_binomial(gsl_rng_bin, 0.5, 10));
    }
    gsl_rng_free(gsl_rng_bin);
    return 0;
}
