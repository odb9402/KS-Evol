/*

*/
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <stddef.h>
#include <pthread.h>
#include <sys/sysinfo.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>

#define HIST_BIN_SIZE 500000
#define FILE_LENGTH_LIMIT 100
#define THREAD_NUM 28

const char *argp_program_version = "binSampling 1.0";
const char *argp_program_bug_address = "<dhehdqls@gmail.com>";
gsl_rng* gsl_rng_bin; //GSL global random number generator
gsl_histogram* gsl_hist[THREAD_NUM]; //GSL global histogram allocator 

struct sampling_args{
    unsigned int average_allele;
    double p_0;
    int time_points;
    long sample_num;
    int thread_num;
};

void print_usage(){
    printf("\t<USAGE of binSampling>\n");
    printf("\tbinSampling will generate p-values based on bootstrapping that samples follow a binomial distribution\n");
    printf("\tbased on a null hypothesis of the time-dependent allele frequency data.\n\n");
    printf("\tbinSampling -n <number of samplings> -t <timepoints> -i <input> -o <output> \n\n");
}

int check_file_len(char* filename){
    /*
    * Check the file length of allele depth file.
    * The number of allele will be divided with its number of timepoints.
    *
    * return: The length of allele depth file.
    */
    FILE *fp;
    fp = fopen(filename, "r");
    if (fp==NULL){
        printf("Could not open file %s. \n", filename);
        return -1;
    }
    int count = 0;
    char c;
    for (c=getc(fp); c!=EOF; c=getc(fp))
        if (c=='\n')
            count++;
    fclose(fp);
    return count;
}

double get_p_0(int* alts, int* refs, int time_points){
    /*
    * Get a null hypothesis from given allele depths for all timepoints.
    * The null hypothesis represents the probability of occuring minor allele(SNP).
    *  
    * return: The null hypothesis of SNP.
    */
    int alt_sum = 0;
    int ref_sum = 0;
    
    for(int i = 0; i < time_points; i++){
        alt_sum += alts[i];
        ref_sum += refs[i];
    }
    
    return alt_sum/(double)(ref_sum + alt_sum);
}

double get_pval(int sampling_num, double ks){
    /*
    * Call an empirical p-value of ks-statistics from pre-constructed
    * gsl-histogram consisting of sampling values.
    *
    * Global variable : gsl_hist[0] (summarized gsl_histogram)
    * return : An empirical p-value of ks distance 
    */
    size_t hist_idx = -1;
    long empirical_cdf_count = 0;
    double p_val = -1.0;
    
    if(gsl_histogram_find(gsl_hist[0], ks, &hist_idx) == GSL_EDOM){
        printf("GSL histogram find error. \n");
        return -1;
    }
    //printf("max ks : %lf\n", ks);
    //printf("ks_hist_idx :  %ld\n", hist_idx);
    
    //for(int i=0; i < HIST_BIN_SIZE; i++)
    //    printf("%lf \t", gsl_histogram_get(gsl_hist,i));
    
    for(int i=0; i <= hist_idx; i++)
        empirical_cdf_count += gsl_histogram_get(gsl_hist[0], i);
    
    p_val = 1.0 - empirical_cdf_count/(double)sampling_num;
    //printf("pval : %.13f\n\n",p_val);
    return p_val;
}

void gsl_ks_distance_sampling(void *args){
    /*
    * Build the sampling histogram to calculate p-value for each SNP.
    *
    * Global variable : gsl_hist (gsl_histogram)
    *                   gsl_rng_bin (gsl random number generator)
    * Macro : THREAD_NUM ( the number of thread )
    * return : None
    */
    struct sampling_args *received_args = (struct sampling_args*)args;
    
    double max_ks_bin;
    unsigned int average_allele = received_args->average_allele;
    double p_0 = received_args->p_0;
    int time_points = received_args->time_points;
    long sample_num = received_args->sample_num;
    int thread_num = received_args->thread_num;
    double* sampled_p;
    clock_t start, end;
    
    sampled_p = malloc(sizeof(double)*time_points);
    start = clock();
    for(int i=0; i<sample_num/THREAD_NUM; i++){
        max_ks_bin = 0;
        
        for(int j=0; j<time_points; j++)
           sampled_p[j] = gsl_ran_binomial(gsl_rng_bin, p_0, average_allele)/(double)(average_allele);
        
        for(int j=0;j<time_points;j++){
            sampled_p[j] = fabs(sampled_p[j] - p_0); // No extra allocations.
            if(sampled_p[j] > max_ks_bin)
                max_ks_bin = sampled_p[j];
        }
        if(gsl_histogram_increment(gsl_hist[thread_num], max_ks_bin) == GSL_EDOM){
            printf("GSL histogram range error. \n");
            exit(-1);
        }
        if(i == 100000){
            end = clock();
            printf("\tETA : %fsec\n", ((sample_num/THREAD_NUM)/100000)*(double)(end-start)/CLOCKS_PER_SEC);
            start = clock();
        }
    }
}

int 
main(int argc, char* argv[]){
    int time_points = 0;
    unsigned long sample_num = 0;
    char* input_name;// = malloc(sizeof(char) * FILE_LENGTH_LIMIT);
    char* output_name;// = malloc(sizeof(char) * FILE_LENGTH_LIMIT);
    FILE *input_fp;
    FILE *output_fp;
    int allele_len = -1;
    
    int opt;
    while((opt=getopt(argc, argv, "ht:n:i:o:")) != -1 ){
        switch(opt){
            case 'h':
                print_usage();
            case 't':
                time_points = atoi(optarg);
                break;
            case 'n':
                sample_num = atol(optarg);
                break;
            case 'i':
                input_name = malloc(sizeof(char)*sizeof(optarg));
                if(optarg==NULL){
                    printf("-i need a second argument as the input file name.\n");
                    print_usage();
                    return -1;
                }
                else
                    strcpy(input_name, optarg);
                break;
            case 'o':
                output_name = malloc(sizeof(char)*sizeof(optarg));
                if(optarg==NULL){
                    printf("-o need a second argument as the output file name.\n");
                    print_usage();
                    return -1;
                }
                else
                    strcpy(output_name, optarg);
                break;
            case '?':
                if(optopt == 'i')
                    printf("option -i requires an input file name.\n");
                else if(optopt == 'o')
                    printf("option -o requires an output file name.\n");
                else if(optopt == 't')
                    printf("option -t requires the timepoints.\n");
                else if(optopt == 'n')
                    printf("option -n requires the number of sampling size.\n");
                else
                    printf("Unknown flag : %c \n", optopt);
                break;
        }
    }
    
    if(input_name == NULL){
        printf("option -i requires an input file name.\n");
        print_usage();
        return -1;
    }
    if(output_name == NULL){
        printf("option -o requires an output file name.\n");
        print_usage();
        return -1;
    }
    
    printf("Time points : %d\n", time_points);
    printf("The input file : %s\n", input_name);
    printf("The output file : %s\n", output_name);
    printf("Sample num : %lu\n", sample_num);
    
    printf("Count the total number of alleles . . . . \n\n");
    allele_len = check_file_len(input_name)/time_points;
    printf("Number of alleles of %s : %d \n\n",input_name,allele_len);
    
    input_fp = fopen(input_name, "r");
    if (input_fp==NULL){
        printf("Could not open file %s. \n", input_name);
        return 0;
    }
    output_fp = fopen(output_name, "w");
    if (output_fp==NULL){
        printf("Could not open file %s. \n", output_name);
        return 0;
    }
    
    ///////////////////// The actual process begin. ///////////////////////
    int* refs = malloc(sizeof(int)*time_points);
    int* alts = malloc(sizeof(int)*time_points);
    unsigned int ref_sum_total = 0;
    unsigned int alt_sum_total = 0;
    unsigned int ref_sum = 0;
    unsigned int alt_sum = 0;
    char* spliter;
    char* lines = NULL;
    int col_idx =0;
    double p_0 = 0;
    double* p_mt = malloc(sizeof(double)*time_points);
    double* ks = malloc(sizeof(double)*time_points);
    double max_ks = 0;
    double* ks_pvals = malloc(sizeof(double)*allele_len);
    size_t len = 0;
    
    const gsl_rng_type *T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    gsl_rng_bin = gsl_rng_alloc(T);
    for(int i=0; i<THREAD_NUM; i++){
        gsl_hist[i] = gsl_histogram_alloc(HIST_BIN_SIZE);
        gsl_histogram_set_ranges_uniform(gsl_hist[i], 0.0, 1.0);
    }
    
    pthread_t sampling_thread[THREAD_NUM];
    
    
    printf("Get p0 for binomial sampling. . . .\n");
    for(int allele_num=0; allele_num<allele_len; allele_num++){
        //File read and extract counts
        for(int i=0; i<time_points; i++){
            if(getline(&lines, &len, input_fp) == -1){
                printf("File read error. %d`th allele\n", allele_num);
                fclose(input_fp);
                fclose(output_fp);
                return -1;
            }
            spliter = strtok(lines,"\t");
            while(spliter != NULL){
                col_idx++;
                spliter = strtok(NULL,"\t");
                if(col_idx==5){
                    ref_sum += atoi(spliter);
                    ref_sum_total += atoi(spliter);
                }
                if(col_idx==6){
                    alt_sum += atoi(spliter);
                    alt_sum_total += atoi(spliter);
                }
            }
            col_idx = 0;
        }
        ref_sum = 0;
        alt_sum = 0;
    }
    rewind(input_fp);
    
    p_0 = alt_sum_total / (double)(alt_sum_total + ref_sum_total);
    printf("p0 for entire alleles : %lf\n", p_0);
    printf("The number of entire alleles : %d\n", ref_sum_total+alt_sum_total);
    printf("The total number of average alleles for each timepoint : %d\n",
            (alt_sum_total + ref_sum_total) / allele_len);
    printf("Binomial sampling based on p0 . . . . \n");
    
    printf("Build the sampling distribution of KS distances . . .");
    //Binomial sample values to ks distances from p_0. 
    struct sampling_args thread_args[THREAD_NUM];
    for(int i=0; i<THREAD_NUM; i++){
        //gsl_ks_distance_sampling((alt_sum_total + ref_sum_total) / allele_len, p_0, time_points, sample_num);
        thread_args[i].average_allele = (alt_sum_total + ref_sum_total) / allele_len;
        thread_args[i].p_0 = p_0;
        thread_args[i].time_points = time_points;
        thread_args[i].sample_num = sample_num;
        thread_args[i].thread_num = i;        
        if(pthread_create(&sampling_thread[i], NULL, gsl_ks_distance_sampling, (void *)&thread_args[i]) != 0){
            printf("Pthread create error. \n");
            return -1;
        }
    }
    for(int i=0; i<THREAD_NUM; i++)
        pthread_join(sampling_thread[i],NULL);
        
    for(int i=1; i<THREAD_NUM; i++){
        gsl_histogram_add(gsl_hist[0], gsl_hist[i]);
        printf("Sampling histogram merged . . . [%lu,%lu]\n",
                  (long)gsl_histogram_sum(gsl_hist[0]),sample_num);
    }
    if((sample_num - gsl_histogram_sum(gsl_hist[0]) > 0)){
        for(int i=1; i < (sample_num - gsl_histogram_sum(gsl_hist[0])); i++)
            continue;
    }
    
    printf("Calculate p-values based on the sampling distribution . . .\n");
    for(int allele_num=0; allele_num<allele_len; allele_num++){
        //File read and extract counts
        for(int i=0; i<time_points; i++){
            if(getline(&lines, &len, input_fp) == -1){
                printf("File read error. %d`th allele\n", allele_num);
                fclose(input_fp);
                fclose(output_fp);
                return -1;
            }
            spliter = strtok(lines,"\t");
            while(spliter != NULL){
                col_idx++;
                spliter = strtok(NULL,"\t");
                if(col_idx==5)
                    refs[i]=atoi(spliter);
                if(col_idx==6)
                    alts[i]=atoi(spliter);
            }
            col_idx = 0;
        }
        max_ks = 0; 
        p_0 = get_p_0(alts, refs, time_points); 
        for(int i=0; i<time_points; i++)
            p_mt[i]=alts[i]/(double)(alts[i]+refs[i]);
        for(int i=0; i<time_points; i++){
            ks[i]=fabs(p_0 - p_mt[i]);
            if(ks[i]>max_ks)
                max_ks = ks[i];
        }
        
        ks_pvals[allele_num] = get_pval(gsl_histogram_sum(gsl_hist[0]), max_ks);
    }
    
    for(int i=0; i<THREAD_NUM; i++){
        gsl_histogram_reset(gsl_hist[i]);
        gsl_histogram_free(gsl_hist[i]);
    }
    for(int i=0; i<allele_len; i++){
        fprintf(output_fp, "%.13lf\n", ks_pvals[i]);
    }
    
    fclose(output_fp);
    fclose(input_fp);
    gsl_rng_free(gsl_rng_bin);
    free(p_mt);
    free(ks);
    free(refs);
    free(alts);
    free(ks_pvals);
    free(spliter);
    return 0;
}
