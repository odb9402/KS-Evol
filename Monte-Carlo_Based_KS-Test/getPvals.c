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
#include <float.h>
#include <sys/sysinfo.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>

#define HIST_BIN_SIZE 1000000
#define HIST_BIAS 0 
#define FILE_LENGTH_LIMIT 100
#define THREAD_NUM 28
#define CLOCK_PER_MIN CLOCK_PER_SEC*60
#define KERN_SIZE 1001 


const char *argp_program_version = "binSampling 1.0";
const char *argp_program_bug_address = "<dhehdqls@gmail.com>";
gsl_rng* gsl_rng_bin; //GSL global random number generator
gsl_histogram* gsl_hist[THREAD_NUM]; //GSL global histogram allocator 
double empirical_cdf[HIST_BIN_SIZE];

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

double get_pval(double sampling_num, double ks){
    /*
    * Call an empirical p-value of ks-statistics from pre-constructed
    * gsl-histogram consisting of sampling values.
    *
    * Global variable : gsl_hist[0] (summarized gsl_histogram)
    *                   empirical_cdf (double array of empirical CDF of the sampling distribution)
    * return : An empirical p-value of ks distance 
    */
    size_t hist_idx = -1;
    double empirical_cdf_count = 0;
    double p_val = -1.0;
    const double shifted_bias = HIST_BIN_SIZE*HIST_BIAS/(double)sampling_num;
    
    if(gsl_histogram_find(gsl_hist[0], ks, &hist_idx) == GSL_EDOM){
        printf("GSL histogram find error. \n");
        return -1;
    }
    //printf("max ks : %lf\n", ks);
    //printf("ks_hist_idx :  %ld\n", hist_idx);
    
    //for(int i=0; i < HIST_BIN_SIZE; i++)
    //    printf("%lf \t", gsl_histogram_get(gsl_hist,i));
    
    //for(int i=0; i <= hist_idx; i++)
    //    empirical_cdf_count += gsl_histogram_get(gsl_hist[0], i);
    
    p_val = (1.0) - empirical_cdf[hist_idx];//empirical_cdf_count/sampling_num;
    if(p_val < 0.0)
        p_val = 1E-30;
    //printf("pval : %.13f\n\n",p_val);
    return p_val;
}

void gsl_ks_distance_sampling(void *args){
    /*
    * Build the sampling histogram to calculate p-value for each SNP.
    * This function will be executed parallelly with POSIX multithreads.
    *
    * Argument : Struct sampling_args (Including timepoints, a sampling number, null hypothesis,
    *            an average number of SNP and the thread number.)
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
    
    sampled_p = malloc(sizeof(double)*time_points);
    for(long i=0; i<sample_num/THREAD_NUM; i++){
        max_ks_bin = 0;
        
        for(int j=0; j<time_points; j++)
           sampled_p[j] = gsl_ran_binomial(gsl_rng_bin, p_0, average_allele)/(double)(average_allele);
        
        for(int j=0;j<time_points;j++){
            sampled_p[j] = fabs(sampled_p[j] - p_0); //No extra allocations.
            if(sampled_p[j] > max_ks_bin)
                max_ks_bin = sampled_p[j];
        }
        if(gsl_histogram_increment(gsl_hist[thread_num], max_ks_bin) == GSL_EDOM){
            printf("GSL histogram range error. \n");
            exit(-1);
        }
        if(i % 1000000 == 0 && thread_num == THREAD_NUM-1){
            if(i == 0)
                continue;
            printf("Sampling : [%ld / %ld]\n", (long)i*THREAD_NUM, sample_num);
        }
    }
}

double* gen_gaussian_filter(double sigma, int filter_size){
    double* gaussian_filter = malloc(sizeof(double)*filter_size);
    double filter_sum = 0.0;
    
    if(filter_size % 2 != 1){
        printf("Filter size should be an odd number.\n");
        exit(-1);
    }
    for(int x= -filter_size/2; x <= filter_size/2; x++){
        //printf("%d:",x+filter_size/2);
        gaussian_filter[x + filter_size/2] = exp(-(x*x)/(2*sigma*sigma))/(sqrt(2*M_PI)*sigma);
        filter_sum += gaussian_filter[x + filter_size/2];
        //printf("<%e,%e> \t", gaussian_filter[x + filter_size/2], filter_sum);
    }
    for(int i=0; i<filter_size; i++)
        gaussian_filter[i] = gaussian_filter[i]/filter_sum;
    return gaussian_filter;
}

void hist_smoothing(double* hist_values, double* kern){
    int nconv = HIST_BIN_SIZE + KERN_SIZE - 1;
    double* expanded = malloc(sizeof(double)*nconv);
    double kern_sum = 0.0;
    
    for(int i=0; i<nconv; i++)
        expanded[i] = 0.0;
    for(int i=0; i<HIST_BIN_SIZE; i++)
        expanded[i+KERN_SIZE/2] = hist_values[i];
    
    for(int i=0; i<HIST_BIN_SIZE; i++){
        kern_sum = 0.0;
        for(int j=0; j<KERN_SIZE; j++){
            kern_sum += expanded[i+j] * kern[j];
        }
        hist_values[i] = kern_sum;
    }
}


int main(int argc, char* argv[]){
    int time_points = 0;
    unsigned long sample_num = 0;
    char* input_name;// = malloc(sizeof(char) * FILE_LENGTH_LIMIT);
    char* output_name;// = malloc(sizeof(char) * FILE_LENGTH_LIMIT);
    FILE *input_fp;
    FILE *output_fp;
    int allele_len = -1;
    FILE *hist_out;
    hist_out = fopen("hist_out.txt","w");
    
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
    long ref_sum_total = 0;
    long alt_sum_total = 0;
    long ref_sum = 0;
    long alt_sum = 0;
    char* spliter;
    char* lines = NULL;
    int col_idx =0;
    double p_0 = 0;
    double* p_mt = malloc(sizeof(double)*time_points);
    double* ks = malloc(sizeof(double)*time_points);
    double max_ks = 0;
    double* ks_pvals = malloc(sizeof(double)*allele_len);
    double cdf_count = 0.0;
    double hist_sum = 0.0;
    double* hist_values = malloc(sizeof(double)*HIST_BIN_SIZE);
    size_t len = 0;
    
    const gsl_rng_type *T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    gsl_rng_bin = gsl_rng_alloc(T);
    for(int i=0; i<THREAD_NUM; i++){
        gsl_hist[i] = gsl_histogram_alloc(HIST_BIN_SIZE);
        gsl_histogram_set_ranges_uniform(gsl_hist[i], 0.0, KERN_SIZE);
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
    printf("The number of entire alleles : %ld\n", ref_sum_total+alt_sum_total);
    printf("The total number of average alleles for each timepoint : %ld\n",
            (alt_sum_total + ref_sum_total) / (allele_len));
    printf("Binomial sampling based on p0 . . . . \n");
    
    printf("Build the sampling distribution of KS distances . . .\n");
    //Binomial sample values to ks distances from p_0. 
    struct sampling_args thread_args[THREAD_NUM];
    for(int i=0; i<THREAD_NUM; i++){
        //gsl_ks_distance_sampling((alt_sum_total + ref_sum_total) / allele_len, p_0, time_points, sample_num);
        thread_args[i].average_allele = (alt_sum_total + ref_sum_total) / (allele_len);
        thread_args[i].p_0 = p_0;
        thread_args[i].time_points = time_points;
        thread_args[i].sample_num = sample_num;
        thread_args[i].thread_num = i;        
        if(pthread_create(&sampling_thread[i],NULL,gsl_ks_distance_sampling,(void *)&thread_args[i])!=0){
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
    gsl_histogram_shift(gsl_hist[0], HIST_BIAS); //
    
    if((sample_num - gsl_histogram_sum(gsl_hist[0]) > 0)){
        for(int i=1; i < (sample_num - gsl_histogram_sum(gsl_hist[0])); i++)
            continue;
    }
    
    printf("Build the empirical CDF of the sampling distribution. . .\n");
    for(int i=0; i<HIST_BIN_SIZE; i++)
        hist_values[i] = gsl_histogram_get(gsl_hist[0], i);
    
    //for(int i=0; i<3000; i++)
    //    printf("%lf\t",hist_values[i]);
    hist_smoothing(hist_values,gen_gaussian_filter(51.0,KERN_SIZE));
    for(int i=0; i<HIST_BIN_SIZE; i++)
        hist_sum += hist_values[i];
    
    //for(int i=0; i<3000; i++)
    //    printf("%lf\t",hist_values[i]);
    //return 4;
    
    for(int i=0; i<HIST_BIN_SIZE; i++){
        cdf_count += hist_values[i]/hist_sum;
        empirical_cdf[i] = cdf_count;
    }
    if(empirical_cdf[HIST_BIN_SIZE-1] != 1.0)
        empirical_cdf[HIST_BIN_SIZE-1] = 1.0;
    cdf_count = -1.0;
    
    printf("Calculate p-values based on the sampling distribution . . .\n");
    for(long allele_num=0; allele_num<allele_len; allele_num++){
        //File read and extract counts
        for(int i=0; i<time_points; i++){
            if(getline(&lines, &len, input_fp) == -1){
                printf("File read error. %ld`th allele\n", allele_num);
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
        
        ks_pvals[allele_num] = get_pval(hist_sum, max_ks);
        if(allele_num % 1000 == 0)
            if(allele_num != 0)
                printf("Get p-values : [%ld / %ld]\n", (long)allele_num, (long)allele_len);
    }
    
    for(int i=0; i < HIST_BIN_SIZE; i++)
        fprintf(hist_out, "%f\n", hist_values[i]);
    
    for(int i=0; i<THREAD_NUM; i++){
        gsl_histogram_reset(gsl_hist[i]);
        gsl_histogram_free(gsl_hist[i]);
    }
    for(int i=0; i<allele_len; i++){
        fprintf(output_fp, "%e\n", ks_pvals[i]);
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
