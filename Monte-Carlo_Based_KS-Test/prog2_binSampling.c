#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <stddef.h>

#define NTHREADS 15

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define FILE_LENGTH_LIMIT 100


const char *argp_program_version = "binSampling 1.0";
const char *argp_program_bug_address = "<dhehdqls@gmail.com>";

gsl_rng* gsl_rng_bin; //GSL global random number generator

void print_usage(){
    printf("\t<USAGE of binSampling>\n");
    printf("\tbinSampling will generate simulation data which follows a binomial distribution\n");
    printf("\tbased on a null hypothesis of the time-dependent allele frequency data.\n\n");
    printf("\tbinSampling -n <number of samplings> -t <timepoints> -i <input> -o <output> \n\n");
}

int check_file_len(char* filename){
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
    int alt_sum = 0;
    int ref_sum = 0;
    
    for(int i = 0; i < time_points; i++){
        alt_sum += alts[i];
        ref_sum += refs[i];
    }
    
    return alt_sum/(double)(ref_sum + alt_sum);
}


int main(int argc, char* argv[]){
    int time_points = 0;
    int sample_num = 0;
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
                sample_num = atoi(optarg);
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
    printf("Sample num : %d\n", sample_num);
    
    printf("Count the total number of alleles . . . . \n\n");
    allele_len = check_file_len(input_name);
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
    char* spliter;
    char* lines = NULL;
    int col_idx =0;
    double p_0 = 0;
    double* p_mt = malloc(sizeof(double)*time_points);
    double* ks = malloc(sizeof(double)*time_points);
    double max_ks = 0;
    double** bin_samples = malloc(sizeof(double*)*time_points);
    double* bin_samples_ks = malloc(sizeof(double)*sample_num);
    int total_n = 0;
    clock_t start, end;
    size_t len = 0;
    
    const gsl_rng_type *T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    gsl_rng_bin = gsl_rng_alloc(T);
    
    for(int i=0; i<time_points; i++){
        bin_samples[i] = malloc(sizeof(double)*sample_num);
    }
    
    ////////////////// KS distance sampling loops /////////////////////////
    start = clock();
    for(int allele_num=0; allele_num<(allele_len/time_points); allele_num++){
        //File read and extract counts
        for(int i=0; i<time_points; i++){
            //if(fgets(lines, 8193, input_fp)==NULL){
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
        
        p_0 = get_p_0(alts, refs, time_points); 
        for(int i=0; i<time_points; i++)
            p_mt[i]=alts[i]/(double)(alts[i]+refs[i]);
        for(int i=0; i<time_points; i++){
            ks[i]=fabs(p_0 - p_mt[i]);
            if(ks[i]>max_ks)
                max_ks = ks[i];
        }
        
        // Binomial sampling for each timepoints
        // bin_samples['time_points']['sample_id']
        for(int i=0; i<time_points; i++){
            total_n = refs[i] + alts[i];
            for(int j=0; j<sample_num; j++){
                bin_samples[i][j] = gsl_ran_binomial(gsl_rng_bin, p_0, total_n)/(double)total_n;
            }
        }
        
        //Binomial sample values to ks distances from p_0. 
        for(int i=0; i<sample_num; i++){
            max_ks = 0;
            for(int j=0;j<time_points;j++){
                bin_samples[j][i] = fabs(bin_samples[j][i] - p_0); // No extra allocations.
                if(bin_samples[j][i] > max_ks)
                    max_ks = bin_samples[j][i];
            }
            bin_samples_ks[i] = max_ks;
        }
        
        
        for(int i=0; i<sample_num; i++){
            fprintf(output_fp, "%f\t", bin_samples_ks[i]);
        }
        fprintf(output_fp, "\n");
        
        if((allele_num % 1000) == 0 && allele_num != 0){
            end = clock();
            printf("Sampling for the %d`th allele . . . .  :%f", allele_num, (double)(end-start)/CLOCKS_PER_SEC);
            printf("\tETA : %fsec\n", (allele_len/(1000*time_points))*(double)(end-start)/CLOCKS_PER_SEC);
            start = clock();
        }
    }
    
    fclose(output_fp);
    fclose(input_fp);
    free(p_mt);
    free(ks);
    free(refs);
    free(alts);
    free(bin_samples_ks);
    for(int i=0; i<time_points; i++)
        free(bin_samples[i]);
    free(bin_samples);
    free(spliter);
    gsl_rng_free(gsl_rng_bin);
    return 0;
}