#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

#define FILE_LENGTH_LIMIT 100

const char *argp_program_version = "binSampling 1.0";
const char *argp_program_bug_address = "<dhehdqls@gmail.com>";

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
    for(int i = 0; i < time_points; i++)
        sum += counts[i];
    int ref_sum = 0;
    for(int i = 0; i < time_points; i++)
        sum += counts[i];
    return alt_sum/double(ref_sum + alt_sum);
}

double binomial_sampling(int n, double p){
    /* 
    *
    *
    */
    double rand_u = 0.0;
    srand(time(NULL));
    rand_u = (double)rand()/RAND_MAX;
    
    
}

int main(int argc, char* argv[]){
    printf("??");
    int time_points = 0;
    int sample_num = 0;
    char* input_name;// = malloc(sizeof(char) * FILE_LENGTH_LIMIT);
    char* output_name;// = malloc(sizeof(char) * FILE_LENGTH_LIMIT);
    FILE *input_fp;
    FILE *output_fp;
    int allele_len = -1;
    
    int opt;
    while((opt=getopt(argc, argv, "t:n:i:o:")) != -1 ){
        switch(opt){
            case 't':
                time_points = atoi(optarg);
                break;
            case 'n':
                sample_num = atoi(optarg);
                break;
            case 'i':
                input_name = malloc(sizeof(char)*sizeof(optarg));
                if(optarg==NULL)
                    printf("-i need a second argument as the input file name.\n");
                else
                    strcpy(input_name, optarg);
                break;
            case 'o':
                output_name = malloc(sizeof(char)*sizeof(optarg));
                if(optarg==NULL)
                    printf("-o need a second argument as the output file name.\n");
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
    }
    if(output_name == NULL){
        printf("option -o requires an output file name.\n");
    }
    
    printf("Time points : %d\n", time_points);
    printf("The input file : %s\n", input_name);
    printf("The output file : %s\n", output_name);
    printf("Sample num : %d\n", sample_num);
    
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
    
    printf("Count the total number of alleles . . . . \n");
    allele_len = check_file_len(input_name);
    printf("Number of alleles of %s : %d \n",input_name,allele_len);
    
    ///////////////////////////////////
    
    int* refs = malloc(sizeof(int)*time_points);
    int* alts = malloc(sizeof(int)*time_points);
    char* spliter;
    char lines[4096];
    int idx = 0;
    int col_idx =0;
    int alt_sum = 0;
    int ref_sum = 0;
    double p_0 = 0;
    double* p_mt = malloc(sizeof(double)*time_points);
    double* ks = malloc(sizeof(double)*time_points);
    double max_ks = 0;
    double** bin_samples = malloc(sizeof(double*)*time_points);
    int total_n = 0;
    
    while(1){
        for(int i=0; i<time_points; i++){
            if(fgets(lines, 4096, input_fp)==NULL){
                printf("File read error.\n");
                return -1;
            }
            spliter = strtok(lines,"\t");
            while(spliter != NULL){
                col_idx++;
//                printf("%s\n", spliter);
                spliter = strtok(NULL,"\t");
                if(col_idx==5)
                    refs[i]=atoi(spliter);
                if(col_idx==6)
                    alts[i]=atoi(spliter);
            }
//          printf("refs: %d, alts: %d\n",refs[i], alts[i]);
            col_idx = 0;
        }
        
        p_0 = get_p_0(alts, refs, time_points); 
        for(int i=0; i<time_points; i++)
            p_mt[i]=alts[i]/double(alts[i]+refs[i]);
        for(int i=0; i<time_points; i++){
            ks[i]=fabs(p_0 - p_mt[i]);
            if(ks[i]>max_ks)
                max_ks = ks[i];
        }
        
        // Binomial sampling for each timepoints
        for(int i=0; i<time_points; i++){
            bin_samples[i] = malloc(sizeof(double)*sample_num);
            total_n = refs[i] + alts[i];
            for(int j=0; j<sample_num; j++){
                bin_samples[i][j] = binomial_sampling(total_n, p_0);
            
            }
        }
        
        return 0;
    }
    
    free(p_mt);
    free(ks);
    free(refs);
    free(alts);
}