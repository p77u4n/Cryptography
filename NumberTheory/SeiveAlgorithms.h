// Implement some algorithms about Eratosthenes Seive
// Author : p77u4n
#ifndef SEIVE_ALGO_H
#define	SEIVE_ALGO_H
#include<stdio.h>
#include<stdlib.h>
#include<gmp.h>

struct erathos_sieve_result{
    int num_prime;
    unsigned long int * list_primes; 
};

struct erathos_sieve_result simple_eratos_sieve(int lim){
	short int is_prime[lim];
    
	int num_prime = 0;
	int iter;
	for(iter = 0 ; iter < lim ; iter++){
		is_prime[iter] = 1;
	}
	int j;
	int i;
	is_prime[0] = 0;
	for(iter = 1 ; iter < lim ; iter++){
		if(is_prime[iter] == 1){
			
			num_prime++;
			i = iter + 1;
			//printf("Detecting : %d \n",i);
			j = 2*i;
			while(j <= lim){
				is_prime[j - 1] = 0;
				j = j + i;
			}
		}
	}
    printf("Num of Primes %d \n" , num_prime);
    unsigned long int *primes = (unsigned long int*)calloc(num_prime,sizeof(unsigned long int));
	//unsigned long int primes[num_prime];
    //if declaring "int primes[num_prime],the result will be zzzzzz" ????
	int iter1 = 0;
	for(iter = 0 ; iter < lim ; iter++){
		if(is_prime[iter] == 1){
            //printf("tim ra %d\n",sizeof(unsigned long int));
			primes[iter1++] = (unsigned long int)iter + 1;
		}
	}

    struct erathos_sieve_result result = {num_prime,primes};
	
	return result;
}

void eratos_sieve_result_to_file(char *name_file,int lim){
    FILE * out = fopen(name_file,"w+");
    struct erathos_sieve_result result = simple_eratos_sieve(lim);
    unsigned long int * result_list = result.list_primes;
    int iter;
    for(iter = 0 ; iter < result.num_prime - 1; iter++){
        fprintf(out,"%d,",result_list[iter]);   
    }
    fprintf(out,"%d",result_list[result.num_prime - 1]);
    fclose(out);
}

#endif
