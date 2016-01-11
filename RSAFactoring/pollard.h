#ifndef POLLARD_H
#define POLLARD_H
#include<stdint.h>
#include<inttypes.h>
#include<gmp.h>
#include<stdio.h>
#include "../NumberTheory/SeiveAlgorithms.h"
#include "../NumberTheory/NumTheoGeneric.h"

struct factoring_result{
    mpz_t q;
    mpz_t p;
};

void F_pollard_rho(mpz_t in,mpz_t out,mpz_t N,int c){
    //mpz_t result;
   // mpz_init(result);
    mpz_powm_ui(out,in,2,N);
    mpz_add_ui(out,out,c);
    mpz_mod(out,out,N);
    
    //return result;
}

struct factoring_result pollard_rho(mpz_t N,int c,int seed,int num_bit_p){
    //printf("done ");
    mpz_t x0,x,x_dot,p,q,subtract;
    mpz_init(subtract);
    mpz_init(x0);
    mpz_init(p);
    mpz_init(q);
    mpz_init(x);
    mpz_init(x_dot);
    
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state,seed);
    mpz_urandomm(x0,state,N);
    mpz_set(x,x0);
    mpz_set(x_dot,x0);
    struct factoring_result result;
    mpz_init(result.q);
    mpz_init(result.p);
   // printf("done ");
    uint64_t iter = 1; 
    iter = iter << (num_bit_p/2);
    while(iter > 0){
        F_pollard_rho(x,x,N,c);
        F_pollard_rho(x_dot,x_dot,N,c);
        F_pollard_rho(x_dot,x_dot,N,c);
        mpz_sub(subtract,x,x_dot); 
        mpz_gcd(p,subtract,N);
        iter--;
        if(mpz_cmp_d(p,1) != 0 && mpz_cmp(p,N) != 0){
            mpz_cdiv_q(q,N,p);
            mpz_set(result.p,p);
            mpz_set(result.q,q);
            return result;
        }
    }
    return result;

}

struct factoring_result pollard_p_1(mpz_t N,int bound,mpz_t a,int bit_p){
    mpz_t gcdN_a,p,q,B,y;
    //mpz_init(N);
    //mpz_init(a);
    mpz_init(gcdN_a);
    mpz_init(p);
    mpz_init(q);
    mpz_init(B);
    mpz_init(y);
    struct factoring_result result;
    mpz_init(result.p);
    mpz_init(result.q);

    char buff[255];
    //FILE *input = fopen("/home/cloudatlas/CodeRespository/ModernCryptography/RSAFactoring/input.txt","r");
    //mpz_inp_str(N,input,10);
    //gmp_randstate_t state;
    //gmp_randinit_default(state);
    
    //mpz_urandomm(a,state,N); 
    mpz_gcd(gcdN_a,N,a);
    if(mpz_cmp_d(gcdN_a,1) != 0){
        mpz_init(q);
        mpz_init(p);
        mpz_set(q,a);
        mpz_cdiv_q(p,N,q);
        mpz_set(result.p,p);
        mpz_set(result.q,q);
        return result;
    }
    
    mpz_set_ui(B,1);
    struct erathos_sieve_result result_sieve = simple_eratos_sieve(bound);
    unsigned long int * list_prime = result_sieve.list_primes;
    int num_primes = result_sieve.num_prime;
    mpz_t temp;
    mpz_init(temp);
    unsigned long int exp;
    int iter;
    //printf("%d\n",logtwo(bound));
    for(iter = 0 ; iter < num_primes ; iter++){
        exp = (unsigned long int)(logtwo(bit_p)/logtwo(list_prime[iter]));
       // exp = (unsigned long int)(logtwo(bound)/logtwo(list_prime[iter]));
        //printf("exp : %d\n",exp);
        mpz_ui_pow_ui(temp,list_prime[iter],exp);
        mpz_mul(B,B,temp);
    }
    /* 
    mpz_out_str(stdout,10,N);
    printf(" is N \n");
    */
    //mpz_out_str(stdout,10,B);
    //printf(" is B \n");
    
    //mpz_out_str(stdout,10,a);
   // printf(" is a \n");
    
    mpz_powm(y,a,B,N);
    mpz_sub_ui(y,y,1);
    // mpz_out_str(stdout,10,y);
    //printf(" is y \n");
    mpz_gcd(p,y,N);
    mpz_cdiv_q(q,N,p);
    //printf("Done pollard \n");
    mpz_set(result.p,p);
    mpz_out_str(stdout,10,p);
    printf("is p \n");
    mpz_set(result.q,q);
    //printf("Done pollard \n");
    return result;
}


#endif
