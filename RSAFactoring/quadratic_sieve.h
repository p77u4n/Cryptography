#ifndef QUAD_SIEVE_H
#define QUAD_SIEVE_H

#include "./pollard.h"
#include "../NumberTheory/SeiveAlgorithms.h"
#include "../LinearAlgebra/gaussian_elimination.h"
void quadratic_sieve_simple(mpz_t N,int B,int num_candidate){
    struct erathos_sieve_result result_erathos = simple_eratos_sieve(B);
    unsigned long int * list_prime = result_erathos.list_primes ;
    int num_prime = result_erathos.num_prime;
    mpz_t list_candidate[num_candidate] ;
    mpz_t list_x[num_candidate];
    mpz_t list_candidate_temp[num_candidate] ; 
    unsigned long int iter = 0 ;
    mpz_t root_N;
    mpz_init(root_N);
    mpz_root(root_N,N,2);
    mpz_add_ui(root_N,root_N,1);
    for(iter = 0 ; iter < num_candidate ; iter++){
        mpz_init(list_candidate[iter]);
        mpz_add_ui(list_candidate[iter],root_N,iter); 
        mpz_set(list_x[iter],list_candidate[iter]);
        mpz_powm_ui(list_candidate[iter],list_candidate[iter],2,N);
        
        mpz_set(list_candidate_temp[iter],list_candidate[iter]);
    }
    
    int matrix[num_candidate][num_prime + 1] ;
    mpz_t s_pow_2,mpz_temp;
    mpz_init(s_pow_2);
    mpz_init(mpz_temp);
    mpz_t mpz_prime;

    //unsigned long int s1,s2;
    struct root_result result_root;
    unsigned long int iter2,num_pass,temp_index,begin_sieve,temp;
    num_pass = 0;
    printf("build matrix \n");
    for(iter = 0 ; iter < num_prime ; iter++){
        printf("kiem tra %d\n",iter);
        mpz_set_ui(mpz_prime,list_prime[iter]);
        mpz_set_ui(s_pow_2,list_prime[iter]*list_prime[iter]);
        
        result_root = shanks_tonelli_alg(s_pow_2,mpz_prime);
        //s1 = mpz_get_ui(result_root.s1);
        //s2 = mpz_get_ui(result_root.s2);
        //m
        // mpz_out_str(stdout,10,result_root.s1);
        // mpz_out_str(stdout,10,result_root.s2);
        if(mpz_cmp(root_N,result_root.s1) < 0){
            mpz_sub(mpz_temp,result_root.s1,root_N);
            begin_sieve = mpz_get_ui(mpz_temp) - 1;
        }else{
            mpz_sub(mpz_temp,root_N,result_root.s1);
            mpz_cdiv_q_ui(mpz_temp,mpz_temp,list_prime[iter]);
            //temp = mpz_get_ui(temp);
            mpz_mul_ui(mpz_temp,mpz_temp,list_prime[iter]);
            mpz_add(mpz_temp,mpz_temp,result_root.s1);
            while(mpz_cmp(mpz_temp,root_N) <= 0){
                mpz_add_ui(mpz_temp,mpz_temp,list_prime[iter]);
            }
            printf("...\n");
            mpz_sub(mpz_temp,mpz_temp,root_N);
            begin_sieve = mpz_get_ui(mpz_temp) - 1;
        }

        //mpz_t root_pass[]
        
        for(iter2 = begin_sieve ; iter2 < num_candidate ;){
            mpz_set(mpz_temp,list_candidate_temp[iter2]);
            while(mpz_divisible_ui_p(mpz_temp,list_prime[iter]) != 0){
                mpz_cdiv_q_ui(mpz_temp,mpz_temp,list_prime[iter]);
                matrix[iter2][iter]++;
            } 
            if(mpz_cmp_ui(mpz_temp,1) == 0){
                matrix[iter2][num_prime] = -1;
                num_pass++;
            }
            mpz_set(list_candidate_temp[iter2],mpz_temp);
            iter2 += list_prime[iter];
        }

        if(mpz_cmp(root_N,result_root.s2) < 0){
            mpz_sub(mpz_temp,result_root.s2,root_N);
            begin_sieve = mpz_get_ui(mpz_temp) - 1;
        }else{
            mpz_sub(mpz_temp,root_N,result_root.s2);
            mpz_cdiv_q_ui(mpz_temp,mpz_temp,list_prime[iter]);
            //temp = mpz_get_ui(temp);
            mpz_mul_ui(mpz_temp,mpz_temp,list_prime[iter]);
            mpz_add(mpz_temp,mpz_temp,result_root.s2);
            while(mpz_cmp(mpz_temp,root_N) <= 0){
                mpz_add_ui(mpz_temp,mpz_temp,list_prime[iter]);
            }
            mpz_sub(mpz_temp,mpz_temp,root_N);
            begin_sieve = mpz_get_ui(mpz_temp) - 1;
        }

        //mpz_t root_pass[]
        
        for(iter2 = begin_sieve ; iter2 < num_candidate ;){
            mpz_set(mpz_temp,list_candidate_temp[iter2]);
            while(mpz_divisible_ui_p(mpz_temp,list_prime[iter]) != 0){
                mpz_cdiv_q_ui(mpz_temp,mpz_temp,list_prime[iter]);
                matrix[iter2][iter]++;
            } 
            if(mpz_cmp_ui(mpz_temp,1) == 0){
                matrix[iter2][num_prime] = -1;
                num_pass++;
            }
            mpz_set(list_candidate_temp[iter2],mpz_temp);
            iter2 += list_prime[iter];
        }
    }
    mpz_t root_pass[num_pass];
    int **matrix_sieved = malloc(num_prime * sizeof(int *) + (num_prime * (num_pass * sizeof(int))));

    printf("num_prime : %d\n",num_prime);
    printf("num_pass : %d\n",num_pass);
    int matrix_save[num_pass][num_prime];
    int index_q[num_pass];
    int row_count = 0;
    for(iter = 0 ; iter < num_candidate ; iter++){
        if(matrix[iter][num_prime] == -1){
            for(iter2 = 0 ; iter2 < num_prime ; iter2++){
                matrix_sieved[iter2][row_count] = matrix[iter][iter2] % 2;
                matrix_save[row_count][iter2] = matrix[iter][iter2];
            }
            index_q[row_count] = iter;
            row_count++;
        }
    }
    int iter3;
    int pow_prime[num_prime];
    mpz_t root2,root1;
    mpz_init(root2);
    mpz_init(root1);
    struct kernel_span span_sol = find_kernel_space(matrix_sieved,num_prime,num_pass);
    int **span_matrix = span_sol.a;
    for(iter = 0 ; iter < span_sol.num_sol ; iter++){
        mpz_set_ui(root1,(unsigned long int)1);
        mpz_set_ui(root2,(unsigned long int) 1);
        for(iter2 = 0 ; iter2 < num_pass ; iter2++){
            if(span_matrix[iter2][iter] == 1){
                mpz_mul(root1,root1,list_x[index_q[iter2]]);
                for(iter3 = 0 ; iter3 < num_prime ;iter3++){
                    pow_prime[iter3] += matrix_save[iter2][iter3];
                }
            }

        }

        
        for(iter2 = 0 ; iter2 < num_prime ; iter2++){
            pow_prime[iter2] = pow_prime[iter2]/2;
            mpz_ui_pow_ui(mpz_temp,list_prime[iter2],pow_prime[iter2]);
            mpz_mul(root2,root2,mpz_temp);
        }

        mpz_sub(mpz_temp,root1,root2);
        mpz_t q;
        mpz_gcd(q,mpz_temp,N);
        if(mpz_cmp_ui(q,1) != 0){
            printf("done\n");
            break;
        }
    }

}
#endif
