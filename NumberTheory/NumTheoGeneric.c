#ifndef NUMTHEOGENERIC_H
#define NUMTHEOGENERIC_H

#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <gmp.h>
#include "../../UsefulFunction/usefulfunc.h"
int euclid_algorithm(int,int);

struct Euclid_Extended_Return {
    int y;
    int alpha;
    int beta;
};

struct root_result{
    mpz_t s1;
    mpz_t s2;
};

struct Euclid_Extended_Return euclid_extended_algorithm(int,int);
/*
int main(void){
    struct Euclid_Extended_Return result;
    result = euclid_extended_algorithm(6,9);
    printf("alpha : %d,beta : %d,gcd : %d",result.alpha,result.beta,result.y);
    return 0;
}
*/
int euclid_algorithm(int a,int b){
    int tmp;
    if(b > a){
        tmp = a;
        a = b;
        b = tmp;
    }
    while(b > 0){
        tmp = b;
        b = a - (a/b)*b;
        a = tmp;
    }

    return a; 
}

struct Euclid_Extended_Return euclid_extended_algorithm(int x,int y){
    int X,Y;
    int alpha_1,alpha_2 = 0;
    int beta_1 = 0,beta_2;
    if(x < 0){
        X = -x;
        alpha_1 = -1;
    }else{
        X = x;
        alpha_1 = 1;
    }

    if(y < 0){
        Y = -y;
        beta_2 = -1;
    }else{
        Y = y;
        beta_2 = 1;
    }
    Y < X?:swap(&X,&Y);
    int quotient,remainder,temp;
    while(Y > 0){
        quotient = X/Y;
        remainder = X - Y*quotient;
        X = Y;
        Y = remainder;
        temp = alpha_1;
        alpha_1 = alpha_2;
        alpha_2 = temp - quotient*alpha_2;
        temp = beta_1;
        beta_1 = beta_2;
        beta_2 = temp - quotient*beta_2;
    }

    struct Euclid_Extended_Return result = {X,alpha_1,beta_1};
    return result;
}

double logtwo(unsigned long int x){
    return log((double)x)/log(2);
}

mp_bitcnt_t find_e(mpz_t p){
   return mpz_scan1(p,(mp_bitcnt_t)0);
}

// SHANK TONELLI ALGRITHM

struct root_result shanks_tonelli_alg(mpz_t a,mpz_t p){
    mpz_t exp,p_subtract_1,s,n,p_sub_1_div_2;
    mpz_init(n);
    mpz_init(p_sub_1_div_2);
    mpz_init(p_subtract_1);
    mpz_init(s);
    mpz_sub_ui(p_subtract_1,p,1);
    mpz_cdiv_q_ui(p_sub_1_div_2,p_subtract_1,2);
    mpz_init(exp);
    mpz_sub_ui(exp,p,1);
    mpz_cdiv_q_ui(exp,exp,2);
    mpz_t a_pow_p_1_div_2;
    mpz_init(a_pow_p_1_div_2);
    mpz_powm(a_pow_p_1_div_2,a,exp,p);
    mpz_add_ui(a_pow_p_1_div_2,a_pow_p_1_div_2,(unsigned long int)1);
    mpz_mod(a_pow_p_1_div_2,a_pow_p_1_div_2,p);
    struct root_result result;
    mpz_init(result.s1);
    mpz_init(result.s2);
    if(mpz_cmp_ui(a_pow_p_1_div_2,0) == 0){
        return result;
    }
    //mpz_set(s,p_subtract_1);
    // Tinh s va e
    mp_bitcnt_t e;
    e = mpz_scan1(p_subtract_1,(mp_bitcnt_t)0);
    mpz_t two_pow_e;
    mpz_init(two_pow_e);
    mpz_ui_pow_ui(two_pow_e,2,e);
    mpz_cdiv_q(s,p_subtract_1,two_pow_e);
    // Tim quadratic non_residue
    unsigned long int n_iter = 2;
    mpz_set_ui(n,n_iter);
    mpz_t is_non_residue_n;
    mpz_init(is_non_residue_n);
    mpz_powm(is_non_residue_n,n,p_sub_1_div_2,p);
    mpz_add_ui(is_non_residue_n,is_non_residue_n,(unsigned long int)1);

    while(1){
        if(mpz_divisible_p(is_non_residue_n,p) != 0){
            break;
        }else{
            n_iter++;
            mpz_set_ui(n,n_iter);
            mpz_powm(is_non_residue_n,n,p_sub_1_div_2,p);
            mpz_add_ui(is_non_residue_n,is_non_residue_n,1);
        }
    }

    mpz_t x,b,g,temp;
    unsigned long int r = e;
    unsigned long int m;
    mpz_init(x);
    mpz_init(b);
    mpz_init(g);
    //mpz_init(r);
    mpz_init(temp);
    mpz_add_ui(temp,s,(unsigned long int)1);
    mpz_cdiv_q_ui(temp,temp,(unsigned long int)2);
    mpz_powm(x,a,temp,p);
    mpz_powm(b,a,s,p);
    mpz_powm(g,n,s,p);
    //mpz_set_ui(r,e);
    //int iter;
    mpz_t two_pow_m,b_pow;
    mpz_init(two_pow_m);
    mpz_init(b_pow);
    for(m = 0;m < r ; m++ ){
        mpz_ui_pow_ui(two_pow_m,2,m);
        mpz_powm(b_pow,b,two_pow_m,p);
        if(mpz_cmp_ui(b_pow,1) == 0){
            //m = iter;
            break;
        }

    }
    while(m != 0){
        mpz_ui_pow_ui(temp,2,r-m-1);
        mpz_powm(temp,g,temp,p);
        mpz_mul(x,x,temp);
        mpz_mod(x,x,p);
        mpz_powm_ui(temp,temp,2,p);
        mpz_mul(b,b,temp);
        mpz_mod(b,b,p);
        mpz_set(g,temp);
        r = m;
        for(m = 0;m < r ; m++ ){
            mpz_ui_pow_ui(two_pow_m,2,m);
            mpz_powm(b_pow,b,two_pow_m,p);

            if(mpz_cmp_ui(b_pow,1) == 0){
                break;
            }
        }

    }
    mpz_set(result.s1,x);
    mpz_sub(result.s2,p,x);
    return result;
    
}   

#endif
