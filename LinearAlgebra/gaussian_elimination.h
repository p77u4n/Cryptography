#ifndef GAUSSIAN_ELIMINATION_H
#define GAUSSIAN_ELIMINATION_H

struct kernel_span{
    int **a;
    int num_sol;
};
struct kernel_span find_kernel_space(int **a ,int row,int col){
    int j,i,count,tmp;
    int boss = 0;
    for(j = 0 ; j < col ; j++){
        for(i = 0 ; i < row ; i++){
            if(a[i][j] == 1){
                if(boss == 0){
                    boss = i;
                }else{
                    for(count = 0 ; count < col ;count++){
                        a[i][count] = (a[i][count] + a[boss][count]) % 2;
                    }
                }
            }
            if(boss != 0 && boss != j){
                for(count = 0 ; count < col ; count++){
                    tmp = a[boss][count];
                    a[boss][count] = a[j][count];
                    a[j][count] = tmp;
                }
            }
            boss = 0;
        }
    }

    int pivot_x[col];
    int row_pivot[col];
    
    int num_pivot = 0;
    for(i = 0 ; i < row ; i++){
        for(j = i ; j < col ; j++){
            if(a[i][j] != 0){
                pivot_x[j] = 1;
                row_pivot[num_pivot] = i;
                num_pivot++;
            }
        }
    }
    int col_pivot[num_pivot];
    int col_free[col - num_pivot];
    j = 0;
    int j2 = 0;
    for(i = 0 ; i < col ; i++){
        if(pivot_x[i] == 1){
            col_pivot[j++] = i;

        }else{
            col_free[j2++] = i;
        }
    }
    int num_free = col - num_pivot;
    //int span[row][num_free] ;
    int **span = malloc(row * sizeof(int *) + (row * (num_free * sizeof(int))));
    for(j = 0 ; j < num_free ; j++){
        span[col_free[j]][j] = 1;
        for(i = 0 ; i < num_pivot ; i++){
            if(a[row_pivot[i]][col_free[j]] == 1){
                span[col_pivot[i]][j] = 1;
            }
        }
    }

    struct kernel_span result = {span,num_free};

    return result;
}

#endif
