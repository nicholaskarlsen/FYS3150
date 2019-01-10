#include <iostream>
#include <stdio.h>

using namespace std;

int main(){

    int N = 2;
    double A[N][N];

    A[0][0] = 123;

    for(int i = 0; i<N; i++){
        printf(" %d ", i)        ;
        for(int j = 0; j<N; j++){
            printf("%d", j);
            printf("%d", A[i][j]);
        }
    }

    return 0;
}
