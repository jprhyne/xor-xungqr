#include<stdlib.h>
#include<string.h>
#include<stdio.h>

int main(int argc, char *argv[])
{
    // Run make to compile the test function if needed
    //system("make");
    char command[2048];
    // Now iterate over several options of m,n,k, and nb
    for (int m = 50; m <= 850; m+=100) {
        for (int n = 50; n <= m; n+=100) {
            for (int k = 10; k <= n; k+=100) {
                for (int nb = 3; nb <= 32 && nb < k; nb+=10) {
                    sprintf(command, "./test_v3.exe -m %d -n %d -k %d -nb %d -e", m,n,k,nb);
                    system(command);
                }
            }
        }
    } 
}
