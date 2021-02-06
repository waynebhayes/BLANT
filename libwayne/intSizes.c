#include <stdio.h>

#define PrintSz(x) printf("%s %ld\n", #x, sizeof(x));

int main(void) {
    PrintSz(char);
    PrintSz(short);
    PrintSz(int);
    PrintSz(long);
    PrintSz(long long);
    PrintSz(__int128);
    return 0;
}
