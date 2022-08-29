#include"main.h"

int main(int argc, char const* argv[]){
    int threadNumber = 1 ;
    for (int i = 0; i < argc; i++) {
    if (!strcmp(argv[i], "-t")) threadNumber = atoi(argv[i + 1]);
    }
    possible_cases(threadNumber);
    return 0;
}


