#include <stdio.h>

void dv_init_disk(double mass, double rate, double radius);
void dv_eval_globals();

float dv_query_f(char* key);

void main() {
    dv_init_disk(10,0.001,10);
    dv_eval_globals();

    printf("mass = %f\n", dv_query_f("mass"));
}

