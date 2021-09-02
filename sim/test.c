#include <stdio.h>
#include <stdlib.h>

typedef struct tracer {
  int tracerid;
  double pos[3];
  int tracertype;
} tracer;

#define ROOT_DIR "."
#define INPUT_DIR "input"

void main(void){
    tracer tracer_array[5];

    char fname[30];
    FILE * fp;
    sprintf(fname, "%s/%s/tracer.csv", ROOT_DIR, INPUT_DIR);   
    fp = fopen (fname, "r");

    for (int i=0;i<5;i++){
        fscanf(fp, "%d, %lf, %lf, %lf, %d", &tracer_array[i].tracerid, 
        &tracer_array[i].pos[0], &tracer_array[i].pos[1], &tracer_array[i].pos[2],
        &tracer_array[i].tracertype);
    }

    fclose(fp);

    for (int i=0;i<5;i++){
      printf("%d, %lf, %lf, %lf, %d\n", tracer_array[i].tracerid, tracer_array[i].pos[0],
       tracer_array[i].pos[1], tracer_array[i].pos[2], tracer_array[i].tracertype);
    }
}