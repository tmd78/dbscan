#include <math.h> // sqrt
#include <stdbool.h> // bool
#include <stdio.h> // FILE
#include <stdlib.h> // RAND_MAX

// The number of points to generate; do not exceed 2,147,483,647 as int is being used.
#define N 100
#define NOISE 0
#define UNDEFINED -1

struct point
{
    int index;
    // The label of the cluster this point belongs to.
    int label;
    float x;
    float y;
};

void create_points(struct point *points);

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        printf("\nThis program requires 2 arguments:\n1. The value of epsilon\n2. The minimum number of points that constitute a cluster\n\n");
        return 0;
    }

    double epsilon;
    FILE *file_pointer;
    // The minimum number of points that constitute a cluster.
    int min_points;
    struct point points[N];

    epsilon = atof(argv[1]);
    min_points = atoi(argv[2]);

    create_points(points);

    // Write points to disk.
    file_pointer = fopen("points.csv", "w");
    fprintf(file_pointer, "x,y\n");
    for (int i = 0; i < N; i++)
    {
        fprintf(file_pointer, "%.2f,%.2f\n", points[i].x, points[i].y);
    }
    fclose(file_pointer);

    // TODO: Implement dbscan.

    // TODO: Write clusters to disk.

    return 0;
}

void create_points(struct point *points)
{
    double max_coordinate = 1000;
    srand(72);

    for (int i = 0; i < N; i++)
    {
        points[i].index = i;
        points[i].label = UNDEFINED;
        points[i].x = max_coordinate * ((double)rand() / (double)RAND_MAX);
        points[i].y = max_coordinate * ((double)rand() / (double)RAND_MAX);
    }
}
