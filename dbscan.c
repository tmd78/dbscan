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

struct args_dbscan
{
    double epsilon;
    int min_points;
    struct point *points;
};

struct args_get_neighbors
{
    double epsilon;
    struct point focal;
    struct point *points;
    int *S;
};

void create_points(struct point *points);
void dbscan(void *args);
int get_neighbors(void *args);

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
        fprintf(file_pointer, "%lf,%lf\n", points[i].x, points[i].y);
    }
    fclose(file_pointer);

    // Package arguments for dbscan.
    struct args_dbscan args = {
        .epsilon = epsilon,
        .min_points = min_points,
        .points = points
    };
    
    // Run dbscan.
    dbscan((void *)&args);

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

/*
Density-based spatial clustering of applications with noise.

Returns: void.
*/
void dbscan(void *args)
{
    double epsilon;
    int label;
    int min_points;
    struct point *points;
    // Seed set.
    int S[N*N];
    int S_count;
    
    // Unpack args.
    struct args_dbscan *packet = (struct args_dbscan *)args;
    epsilon = packet->epsilon;
    min_points = packet->min_points;
    points = packet->points;

    // Loop through all points.
    for (int i = 0; i < N; i++)
    {
        if (points[i].label != UNDEFINED)
        {
            // Don't revisit processed points.
            continue;
        }

        S_count = 0;

        // Package arguments for get_neighbors.
        struct args_get_neighbors args = {
            .epsilon = epsilon,
            .focal = points[i],
            .points = points,
            .S = S
        };

        // Find neighbors of points[i].
        int S_count = get_neighbors((void *)&args);
    }
}

/*
Find neighbors of focal and record their indices in S.

Returns: An int representing the number of neighbors found.
*/
int get_neighbors(void *args)
{
    double delta_x;
    double delta_y;
    double distance;
    double epsilon;
    struct point focal;
    int next_index;
    struct point *points;
    int *S;

    // Unpack args.
    struct args_get_neighbors *packet = (struct args_get_neighbors *)args;
    epsilon = packet->epsilon;
    focal = packet->focal;
    // The next available index in S.
    next_index = 0;
    points = packet->points;
    S = packet->S;
    
    // Loop through points to find neighbors--points within epsilon distance of focal.
    for (int i = 0; i < N; i++)
    {
        // Calculate distance.
        delta_x = focal.x - points[i].x;
        delta_y = focal.y - points[i].y;
        distance = sqrt(delta_x*delta_x + delta_y*delta_y);

        if (distance <= epsilon)
        {
            // We found a neighbor; record its index in S.
            S[next_index] = i;
            next_index += 1;
        }
    }

    return next_index;
}
