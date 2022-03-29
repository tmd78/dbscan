#include <math.h> // sqrt
#include <stdbool.h> // bool
#include <stdio.h> // FILE
#include <stdlib.h> // RAND_MAX

// The number of points to generate; do not exceed 2,147,483,647 as int is being used.
#define N 1000
#define NOISE 0
#define UNDEFINED -1

struct point
{
    int index;
    // The label of the cluster this point belongs to.
    int label;
    double x;
    double y;
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
    int *neighbors;
    struct point *points;
};

void create_grid(struct point *points);
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

    create_grid(points);

    // TODO: Run dense box algorithm.

    // Package arguments for dbscan.
    struct args_dbscan args = {
        .epsilon = epsilon,
        .min_points = min_points,
        .points = points
    };
    
    // Run dbscan.
    dbscan((void *)&args);

    // Write clusters to disk.
    file_pointer = fopen("clusters.csv", "w");
    fprintf(file_pointer, "x,y,cluster\n");
    for (int i = 0; i < N; i++)
    {
        fprintf(file_pointer, "%lf,%lf,%d\n", points[i].x, points[i].y, points[i].label);
    }
    fclose(file_pointer);

    return 0;
}

void create_grid(struct point *points)
{
    double min_x;
    double min_y;
    double max_x;
    double max_y;

    // Find bounds for grid.
    min_x = points[0].x;
    min_y = points[0].y;
    max_x = points[0].x;
    max_y = points[0].y;
    for (int i = 1; i < N; i++)
    {
        // Update min bounds if needed.
        if (points[i].x < min_x)
        {
            min_x = points[i].x;
        }
        if (points[i].y < min_y)
        {
            min_y = points[i].y;
        }
        // Update max bounds if needed.
        if (points[i].x > max_x)
        {
            max_x = points[i].x;
        }
        if (points[i].y > max_y)
        {
            max_y = points[i].y;
        }
    }
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
    // Reusable; holds indices from get_neighbors.
    int neighbors[N];
    struct point *points;
    // Seed set.
    int S[N*N];
    int S_count;
    
    // Unpack args.
    struct args_dbscan *packet = (struct args_dbscan *)args;
    epsilon = packet->epsilon;
    label = 0;
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
            .neighbors = S,
            .points = points
        };

        // Find neighbors of points[i].
        int S_count = get_neighbors((void *)&args);

        if (S_count < min_points)
        {
            // min_points condition does not hold.
            points[i].label = NOISE;
            continue;
        }

        // Add points[i] to cluster whose ID = label.
        label += 1;
        points[i].label = label;

        // Loop through S; note that S_count may increase during iteration.
        int j = 0;
        while (j < S_count)
        {
            // The points index of the seed in this iteration.
            int s = S[j];
            
            // Increment j to show we've visited this seed.
            j += 1;

            if (s == i)
            {
                // Don't revisit points[i].
                continue;
            }

            if (points[s].label == NOISE)
            {
                // This is a border point.
                points[s].label = label;
                continue;
            }

            if (points[s].label != UNDEFINED)
            {
                // Don't revisit processed points.
                continue;
            }

            // Add this seed to cluster whose ID = label.
            points[s].label = label;

            // Package arguments for get_neighbors.
            struct args_get_neighbors args = {
                .epsilon = epsilon,
                .focal = points[s],
                .neighbors = neighbors,
                .points = points
            };

            int neighbor_count = get_neighbors((void *)&args);

            if (neighbor_count < min_points)
            {
                // min_points condition does not hold.
                continue;
            }

            // Add this seed's neighbor indices to S.
            for (int k = 0; k < neighbor_count; k++)
            {
                S[S_count] = neighbors[k];
                S_count += 1;
            }
        }
    }
}

/*
Find neighbors of focal and record their indices in neighbors.

Returns: An int representing the number of neighbors found.
*/
int get_neighbors(void *args)
{
    double delta_x;
    double delta_y;
    double distance;
    double epsilon;
    struct point focal;
    int *neighbors;
    // The next available index in neighbors.
    int next_index;
    struct point *points;

    // Unpack args.
    struct args_get_neighbors *packet = (struct args_get_neighbors *)args;
    epsilon = packet->epsilon;
    focal = packet->focal;
    next_index = 0;
    neighbors = packet->neighbors;
    points = packet->points;
    
    // Loop through points to find neighbors--points within epsilon distance of focal.
    for (int i = 0; i < N; i++)
    {
        // Calculate distance.
        delta_x = focal.x - points[i].x;
        delta_y = focal.y - points[i].y;
        distance = sqrt(delta_x*delta_x + delta_y*delta_y);

        if (distance <= epsilon)
        {
            // We found a neighbor; record its index.
            neighbors[next_index] = i;
            next_index += 1;
        }
    }

    return next_index;
}
