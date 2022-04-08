#include "cvector.h"
#include <math.h> // sqrt
#include <stdbool.h> // bool
#include <stdio.h> // FILE
#include <stdlib.h> // RAND_MAX

#define CVECTOR_LOGARITHMIC_GROWTH
// The number of points to generate; do not exceed 2,147,483,647 as int is being used.
#define N 1000
#define NOISE 0
#define UNDEFINED -1

struct cell
{
    // Holds the indices of points contained in this cell.
    cvector_vector_type(int) points;
    int x;
    int y;
};

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

struct args_dense_box
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

struct args_label_points
{
    cvector_vector_type(int) indices;
    int label;
    struct point *points;
};

struct args_merge_dense_boxes
{
    int center;
    struct cell *G;
    int G_x;
    int G_y;
    int min_points;
    struct point *points;
    int *queue;
};

void create_points(struct point *points);
void dbscan(void *args);
void dense_box(void *args);
int get_neighbors(void *args);
void label_points(void *args);
void merge_dense_boxes(void *args);

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

    // Package arguments for dense_box.
    struct args_dense_box args_dense_box = {
        .epsilon = epsilon,
        .min_points = min_points,
        .points = points
    };

    // Avoid distance calculations with dense box.
    dense_box((void *)&args_dense_box);

    // Package arguments for dbscan.
    struct args_dbscan args_dbscan = {
        .epsilon = epsilon,
        .min_points = min_points,
        .points = points
    };
    
    // Run dbscan.
    dbscan((void *)&args_dbscan);

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
    label = 0;
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

        // Iterator.
        int j = 0;
        int neighbor_count;
        // Holds indices from get_neighbors.
        int neighbors[N];
        // The points index of seed in iteration below.
        int s;
        // Loop through S; note that S_count may increase during iteration.
        while (j < S_count)
        {
            s = S[j];
            
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

            neighbor_count = get_neighbors((void *)&args);

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
Avoid distance calculations by assigning points to clusters with the dense box algorithm.

Returns: void.
*/
void dense_box(void *args)
{
    struct args_label_points args_label_points;
    struct args_merge_dense_boxes args_merge_dense_boxes;
    double cell_length;
    double epsilon;
    // Number of cells in grid's x-dimension.
    int G_x;
    // Number of cells in grid's y-dimension.
    int G_y;
    int i;
    int label;
    int min_points;
    double min_x;
    double min_y;
    double max_x;
    double max_y;
    struct point *points;
    cvector_vector_type(int) queue;
    int x;
    int y;

    // Unpack args.
    struct args_dense_box *packet = (struct args_dense_box *)args;
    epsilon = packet->epsilon;
    min_points = packet->min_points;
    points = packet->points;

    // Find bounds for grid.
    min_x = points[0].x;
    min_y = points[0].y;
    max_x = points[0].x;
    max_y = points[0].y;
    for (i = 1; i < N; i++)
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

    // Find grid dimensions.
    cell_length = epsilon / (2 * sqrt(2));
    G_x = (max_x - min_x) / cell_length + 1;
    G_y = (max_y - min_y) / cell_length + 1;
    
    // Grid whose cells are potential dense boxes.
    struct cell *G = malloc((G_x * G_y) * sizeof(struct cell));
    for (x = 0; x < G_x; x++)
    {
        for (y = 0; y < G_y; y++)
        {
            G[x * G_y + y].x = x;
            G[x * G_y + y].y = y;
            
            // Allocate this cell's points vector.
            G[x * G_y + y].points = NULL;
            cvector_push_back(G[x * G_y + y].points, 0);
            cvector_pop_back(G[x * G_y + y].points);
        }
    }

    // Assign each point to a cell in grid.
    for (i = 0; i < N; i++)
    {
        x = (points[i].x - min_x) / cell_length;
        y = (points[i].y - min_y) / cell_length;
        cvector_push_back(G[x * G_y + y].points, i);
    }

    // Allocate queue.
    queue = NULL;
    cvector_push_back(queue, 0);
    cvector_pop_back(queue);
    
    // Update arguments.
    args_label_points.points = points;
    args_merge_dense_boxes.G = G;
    args_merge_dense_boxes.G_x = G_x;
    args_merge_dense_boxes.G_y = G_y;
    args_merge_dense_boxes.min_points = min_points;
    args_merge_dense_boxes.points = points;
    args_merge_dense_boxes.queue = queue;

    // Merge dense boxes.
    label = 0;
    for (i = 0; i < (G_x * G_y); i++)
    {
        if (cvector_size(G[i].points) < min_points)
        {
            // Skip this cell if not a dense box.
            continue;
        }

        if (cvector_size(G[i].points) == 0)
        {
            // Skip this dense box if there's no points to process.
            continue;
        }

        if (points[G[i].points[0]].label > 0)
        {
            // Skip this dense box if already merged.
            continue;
        }

        // Label all points in this cell with label's current value.
        label += 1;
        args_label_points.indices = G[i].points;
        args_label_points.label = label;
        label_points((void *)&args_label_points);

        // Update center in arguments.
        args_merge_dense_boxes.center = i;

        // TODO: Debug merge_dense_boxes().
        //merge_dense_boxes((void *)&args_merge_dense_boxes);

        // TODO: Continue merge with current label value and cells on queue.
    }

    // Free heap memory.
    cvector_free(queue);
    for (i = 0; i < (G_x * G_y); i++)
    {
        cvector_free(G[i].points);
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
    neighbors = packet->neighbors;
    points = packet->points;
    
    // Loop through points to find neighbors--points within epsilon distance of focal.
    next_index = 0;
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

/*
Label each point whose index is in indices.

Returns: void.
*/
void label_points(void *args)
{
    int index;
    // The indices (to points array) of points to label.
    int *indices;
    int label;
    struct point *points;

    // Unpack args.
    struct args_label_points *packet = (struct args_label_points *)args;
    indices = packet->indices;
    label = packet->label;
    points = packet->points;

    for (int i = 0; i < cvector_size(indices); i++)
    {
        index = indices[i];
        points[index].label = label;
    }
}

/*
A description.

Returns: void.
*/
void merge_dense_boxes(void *args)
{
    struct args_label_points args_label_points;
    // An index to G.
    int center;
    struct cell *G;
    int G_x;
    int G_y;
    bool in_bounds;
    int min_points;
    struct point *points;
    int *queue;
    int x;
    int y;
    
    // Unpack args.
    struct args_merge_dense_boxes *packet = (struct args_merge_dense_boxes *)args;
    center = packet->center;
    G = packet->G;
    G_x = packet->G_x;
    G_y = packet->G_y;
    min_points = packet->min_points;
    points = packet->points;
    queue = packet->queue;
    
    // Update arguments.
    args_label_points.label = points[G[center].points[0]].label;
    args_label_points.points = points;
    
    // Check top.
    x = G[center].x;
    y = G[center].y + 1;
    in_bounds = (x > -1) && (x < G_x) && (y > -1) && (y < G_y);
    if (in_bounds && cvector_size(G[x * G_y + y].points) >= min_points)
    {
        // Label all points in this dense box.
        args_label_points.indices = G[x * G_y + y].points;
        label_points((void *)&args_label_points);

        // Queue cell.
        cvector_push_back(queue, x * G_y + y);
    }
    else
    {
        x = G[center].x;
        y = G[center].y + 2;
        in_bounds = (x > -1) && (x < G_x) && (y > -1) && (y < G_y);
        if (in_bounds && cvector_size(G[x * G_y + y].points) >= min_points)
        {
            // TODO: Attempt to merge non-adjacent dense box.
        }
    }

    // Check top-right.
    x = G[center].x + 1;
    y = G[center].y + 1;
    in_bounds = (x > -1) && (x < G_x) && (y > -1) && (y < G_y);
    if (in_bounds && cvector_size(G[x * G_y + y].points) >= min_points)
    {
        // Label all points in this dense box.
        args_label_points.indices = G[x * G_y + y].points;
        label_points((void *)&args_label_points);

        // Queue cell.
        cvector_push_back(queue, x * G_y + y);
    }

    // Check right.
    x = G[center].x + 1;
    y = G[center].y;
    in_bounds = (x > -1) && (x < G_x) && (y > -1) && (y < G_y);
    if (in_bounds && cvector_size(G[x * G_y + y].points) >= min_points)
    {
        // Label all points in this dense box.
        args_label_points.indices = G[x * G_y + y].points;
        label_points((void *)&args_label_points);

        // Queue cell.
        cvector_push_back(queue, x * G_y + y);
    }
    else
    {
        x = G[center].x + 2;
        y = G[center].y;
        in_bounds = (x > -1) && (x < G_x) && (y > -1) && (y < G_y);
        if (in_bounds && cvector_size(G[x * G_y + y].points) >= min_points)
        {
            // TODO: Attempt to merge non-adjacent dense box.
        }
    }

    // Check bottom-right.
    x = G[center].x + 1;
    y = G[center].y - 1;
    in_bounds = (x > -1) && (x < G_x) && (y > -1) && (y < G_y);
    if (in_bounds && cvector_size(G[x * G_y + y].points) >= min_points)
    {
        // Label all points in this dense box.
        args_label_points.indices = G[x * G_y + y].points;
        label_points((void *)&args_label_points);

        // Queue cell.
        cvector_push_back(queue, x * G_y + y);
    }

    // Check bottom.
    x = G[center].x;
    y = G[center].y - 1;
    in_bounds = (x > -1) && (x < G_x) && (y > -1) && (y < G_y);
    if (in_bounds && cvector_size(G[x * G_y + y].points) >= min_points)
    {
        // Label all points in this dense box.
        args_label_points.indices = G[x * G_y + y].points;
        label_points((void *)&args_label_points);

        // Queue cell.
        cvector_push_back(queue, x * G_y + y);
    }
    else
    {
        x = G[center].x;
        y = G[center].y - 2;
        in_bounds = (x > -1) && (x < G_x) && (y > -1) && (y < G_y);
        if (in_bounds && cvector_size(G[x * G_y + y].points) >= min_points)
        {
            // TODO: Attempt to merge non-adjacent dense box.
        }
    }

    // Check bottom-left.
    x = G[center].x - 1;
    y = G[center].y - 1;
    in_bounds = (x > -1) && (x < G_x) && (y > -1) && (y < G_y);
    if (in_bounds && cvector_size(G[x * G_y + y].points) >= min_points)
    {
        // Label all points in this dense box.
        args_label_points.indices = G[x * G_y + y].points;
        label_points((void *)&args_label_points);

        // Queue cell.
        cvector_push_back(queue, x * G_y + y);
    }

    // Check left.
    x = G[center].x - 1;
    y = G[center].y;
    in_bounds = (x > -1) && (x < G_x) && (y > -1) && (y < G_y);
    if (in_bounds && cvector_size(G[x * G_y + y].points) >= min_points)
    {
        // Label all points in this dense box.
        args_label_points.indices = G[x * G_y + y].points;
        label_points((void *)&args_label_points);

        // Queue cell.
        cvector_push_back(queue, x * G_y + y);
    }
    else
    {
        x = G[center].x - 2;
        y = G[center].y;
        in_bounds = (x > -1) && (x < G_x) && (y > -1) && (y < G_y);
        if (in_bounds && cvector_size(G[x * G_y + y].points) >= min_points)
        {
            // TODO: Attempt to merge non-adjacent dense box.
        }
    }

    // Check top-left.
    x = G[center].x - 1;
    y = G[center].y + 1;
    in_bounds = (x > -1) && (x < G_x) && (y > -1) && (y < G_y);
    if (in_bounds && cvector_size(G[x * G_y + y].points) >= min_points)
    {
        // Label all points in this dense box.
        args_label_points.indices = G[x * G_y + y].points;
        label_points((void *)&args_label_points);

        // Queue cell.
        cvector_push_back(queue, x * G_y + y);
    }
}
