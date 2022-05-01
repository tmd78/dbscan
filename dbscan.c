#include "cvector.h"
#include <math.h> // sqrt
#include "rtree.h"
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
    double epsilon;
    int focal;
    struct cell *G;
    int G_x;
    int G_y;
    int label;
    int *merged;
    int min_points;
    struct point *points;
};

struct args_rtree_search
{
    int *candidates;
    int *candidates_count;
};

struct args_should_merge
{
    struct cell a;
    struct cell b;
    double epsilon;
    struct point *points;
};

void create_points(struct point *points);
void dbscan(void *args);
void dense_box(void *args);
int get_neighbors(void *args);
void label_points(void *args);
int merge_dense_boxes(void *args);
bool rtree_search_item_handler(const double *mbr, const void *item, void *args);
bool should_merge(void *args);

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
    // Allocate on heap.
    struct point *points = malloc(N * sizeof(struct point));

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
    double max_coordinate = N;
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
    struct args_get_neighbors args_get_neighbors;
    struct args_rtree_search args_rtree_search;
    int *candidates = malloc(N * sizeof(int));
    int candidates_count;
    double epsilon;
    int label;
    double mbr[4];
    int min_points;
    int *neighbors = malloc(N * sizeof(int));
    int neighbors_count;
    struct point *points;
    struct rtree *r_tree = rtree_new(sizeof(struct point *), 2);
    // Seed set.
    cvector_vector_type(int) S = NULL;
    cvector_push_back(S, 0);
    cvector_pop_back(S);
    
    // Unpack args.
    struct args_dbscan *packet = (struct args_dbscan *)args;
    epsilon = packet->epsilon;
    min_points = packet->min_points;
    points = packet->points;

    // Populate r_tree.
    for (int i = 0; i < N; i++)
    {
        mbr[0] = points[i].x;
        mbr[1] = points[i].y;
        mbr[2] = points[i].x;
        mbr[3] = points[i].y;
        rtree_insert(r_tree, mbr, &points[i]);
    }

    // Update arguments for get_neighbors.
    args_get_neighbors.epsilon = epsilon;
    args_get_neighbors.neighbors = neighbors;
    args_get_neighbors.points = points;

    // Update arguments for rtree_search.
    args_rtree_search.candidates = candidates;

    // Loop through all points.
    label = 0;
    for (int i = 0; i < N; i++)
    {
        if (points[i].label != UNDEFINED)
        {
            // Don't revisit processed points.
            continue;
        }

        // Update arguments for rtree_search.
        candidates_count = 0;
        args_rtree_search.candidates_count = &candidates_count;
        mbr[0] = points[i].x - epsilon;
        mbr[1] = points[i].y - epsilon;
        mbr[2] = points[i].x + epsilon;
        mbr[3] = points[i].y + epsilon;

        // Search r_tree for neighbor candidates.
        rtree_search(r_tree, mbr, rtree_search_item_handler, (void *)&args_rtree_search);
        
        // Update arguments for get_neighbors.
        args_get_neighbors.focal = points[i];

        // Find neighbors of points[i].
        neighbors_count = get_neighbors((void *)&args_get_neighbors);

        if (neighbors_count < min_points)
        {
            // min_points condition does not hold.
            points[i].label = NOISE;
            continue;
        }

        // Add points[i] to current cluster.
        label += 1;
        points[i].label = label;
        
        // Add neighbors to seed set.
        for (int j = 0; j < neighbors_count; j++)
        {
            cvector_push_back(S, neighbors[j]);
        }

        // Exhaust seed set.
        while (cvector_size(S) > 0)
        {
            // An index to points array.
            int s = S[0];
            cvector_erase(S, 0);

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

            // Add points[s] to current cluster.
            points[s].label = label;

            // Update arguments for get_neighbors.
            args_get_neighbors.focal = points[s];

            neighbors_count = get_neighbors((void *)&args_get_neighbors);

            if (neighbors_count < min_points)
            {
                // min_points condition does not hold.
                continue;
            }

            // Add neighbors to seed set.
            for (int k = 0; k < neighbors_count; k++)
            {
                cvector_push_back(S, neighbors[k]);
            }
        }
    }

    // Free heap.
    cvector_free(S);
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
    int j;
    int label;
    int merged[8];
    int merged_count;
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
    
    // Update arguments for label_points.
    args_label_points.points = points;
    
    // Update arguments for merge_dense_boxes.
    args_merge_dense_boxes.epsilon = epsilon;
    args_merge_dense_boxes.G = G;
    args_merge_dense_boxes.G_x = G_x;
    args_merge_dense_boxes.G_y = G_y;
    args_merge_dense_boxes.merged = merged;
    args_merge_dense_boxes.min_points = min_points;
    args_merge_dense_boxes.points = points;

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

        // Label all points in this dense box with label.
        label += 1;
        args_label_points.indices = G[i].points;
        args_label_points.label = label;
        label_points((void *)&args_label_points);

        // Update arguments for merge_dense_boxes.
        args_merge_dense_boxes.focal = i;
        args_merge_dense_boxes.label = label;
        
        // Merge surrounding dense boxes.
        merged_count = merge_dense_boxes((void *)&args_merge_dense_boxes);

        if (merged_count == 0)
        {
            continue;
        }

        // Queue merged dense boxes.
        for (j = 0; j < merged_count; j++)
        {
            cvector_push_back(queue, merged[j]);
        }

        // Continue merge using queue.
        j = cvector_size(queue);
        while (j > 0)
        {
            // Dequeue.
            int focal = queue[0];
            cvector_erase(queue, 0);
            j -= 1;

            // Update arguments for merge_dense_boxes.
            args_merge_dense_boxes.focal = focal;
            args_merge_dense_boxes.label = label;

            // Merge surrounding dense boxes.
            merged_count = merge_dense_boxes((void *)&args_merge_dense_boxes);

            if (merged_count == 0)
            {
                continue;
            }

            // Queue merged dense boxes.
            for (int k = 0; k < merged_count; k++)
            {
                cvector_push_back(queue, merged[k]);
                j += 1;
            }
        }
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
Merge surrounding dense boxes.

Returns: An int representing the number of dense boxes merged.
*/
int merge_dense_boxes(void *args)
{
    int adjacent_count;
    struct args_label_points args_label_points;
    struct args_should_merge args_should_merge;
    double epsilon;
    // An index to G.
    int focal;
    struct cell *G;
    int G_x;
    int G_y;
    int i;
    bool in_bounds;
    int j;
    int label;
    // Holds indices of dense boxes merged to focal.
    int *merged;
    int merged_count;
    int min_points;
    int non_adjacent[4];
    int non_adjacent_count;
    struct point *points;
    int *queue;
    int x;
    int y;
    
    // Unpack args.
    struct args_merge_dense_boxes *packet = (struct args_merge_dense_boxes *)args;
    epsilon = packet->epsilon;
    focal = packet->focal;
    G = packet->G;
    G_x = packet->G_x;
    G_y = packet->G_y;
    label = packet->label;
    merged = packet->merged;
    min_points = packet->min_points;
    points = packet->points;

    // Get diagonally-adjacent cells.
    int adjacent[8] = {
        (G[focal].x + 1) * G_y + (G[focal].y + 1),
        (G[focal].x + 1) * G_y + (G[focal].y - 1),
        (G[focal].x - 1) * G_y + (G[focal].y - 1),
        (G[focal].x - 1) * G_y + (G[focal].y + 1)
    };
    adjacent_count = 4;

    // Get horizontal/vertical cells.
    non_adjacent_count = 0;
    // Top.
    x = G[focal].x;
    y = G[focal].y + 1;
    in_bounds = (x > -1) && (x < G_x) && (y > -1) && (y < G_y);
    if (in_bounds && cvector_size(G[x * G_y + y].points) < min_points)
    {
        non_adjacent[non_adjacent_count] = x * G_y + (y + 1);
        non_adjacent_count += 1;
    }
    else
    {
        adjacent[adjacent_count] = x * G_y + y;
        adjacent_count += 1;
    }
    // Right.
    x = G[focal].x + 1;
    y = G[focal].y;
    in_bounds = (x > -1) && (x < G_x) && (y > -1) && (y < G_y);
    if (in_bounds && cvector_size(G[x * G_y + y].points) < min_points)
    {
        non_adjacent[non_adjacent_count] = (x + 1) * G_y + y;
        non_adjacent_count += 1;
    }
    else
    {
        adjacent[adjacent_count] = x * G_y + y;
        adjacent_count += 1;
    }
    // Bottom.
    x = G[focal].x;
    y = G[focal].y - 1;
    in_bounds = (x > -1) && (x < G_x) && (y > -1) && (y < G_y);
    if (in_bounds && cvector_size(G[x * G_y + y].points) < min_points)
    {
        non_adjacent[non_adjacent_count] = x * G_y + (y - 1);
        non_adjacent_count += 1;
    }
    else
    {
        adjacent[adjacent_count] = x * G_y + y;
        adjacent_count += 1;
    }
    // Left.
    x = G[focal].x - 1;
    y = G[focal].y;
    in_bounds = (x > -1) && (x < G_x) && (y > -1) && (y < G_y);
    if (in_bounds && cvector_size(G[x * G_y + y].points) < min_points)
    {
        non_adjacent[non_adjacent_count] = (x - 1) * G_y + y;
        non_adjacent_count += 1;
    }
    else
    {
        adjacent[adjacent_count] = x * G_y + y;
        adjacent_count += 1;
    }

    // Update arguments for label_points.
    args_label_points.label = label;
    args_label_points.points = points;

    // Update arguments for should_merge.
    args_should_merge.epsilon = epsilon;
    args_should_merge.points = points;
    
    // Attempt to merge adjacent cells.
    merged_count = 0;
    for (i = 0; i < adjacent_count; i++)
    {
        if (adjacent[i] < 0 || adjacent[i] >= G_x * G_y)
        {
            // Index out of bounds of G.
            continue;
        }

        if (cvector_size(G[adjacent[i]].points) < min_points)
        {
            // Not a dense box.
            continue;
        }

        if (points[G[adjacent[i]].points[0]].label > 0)
        {
            // Already merged.
            continue;
        }

        // Label all points in this dense box.
        args_label_points.indices = G[adjacent[i]].points;
        label_points((void *)&args_label_points);

        // Add this dense box to merged.
        merged[merged_count] = adjacent[i];
        merged_count += 1;
    }

    // Attempt to merge non-adjacent cells.
    for (i = 0; i < non_adjacent_count; i++)
    {
        if (non_adjacent[i] < 0 || non_adjacent[i] >= G_x * G_y)
        {
            // Index out of bounds of G.
            continue;
        }

        if (cvector_size(G[non_adjacent[i]].points) < min_points)
        {
            // Not a dense box.
            continue;
        }

        if (points[G[non_adjacent[i]].points[0]].label > 0)
        {
            // Already merged.
            continue;
        }

        // Update arguments for should_merge.
        args_should_merge.a = G[focal];
        args_should_merge.b = G[non_adjacent[i]];
        
        if (!should_merge((void *)&args_should_merge))
        {
            continue;
        }
        
        // Label all points in this dense box.
        args_label_points.indices = G[non_adjacent[i]].points;
        label_points((void *)&args_label_points);

        // Add this dense box to merged.
        merged[merged_count] = non_adjacent[i];
        merged_count += 1;
    }

    return merged_count;
}

/*
A description.

Returns: A bool.
*/
bool rtree_search_item_handler(const double *mbr, const void *item, void *args)
{
    int *candidates;
    int *candidates_count;
    const struct point *point;
    
    // Unpack args.
    struct args_rtree_search *packet = (struct args_rtree_search *)args;
    candidates = packet->candidates;
    candidates_count = packet->candidates_count;
    point = (const struct point *)item;

    // Add point's index to candidates.
    candidates[*candidates_count] = point->index;
    *candidates_count += 1;

    return true;
}

/*
A description.

Returns: A bool.
*/
bool should_merge(void *args)
{
    struct cell a;
    struct cell b;
    double delta_x;
    double delta_y;
    double distance;
    double epsilon;
    int i;
    int j;
    struct point *points;

    // Unpack args.
    struct args_should_merge *packet = (struct args_should_merge *)args;
    a = packet->a;
    b = packet->b;
    epsilon = packet->epsilon;
    points = packet->points;

    for (i = 0; i < cvector_size(a.points); i++)
    {
        for (j = 0; j < cvector_size(b.points); j++)
        {
            // Calculate distance.
            delta_x = points[a.points[i]].x - points[b.points[j]].x;
            delta_y = points[a.points[i]].y - points[b.points[j]].y;
            distance = sqrt(delta_x * delta_x + delta_y * delta_y);

            if (distance <= epsilon)
            {
                return true;
            }
        }
    }

    return false;
}
