#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct graph
{
    int* edges;
    int node_counter;
};

struct graph graph_a;
struct graph graph_b;

void init_graph(int graph_id, int node_count)
{
    //printf("\ninit_graph called");
    if (graph_id == 0)
    {
        graph_a.edges = (int*) malloc(node_count * node_count * sizeof(int));
        memset(graph_a.edges, 0, node_count * node_count * sizeof(int));
        graph_a.node_counter = node_count;
    }
    else if (graph_id == 1)
    {
        graph_b.edges = (int*) malloc(node_count * node_count * sizeof(int));
        memset(graph_b.edges, 0, node_count * node_count * sizeof(int));
        graph_b.node_counter = node_count;
    }
    else
    {
        printf("error, no such graph\n");
    }
}

void add_edge(int graph_id, int from_node, int to_node, int edge_weight)
{
    //printf("\nadd_edge called");
    if (graph_id == 0)
    {
        graph_a.edges[from_node * graph_a.node_counter + to_node] = edge_weight;
        graph_a.edges[to_node * graph_a.node_counter + from_node] = edge_weight;
    }
    else if (graph_id == 1)
    {
        graph_b.edges[from_node * graph_b.node_counter + to_node] = edge_weight;
        graph_b.edges[to_node * graph_b.node_counter + from_node] = edge_weight;
    }
    else
    {
        printf("error, no such graph\n");
    }
}

int has_edge(int graph_id, int from_node, int to_node)
{
    //printf("\nhas_edge called");
    if (graph_id == 0)
    {
        return graph_a.edges[(from_node * graph_a.node_counter) + to_node];
    }
    else if (graph_id == 1)
    {
        return graph_b.edges[(from_node * graph_b.node_counter) + to_node];
    }
    else
    {
        printf("error, no such graph\n");
        return -1;
    }
}

int _len_(int graph_id)
{
    //printf("\n_len_ called");
    if (graph_id == 0)
    {
        int result = 0;
        int count = 0;
        while (count < graph_a.node_counter * graph_a.node_counter)
        {
            if (graph_a.edges[count] > 0)
                result += 1;
            count += 1;
        }
        return result;
    }
    else if (graph_id == 1)
    {
        int result = 0;
        int count = 0;
        while (count < graph_b.node_counter * graph_b.node_counter)
        {
            if (graph_b.edges[count] > 0)
                result += 1;
            count += 1;
        }
        return result;
    }
    else
    {
        printf("error, no such graph\n");
        return -1;
    }
}
