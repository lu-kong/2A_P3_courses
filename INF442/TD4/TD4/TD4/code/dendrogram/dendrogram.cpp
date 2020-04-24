#include "dendrogram.hpp" // This is the header for the class implemented here

#include "cloud.hpp"      // Variable c of type cloud& is used in the constructor and several other methods
#include "edge.hpp"       // Used in several methods, e.g. dendrogram::merge(edge *e)
#include "graph.hpp"      // Variable g of type graph* is used in the constructor and several other methods

#include <iostream> // This provides the input/output functionality
#include <cassert>  // This provides the assert() method

using namespace std;

dendrogram::dendrogram(cloud &_c)
{
    c = &_c;
    g = new graph(_c);

    int n = c->get_n();

    parent = new int[n];
    rank = new int[n];

    left = new int[n];
    down = new int[n];

    height = new double[n];

    clusters = new int[n];

    for (int i = 0; i < n; i++)
    {
        parent[i] = -1;
        rank[i] = 0;

        left[i] = -1;
        down[i] = -1;

        height[i] = -1;

        clusters[i] = -1;
    }

    cut_height = -1;
    nb_clusters = 0;
}

dendrogram::~dendrogram()
{
    delete g;

    delete[] parent;
    delete[] left;
    delete[] down;
    delete[] rank;
    delete[] height;
    delete[] clusters;

    if (sign_heights != nullptr)
        delete[] sign_heights;
}

int dendrogram::find(int i)
{
    assert(0 <= i && i < c->get_n());
    int result = 0;

    // TODO: Exercise 2.2
    return result;
}

void dendrogram::merge(edge *e)
{
    // TODO: Exercise 2.2
}

void dendrogram::build()
{
    g->start_iteration();

    // TODO: Exercise 2.2
}

double dendrogram::get_dendro_height()
{
    return height[down[find(0)]];
}

void dendrogram::set_clusters(int i, double h)
{
    assert(0 <= i && i < c->get_n());
    point &p = c->get_point(i);

    // TODO: Exercise 3.1
}

void dendrogram::set_clusters(double h)
{
    if (cut_height != h) {
        cut_height = h;
        set_clusters(find(0), h);
    }
}

int dendrogram::count_clusters(int i) {
    int count = 0;
    
    // TODO: Exercise 3.2

    return count;
}

int dendrogram::count_clusters() {
    if (nb_clusters == 0)
        nb_clusters = count_clusters(find(0));
    return nb_clusters;
}

void dendrogram::clear_clusters()
{
    for (int i = 0; i < c->get_n(); i++)
    {
        c->get_point(i).label = -1;
        clusters[i] = -1;
    }
    nb_clusters = 0;
    cut_height = -1;
}

double dendrogram::get_cluster_height(int cluster)
{
    assert(0 <= cluster && cluster < c->get_n());

    // TODO: Exercise 3.3

    return 0; // Singleton cluster at a node
}

/******** Significant heights ********/

int count_non_zero(double *unfiltered, int countu)
{
    int countf = 0;
    for (int i = 0; i < countu; i++)
        if (unfiltered[i] > 0)
            countf++;

    return countf;
}

double *filter_double_array(double *unfiltered, int countu, int countf)
{
    double *filtered = new double[countf];
    int f = 0;
    for (int i = 0; i < countu; i++)
        if (unfiltered[i] > 0)
            filtered[f++] = unfiltered[i];
    return filtered;
}

void dendrogram::find_heights(double eps)
{
    double h = get_dendro_height();
    int slots = 1 / eps + 1;
    double *buckets = new double[slots];
    for (int i = 0; i < slots; i++)
        buckets[i] = 0;

    for (int i = 0; i < c->get_n(); i++)
        if (height[i] > buckets[(int)(height[i] / eps / h)])
            buckets[(int)(height[i] / eps / h)] = height[i];

    count_sign_heights = count_non_zero(buckets, slots);

    if (sign_heights != nullptr)
        delete[] sign_heights;

    sign_heights = filter_double_array(buckets, slots, count_sign_heights);
    delete[] buckets;
}

/***** GETTERS *****/

cloud *dendrogram::get_cloud()
{
    return c;
}

int dendrogram::get_parent(int i)
{
    assert(0 <= i && i < c->get_n());
    return parent[i];
}

int dendrogram::get_down(int i)
{
    assert(0 <= i && i < c->get_n());
    return down[i];
}

int dendrogram::get_left(int i)
{
    assert(0 <= i && i < c->get_n());
    return left[i];
}

int dendrogram::get_rank(int i)
{
    assert(0 <= i && i < c->get_n());
    return rank[i];
}

double dendrogram::get_height(int i)
{
    assert(0 <= i && i < c->get_n());
    return height[i];
}

int dendrogram::get_count_sign_heights()
{
    return count_sign_heights;
}

double dendrogram::get_sign_height(int i)
{
    assert(0 <= i && i < count_sign_heights);
    return sign_heights[i];
}

int dendrogram::get_count_clusters()
{
    return nb_clusters;
}

double dendrogram::get_cut_height()
{
    return cut_height;
}

/****** Methods for testing ******/

void dendrogram::print_node(int i)
{
    cout << "Node " << i
         << "(parent =  " << parent[i]
         << ", down = " << down[i]
         << ", left = " << left[i]
         << ", rank = " << rank[i]
         << ", height = " << height[i]
         << ", clusters = " << clusters[i]
         << ")";
}

void dendrogram::print_dendrogram()
{
    cout << "node\tpoint\tparent\trank\tleft\tdown\theight\tclusters" << endl;
    for (int i = 0; i < c->get_n(); i++)
        cout << i << "\t"
             << c->get_point(i).etiquette << "\t"
             << parent[i] << "\t"
             << rank[i] << "\t"
             << left[i] << "\t"
             << down[i] << "\t"
             << height[i] << "\t"
             << clusters[i] << "\t"
             << endl;
}

void dendrogram::print_clusters()
{
    // Go through all nodes
    for (int i = 0; i < c->get_n(); i++)
        // For each encountered cluster representative
        if (c->get_point(i).label == i && clusters[i] != -1)
        {
            // Print the corresponding cluster
            cout << "Cluster \"" << c->get_point(i).etiquette
                 << "\" (node: " << i << ";"
                 << " height: " << get_cluster_height(i)
                 << ")" << endl;

            // Print all points in the cluster
            int point = i;
            do
            {
                cout << point;
                if (c != nullptr)
                {
                    cout << " (";
                    c->get_point(point).print();
                    cout << ")" << endl;
                }
                point = clusters[point];
            } while (point != -1);
        }
}

void dendrogram::iterate_sign_heights()
{
    static int i = 0;
    double cut = get_sign_height(i);
    cout << "Setting clusters at height " << cut
         << "..." << endl;
    set_clusters(cut);
    count_clusters();
    cout << "\t" << get_count_clusters() << " non-singleton cluster"
         << (get_count_clusters() > 1 ? "s" : "")
         << " below " << cut << endl;
    i = (i + 1) % get_count_sign_heights();
}

void dendrogram::set_parent(int *_parent) {
    for (int i = 0; i < c->get_n(); i++)
        parent[i] = _parent[i];
}
