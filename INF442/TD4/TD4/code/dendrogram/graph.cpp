#include "graph.hpp" // This is the header for the class implemented here

#include "cloud.hpp" // Used in the constructor
#include "edge.hpp"  // Used in almost all methods

#include <algorithm> // This provides the sort() method for the constructor

/* graph -- method implementations */

graph::graph(cloud &_c)
{
    c = &_c;
    edge::set_cloud(c);
    iterator_pos = 0;
    size = c->get_n() * (c->get_n() - 1) / 2;
    // TODO: Exercise 2.1
    int n = c->get_n();
    edges = new edge*[size];
    for (int i = 0; i < c->get_n() - 1; i++)
    { //i<n-1
        //*edges[i] = new edge [10];
        for (int j = i + 1; j < c->get_n(); j++)
        {
            edges[iterator_pos++] = new edge(j,i);
        }
    }
    std::sort(edges,edges+size,edge::compare);
    
}

graph::~graph()
{
    // TODO: Exercise 2.1
    // graph::c = NULL;
    // graph::size = 0;
    // graph::iterator_pos = 0;
    for(int i = 0; i<size;i++)
    {
        delete edges[i];
    }
    delete[] edges;
}

long graph::get_size()
{
    // TODO: Exercise 2.1
    return size;
}

void graph::start_iteration()
{
    iterator_pos = 0;
}

edge *graph::get_next()
{
    // TODO: Exercise 2.1
    if(iterator_pos<size) return edges[iterator_pos++];
    else return NULL;
}
