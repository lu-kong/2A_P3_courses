#pragma once // ensures that this header file is only included once

class point  // this is just a declaration of the class 
             // implementation is in a separate file: point.cpp
{
private:
    // @OFF
    // static int d; // Must be moved here to make private 
    static int d; 
    // @ON
  
public:
    // @OFF
    // Move above to make private
    // @ON 
    double *coords;
    int label;

    static bool set_dim(int _d);
    static int get_dim();
    
    point();
    ~point();

    void print();
    double dist(point &q);
};
