#include <iostream>  // for cout
#include <fstream>  // for ifstream
#include <cfloat>  // for DBL_MAX
#include <cmath>  // for sqrt
#include <cstdlib>  // for rand, srand
#include <ctime>  // for rand seed
#include <cassert>  // for assertions
#include <algorithm>  // for sort
#include <ctime>  // for clock

const bool debug = false;  // debug flag, r
// remember to turn it back to false before submitting your code for automatic validation.

typedef double* point;  // point = array of coordinates (doubles)

void print(point p, int dim) {
  std::cout << p[0];
  for ( int j = 1; j < dim; j++)
    std::cout << " " << p[j];
  std::cout << "\n";
}

void swap (point* P, int i, int j) {  // swap 2 points in data set
  point temp = P[i];
  P[i] = P[j];
  P[j] = temp;
}
void test_swap() {
  int n = 5;  // n points in R^{dim}
  int dim = 2;
  // random input data points (uniformly sampled in unit cube)
  srand (time(NULL));
  point P[n];
  for (int i=0; i<n; i++) {
    P[i] = new double [dim];
    for (int j=0; j<dim; j++)
      P[i][j] = (double)rand() / RAND_MAX;
  }
  for (int i=0; i<n; i++)
    print(P[i], dim);
  swap (P, 1, 3);  // swap 2 points in data set
  std::cout << "After swap" << std::endl;
  for (int i=0; i<n; i++)
    print(P[i], dim);
}


double dist (point p, point q, int dim) {
  // complete at ex. 1.1
  double res = 0;
  for(int i=0;i<dim;i++){
    res += (p[i] - q[i]) * (p[i] - q[i]);
  }
  return sqrt(res);
}

int linear_scan (point q, int dim, point* P, int n) {
  // complete at ex. 1.2
  int res = 0;
  double dmin = DBL_MAX;
  while(n --> 0) {
    double dis = dist(q,P[n],dim);
    if(dis < dmin){
      dmin = dis;
      res = n;
    }
  }
  return res;
}

double computeMedian(point* P, int start, int end, int c) {
  // complete at ex. 2.1
  
  double d[end-start];
  for(int i = 0; i<(end-start);i++){
    d[i] = P[start+i][c];
  }
  std::sort(d,d+(end-start));

  int mid = (end-start)/2;
  return d[mid];
}


int partition(point* P, int start, int end, int c, int dim) {
  double m = computeMedian (P, start, end, c);
  // complete at ex. 2.2
  // for(int i = start; i < end; i++){
  //   int tmp = i;
  //   while(P[tmp][c]<) swap(P,)
  // }
  
  c = c%dim;
  int p = (start + end) / 2;
  end--;
  while (end-start>0)
  {
    swap(P,start,end);
    while (P[start][c] < m) {start++;}
    while (P[end][c] > m) {end--;}
  }
  // if(P[p][c]!=m) {
  //   if(P[p+1][c]== m) swap(P,p,p+1);;
  //   if(P[p-1][c] ==m) swap(P,p,p-1);
  // }
  assert(P[p][c]==m);
  //if(P[p+1][c]<P[p-1][c]) swap(P,p+1,p-1);
  return p;
// for (int i = start; i < p; i++)
// {
//   if (P[i][c] = m)
//   {
//     swap[P, i, m];
//     i--;
//   }
//   else if (P[i][c] > m)
//   {
//     for (int j = p; j < end; j++)
//     {
//       if (P[j][c] < m)
//         swap(P, i, j);
//     }
//   }

}

typedef struct node {  // node for the kd-tree
  // complete at ex. 2.3
  int c; //coordinate for split
  double m;// split value
  int p;// index of data point srored at the node
  node *left, *right; //children

} node;

node* create_node (int _p) {  // creates a leaf node
  // complete at ex. 2.3
  node* n = new node;
  n->p = _p;
  n->left = NULL;
  n->right = NULL;
  //return n;
}

node* create_node (int _c, double _m, int _p, node* _left, node* _right) {  // creates an internal node
  // complete at ex. 2.3
  node* n = new node;
  n->c = _c;
  n->m = _m;
  n->p = _p;
  n->left = _left; n->right = _right;
  //return n;
}

node* build (point* P, int start, int end, int c, int dim) {
  
  // builds tree for sub-cloud P[start -> end-1]
  
  assert (end-start >= 0);
  if (debug)
    std::cout << "start=" << start << ", end=" << end << ", c=" << c << std::endl;
  
  if (end-start == 0)  // no data point left to process
    return NULL;  
  else if (end-start == 1)  // leaf node
    return create_node (start);
  
  else {  // internal node
    if (debug) {
      std::cout << "array:\n";
      for (int i=start; i<end; i++)
        print(P[i],dim);
	// std::cout << P[i] << ((i==end-1)?"\n":" ");
    }
    
    // compute partition
    // rearrange subarray (less-than-median first, more-than-median last)
    
    int p = partition (P, start, end, c, dim);
    double m = P[p][c];
    
    // prepare for recursive calls
    
    int cc = (c+1)%dim;  // next coordinate
    return create_node (c, m, p,
		     build (P, start, p, cc, dim),
		     build (P, p+1, end, cc, dim));
  }
}
  
void dsearch (node* n, point q, int dim, point* P, double& res, int& nnp) {
  // complete at ex. 2.4
  if (n != NULL) {
    double dis = dist(q,P[n->p],dim);
    if(res>dis) {
      res = dis;
      nnp = n->p;
    }
    if (n->left != NULL || n->right != NULL)  // internal node
      dsearch ( (q[n->c] <= n->m)?n->left:n->right, q, dim, P, res,nnp);
  }
}

void bsearch (node* n, point q, int dim, point* P, double& res, int& nnp) {
  // complete at ex. 2.5
    if (n != NULL) {
    double dis = dist(q,P[n->p],dim);
    if(res>dis) {
      res = dis;
      nnp = n->p;
    }
    if (n->left != NULL || n->right != NULL) {  // internal node
      if (q[n->c] - res <= n->m)  // ball intersects left side of H
	      bsearch ( n->left, q, dim, P, res,nnp);
      if (q[n->c] + res > n->m)  // ball intersects right side of H
	      bsearch ( n->right, q, dim, P, res,nnp);
    }
  }
}

void computeM_test(){
  std::cout <<"----------------------------begin computeMedium Test------------------------------"<<std::endl;
  int n=10;
  int dim = 2;
  srand (time(NULL));
  point P[n];
  for (int i=0; i<n; i++) {
    P[i] = new double [dim];
    for (int j=0; j<dim; j++)
      P[i][j] = (double)rand() / RAND_MAX;
  }
  for (int i=0; i<n; i++)
    print(P[i], dim);
  double m_3_7,m_4_8,m_2_9;
  m_2_9 = computeMedian(P,2,9,0);
  m_3_7 = computeMedian(P,3,7,1);
  m_4_8 = computeMedian(P,4,8,0);
  std::cout << "midium de [2,9[, en dim=0" << std::endl;
  std::cout << m_2_9 << std::endl;
  std::cout << "midium de [3,7[, en dim=1" << std::endl;
  std::cout << m_3_7 << std::endl;
  std::cout << "midium de [4,8[, en dim=0" << std::endl;
  std::cout << m_4_8 << std::endl;

  std::cout <<"----------------------------begin partition Test------------------------------"<<std::endl;

  
}

int main () {
  // test_swap();
   computeM_test();
 
  const int dim = 4;  // dimension (hard-coded)
  int n = 20000;  // upper bound on number of data points in R^{dim}
  int nt = 1000;  // nt query points

  point P[n];
  char names[n][255];
  
  srand (time(NULL));
  
  for (int i=0; i<n; i++)
    P[i] = new double[dim];
  
  std::ifstream is("iris2.data");
  
  assert(is.is_open());
  for (n=0; is >> P[n][0]; n++) {
    for (int i=1; i<dim; i++)
      is >> P[n][i];
    is >> names[n];
  }
  
  std::cout << n << " observations" << std::endl;


/********* section 1 ***************/
  
  // Some inputs / outputs
  // while (true){
    
  //   std::cout << "\nEnter your own measurements (4 space-separated real numbers in [0,10]):\n(type CTRL-D to exit): " << std::endl;
    
  //   if (std::cin.peek() == EOF) break;
    
  //   point query = new double [dim];
  //   std::cin >> query[0] >> query[1] >> query[2] >> query[3];
    
  //   double distb = DBL_MAX;
  //   int nn;
  //   nn = linear_scan (query, dim, P, n);
  //   // bsearch(tree, query, dim, P, distb, nn);
  //   std::cout << "\n  -> Your iris is of type " << names[nn] << std::endl;
  // }

  
/******* section 2 ********************************************/    


  // build kd-tree (warning: rearranges points in data set)
  std::cout << "building kd-tree..." << std::flush;
  node* tree = build (P, 0, n, 0, dim);
  std::cout << " done" << std::endl;
  

  double* d = new double[nt];  // recorded distances
  point* q = new point[nt];  // query points
  int* nnp = new int[nt];  // nearest neighbors

  // random query points
  for (int i=0; i<nt; i++) {
    q[i] = new double [dim];
    for (int j=0; j<dim; j++)
      q[i][j] = 10*(double)rand() / RAND_MAX;
  }


  // benchmarking, start and end time for chronos
  std::clock_t start, end;


  // linear scan
  std::cout << "Benchmarking linear scan..." << std::flush;
  start = std::clock();  
  for (int i=0; i<nt; i++) {
    int nn = linear_scan (q[i], dim, P, n);
    d[i] = dist(q[i], P[nn], dim);
    nnp[i] = nn;
    // std::cout << "nearest neighbor (linear scan) = P[" << nn
    // 	      << "] / distance = " << d[i] << std::endl; 
  }
  end = std::clock();  
  std::cout << " done (avg time = "
	    << (end - start)/nt
	    << " us)" << std::endl;


  // defeatist search
  std::cout << "Benchmarking defeatist search..." << std::flush;
  int nf = 0;
  start = std::clock();  
  for (int i=0; i<nt; i++) {
    double distd = DBL_MAX;
    int nn;
    dsearch(tree, q[i], dim, P, distd, nn);
    assert (distd >= d[i]);
    if (nn != nnp[i]) nf++;
    // std::cout << "nearest neighbor (defeatist search): " << nnp
    // 	      << " | distance = " << distd << std::endl; 
  }
  end = std::clock();  
  std::cout << " done (avg time = "
	    << (end - start)/nt
	    << " us, failure rate = " << (nf*100.0/nt)
	    << "\%)" << std::endl;

  // backtracking search
  std::cout << "Benchmarking backtracking search..." << std::flush;
  nf = 0;
  start = std::clock();  
  for (int i=0; i<nt; i++) {
    // std::cerr << i << " ";
    double distb = DBL_MAX;
    int nn;
    bsearch(tree, q[i], dim, P, distb, nn);
    // assert (nn == nnp[i]); 
    if (distb != d[i]) nf++;
  }
  end = std::clock();  
  std::cout << " done (avg time = "
	    << (end - start)/nt
	    << " us, failure rate = " << (nf*100.0/nt)
	    << "\%)" << std::endl;



  return 0;
}
