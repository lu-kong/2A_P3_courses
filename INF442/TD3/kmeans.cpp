#include <iostream>
#include <cassert>
#include <cmath>	// for sqrt, fabs
#include <cfloat>	// for DBL_MAX
#include <cstdlib>	// for rand, srand
#include <ctime>	// for rand seed
#include <fstream>
#include <cstdio>	// for EOF
#include <string>
#include <algorithm>	// for count
#include <vector>

#include <gtkmm/application.h>
#include <gtkmm/window.h>
#include <gtkmm/drawingarea.h>

struct point
{
        static int d;
        double *coords;
        int label;
		point () {coords = new double[d];label = 0;
		for(int j =0;j<d;j++){
			coords[j] = 0.0;
		}}
		~point () {delete[] coords;}

		void print() {
			for(int k=0;k<d-1;k++){
				std::cout << coords[k] << '\t';
			}
			std::cout << coords[d-1] << std::endl;
		}

		double dist(point &q){
			double res=0;
			int ind = d;
			while (ind-->0)
			{
				res += (coords[ind]-q.coords[ind])*(coords[ind]-q.coords[ind]);
			}
			return sqrt(res);
			
		}
};

int point::d;

class cloud
{
	private:
	int d;
	int n;
	int k;

	// maximum possible number of points
	int nmax;

	point *points;
	point *centers;

	public:
	cloud(int _d, int _nmax, int _k)
	{
		d = _d;
		point::d = _d;
		n = 0;
		k = _k;

		nmax = _nmax;

		points = new point[nmax];
		centers = new point[k];
	}

	~cloud()
	{
		delete[] centers;
		delete[] points;
	}

	void add_point(point &p, int label)
	{
		assert(n < nmax);

		for(int m = 0; m < d; m++)
		{
			points[n].coords[m] = p.coords[m];
		}

		points[n].label = label;

		n++;
	}

	int get_d()
	{
		return d;
	}

	int get_n()
	{
		return n;
	}

	int get_k()
	{
		return k;
	}

	point& get_point(int i)
	{
		return points[i];
	}

	point& get_center(int j)
	{
		return centers[j];
	}
	void summury(){
		std::cout<<"cloud_summury-------------\n"<<"number of points : "<<n;
		std::cout<<"\t number of clusers : "<<k<<"\t dim : "<< d<<std::endl;
		std::cout<<"centers : \n";
		int m = std::min(10,k);
		int mm = std::min(10,n);
		for(int j =0;j<m;j++){
			std::cout<<j<<"\t\t";
			centers[j].print();
		}
		std::cout<<"points : \n";
		for(int i=0;i<mm;i++){
			std::cout<<"of label : "<<points[i].label<<"\t";
			points[i].print();
		}
		std::cout<<"\n";
	}
	void set_center(point &p, int j)
	{
		for(int m = 0; m < d; m++)
			centers[j].coords[m] = p.coords[m];
	}

	double intracluster_variance()
	{
		double res=0.0;
		for(int i=0;i<n;i++){
			double dis = points[i].dist(centers[(points[i].label)]);
			res += dis*dis;
		}
		return res/n;
	}

	void set_voronoi_labels()
	{
		for(int ip=0; ip<n; ip++){
			double dis = DBL_MAX;
			int clu;
			for(int i=0;i<k;i++){
				if(points[ip].dist(centers[i])<dis){clu = i;dis = points[ip].dist(centers[i]);}
			}
			points[ip].label = clu;
		}
	}

	void set_centroid_centers(){
		double** res = new double* [k];
		int* num = new int[k]();
		for(int j=0;j<k;j++){ res[j] = new double[d]();}//initialize the matrix to store values
								//res[i][d] = 1} d is the label bool for detect the changement
		//assign to each center
		
		for(int i=0;i<n;i++){
			int lab = points[i].label;
			// if(res[lab][d]!=1){res[lab][d]++;}
			num[lab]++;
			for(int m=0;m<d;m++){res[lab][m]+=points[i].coords[m];}
		}

		//final set
		for(int l=0;l<k;l++){
			if(num[l]==0){break;}
			for(int dim=0;dim<d;dim++){
				res[l][dim]/=num[l]; //divided by number of points
			}
			centers[l].coords = res[l];
		}
	}
	
	void kmeans()
	{
		int times = 1000;
		double icv_old = DBL_MAX;
		double icv_new = DBL_MAX/2;
		const double eps = 0.000001;
		set_centroid_centers();
		while(times-->0 && icv_old-icv_new>eps){
			set_voronoi_labels();
			set_centroid_centers();
			icv_old=icv_new;
			icv_new=intracluster_variance();
		}
		std::cout<<1000-times<<std::endl;
	}

	void init_forgy()
	{
	}

	void init_plusplus()
	{
	}

	void init_random_partition()
	{
	}
};

double dist(point& p, point& q, int d) {
	double res = 0.0;
	while(d-->0) {
		res += (p.coords[d]-q.coords[d])*(p.coords[d]-q.coords[d]);
	}
	return sqrt(res);
}
// test functions
void test_point(){

	point::d = 5; //set the static value of dimension
	point p,origin; // call point p where all the coords are 0
	std::cout<< "Testing class point..."<<std::endl;
	p.coords = new double[5] {1.0, 2.0, 3.0, 4.0, 5.0};

	std::cout<< "coordinates of point p..."<<std::endl;
	p.print();
	std::cout<< "coordinates of point origin..."<<std::endl;
	origin.print();
	std::cout<< "dist() test..."<<std::endl;
	std::cout<< "p.dist(origin) = "<<p.dist(origin)<< std::endl;
	std::cout<< "dist(p,origin,d) = "<<dist(p,origin,5)<< std::endl;
	assert(p.dist(origin)==dist(p,origin,5));
	std::cout<< "Testing class point... done."<<std::endl;
}
void test_set_voronoi_labels(){
	std::cout << "Testing function set_voronoi_labels()...\n";

	const double eps = 0.0001;
	// dimension used for tests
	point::d = 1;
	// temporary container
	point p,c,b,a;

	// test case 4
	const double dist_twoclusters = 6.8125;
	cloud twoclusters(1, 4, 2);
	p.coords[0] = -1.0;
	twoclusters.add_point(p, 0);
	p.coords[0] = 0.5;
	twoclusters.add_point(p, 0);
	p.coords[0] = -0.5;
	twoclusters.set_center(p, 0);
	p.coords[0] = -2.0;
	twoclusters.add_point(p, 1);
	p.coords[0] = 2.0;
	twoclusters.add_point(p, 1);
	p.coords[0] = -3.0;
	twoclusters.set_center(p, 1);
	std::cout<<"----------before set labels: -----------\n";
	twoclusters.summury();
	twoclusters.set_voronoi_labels();
	std::cout<<"----------after set labels: -----------\n";
	twoclusters.summury();
	std::cout << "\t[OK]" << std::endl;
	
}

void test_set_centroid_centers(){
	std::cout << "Testing function set_voronoi_labels()...\n";

	const double eps = 0.0001;
	// dimension used for tests
	point::d = 1;
	// temporary container
	point p,c,b,a;

	// test case 4
	const double dist_twoclusters = 6.8125;
	cloud twoclusters(1, 4, 2);
	p.coords[0] = -1.0;
	twoclusters.add_point(p, 0);
	p.coords[0] = 0.5;
	twoclusters.add_point(p, 0);
	p.coords[0] = -0.5;
	twoclusters.set_center(p, 0);
	p.coords[0] = -2.0;
	twoclusters.add_point(p, 1);
	p.coords[0] = 2.0;
	twoclusters.add_point(p, 1);
	p.coords[0] = -3.0;
	twoclusters.set_center(p, 1);
	std::cout<<"----------original: -----------\n";
	twoclusters.summury();
	twoclusters.set_voronoi_labels();
	std::cout<<"----------after set labels: -----------\n";
	twoclusters.summury();
	twoclusters.set_centroid_centers();
	std::cout<<"----------after set centers: -----------\n";
	twoclusters.summury();
	std::cout << "\t[OK]" << std::endl;
}

void test_intracluster_variance()
{
	// tolerance for comparison of doubles
	const double eps = 0.0001;

	// dimension used for tests
	point::d = 1;

	// temporary container
	point p;

	std::cout << "Testing function intracluster_variance()...";

	// test case 1
	const double dist_onepoint_zerodist = 0.0;
	cloud onepoint_zerodist(1, 1, 1);
	p.coords[0] = 0.0;
	onepoint_zerodist.add_point(p, 0);
	assert(std::fabs(onepoint_zerodist.intracluster_variance() - dist_onepoint_zerodist) < eps);

	// test case 2
	const double dist_onepoint_posdist = 0.25;
	cloud onepoint_posdist(1, 1, 1);
	p.coords[0] = 0.5;
	onepoint_posdist.add_point(p, 0);
	assert(std::fabs(onepoint_posdist.intracluster_variance() - dist_onepoint_posdist) < eps);

	// test case 3
	const double dist_twopoints = 0.625;
	cloud twopoints(1, 2, 1);
	p.coords[0] = -1.0;
	twopoints.add_point(p, 0);
	p.coords[0] = 0.5;
	twopoints.add_point(p, 0);
	p.coords[0] = -0.5;
	twopoints.set_center(p, 0);
	assert(std::fabs(twopoints.intracluster_variance() - dist_twopoints) < eps);

	// test case 4
	const double dist_twoclusters = 6.8125;
	cloud twoclusters(1, 4, 2);
	p.coords[0] = -1.0;
	twoclusters.add_point(p, 0);
	p.coords[0] = 0.5;
	twoclusters.add_point(p, 0);
	p.coords[0] = -0.5;
	twoclusters.set_center(p, 0);
	p.coords[0] = -2.0;
	twoclusters.add_point(p, 1);
	p.coords[0] = 2.0;
	twoclusters.add_point(p, 1);
	p.coords[0] = -3.0;
	twoclusters.set_center(p, 1);
	assert(std::fabs(twoclusters.intracluster_variance() - dist_twoclusters) < eps);

	std::cout << "\t[OK]" << std::endl;
}

void test_kmeans()
{
	// TODO
}

void test_init_forgy()
{
	// number of random experiments
	const int K = 10000;
	// tolerance in probability
	const double delta = 0.0625;

	// dimenstion used for tests
	point::d = 1;

	// temporary container
	point p;

	std::cout << "Testing function init_forgy()...";

	const double prob_threepoints = 0.3333;
	cloud threepoints(1, 3, 1);
	p.coords[0] = 0.0;
	threepoints.add_point(p, 0);
	p.coords[0] = 1.0;
	threepoints.add_point(p, 0);
	p.coords[0] = 2.0;
	threepoints.add_point(p, 0);
	int cnt = 0;
	for(int k = 0; k < K; k++)
	{
		threepoints.init_forgy();
		if(threepoints.get_center(0).coords[0] == 1.0)
			cnt++;
	}
	assert(std::fabs(cnt/(double)K - prob_threepoints) < delta);

	std::cout << "\t\t[OK]" << std::endl;
}

void test_init_plusplus()
{
	// number of random experiments
	const int K = 10000;
	// tolerance in probability
	const double delta = 0.0625;

	// dimenstion used for tests
	point::d = 1;

	// temporary container
	point p;

	std::cout << "Testing function init_plusplus()...";

	// test case 1
	const double prob_threepoints = 0.3333;
	cloud threepoints(1, 3, 1);
	p.coords[0] = 0.0;
	threepoints.add_point(p, 0);
	p.coords[0] = 1.0;
	threepoints.add_point(p, 0);
	p.coords[0] = 2.0;
	threepoints.add_point(p, 0);
	int cnt = 0;
	for(int k = 0; k < K; k++)
	{
		threepoints.init_plusplus();
		if(threepoints.get_center(0).coords[0] == 1.0)
			cnt++;
	}
	assert(std::fabs(cnt/(double)K - prob_threepoints) < delta);

	// test case 2
	const double prob_twoclusters = 0.125;
	cloud twoclusters(1, 4, 2);
	p.coords[0] = 0.0;
	twoclusters.add_point(p, 0);
	p.coords[0] = 0.0;
	twoclusters.add_point(p, 0);
	p.coords[0] = 1.0;
	twoclusters.add_point(p, 0);
	p.coords[0] = 2.0;
	twoclusters.add_point(p, 0);
	cnt = 0;
	for(int k = 0; k < K; k++)
	{
		twoclusters.init_plusplus();
		if(twoclusters.get_center(1).coords[0] == 1.0)
			cnt++;
	}
	assert(std::fabs(cnt/(double)K - prob_twoclusters) < delta);

	std::cout << "\t\t[OK]" << std::endl;
}

void test_init_random_partition()
{
	// number of random experiments
	const int K = 10000;
	// tolerance in probability
	const double delta = 0.0625;

	// dimenstion used for tests
	point::d = 1;

	// temporary container
	point p;

	std::cout << "Testing function init_random_parition()...";

	const double prob_threepoints = 0.3333;
	cloud threepoints(1, 3, 3);
	p.coords[0] = 0.0;
	threepoints.add_point(p, 0);
	p.coords[0] = 1.0;
	threepoints.add_point(p, 0);
	p.coords[0] = 2.0;
	threepoints.add_point(p, 0);
	int cnt = 0;
	for(int k = 0; k < K; k++)
	{
		threepoints.init_random_partition();
		if(threepoints.get_point(2).label == 1)
			cnt++;
	}
	assert(std::fabs(cnt/(double)K - prob_threepoints) < delta);

	std::cout << "\t[OK]" << std::endl;
}


// for graphical interface
class MyArea : public Gtk::DrawingArea
{
	private:
	cloud *c;
	int x_column;
	int y_column;

	public:
	MyArea(cloud *_c, int _x_column, int _y_column)
	{
		c = _c;
		x_column = _x_column;
		y_column = _y_column;
	}

	virtual ~MyArea() {}

protected:
	//Override default signal handler:
	bool on_draw(const Cairo::RefPtr<Cairo::Context> &cr) override;
};

/**
 * Counts the number of tab-separated columns in the given line.
 */
int nb_columns(const std::string &line)
{
	return std::count(line.begin(), line.end(), '\t') + 1;
}

int main(int argc, char **argv)
{
	bool run_tests = false;

	if(argc < 5 || argc > 6)
	{
		std::cerr << "Usage: " << argv[0] << " csv nb_clusters x_column y_column [ test ]" << std::endl;
		std::exit(1);
	}
	std::string csv_filename = argv[1];
	int nb_clusters = std::stoi(argv[2]);
	int x_column = std::stoi(argv[3]);
	int y_column = std::stoi(argv[4]);

	if(argc >= 6)
		run_tests = true;

	srand(time(NULL));

	if(run_tests)
	{
		// test_point();
		//test_set_voronoi_labels();
		test_set_centroid_centers();
		test_intracluster_variance();
		//test_kmeans();
		// test_init_forgy();
		// test_init_plusplus();
		// test_init_random_partition();
	}

	// open data file
	std::ifstream is(csv_filename);
	assert(is.is_open());

	// read header line
	std::string header_line;
	std::getline(is, header_line);
	std::vector<std::string> names;

	const int d = nb_columns(header_line) - 1;
	const int nmax = 150;
	const int k = nb_clusters;

	// construct point cloud
	cloud c(d, nmax, k);

	// point to read into
	point p;

	// labels to cycle through
	int label = 0;

	// while not at end of file
	while(is.peek() != EOF)
	{
		// read new points
		for(int m = 0; m < d; m++)
		{
			is >> p.coords[m];
		}

		c.add_point(p, label);

		label = (label + 1) % k;

		// read ground-truth labels
		// unused in normal operation
		std::string next_name;
		is >> next_name;
		names.push_back(next_name);

		// consume \n
		is.get();
	}

	// execute k-means algorithm
	std::cout << "Intracluster variance before k-means: " << c.intracluster_variance() << std::endl;
	c.kmeans();
	std::cout << "Intracluster variance after k-means: " << c.intracluster_variance() << std::endl;

	std::cout << "Saving clustering into \"output.csv\"... "; 
	std::ofstream os("output.csv");
	assert(os.is_open());
	os << header_line << '\n';
	for(int i = 0; i < c.get_n(); ++i)
	{
		for(int j = 0; j < c.get_d(); ++j)
		{
			os << c.get_point(i).coords[j] << '\t';
		}
		os << names[i] << "_Label_" << c.get_point(i).label;
		if(i != c.get_n()-1)
			os << '\n';
	}
	std::cout << "done" << std::endl;

	// launch graphical interface
	int gtk_argc = 1;
	Glib::RefPtr<Gtk::Application> app = Gtk::Application::create(gtk_argc, argv, "inf442.td3");

	Gtk::Window win;
	win.set_title("TD 3");
	win.set_default_size(400, 400);

	MyArea area(&c, x_column, y_column);
	win.add(area);
	area.show();

	return app->run(win);
}


bool MyArea::on_draw(const Cairo::RefPtr<Cairo::Context> &cr)
{
	Gtk::Allocation allocation = get_allocation();
	const int width = allocation.get_width();
	const int height = allocation.get_height();

	// find min and max on each axis
	double x_min = DBL_MAX;
	double x_max = DBL_MIN;
	double y_min = DBL_MAX;
	double y_max = DBL_MIN;
	for(int i = 0; i < c->get_n(); i++)
	{
		if(c->get_point(i).coords[x_column] < x_min)
			x_min = c->get_point(i).coords[x_column];

		if(c->get_point(i).coords[x_column] > x_max)
			x_max = c->get_point(i).coords[x_column];

		if(c->get_point(i).coords[y_column] < y_min)
			y_min = c->get_point(i).coords[y_column];

		if(c->get_point(i).coords[y_column] > y_max)
			y_max = c->get_point(i).coords[y_column];
	}

	// plot all points
	for(int i = 0; i < c->get_n(); i++)
	{
		cr->save(); // save current drawing context (opaque black)
		cr->arc((c->get_point(i).coords[x_column]-x_min)*width/(x_max-x_min), (c->get_point(i).coords[y_column]-y_min)*height/(y_max-y_min), 10.0, 0.0, 2.0 * M_PI); // full circle

		// choose color depending on label
		switch(c->get_point(i).label)
		{
			case 0:
			cr->set_source_rgba(1.0, 0.0, 0.0, 0.6); // red, partially translucent
			break;

			case 1:
			cr->set_source_rgba(0.0, 0.0, 0.8, 0.6); // 0.8 blue, partially translucent
			break;

			case 2:
			cr->set_source_rgba(0.0, 1.0, 0.0, 0.6); // green, partially translucent
			break;

			default:
			double shade = c->get_point(i).label/(double)c->get_k();
			cr->set_source_rgba(shade, shade, shade, 1.0); // gray
			break;
		}

		cr->fill_preserve();
		cr->restore();  // restore drawing context to opaque black
		cr->stroke();
	}

	return true;
}
