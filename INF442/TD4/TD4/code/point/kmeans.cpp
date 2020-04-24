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

	point()
	{
		coords = new double[d];
		for(int m = 0; m < d; m++)
			coords[m] = 0.0;
		label = 0;
	}

	~point()
	{
		delete[] coords;
	}

	void print()
	{
		std::cout << coords[0];

		for (int j = 1; j < d; j++)
			std::cout << '\t' << coords[j];

		std::cout << std::endl;
	}

	double dist(point &q)
	{
		double sqd = 0.0;

		for (int m = 0; m < d; m++)
			sqd += (coords[m]-q.coords[m]) * (coords[m]-q.coords[m]);

		return std::sqrt(sqd);
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

	void set_center(point &p, int j)
	{
		for(int m = 0; m < d; m++)
			centers[j].coords[m] = p.coords[m];
	}

	double intracluster_variance()
	{
		
		double sum = 0.0;
		for(int i = 0; i < n; i++)
		{
			sum += points[i].dist(centers[points[i].label]) * points[i].dist(centers[points[i].label]);
		}

		return sum / n;
		
	}

	void set_voronoi_labels()
	{
		
		for(int i = 0; i < n; i++)
		{
			double min_dist = DBL_MAX;
			int min_ind = -1;

			for(int j = 0; j < k; j++)
			{
				if(points[i].dist(centers[j]) < min_dist)
				{
					min_dist = points[i].dist(centers[j]);
					min_ind = j;
				}
			}

			points[i].label = min_ind;
		}
		
	}

	void set_centroid_centers()
	{
		
		int *counts = new int[k];
		for(int j = 0; j < k; j++)
			counts[j] = 0;
		for(int i = 0; i < n; i++)
			counts[points[i].label]++;


		for(int j = 0; j < k; j++)
			if(counts[j] != 0)
				for(int m = 0; m < d; m++)
					centers[j].coords[m] = 0.0;

		for(int i = 0; i < n; i++)
		{
			for(int m = 0; m < d; m++)
				centers[points[i].label].coords[m] += points[i].coords[m];
		}

		for(int j = 0; j < k; j++)
			if(counts[j] != 0)
				for(int m = 0; m < d; m++)
					centers[j].coords[m] /= counts[j];

		delete[] counts;
		
	}

	void kmeans()
	{
		
		// set_centroid_centers(); // trivial initialization
		// init_forgy();
		init_plusplus();
		// init_random_partition();

		bool stabilized = false;
		int *old_labels = new int[n];

		while(!stabilized)
		{
			for(int i = 0; i < n; i++)
				old_labels[i] = points[i].label;

			set_voronoi_labels();
			set_centroid_centers();

			bool equal = true;
			for(int i = 0; i < n; i++)
			{
				if(points[i].label != old_labels[i])
				{
					equal = false;
					break;
				}
			}

			stabilized = equal;
		}

		delete[] old_labels;
		
	}

	void init_forgy()
	{
		
		int *already_chosen = new int[k];

		for(int j = 0; j < k; j++)
		{
			// choose index i different from those already chosen
			int i;
			bool new_index = false;
			while(!new_index)
			{
				i = rand() % n;

				// check whether i was already chosen
				new_index = true;
				for(int r = 0; r < j; r++)
					if(already_chosen[r] == i)
					{
						new_index = false;
						break;
					}
			}

			already_chosen[j] = i;

			for(int m = 0; m < d; m++)
				centers[j].coords[m] = points[i].coords[m];
		}

		delete[] already_chosen;
		
	}

	void init_plusplus()
	{
		
		// choose first center
		int i = rand() % n;
		for(int m = 0; m < d; m++)
			centers[0].coords[m] = points[i].coords[m];

		// number of centers already chosen
		int nb_already_chosen = 1;

		// will hold squared distances to nearest already chosen center
		double *distances = new double[n];

		while(nb_already_chosen < k)
		{
			double sum_distances = 0.0;

			// calculate squared distance to nearest already chosen center
			for(i = 0; i < n; i++)
			{
				distances[i] = DBL_MAX;

				for(int j = 0; j < nb_already_chosen; j++)
					if(points[i].dist(centers[j])*points[i].dist(centers[j]) < distances[i])
						distances[i] = points[i].dist(centers[j])*points[i].dist(centers[j]);

				sum_distances += distances[i];
			}


			// choose random point proportional to square distance
			double random = ((double)rand() / RAND_MAX) * sum_distances;
			double prob_sum = 0.0;
			i = 0;
			while(prob_sum <= random && i < n)
			{
				prob_sum += distances[i];
				i++;
			}

			// i-1 is now the index of the chosen point
			for(int m = 0; m < d; m++)
				centers[nb_already_chosen].coords[m] = points[i-1].coords[m];

			nb_already_chosen++;
		}

		delete[] distances;
		
	}

	void init_random_partition()
	{
		
		for(int i = 0; i < n; i++)
		{
			points[i].label = rand() % k;
		}

		set_centroid_centers();
		
	}
};


// test functions
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
		test_intracluster_variance();
		test_kmeans();
		test_init_forgy();
		test_init_plusplus();
		test_init_random_partition();
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
