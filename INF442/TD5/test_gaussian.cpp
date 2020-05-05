#include <iostream>
#include <cmath>
#include <cassert>

#include "point.hpp"
#include "cloud.hpp"
#include "gaussian.hpp"


int main()
{
	const double eps = 0.01;

	cloud c(7, 5);
	point p;

	gaussian ker(&c, 5.00);

	assert( std::abs(ker.volume() - 621.77) < eps );

	assert( std::abs(ker.profile(0.5) - 0.778801) < eps );
	assert( std::abs(ker.profile(1.5) - 0.472367) < eps );
	assert( std::abs(ker.profile(2.5) - 0.286505) < eps );

	std::cout << "tests for gaussian passed" << std::endl;

	return 0;
}
