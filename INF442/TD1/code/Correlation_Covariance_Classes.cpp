#include <stdio.h>
#include <math.h>
 
class StdDeviation

{

private:

    int max;

    double value[100];

    double mean;

public:

    double ComputeMean()

    {

        // Complete the function

    }

   double ComputeVariance()

    {
        mean = ComputeMean();

        // Complete the function

    }

    double ComputeSampleVariane()

    {
        mean = ComputeMean();

        // Complete the function
    }

    int SetValues(double *p, int count)

    {

        if(count > 100)

            return -1;

        max = count;

        for(int i = 0; i < count; i++)

            value[i] = p[i];

        return 0;

    }

  double ComputeStandardDeviation()
    {
        // Complete the function
    }

 double ComputeSampleStandardDeviation()
    {
        // Complete the function
    }
};

class StatCalculator

{

  private:

    double XSeries[100];

    double YSeries[100];

    int max; 

    StdDeviation x;

    StdDeviation y;

public:

    void SetValues(double *xvalues, double *yvalues, int count)

    {
        for(int i = 0; i < count; i++)

        {

            XSeries[i] = xvalues[i];

            YSeries[i] = yvalues[i];

        }

        x.SetValues(xvalues, count);

        y.SetValues(yvalues, count);

        max = count;

    }

  double ComputeCovariance()

    {

        double xmean = x.ComputeMean();

        double ymean = y.ComputeMean();

        double total = 0;

        for(int i = 0; i < max; i++)

        {

            // Complete the function

        }

        return total / max;

    }

    double ComputeCorrelation()

    {

        // Complete the function

    }
  
  void printStats(StatCalculator &sd)
    {
	// Complete the function to output the covariance and the correlation of the input data
			   
    }

};
 
int  main()

{
    StatCalculator calc;
    {

        printf("\n\nZero Correlation and Covariance Data Set\n");

        double xarr[] = { 8, 6, 4, 6, 8 };

        double yarr[] = { 10, 12, 14, 16, 18 };

        calc.SetValues(xarr,yarr,sizeof(xarr) / sizeof(xarr[0]));

       // Complete the function to output the covariance and the correlation 

    }


    {

        printf("\n\nPositive Correlation and Low Covariance Data Set\n");

        double xarr[] = { 0, 2, 4, 6, 8 };

        double yarr[] = { 6, 13, 15, 16, 20 };

        // Complete the function to output the covariance and the correlation 

    }

    {

        printf("\n\nNegative Correlation and Low Covariance Data Set\n");

        double xarr[] = { 8, 6, 4, 2, 0 };

        double yarr[] = { 6, 13, 15, 16, 20 };

        // Complete the function to output the covariance and the correlation 
    }

    {

        printf("\n\nPositive Correlation and High Covariance Data Set\n");

        double xarr[] = { 8, 6, 4, 2, 0 };

        double yarr[] = { 1006, 513, 315, 216, 120 };

        // Complete the function to output the covariance and the correlation 
    }

    {

        printf("\n\nNegative Correlation and High Covariance Data Set\n");

        double xarr[] = { 8, 6, 4, 2, 0 };

        double yarr[] = { 120, 216, 315, 513, 1006 };

       // Complete the function to output the covariance and the correlation 

    }
    return 0;

}

