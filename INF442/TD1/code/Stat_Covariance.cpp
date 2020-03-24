#include <stdio.h>

#include <math.h>

 #include <iostream>

class StdDeviation

{

private:

    int max;

    double value[100];

    double mean;

   

public:

  double ComputeMean()

    {

        double sum = 0;

        for(int i = 0; i < max; i++)

            sum += value[i];

        return (sum / max);

    }

 

    double ComputeVariance()

    {

        mean = ComputeMean();

 

        double temp = 0;

        for(int i = 0; i < max; i++)

        {

             temp += (value[i] - mean) * (value[i] - mean) ;

        }

        return temp / max;

    }

 

    double ComputeSampleVariance()

    {

        mean = ComputeMean();

 

        double temp = 0;

        for(int i = 0; i < max; i++)

        {

             temp += (value[i] - mean) * (value[i] - mean) ;

        }

        return temp / (max - 1);

    }

 

    int SetValues(int p[], int count)

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

        return sqrt(ComputeVariance());

    }

 

    double ComputeSampleStandardDeviation()

    {

        return sqrt(ComputeSampleVariance());

    }

 

};

 

class StatCalculator

{

private:

  int _XSeries[100];
  
  int _YSeries[100];

  int _max;

  int _M_WINES;

  int _N_REGIONS;

  StdDeviation x;

  StdDeviation y;
 

public:
  
  StatCalculator( int M_WINES, int N_REGIONS): _M_WINES(M_WINES), _N_REGIONS(N_REGIONS) {};
  ~StatCalculator() {};
  
  void SetValues(int xvalues[], int yvalues[], int count)
  {
    for(int i = 0; i < count; i++)
      {
	_XSeries[i] = xvalues[i];

	_YSeries[i] = yvalues[i];
      }

    x.SetValues(xvalues, count);
    
    y.SetValues(yvalues, count);
    
    _max = count;
  }

  double ComputeCovariance()
  {
    double xmean = x.ComputeMean();
    double ymean = y.ComputeMean();
    double total = 0;

    for(int i = 0; i < _max; i++)
      {
	total += (_XSeries[i] - xmean) * (_YSeries[i] - ymean);	
      }

    return total / _max;
  }

  double ComputeCorrelation()
  {
    double cov = ComputeCovariance();
    double correlation = cov / (x.ComputeStandardDeviation() * y.ComputeStandardDeviation());

    return correlation;
  }

  /*Read data matrix from the file "wines.txt" given as a parameter while runing the program : ./Stat_Covariance < wines.txt*/

  void readData(int** DataMatrix)
  {  
    for (int i = 0; i < _N_REGIONS; i++)
      {
	DataMatrix[i] = new int[_M_WINES];
	for (int j = 0; j < _M_WINES; j++)
	  std::cin >> DataMatrix[i][j];
      }
  }

  /*Print the readed data matrix on the screen*/
  
  void printData(int**  DataMatrix){
    std::cout << "Print data matrix" << std::endl;  
    
        // *************   Complete the function  ***************
    
  }


  /* Print Means and Variance */

  void printSatResults(int** DataMatrix){

    int row1[_N_REGIONS];
    int column1[_M_WINES];
     int row2[_N_REGIONS];
    int column2[_M_WINES];
        
 
    /*for (int i =0; i < _N_REGIONS; i++)
    
	  {
	    getRow(DataMatrix, i, row1);
	    getRow(DataMatrix, i+1, row2);
	    SetValues(row1, row2, _M_WINES);
	    std::cout << "MeanCountry[" << i << "] = " << x.ComputeMean() << " " << std::endl;
	    }*/
    
  for (int j =0; j < _M_WINES; j++)
    
	  {
	    getColumn(DataMatrix, j, column1);
	    getColumn(DataMatrix, j+1, column2);
	    SetValues(column1, column2,_N_REGIONS);
	    std::cout << "MeanWine[" << j << "] = " << x.ComputeMean() << " " << std::endl;
	  }

  
  /* for (int i =0; i < _N_REGIONS; i++)
    
	  {
	    getRow(DataMatrix, i, row1);
	    getRow(DataMatrix, i+1, row2);
	    SetValues(row1, row2, _M_WINES);
	    std::cout << "VarianceCountry[" << i << "] = " << x.ComputeVariance() << " " << std::endl;
	    }*/
    

  for (int j =0; j < _M_WINES; j++)
    
	  {
	    getColumn(DataMatrix, j, column1);
	    getColumn(DataMatrix, j+1, column2);
	    SetValues(column1, column2,_N_REGIONS);
	    std::cout << "VarianceWine[" << j << "] = " << x.ComputeVariance() << " " << std::endl;
	  }

  for (int j =0; j < _M_WINES; j++)
    
	  {
	    getColumn(DataMatrix, j, column1);
	    getColumn(DataMatrix, j+1, column2);
	    SetValues(column1, column2,_N_REGIONS);
	    std::cout << "CorrelationWine[" << j << " ," << j+1 <<" ] = " << ComputeCorrelation() << " " << std::endl;
	  }
  
  } /*end function print Mean Variance*/
  

  /* get the row of index "row_index" from the data matrix*/
  
  void getRow(int** DataMatrix, int row_index, int row[]){
    for (int j = 0; j < _M_WINES; j++)
      
      {
	 // *************   Complete the function  ***************
	
      }
  }

  /* get the column of index "column_index" from the data matrix*/
  
  void getColumn(int** DataMatrix, int column_index, int column[]){
    for (int i = 0; i < _N_REGIONS; i++)
      {
	 // *************   Complete the function  ***************
      }
  }
  
  
  /* Compute the Covariance Matrix */
 
  
  void ComputeCovarianceMatrix(int** DataMatrix, double** CovarianceMatrix){
    int column1[_N_REGIONS];
    int column2[_M_WINES];
    
    for (int j =0; j < _M_WINES; j++)
      {
	CovarianceMatrix[j] = new double[_M_WINES];
	getColumn(DataMatrix, j, column1);
	CovarianceMatrix[j][j] = x.ComputeVariance();
	/*std::cout << "in Covariance: CovarianceMatrix[" << j << " , ] = " << CovarianceMatrix[j][j] << " " << std::endl;*/
	
	for (int k = j+1; k < _M_WINES ; k++)
	  {
	     // *************   Complete the function  ***************    
	  }
      }

    /*Symetric values*/
    for (int j =0; j < _M_WINES; j++)
      {
	for (int k = 0; k < j ; k++)
	  {	    
	     // *************   Complete the function  ***************	    
	  }
      }
    
  }

  /* Output the obtained Covariance Matrix*/

void printCovarianceMatrix(double** CovarianceMatrix){

    // *************   Complete the function  ***************
  } 
  
  
}; /*end class StatCalculator*/

 
int main(int argc, char** argv)

{
  printf("\n\n Correlation and Covariance Data Set\n");
 
  int M_WINES;
  int N_REGIONS;  
  
  std::cin >> M_WINES;
  std::cin >> N_REGIONS;
  
  
  StatCalculator calc(M_WINES, N_REGIONS);

  {
  std::cout << "DataMatrix [ "<< M_WINES  <<" Wines , " << N_REGIONS << " Countries]"<< std::endl;
  
  int** DataMatrix = 0;
  DataMatrix = new int*[N_REGIONS];

  double** CovarianceMatrix = 0;
  CovarianceMatrix = new double*[M_WINES];
 
  int row1[M_WINES];
  int column1[N_REGIONS];
  int row2[M_WINES];
  int column2[N_REGIONS];

  /* Read data matrix from wine.txt file*/
  
  calc.readData(DataMatrix);
  
  /* Output data matrix*/
  
  calc.printData(DataMatrix);


  calc.printSatResults(DataMatrix);
  
  /* Extract from the data matrix the vectors: row1 and row2, column1 and column2*/
    
  for (int j = 0; j < N_REGIONS; j++)
    {
      row1[j]=0;
      row2[j]=0;
    }

  for (int i = 0; i < M_WINES; i++)
    {
      column1[i]=0;
      column2[i]=0;
    }
  
  int row_index1 = 0;
  int row_index2 = row_index1 + 1;
  int column_index1 = 0;
  int column_index2 = column_index1 + 1;
  
  calc.getRow(DataMatrix, row_index1, row1);
  
  /*for (int j = 0; j < M_WINES; j++)
    {
    std::cout << "in main getRow results : row1[" << j << "]="<< row1[j] << " " << std::endl;
    }*/
  
  calc.getRow(DataMatrix, row_index2, row2);
  
  /*for (int j = 0; j < M_WINES; j++)
    {
    std::cout << "in main getRow results : row1[" << j << "]="<< row2[j] << " " << std::endl;
    }*/
  
  calc.getColumn(DataMatrix, column_index1, column1);
  
  /*for (int i = 0; i < N_REGIONS; i++)
    {
    std::cout << "in main getColoumn results : column1[" << i << "]="<< column1[i] << " " << std::endl;
    }*/
  
  calc.getColumn(DataMatrix, column_index2, column2);
  
  /*for (int i = 0; i < N_REGIONS; i++)
    {
    std::cout << "in main getColoumn results : column2[" << i << "]="<< column2[i] << " " << std::endl;
    }*/
   
  /* Compute and print covaraince and correlation of row1 and row2 */
   
  calc.SetValues(row1,row2,sizeof(row1) / sizeof(row1[0])); 
  
  printf("Covariance = %.10lf\n",  calc.ComputeCovariance());
  
  printf("Correlation = %.10lf\n", // *************   Complete the function  ***************

   /* Compute and print covaraince and correlation of column1 and column2 */
  
  calc.SetValues(column1, column2,sizeof(column1) / sizeof(column1[0])); 
	
  printf("Covariance = %.10lf\n",  // *************   Complete the function  ***************);
  
  printf("Correlation = %.10lf\n",  // *************   Complete the function  ***************);
  
  
  /*Compute Covariance Matrix*/

  /*std::cout << "In main function Covariance Matrix:" << std::endl;*/
  
   
  calc.ComputeCovarianceMatrix(DataMatrix, CovarianceMatrix);


  /*Print the Covariance Matrix*/
  
    std::cout << "Ouput Covariance Matrix" << std::endl;
    
    calc.printCovarianceMatrix(CovarianceMatrix);  
 	
  
  
  }
  return 0;/*The end of main*/
}

