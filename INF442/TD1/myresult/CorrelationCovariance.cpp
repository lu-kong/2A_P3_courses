#include <iostream>
#include <string>
#include <math.h>


/*Read data matrix from the file "wines.txt" given as a parameter while runing the program : "./Stat_Covariance < wines.txt" */
/* You Should output to a file "./Stat_Covariance < wines.txt >> Results_Stat.txt" */
/* Open "Results_Stat.txt" to read the runing results*/

void readData(double** DataMatrix, int N_REGIONS,  int M_WINES)
{  
    for (int i = 0; i < N_REGIONS; i++)
      {
	DataMatrix[i] = new double[M_WINES];
	for (int j = 0; j < M_WINES; j++)
	  std::cin >> DataMatrix[i][j];
      }
  }


double CalculateMean(double values[], int length) {
  double sum = 0;
  for(int i=0; i<length;i++){
    sum += values[i];
  }
  // Complete the function
  return (sum / length);
}

double CalculateVariance(double values[], int length) {
  double mean = CalculateMean(values, length);
  // Complete the function
  double temp = 0;
  for(int i=0;i<length;i++){
    temp += (values[i]-mean)*(values[i]-mean);
  }
  return temp / length;
}


/*Print the readed data matrix on the screen */

void printData(double**  DataMatrix, int N_REGIONS, int M_WINES){
  std::cout << "Let's look at our dataset: " <<std::endl;
  std::cout << "----------------------------------------------------------------------" << std::endl;
  //Complete the function
  for(int i=0;i<N_REGIONS;i++){
    for(int j=0;j<M_WINES;j++){
      std::cout << DataMatrix[i][j] << "\t";
    }
    std::cout << '\n'; 
  }
    
}

double CalculateSampleVariance(double values[], int length) {
  double mean = CalculateMean(values, length);
  // Complete the function
  double temp = 0;
  for(int i=0;i<length;i++){
    temp += (values[i]-mean)*(values[i]-mean);
  }
  return temp / (length - 1);
}

double GetStandardDeviation(double values[], int length) {
  return sqrt(CalculateVariance(values,length));// Complete the function;
}

double GetSampleStandardDeviation(double values[], int length) {
  return sqrt(CalculateSampleVariance(values,length));// Complete the function;
}

void StatCalcul_vector(double sample[], int sample_length) {

  double mean = CalculateMean(sample, sample_length); 
  double variance = CalculateVariance(sample, sample_length); 
  double samplevariance = CalculateSampleVariance(sample, sample_length); 
  double sampledevi = GetSampleStandardDeviation(sample, sample_length); 
  double devi = GetStandardDeviation(sample, sample_length);

  std::cout << "Mean \t\t = " << mean << ", \t\t Variance = \t\t" << variance<< std::endl;

  /*std::cout << "Total Numbers\t\t\t: " << sample_length << "\n"
            << "Mean\t\t\t\t: " << mean << "\n"
            << "Population Variance\t\t: " << variance << "\n"
            << "Sample variance\t\t\t: " << samplevariance << "\n"
            << "Population Standard Deviation\t: " << devi << "\n"
            << "Sample Standard Deviation\t: " << sampledevi << "\n";*/
}

double CalculateCovariance(double values1[], double values2[], int length) {

  double xmean = CalculateMean(values1, length);
  double ymean = CalculateMean(values2, length);
  double total = 0;
  for(int i=0; i<length;i++){
    total+=(values1[i]-xmean)*(values2[i] - ymean);
  }
  return total/length;
  //Complete the function
}

/*Compute Correlation*/


double CalculateCorrelation(double values1[], double values2[], int length)
  {
    double cov = CalculateCovariance(values1, values2, length);
    double sd1 = GetStandardDeviation(values1, length);
    double sd2 = GetStandardDeviation(values2, length);
    return cov/sd1/sd2;
    //Complete the function
    
  }

/* get the row of index "row_index" from the data matrix*/
  
void getRow(double** DataMatrix,  int row_index, int M_WINES, double row[]){
  
    for (int j = 0; j < M_WINES ; j++)
      
      {
	     //Complete the function
      row[j] = DataMatrix[row_index][j];
	
	/*std::cout << "in function getRow: row[" << j << "]="<< row[j] << " " << std::endl;*/
      }
  }

/* get the column of index "column_index" from the data matrix*/
  
void getColumn(double** DataMatrix, int N_REGIONS, int column_index, double column[]){
    for (int i = 0; i < N_REGIONS; i++)
      {
      
	      //Complete the function
	    column[i] = DataMatrix[i][column_index];
	/*std::cout << "in function getColumn: column[" << i << "]="<< column[i] << " " << std::endl;*/
      }
  }

/* Covariance matrix */

void ComputeCovarianceMatrix(double** DataMatrix, double** CovarianceMatrix, int M_WINES, int N_REGIONS){
    
  int column_index1 = 0;
  int column_index2 = 0;
  double column1[N_REGIONS];
  double column2[N_REGIONS];

 for (int j = 0; j <  M_WINES; j++)
      {
        column_index1 = j;
        CovarianceMatrix[column_index1] = new double [M_WINES];
        
        getColumn(DataMatrix,N_REGIONS,column_index1,column1);
        for(int i=j+1;i<M_WINES;i++)
        {
          column_index2 = i;
          getColumn(DataMatrix,N_REGIONS,column_index2,column2);
          CovarianceMatrix[column_index1][column_index2] = CalculateCovariance(column1,column2,N_REGIONS);
        }
	//Complete the function
	
      }
 
 for (int j =0; j < M_WINES; j++) /* Symetric Values*/
      {
        double column[N_REGIONS];
        getColumn(DataMatrix,N_REGIONS,j,column);
        CovarianceMatrix[j][j] = CalculateVariance(column,N_REGIONS);
      for(int i=0;i<j;i++)
      {
        CovarianceMatrix[j][i] = CovarianceMatrix[i][j];
      }
	//Complete the function
	
      }
}

/* Print Covariance Matrix*/
void printCovarianceMatrix(double** CovarianceMatrix, int M_WINES){

  //  std::cout << "Let's look at our dataset: " << std::endl;
  //  std::cout << "----------------------------------------------------------------------" << std::endl;
   //Complete the function
   for (int i = 0; i < M_WINES; i++)
   {
     for (int j = 0; j < M_WINES; j++)
     {
       std::cout << CovarianceMatrix[i][j] << "\t";
     }
     std::cout << '\n';
   }
  } 


int main(int argc, char** argv) {
  
  int M_WINES=0;
  int N_REGIONS=0;  
  
  std::cin >> M_WINES;
  std::cin >> N_REGIONS;

  double** DataMatrix = 0;
  DataMatrix = new double*[N_REGIONS];
  
  double sample1[M_WINES];  /*for contry statistical calculation*/
  double sample2[N_REGIONS]; /* for wine statistical calculation*/
  
  readData(DataMatrix, N_REGIONS,  M_WINES);
  
std::cout << "----------------------------------------------------------------------" << std::endl;


  /* Print Data Matrix */

  printData(DataMatrix, N_REGIONS, M_WINES);

  std::cout << "----------------------------------------------------------------------" << std::endl;

  
    /*int DataMatrix_length = sizeof(DataMatrix) / sizeof(DataMatrix[0]); // usual trick*/


  /*Statiscal Calculation for the country "row_index" */
  
std::cout << "----------------------------------------------------------------------" << std::endl;
 
  std::cout << "Countries Mean and Variance:" << std::endl;
 std::cout << "----------------------------------------------------------------------" << std::endl;
 
for (int i = 0; i < N_REGIONS; i++)
  {
    int row_length = M_WINES;
    std::cout << i<< "  " ;
    for (int j = 0; j < M_WINES; j++)
      {
	sample1[j] =  DataMatrix[i][j];
	
      }
    StatCalcul_vector(sample1, row_length);
  }   


 
std::cout << "----------------------------------------------------------------------" << std::endl;
/*Statiscal Calculation for the wine "column_index" */

 std::cout << "Wines Mean and Variance:" << std::endl;
std::cout << "----------------------------------------------------------------------" << std::endl;

 for (int j = 0; j < M_WINES; j++)
    {

      int column_length = N_REGIONS;
       std::cout << j << "  " ;
      for (int i = 0; i < N_REGIONS; i++)
      {
	sample2[i] =  DataMatrix[i][j];
      }
    StatCalcul_vector(sample2, column_length);
  }
 
 /*Countries Correlation*/
std::cout << "----------------------------------------------------------------------" << std::endl;

 std::cout << " Countires Correlation:" << std::endl;
std::cout << "----------------------------------------------------------------------" << std::endl;

 int row_index1 = 0;
 int row_index2 = 0;
 double row1[M_WINES];
 double row2[M_WINES];
 
 
  for (int i = 0; i <  N_REGIONS; i++)
      {
	row_index1 = i;
	getRow(DataMatrix, row_index1, M_WINES, row1);
	for (int k = i+1; k <  N_REGIONS; k++)
	  {
	    row_index2 = k;
	    getRow(DataMatrix, row_index2, M_WINES, row2);
	     std::cout << "Correlation[" << i << ", " << k << "] = " << CalculateCorrelation(row1, row2, M_WINES) << std::endl;
	  }
  }

  /*Countries Covariance*/
std::cout << "----------------------------------------------------------------------" << std::endl;

 std::cout << " Countires Covariance:" << std::endl;
std::cout << "----------------------------------------------------------------------" << std::endl;

 for (int i = 0; i <  N_REGIONS; i++)
      {
	row_index1 = i;
	getRow(DataMatrix, row_index1, M_WINES, row1);
	for (int k = i+1; k <  N_REGIONS; k++)
	  {
	    row_index2 = k;
	    getRow(DataMatrix, row_index2, M_WINES, row2);
	     std::cout << "Correlation[" << i << ", " << k << "] = " << CalculateCovariance(row1, row2, M_WINES) << std::endl;
	  }
  }

 
 
 /* Wines Correlation*/

 std::cout << "----------------------------------------------------------------------" << std::endl;

 std::cout << "Wines Correlation:" << std::endl;
std::cout << "----------------------------------------------------------------------" << std::endl;

 int column_index1 = 0;
 int column_index2 = 0;
 double column1[N_REGIONS];
 double column2[N_REGIONS];
 
 
  for (int j = 0; j <  M_WINES; j++)
      {
	column_index1 = j;
	getColumn(DataMatrix, N_REGIONS, column_index1, column1);
	for (int k = j+1; k <  M_WINES; k++)
	  {
	    column_index2 = k;
	    getColumn(DataMatrix, N_REGIONS, column_index2, column2);
	     std::cout << "Correlation[" << j << ", " << k << "] = " << CalculateCorrelation(column1, column2, N_REGIONS) << std::endl;
	  }
  }

  /* Wines Covariance*/

 std::cout << "----------------------------------------------------------------------" << std::endl;

 std::cout << "Wines Covariance:" << std::endl;
std::cout << "----------------------------------------------------------------------" << std::endl;



 for (int j = 0; j <  M_WINES; j++)
      {
	column_index1 = j;
	getColumn(DataMatrix, N_REGIONS, column_index1, column1);
	for (int k = j+1; k <  M_WINES; k++)
	  {
	    column_index2 = k;
	    getColumn(DataMatrix, N_REGIONS, column_index2, column2);
	     std::cout << "Covariance[" << j << ", " << k << "] = " << CalculateCovariance(column1, column2, N_REGIONS) << std::endl;
	  }
  }

  /* Covaraince Matrix */
  double** CovarianceMatrix = 0;
  
  CovarianceMatrix = new double*[M_WINES];
  ComputeCovarianceMatrix(DataMatrix, CovarianceMatrix, M_WINES, N_REGIONS);

  /*Print Covariance Matrix */

  std::cout << "----------------------------------------------------------------------" << std::endl;

 std::cout << "Wines Covariance Matrix:" << std::endl;
std::cout << "----------------------------------------------------------------------" << std::endl;
 
 printCovarianceMatrix(CovarianceMatrix, M_WINES);
  
  return 0;
}

