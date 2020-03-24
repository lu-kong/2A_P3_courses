#include <iostream>
using namespace std;
#include <string>
        
int main(int argc, char **argv)
{
  std::string name;
   std::cout << "Thank you for running my first C++ programm !" << std::endl;
  std::cout << "Please enter your name: " << std::endl;
  std::cin >> name;         
  std::cout << "Hello " << name << std::endl;   
    
  return 0;
}