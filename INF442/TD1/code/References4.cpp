#include <iostream>
 void swap(int a, int b)
{
  int s=a;
  a=b;
  b=s;
}

int main(int argc, char **argv)
{
 int a = 0;
 int b = 0;
 
 std::cin >> a;
 std::cin >> b;

  swap(a,b);
 
  std::cout << a << " " << b << std::endl;
  return 0;
}