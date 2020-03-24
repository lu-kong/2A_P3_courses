#include <iostream> 

void incr (int& r) {
  r++;
}

int main(){
    int a = 1;
    std::cout<<"a origin : "<< a << std::endl;
}