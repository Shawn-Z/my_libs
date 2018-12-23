#include "SVector.hpp"
#include "STime.hpp"
#include "SJerk.hpp"


#include <iostream>
#include <cmath>
#include <cfloat>

void fun(std::vector<double> * p) {
    std::cout<<(*p).size()<<std::endl;
    auto begin = (*p).begin();
    auto end = (*p).end();
    std::cout<<*begin<<std::endl;
    std::cout<<*(end-1)<<std::endl;

}


int main() {
//    std::vector<double> a{1, 2, 3, 4};
//    std::cout<<a.size()<<std::endl;
//    fun(&a);

//    shawn::SVector sVector;
//    std::vector<char> b{'a'};
//    std::vector<double> c{1};
//    std::vector<double> c1{1};
//    double aa = sVector.get_max_absolute_diff<double_t>(c1, c);
//    std::cout<<DBL_MAX*2<<std::endl;
//    shawn::SJerk sJerk;


    shawn::STime sTime;
    sTime.updateTimestampsTillNow();
    sTime.updateTimestampsTillNow();
    sTime.updateTimestampsTillNow();
    sTime.updateTimestampsTillNow();
    sTime.updateTimestampsTillNow();
    sTime.updateTimestampsTillNow();
    sTime.updateTimestampsTillNow();
    sTime.updateTimestampsTillNow();
    sTime.updateTimestampsTillNow();

    return 0;
}