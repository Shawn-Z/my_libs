#include "SVector.hpp"
#include "STime.hpp"
#include "SJerk.hpp"
#include "SHandle.hpp"
#include "SLog.hpp"


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
    shawn::SHandle sHandle;
    shawn::SLog sLog;
    sLog.init("haha");
    uint8_t dta[3] = {1,2,3};
    std::cout<< sizeof(dta)<<std::endl;
    sLog.logUint8Array((char *) dta, sizeof(dta), google::INFO);
    sLog.logUint8Array((char *) dta, sizeof(dta), google::INFO);
    sLog.logUint8Array((char *) dta, sizeof(dta), google::INFO);
    sLog.logUint8Array((char *) dta, sizeof(dta), google::INFO);
    sLog.logUint8Array((char *) dta, sizeof(dta), google::INFO);
    sLog.logUint8Array((char *) dta, sizeof(dta), google::INFO);

//    sHandle.setID(100);


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


//    shawn::STime sTime;
//    sTime.updateTimestampsTillNowMs();
//    sTime.updateTimestampsTillNowMs();
//    sTime.updateTimestampsTillNowMs();
//    sTime.updateTimestampsTillNowMs();
//    sTime.updateTimestampsTillNowMs();
//    sTime.updateTimestampsTillNowMs();
//    sTime.updateTimestampsTillNowMs();
//    sTime.updateTimestampsTillNowMs();
//    sTime.updateTimestampsTillNowMs();

    return 0;
}