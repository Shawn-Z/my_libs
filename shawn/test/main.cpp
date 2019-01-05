#include "SVector.hpp"
#include "STime.hpp"
#include "SJerk.hpp"
#include "SHandle.hpp"
#include "SLog.hpp"
#include "SProportion.hpp"


#include <iostream>
#include <cmath>
#include <cfloat>

void fun(int i){
    static std::vector<int> vec(20, 0);
    vec.erase(vec.begin());
    vec.emplace_back(i);
    for (int j = 0; j < vec.size(); ++j) {
        std::cout<<vec[j]<<" ";
    }
    std::cout<<std::endl;
}

int main() {
    fun(1);
    fun(2);
    fun(3);
    fun(4);
    fun(5);
    fun(6);
    fun(7);
    fun(8);
    fun(9);
    fun(10);
    fun(11);
    fun(12);
    fun(13);
    fun(14);
    fun(15);
    fun(16);
    fun(17);
    fun(18);
    fun(19);
    fun(20);
    fun(10);
    fun(10);
    fun(10);
    return 0;
}