#include "AreaChange.h"
#include <iomanip>
#include <iostream>



AreaChange::AreaChange(int a, double p) {

    area = a;
    position = p;
}

AreaChange::AreaChange(AreaChange& a) {
    
    area = a.area;
    position = a.position;
}

AreaChange::~AreaChange(void) {
    
}

bool AreaChange::operator<(const AreaChange& a) const {

    if ( position < a.position )
        return true;
    return false;
}

void AreaChange::print(void) {

    std::cout << std::setw(4) << area << " " << std::fixed << std::setprecision(6) << position << std::endl;
}