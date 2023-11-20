/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

///\cond
#include <iostream> 
///\endcond
#include "IntRuleTriangle.h"


IntRuleTriangle::IntRuleTriangle(){
    SetOrder(0);

}

IntRuleTriangle::IntRuleTriangle(int order) :IntRule(order){
    SetOrder(order);
    
}

void IntRuleTriangle::SetOrder(int order) {
    std::cout << __PRETTY_FUNCTION__ << " needs to be implemented\n";
    DebugStop();
}
