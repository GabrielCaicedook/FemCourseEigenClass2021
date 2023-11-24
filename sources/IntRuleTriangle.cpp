/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

///\cond
#include <iostream> 
///\endcond
#include "IntRuleTriangle.h"
#include "IntRule1d.h"
#include "IntRuleQuad.h"


IntRuleTriangle::IntRuleTriangle(){
    SetOrder(0);

}

IntRuleTriangle::IntRuleTriangle(int order) :IntRule(order){
    SetOrder(order);
    
}

void IntRuleTriangle::SetOrder(int order) {
    VecDouble co, w;
    
    fOrder = order;
    if (order < 0 || order > MaxOrder()) DebugStop();

    if (order == 0 || order == 1) {
        fPoints.resize(1, 2);
        fWeights.resize(1);

        fPoints(0, 0) = 0.333333333333333;
        fPoints(0, 1) = 0.333333333333333;
        fWeights[0] = 1/2.;

    } else if (order == 2 ) {
        fPoints.resize(3, 2);
        fWeights.resize(3);

        fPoints(0, 0) = 0.166666666666667;
        fPoints(0, 1) = 0.166666666666667;
        fWeights[0] = 0.333333333333333/2.;

        fPoints(1, 0) = 0.166666666666667;
        fPoints(1, 1) = 0.166666666666667;
        fWeights[1] = 0.333333333333333/2.;

        fPoints(2, 0) = 0.166666666666667;
        fPoints(2, 1) = 0.166666666666667;
        fWeights[2] = 0.333333333333333/2;

    } else if (order == 3) {
        fPoints.resize(4, 2);
        fWeights.resize(4);

        fPoints(0, 0) = 0.333333333333333;
        fPoints(0, 1) = 0.333333333333333;
        fWeights[0] = -0.5625/2.;

        fPoints(1, 0) = 0.2;
        fPoints(1, 1) = 0.6;
        fWeights[1] = 0.520833333333333/2.;

        fPoints(2, 0) = 0.2;
        fPoints(2, 1) = 0.2;
        fWeights[2] = 0.520833333333333/2;

        fPoints(3, 0) = 0.6;
        fPoints(3, 1) = 0.2;
        fWeights[3] = 0.520833333333333/2;
    } else if (order == 4) { 
        fPoints.resize(6, 2);
        fWeights.resize(6);
        fPoints(0, 0) = 0.445948490915965; 
        fPoints(0, 1) = 0.108103018168070; 
        fWeights[0] = 0.223381589678011/2.;
        
        fPoints(1, 0) = 0.445948490915965; 
        fPoints(1, 1) = 0.445948490915965; 
        fWeights[1] = 0.223381589678011/2.;
        
        fPoints(2, 0) = 0.108103018168070; 
        fPoints(2, 1) = 0.445948490915965; 
        fWeights[2] = 0.223381589678011/2.;
        
        fPoints(3, 0) = 0.091576213509771; 
        fPoints(3, 1) = 0.816847572980459; 
        fWeights[3] = 0.109951743655322/2.;
        
        fPoints(4, 0) = 0.091576213509771; 
        fPoints(4, 1) = 0.091576213509771; 
        fWeights[4] = 0.109951743655322/2.;
        
        fPoints(5, 0) = 0.816847572980459; 
        fPoints(5, 1) = 0.091576213509771; 
        fWeights[5] = 0.109951743655322/2.;
    }
    
}


