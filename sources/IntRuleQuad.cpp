/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

///\cond
#include <iostream> 
///\endcond
#include "IntRule1d.h"
#include "IntRuleQuad.h"

IntRuleQuad::IntRuleQuad(){
     SetOrder(0);
}

IntRuleQuad::IntRuleQuad(int order) :IntRule(order) {
    SetOrder(order);
}


void IntRuleQuad::SetOrder(int order) {
    if(order<0 || order>MaxOrder()){
        DebugStop();
    };
    fOrder = order;
    int np =(order/2 +1);
    VecDouble co(2*np);
    fWeights.resize(np);
    gaulegQuad(-1., 1., co, fWeights);
    np = fWeights.size();
    this->fPoints.resize(np,2);
    for(int ip=0;ip<np;ip++){
        fPoints(ip,0)=co[ip];
        fPoints(ip,1)=co[ip+np];
    }
}

void IntRuleQuad::gaulegQuad(const double x1, const double x2, VecDouble &co, VecDouble &w) {
    IntRule1d x;
    IntRule1d y;
    
    int n = w.size();   
    VecDouble cox(n);
    VecDouble coy(n);
    VecDouble wx(n);
    VecDouble wy(n);


    x.gauleg(x1, x2, cox, wx);
    y.gauleg(x1, x2, coy, wy);
    
    // co debería tener una "," porque ahora será un par ordenado???
    co.resize(2*n*n);
    w.resize(n * n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            co[j + i * n] = cox[j];
            co[j + i * n + n * n] = coy[i];
            w[n * i + j] = wx[i] * wy[j];
        }
    }
}
