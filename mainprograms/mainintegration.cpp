

#include <iostream>
#include <math.h>
#include "IntRule.h"
#include "IntRule1d.h"
#include "IntRuleQuad.h"
#include "IntRuleTetrahedron.h"
#include "IntRuleTriangle.h"
#include "Topology1d.h"
#include "TopologyTriangle.h"
#include "TopologyQuad.h"
#include "TopologyTetrahedron.h"
#include "DataTypes.h"
#include "Analysis.h"
#include "VTKGeoMesh.h"
#include "ReadGmsh.h"

using std::cout;
using std::endl;
using std::cin;

//void Integrate1D();
void IntegrateTriangle();
//void Integrate2D();


int main() {
    //Integrate1D();
    //Integrate2D();
    IntegrateTriangle();
   
    return 0;
}

void IntegrateTriangle(){
    
    auto func = [](VecDouble x){ return x[0]*x[0]*x[1]*x[1];};
    IntRuleTriangle quadrule(2);
    const int np = quadrule.NPoints();
    double integral = 0.0;
    VecDouble co(2);
    double weight;
    for (int ip=0; ip < np; ip++) {
        quadrule.Point(ip, co, weight);
        //std::cout << "ip  = "<< ip << " co ="<<co << " weight = "<<weight << std::endl;
        double val = func (co);
        //std::cout << "func value = "<< val << std::endl;
        integral += val*weight;

    }
    std::cout << "se espera 0.0007716, se obtiene " << integral << std::endl;
}




void Integrate2D(){
    
    //func: func lamba, es la función que vamos a integrar, en este caso x2 y2
    auto func = [](VecDouble x){ return x[0]*x[0]*x[1]*x[1];};
    IntRuleQuad quadrule(2 /*Order*/);
    const int np = quadrule.NPoints();
    double integral = 0.0;
    VecDouble co(2);
    double weight;
    for (int ip=0; ip < np; ip++) {
        quadrule.Point(ip, co, weight);
        std::cout << "ip  = "<< ip << " co ="<<co << " weight = "<<weight << std::endl;
        double val = func (co);
        integral += val*weight;
    }
    std::cout << "se espera 4./9, se obtiene" << integral << std::endl;
}

void Integrate1D (){

    // test an integration rule
    // lambda expression
    auto func = [](double x){return x*x;};
    //int val=3;
    IntRule1d gab;

   //Creando un objeto tipo Intrule y seteandolo el orden 2    
    IntRule1d oned(4);
   // auto val2 = oned.gabriel(3);
   //int np = oned.NPoints(); (Entrega el num de puntos o la fila que se va a emplear en la integración Gaussiana. Para el orden especificado en oned (2);)

    int np = oned.NPoints();
    double integral = 0.;
    VecDouble co(1);
    double weight;
    for (int ip = 0; ip < np; ip++) {
        
        oned.Point(ip, co, weight);
        double val = func(co[0]);
        integral += val*weight;
        double valoresperado=2/3;
        
        std::cout << "punto de integración = " << co[0] << "  pesos = " << weight << " np  = " << np << std::endl;
    }
    std::cout << "espera se "<< 2./3 <<" obtem se " << integral << std::endl; 
       
}
