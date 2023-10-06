#include <iostream>
#include <math.h>
#include "VTKGeoMesh.h"
#include "Geom1d.h"
#include "CompElementTemplate.h"
#include "Shape1d.h"
#include "CompMesh.h"
#include "GeoElementTemplate.h"
#include "GeoMesh.h"
#include "Poisson.h"

using namespace std;


void CallStiffFunc();

int main() {
    CallStiffFunc();
    return 0;
}

void CallStiffFunc(){
    
    int ordenP=1;
    const int dim = 1;
    const int numnodes = 2;
    
    //const int numnodes = 3;

    GeoMesh gmesh;
    gmesh.SetDimension(dim);
    gmesh.SetNumNodes(numnodes);
    GeoNode gnod0, gnod1, gnod2;
    VecDouble co0(3), co1(3), co2(3);
    co0[0]=0.0;
    co1[0]=16.0;
    //co2[0]=32.0;

//    co2[0]=6.0;
   // co0 << 0.;
   // co1 << 4.;
    gnod0.SetCo(co0);
    gnod1.SetCo(co1);
//    gnod2.SetCo(co2);
    
    
    gmesh.Node(0) = gnod0;
    gmesh.Node(1) = gnod1;
  //  gmesh.Node(2) = gnod2;
    
    int materialid = 0;
    //int materialid1 = 1;

    VecInt nodeindices(2), nodeindices1(2);
    nodeindices[0]=0;
    nodeindices[1]=1;

    // nodeindices1[0]=1;
    //nodeindices1[1]=2;


    
    int index = 0;
    //int index = 1; 

//  const VecInt &nodeindices, int materialid, GeoMesh *gmesh, int index

//GeoElementTemplate Crea un elemento geometrico del tipo < > (que se especifique dentro de < >)
//"geo" es el objeto (un elemento) que estamos creando que va a tomar como argumento (nodeindices, material, &gmesh, index)
//nodeindices = son los nodos asociados al elemento geo
//materialid = es el material de ese elemento
//&gmesh = es la malla a la cual estará asociada ese elemento
//index = (identificación del elemento) 

    GeoElementTemplate<Geom1d> geo(nodeindices,materialid,&gmesh,index);
//GeoElementTemplate<Geom1d> geo1(nodeindices1, materialid1, &gmesh, index1);

// Crear un segundo elemento en la malla geometrica

//
    CompMesh cmesh(&gmesh);
    cmesh.SetDefaultOrder(ordenP);
    // Criando material
    MatrixDouble perm(3,3);
    perm.setIdentity();
    perm(0,0) = 1.;
    Poisson poi(materialid, perm);
    cmesh.SetMathStatement(materialid, &poi);
    CompElementTemplate<Shape1d> cel(index,&cmesh,&geo);
    
    MatrixDouble ek(2,2),ef(2,1);
    IntRule1d intrule(ordenP);







    cel.CalcStiff(ek, ef);

    std::cout<<"********************Element Matrix********************"<<std::endl;
    cout << ek << endl;
    
}
