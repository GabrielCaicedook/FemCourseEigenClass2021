//
//  TestOneDProblem.cpp
//  FemSC
//
//  Created by Eduardo Ferri on 08/17/15.
//
//
//TestOneDProblem cpp
/*
        Os testes foram preparados com um proposito educacional,
        recomenda-se que o aluno entenda a funcionalidade de cada
        teste e posteriormente use com seu código caso a caso
*/
//      Obs: O xmax e xmin estao tomados como 4 e 0, respectivamente,
//      caso estes valores sejam alterados, editar o teste TestNodes.
//
//
#include <iostream>
#include <math.h>
#include "GeoMesh.h"
#include "ReadGmsh.h"
#include "CompMesh.h"
#include "GeoElement.h"
#include "GeoElementTemplate.h"
#include "MathStatement.h"
#include "L2Projection.h"
#include "Analysis.h"
#include "IntRule.h"
#include "PostProcessTemplate.h"
#include "Poisson.h"
#include "VTKGeoMesh.h"

using std::cout;
using std::endl;
using std::cin;

void exact(const VecDouble &point,VecDouble &val, MatrixDouble &deriv);

int main ()
{
    //Criando da malha geometrica 
    //Basado en OneD.msh
    GeoMesh gmesh;
    ReadGmsh read;
    //Se agrega la localización del archivo OneD.msh 
    std::string filename("oneD1.msh");
#ifdef MACOSX
    filename = "../"+filename;
#endif
    read.Read(gmesh,filename);
    const std::string filenamegeo("GeomeshGabriel.vtk");
    
    VTKGeoMesh::PrintGMeshVTK(&gmesh,filenamegeo);

//Criando de malha computaciononal e de mathstatement de Poisoon com materialid 1

   int orderp = 1;
    CompMesh cmesh(&gmesh);
    MatrixDouble perm(3,3);
    perm.setZero();
    perm(0,0) = 1.;
    perm(1,1) = 1.;
    perm(2,2) = 1.;
    Poisson *mat1 = new Poisson(1/*material id*/,perm);
    mat1->SetDimension(1);
    
    //Cordenada &x y devuelve res
    //Qué es res?
    auto force = [](const VecDouble &x, VecDouble &res)
    {
       //res[0] = x[0]*x[0];
        res[0] = 1.;
    };
    
    mat1->SetForceFunction(force);

    //Por qué los setea en 1

    MatrixDouble proj(1,1),val1(1,1),val2(1,1);
    proj.setZero();
    val1.setZero();
    val2.setZero();

    //Boundary condition = 0

    L2Projection *bc_linha = new L2Projection(0,2/*matid*/,proj,val1,val2);
    L2Projection *bc_point = new L2Projection(0,3/*matid*/,proj,val1,val2);



//Setando todos los mathstatements na malha computacional incluindo
//as condiciones contorno
    std::vector<MathStatement *> mathvec = {0,mat1,bc_point,bc_linha};
    cmesh.SetMathVec(mathvec);
    cmesh.SetDefaultOrder(orderp);
    cmesh.AutoBuild();  
    cmesh.Resequence();

//Creación de analisis e rodando a simulacao que envolve Assemble() e Solve()    
    Analysis AnalysisLoc(&cmesh);
    AnalysisLoc.RunSimulation();
    
    //Postprocesamiento do error baseado na solucao exata
    PostProcessTemplate<Poisson> postprocess;
    postprocess.SetExact(exact);
    
    VecDouble errvec;
    errvec = AnalysisLoc.PostProcessError(std::cout, postprocess);

    postprocess.AppendVariable("Sol");


    const std::string filename2("Solution.vtk");
    AnalysisLoc.PostProcessSolution(filename2, postprocess);
    
    
    return 0;
}

void exact(const VecDouble &point,VecDouble &val, MatrixDouble &deriv){

    deriv(0,0) = 4-point[0];
    val[0]=point[0]*(8.-point[0])/2.;
    return;
}


