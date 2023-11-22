//
//  TestOneDProblem.cpp MODIFICADO DO ORIGINAL
//  FemSC
//
//  Created by Eduardo Ferri on 08/17/15.
//
//
//TestOneDProblem cpp
/*
 Os testes foram preparados com um proposito educacional,
 recomenda-se que o aluno entenda a funcionalidade de cada
 teste e posteriormente use com seu cÛdigo caso a caso
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
#include "Poisson.h"
#include "L2Projection.h"
#include "Analysis.h"
#include "PostProcessTemplate.h"
 #include "VTKGeoMesh.h"

int main ()
{
    GeoMesh gmesh;
    ReadGmsh read;
    std::string filename("test 222.msh");


#ifdef MACOSX
    filename = "../"+filename;
#endif

    
    read.Read(gmesh,filename);
    const std::string filenamevtk("geomesh.vtk");
    VTKGeoMesh::PrintGMeshVTK(&gmesh, filenamevtk);

    CompMesh cmesh(&gmesh);
    MatrixDouble perm(3,3);
    perm.setZero();
    perm(0,0) = 1.;
    perm(1,1) = 1.;
    perm(2,2) = 1.;
    Poisson *mat1 = new Poisson(1,perm);
    mat1->SetDimension(2);

    auto force = [](const VecDouble &x, VecDouble &res)
    {
        res[0] = 2.*(1.-x[0])*x[0]+2.*(1-x[1])*x[1];
    };
    mat1->SetForceFunction(force);
    MatrixDouble proj(1,1),val1(1,1),val2(1,1);
    proj.setZero();
    val1.setZero();
    val2.setZero();

    int BCD = 0;
    int BCN = 0;
    int matid0 = 2;
    int matid1 = 3;
    int matid2 = 4;
    int matid3 = 5;
    
    L2Projection *bc_linha0 = new L2Projection(BCD,matid0,proj,val1,val2);
    L2Projection *bc_linha1 = new L2Projection(BCD,matid1,proj,val1,val2);
    L2Projection *bc_linha2 = new L2Projection(BCD,matid2,proj,val1,val2);
    L2Projection *bc_linha3 = new L2Projection(BCD,matid3,proj,val1,val2);
   


    std::vector<MathStatement *> mathvec = {0,mat1,bc_linha0,bc_linha1,bc_linha2,bc_linha3};
    cmesh.SetMathVec(mathvec);
    cmesh.SetDefaultOrder(1);
    cmesh.AutoBuild();
    cmesh.Resequence();

    Analysis locAnalysis(&cmesh);
    locAnalysis.RunSimulation();
    PostProcessTemplate<Poisson> postprocess;
    auto exact = [](const VecDouble &x, VecDouble &val, MatrixDouble &deriv)
    {
        val[0] = (1.-x[0])*x[0]*(1-x[1])*x[1];
        deriv(0,0) = (1.-2.*x[0])*(1-x[1])*x[1];
        deriv(1,0) = (1-2.*x[1])*(1-x[0])*x[0];
    };

    mat1->SetExactSolution(exact);
    const std::string filenameSol("solutionQuad3.vtk");
    const std::string namevar("Sol");
    const std::string namevar2("SolExact");

    postprocess.AppendVariable(namevar);
    postprocess.AppendVariable(namevar2);
    locAnalysis.PostProcessSolution(filenameSol, postprocess);

//    if (!strcmp("Sol", name.c_str())) return ESol;
//    if (!strcmp("DSol", name.c_str())) return EDSol;
//    if (!strcmp("Flux", name.c_str())) return EFlux;
//    if (!strcmp("Force", name.c_str())) return EForce;
//    if (!strcmp("SolExact", name.c_str())) return ESolExact;
//    if (!strcmp("DSolExact", name.c_str())) return EDSolExact;


    postprocess.AppendVariable("Sol");
    postprocess.AppendVariable("SolExact");
    
    // postprocess.AppendVariable("DSol");
    // postprocess.AppendVariable("Flux");
    // postprocess.AppendVariable("Force");
    
    // postprocess.AppendVariable("DSolExact");
    postprocess.SetExact(exact);
    mat1->SetExactSolution(exact);
    locAnalysis.PostProcessSolution("triangle.vtk", postprocess);

    //locAnalysis.PostProcessSolution("quads.vtk", postprocess);


    VecDouble errvec;
    errvec = locAnalysis.PostProcessError(std::cout, postprocess);
    
    return 0;
}
