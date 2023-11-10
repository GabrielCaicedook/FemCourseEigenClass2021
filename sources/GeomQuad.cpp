/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "GeomQuad.h"

GeomQuad::GeomQuad() {
}

GeomQuad::~GeomQuad() {
}

GeomQuad::GeomQuad(const GeomQuad &copy) {
    fNodeIndices = copy.fNodeIndices;
}

GeomQuad& GeomQuad::operator=(const GeomQuad& copy) {
    fNodeIndices = copy.fNodeIndices;
    return *this;
}

void GeomQuad::Shape(const VecDouble &xi, VecDouble &phi, MatrixDouble &dphi) {
    
    if(xi.size() != Dimension || phi.size() != nCorners || dphi.rows() != Dimension || dphi.cols() != nCorners) DebugStop();
    
    double ksi = xi[0];
    double eta = xi[1];

    
    phi[0]=0.25*(1-ksi)*(1-eta);
    phi[1]=0.25*(1+ksi)*(1-eta);
    phi[2]=0.25*(1+ksi)*(1+eta);
    phi[3]=0.25*(1-ksi)*(1+eta);

    dphi(0,0)=-0.25*(1-eta);
    dphi(0,1)=-0.25*(1+ksi);
    dphi(0,2)=0.25*(1-eta);
    dphi(0,3)=-0.25*(1+ksi);
    dphi(1,0)=0.25*(1+eta);
    dphi(1,1)=0.25*(1+ksi);
    dphi(1,2)=-0.25*(1+eta);
    dphi(1,3)=0.25*(1-ksi);
    std::cout << "pasé por aquí" << std::endl;
}

void GeomQuad::X(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x) {
    
    
    if(xi.size() != Dimension) DebugStop();
    if(x.size() < NodeCo.rows()) DebugStop();
    if(NodeCo.cols() != nCorners) DebugStop();

    int nrow = NodeCo.rows();
    int ncol = NodeCo.cols();
    if (x.size() < nrow) {
        x.resize(2);
    }
   
    x.setZero();
    
    VecDouble phi(nCorners);
    MatrixDouble dphi(Dimension, nCorners);

    Shape(xi, phi, dphi);

    int space = NodeCo.rows();

    for(int i = 0; i < space; i++) {
        x[i] = 0.0;
        for(int j = 0; j < 4; j++) {
            x[i] = phi(j,0)*NodeCo(i,j);
            std::cout << "funcion map " <<x[i]<< std::endl; 
            std::cout << "NodeCo " << NodeCo(i,j)<< std::endl;
            std::cout << "phi " <<phi[j]<< std::endl; 
           
        }
    }

;

    // // Number of coordinates of the resulting x needs to be 1, 2 or 3d according to the NodeCo provided
    // if (x.size() < nrow) x.resize(nrow);
    // x.setZero();

    // // Could have called Shape and used phi, but the mapping is quite
    // // simple in this case
    // for (int i = 0; i < nrow; i++) {
    //     x[i] = 0.25 * NodeCo(i, 0) * (1. - xi[0]) + 0.25 * NodeCo(i, 1) * (1. + xi[0]) + 0.25 * NodeCo(i, 1) * (1. + xi[0])+ 0.25 * NodeCo(i, 1) * (1. + xi[0]);
    //     x[i+1] = 0.25 * NodeCo(i, 1) * (1. - xi[1]) + 0.25 * NodeCo(i, 1) * (1. + xi[1]) + 0.25 * NodeCo(i, 2) * (1. + xi[1])+ 0.25 * NodeCo(i, 2) * (1. + xi[1]);
    //     std::cout <<"NodeCo x "<< NodeCo(i, 0) << " NodeCo y " << NodeCo(i, 1)<< std::endl;
    //     std::cout <<"NodeCo x "<< NodeCo(1,0) << " NodeCo y " << NodeCo(1,1)<< std::endl;
    

    



}

void GeomQuad::GradX(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x, MatrixDouble &gradx) {
    std::cout << "\nPLEASE IMPLEMENT ME\n" << __PRETTY_FUNCTION__ << std::endl;
    if(xi.size() != Dimension) DebugStop();
    if(x.size() != NodeCo.rows()) DebugStop();
    if(NodeCo.cols() != nCorners) DebugStop();
    DebugStop();
}

void GeomQuad::SetNodes(const VecInt &nodes) {
    if(nodes.size() != nCorners) DebugStop();
    fNodeIndices = nodes;
}

void GeomQuad::GetNodes(VecInt &nodes) const{
    nodes = fNodeIndices;
}

int GeomQuad::NodeIndex(int node) const {
    return fNodeIndices[node];
}

int GeomQuad::NumNodes() {
    return nCorners;
}

GeoElementSide GeomQuad::Neighbour(int side) const {
    return fNeighbours[side];
}

void GeomQuad::SetNeighbour(int side, const GeoElementSide &neighbour) {
    fNeighbours[side] = neighbour;
}
