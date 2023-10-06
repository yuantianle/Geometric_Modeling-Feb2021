//**************************************************************//
/*  The algorithm of the codes is refered from:                 */
/*  https://www.cnblogs.com/TooManyWayToBe/default.html?page=2  */
/*  Edited by Tianle Yuan, Oregon State University for CS554    */
/*  Last updated date: 3/17/2021                                */
//**************************************************************//


#pragma once

#include "glm/glm/glm.hpp"
#include "glm/glm/gtc/matrix_transform.hpp"
#include "glm/glm/gtc/type_ptr.hpp"
#include <vector>

#include <iostream>
#include "learnply.h"



void Polyhedron::InitAxesAndMatrixSelf()
{
    //Firstly construct Local Frame
    for (int i = 0; i < nverts; i++)
    {
        glm::vec3 vertexNormal, E1, E2;

        Corner* c;
        for (int k = 0; k < 3; k++)
        {
            if (vlist[i] == vlist[i]->tris[0]->corners[k]->v)
            {
                c = vlist[i]->tris[0]->corners[k];   // record the first triangle's corner of Vi
                break;
            }
        }

        for (int j = 0; j < vlist[i]->ntris; j++)
        {
            glm::vec3 tj_v0 = glm::vec3(c->v->x, c->v->y, c->v->z);
            glm::vec3 tj_v1 = glm::vec3(c->n->v->x, c->n->v->y, c->n->v->z);
            glm::vec3 tj_v2 = glm::vec3(c->n->n->v->x, c->n->n->v->y, c->n->n->v->z);

            glm::vec3 BA = tj_v0 - tj_v1;
            glm::vec3 BC = tj_v2 - tj_v1;
            glm::vec3 normal = cross(BA,BC);
            normal = normalize(normal);
            vertexNormal += normal;
            
            c = c->p->o->p;
        }
        vertexNormal /= vlist[i]->ntris;
        vertexNormal = normalize(vertexNormal);
        if (vertexNormal != glm::vec3(0.,0.,1.) && vertexNormal != glm::vec3(0., 0., -1.))
        {
            E1 = glm::vec3(-vertexNormal.y, vertexNormal.x, 0.0);
            E2 = cross(vertexNormal, E1);
        }
        else if (vertexNormal == glm::vec3(0., 0., 1.))
        {
            E1 = glm::vec3(1, 0, 0);
            E2 = glm::vec3(0, 1, 0);
        }
        else
        {
            E1 = glm::vec3(0, 1, 0);
            E2 = glm::vec3(1, 0, 0);
        }
        
        E1 = normalize(E1);
        E2 = normalize(E2);

        vlist[i]->E1 = E1;
        vlist[i]->E2 = E2;
        vlist[i]->N = vertexNormal;
        vlist[i]->axesSelf.resize(3, 3);

        std::vector<Eigen::Triplet<float>> vectorTriplet;
        vectorTriplet.push_back(Eigen::Triplet<float>(0, 0, E1.x));
        vectorTriplet.push_back(Eigen::Triplet<float>(1, 0, E1.y));
        vectorTriplet.push_back(Eigen::Triplet<float>(2, 0, E1.z));

        vectorTriplet.push_back(Eigen::Triplet<float>(0, 1, E2.x));
        vectorTriplet.push_back(Eigen::Triplet<float>(1, 1, E2.y));
        vectorTriplet.push_back(Eigen::Triplet<float>(2, 1, E2.z));

        vectorTriplet.push_back(Eigen::Triplet<float>(0, 2, vertexNormal.x));
        vectorTriplet.push_back(Eigen::Triplet<float>(1, 2, vertexNormal.y));
        vectorTriplet.push_back(Eigen::Triplet<float>(2, 2, vertexNormal.z));

        Eigen::SparseMatrix<float> matrixTemp(3, 3);
        matrixTemp.setFromTriplets(vectorTriplet.begin(), vectorTriplet.end());
        vlist[i]->axesSelf = matrixTemp;

        // Construct the Laplacian coordinate length coeffient 
        // of the relative local coordinate system.
        vlist[i]->relativeX = dot(vlist[i]->lplsPos, E1);
        vlist[i]->relativeY = dot(vlist[i]->lplsPos, E2);
        vlist[i]->relativeZ = dot(vlist[i]->lplsPos, vertexNormal);
    }

    //Construct the transformation matrix from the surrounding first-order vertex to this vertex
    for (int i = 0; i < nverts; i++)
    {
        Corner* c;
        for (int k = 0; k < 3; k++)
        {
            if (vlist[i] == vlist[i]->tris[0]->corners[k]->v)
            {
                c = vlist[i]->tris[0]->corners[k];   // record the first triangle's corner of Vi
                break;
            }
        }

        Eigen::Matrix3f axesSelf = vlist[i]->axesSelf;
        //Eigen::Matrix3f matrixSelf = axesSelf;
       
        for (int j = 0; j < vlist[i]->ntris; j++)
        {
          
            Eigen::Matrix3f axesOther = c->n->v->axesSelf;
            //std::cout << axesOther << std::endl;
            Eigen::Matrix3f axesOtherInverse = axesOther.inverse();
            //std::cout << axesOtherInverse << std::endl;
            // Ax = b; T * axesOther = axesSelf;
            Eigen::Matrix3f  matrixT = axesSelf * axesOtherInverse;
            vlist[i]->mapMatrixOtherToSelf.insert(std::pair<int, Eigen::Matrix3f>(c->n->v->index, matrixT));
            //test
            Eigen::Matrix3f matrixTemp = matrixT * axesOther;
            //std::cout << matrixTemp << std::endl;
            //std::cout << axesSelf << std::endl;

            c = c->p->o->p;
        }
    }
}


void Polyhedron::InitNewAxesMatrix()
{
    int sizeOfVertex = nverts;
    std::vector<Eigen::Triplet<float> > vectorTriplet;
    for (int i = 0; i < sizeOfVertex; i++)
    {
        
        int vertexID = vlist[i]->index;
        float sizeOfNeighborVertex = vlist[i]->ntris;

        vectorTriplet.push_back(Eigen::Triplet<float>(vertexID * 3, vertexID * 3, -sizeOfNeighborVertex));
        vectorTriplet.push_back(Eigen::Triplet<float>(vertexID * 3 + 1, vertexID * 3 + 1, -sizeOfNeighborVertex));
        vectorTriplet.push_back(Eigen::Triplet<float>(vertexID * 3 + 2, vertexID * 3 + 2, -sizeOfNeighborVertex));

        std::map<int, Eigen::Matrix3f>::iterator itorOfT;
        for (itorOfT = vlist[i]->mapMatrixOtherToSelf.begin(); itorOfT != vlist[i]->mapMatrixOtherToSelf.end(); itorOfT++)
        {
            int vertexNeighborID = itorOfT->first;
            Eigen::Matrix3f T = itorOfT->second;
            vectorTriplet.push_back(Eigen::Triplet<float>(vertexID * 3, vertexNeighborID * 3, T(0, 0)));
            vectorTriplet.push_back(Eigen::Triplet<float>(vertexID * 3, vertexNeighborID * 3 + 1, T(0, 1)));
            vectorTriplet.push_back(Eigen::Triplet<float>(vertexID * 3, vertexNeighborID * 3 + 2, T(0, 2)));

            vectorTriplet.push_back(Eigen::Triplet<float>(vertexID * 3 + 1, vertexNeighborID * 3, T(1, 0)));
            vectorTriplet.push_back(Eigen::Triplet<float>(vertexID * 3 + 1, vertexNeighborID * 3 + 1, T(1, 1)));
            vectorTriplet.push_back(Eigen::Triplet<float>(vertexID * 3 + 1, vertexNeighborID * 3 + 2, T(1, 2)));

            vectorTriplet.push_back(Eigen::Triplet<float>(vertexID * 3 + 2, vertexNeighborID * 3, T(2, 0)));
            vectorTriplet.push_back(Eigen::Triplet<float>(vertexID * 3 + 2, vertexNeighborID * 3 + 1, T(2, 1)));
            vectorTriplet.push_back(Eigen::Triplet<float>(vertexID * 3 + 2, vertexNeighborID * 3 + 2, T(2, 2)));
        }
    }
    //Add constrained vertices; 
    // the first circle and the last circle can be adjusted by themselves
    for (int i = 0; i < (m_iAroundNumber + m_iHeightNumber); i++)
    {
        int vertexID;
        float value = 1.0;
        if (i < m_iAroundNumber)
        {
            int p =  m_fixedindex[i];
            vertexID = vlist[p]->index;
        }
        else
        {
            //int q = i + m_iAroundNumber * (m_iHeightNumber - 2);
            int q = m_editedintex[i - m_iAroundNumber];
            vertexID = vlist[q]->index;
        }
        vectorTriplet.push_back(Eigen::Triplet<float>(sizeOfVertex * 3 + i * 3, vertexID * 3, value));
        vectorTriplet.push_back(Eigen::Triplet<float>(sizeOfVertex * 3 + i * 3 + 1, vertexID * 3 + 1, value));
        vectorTriplet.push_back(Eigen::Triplet<float>(sizeOfVertex * 3 + i * 3 + 2, vertexID * 3 + 2, value));
    }

    m_NewAxesMatrix.resize(sizeOfVertex * 3 + ((m_iAroundNumber + m_iHeightNumber) * 3), sizeOfVertex * 3);
    m_NewAxesMatrix.setFromTriplets(vectorTriplet.begin(), vectorTriplet.end());
}

void Polyhedron::CalculateNewAxes()
{
    Eigen::SparseMatrix<float> matrixNewAxesTrans = m_NewAxesMatrix.transpose();
    Eigen::SparseMatrix<float> matrixRes = matrixNewAxesTrans * m_NewAxesMatrix;
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<float> > matrixCholesky(matrixRes);

    std::vector<Eigen::Triplet<float> > vectorTriplet;
    int sizeOfVertex = nverts;

    //The right side of the equation system (b vector)
    Eigen::VectorXf vectorX, vectorY, vectorZ;
    vectorX.resize(sizeOfVertex * 3 + (m_iAroundNumber + m_iHeightNumber) * 3);
    vectorX.setZero();
    vectorY.resize(sizeOfVertex * 3 + (m_iAroundNumber + m_iHeightNumber) * 3);
    vectorY.setZero();
    vectorZ.resize(sizeOfVertex * 3 + (m_iAroundNumber + m_iHeightNumber) * 3);
    vectorZ.setZero();
    for (int i = 0; i < (m_iAroundNumber + m_iHeightNumber); i++)
    {
        //int vertexID;
        Eigen::Matrix3f matrixSelf;
        float value = 1.0;
        if (i < m_iAroundNumber)
        {
            //vertexID = vlist[q]->index;

            //int q = i + m_iAroundNumber * (m_iHeightNumber - 2);
            int q = m_fixedindex[i];
            matrixSelf = vlist[q]->axesSelf;
        }
        else
        {   
            //vertexID = vlist[i]->index;
            int p = m_editedintex[i - m_iAroundNumber];

            Eigen::Matrix3f rotate;
            //m_rotate_theta = -PI / 6.;  //for Z
            //m_rotate_theta = PI / 6.;
            //m_rotate_theta = -PI / 2.;     //for X
            m_rotate_theta = PI / 2.;     //for X
            //rotate << 1, 0, 0,
            //    0, 1, 0,
            //    0, 0, 1;
            //rotate << 0.5, 0, 0,
            //    0, 0.5, 0,
            //    0, 0, 0.5;

            rotate << 1, 0, 0,
                0, cos(m_rotate_theta), -sin(m_rotate_theta),
                0, sin(m_rotate_theta), cos(m_rotate_theta);

            //rotate << cos(m_rotate_theta), 0, sin(m_rotate_theta),
            //    0, 1, 0,
            //    -sin(m_rotate_theta), 0, cos(m_rotate_theta);

            //rotate << cos(m_rotate_theta), -sin(m_rotate_theta), 0,
            //          sin(m_rotate_theta), cos(m_rotate_theta),  0,
            //          0,          0,           1;
            //matrixSelf = rotate * vlist[p]->axesSelf;
        }

        if (i == m_iAroundNumber + m_iHeightNumber - 1)
        {
            Eigen::Matrix3f rotate2;
            rotate2 << 1, 0, 0,
                0, cos(m_rotate_theta), -sin(m_rotate_theta),
                0, sin(m_rotate_theta), cos(m_rotate_theta);
            matrixSelf *= rotate2;
        }
        

        vectorX[sizeOfVertex * 3 + i * 3] = matrixSelf(0, 0);
        vectorX[sizeOfVertex * 3 + i * 3 + 1] = matrixSelf(1, 0);
        vectorX[sizeOfVertex * 3 + i * 3 + 2] = matrixSelf(2, 0);

        vectorY[sizeOfVertex * 3 + i * 3] = matrixSelf(0, 1);
        vectorY[sizeOfVertex * 3 + i * 3 + 1] = matrixSelf(1, 1);
        vectorY[sizeOfVertex * 3 + i * 3 + 2] = matrixSelf(2, 1);

        vectorZ[sizeOfVertex * 3 + i * 3] = matrixSelf(0, 2);
        vectorZ[sizeOfVertex * 3 + i * 3 + 1] = matrixSelf(1, 2);
        vectorZ[sizeOfVertex * 3 + i * 3 + 2] = matrixSelf(2, 2);
    }

    vectorX = matrixNewAxesTrans * vectorX;
    vectorY = matrixNewAxesTrans * vectorY;
    vectorZ = matrixNewAxesTrans * vectorZ;

    Eigen::VectorXf vectorResX, vectorResY, vectorResZ;
    vectorResX = matrixCholesky.solve(vectorX);
    vectorResY = matrixCholesky.solve(vectorY);
    vectorResZ = matrixCholesky.solve(vectorZ);

    //Get the new frame
    for (int i = 0; i < nverts; i++)
    {
        Eigen::Matrix3f matrixSelf;
        float x1, x2, x3, y1, y2, y3, z1, z2, z3;
        x1 = vectorResX[i * 3];
        x2 = vectorResX[i * 3 + 1];
        x3 = vectorResX[i * 3 + 2];

        y1 = vectorResY[i * 3];
        y2 = vectorResY[i * 3 + 1];
        y3 = vectorResY[i * 3 + 2];

        z1 = vectorResZ[i * 3];
        z2 = vectorResZ[i * 3 + 1];
        z3 = vectorResZ[i * 3 + 2];

        vlist[i]->E1 = glm::vec3(x1, x2, x3);
        vlist[i]->E2 = glm::vec3(y1, y2, y3);
        vlist[i]->N =  glm::vec3(z1, z2, z3);

        //Get the new laplasian coordinates
        vlist[i]->lplsPos = (vlist[i]->E1 * vlist[i]->relativeX) + (vlist[i]->E2 * vlist[i]->relativeY) + (vlist[i]->N * vlist[i]->relativeZ);
    }
}

void Polyhedron::InitLpls()
{
    //Method One: Calculate the weight of the vertex first order adjacent to the vertex.
    for (int i = 0; i < nverts; i++)
    {
        Corner* c;
        for (int k = 0; k < 3; k++)
        {
            if (vlist[i] == vlist[i]->tris[0]->corners[k]->v)
            {
                c = vlist[i]->tris[0]->corners[k];   // record the first triangle's corner of Vi
                break;
            }
        }

        for (int j = 0; j < vlist[i]->ntris; j++)
        {

            float weight = float(1) / float(vlist[i]->ntris);  // Uniform method
            vlist[i]->mapWeightForOther.insert(std::pair<int, float>(c->n->v->index, weight));
            vlist[i]->fTotalWeight = 1.0;

            c = c->p->o->p;
        }
    }

    /*
    //Method Two: Construct the weights of the Laplace matrix.
    for (int i = 0; i < vertexNumber; i++)
    {
        pVERTEX pVertex = m_vecAllVertex.at(i);
        float totalWeight = 0.0;
        for (int j = 0; j < pVertex->vecNeighborEdge.size(); j++)
        {
            pEDGE pEdge = pVertex->vecNeighborEdge.at(j);
            pVERTEX pA = pEdge->pA;
            pVERTEX pB = pEdge->pB;

            pVERTEX pTarget;
            if (pA->id == pVertex->id)
                pTarget = pB;
            else
                pTarget = pA;

            std::vector<pTRIANGLE> vecTri = pEdge->vecNeighborTri;
            pVERTEX pC = NULL;
            float weight = 0.0;
            for (int k = 0; k < vecTri.size(); k++)
            {
                pTRIANGLE pTri = vecTri.at(k);
                pVERTEX p1 = pTri->pA;
                pVERTEX p2 = pTri->pB;
                pVERTEX p3 = pTri->pC;
                if ((pA->id == p1->id && pB->id == p2->id) || (pA->id == p2->id && pB->id == p1->id))
                {
                    pC = p3;
                }
                else if ((pA->id == p1->id && pB->id == p3->id) || (pA->id == p3->id && pB->id == p1->id))
                {
                    pC = p2;
                }
                else if ((pA->id == p2->id && pB->id == p3->id) || (pA->id == p3->id && pB->id == p2->id))
                {
                    pC = p1;
                }
                //Start to get the weight value
                float cotAngle = 0.0;
                GetCotAngle(pA->pos, pB->pos, pC->pos, cotAngle);
                weight += cotAngle;
            }
            weight /= 2.0;
            pVertex->mapWeightForOther.insert(std::pair<int, float>(pTarget->id, weight) );
            totalWeight += weight;
        }
        pVertex->fTotalWeight = totalWeight;
    }
    */
    
    //Calculate the Laplasian coordinates
    for (int i = 0; i < nverts; i++)
    {
        Corner* c;
        for (int k = 0; k < 3; k++)
        {
            if (vlist[i] == vlist[i]->tris[0]->corners[k]->v)
            {
                c = vlist[i]->tris[0]->corners[k];   // record the first triangle's corner of Vi
                break;
            }
        }

        glm::vec3 pos = glm::vec3(vlist[i]->x, vlist[i]->y, vlist[i]->z);
        pos = pos * vlist[i]->fTotalWeight;
        glm::vec3 otherPos = glm::vec3(0.0, 0.0, 0.0);
        for (int j = 0; j < vlist[i]->ntris; j++)
        {
            
            std::map<int, float>::iterator itor = vlist[i]->mapWeightForOther.find(c->n->v->index);
            float weight = itor->second;
            otherPos += glm::vec3(c->n->v->x, c->n->v->y, c->n->v->z) * weight;

            c = c->p->o->p;
        }
        vlist[i]->lplsPos = pos - otherPos;
    }

    int count = 0;
    std::vector<int> beginNumber(nverts);
    for (int i = 0; i < nverts; i++)
    {
        beginNumber[i] = count;
        count += vlist[i]->ntris + 1;
    }

    std::vector<Eigen::Triplet<float> > vectorTriplet(count);
    for (int i = 0; i < nverts; i++)
    {
        //Get the original Laplasian matrix
        vectorTriplet[beginNumber[i]] = Eigen::Triplet<float>(i, i, vlist[i]->fTotalWeight);
        int j = 0;
        std::map<int, float>::iterator itor;
        for (itor = vlist[i]->mapWeightForOther.begin(); itor != vlist[i]->mapWeightForOther.end(); itor++)
        {
            int neighborID = itor->first;
            float weight = itor->second;
            vectorTriplet[beginNumber[i] + j + 1] = Eigen::Triplet<float>(i, neighborID, -weight);
            j++;
        }
    }

    //Top a circle, Bottom a circle
    for (int i = 0; i < (m_iAroundNumber + m_iHeightNumber); i++)
    {
        float weight = 1.0;

        if (i < m_iAroundNumber)
        {      
            int p = m_fixedindex[i];
            int row = i + nverts;
            vectorTriplet.push_back(Eigen::Triplet<float>(row, vlist[p]->index, weight));
        }
        else
        {
            //int q = i + m_iAroundNumber * 1; //14
            int q = m_editedintex[i- m_iAroundNumber];
            int row = i + nverts;
            vectorTriplet.push_back(Eigen::Triplet<float>(row, vlist[q]->index, weight));
        }


    }

    m_Matrix.resize(nverts + (m_iAroundNumber + m_iHeightNumber), nverts);
    //m_Matrix.resize(vertexNumber, vertexNumber);
    m_Matrix.setFromTriplets(vectorTriplet.begin(), vectorTriplet.end());
    //std::cout << m_Matrix << std::endl;
}

void Polyhedron::CalculateNewMatrixLpls()
{
    Eigen::SparseMatrix<float> Matrix = m_Matrix;

    Eigen::SparseMatrix<float> MatrixTranspose = Matrix.transpose();
    Eigen::SparseMatrix<float> MatrixLpls = MatrixTranspose * Matrix;

    Eigen::VectorXf targetX, targetY, targetZ;
    int vertexNumber = nverts;

    targetX.resize(vertexNumber + (m_iAroundNumber + m_iHeightNumber));
    targetY.resize(vertexNumber + (m_iAroundNumber + m_iHeightNumber));
    targetZ.resize(vertexNumber + (m_iAroundNumber + m_iHeightNumber));
    
    //Right side of the euqation system
    for (int i = 0; i < vertexNumber + (m_iAroundNumber + m_iHeightNumber); i++)
    {
        if (i < vertexNumber)
        {
            targetX[i] = vlist[i]->lplsPos[0];
            targetY[i] = vlist[i]->lplsPos[1];
            targetZ[i] = vlist[i]->lplsPos[2];
        }
        else if (i < vertexNumber + m_iAroundNumber)
        {
            //The First level's vertices coordinates
            //targetX[i] = vlist[i - vertexNumber]->x + 0;
            //targetY[i] = vlist[i - vertexNumber]->y + 0;
            //targetZ[i] = vlist[i - vertexNumber]->z + 0.2;
            int p = m_fixedindex[i-vertexNumber];
            targetX[i] = vlist[p]->x ;
            targetY[i] = vlist[p]->y ;
            targetZ[i] = vlist[p]->z ;
        }

        //------------------------------------Z negative PI/6------------------------------------
        //else if (i == vertexNumber + m_iAroundNumber)
        //{
        //    //The Topest level's vertices coordinates
        //    //int q = i - (vertexNumber + m_iAroundNumber) + (vertexNumber - m_iAroundNumber);
        //    int q = m_editedintex[i - vertexNumber - m_iAroundNumber];
        //    targetX[i] = vlist[q]->x - 0.1;//- 0.75;//;
        //    targetY[i] = vlist[q]->y - 0.3;//-0.5;//+ sin(m_rotate_theta);//- 1.2;//;
        //    targetZ[i] = vlist[q]->z ;//+0.5 ;//+ cos(m_rotate_theta);//+ 0.;//;
        //}
        //else if (i == vertexNumber + m_iAroundNumber+1)
        //{
        //    //The Topest level's vertices coordinates
        //    //int q = i - (vertexNumber + m_iAroundNumber) + (vertexNumber - m_iAroundNumber);
        //    int q = m_editedintex[i - vertexNumber - m_iAroundNumber];
        //    targetX[i] = vlist[q]->x - 0.1;//- 0.75;//;
        //    targetY[i] = vlist[q]->y - 0.3;//-0.5;//+ sin(m_rotate_theta);//- 1.2;//;
        //    targetZ[i] = vlist[q]->z;//+0.5 ;//+ cos(m_rotate_theta);//+ 0.;//;
        //}
        //else if (i == vertexNumber + m_iAroundNumber+2)
        //{
        //    //The Topest level's vertices coordinates
        //    //int q = i - (vertexNumber + m_iAroundNumber) + (vertexNumber - m_iAroundNumber);
        //    int q = m_editedintex[i - vertexNumber - m_iAroundNumber];
        //    targetX[i] = vlist[q]->x - 0.5;//- 0.75;//;
        //    targetY[i] = vlist[q]->y - 0.2;//-0.5;//+ sin(m_rotate_theta);//- 1.2;//;
        //    targetZ[i] = vlist[q]->z;//+0.5 ;//+ cos(m_rotate_theta);//+ 0.;//;
        //}
        //else
        //{
        //    //The Topest level's vertices coordinates
        //    //int q = i - (vertexNumber + m_iAroundNumber) + (vertexNumber - m_iAroundNumber);
        //    int q = m_editedintex[i - vertexNumber - m_iAroundNumber];
        //    targetX[i] = vlist[q]->x - 0.5;//- 0.75;//;
        //    targetY[i] = vlist[q]->y - 0.2;//-0.5;//+ sin(m_rotate_theta);//- 1.2;//;
        //    targetZ[i] = vlist[q]->z;//+0.5 ;//+ cos(m_rotate_theta);//+ 0.;//;
        //}
        //------------------------------------Z positive PI/6------------------------------------
        //else if (i == vertexNumber + m_iAroundNumber)
        //{
        ////The Topest level's vertices coordinates
        ////int q = i - (vertexNumber + m_iAroundNumber) + (vertexNumber - m_iAroundNumber);
        //int q = m_editedintex[i - vertexNumber - m_iAroundNumber];
        //targetX[i] = vlist[q]->x - 0.5;//- 0.75;//;
        //targetY[i] = vlist[q]->y + 0.2;//-0.5;//+ sin(m_rotate_theta);//- 1.2;//;
        //targetZ[i] = vlist[q]->z;//+0.5 ;//+ cos(m_rotate_theta);//+ 0.;//;
        //}
        //else if (i == vertexNumber + m_iAroundNumber + 1)
        //{
        ////The Topest level's vertices coordinates
        ////int q = i - (vertexNumber + m_iAroundNumber) + (vertexNumber - m_iAroundNumber);
        //int q = m_editedintex[i - vertexNumber - m_iAroundNumber];
        //targetX[i] = vlist[q]->x - 0.5;//- 0.75;//;
        //targetY[i] = vlist[q]->y + 0.2;//-0.5;//+ sin(m_rotate_theta);//- 1.2;//;
        //targetZ[i] = vlist[q]->z;//+0.5 ;//+ cos(m_rotate_theta);//+ 0.;//;
        //}
        //else if (i == vertexNumber + m_iAroundNumber + 2)
        //{
        ////The Topest level's vertices coordinates
        ////int q = i - (vertexNumber + m_iAroundNumber) + (vertexNumber - m_iAroundNumber);
        //int q = m_editedintex[i - vertexNumber - m_iAroundNumber];
        //targetX[i] = vlist[q]->x - 0.1;//- 0.75;//;
        //targetY[i] = vlist[q]->y + 0.3;//-0.5;//+ sin(m_rotate_theta);//- 1.2;//;
        //targetZ[i] = vlist[q]->z;//+0.5 ;//+ cos(m_rotate_theta);//+ 0.;//;
        //}
        //else
        //{
        ////The Topest level's vertices coordinates
        ////int q = i - (vertexNumber + m_iAroundNumber) + (vertexNumber - m_iAroundNumber);
        //int q = m_editedintex[i - vertexNumber - m_iAroundNumber];
        //targetX[i] = vlist[q]->x - 0.1;//- 0.75;//;
        //targetY[i] = vlist[q]->y + 0.3;//-0.5;//+ sin(m_rotate_theta);//- 1.2;//;
        //targetZ[i] = vlist[q]->z;//+0.5 ;//+ cos(m_rotate_theta);//+ 0.;//;
        //}
        //------------------------------------X negative PI/2------------------------------------
        //else if (i == vertexNumber + m_iAroundNumber)
        //{
        ////The Topest level's vertices coordinates
        ////int q = i - (vertexNumber + m_iAroundNumber) + (vertexNumber - m_iAroundNumber);
        //int q = m_editedintex[i - vertexNumber - m_iAroundNumber];
        //targetX[i] = vlist[q]->x - 0.4;//- 0.75;//;
        //targetY[i] = vlist[q]->y;//-0.5;//+ sin(m_rotate_theta);//- 1.2;//;
        //targetZ[i] = vlist[q]->z - 1;//+0.5 ;//+ cos(m_rotate_theta);//+ 0.;//;
        //}
        //else if (i == vertexNumber + m_iAroundNumber + 1)
        //{
        ////The Topest level's vertices coordinates
        ////int q = i - (vertexNumber + m_iAroundNumber) + (vertexNumber - m_iAroundNumber);
        //int q = m_editedintex[i - vertexNumber - m_iAroundNumber];
        //targetX[i] = vlist[q]->x - 0.4;//- 0.75;//;
        //targetY[i] = vlist[q]->y - 1;//-0.5;//+ sin(m_rotate_theta);//- 1.2;//;
        //targetZ[i] = vlist[q]->z;//+0.5 ;//+ cos(m_rotate_theta);//+ 0.;//;
        //}
        //else if (i == vertexNumber + m_iAroundNumber + 2)
        //{
        ////The Topest level's vertices coordinates
        ////int q = i - (vertexNumber + m_iAroundNumber) + (vertexNumber - m_iAroundNumber);
        //int q = m_editedintex[i - vertexNumber - m_iAroundNumber];
        //targetX[i] = vlist[q]->x - 0.4;//- 0.75;//;
        //targetY[i] = vlist[q]->y + 1;//-0.5;//+ sin(m_rotate_theta);//- 1.2;//;
        //targetZ[i] = vlist[q]->z;//+0.5 ;//+ cos(m_rotate_theta);//+ 0.;//;
        //}
        //else
        //{
        ////The Topest level's vertices coordinates
        ////int q = i - (vertexNumber + m_iAroundNumber) + (vertexNumber - m_iAroundNumber);
        //int q = m_editedintex[i - vertexNumber - m_iAroundNumber];
        //targetX[i] = vlist[q]->x - 0.4;//- 0.75;//;
        //targetY[i] = vlist[q]->y;//-0.5;//+ sin(m_rotate_theta);//- 1.2;//;
        //targetZ[i] = vlist[q]->z + 1;//+0.5 ;//+ cos(m_rotate_theta);//+ 0.;//;
        //}
        //------------------------------------X positive PI/2------------------------------------
        else if (i == vertexNumber + m_iAroundNumber)
        {
        //The Topest level's vertices coordinates
        //int q = i - (vertexNumber + m_iAroundNumber) + (vertexNumber - m_iAroundNumber);
        int q = m_editedintex[i - vertexNumber - m_iAroundNumber];
        targetX[i] = vlist[q]->x - 0.4;//- 0.75;//;
        targetY[i] = vlist[q]->y - 1;//-0.5;//+ sin(m_rotate_theta);//- 1.2;//;
        targetZ[i] = vlist[q]->z;//+0.5 ;//+ cos(m_rotate_theta);//+ 0.;//;
        }
        else if (i == vertexNumber + m_iAroundNumber + 1)
        {
        //The Topest level's vertices coordinates
        //int q = i - (vertexNumber + m_iAroundNumber) + (vertexNumber - m_iAroundNumber);
        int q = m_editedintex[i - vertexNumber - m_iAroundNumber];
        targetX[i] = vlist[q]->x - 0.4;//- 0.75;//;
        targetY[i] = vlist[q]->y;//-0.5;//+ sin(m_rotate_theta);//- 1.2;//;
        targetZ[i] = vlist[q]->z + 1;//+0.5 ;//+ cos(m_rotate_theta);//+ 0.;//;
        }
        else if (i == vertexNumber + m_iAroundNumber + 2)
        {
        //The Topest level's vertices coordinates
        //int q = i - (vertexNumber + m_iAroundNumber) + (vertexNumber - m_iAroundNumber);
        int q = m_editedintex[i - vertexNumber - m_iAroundNumber];
        targetX[i] = vlist[q]->x - 0.4;//- 0.75;//;
        targetY[i] = vlist[q]->y;//-0.5;//+ sin(m_rotate_theta);//- 1.2;//;
        targetZ[i] = vlist[q]->z - 1;//+0.5 ;//+ cos(m_rotate_theta);//+ 0.;//;
        }
        else
        {
        //The Topest level's vertices coordinates
        //int q = i - (vertexNumber + m_iAroundNumber) + (vertexNumber - m_iAroundNumber);
        int q = m_editedintex[i - vertexNumber - m_iAroundNumber];
        targetX[i] = vlist[q]->x - 0.4;//- 0.75;//;
        targetY[i] = vlist[q]->y + 1;//-0.5;//+ sin(m_rotate_theta);//- 1.2;//;
        targetZ[i] = vlist[q]->z;//+0.5 ;//+ cos(m_rotate_theta);//+ 0.;//;
        }
    }

    Eigen::SimplicialCholesky<Eigen::SparseMatrix<float> > MatrixCholesky(MatrixLpls);
    Eigen::VectorXf resX, resY, resZ;
    resX = MatrixTranspose * targetX;
    resY = MatrixTranspose * targetY;
    resZ = MatrixTranspose * targetZ;

    Eigen::VectorXf XX, YY, ZZ;
    XX = MatrixCholesky.solve(resX);
    //std::cout << X << std::endl;
    YY = MatrixCholesky.solve(resY);
    //std::cout << Y << std::endl;
    ZZ = MatrixCholesky.solve(resZ);
    //std::cout << Z << std::endl;


    for (int i = 0; i < nverts; i++)
    {
        int apointflag = false;
        for (int p = 0; p < m_fixedindex.size(); p++)
        {
            if (i == m_fixedindex[p])
            {
                apointflag = true;
                break;
            }
        }
        for (int p = 0; p < m_editedintex.size(); p++)
        {
            if (i == m_editedintex[p])
            {
                apointflag = true;
                break;
            }
        }

        if (apointflag == false)
        {
            vlist[i]->x = XX[i];
            vlist[i]->y = YY[i];
            vlist[i]->z = ZZ[i];
        }

        //for (int p = 0; p < m_editedintex.size(); p++)
        //{
        //    if (i == m_editedintex[p])
        //    {
        //        vlist[i]->x;//-= 0.75;//;
        //        vlist[i]->y ;//-= 0.5;//+= sin(m_rotate_theta);//-= 1.2;//;
        //        vlist[i]->z ;//+= 0.5;//+= cos(m_rotate_theta);//+= 0.;//;  
        //    }
        //}

        //------------------------------------Z negative PI/6------------------------------------
        //if (i == m_editedintex[0])
        //{
        //    vlist[i]->x -= 0.1;//-= 0.75;//;
        //    vlist[i]->y -= 0.3;//-= 0.5;//+= sin(m_rotate_theta);//-= 1.2;//;
        //    vlist[i]->z ;//+= 0.5;//+= cos(m_rotate_theta);//+= 0.;//;  
        //}
        //if (i == m_editedintex[1])
        //{
        //    vlist[i]->x -= 0.1;//-= 0.75;//;
        //    vlist[i]->y -= 0.3;//-= 0.5;//+= sin(m_rotate_theta);//-= 1.2;//;
        //    vlist[i]->z;//+= 0.5;//+= cos(m_rotate_theta);//+= 0.;//;  
        //}
        //if (i == m_editedintex[2])
        //{
        //    vlist[i]->x -= 0.5;//-= 0.75;//;
        //    vlist[i]->y -= 0.2;//+= sin(m_rotate_theta);//-= 1.2;//;
        //    vlist[i]->z;//+= 0.5;//+= cos(m_rotate_theta);//+= 0.;//;  
        //}
        //if (i == m_editedintex[3])
        //{
        //    vlist[i]->x -= 0.5;//-= 0.75;//;
        //    vlist[i]->y -= 0.2;//+= sin(m_rotate_theta);//-= 1.2;//;
        //    vlist[i]->z;//+= 0.5;//+= cos(m_rotate_theta);//+= 0.;//;  
        //}
        //------------------------------------Z positive PI/6------------------------------------
        //if (i == m_editedintex[0])
        //{
        //    vlist[i]->x -= 0.5;//-= 0.75;//;
        //    vlist[i]->y += 0.2;//-= 0.5;//+= sin(m_rotate_theta);//-= 1.2;//;
        //    vlist[i]->z;//+= 0.5;//+= cos(m_rotate_theta);//+= 0.;//;  
        //}
        //if (i == m_editedintex[1])
        //{
        //    vlist[i]->x -= 0.5;//-= 0.75;//;
        //    vlist[i]->y += 0.2;//-= 0.5;//+= sin(m_rotate_theta);//-= 1.2;//;
        //    vlist[i]->z;//+= 0.5;//+= cos(m_rotate_theta);//+= 0.;//;  
        //}
        //if (i == m_editedintex[2])
        //{
        //    vlist[i]->x -= 0.1;//-= 0.75;//;
        //    vlist[i]->y += 0.3;//+= sin(m_rotate_theta);//-= 1.2;//;
        //    vlist[i]->z;//+= 0.5;//+= cos(m_rotate_theta);//+= 0.;//;  
        //}
        //if (i == m_editedintex[3])
        //{
        //    vlist[i]->x -= 0.1;//-= 0.75;//;
        //    vlist[i]->y += 0.3;//+= sin(m_rotate_theta);//-= 1.2;//;
        //    vlist[i]->z;//+= 0.5;//+= cos(m_rotate_theta);//+= 0.;//;  
        //}
        //------------------------------------X negative PI/2------------------------------------
        //if (i == m_editedintex[0])
        //{
        //    vlist[i]->x -= 0.4;//-= 0.75;//;
        //    vlist[i]->y;//-= 0.5;//+= sin(m_rotate_theta);//-= 1.2;//;
        //    vlist[i]->z -= 1;//+= 0.5;//+= cos(m_rotate_theta);//+= 0.;//;  
        //}
        //if (i == m_editedintex[1])
        //{
        //    vlist[i]->x -= 0.4;//-= 0.75;//;
        //    vlist[i]->y -= 1;//-= 0.5;//+= sin(m_rotate_theta);//-= 1.2;//;
        //    vlist[i]->z;//+= 0.5;//+= cos(m_rotate_theta);//+= 0.;//;  
        //}
        //if (i == m_editedintex[2])
        //{
        //    vlist[i]->x -= 0.4;//-= 0.75;//;
        //    vlist[i]->y += 1;//+= sin(m_rotate_theta);//-= 1.2;//;
        //    vlist[i]->z;//+= 0.5;//+= cos(m_rotate_theta);//+= 0.;//;  
        //}
        //if (i == m_editedintex[3])
        //{
        //    vlist[i]->x -= 0.4;//-= 0.75;//;
        //    vlist[i]->y;//+= sin(m_rotate_theta);//-= 1.2;//;
        //    vlist[i]->z += 1;//+= 0.5;//+= cos(m_rotate_theta);//+= 0.;//;  
        //}
        ////------------------------------------X positive PI/2------------------------------------
        if (i == m_editedintex[0])
        {
            vlist[i]->x -= 0.4;//-= 0.75;//;
            vlist[i]->y -= 1;//-= 0.5;//+= sin(m_rotate_theta);//-= 1.2;//;
            vlist[i]->z;//+= 0.5;//+= cos(m_rotate_theta);//+= 0.;//;  
        }
        if (i == m_editedintex[1])
        {
            vlist[i]->x -= 0.4;//-= 0.75;//;
            vlist[i]->y;//-= 0.5;//+= sin(m_rotate_theta);//-= 1.2;//;
            vlist[i]->z += 1;//+= 0.5;//+= cos(m_rotate_theta);//+= 0.;//;  
        }
        if (i == m_editedintex[2])
        {
            vlist[i]->x -= 0.4;//-= 0.75;//;
            vlist[i]->y;//+= sin(m_rotate_theta);//-= 1.2;//;
            vlist[i]->z -= 1;//+= 0.5;//+= cos(m_rotate_theta);//+= 0.;//;  
        }
        if (i == m_editedintex[3])
        {
            vlist[i]->x -= 0.4;//-= 0.75;//;
            vlist[i]->y += 1;//+= sin(m_rotate_theta);//-= 1.2;//;
            vlist[i]->z;//+= 0.5;//+= cos(m_rotate_theta);//+= 0.;//;  
        }
    }

}

//fixed red 3
//move blue 0