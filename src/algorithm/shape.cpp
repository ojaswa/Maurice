/***************************************************************************
**                                                                        **
**  Maurice - a 2D deformable mapping/animation program                   **
**  Copyright (C) 2014-2016 Graphics Research Group, IIIT Delhi           **
**                                                                        **
**  This program is free software: you can redistribute it and/or modify  **
**  it under the terms of the GNU General Public License as published by  **
**  the Free Software Foundation, either version 3 of the License, or     **
**  (at your option) any later version.                                   **
**                                                                        **
**  This program is distributed in the hope that it will be useful,       **
**  but WITHOUT ANY WARRANTY; without even the implied warranty of        **
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         **
**  GNU General Public License for more details.                          **
**                                                                        **
**  You should have received a copy of the GNU General Public License     **
**  along with this program.  If not, see http://www.gnu.org/licenses/.   **
**                                                                        **
****************************************************************************
**           Author: Ojaswa Sharma                                        **
**           E-mail: ojaswa@iiitd.ac.in                                   **
**           Date  : 01.11.2016                                           **
****************************************************************************/
#include "shape.h"
#include <vector>
#include <limits.h>

#include <Eigen/Dense>
//#include <omp.h>

using namespace cv;
using namespace std;
using namespace Eigen;

Shape::Shape()
{
    srand (time(NULL));
    m_useEuclideanWeights = false;
    m_cageDilation = 2;
    m_cageSimplification = 10;
    m_baryVertices = NULL;
    m_baryHandles_p = NULL;
}

Shape::~Shape()
{
    m_image.release();
    m_alpha.release();
    if(m_baryVertices) delete []m_baryVertices;
    if(m_baryHandles_p) delete []m_baryHandles_p;
}

void Shape::loadImage(const char *imagefile)
{
    //Read image
    m_image = imread(imagefile,IMREAD_UNCHANGED);
    flip(m_image, m_image, 0);
    copyMakeBorder(m_image, m_image, IMAGE_PADDING, IMAGE_PADDING, IMAGE_PADDING, IMAGE_PADDING, BORDER_CONSTANT, Scalar::all(0)); //padding the image

    //Get alpha channel
    vector<Mat> layers;
    split(m_image, layers);
    threshold( layers[3], m_alpha, 128, 255,THRESH_BINARY );

    //Extract border
    vector<vector<Point> > contours;
    dilate(m_alpha, m_alpha, Mat(), Point(-1, -1), 1); //Dilate by one-pixel to get rid of skinny features
    findContours(m_alpha, contours, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_TC89_L1);
    //findContours(m_alpha, contours, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_NONE);
    //approxPolyDP(contours[0], m_shape, 2.0f, true);
    m_shape = contours[0];

    //Compute cage
    computeCage();

    // Triangulate shape
    triangulate();

    // Compute Gram matrix
    computeGramMatrix();

    //Compute barycentric coordinates
    precomputeBaryVertices();
}

void Shape::precomputeBaryVertices()
{
    int numPoints = m_points.size();
    int dimBary = m_cage.size();
    if(m_baryVertices) delete []m_baryVertices;
    m_baryVertices = new float[numPoints*dimBary];
    memset(m_baryVertices, numPoints*dimBary, sizeof(float));

#pragma omp parallel for
    for(int i=0; i<numPoints; i++) {
        mvc(m_points[i], m_baryVertices + i*dimBary);
    }
}

//Must be called just after any handle is created or destroyed.
void Shape::precomputeBaryHandles()
{
    int dimBary = m_cage.size();
    int nHandles = m_p.size()/4;
    if(m_baryHandles_p) delete []m_baryHandles_p;
    m_baryHandles_p = new float[2*nHandles*dimBary]; //Each handle has two points
    memset(m_baryHandles_p, nHandles*dimBary, sizeof(float));

    p2t::Point pt;
    for(int i=0; i<nHandles; i++) {
        pt.x = m_p[i*4]; pt.y = m_p[i*4+1];
        mvc(&pt, m_baryHandles_p + 2*i*dimBary);
        pt.x = m_p[i*4+2]; pt.y = m_p[i*4+3];
        mvc(&pt, m_baryHandles_p + (2*i + 1)*dimBary);
    }
}

void Shape::computeCage()
{
    vector<Mat> layers;
    split(m_image, layers);
    threshold( layers[3], m_alpha, 128, 255,THRESH_BINARY );
    dilate(m_alpha, m_alpha, Mat(), Point(-1, -1), m_cageDilation);
    vector<vector<Point> > m_dilatedContour;
    findContours(m_alpha, m_dilatedContour, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_NONE);
    approxPolyDP(m_dilatedContour[0], m_cage, (float)m_cageSimplification, true);
}

void Shape::triangulate()
{
    //Create boundary list
    long idx = 0;
    for(int i=0; i<m_shape.size(); i++)
    {
        Point &p = m_shape[i];
        m_points.push_back(new p2t::Point(p.x, p.y, idx++));
    }

    m_cdt = new p2t::CDT(m_points); // Add primay polygon

    //Add steiner points
    Rect bbox  = boundingRect(m_shape);

    float px, py;
    for(float y = bbox.y; y<(bbox.y + bbox.height); y+=TRIANGULATION_DENSITY)
        for(float x = bbox.x; x<(bbox.x + bbox.width); x+=TRIANGULATION_DENSITY)
        {
            //Generate a random stratified point inside the polygon
            px = x + (float(rand())/float(RAND_MAX) - 0.5)*TRIANGULATION_RANDOMIZE_DIST;
            py = y + (float(rand())/float(RAND_MAX) - 0.5)*TRIANGULATION_RANDOMIZE_DIST;
            if(pointPolygonTest(m_shape, Point2f(px, py), false) > 0.0) {
                p2t::Point *pt = new p2t::Point(px, py, idx++);
                m_cdt->AddPoint(pt);
                m_points.push_back(pt);
            }
        }

    m_cdt->Triangulate();
}

//Code to Compute the barycentric cordinates for a cage and a query cordinate
//The caller is responsible to allocate 'bary'.
int Shape::mvc(p2t::Point *p, float *bary)
{
    int nSize = m_cage.size();
    assert( nSize );

    double dx, dy;

    vector<Point2d> s(nSize);
    for( int i = 0; i < nSize; i++)
    {
        dx  =   m_cage[i].x - p->x;
        dy  =   m_cage[i].y - p->y;
        s[i].x  =   dx;
        s[i].y  =   dy;
    }

    double epsilon = 10.0*std::numeric_limits<double>::min();
    int ip, im;      // (i+1) and (i-1)
    double rp, dl, mu;  // Distance

    double ri[nSize], Ai[nSize], Di[nSize];

    // First check if any coordinates close to the cage point or
    // lie on the cage boundary. These are special cases.
    for( int i = 0; i < nSize; i++)
    {
        ip = (i+1)%nSize;
        ri[i] = sqrt( s[i].x*s[i].x + s[i].y*s[i].y );
        Ai[i] = 0.5*(s[i].x*s[ip].y - s[ip].x*s[i].y);
        Di[i] = s[ip].x*s[i].x + s[ip].y*s[i].y;
        if( ri[i] <= epsilon)
        {
            bary[i] = 1.0;
            return 0;
        }
        if( fabs(Ai[i]) <= epsilon && Di[i] < 0.0)
        {
            dx = m_cage[ip].x - m_cage[i].x;
            dy = m_cage[ip].y - m_cage[i].y;
            dl = sqrt(dx*dx + dy*dy);
            assert(dl > epsilon);
            dx = p->x - m_cage[i].x;
            dy = p->y - m_cage[i].y;
            mu = sqrt(dx*dx + dy*dy)/dl;
            assert( mu >= 0.0 && mu <= 1.0);
            bary[i]  = 1.0-mu;
            bary[ip] = mu;
            return 0;
        }
    }

    double w[nSize], wsum = 0.0;
    for( int i = 0; i < nSize; i++)
    {
        w[i] = 0;
        ip = (i+1)%nSize;
        im = (nSize-1+i)%nSize;
        if(fabs(Ai[im]) > epsilon)
            w[i] = w[i] + (ri[im] - Di[im]/ri[i]) / Ai[im];
        if (fabs(Ai[i]) > epsilon)
            w[i] = w[i] + (ri[ip]  - Di[i]/ri[i]) / Ai[i];
    }

    for (int i = 0; i < nSize; i++)
        wsum += w[i];

    for (int i = 0; i < nSize; i++)
    {
        w[i] /= wsum;
        bary[i] = w[i];
    }

    return 0;
}

int Shape::mvc_omp(const Point &p, vector<double> &bary)
{
    int nSize = m_cage.size();
    assert( nSize );
    double dx, dy;
    vector<Point2d> s(nSize);
    double epsilon = 10.0*std::numeric_limits<double>::min();
    int ip, im;      // (i+1) and (i-1)
    double dl, mu;  // Distance
    double ri[nSize], Ai[nSize], Di[nSize];
    bary.resize(nSize);
    double w[nSize], wsum = 0.0;
    bool break_loop = false;
    int break_loop_i = 0, break_loop_ip = 0;
    float break_loop_i_val, break_loop_ip_val;

#pragma omp parallel
{
#pragma omp for
    for( int i = 0; i < nSize; i++)
    {
        s[i].x  =   m_cage[i].x - p.x;
        s[i].y  =   m_cage[i].y - p.y;
    }

#pragma omp for
    for( int i = 0; i < nSize; i++) {
        bary[i] = 0.0;}

    // First check if any coordinates close to the cage point or
    // lie on the cage boundary. These are special cases.
#pragma omp for private(ip, im, dx, dy, dl, mu)
    for( int i = 0; i < nSize; i++)
    {
        ip = (i+1)%nSize;
        ri[i] = sqrt( s[i].x*s[i].x + s[i].y*s[i].y );
        Ai[i] = 0.5*(s[i].x*s[ip].y - s[ip].x*s[i].y);
        Di[i] = s[ip].x*s[i].x + s[ip].y*s[i].y;
        if( ri[i] <= epsilon)
        {
            //bary[i]  = 1.0;
            break_loop = true;
            break_loop_i = i;
            break_loop_ip = -1;
            break_loop_i_val = 1.0;
        }
        if( fabs(Ai[i]) <= epsilon && Di[i] < 0.0)
        {
            dx = m_cage[ip].x - m_cage[i].x;
            dy = m_cage[ip].y - m_cage[i].y;
            dl = sqrt(dx*dx + dy*dy);
            assert(dl > epsilon);
            dx = p.x - m_cage[i].x;
            dy = p.y - m_cage[i].y;
            mu = sqrt(dx*dx + dy*dy)/dl;
            assert( mu >= 0.0 && mu <= 1.0);
            //bary[i]  = 1.0-mu;
            //bary[ip] = mu;

            break_loop = true;
            break_loop_i = i;
            break_loop_ip = ip;
            break_loop_i_val = 1.0-mu;
            break_loop_ip_val = mu;
        }
    }

#pragma omp master
    {
        if(break_loop) {
            for( int i = 0; i < nSize; i++)
                bary[i] = 0.0;
            bary[break_loop_i] = break_loop_i_val;
            if(break_loop_ip > -1)
                bary[break_loop_ip] = break_loop_ip_val;
        }
    }

#pragma omp for private(ip, im)
    for( int i = 0; i < nSize; i++)
    {
        w[i] = 0;
        ip = (i+1)%nSize;
        im = (nSize-1+i)%nSize;
        if(fabs(Ai[im]) > epsilon)
            w[i] = w[i] + (ri[im] - Di[im]/ri[i]) / Ai[im];
        if (fabs(Ai[i]) > epsilon)
            w[i] = w[i] + (ri[ip]  - Di[i]/ri[i]) / Ai[i];
    }

#pragma omp for reduction(+:wsum)
    for (int i = 0; i < nSize; i++) {
        wsum += w[i];}

#pragma omp for
    for (int i = 0; i < nSize; i++)
    {
        w[i] /= wsum;
        bary[i] = w[i];
    }
}
    return 0;
}

void Shape::computeGramMatrix()
{
    int size = m_cage.size();
    //Calculate double centering matrix J
    MatrixXd J = MatrixXd::Identity(size, size) - MatrixXd::Ones(size, size)/float(size);

    //Compute Distance matrix (pairwise shortest distance)
    MatrixXd D = MatrixXd(size, size);

    if(!m_useEuclideanWeights) {
        MatrixXi V = MatrixXi(size, size);
        computeVisibilityMatrix(V);
        computeDistanceMatrix(V, D);

        //Square the D-matrix
        for(int i = 0; i < size; i++)
            for(int j = 0; j < size; j++)
                D(i,j) = D(i, j) * D(i,j);
    } else {
        double dx, dy;
        for(int i=0; i<size; i++)
            D(i, i) = 0.0;
        for(int i = 0; i < (size-1); i++)
            for(int j = (i+1); j < size; j++) {
                dx = m_cage[i].x - m_cage[j].x;
                dy = m_cage[i].y - m_cage[j].y;
                D(i, j) = dx*dx + dy*dy;
                D(j, i) = D(i, j);
            }
    }
    // Compute Gram matrix A
    m_A = -0.5 * J*D*J;

    //Compute nearest positive definite approximation
    approxPD();
}

void Shape::approxPD()
{
    int size = m_cage.size();

    SelfAdjointEigenSolver<MatrixXd> eigensolver(m_A);
    VectorXd lambda = eigensolver.eigenvalues();
    MatrixXd V = eigensolver.eigenvectors();
    for(int i=0; i<size; i++)
        if(lambda(i) < 0.0) lambda(i) = 1.0e-8;
    MatrixXd D = lambda.asDiagonal();
    m_A = V * D * V.inverse();
}

// Use the cage polygon to compute visibility matrix
void Shape::computeVisibilityMatrix(MatrixXi &V)
{
    Mat polypoints = Mat(m_cage); //get points extracted from cage
    // polypoints : polygon points
    // polypoints stores the vector<Point> of the contour
    // polypoints stores points in contigous locations (x1,y1,x2,y2,x3,y3.....)
    // points that are adjacent to each other are connected
    for (int ii=0;ii<polypoints.rows;ii++)
    {
        for(int jj=0;jj<polypoints.rows;jj++)
        {
            V(ii, jj)=0; // intialize the adjacency matrix with zero
            if(ii==jj)
                V(ii, jj)=1; // diagonal vertices are connected
            if(ii==0)
                V(ii, polypoints.rows-1)=1; // connect first point to last
            if(ii==polypoints.rows-1)
                V(ii, 0)=1; // connect last point to first (to maintain symmetric nature)
            if(abs(ii-jj)==1) // if points are adjacent
                V(ii, jj)=1; // connect the adjacent points
        }
    }

    int numVertices=polypoints.rows;//number of vertices

    //*************partial outside test begins ************************//
    //****this test takes the adjacency matrix , and checks for each edge if it lies outside the bpundary
    double boundaryX1,boundaryX2,boundaryY1,boundaryY2; // stores a boundary edge (which is connected in adjmat.
    //Edge is represented by (boundaryX1,boundaryY1) and (boundaryX2,boundaryY2))
    int intersect=0;//find intersection between two edges

    for (int i=0;i<numVertices;i++)
        for (int j=0;j<numVertices;j++)
        {
            //getting non connected vertices and rejecting edges
            if((V(i, j)!=1)&&(i!=j)&&(i>=j)) // getting candidate edges
            {
                //now we have two vertices and we need to check partial outside
                double x1=polypoints.at<int>(2*i); //getting the cordinates of vertices of candidate edge
                double y1=polypoints.at<int>(2*i+1);
                double x2=polypoints.at<int>(2*j);
                double y2=polypoints.at<int>(2*j+1);

                for (int k=0; k<numVertices-1;k++)
                {
                    if(k!=numVertices-1)
                    {
                        boundaryX1=polypoints.at<int>(2*k); //getting the cordinates of boundary points
                        boundaryY1=polypoints.at<int>(2*k+1);
                        boundaryX2=polypoints.at<int>(2*(k+1));
                        boundaryY2=polypoints.at<int>(2*(k+1)+1);
                    }
                    else
                    {   boundaryX1=polypoints.at<int>(2*k); //getting the cordinates of boundary points
                        boundaryY1=polypoints.at<int>(2*k+1);
                        boundaryX2=polypoints.at<int>(0);
                        boundaryY2=polypoints.at<int>(1);
                    }
                    int indpi=i-1;
                    if((i-1)<0)
                        indpi=numVertices-1; //predeccesor vertex of i
                    int indsi=i+1;
                    if((i+1)>numVertices-1)
                        indsi=0; //successor of vertex i


                    int indpj=j-1;
                    if((j-1)<0)
                        indpj=numVertices-1; //predecessor vertex of j
                    int indsj=j+1;
                    if((j+1)>numVertices-1)
                        indsj=0; //successor vertex of j

                    if(!((indpi==k && i==(k+1)) || (i==k && indsi==(k+1)) || (indpj==k && j==(k+1))||(j==k && indsj==(k+1)))) //if current cordinates for k are  boundary edges with respect to i and j
                    {
                        // to test intersection : check if points of one line lie on the same side of other line
                        // test 1: if no
                        // test 2 : if no
                        // if both no then they intersect
                        // if even one intersects remove it
                        intersect = edgeIntersect(boundaryX1,boundaryY1,boundaryX2,boundaryY2,x1,y1,x2,y2); // tells if two edges will intersect
                        if(intersect==1)
                        {
                            V(i, j)=INT_MAX; //remove the edge if it intersects with non neighbour boundary
                            V(j, i)=INT_MAX;

                            intersect=0;
                            break;
                        }
                    }
                }
            }
        }
    //*************************polygon point code starts to exclude edges that are totally outside the boundary ***************************//
    //checks if mid point is outside the polygon
    float distpointpoly; // stores the distance of point from polygon
    float ptx,pty;

    for (int a=0;a<numVertices;a++)
    {
        for(int b=0;b<numVertices;b++)
        {
            if((a!=b)&&V(a, b)!=1 && V(a, b)!=INT_MAX) // for  points in the adjacency matrix that are still not connected
            {
                double x1=polypoints.at<int>(2*a); // get a line from the points
                double y1=polypoints.at<int>(2*a+1);
                double x2=polypoints.at<int>(2*b);
                double y2=polypoints.at<int>(2*b+1);

                ptx=(x1+x2)/2; // get the mid point cordinates
                pty=(y1+y2)/2;

                distpointpoly=pointPolygonTest(m_cage,Point2f(ptx,pty),true);// check distance of mid point and polygon

                if(distpointpoly<=2)
                {
                    V(a, b)=INT_MAX; // remove (mark distance infinity) if poutside the polygon
                    V(b, a)=INT_MAX;
                }
            }
        }
    }

    //Sanity -- required??
    for (int i=0; i<m_cage.size();i++)
        for (int j=0;j<m_cage.size();j++)
            if(V(i,j)!=INT_MAX && i!=j)
                V(i,j)=1; //if the pair of vertices has passed all tests make it 1 and add it to the drawinf
}

void Shape::computeDistanceMatrix(MatrixXi &V, MatrixXd &D)
{
    int numVertices  = m_cage.size();
    MatrixXd cost = MatrixXd(numVertices, numVertices);
    Mat polypoints = Mat (m_cage); // The points of contour

    for(int i=0;i<numVertices;i++)
        for(int j=0;j<numVertices;j++)
        {
            if(V(i,j)==1) // If there is a path (edge) get the distance
            {
                double pt1x=polypoints.at<int>(2*i);
                double pt1y=polypoints.at<int>(2*i+1);
                double pt2x=polypoints.at<int>(2*j);
                double pt2y=polypoints.at<int>(2*j+1);
                cost(i,j) =  sqrt((pt1x - pt2x)*(pt1x - pt2x) + (pt1y - pt2y)*(pt1y - pt2y)); //  Distance between the vertex
            }
            else
                cost(i,j)=INT_MAX;
        }
    floydWarshall(cost, D, numVertices);
}

//Floyd-Warshall shortest path
void Shape::floydWarshall(MatrixXd &cost, MatrixXd &dist,int numVertices)
{
    int i, j, k;

    for( i = 0; i < numVertices; i++)
        for( j = 0; j < numVertices; j++)
            dist(i,j) = cost(i,j); // Intialize the distances

    for(k = 0; k < numVertices; k++)
        for (i = 0; i < numVertices; i++)
            for( j = 0; j < numVertices; j++)
                if(dist(i,k) + dist(k,j) < dist(i,j))
                    dist(i,j) = dist(i,k) + dist(k,j); //Floyd Warshall
}

// Checks if two edges of polygon intersect
int Shape::edgeIntersect(double boundaryX1,double boundaryY1,double boundaryX2,double boundaryY2,double x1,double y1,double x2,double y2) //returns if two edges intersect
{
    int result1=0,result2=0; //result of test 1 and test 2
    if( ((y1-y2)*(boundaryX1-x1)+ (x2-x1)*(boundaryY1-y1))*((y1-y2)*(boundaryX2-x1)+(x2-x1)*(boundaryY2-y1)  ) > 0) // check if boundary edge lies on same side of candidate edge
        result1=1;
    else
        result1=0;

    if( ((boundaryY1-boundaryY2)*(x1-boundaryX1) + (boundaryX2-boundaryX1)*(y1-boundaryY1))*((boundaryY1-boundaryY2)*(x2-boundaryX1)+(boundaryX2-boundaryX1)*(y2-boundaryY1)  ) > 0)
        result2=1; // check if candidate edge lies on same side of boundary edge
    else
        result2=0;
    if(result1==0 && result2==0) // if both false return psoitive intersection
        return 1;
    else
        return 0;
}

void Shape::deform(int index, float *deformedPoint)
{
    //int verticesCount = approx_contour.size();
    //Matrix4d weight[lineCount];
    int nHandles = m_q.size()/4;
    Matrix4d *weight = new Matrix4d[nHandles];
    Vector2d pStar, qStar;

    double delta_00[nHandles], delta_01[nHandles], delta_11[nHandles];
    Matrix2d aCapMatrix;
    aCapMatrix.setZero(2, 2);

    int dimBary = m_cage.size();
    float *baryA, *baryB;
    float *baryV = m_baryVertices + index*dimBary;
    int x1, x2, y1, y2;

    p2t::Point a, b, v;      //a, b representative of end points of a line in initial line handles. v representative of pointToDeform

    Vector2d vVector;
    vVector<<m_points[index]->x, m_points[index]->y;

    v.x = m_points[index]->x;
    v.y = m_points[index]->y;

    MatrixXd phi(dimBary, 2);
    double d,e,f,g, g_square, gamma, lineDerivativeMagnitude;

    //Calculating 3 integrals
    for(int i = 0; i < nHandles; i++)
    {
        x1 = m_p[4*i];
        y1 = m_p[4*i + 1];
        x2 = m_p[4*i + 2];
        y2 = m_p[4*i + 3];

        a.x = x1;
        a.y = y1;

        b.x = x2;
        b.y = y2;

        lineDerivativeMagnitude = sqrt( (x2 - x1) * (x2 - x1) + (y2 - y1)*(y2 - y1) );
        baryA = m_baryHandles_p + 2*i*dimBary;
        baryB = m_baryHandles_p + (2*i+1)*dimBary;

        for(int j=0; j<dimBary; j++) {
            phi(j, 0) = baryB[j] - baryA[j];
            phi(j, 1) = baryA[j] - baryV[j];
        }

        aCapMatrix = phi.transpose() * m_A * phi;

        d = aCapMatrix(0,0);
        e = aCapMatrix(0,1);
        f = aCapMatrix(1,1);

        g_square = aCapMatrix.determinant();
        g = sqrt(g_square);

        gamma = atan((d+e)/ g ) - atan(e/g);

        delta_00[i] = (lineDerivativeMagnitude / (2 * g_square)) * ( ( -1 * (e + f) / f) + ( ((d + 2*e + f) * gamma) / g));
        delta_01[i] = (lineDerivativeMagnitude / (2 * g_square)) * ( 1 - ( ((e + f) * gamma) / g ) );
        delta_11[i] = (lineDerivativeMagnitude / (2 * g_square)) * ( ( -1 * (e + f) / (d + 2*e + f) ) + (f * gamma / g ));

        weight[i].setZero(4 , 4);
        weight[i](0,0) = delta_00[i];
        weight[i](1,1) = delta_00[i];
        weight[i](0,2) = delta_01[i];
        weight[i](1,3) = delta_01[i];
        weight[i](2,0) = delta_01[i];
        weight[i](3,1) = delta_01[i];
        weight[i](2,2) = delta_11[i];
        weight[i](3,3) = delta_11[i];
    }

    //Vector2d pStar;
    pStar<<0,0;

    //Vector2d qStar;
    qStar<<0,0;

    Vector2d ai, bi, ci , di;
    double denominator = 0.0;

    for(int i = 0; i < nHandles; i++)
    {
        ai<<m_p[4*i], m_p[4*i+1];
        bi<<m_p[4*i + 2], m_p[4*i + 3];

        ci<<m_q[4*i], m_q[4*i+1];
        di<<m_q[4*i + 2], m_q[4*i + 3];

        ai = ai * (delta_00[i] + delta_01[i]);
        bi = bi * (delta_01[i] + delta_11[i]);

        ci = ci * (delta_00[i] + delta_01[i]);
        di = di * (delta_01[i] + delta_11[i]);

        pStar += ai + bi;
        qStar += ci + di;

        denominator += delta_00[i] + 2 * delta_01[i] + delta_11[i];
    }

    pStar = pStar / denominator;
    qStar = qStar / denominator;

    Vector2d aCap, bCap, cCap, dCap;
    MatrixXd abMat(4,2);
    MatrixXd vpMat(2,2);
    MatrixXd cdMat(1,4);
    MatrixXd rotationVector(1,2);
    MatrixXd result(4,2);

    rotationVector.setZero(1,2);
    for(int i = 0; i < nHandles; i++)
    {
        ai<<m_p[4*i], m_p[4*i+1];
        bi<<m_p[4*i + 2], m_p[4*i + 3];

        ci<<m_q[4*i], m_q[4*i+1];
        di<<m_q[4*i + 2], m_q[4*i + 3];

        aCap = ai - pStar;
        bCap = bi - pStar;
        cCap = ci - qStar;
        dCap = di - qStar;

        abMat.row(0) = aCap;
        abMat(1,0) = aCap(1);
        abMat(1,1) = -1 * aCap(0);
        abMat.row(2) = bCap;
        abMat(3,0) = bCap(1);
        abMat(3,1) = -1 * bCap(0);

        vpMat.row(0) = vVector - pStar;
        vpMat(1,0) = vpMat(0,1);
        vpMat(1,1) = -1 * vpMat(0,0);

        result = weight[i] * abMat * vpMat;

        cdMat<<cCap.transpose(), dCap.transpose();
        rotationVector += cdMat * result;
    }
     Vector2d final, temp;
    temp = rotationVector.transpose();

    final = ((vVector - pStar).norm() / temp.norm()) * temp + qStar;

    deformedPoint[0] = final(0);
    deformedPoint[1] = final(1);
    delete []weight;
}

// Add line/polyline handle
void Shape::addLineHandle(vector<float> handle)
{
    vector<double> v(handle.begin(), handle.end()); //Convert to double

    //Insert handle
    int start = m_p.size();
    m_p.insert(m_p.end(), v.begin(), v.end()); // Add to p
    m_q.insert(m_q.end(), v.begin(), v.end());// Add to q
    int end = m_p.size() - 4; //Starting index of the last element of this line/polyline handle.
    m_handleMap.push_back(glm::ivec2(start, end));

    //Recompute barycentric coordinates for the handles.
    precomputeBaryHandles();
}

void Shape::inverseKinematicsFABRIK(int handleID, int endEffector, glm::vec2 target)
{
    //Copy chain
    int sz = (m_handleMap[handleID].y - m_handleMap[handleID].x)/4 + 2; //No. of  points in the kinematic chain.
    glm::vec2 *chain = new glm::vec2[sz];
    int i, j = 0;
    if(m_handleMap[handleID].x == endEffector) {
        for(i= m_handleMap[handleID].x; i<= m_handleMap[handleID].y; i+=4) {
            chain[j].x = m_q[i];
            chain[j++].y = m_q[i+1];
        }
        chain[j].x = m_q[m_handleMap[handleID].y+2];
        chain[j].y = m_q[m_handleMap[handleID].y+3];
    } else {
        for(i=(m_handleMap[handleID].y + 2); i>=(m_handleMap[handleID].x + 2); i-=4) {
            chain[j].x = m_q[i];
            chain[j++].y = m_q[i+1];
        }
        chain[j].x = m_q[m_handleMap[handleID].x];
        chain[j].y = m_q[m_handleMap[handleID].x + 1];
    }

    //Compute distances
    float *distances = new float[sz-1];
    for(int i=0; i<(sz-1); i++)
        distances[i] = glm::distance(chain[i], chain[i+1]);

    //Run FABRIK Cycles
    glm::vec2 newp;
    glm::vec2 fixed = chain[sz - 1];
    for(int i=0; i<FABRIK_CYCLES; i++)
    {
        //Forward pass, endEffector -> fixedEnd
        chain[0] = target;
        for(int j=1; j<sz; j++)
            chain[j] = chain[j-1] + distances[j-1] * glm::normalize(chain[j] - chain[j-1]);

        //Reverese pass, fixedEnd -> endEffector
        chain[sz - 1] = fixed;
        for(int j=(sz-2); j>=0; j--)
            chain[j] = chain[j+1] + distances[j] * glm::normalize(chain[j] - chain[j+1]);
    }

    //Copy back modified chain to handle storage
    if(m_handleMap[handleID].x == endEffector) {
        int j = m_handleMap[handleID].x;
        for(int i=0; i<(sz-1); i++) {
            m_q[j++] = chain[i].x;
            m_q[j++] = chain[i].y;
            m_q[j++] = chain[i+1].x;
            m_q[j++] = chain[i+1].y;
        }
    } else {
        int j = m_handleMap[handleID].y + 3;
        for(int i=0; i<(sz-1); i++) {
            m_q[j--] = chain[i].y;
            m_q[j--] = chain[i].x;
            m_q[j--] = chain[i+1].y;
            m_q[j--] = chain[i+1].x;
        }
    }
    delete []chain;
}

void Shape::dilateCage(int value)
{
    m_cageDilation = value;

    //Compute cage
    computeCage();
}

void Shape::simplifyCage(int value)
{
    m_cageSimplification = value;

    //Compute cage
    computeCage();
}

//To be called once the cage modification is over.
void Shape::postCageModify()
{
    computeGramMatrix();
    precomputeBaryVertices();
    precomputeBaryHandles();

}
