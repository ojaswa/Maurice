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
#ifndef SHAPE_H
#define SHAPE_H

#include "opencv2/highgui/highgui.hpp"
#include "opencv2/opencv.hpp"
#include "poly2tri.h"
#include <Eigen/Dense>

#define GLM_COMPILER 0
#define GLM_FORCE_RADIANS
#include <glm/gtc/type_ptr.hpp>
#include <glm/glm.hpp>

#define TRIANGULATION_DENSITY 10 //triangle side of 10 pixels
#define TRIANGULATION_RANDOMIZE_DIST 5
#define IMAGE_PADDING 20
#define FABRIK_CYCLES 2

// Shape Declarations
class Shape {
public:
    Shape();
    ~Shape();
    void loadImage(const char *imagefile);
    void triangulate();
    int mvc(p2t::Point *p, float *bary);
    int mvc_omp(const cv::Point &m_p, std::vector<double> &bary);
    void computeGramMatrix();
    void approxPD();
    void deform(int index, float *deformedPoint);
    void addLineHandle(std::vector<float>);
    void inverseKinematicsFABRIK(int handleID, int endEffector, glm::vec2 target);
    void dilateCage(int);
    void simplifyCage(int);
    void postCageModify();

    //Public variables
    cv::Mat m_image;
    cv::Mat m_alpha;
    std::vector<cv::Point> m_shape;
    std::vector<cv::Point> m_cage;
    p2t::CDT *m_cdt;
    std::vector<p2t::Point *>m_points;
    Eigen::MatrixXd m_A;
    bool m_useEuclideanWeights;
    int m_cageDilation, m_cageSimplification;

    //Handles
    std::vector<double> m_p; // Line handles. Each line handle is four floats: [v0x, v0y, v1x, v1y]
    std::vector<double> m_q; // Modified line handles
    std::vector<glm::ivec2> m_handleMap; // Stores a map of line/polyline handles. Each entry is a range <i,j> indexing into p or q.
    //i points to (starting) index of start of a handle; j points to (starting) index of end of a handle
    // for a line handle i == j.
    float *m_baryVertices; //barycentric coordinates for all points of triangulation --precomputed
    float *m_baryHandles_p; //barycentric coordinates for all line handles --precomputed

private:
    //Private functions
    void computeVisibilityMatrix(Eigen::MatrixXi &V);
    void computeDistanceMatrix(Eigen::MatrixXi &V, Eigen::MatrixXd &D);
    int edgeIntersect(double ,double ,double ,double ,double ,double ,double ,double );
    void floydWarshall(Eigen::MatrixXd &cost, Eigen::MatrixXd &dist,int numVertices);
    void computeCage();
    void precomputeBaryVertices();
    void precomputeBaryHandles();
};

#endif // SHAPE_H
