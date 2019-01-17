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
**           Date  : 28.10.2015                                           **
****************************************************************************/
#include "shader_utils.h"
#include "openglwidget.h"
#include "gl_utils.h"

#include <QMouseEvent>
#include <QKeyEvent>
#include <QFileDialog>
#include <QMessageBox>
#include <QCoreApplication>
//#include <QProcess>

//#include <omp.h>

using namespace std;

OpenGLWidget::OpenGLWidget(QWidget *parent) : QOpenGLWidget(parent)
{
    //Widget specific
    setMouseTracking(true);
    setFocusPolicy(Qt::StrongFocus);

    //OpenGL specific
    m_shape = new Shape();
    m_isShapeTriangulated = false;
    m_drawWireframe = false;
    m_isLineDrawing = false;
    m_selectedNodeIndex = -1; // None selected
    m_hideHandles = false;
    m_showCage = false;

    m_shapeVertices = NULL;

    m_recordVideo = false;

    m_videoPath = QCoreApplication::applicationDirPath();
    m_videoPath.append("/video/");
    fprintf(stderr, "Video screenshots save path: %s\n", m_videoPath.toLocal8Bit().constData());
}

OpenGLWidget::~OpenGLWidget()
{
    glDeleteProgram(m_shapeProgram);
    glDeleteBuffers(1, &m_shapeVBO);
    glDeleteBuffers(1, &m_shapeIBO);
    glDeleteTextures(1, &m_shapeTexture);
    glDeleteBuffers(1, &m_cageVBO);
    if(m_shapeVertices) delete []m_shapeVertices;
}

void OpenGLWidget::initializeGL()
{
    initializeOpenGLFunctions();
    glEnable(GL_MULTISAMPLE); //Draw smoothed polygons
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_POINT_SMOOTH);

    //Create program
    m_shapeProgram = createProgram("./shaders/shapeTextured.vs", "./shaders/shapeTextured.fs");
    m_lineProgram = createProgram("./shaders/lineHandle.vs", "./shaders/lineHandle.fs");
    setupProjectionTransformation();
}

Shape* OpenGLWidget::getShape()
{
    return m_shape;
}

void OpenGLWidget::resizeGL(int w, int h)
{
    m_screenWidth = w;
    m_screenHeight = h;
    glViewport(0, 0, m_screenWidth, m_screenHeight);
    setupProjectionTransformation();
}

void OpenGLWidget::paintGL()
{
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);

    if(m_isShapeTriangulated) {
        glUseProgram(m_shapeProgram);
        glUniform2f(m_shapeImageSizeUniform, float(m_shape->m_image.cols), float(m_shape->m_image.rows));
        if (m_drawWireframe) {
            glColor3f(0.0, 0.0, 0.0);
            glLineWidth(1.0);
            glDisable(GL_TEXTURE_2D);
            glPolygonMode ( GL_FRONT_AND_BACK, GL_LINE ) ;
            glUniform1f(m_shapeUniformWireframe, 1.0);
        } else
        {
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glUniform1f(m_shapeUniformWireframe, 0.0);

            glEnable(GL_TEXTURE_2D);
            glActiveTexture(GL_TEXTURE0);
            glUniform1i(m_shapeUniformTexture, /*GL_TEXTURE*/0);
            glBindTexture(GL_TEXTURE_2D, m_shapeTexture);
        }

        glEnableVertexAttribArray(m_shapeVertexAttrib);
        glBindBuffer(GL_ARRAY_BUFFER, m_shapeVBO);
        glVertexAttribPointer(m_shapeVertexAttrib, 2, GL_FLOAT, GL_FALSE, 0, 0);

        glEnableVertexAttribArray(m_shapeTextureAttrib);
        glBindBuffer(GL_ARRAY_BUFFER, m_texVBO);
        glVertexAttribPointer(m_shapeTextureAttrib, 2, GL_FLOAT, GL_FALSE, 0, 0);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_shapeIBO);
        int size;  glGetBufferParameteriv(GL_ELEMENT_ARRAY_BUFFER, GL_BUFFER_SIZE, &size);
        glDrawElements(GL_TRIANGLES, size/sizeof(GLuint), GL_UNSIGNED_INT, 0);

        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        glDisableVertexAttribArray(m_shapeVertexAttrib);
    }

    //Display line handles
    glUseProgram(m_lineProgram);
    glUniform2f(m_lineImageSizeUniform, float(m_shape->m_image.cols), float(m_shape->m_image.rows));
    glEnableVertexAttribArray(m_lineVertexAttrib);

    //Transient rendering for handles being drawn currently
    if(!m_hideHandles) {
        glLineWidth(8.0);
        if((m_modeAddLineHandle || m_modeAddPolylineHandle) && m_isLineDrawing) {
            glVertexAttribPointer(m_lineVertexAttrib, 2, GL_FLOAT, GL_FALSE, 0, m_linePoints.data());
            glColor4f(0.0, 1.0, 1.0, 0.4);
            glDrawArrays(GL_LINES, 0, m_linePoints.size()/2);

            glPointSize(32.0);
            glVertexAttribPointer(m_lineVertexAttrib, 2, GL_FLOAT, GL_FALSE, 0, m_linePoints.data());
            glColor4f(1.0, 1.0, 1.0, 0.4);
            glDrawArrays(GL_POINTS, 0, m_linePoints.size()/2);

            glPointSize(20.0);
            glVertexAttribPointer(m_lineVertexAttrib, 2, GL_FLOAT, GL_FALSE, 0, m_linePoints.data());
            glColor4f(0.0, 0.0, 0.0, 0.4);
            glDrawArrays(GL_POINTS, 0, m_linePoints.size()/2);
        }

        //Display all handles (in their current state; i.e., display q handles)
        if(m_shape->m_q.size() > 0) {
            glVertexAttribPointer(m_lineVertexAttrib, 2, GL_DOUBLE, GL_FALSE, 0, m_shape->m_q.data());
            glColor4f(0.0, 1.0, 1.0, 0.4);
            glDrawArrays(GL_LINES, 0, m_shape->m_q.size()/2);

            glPointSize(32.0);
            glVertexAttribPointer(m_lineVertexAttrib, 2, GL_DOUBLE, GL_FALSE, 0, m_shape->m_q.data());
            glColor4f(1.0, 1.0, 1.0, 0.4);
            glDrawArrays(GL_POINTS, 0, m_shape->m_q.size()/2);

            glPointSize(20.0);
            glVertexAttribPointer(m_lineVertexAttrib, 2, GL_DOUBLE, GL_FALSE, 0, m_shape->m_q.data());
            glColor4f(0.0, 0.0, 0.0, 0.4);
            glDrawArrays(GL_POINTS, 0, m_shape->m_q.size()/2);

        }
    }
    //Display selected node
//    if(m_selectedNodeIndex != -1) {
//        glColor3f(0.0, 0.6, 0.0);
//        glLineWidth(1.0);
//        float x = m_shape->m_q[m_selectedNodeIndex];
//        float y = m_shape->m_q[m_selectedNodeIndex+1];
//        float box[8] = {x - NODE_SELECT_SIZE, y - NODE_SELECT_SIZE,
//                       x + NODE_SELECT_SIZE, y - NODE_SELECT_SIZE,
//                       x + NODE_SELECT_SIZE, y + NODE_SELECT_SIZE,
//                       x - NODE_SELECT_SIZE, y + NODE_SELECT_SIZE};
//        glVertexAttribPointer(m_lineVertexAttrib, 2, GL_FLOAT, GL_FALSE, 0, box);
//        glDrawArrays(GL_LINE_LOOP, 0, 4);
//    }

    if(m_showCage) {
        glBindBuffer(GL_ARRAY_BUFFER, m_cageVBO);
        glVertexAttribPointer(m_lineVertexAttrib, 2, GL_FLOAT, GL_FALSE, 0, 0);
        glColor4f(1.0, 0.0, 0.0, 1.0);
        glLineWidth(1.0);
        glDrawArrays(GL_LINE_LOOP, 0, m_shape->m_cage.size());
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
    glDisableVertexAttribArray(m_lineVertexAttrib);
}


void OpenGLWidget::setupProjectionTransformation()
{
    //Projection transformation (Orthographic projection)
    float aspect = (float)m_screenWidth/(float)m_screenHeight;
    float image_aspect = 1.0f;
    if(m_shape) (float)m_shape->m_image.cols / (float)m_shape->m_image.rows;
    aspect /=image_aspect;

    if(image_aspect > 1.0)
        m_projection = glm::ortho(-2.0f, 2.0f, -2.0f/aspect, 2.0f/aspect);
    else
        m_projection = glm::ortho(-2.0f*aspect, 2.0f*aspect, -2.0f, 2.0f);

    //Pass on the projection matrix to the vertex shader of Shape program
    glUseProgram(m_shapeProgram);
    m_shapeProjectionUniform = glGetUniformLocation(m_shapeProgram, "vProjection");
    if(m_shapeProjectionUniform == -1){
        fprintf(stderr, "Could not bind location: vProjection\n");
        exit(0);
    }
    glUniformMatrix4fv(m_shapeProjectionUniform, 1, GL_FALSE, glm::value_ptr(m_projection));

    //Pass on the projection matrix to the vertex shader of Line program
    glUseProgram(m_lineProgram);
    m_lineProjectionUniform = glGetUniformLocation(m_lineProgram, "vProjection");
    if(m_lineProjectionUniform == -1){
        fprintf(stderr, "Could not bind location: vProjection\n");
        exit(0);
    }
    glUniformMatrix4fv(m_lineProjectionUniform, 1, GL_FALSE, glm::value_ptr(m_projection));
}

void OpenGLWidget::setupLineHandleDisplay()
{
    glUseProgram(m_lineProgram);

    //Bind shader variables
    m_lineVertexAttrib = glGetAttribLocation(m_lineProgram, "vVertex");
    if(m_lineVertexAttrib == -1) {
        fprintf(stderr, "Could not bind location: vVertex\n");
        exit(0);
    }

    m_lineImageSizeUniform = glGetUniformLocation(m_lineProgram, "vImageSize");
    if(m_lineImageSizeUniform == -1) {
        fprintf(stderr, "Could not bind location: vImageSize\n");
        exit(0);
    }

    //Setup Line style

}

void OpenGLWidget::createShape()
{
    glUseProgram(m_shapeProgram);

    //Bind shader variables
    m_shapeVertexAttrib = glGetAttribLocation(m_shapeProgram, "vVertex");
    if(m_shapeVertexAttrib == -1) {
        fprintf(stderr, "Could not bind location: vVertex\n");
        exit(0);
    }

    m_shapeTextureAttrib = glGetAttribLocation(m_shapeProgram, "vTexCoord");
    if(m_shapeTextureAttrib == -1) {
        fprintf(stderr, "Could not bind location: vTexCoord\n");
        exit(0);
    }

    m_shapeUniformWireframe = glGetUniformLocation(m_shapeProgram, "fDrawWireframe");
    if(m_shapeUniformWireframe == -1) {
        fprintf(stderr, "Could not bind location: fDrawWireframe\n");
        exit(0);
    }

    m_shapeImageSizeUniform = glGetUniformLocation(m_shapeProgram, "vImageSize");
    if(m_shapeImageSizeUniform == -1) {
        fprintf(stderr, "Could not bind location: vImageSize\n");
        exit(0);
    }

    glUseProgram(0);

    // Create vertex buffer
    std::vector<p2t::Point*> points = m_shape->m_points;
    m_shapeVertices = new GLfloat[points.size() * 2];
    for(int i=0; i<points.size(); i++) {
        m_shapeVertices[2*i + 0] = points[i]->x;
        m_shapeVertices[2*i + 1] = points[i]->y;
    }
    glGenBuffers(1, &m_shapeVBO);
    glBindBuffer(GL_ARRAY_BUFFER, m_shapeVBO);
    glBufferData(GL_ARRAY_BUFFER, points.size() * 2 * sizeof(GLfloat), m_shapeVertices, GL_DYNAMIC_DRAW);

    //Create texture buffer
    GLfloat width = m_shape->m_image.cols;
    GLfloat height = m_shape->m_image.rows;
    GLfloat *texCoords = new GLfloat[points.size() * 2];
    for(int i=0; i<points.size(); i++) {
        texCoords[2*i + 0] = m_shapeVertices[2*i + 0] / width;
        texCoords[2*i + 1] = m_shapeVertices[2*i + 1] / height;
    }
    glGenBuffers(1, &m_texVBO);
    glBindBuffer(GL_ARRAY_BUFFER, m_texVBO);
    glBufferData(GL_ARRAY_BUFFER, points.size() * 2 * sizeof(GLfloat), texCoords, GL_STATIC_DRAW);
    delete []texCoords;

    //Create index buffer
    vector<p2t::Triangle *> triangles = m_shape->m_cdt->GetTriangles();
    GLuint *shape_indices = new GLuint[triangles.size() * 3];
    for (int i=0; i<triangles.size(); i++)
    {
        p2t::Triangle& t = *triangles[i];
        shape_indices[3*i + 0] = t.GetPoint(0)->index;
        shape_indices[3*i + 1] = t.GetPoint(1)->index;
        shape_indices[3*i + 2] = t.GetPoint(2)->index;
    }
    glGenBuffers(1, &m_shapeIBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_shapeIBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, triangles.size() * 3 * sizeof(GLuint), shape_indices, GL_STATIC_DRAW);
    delete []shape_indices;

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    //Create texture
    //use fast 4-byte alignment (default anyway) if possible
    glPixelStorei(GL_UNPACK_ALIGNMENT, (m_shape->m_image.step & 3) ? 1 : 4);
    //set length of one complete row in data (doesn't need to equal image.cols)
    glPixelStorei(GL_UNPACK_ROW_LENGTH, m_shape->m_image.step/m_shape->m_image.elemSize());

    glGenTextures(1, &m_shapeTexture);
    glBindTexture(GL_TEXTURE_2D, m_shapeTexture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, // target
                 0,  // level, 0 = base, no minimap,
                 GL_RGBA, // internalformat
                 m_shape->m_image.cols,  // width
                 m_shape->m_image.rows,  // height
                 0,  // border, always 0 in OpenGL ES
                 GL_BGRA,  // format
                 GL_UNSIGNED_BYTE, // type
                 m_shape->m_image.ptr());

    //Create cage VBO
    createCage();
}

void OpenGLWidget::createCage()
{
    std::vector<cv::Point> &cageVec = m_shape->m_cage;
    float *cageArr = new float[2*cageVec.size()];
    for(int i=0; i<cageVec.size(); i++) {
        cageArr[2*i] = cageVec[i].x;
        cageArr[2*i+1] = cageVec[i].y;
    }
    glGenBuffers(1, &m_cageVBO);
    glBindBuffer(GL_ARRAY_BUFFER, m_cageVBO);
    glBufferData(GL_ARRAY_BUFFER, cageVec.size()*2*sizeof(GL_FLOAT), cageArr, GL_STREAM_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    delete []cageArr;
}

// Convert mouse coordinates from screen space to image space
glm::vec2 OpenGLWidget::screenToImage(int x, int y)
{
    float nx = 2.0 * float(x)/ float(m_screenWidth) - 1.0;
    float ny = -2.0 * float(y)/ float(m_screenHeight) + 1.0;
    glm::vec4 in = glm::vec4(nx, ny, 0.0, 1.0);
    glm::vec4 out = glm::inverse(m_projection) * in;

    //Convert to image space
    GLfloat image_w = m_shape->m_alpha.cols;
    GLfloat image_h = m_shape->m_alpha.rows;
    out.x = 0.5*image_w*(out.x + 1.0);
    out.y = 0.5*image_h*(out.y + 1.0);
    return glm::vec2(out);
}

void OpenGLWidget::deformVertices()
{
    //Deform vertices
    int numPoints = m_shape->m_points.size();
    //int nthreads = omp_get_num_threads();

#pragma omp parallel for
    for(int i=0; i<numPoints; i++) {
        m_shape->deform(i, m_shapeVertices + 2*i);
    }

    //Update VBO
    glUseProgram(m_shapeProgram);
    glBindBuffer(GL_ARRAY_BUFFER, m_shapeVBO);
    glBufferData(GL_ARRAY_BUFFER, numPoints * 2 * sizeof(GLfloat), m_shapeVertices, GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    update();
}

//Mouse events
void OpenGLWidget::mouseMoveEvent(QMouseEvent *event)
{
    glm::vec2 p = screenToImage(event->x(), event->y());
    if(m_modeAddLineHandle || m_modeAddPolylineHandle) { // Draw rubberband handle to add new
        if(m_isLineDrawing) {
            m_linePoints[m_linePoints.size() - 2] = p.x;
            m_linePoints[m_linePoints.size() - 1] = p.y;
        }
    } else {
        if ((event->buttons() & Qt::LeftButton) && m_selectedNodeIndex > -1) { // Move line handle
            vector<double> &handles = m_shape->m_q;
            vector<glm::ivec2> &map = m_shape->m_handleMap;
            if(map[m_selectedHandleIndex].x == map[m_selectedHandleIndex].y) { //A line handle --constrain length of handle
                int fixed = map[m_selectedHandleIndex].x;
                fixed += (m_selectedNodeIndex == fixed) ? 2 : 0;
                glm::vec2 p0(handles[fixed], handles[fixed+1]);
                glm::vec2 p1(handles[m_selectedNodeIndex], handles[m_selectedNodeIndex+1]);
                float dist = glm::distance(p0, p1);
                glm::vec2 dir(p - p0);
                p1 = p0 + dist * glm::normalize(dir);
                handles[m_selectedNodeIndex] = p1.x;
                handles[m_selectedNodeIndex+1] = p1.y;
            } else { // A polyline handle -- Perform IK
                m_shape->inverseKinematicsFABRIK(m_selectedHandleIndex, m_selectedNodeIndex, p);
            }
            deformVertices();
        } else {//Search for a nearest node and select
            vector<glm::ivec2> &map = m_shape->m_handleMap;
            if (map.size() == 0) return;
            vector<double> handles = m_shape->m_q;
            float dist = (float)INT_MAX;
            m_selectedNodeIndex = -1;
            m_selectedHandleIndex = -1;
            for(int i=0; i<map.size(); i++) {
                int idx1 = map[i].x; //First handle - first coordinate
                glm::vec2 p1(handles[idx1], handles[idx1+1]);
                float dist1 = glm::distance(p, p1);
                if(dist1 < HANDLE_SELECT_THRESHOLD)
                    if (dist1 < dist) {
                        dist = dist1;
                        m_selectedNodeIndex = idx1;
                        m_selectedHandleIndex = i;
                    }

                int idx2 = map[i].y + 2; // Last handle - last coordinate
                glm::vec2 p2(handles[idx2], handles[idx2+1]);
                float dist2 = glm::distance(p, p2);
                if(dist2 < HANDLE_SELECT_THRESHOLD)
                    if (dist2 < dist) {
                        dist = dist2;
                        m_selectedNodeIndex = idx2;
                        m_selectedHandleIndex = i;
                    }
            }

        }
    }
    update();
    if(m_recordVideo) {
        QString path = QString("%1%2.png").arg(m_videoPath).arg(m_videoFrame++, 4, 10, QLatin1Char('0'));
        takeSnapshot(path);
    }
}

void OpenGLWidget::mousePressEvent(QMouseEvent *event)
{
    if(m_modeAddLineHandle) {
        if (event->button() == Qt::LeftButton) {
            if(!m_isLineDrawing) { // Start drawing line handle
                glm::vec2 ndcs = screenToImage(event->x(), event->y());
                m_linePoints.clear();
                m_linePoints.push_back(ndcs.x);
                m_linePoints.push_back(ndcs.y);
                m_linePoints.push_back(ndcs.x);//Current point x
                m_linePoints.push_back(ndcs.y);//Current point y
                m_isLineDrawing = true;
            } else { //Done with this handle. Add to Handle list
                glm::vec2 ndcs = screenToImage(event->x(), event->y());
                m_linePoints[m_linePoints.size() - 2] = ndcs.x;
                m_linePoints[m_linePoints.size() - 1] = ndcs.y;
                m_shape->addLineHandle(m_linePoints);
                m_isLineDrawing = false;
            }
        }
    } else if(m_modeAddPolylineHandle)
    {
        if (event->button() == Qt::LeftButton) { // Start drawing polyline handle
            if(!m_isLineDrawing) m_isLineDrawing = true;
            glm::vec2 ndcs = screenToImage(event->x(), event->y());
            // Begin new segment
            m_linePoints.push_back(ndcs.x);
            m_linePoints.push_back(ndcs.y);
            m_linePoints.push_back(ndcs.x);//Current point x
            m_linePoints.push_back(ndcs.y);//Current point y
        }
    }
    update();
}

void OpenGLWidget::mouseReleaseEvent(QMouseEvent *event)
{
    if(m_modeAddPolylineHandle) {
        if(event->button() == Qt::RightButton) { //Done with this polyLine handle. Add to Handle list
            glm::vec2 ndcs = screenToImage(event->x(), event->y());
            m_linePoints[m_linePoints.size() - 2] = ndcs.x;
            m_linePoints[m_linePoints.size() - 1] = ndcs.y;
            m_shape->addLineHandle(m_linePoints);
            m_isLineDrawing = false;
            m_linePoints.clear();
        }
    }
    update();
}

void OpenGLWidget::keyPressEvent( QKeyEvent *event )
{
    //event->modifiers().testFlag(Qt::ControlModifier)
    if (event->key() == Qt::Key_V) {
        m_recordVideo = !m_recordVideo; //Toggle
        if(m_recordVideo) {
            fprintf(stderr, "Begin video frame recording...");

            m_videoFrame = 0L;
            if(!QDir(m_videoPath).exists()) //Create
                QDir().mkdir(m_videoPath);
            else{ // Empty it!
                QDir dir(m_videoPath);
                dir.setNameFilters(QStringList() << "*.*");
                dir.setFilter(QDir::Files);
                foreach(QString dirFile, dir.entryList())
                    dir.remove(dirFile);
            }
            m_timerId = startTimer(10, Qt::PreciseTimer);
            m_elapsedTimer.start();
        } else {
            killTimer(m_timerId);
            //ffmpeg -f image2 -r 60 -i %04d.png -c:v libx264 -pix_fmt yuv420p out.mp4
            //QProcess::execute ("ffmpeg -f image2 -r 60 -i video/%04d.png -c:v libx264 -pix_fmt yuv420p out.mp4");
            long duration = m_elapsedTimer.elapsed();
            fprintf(stderr, "done.\nnFrames: %d, time: %d ms, FPS: %f\n", m_videoFrame+1, duration, double(m_videoFrame+1)*1000.0/double(duration));
        }
    }
}

void OpenGLWidget::timerEvent(QTimerEvent *event)
{
    if(m_recordVideo) {
        QString path = QString("%1%2.png").arg(m_videoPath).arg(m_videoFrame++, 4, 10, QLatin1Char('0'));
        takeSnapshot(path);
    }
}

// Slots
void OpenGLWidget::imageLoaded()
{
    createShape();
    setupLineHandleDisplay();
    m_isShapeTriangulated = true;
    update();
}

void OpenGLWidget::toggleDisplayMesh(bool state)
{
    m_drawWireframe = state;
    update();
}

void OpenGLWidget::toggleDisplayCage(bool state)
{
    m_showCage = state;
    update();
}

void OpenGLWidget::toggleAddLineHandle(bool state)
{
    m_modeAddLineHandle = state;
    m_modeAddPolylineHandle = false;
    m_linePoints.clear();
}

void OpenGLWidget::toggleAddPolylineHandle(bool state)
{
    m_modeAddPolylineHandle = state;
    m_modeAddLineHandle = false;
    m_linePoints.clear();
}

void OpenGLWidget::toggleUseEuclideanDistanceWeights(bool state)
{
    m_shape->m_useEuclideanWeights = state;
    m_shape->computeGramMatrix();
    deformVertices();
    update();
}

void OpenGLWidget::resetDeformation()
{
    //Reset vertices
    std::vector<p2t::Point*> points = m_shape->m_points;
    for(int i=0; i<points.size(); i++) {
        m_shapeVertices[2*i + 0] = points[i]->x;
        m_shapeVertices[2*i + 1] = points[i]->y;
    }

    //Update VBO
    glUseProgram(m_shapeProgram);
    glBindBuffer(GL_ARRAY_BUFFER, m_shapeVBO);
    glBufferData(GL_ARRAY_BUFFER, points.size() * 2 * sizeof(GLfloat), m_shapeVertices, GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    //Reset Handles
    m_shape->m_q = m_shape->m_p;

    //Update display
    update();
}

void OpenGLWidget::deleteAllHandles()
{
    resetDeformation();

    m_shape->m_q.clear();
    m_shape->m_p.clear();
    m_shape->m_handleMap.clear();

    update();
}

void OpenGLWidget::takeSnapshot()
{
    QString filename = QFileDialog::getSaveFileName( this, "Save File", getenv( "HOME" ), "PNG Image (*.png)" );
    if( !filename.endsWith( ".png" ) && !filename.endsWith( ".png" ) )
        filename.append( ".png" );
    takeSnapshot(filename);
}

void OpenGLWidget::takeSnapshot(QString filename)
{
    QImage image = grabFramebuffer( );
    if( !image.save( filename, "PNG", 100) )
        QMessageBox::warning( this, "Save Image", "Error saving image." );
}

void OpenGLWidget::hideHandles(bool state)
{
    m_hideHandles = state;
    update();
}

void OpenGLWidget::dilateCage(int value)
{
    m_shape->dilateCage(value);
    createCage();
    update();
}

void OpenGLWidget::simplifyCage(int value)
{
    m_shape->simplifyCage(value);
    createCage();
    update();
}
