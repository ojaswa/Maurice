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
**           Date  : 10.10.2015                                           **
****************************************************************************/
#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    openGLWidget = ui->openGLWidget;
    m_shape = openGLWidget->getShape();
    cageDialog = new CageDialog(this);
}

MainWindow::~MainWindow()
{
    delete ui;
    delete m_shape;
}

void MainWindow::on_action_Open_Image_triggered()
{

    QFileDialog dialog(this); //declare a dialog box to help assist selection
    dialog.setNameFilter(tr("Images (*.png)")); // set type of images that are compatible

    dialog.setViewMode(QFileDialog::Detail); // set the view mode
    QStringList list; // List to store the selections
    if(dialog.exec()) //execute the dialog box
    {
        list = dialog.selectedFiles(); // capture list
    }
    m_shape->loadImage(list.first().toStdString().c_str());
    emit imageLoaded();
}

void MainWindow::on_actionAdd_Line_handles_triggered(bool checked)
{
    ui->actionAdd_Polyline_handle->setChecked(false);
}

void MainWindow::on_actionAdd_Polyline_handle_triggered(bool checked)
{
    ui->actionAdd_Line_handles->setChecked(false);
}

void MainWindow::on_actionModify_cage_triggered()
{
    cageDialog->exec();
    m_shape->postCageModify();
}
