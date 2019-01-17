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
**           Date  : 23.09.2015                                           **
****************************************************************************/
#include "cagedialog.h"
#include "ui_cagedialog.h"
#include "mainwindow.h"

CageDialog::CageDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::CageDialog)
{
    ui->setupUi(this);
    setFixedSize(sizeHint());
}

CageDialog::~CageDialog()
{
    delete ui;
}

void CageDialog::on_spinBoxCageDilation_valueChanged(int arg1)
{
    ((MainWindow*)parent())->openGLWidget->dilateCage(arg1);
}

void CageDialog::on_spinBoxCageThreshold_valueChanged(int arg1)
{
    ((MainWindow*)parent())->openGLWidget->simplifyCage(arg1);
}
