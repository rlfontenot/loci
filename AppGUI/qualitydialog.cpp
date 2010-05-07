/////////////////////////////////////////////////
//  Filename: cutdialog.cpp
//
//  Contains: Implementation of CutDialog class
/////////////////////////////////////////////////

#include <QtGui>
#include <QObject>
#include <stdlib.h>
#include <unistd.h>
#include <QSignalMapper> 
#include <QFile>
#include <QTextStream>
#include <QString>
#include "qualitydialog.h"
#include <QDialog>
#include <QErrorMessage>

bool QualityDialog::processQualityReport(QString filename){

  
  QFile file(filename);
  if (!file.open(QFile::ReadOnly | QFile::Text)) {
      QMessageBox::information(window(), "quality dialog",
                                 tr("Cannot open ") + filename + tr(" for reading!"));
    
    return false;
  }
  
  QTextStream in(&file);
  QString quality = in.readAll();
  
  QString head= "<?xml version='1.0' encoding='ISO-8859-1'?>\n" ;
  head += "<main>\n";
  QString tail = "</main>\n";
  quality = head+quality+tail;
   
  QString errorStr;
  int errorLine;
  int errorColumn;
  if (!doc.setContent(quality, true, &errorStr, &errorLine,
                      &errorColumn)) {
    QMessageBox::information(window(), tr("DOM Bookmarks"),
                             tr("Parse error at line %1, column %2:\n%3")
                             .arg(errorLine)
                             .arg(errorColumn)
                             .arg(errorStr));
    return false;
  }
  
  QDomElement root = doc.documentElement();
  if(root.isNull()){
    QMessageBox::warning(window(), tr("quality report"),
                         tr("root is NUll")
                         );
    return false;
    
  }
  QDomElement minVol_elem = root.firstChildElement("minVol");
 if(minVol_elem.isNull()){
    QMessageBox::warning(window(), tr("quality report"),
                         tr("minVol is NUll")
                         );
    return false;
    
  }
  
  if(minVol_elem.text().toDouble()<=0){
    minVol_elem.setAttribute("warning", "Negative or zero volume cell!\nThe grid quality is too poor to use!");
    
  }
  minVol_elem.setAttribute("detail", "Minimum cell volume should be greater than zero");

  QDomElement volRatio_elem = root.firstChildElement("volRatio");
   if(volRatio_elem.isNull()){
    QMessageBox::warning(window(), tr("Quality report"),
                         tr("volRatio is NUll")
                         );
    return false;
    
  }
  
  if(volRatio_elem.text().toDouble() > 1000){
    volRatio_elem.setAttribute("warning", "The volume ratio is greater than 1000!\nThe grid quality is too poor to use!"); 
    
  }
  volRatio_elem.setAttribute("detail", "The volume ratio is the ratio between maximum cell volume and the minimum cell volume.\nThe lower this number the better the mesh quality.\nValues below 10 are desirable.\nValues between 10 and 50 indicate good mesh quality,\nValues between 50 and 1000 indicate  poor mesh quality.");
  
  QDomElement angle_elem = root.firstChildElement("maxCellAngle");
 if(angle_elem.isNull()){
    QMessageBox::warning(window(), tr("quality report"),
                         tr("angle is NUll")
                         );
    return false;
    
  }
  

  if(angle_elem.text().toDouble()>179){
     angle_elem.setAttribute("warning", "The angle between face normal and the vector from cell\n centroid and face center is greater than 179 degree !\nThe grid quality is too poor to use!");
    
   }
  angle_elem.setAttribute("detail", QString("The angle between the face normal and cell centroids provides an indication of mesh isotropy.")+
                          QString( "\nThe lower this number the better the mesh quality.")+
                          QString("\nValues below 100 are desirable.")+
                          QString("\nValues above 150 indicate very poor mesh quality."));
  
  

  QDomElement twist_elem = root.firstChildElement("maxTwist");
   if(twist_elem.isNull()){
    QMessageBox::warning(window(), tr("Quality report"),
                         tr("twist is NUll")
                         );
    return false;
    
  }
  
  if(twist_elem.text().toDouble()>0.8){
     twist_elem.setAttribute("warning", "The maximum twist  greater than 0.8 !\nThe grid quality is too poor to use!");
    
  }
 twist_elem.setAttribute("detail", QString("For non-triangular faces it is possible for the face to be non-planar, i.e. twisted.\n")+
  QString("The twist metric measure the  non-planar component of the face geometry. \n")+
 QString("In other words a value of 0.1 indicates that the face geometry deviates from the planar description by 10 percent.\n")+
 QString("A value below 0.1 is desirable"));
 
 
 QDomElement sheartwist_elem = root.firstChildElement("maxShearTwist");
  if(sheartwist_elem.isNull()){
    QMessageBox::warning(window(), tr("quality report"),
                         tr("sheartwist is NUll")
                         );
    return false;
    
  }
  
 if(sheartwist_elem.text().toDouble()>0.8){
   sheartwist_elem.setAttribute("warning", "The maximum sheartwist  greater than 0.8 !\nThe grid quality is too poor to use!");
  
   }
 sheartwist_elem.setAttribute("detail", QString("For non-triangular faces it is possible for the face to be non-planar, i.e. twisted.\n")+
                              QString("The sheartwist metric measure the  non-planar component of the face geometry. \n")+
                              QString("In other words a value of 0.1 indicates that the face geometry deviates from the planar description by 10 percent.\n")+
                              QString("A value below 0.1 is desirable"));
 
 
 QDomElement quality_elem = root.firstChildElement("quality");
 if(quality_elem.isNull()){
   QMessageBox::warning(window(), tr("Quality report"),
                        tr("quality is NUll")
                        );
   return false;
   
 }
  
 root.removeChild(quality_elem);
 
 
 QDomElement report_elem = root.firstChildElement("qualityReport");
 if(report_elem.isNull()){
    QMessageBox::warning(window(), tr("Quality report"),
                         tr("report is NUll")
                         );
    return false;
    
  }
  
 if(report_elem.text()=="UNUSABLE"){
   report_elem.setAttribute("warning",  "The grid quality is too poor to use!\nSuggest to regenerate grid");

 }
 report_elem.setAttribute("detail", QString("Mesh quality can be defined as: excellent, good, poor, marginal or UNUSABLE. \n")+
                          QString("excellent: volRatio <= 10 && maxCellAngle < 90 || maxTwist <= 0.1 && maxShearTwist <= 0.1 \n")+
                          QString("good: volRatio > 10 || maxCellAngle > 90 || maxTwist > 0.1 || maxShearTwist <= 0.1 \n")+
                          QString("poor: volRatio > 50 || maxCellAngle > 150 || maxTwist > 0.2 || maxShearTwist <= 0.2 \n")+
                          QString("marginal: volRatio > 100 || maxCellAngle > 170 || maxTwist > 0.45 || maxShearTwist <= 0.45 \n")+
                          QString("UNUSABLE: volRatio > 1000 || maxCellAngle > 179 || maxTwist > 0.8 || maxShearTwist <= 0.8 || minVol <0 || convexCell >0 "));

 return true;
}









  
//////////////////////////////////////////////////////////
//  public:
//    QualityDialog(dialog_info *info, QWidget *parent = 0);
//
//  Assembles the dialog shown by the "quality grid " in 'file' menu.
//////////////////////////////////////////////////////////

QualityDialog::QualityDialog(QString filename, QWidget *parent)
  : QDialog(parent)
{
 
  if(! processQualityReport(filename)) return;  
  signalMapper = new QSignalMapper(this);


  QDomElement root = doc.documentElement();
  QDomElement elem = root.firstChildElement();
  if(root.isNull()){
    QMessageBox::warning(window(), tr("quality dialog"),
                         tr("root is NUll")
                         );
    return;
    
 }
  
  QVBoxLayout *qualityLayout = new QVBoxLayout;

  
  
  for(; !elem.isNull(); elem= elem.nextSiblingElement()){
  
    QHBoxLayout* itemLayout = new QHBoxLayout;
    QLabel*   variableLabel = new QLabel(elem.tagName() + tr(": "));
    QLabel*   valueLabel = new QLabel(elem.text());
    itemLayout->addWidget(variableLabel);
    itemLayout->addWidget(valueLabel);
    itemLayout->addSpacing(8);
    if(elem.hasAttribute("warning")){
      QPalette pal(QColor(255,255,255));
      pal.setColor( QPalette::Text, Qt::red );
      pal.setColor( QPalette::Foreground, Qt::red );
      valueLabel->setPalette(pal);
    }
    if(elem.hasAttribute("warning")){
      QLabel* warningLabel = new QLabel(elem.attribute(tr("warning")));
      itemLayout->addWidget(warningLabel);
      QPalette pal(QColor(255,255,255));
      pal.setColor( QPalette::Text, Qt::darkRed );
      pal.setColor( QPalette::Foreground, Qt::darkRed );
      warningLabel->setPalette(pal);
      
    }
    if(elem.hasAttribute("detail")){
      QLabel*   detailLabel = new QLabel(elem.attribute("detail"));
      
      
      QPushButton* detailButton = new QPushButton(tr("Show detail"));
      itemLayout->addWidget(detailButton);
      detailButton->setCheckable(true);
      
      itemLayout->addWidget(detailButton);
      itemLayout->addWidget(detailLabel);
      detailLabel->hide();
      
      buttonList<<detailButton;
      labelList<<detailLabel;
      
      connect(detailButton,SIGNAL(toggled(bool)), signalMapper, SLOT(map())); 
      signalMapper->setMapping(detailButton, int(buttonList.size())-1);
      
    }
    qualityLayout->addLayout(itemLayout);
  }
  connect(signalMapper, SIGNAL(mapped(int)), this, SLOT(buttonToggled(int)));     
  
  setLayout(qualityLayout);
  
}



//////////////////////////////////////////////////////////////////////////////
//  private slots:
//    void buttonToggled();
//
//////////////////////////////////////////////////////////////////////////////

void QualityDialog::buttonToggled(int id)
  {
    if(buttonList[id]->isChecked()){
      buttonList[id]->setText(tr(" Hide detail"));
      labelList[id]->show();
    }else{
      buttonList[id]->setText(tr("Show detail"));
       labelList[id]->hide();
    }
  }
 

