#ifndef QUALITYDIALOG_H
#define QUALITYDIALOG_H

#include <QDialog>
#include "grid.h"
#include <QDomDocument>
#include <QWidget>
#include <QObject>
#include <QList>
#include <QSignalMapper>

class QTextEdit;
class QLabel;
class QPushButton;



class QualityDialog: public QDialog
{
  Q_OBJECT

public:
  QualityDialog(QString filename, QWidget *parent = 0);
 
private slots:
void buttonToggled(int);
private:
  bool processQualityReport(QString filename);
   
  QDomDocument doc;
  QSignalMapper* signalMapper;
  QList<QPushButton*> buttonList;
  QList<QLabel*> labelList;  
};

#endif
