#ifndef STATEREGION_H
#define STATEREGION_H

#include <QWidget>
#include <QStringList>
#include <QDomDocument>
#include <QModelIndex>
#include <QGroupBox>
#include <QModelIndex>
#include <QStandardItemModel>
#include <QTableView>
#include <QList>
#include <QHBoxLayout>

#include <utility>
#include "pages.h"
#include "grid.h"
#include "fvmadapt.h"

class QListWidget;
class QListWidgetItem;
class QStackedWidget;
class QCheckBox;
class QButtonGroup;
class QStackedLayout;
using namespace std;

//typedef affineMapping affineMapping2;




class RegionWindow : public GeneralWindow
{
  Q_OBJECT
  
public:
  RegionWindow( QDomElement& elem,  QDomElement& theroot,  QWidget *parent=0);
  QTreeWidgetItem* getRoot();
  QString currentText();
public slots:
  void addStateClicked();
  void setItemText(const QString&);
  
  void addShape();
  void updateShape();
  void showData(QTreeWidgetItem* item);
  void assignState();
 
  void next();
  void previous();
  void updateText();
signals:
  void valueChanged(const QTreeWidgetItem*);
private:
  
private:
  QGroupBox *stateGroup;
   
  QListWidget *typesWidget;
  QStackedWidget *pagesWidget;

  QGroupBox *regionGroup;
  QStackedWidget *paraPages;
  QTreeWidget* tree;
  QTreeWidgetItem* root;
  vector<Shape*> defaultShapes;
  
   QStackedWidget *myPages;
  QStringList stateFunctions;
  QTextEdit* textEdit;
};







#endif
