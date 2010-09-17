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




class RegionWindow : public GeneralGroup
{
  Q_OBJECT
  
public:
  RegionWindow( QDomElement& elem, QWidget *parent=0);
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
  void updateShowStatus(const bool& show);
signals:
  void valueChanged(const QTreeWidgetItem*);
  void textChanged();
private:
  
private:
  QGroupBox *stateGroup;
  QListWidget *stateList;
  QStackedWidget *statePages;

  QGroupBox *regionGroup;
  QStackedWidget *regionPages;
  QTreeWidget* tree;
  // QTreeWidgetItem* root;
  vector<Shape*> defaultShapes;
  
  QStackedWidget *myPages;
  QStringList stateFunctions;
  QTextEdit* textEdit;
};







#endif
