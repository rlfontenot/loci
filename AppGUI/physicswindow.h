#ifndef PHYSICSWINDOW_H
#define PHYSICSWINDOW_H

#include <QWidget>
#include <QStringList>
#include <QDomDocument>
#include <QDomElement>
#include <QModelIndex>
#include <QItemDelegate>
#include <QModelIndex>
#include <QObject>
#include <QSize>
#include <QHeaderView>
#include <QItemSelectionModel>
#include <QStandardItemModel>
#include <QTableView>
#include <QObject>
#include "pages.h"

class QButtonGroup;
class QStackedWidget;
class QRadioButton;
class QSpinBox;
class QStandardItemModel;






class FloatEditDelegate : public QItemDelegate
 {
     Q_OBJECT

 public:
     FloatEditDelegate(QObject *parent = 0);

     QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &option,
                           const QModelIndex &index) const;

     void setEditorData(QWidget *editor, const QModelIndex &index) const;
     void setModelData(QWidget *editor, QAbstractItemModel *model,
                       const QModelIndex &index) const;

     void updateEditorGeometry(QWidget *editor,
         const QStyleOptionViewItem &option, const QModelIndex &index) const;
 };

class CpWindow : public QWidget
{
  Q_OBJECT
  
public:
  CpWindow( QDomElement& elem , QWidget* parent = 0);
  
public slots:
void setNumberOfIntervals(int);
  bool save();
  void replaceButtonToggled(bool);
private:
 
  QRadioButton* replaceButton;
  QRadioButton* augmentButton;

  QRadioButton* shomateButton;
  QRadioButton* polyButton;

  QGroupBox *mGroup;
  FloatEdit* mEdit;
  IntEdit* nEdit;
  QLineEdit* nameEdit;
  QStandardItemModel* model; 
  int numberOfIntervals;
  QTableView* tableView;
  std::vector<double> temperatues;
  //QString unitOfTemperature;
  
  QDomElement myelem;
  
 
  
};










class PhysicsWindow : public GeneralWindow
{
  Q_OBJECT
  
public:
  PhysicsWindow(QDomElement& theelem, QDomElement& myroot, QWidget* parent = 0);
  // QStringList currentState();
  signals:
  //  void componentsChanged();
  void parentStateChanged(QString);
  void stateUpdated(const QStringList&);
private slots:
void updateState(QString);
  void checkStatus();
private:
  QStringList state; 
  QList<StackGroup2*> stacks;
 };

#endif
