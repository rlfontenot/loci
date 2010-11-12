#ifndef PAGES_H
#define PAGES_H
#include <QWidget>
#include <QDomDocument>
#include <QString>
#include <QObject>
#include <QDomElement>
#include <QLineEdit>
#include <QGroupBox>
#include <QStringList>
#include <QDoubleValidator>
#include <QSignalMapper>
#include <QList>
#include <QComboBox>
#include <QPointer>
#include <iostream>
#include <QString>
#include <QItemDelegate>
#include <QDoubleSpinBox>
#include "grid.h"
class QPushButton;
class QButtonGroup;
class QLabel;
class QListWidget;
class QStackedWidget;
class QListWidgetItem;
class QCheckBox;





bool conditionIsSatisfied(const QDomElement& theroot, const QString& condition);






class DoubleEdit : public QLineEdit{
  
  Q_OBJECT
  
  public:
  DoubleEdit(QWidget *parent = 0);
  DoubleEdit(double d, QWidget *parent = 0);
  void setBottom(double);
  void setTop(double);
  void setRange(double minn, double maxx);
  double value();
public slots:
  void mapValue(int);
  void setValue(double);
private slots:
  void changeValue(const QString&);
signals:
  void valueChanged(double);
private:
  QDoubleValidator* validator;
};



class DoubleSpinBox:public QDoubleSpinBox{
  Q_OBJECT
  public:
  DoubleSpinBox(QWidget* parent = 0);
signals:
  void paraChanged();//this signal is emitted when the step value or decimal value changed
protected:
  void keyPressEvent(QKeyEvent *event);
};



class LabeledDoubleSpBox : public QWidget{
  
  Q_OBJECT
  
  public:
  
  LabeledDoubleSpBox(const QString& title,  QWidget *parent = 0);
  void setRange(double minn, double maxx);
  double value();
  inline double singleStep(){return edit->singleStep();}
  inline int decimals(){return edit->decimals();}
  inline void setSingleStep(double d){edit->setSingleStep(d);}
  inline void setDecimals(int n){edit->setDecimals(n);} 
public slots:
  void setValue(double);
  void display();
  void undisplay();
private slots:
  
signals:
  void valueChanged(double);
  void paraChanged();
protected:
  DoubleSpinBox* edit;
  QLabel* valueLabel;
};


  
class VectSpBox:public QGroupBox{
  Q_OBJECT
  
  public:
    
  VectSpBox( const QString& title, QWidget *parent = 0);
  void setRange(double d1, double d1);
  void setXRange(double minn, double maxx);
  void setYRange(double minn, double maxx);
  void setZRange(double minn, double maxx);
  positions3d value();
public slots:
  void setValue(const positions3d& p);
private slots:
  void setInfo();
  void setPara();
signals:
  void valueChanged(const positions3d&);
private:
  LabeledDoubleSpBox* xedit;
  LabeledDoubleSpBox* yedit;
  LabeledDoubleSpBox* zedit;
 
};


class IntEdit : public QLineEdit{
  
  Q_OBJECT
  
  public:
  IntEdit(QWidget *parent = 0);
  IntEdit(int d, QWidget *parent = 0);
  void setBottom(int);
  void setTop(int);
  void setRange(int minn, int maxx);
 
  int value();
public slots:
  void setValue(int);
private slots:
  void changeValue(const QString&);
signals:
  void valueChanged(int);
private:
  QIntValidator* validator;
};


class GeneralGroup : public QGroupBox
{
  Q_OBJECT
  
  public:
  GeneralGroup(QDomElement& elem,QWidget* parent = 0);
  QString currentText();
public slots:
  void updateComponents(); //for dvector
  void changeState();
  void updateChecked();
  void updateShowStatus(const bool&);
  void gotoUnfinishedPage(){};
signals:
  void updateStatus(const QString&);
  void updateStatusTip(int);
  void stateChanged();
  void componentsChanged();
  void textChanged(const QString &);
  void showStatus(const bool &);
protected:
  QDomElement myelem;

};




class VarGBox : public GeneralGroup{
  
  Q_OBJECT
  
  public:
  VarGBox(QDomElement& elem, QWidget *parent = 0);
  QString currentText();
signals:
 
public slots:
  void updateComponents(); //for dvector
  void changeState();
private slots:
  void updateCurrentX(const QString&);
  void updateCurrentY(const QString&);
  void updateCurrentZ(const QString&);
  void updateCurrentUnit(const QString&);
  void updateUnitList(const QString&);
  void updateLabels(int );
  void updateCurrent(const QString&);
  void updateShowStatus(const bool&);
  void updateSelection(int); //for flags
  void setDefault(bool satisfied);//for default condition, when the condition is satisfied

private:
 /*!
   for vectors and dvectors 
 */ 
  QList<QLabel*> labels;
  QList<QWidget*> mfs;
 
};

class AllGroup : public GeneralGroup{
  
  Q_OBJECT
  
  public:
  AllGroup(  QDomElement& elem,  bool isVertical = true, QWidget *parent = 0);
  QString currentText();
  
signals:
public slots:
  void childrenClicked(int i);//for 2of3 element
  void updateCurrentText();
 
private:
  //for 2of3 element
  QList<QWidget*> objs;
  QSignalMapper *signalMapper;
};


  

class StackGroup : public GeneralGroup{
  
  Q_OBJECT
  
  public:
  StackGroup(  QDomElement& elem, QWidget *parent = 0);
  QString currentText();
  void clearText();
signals:
  
public slots:
  void updateCurrentText();
  void changePage(int);
private:
  
  QComboBox *typesWidget;
  QStackedWidget *pagesWidget;
  QStringList whatsThisList;
  QStringList toolTipList;

};






//class for a panel
class VarPanel : public GeneralGroup{
  Q_OBJECT
  public:
  VarPanel(QDomElement &elem,
          QWidget *parent=0);
  ~VarPanel();
  QString currentText();
  void clearText();
  
public slots:
  void advancedButtonClicked(int);
  void updateCurrentText();
  void updateShowStatus(const bool&);
signals:

  void updateStatusTip(int);
  void updateStatus( QString);

protected:
  QList<QWidget*> myAdvancedGroup;
  QButtonGroup* buttonGroup;
  QLabel* currentLabel;
  bool showStatus;
};



//for physicswindow and solverwindow


class AllVBWindow : public GeneralGroup{
  
  Q_OBJECT
   
  public:
  AllVBWindow(  QDomElement& elem, QWidget *parent = 0);
  
public slots:
  bool save();
private:

};

class ChoiceGroup : public GeneralGroup{
  Q_OBJECT
  public:
  ChoiceGroup( QDomElement &elem,
               QWidget *parent=0);
  ~ChoiceGroup();
  QString currentText();
public slots:
  void editButtonPressed();
  void changeState();
  void update(int);
  void updateCurrentText();
  void updateShowStatus(const bool&);
signals:
private:
  
  QPushButton* editButton;
  QWidget* editGroup;
  QComboBox* comboBox;
  QStringList editItems;
  QStringList whatsThisList;
  QStringList toolTipList;
  bool showStatus;
 
};


class StateStackGroup : public GeneralGroup{
  
  Q_OBJECT
  
  public:
  StateStackGroup(  QDomElement& elem, 
                QWidget *parent = 0);
  QString currentText();
  
signals:
   /*!
     This signal is emited when my state changed. To tell \a parent update its state
   */ 
  void stateChanged(QString);

  
public slots:
  void changePage(int i);//when buttonGroup clicked
  QString myState();
  void add();
  void updateCurrentText();
  

private:
  QButtonGroup *buttonGroup;
  QStackedWidget *pagesWidget;
  
};






class VarPage : public GeneralGroup{
  Q_OBJECT
  public:
  VarPage(QDomElement &elem,
       QWidget *parent=0);
   
private slots:
  void updateCurrentText();
signals:
  void textChanged(const QString &);
protected:
};

// Class that delegates for the boundary visibility model
class showDelegate: public QItemDelegate
{
  Q_OBJECT
  public:
  QWidget* createEditor(QWidget* parent, const QStyleOptionViewItem &option,
			const QModelIndex &index) const;
  showDelegate(QObject* parent=0):QItemDelegate(parent){};
  
};

class colorDelegate: public QItemDelegate
{
  
  Q_OBJECT
  public:
  bool  editorEvent(QEvent *event, QAbstractItemModel *model,
                    const QStyleOptionViewItem &option,
                    const QModelIndex &index);
  colorDelegate(QObject* parent=0):QItemDelegate(parent){};
  
};





#endif

