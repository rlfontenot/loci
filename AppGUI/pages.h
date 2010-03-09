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


class GeneralWindow : public QWidget
{
  Q_OBJECT
  
public:
  GeneralWindow(QDomElement& elem, QDomElement& root,  QWidget* parent = 0);
 
 public slots:
 void changeState();
  void updateShowStatus(const bool &);
  signals:
 void updateStatus(const QString&);
  void updateStatusTip(int);
  void stateChanged();
  void componentsChanged();
   void showStatus(const bool &);
  
protected:
  QDomElement myelem;
  QDomElement myroot;
  
  
};



class FloatEdit : public QLineEdit{
  
    Q_OBJECT
  
public:
    FloatEdit(QWidget *parent = 0);
  FloatEdit(double d, QWidget *parent = 0);
  void setBottom(double);
  void setTop(double);
  //   void setDecimals(int);
  void setRange(double minn, double maxx);
   double value();
public slots:
void setValue(int);
  void setValue(double);
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
  void paraChanged();
   
protected:
  
  void keyPressEvent(QKeyEvent *event);
};
  
class FloatSlider : public QWidget{
  
  Q_OBJECT
  
public:
    
  FloatSlider(const QString& title,  QWidget *parent = 0);
  void setRange(double minn, double maxx);
  double value();
public slots:
  void setValue(double);
  void display();
  void undisplay();
private slots:

  signals:
  void valueChanged(double);
protected:
  DoubleSpinBox* edit;
  QLabel* valueLabel;
};

// class FloatSpinBox: public FloatSlider{
  
//   Q_OBJECT
  
// public:
  
//   FloatSpinBox(QDomElement& elem,  QWidget *parent = 0);
//  public slots:
//  void updateCurrent();
// private:
//   QDomElement myelem;
// };

  
class VectSlider:public QGroupBox{
   Q_OBJECT
  
public:
    
   VectSlider( const QString& title, QWidget *parent = 0);
  void setRange(double d1, double d1);
  void setXRange(double minn, double maxx);
  void setYRange(double minn, double maxx);
  void setZRange(double minn, double maxx);
  positions3d value();
 public slots:
 void setValue(const positions3d& p);
private slots:
 void setInfo();
  signals:
void valueChanged(const positions3d&);
private:
  FloatSlider* xedit;
  FloatSlider* yedit;
  FloatSlider* zedit;
 
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
  GeneralGroup(QDomElement& elem, QDomElement& root,  QWidget* parent = 0);
    QString currentText();
public slots:
void changeState();
  void updateChecked();
  void updateShowStatus(const bool&);
  signals:
 void stateChanged();
  void componentsChanged();
  void textChanged(const QString &);
  void showStatus(const bool &);
protected:
  QDomElement myelem;
  QDomElement myroot;

  
};


class OpGroup : public GeneralGroup{
  
  Q_OBJECT
  
public:
  OpGroup(QDomElement& elem,  QDomElement& root,  QWidget *parent = 0);
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
  //void updateChecked();
  void updateSelection(int); //for flags

 
 void setDefault(bool satisfied);//for default condition, when the condition is satisfied

 private:
  
  QList<QLabel*> labels;//for vectors and dvectors
  QList<QWidget*> mfs;//for dvectors
 
};

class AllVBGroup : public GeneralGroup{
  
  Q_OBJECT
  
public:
  AllVBGroup(  QDomElement& elem,  QDomElement& root, bool isVertical = true, QWidget *parent = 0);
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
  StackGroup(  QDomElement& elem, QDomElement& root,  QWidget *parent = 0);
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
class OptionPage : public GeneralGroup{
   Q_OBJECT
public:
   OptionPage(QDomElement &elem,
              QDomElement &root,
              QWidget *parent=0);
  ~OptionPage();
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


class AllVBWindow : public GeneralWindow{
  
  Q_OBJECT
   
public:
  AllVBWindow(  QDomElement& elem,  QDomElement& root,  QWidget *parent = 0);
  
public slots:
bool save();
private:

};

class ChoiceGroup : public GeneralGroup{
  Q_OBJECT
public:
  ChoiceGroup( QDomElement &elem,
                            QDomElement &root,
                           QWidget *parent=0);
  ~ChoiceGroup();
  QString currentText();
 public slots:
 void editButtonPressed();
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


class StackGroup2 : public GeneralGroup{
  
   Q_OBJECT
  
public:
   StackGroup2(  QDomElement& elem,  QDomElement& root,
                 const QStringList& pState,  QWidget *parent = 0);
  QString currentText();
  
  signals:
  void stateChanged(QString);

  
public slots:
void parentStateChanged(QString);
  void changePage(int i);//when buttonGroup clicked
  void updateState(QString);//when children's state changed
  QString myState();
  void  setParentState(const QStringList&);//when parent finish construction
  void add();
  void updateCurrentText();
  void updateComponents();
  void changeState();//try to disable the slot of GeneralGroup
private:
 
  QString name;
  QStringList values;
  int current;
  QButtonGroup *buttonGroup;
  QStackedWidget *pagesWidget;
  QList<StackGroup2*> stacks;
  QList<ChoiceGroup*> conditionalModels;
  QList<QDomElement> conditionalElems;
  QList<ChoiceGroup*> unconditionalModels;
  QStringList parentState;
};






class Page : public GeneralWindow{
  Q_OBJECT
public:
  Page(QDomElement &elem,
       QDomElement &root,
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
