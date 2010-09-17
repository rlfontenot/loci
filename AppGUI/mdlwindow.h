#ifndef MDLWINDOW_H
#define MDLWINDOW_H

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
#include "pages.h"

class QButtonGroup;
class QStackedWidget;
class QRadioButton;
class QSpinBox;
class QStandardItemModel;
class QHBoxLayout;

class DoubleEditDelegate : public QItemDelegate
 {
     Q_OBJECT

 public:
     DoubleEditDelegate(QObject *parent = 0);

     QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &option,
                           const QModelIndex &index) const;

     void setEditorData(QWidget *editor, const QModelIndex &index) const;
     void setModelData(QWidget *editor, QAbstractItemModel *model,
                       const QModelIndex &index) const;

     void updateEditorGeometry(QWidget *editor,
         const QStyleOptionViewItem &option, const QModelIndex &index) const;
 };




class CpGroup : public QGroupBox
{
  Q_OBJECT
  
public:
  CpGroup(QWidget* parent = 0);
  QString text();
  void clear();
public slots:
void setNumberOfIntervals(int);
  signals:
void textChanged();
private:
  IntEdit* nEdit;
  QRadioButton* shomateButton;
  QRadioButton* polyButton;
  QStandardItemModel* model; 
  int numberOfIntervals;
  QTableView* tableView;
  // std::vector<double> temperatues;
    
};


class ListGroup:public  QGroupBox
{
  Q_OBJECT
  
public:
  ListGroup(const QString& title=tr("please specify theta_v:"), QWidget* parent = 0);
  
public slots:
void remove();
  void add();
  void clear();
  QString text();
  signals:
  void textChanged();
private:
  QLabel *display;
  DoubleEdit* fEdit;
  QPushButton *removeButton;
};


class SpeciesGroup : public QGroupBox
{
  Q_OBJECT
  
public:
  SpeciesGroup(const QString& title=tr("please specify species:"), QWidget* parent = 0);
  QString text();
  QString getSpecies();
public slots:
//  bool save();
  void replaceBGroupClicked(int);
  void changePage(int);
  void nameBGroupClicked(int);
  void setText();
  void clear();
  signals:
  void accept();
private:
  QButtonGroup* nameBGroup;
  QButtonGroup* replaceBGroup;
  QButtonGroup* thetaCpBGroup;

  QGroupBox *mGroup;
  QGroupBox *mfGroup;
  QGroupBox* sGroup;
  QGroupBox* tGroup ;
  QGroupBox* pGroup;
  QGroupBox*  hGroup ;
  QGroupBox  *thetaCpGroup;
  DoubleEdit* mEdit;
  DoubleEdit* hrefEdit;
  DoubleEdit* srefEdit;
  DoubleEdit* trefEdit;
  DoubleEdit* prefEdit;
  DoubleEdit* mfEdit;
  
  ListGroup* thetavGroup;
  CpGroup* cpGroup;

  QStackedWidget* pagesWidget;
  QLineEdit* nameEdit;

  QLabel* display;
    
};

class CpWindow:public QGroupBox
{
  Q_OBJECT
  public:
  CpWindow(const QString& title, const QString&, QWidget* parent = 0);
  QString text();
public slots:
  bool save();
  void replaceBGroupClicked(int);
  void nameBGroupClicked(int);
  void setText();
private:
  QButtonGroup* nameBGroup;
  QButtonGroup* replaceBGroup;
  
  QGroupBox *mGroup;
  DoubleEdit* mEdit;
  CpGroup* cpGroup;
  QLineEdit* nameEdit;
  QString prefix;
  QLabel* display;




  
};



class KFKC: public QGroupBox
{
  Q_OBJECT
  
public:
  KFKC(const QString& title, QWidget* parent = 0);
  public:
  QString text();
  void checkTher();
   void clear();
  signals:
  void textChanged();
private slots:
void arrToggled(bool);
 
private:
  QRadioButton* ther;
  QRadioButton* arr;
  DoubleEdit* f1;
  DoubleEdit* f2;
  DoubleEdit* f3;
};

class NumberString: public QGroupBox
{
  Q_OBJECT
  
public:
  NumberString(const QString& title, bool floatAllowed = false, QWidget* parent = 0);
  QString text();
  QString cleanText();
  signals:
  void textChanged();
public slots:
  void clear();
private:
  QWidget* edit;
  bool floatAllowed;
  
};
class ExpNu:public QGroupBox
{
  Q_OBJECT
  
  public:
  ExpNu(const QStringList& sp, QWidget* parent = 0);
  
  QString text();
  void clear();
signals:
  void textChanged();
private slots:
private:
  QStringList species;
  QList<NumberString*> exp_nu;
 
};

class RateModifier:public QGroupBox
{
  Q_OBJECT
  
  public:
  RateModifier(const QString& vname, int n, QWidget* parent = 0);
  
  QString text();
  void clear();
signals:
  void textChanged();
private slots:
private:
  QString name;
  int numValues;
  QList<DoubleEdit*> edits;
 
};  


  
class Reactants: public QGroupBox
{
  Q_OBJECT
  
public:
  Reactants(const QStringList& sp, const QString& title, QWidget* parent = 0);
  QString text();
  signals:
  void textChanged();
 public slots:
 void clear();
private:
  QList<NumberString*> species;
};



class ReactionGroup: public QGroupBox
{
  Q_OBJECT
  
public:
  ReactionGroup(const QStringList& sps, const QString& title=tr("please specify reaction:"), QWidget* parent = 0);
  QString text();
 
  
   signals:
  void accept();
public slots:
void setText();
  void directionChanged();
  void clear();
private:
 
  QCheckBox* mbody;
  Reactants* reactants;
  Reactants* products;
  QComboBox* direction;
  QLabel* display;
  KFKC* kf;
  KFKC* kc;
  ExpNu* expnu;
  NumberString* minMf;
  RateModifier* mdf;
    
};

class ReactionWindow: public QGroupBox
{

  Q_OBJECT
  
public:
  ReactionWindow(const QStringList& species, const QString& title=tr("please specify reaction:"), QWidget* parent = 0);
  QStringList getSpecies();
  QString text();
  public slots:
  void add();
  void clear();
  void resetSpecies(const QStringList& sps);
private:
  QStringList species;
  ReactionGroup* reactionPage;
  QListWidget *reactionList;
  QHBoxLayout* mainLayout ;
};



class SpeciesWindow:public QGroupBox
{
  Q_OBJECT
public:
  SpeciesWindow(const QString& title ="please specify species: ", QWidget* parent = 0);
  QStringList getSpecies();
  QString text();
public slots:
void add();
  void clear();

  
private:
  QListWidget *speciesList;
  SpeciesGroup* speciesPage; 
  QStringList species;
};




  
class ChemistryMdl: public QWidget
{
  Q_OBJECT
  public:
  ChemistryMdl(bool noReaction = false, QWidget* parent = 0);
  bool save();
  void hideReaction();
public slots:
  void changePage(int);
  
private:
  QStackedWidget* center;
  QRadioButton* reactionButton;
  SpeciesWindow* speciesWindow;
  ReactionWindow* reactionWindow;
  
};

#endif
