/****************************************************************************
** Meta object code from reading C++ file 'bdcndwindow.h'
**
** Created: Mon Sep 14 09:59:15 2009
**      by: The Qt Meta Object Compiler version 59 (Qt 4.4.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "bdcndwindow.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'bdcndwindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 59
#error "This file was generated using the moc from 4.4.3. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_BdCndWindow[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       6,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // signals: signature, parameters, type, tag, flags
      13,   12,   12,   12, 0x05,

 // slots: signature, parameters, type, tag, flags
      22,   12,   12,   12, 0x0a,
      46,   12,   12,   12, 0x08,
      75,   12,   12,   12, 0x08,
      97,   12,   12,   12, 0x08,
     113,   12,   12,   12, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_BdCndWindow[] = {
    "BdCndWindow\0\0closed()\0setCurrent(QModelIndex)\0"
    "updateConditionView(QString)\0"
    "updateConditionView()\0changePage(int)\0"
    "checkStatus()\0"
};

const QMetaObject BdCndWindow::staticMetaObject = {
    { &GeneralWindow::staticMetaObject, qt_meta_stringdata_BdCndWindow,
      qt_meta_data_BdCndWindow, 0 }
};

const QMetaObject *BdCndWindow::metaObject() const
{
    return &staticMetaObject;
}

void *BdCndWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_BdCndWindow))
        return static_cast<void*>(const_cast< BdCndWindow*>(this));
    return GeneralWindow::qt_metacast(_clname);
}

int BdCndWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = GeneralWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: closed(); break;
        case 1: setCurrent((*reinterpret_cast< QModelIndex(*)>(_a[1]))); break;
        case 2: updateConditionView((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 3: updateConditionView(); break;
        case 4: changePage((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 5: checkStatus(); break;
        }
        _id -= 6;
    }
    return _id;
}

// SIGNAL 0
void BdCndWindow::closed()
{
    QMetaObject::activate(this, &staticMetaObject, 0, 0);
}
QT_END_MOC_NAMESPACE
