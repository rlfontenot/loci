/****************************************************************************
** Meta object code from reading C++ file 'initcndwindow.h'
**
** Created: Mon Sep 14 09:59:16 2009
**      by: The Qt Meta Object Compiler version 59 (Qt 4.4.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "initcndwindow.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'initcndwindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 59
#error "This file was generated using the moc from 4.4.3. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_InitCndWindow[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       2,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
      32,   15,   14,   14, 0x08,
      78,   14,   14,   14, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_InitCndWindow[] = {
    "InitCndWindow\0\0current,previous\0"
    "changePage(QListWidgetItem*,QListWidgetItem*)\0"
    "checkStatus()\0"
};

const QMetaObject InitCndWindow::staticMetaObject = {
    { &GeneralWindow::staticMetaObject, qt_meta_stringdata_InitCndWindow,
      qt_meta_data_InitCndWindow, 0 }
};

const QMetaObject *InitCndWindow::metaObject() const
{
    return &staticMetaObject;
}

void *InitCndWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_InitCndWindow))
        return static_cast<void*>(const_cast< InitCndWindow*>(this));
    return GeneralWindow::qt_metacast(_clname);
}

int InitCndWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = GeneralWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: changePage((*reinterpret_cast< QListWidgetItem*(*)>(_a[1])),(*reinterpret_cast< QListWidgetItem*(*)>(_a[2]))); break;
        case 1: checkStatus(); break;
        }
        _id -= 2;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
