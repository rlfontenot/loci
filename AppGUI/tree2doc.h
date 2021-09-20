#ifndef TREE2DOC_H
#define TREE2DOC_H
#include <QDomDocument>
#include <QTreeWidget>
#include <QTreeWidgetItem>
#include "defines.h"

QDomNode makeElement(QDomDocument& doc, const QTreeWidgetItem* item);
void graft_tree(QDomDocument& doc);
QDomDocument tree2dom(const QTreeWidgetItem* root);
bool anglesBetween(positions3d v1,positions3d v2, double& heading, double& attitude, double& bank);
#endif
