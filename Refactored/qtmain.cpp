#include "qtwidget.h"
#include <QApplication>

int main(int argc, char *argv[])
{
  QApplication a(argc, argv);
  QtWidget w;
  w.show();
  return a.exec();
}
