#ifndef QTWIDGET_H
#define QTWIDGET_H

#include <QWidget>
#include <QImage>

#include "ntsim.h"

namespace Ui {
  class QtWidget;
}

struct QImage;

class QtWidget : public QWidget
{
  Q_OBJECT

public:
  explicit QtWidget(
    const int width = 600,
    const int height = 400,
    QWidget *parent = 0
  );
  QtWidget(const QtWidget&) = delete;
  QtWidget& operator=(const QtWidget&) = delete;
  ~QtWidget();

  void SetPixel(const int x, const int y, const QColor color);
  void SetResolution(const int width, const int height);

protected:
  void paintEvent(QPaintEvent *);
private:
  Ui::QtWidget *ui;
  QImage m_image;

  NTsim m_sim;
  int m_t; //Timestep = x coordinat
  void Display();

private slots:
  //Increments m_t
  void OnTimer();
};

#endif // QTWIDGET_H
