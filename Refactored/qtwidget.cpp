#include "qtwidget.h"

#include <QImage>
#include <QPainter>
#include <QPixmap>
#include <QTimer>

#include "ui_qtwidget.h"

QtWidget::QtWidget(
  const int width,
  const int height,
  QWidget *parent
)
  : QWidget(parent),
    ui(new Ui::QtWidget),
    m_image(width,height,QImage::Format_RGB32),
    m_sim{},
    m_t{0}
{
  ui->setupUi(this);

  const int seed = 42;
  const int metacommunity_size = 100;
  const double mutation_rate = 0.1;
  const double selection_strength = 0.1;
  m_sim.setup(
    seed,
    metacommunity_size,
    mutation_rate,
    selection_strength
  );

  Display();

  {
    QTimer * const timer{new QTimer(this)};
    timer->setInterval(100);
    QObject::connect(timer,SIGNAL(timeout()),this,SLOT(OnTimer()));
    timer->start();
  }
}

QtWidget::~QtWidget()
{
  delete ui;
}

void QtWidget::Display()
{
  const auto& cs = m_sim.GetFitnessCategories();
  const int sz{static_cast<int>(cs.size())};
  for (int i=0; i!=sz; ++i)
  {
    const int c = cs[i];
    //if (c < 0) continue;
    //const int r = qRed(m_image.pixel(m_t,i));
    const int w = 32;
    SetPixel(m_t,i,qRgb(128+(c*w),128+(c*w),128+(c*w)));
  }
  repaint();
}

void QtWidget::OnTimer()
{
  //std::clog << ".";
  m_sim.sim_step();
  ++m_t;
  Display();
}

void QtWidget::SetPixel(const int x, const int y, const QColor color)
{
  m_image.setPixel(x,y,color.rgb());
}

void QtWidget::SetResolution(const int width, const int height)
{
  m_image = QImage(width,height,QImage::Format_RGB32);
}

void QtWidget::paintEvent(QPaintEvent *)
{
  QPainter painter(this);
  painter.drawPixmap(
    this->rect(),
    QPixmap::fromImage(m_image)
  );
}

/*

  NTsim sim;

  {
    const auto& cs = sim.GetFitnessCategories();
    std::copy(std::begin(cs),std::end(cs),
      std::ostream_iterator<long>(std::cout," ")
    );
    std::cout << std::endl;
  }

  for (int i=0; i!=100; ++i)
  {
    for (int j=0; j!=100; ++j)
    {
      sim.sim_step();
    }
    {
      const auto& cs = sim.GetFitnessCategories();
      std::copy(std::begin(cs),std::end(cs),std::ostream_iterator<long>(std::cout," "));
      std::cout << std::endl;
    }
  }

*/
