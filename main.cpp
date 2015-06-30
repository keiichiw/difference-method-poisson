#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <Eigen/Sparse>
using namespace std;

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> Tri;
typedef Eigen::VectorXd VecXd;

struct Poisson {

  const int dirx[4] = {0, 0, 1, -1};
  const int diry[4] = {1, -1, 0, 0};

  int N;
  int M;
  int SIZE;
  double dx;
  double dy;

  Poisson(int N, int M) : N(N), M(M) {
    SIZE = N * (M - 1);
    dx = 1.0 / N;
    dy = 1.0 / M;
  }

  double f(double x, double y) {
    return x + y;
  }
  double u_xy(double x, double y) {
    return (pow(x, 3.0) + pow(y, 3.0)) / 6.0;
  }

  bool has_Dirichlet(int i, int j) {
    return i == 0 || i == N || j == 0;
  }

  double dirichlet(double x, double y) {
    return (pow(x, 3) + pow(y, 3)) / 6.0;
  }

  bool has_Neumann(int , int j) {
    return j == M;
  }

  double neumann(double ) {
    return 0.5;
  }

  int vec_nth(int i, int j) {
    assert(!has_Dirichlet(i, j));
    return (N-1) * (j - 1) + i - 1;
  }

  bool valid(int x, int y) {
    return 0 <= x && x <= N && 0 <= y && y <= M;
  }

  VecXd expect() {
    VecXd x = VecXd::Zero(SIZE);
    int row = 0;

    for (int j = 1; j <= M; ++j) {
      for (int i = 1; i < N; ++i) {
        x(row) = u_xy(i * dx, j * dy);
        row++;
      }
    }

    return x;
  }

  tuple<SpMat, VecXd> fdm() {

    vector<Tri> vList;
    SpMat A(SIZE, SIZE);
    VecXd b = VecXd::Zero(SIZE);

    int row = 0;

    for (int j=1; j <= M; ++j) {
      for (int i = 1; i < N; ++i) {

        if (has_Neumann(i, j)) {

          vList.push_back(Tri(row, row, 1.0));
          vList.push_back(Tri(row, vec_nth(i, j-1), -1.0));
          b(row) = dy * neumann(i * dx);

        } else {

          vList.push_back(Tri(row, row, 4.0));
          b(row) += dx * dx * f(i * dx, j * dy);

          for (int k = 0; k < 4; ++k) {
            int nx = i + dirx[k];
            int ny = j + diry[k];
            if (!valid(nx, ny)) continue;

            if (has_Dirichlet(nx, ny)) {
              b(row) += dirichlet(nx * dx, ny * dy);
            } else {
              vList.push_back(Tri(row, vec_nth(nx, ny), -1.0));
            }
          }
        }
        row++;
      }
    }

    A.setFromTriplets(vList.begin(), vList.end());

    return make_tuple(A, b);
  }

  // return (x, iteration count)
  pair<VecXd, int> cg(SpMat A, VecXd b) {

    VecXd x = VecXd::Zero(SIZE);
    VecXd r = b - A * x;
    VecXd p = r;

    int k = 0;


    while (1) {
      auto ax_b = A*x - b;
      if (ax_b.norm() / b.norm() < 1e-12)
        break;

      double alpha = (r.transpose() * r);
      alpha /= (p.transpose() * (A * p));

      x += alpha * p;
      VecXd r1 = r - ((alpha * A) * p);

      double beta = r1.transpose() * r1;
      beta /= r.transpose() * r;

      p = r1 + beta * p;
      r = r1;
      ++k;
    }

    return make_pair(x, k);
  }
};

// 反復回数
void iter_count() {
  ofstream ofs("./data/iter.dat");
  for (int n = 4; n <= 100; n += 4) {
    Poisson p = Poisson(n, n);
    auto t = p.fdm();
    auto result = p.cg(get<0>(t), get<1>(t));

    ofs << 1.0 / n << " " << result.second << endl;
  }
}

// 適合性
void compatibility() {
  ofstream ofs("./data/comp.dat");
  for (int n = 4; n <= 100; n += 4) {
    Poisson p = Poisson(n, n);
    auto t = p.fdm();
    auto A = get<0>(t);
    auto b = get<1>(t);
    auto x = p.expect();
    auto r = (A * x) - b;
    double v = r.cwiseAbs().maxCoeff();
    ofs << n << " " << v << endl;
  }
}

// 収束性
void convergence() {
  ofstream ofs("./data/conv.dat");
  for (int n = 4; n <= 100; n += 4) {
    Poisson p = Poisson(n, n);
    auto t = p.fdm();
    auto A = get<0>(t);
    auto b = get<1>(t);
    auto ex = p.expect();
    auto x = p.cg(A, b).first;
    double v = (ex - x).cwiseAbs().maxCoeff();
    ofs << n << " " << v << endl;
  }
}

// 誤差ノルム
void error_norm() {
  ofstream ofs("./data/error.dat");
  for (int n = 4; n <= 100; n += 4) {
    Poisson p = Poisson(n, n);
    auto t = p.fdm();
    auto A = get<0>(t);
    auto b = get<1>(t);
    auto ex = p.expect();
    auto x = p.cg(A, b).first;
    double norm = (ex - x).norm();
    ofs << 1.0 / n << " " << norm << endl;
  }
}


int main(){

  cerr << "iteration count" << endl;
  iter_count();
  cerr << "compatibility" << endl;
  compatibility();
  cerr << "convergence" << endl;
  convergence();
  cerr << "error norm" << endl;
  error_norm();

}
