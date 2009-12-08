#include <math.h>
#include <gsl/gsl_errno.h>
#include "bicycle.h"
//#include <stdio.h>
//#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_odeiv.h>


void eval_qdots(double qd[], const double q[], const double u[], const Bicycle_t *bike)
{
  qd[0] = (-sin(q[2])*u[0] + cos(q[2])*u[2]) / cos(q[1]);
  qd[1] = cos(q[2])*u[0] + sin(q[2])*u[2];
  qd[2] = sin(q[2])*tan(q[1])*u[0] + u[1] - cos(q[2])*tan(q[1])*u[2];
  qd[3] = u[3] - u[1];
  qd[4] = u[4] - u[2];
  qd[5] = sin(q[4])*u[0] - cos(q[4])*u[1] + u[5];
  qd[6] = -bike->rr*cos(q[0])*(qd[2] + qd[3]) - bike->rrt*(sin(q[0])*qd[1] +
          cos(q[0])*cos(q[1])*(qd[2] + qd[3]));
  qd[7] = -bike->rr*sin(q[0])*(qd[2] + qd[3]) + bike->rrt*(cos(q[0])*qd[1] +
          sin(q[0])*cos(q[1])*(qd[2] + qd[3]));
} // eval_qdots


void Bicycle_to_Bicycle_benchmark(Bicycle_benchmark_t *bike_b, const Bicycle_t *bike)
{
} // Bicycle_to_Bicycle_benchmark

void Bicycle_benchmark_to_Bicycle(Bicycle_t *bike, const Bicycle_benchmark_t *bike_b)
{
} // Bicycle_benchmark_to_Bicycle


void eval_dependent_speeds(double u[], double B[3][6], const int ind[], const int dep[])
{
  double Bi[3][3], Bd[3][3], Bd_adj[3][3], Bd_det, ui[3], ud[3];

  // Independent speeds
  ui[0] = u[ind[0]];
  ui[1] = u[ind[1]];
  ui[2] = u[ind[2]];

  // Dependent speeds
  ud[0] = u[dep[0]];
  ud[1] = u[dep[1]];
  ud[2] = u[dep[2]];

  // Form -Bi
  Bi[0][0] = -B[0][ind[0]];
  Bi[0][1] = -B[0][ind[1]];
  Bi[0][2] = -B[0][ind[2]];

  Bi[1][0] = -B[1][ind[0]];
  Bi[1][1] = -B[1][ind[1]];
  Bi[1][2] = -B[1][ind[2]];

  Bi[2][0] = -B[2][ind[0]];
  Bi[2][1] = -B[2][ind[1]];
  Bi[2][2] = -B[2][ind[2]];

  // Form Bd
  Bd[0][0] = B[0][dep[0]];
  Bd[0][1] = B[0][dep[1]];
  Bd[0][2] = B[0][dep[2]];

  Bd[1][0] = B[1][dep[0]];
  Bd[1][1] = B[1][dep[1]];
  Bd[1][2] = B[1][dep[2]];

  Bd[2][0] = B[2][dep[0]];
  Bd[2][1] = B[2][dep[1]];
  Bd[2][2] = B[2][dep[2]];

  // Forming Bd_adj and Bd_det
  madj3(Bd_adj, Bd);
  mdet3(&Bd_det, Bd);

  // Bd_adj * (-Bi), store the result in Bd instead of creating a new matrix to
  // store the result
  mm3(Bd, Bd_adj, Bi);
  mv3(ud, Bd, ui);

  // Form the dependent speeds
  u[dep[0]] = ud[0] / Bd_det;
  u[dep[1]] = ud[1] / Bd_det;
  u[dep[2]] = ud[2] / Bd_det;
} // eval_dependent_speeds

/*
 * 3x3 Matrix - Matrix multiplication
 */
void mm3(double AB[3][3], double A[][3], double B[][3])
{
  unsigned short i, j, k;
  double sum;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      sum = 0.0;
      for (k = 0; k < 3; k++) {
        sum += A[i][k] * B[k][j];
      } // for k
      AB[i][j] = sum;
    } // for j
  } // for i
} // mm3

/*
 * 3x3 Matrix - 3x1 Vector multiplication
 */
void mv3(double Ab[3], double A[3][3], double b[3])
{
  unsigned short i, j;
  double sum;
  for (i = 0; i < 3; i++) {
    sum = 0.0;
    for (j = 0; j < 3; j++) {
      sum += A[i][j]*b[j];
    } // for j
    Ab[i] = sum;
  } // for i
} // mv 3

/*
 * 3x3 Matrix adjugate
 */
void madj3(double adj[3][3], double m[3][3])
{
  adj[0][0] = m[1][1]*m[2][2] - m[1][2]*m[2][1];
  adj[0][1] = m[0][2]*m[2][1] - m[0][1]*m[2][2];
  adj[0][2] = m[0][1]*m[1][2] - m[0][2]*m[1][1];

  adj[1][0] = m[1][2]*m[2][0] - m[1][0]*m[2][2];
  adj[1][1] = m[0][0]*m[2][2] - m[0][2]*m[2][0];
  adj[1][2] = m[0][2]*m[1][0] - m[0][0]*m[1][2];

  adj[2][0] = m[1][0]*m[2][1] - m[1][1]*m[2][0];
  adj[2][1] = m[0][1]*m[2][0] - m[0][0]*m[2][1];
  adj[2][2] = m[0][0]*m[1][1] - m[0][1]*m[1][0];
} // madj3

/*
 * Determinant of a 3x3 matrix
 */
void mdet3(double *det, double m[3][3])
{
  *det = m[0][0]*(m[1][1]*m[2][2] - m[1][2]*m[2][1]) +
         m[0][1]*(m[1][2]*m[2][0] - m[1][0]*m[2][2]) +
         m[0][2]*(m[1][0]*m[2][1] - m[1][1]*m[2][0]);
} // mdet3

void form_constraint_matrix(double B[3][6], const double q[], const Bicycle_t *bike)
{
  double den = sqrt(pow(sin(q[1])*sin(q[4]) - cos(q[1])*sin(q[2])*cos(q[4]), 2)
               + pow(cos(q[1])*cos(q[3]), 2)),
         g31 = (sin(q[1])*sin(q[4]) - cos(q[1])*sin(q[2])*cos(q[5])) / den,
         g33 = cos(q[1])*cos(q[2]) / den;

  B[0][0] = cos(q[4])*sin(q[4])*(g33*bike->rf + cos(q[1])*cos(q[2])*bike->rft);
  B[0][1] = bike->ls +
            pow(sin(q[4]), 2)*(g33*bike->rf + cos(q[1])*cos(q[2])*bike->rft);
  B[0][2] = bike->rrt*sin(q[2]);
  B[0][3] = - cos(q[2])*(bike->rr + cos(q[1])*bike->rrt);
  B[0][4] = - (sin(q[1])*bike->rft + sin(q[4])*(g31*bike->rf + bike->lf));
  B[0][5] = cos(q[4])*(g33*bike->rf + cos(q[1])*cos(q[2])*bike->rft);

  B[1][0] = - bike->ls + cos(q[2])*(bike->rr + cos(q[1])*bike->rrt)
            - pow(cos(q[4]), 2)*(g33*bike->rf + cos(q[1])*cos(q[2])*bike->rft);
  B[1][1] = -cos(q[4])*sin(q[4])*(g33*bike->rf + cos(q[1])*cos(q[2])*bike->rft);
  B[1][2] = bike->lr + sin(q[2])*(bike->rr + cos(q[1])*bike->rrt);
  B[1][3] = 0.0;
  B[1][4] = cos(q[4])*(bike->lf + g31*bike->rf) + cos(q[1])*sin(q[2])*bike->rft;
  B[1][5] = 0.0;

  B[2][0] = 0.0;
  B[2][1] = 0.0;
  B[2][2] = 0.0;
  B[2][3] = 0.0;
  B[2][4] = 0.0;
  B[2][5] = 0.0;
} // form_constraint_matrix
