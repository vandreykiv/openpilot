#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with sympy 1.9                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_3484915228138340995) {
   out_3484915228138340995[0] = delta_x[0] + nom_x[0];
   out_3484915228138340995[1] = delta_x[1] + nom_x[1];
   out_3484915228138340995[2] = delta_x[2] + nom_x[2];
   out_3484915228138340995[3] = delta_x[3] + nom_x[3];
   out_3484915228138340995[4] = delta_x[4] + nom_x[4];
   out_3484915228138340995[5] = delta_x[5] + nom_x[5];
   out_3484915228138340995[6] = delta_x[6] + nom_x[6];
   out_3484915228138340995[7] = delta_x[7] + nom_x[7];
   out_3484915228138340995[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_101761715903343570) {
   out_101761715903343570[0] = -nom_x[0] + true_x[0];
   out_101761715903343570[1] = -nom_x[1] + true_x[1];
   out_101761715903343570[2] = -nom_x[2] + true_x[2];
   out_101761715903343570[3] = -nom_x[3] + true_x[3];
   out_101761715903343570[4] = -nom_x[4] + true_x[4];
   out_101761715903343570[5] = -nom_x[5] + true_x[5];
   out_101761715903343570[6] = -nom_x[6] + true_x[6];
   out_101761715903343570[7] = -nom_x[7] + true_x[7];
   out_101761715903343570[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_7962091599851903711) {
   out_7962091599851903711[0] = 1.0;
   out_7962091599851903711[1] = 0;
   out_7962091599851903711[2] = 0;
   out_7962091599851903711[3] = 0;
   out_7962091599851903711[4] = 0;
   out_7962091599851903711[5] = 0;
   out_7962091599851903711[6] = 0;
   out_7962091599851903711[7] = 0;
   out_7962091599851903711[8] = 0;
   out_7962091599851903711[9] = 0;
   out_7962091599851903711[10] = 1.0;
   out_7962091599851903711[11] = 0;
   out_7962091599851903711[12] = 0;
   out_7962091599851903711[13] = 0;
   out_7962091599851903711[14] = 0;
   out_7962091599851903711[15] = 0;
   out_7962091599851903711[16] = 0;
   out_7962091599851903711[17] = 0;
   out_7962091599851903711[18] = 0;
   out_7962091599851903711[19] = 0;
   out_7962091599851903711[20] = 1.0;
   out_7962091599851903711[21] = 0;
   out_7962091599851903711[22] = 0;
   out_7962091599851903711[23] = 0;
   out_7962091599851903711[24] = 0;
   out_7962091599851903711[25] = 0;
   out_7962091599851903711[26] = 0;
   out_7962091599851903711[27] = 0;
   out_7962091599851903711[28] = 0;
   out_7962091599851903711[29] = 0;
   out_7962091599851903711[30] = 1.0;
   out_7962091599851903711[31] = 0;
   out_7962091599851903711[32] = 0;
   out_7962091599851903711[33] = 0;
   out_7962091599851903711[34] = 0;
   out_7962091599851903711[35] = 0;
   out_7962091599851903711[36] = 0;
   out_7962091599851903711[37] = 0;
   out_7962091599851903711[38] = 0;
   out_7962091599851903711[39] = 0;
   out_7962091599851903711[40] = 1.0;
   out_7962091599851903711[41] = 0;
   out_7962091599851903711[42] = 0;
   out_7962091599851903711[43] = 0;
   out_7962091599851903711[44] = 0;
   out_7962091599851903711[45] = 0;
   out_7962091599851903711[46] = 0;
   out_7962091599851903711[47] = 0;
   out_7962091599851903711[48] = 0;
   out_7962091599851903711[49] = 0;
   out_7962091599851903711[50] = 1.0;
   out_7962091599851903711[51] = 0;
   out_7962091599851903711[52] = 0;
   out_7962091599851903711[53] = 0;
   out_7962091599851903711[54] = 0;
   out_7962091599851903711[55] = 0;
   out_7962091599851903711[56] = 0;
   out_7962091599851903711[57] = 0;
   out_7962091599851903711[58] = 0;
   out_7962091599851903711[59] = 0;
   out_7962091599851903711[60] = 1.0;
   out_7962091599851903711[61] = 0;
   out_7962091599851903711[62] = 0;
   out_7962091599851903711[63] = 0;
   out_7962091599851903711[64] = 0;
   out_7962091599851903711[65] = 0;
   out_7962091599851903711[66] = 0;
   out_7962091599851903711[67] = 0;
   out_7962091599851903711[68] = 0;
   out_7962091599851903711[69] = 0;
   out_7962091599851903711[70] = 1.0;
   out_7962091599851903711[71] = 0;
   out_7962091599851903711[72] = 0;
   out_7962091599851903711[73] = 0;
   out_7962091599851903711[74] = 0;
   out_7962091599851903711[75] = 0;
   out_7962091599851903711[76] = 0;
   out_7962091599851903711[77] = 0;
   out_7962091599851903711[78] = 0;
   out_7962091599851903711[79] = 0;
   out_7962091599851903711[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_3175192585520858313) {
   out_3175192585520858313[0] = state[0];
   out_3175192585520858313[1] = state[1];
   out_3175192585520858313[2] = state[2];
   out_3175192585520858313[3] = state[3];
   out_3175192585520858313[4] = state[4];
   out_3175192585520858313[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_3175192585520858313[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_3175192585520858313[7] = state[7];
   out_3175192585520858313[8] = state[8];
}
void F_fun(double *state, double dt, double *out_9023148077527241608) {
   out_9023148077527241608[0] = 1;
   out_9023148077527241608[1] = 0;
   out_9023148077527241608[2] = 0;
   out_9023148077527241608[3] = 0;
   out_9023148077527241608[4] = 0;
   out_9023148077527241608[5] = 0;
   out_9023148077527241608[6] = 0;
   out_9023148077527241608[7] = 0;
   out_9023148077527241608[8] = 0;
   out_9023148077527241608[9] = 0;
   out_9023148077527241608[10] = 1;
   out_9023148077527241608[11] = 0;
   out_9023148077527241608[12] = 0;
   out_9023148077527241608[13] = 0;
   out_9023148077527241608[14] = 0;
   out_9023148077527241608[15] = 0;
   out_9023148077527241608[16] = 0;
   out_9023148077527241608[17] = 0;
   out_9023148077527241608[18] = 0;
   out_9023148077527241608[19] = 0;
   out_9023148077527241608[20] = 1;
   out_9023148077527241608[21] = 0;
   out_9023148077527241608[22] = 0;
   out_9023148077527241608[23] = 0;
   out_9023148077527241608[24] = 0;
   out_9023148077527241608[25] = 0;
   out_9023148077527241608[26] = 0;
   out_9023148077527241608[27] = 0;
   out_9023148077527241608[28] = 0;
   out_9023148077527241608[29] = 0;
   out_9023148077527241608[30] = 1;
   out_9023148077527241608[31] = 0;
   out_9023148077527241608[32] = 0;
   out_9023148077527241608[33] = 0;
   out_9023148077527241608[34] = 0;
   out_9023148077527241608[35] = 0;
   out_9023148077527241608[36] = 0;
   out_9023148077527241608[37] = 0;
   out_9023148077527241608[38] = 0;
   out_9023148077527241608[39] = 0;
   out_9023148077527241608[40] = 1;
   out_9023148077527241608[41] = 0;
   out_9023148077527241608[42] = 0;
   out_9023148077527241608[43] = 0;
   out_9023148077527241608[44] = 0;
   out_9023148077527241608[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_9023148077527241608[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_9023148077527241608[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_9023148077527241608[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_9023148077527241608[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_9023148077527241608[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_9023148077527241608[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_9023148077527241608[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_9023148077527241608[53] = -9.8000000000000007*dt;
   out_9023148077527241608[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_9023148077527241608[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_9023148077527241608[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_9023148077527241608[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_9023148077527241608[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_9023148077527241608[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_9023148077527241608[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_9023148077527241608[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_9023148077527241608[62] = 0;
   out_9023148077527241608[63] = 0;
   out_9023148077527241608[64] = 0;
   out_9023148077527241608[65] = 0;
   out_9023148077527241608[66] = 0;
   out_9023148077527241608[67] = 0;
   out_9023148077527241608[68] = 0;
   out_9023148077527241608[69] = 0;
   out_9023148077527241608[70] = 1;
   out_9023148077527241608[71] = 0;
   out_9023148077527241608[72] = 0;
   out_9023148077527241608[73] = 0;
   out_9023148077527241608[74] = 0;
   out_9023148077527241608[75] = 0;
   out_9023148077527241608[76] = 0;
   out_9023148077527241608[77] = 0;
   out_9023148077527241608[78] = 0;
   out_9023148077527241608[79] = 0;
   out_9023148077527241608[80] = 1;
}
void h_25(double *state, double *unused, double *out_2256446242627236920) {
   out_2256446242627236920[0] = state[6];
}
void H_25(double *state, double *unused, double *out_7575440248946220157) {
   out_7575440248946220157[0] = 0;
   out_7575440248946220157[1] = 0;
   out_7575440248946220157[2] = 0;
   out_7575440248946220157[3] = 0;
   out_7575440248946220157[4] = 0;
   out_7575440248946220157[5] = 0;
   out_7575440248946220157[6] = 1;
   out_7575440248946220157[7] = 0;
   out_7575440248946220157[8] = 0;
}
void h_24(double *state, double *unused, double *out_986700759327517878) {
   out_986700759327517878[0] = state[4];
   out_986700759327517878[1] = state[5];
}
void H_24(double *state, double *unused, double *out_6124794497796822425) {
   out_6124794497796822425[0] = 0;
   out_6124794497796822425[1] = 0;
   out_6124794497796822425[2] = 0;
   out_6124794497796822425[3] = 0;
   out_6124794497796822425[4] = 1;
   out_6124794497796822425[5] = 0;
   out_6124794497796822425[6] = 0;
   out_6124794497796822425[7] = 0;
   out_6124794497796822425[8] = 0;
   out_6124794497796822425[9] = 0;
   out_6124794497796822425[10] = 0;
   out_6124794497796822425[11] = 0;
   out_6124794497796822425[12] = 0;
   out_6124794497796822425[13] = 0;
   out_6124794497796822425[14] = 1;
   out_6124794497796822425[15] = 0;
   out_6124794497796822425[16] = 0;
   out_6124794497796822425[17] = 0;
}
void h_30(double *state, double *unused, double *out_2413632910978366065) {
   out_2413632910978366065[0] = state[4];
}
void H_30(double *state, double *unused, double *out_658749907454603402) {
   out_658749907454603402[0] = 0;
   out_658749907454603402[1] = 0;
   out_658749907454603402[2] = 0;
   out_658749907454603402[3] = 0;
   out_658749907454603402[4] = 1;
   out_658749907454603402[5] = 0;
   out_658749907454603402[6] = 0;
   out_658749907454603402[7] = 0;
   out_658749907454603402[8] = 0;
}
void h_26(double *state, double *unused, double *out_1037742327441068596) {
   out_1037742327441068596[0] = state[7];
}
void H_26(double *state, double *unused, double *out_4270914279185419556) {
   out_4270914279185419556[0] = 0;
   out_4270914279185419556[1] = 0;
   out_4270914279185419556[2] = 0;
   out_4270914279185419556[3] = 0;
   out_4270914279185419556[4] = 0;
   out_4270914279185419556[5] = 0;
   out_4270914279185419556[6] = 0;
   out_4270914279185419556[7] = 1;
   out_4270914279185419556[8] = 0;
}
void h_27(double *state, double *unused, double *out_2745227025033608725) {
   out_2745227025033608725[0] = state[3];
}
void H_27(double *state, double *unused, double *out_2833513219255028313) {
   out_2833513219255028313[0] = 0;
   out_2833513219255028313[1] = 0;
   out_2833513219255028313[2] = 0;
   out_2833513219255028313[3] = 1;
   out_2833513219255028313[4] = 0;
   out_2833513219255028313[5] = 0;
   out_2833513219255028313[6] = 0;
   out_2833513219255028313[7] = 0;
   out_2833513219255028313[8] = 0;
}
void h_29(double *state, double *unused, double *out_5018038231787063716) {
   out_5018038231787063716[0] = state[1];
}
void H_29(double *state, double *unused, double *out_148518563140211218) {
   out_148518563140211218[0] = 0;
   out_148518563140211218[1] = 1;
   out_148518563140211218[2] = 0;
   out_148518563140211218[3] = 0;
   out_148518563140211218[4] = 0;
   out_148518563140211218[5] = 0;
   out_148518563140211218[6] = 0;
   out_148518563140211218[7] = 0;
   out_148518563140211218[8] = 0;
}
void h_28(double *state, double *unused, double *out_67325078236244658) {
   out_67325078236244658[0] = state[0];
}
void H_28(double *state, double *unused, double *out_5230917580209741792) {
   out_5230917580209741792[0] = 1;
   out_5230917580209741792[1] = 0;
   out_5230917580209741792[2] = 0;
   out_5230917580209741792[3] = 0;
   out_5230917580209741792[4] = 0;
   out_5230917580209741792[5] = 0;
   out_5230917580209741792[6] = 0;
   out_5230917580209741792[7] = 0;
   out_5230917580209741792[8] = 0;
}
void h_31(double *state, double *unused, double *out_5893895770881259270) {
   out_5893895770881259270[0] = state[8];
}
void H_31(double *state, double *unused, double *out_498764998434402904) {
   out_498764998434402904[0] = 0;
   out_498764998434402904[1] = 0;
   out_498764998434402904[2] = 0;
   out_498764998434402904[3] = 0;
   out_498764998434402904[4] = 0;
   out_498764998434402904[5] = 0;
   out_498764998434402904[6] = 0;
   out_498764998434402904[7] = 0;
   out_498764998434402904[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_3484915228138340995) {
  err_fun(nom_x, delta_x, out_3484915228138340995);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_101761715903343570) {
  inv_err_fun(nom_x, true_x, out_101761715903343570);
}
void car_H_mod_fun(double *state, double *out_7962091599851903711) {
  H_mod_fun(state, out_7962091599851903711);
}
void car_f_fun(double *state, double dt, double *out_3175192585520858313) {
  f_fun(state,  dt, out_3175192585520858313);
}
void car_F_fun(double *state, double dt, double *out_9023148077527241608) {
  F_fun(state,  dt, out_9023148077527241608);
}
void car_h_25(double *state, double *unused, double *out_2256446242627236920) {
  h_25(state, unused, out_2256446242627236920);
}
void car_H_25(double *state, double *unused, double *out_7575440248946220157) {
  H_25(state, unused, out_7575440248946220157);
}
void car_h_24(double *state, double *unused, double *out_986700759327517878) {
  h_24(state, unused, out_986700759327517878);
}
void car_H_24(double *state, double *unused, double *out_6124794497796822425) {
  H_24(state, unused, out_6124794497796822425);
}
void car_h_30(double *state, double *unused, double *out_2413632910978366065) {
  h_30(state, unused, out_2413632910978366065);
}
void car_H_30(double *state, double *unused, double *out_658749907454603402) {
  H_30(state, unused, out_658749907454603402);
}
void car_h_26(double *state, double *unused, double *out_1037742327441068596) {
  h_26(state, unused, out_1037742327441068596);
}
void car_H_26(double *state, double *unused, double *out_4270914279185419556) {
  H_26(state, unused, out_4270914279185419556);
}
void car_h_27(double *state, double *unused, double *out_2745227025033608725) {
  h_27(state, unused, out_2745227025033608725);
}
void car_H_27(double *state, double *unused, double *out_2833513219255028313) {
  H_27(state, unused, out_2833513219255028313);
}
void car_h_29(double *state, double *unused, double *out_5018038231787063716) {
  h_29(state, unused, out_5018038231787063716);
}
void car_H_29(double *state, double *unused, double *out_148518563140211218) {
  H_29(state, unused, out_148518563140211218);
}
void car_h_28(double *state, double *unused, double *out_67325078236244658) {
  h_28(state, unused, out_67325078236244658);
}
void car_H_28(double *state, double *unused, double *out_5230917580209741792) {
  H_28(state, unused, out_5230917580209741792);
}
void car_h_31(double *state, double *unused, double *out_5893895770881259270) {
  h_31(state, unused, out_5893895770881259270);
}
void car_H_31(double *state, double *unused, double *out_498764998434402904) {
  H_31(state, unused, out_498764998434402904);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_init(car);
