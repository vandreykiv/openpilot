#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_3484915228138340995);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_101761715903343570);
void car_H_mod_fun(double *state, double *out_7962091599851903711);
void car_f_fun(double *state, double dt, double *out_3175192585520858313);
void car_F_fun(double *state, double dt, double *out_9023148077527241608);
void car_h_25(double *state, double *unused, double *out_2256446242627236920);
void car_H_25(double *state, double *unused, double *out_7575440248946220157);
void car_h_24(double *state, double *unused, double *out_986700759327517878);
void car_H_24(double *state, double *unused, double *out_6124794497796822425);
void car_h_30(double *state, double *unused, double *out_2413632910978366065);
void car_H_30(double *state, double *unused, double *out_658749907454603402);
void car_h_26(double *state, double *unused, double *out_1037742327441068596);
void car_H_26(double *state, double *unused, double *out_4270914279185419556);
void car_h_27(double *state, double *unused, double *out_2745227025033608725);
void car_H_27(double *state, double *unused, double *out_2833513219255028313);
void car_h_29(double *state, double *unused, double *out_5018038231787063716);
void car_H_29(double *state, double *unused, double *out_148518563140211218);
void car_h_28(double *state, double *unused, double *out_67325078236244658);
void car_H_28(double *state, double *unused, double *out_5230917580209741792);
void car_h_31(double *state, double *unused, double *out_5893895770881259270);
void car_H_31(double *state, double *unused, double *out_498764998434402904);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}