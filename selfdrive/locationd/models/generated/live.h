#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_4341031410396474634);
void live_err_fun(double *nom_x, double *delta_x, double *out_4797538582936705820);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_8826134262153630586);
void live_H_mod_fun(double *state, double *out_6756890915510285479);
void live_f_fun(double *state, double dt, double *out_1045568842516687553);
void live_F_fun(double *state, double dt, double *out_1217235508240824483);
void live_h_4(double *state, double *unused, double *out_1187195172086222192);
void live_H_4(double *state, double *unused, double *out_2881209398931606613);
void live_h_9(double *state, double *unused, double *out_468430767723912667);
void live_H_9(double *state, double *unused, double *out_8278315739513497533);
void live_h_10(double *state, double *unused, double *out_7900897131025096528);
void live_H_10(double *state, double *unused, double *out_9004649023905477017);
void live_h_12(double *state, double *unused, double *out_8041879091755035437);
void live_H_12(double *state, double *unused, double *out_7900665806963568408);
void live_h_31(double *state, double *unused, double *out_7319727340276679641);
void live_H_31(double *state, double *unused, double *out_6247871456304213989);
void live_h_32(double *state, double *unused, double *out_4286017220412784756);
void live_H_32(double *state, double *unused, double *out_7759747416440127621);
void live_h_13(double *state, double *unused, double *out_7905065231247966384);
void live_H_13(double *state, double *unused, double *out_3633141697491679160);
void live_h_14(double *state, double *unused, double *out_468430767723912667);
void live_H_14(double *state, double *unused, double *out_8278315739513497533);
void live_h_33(double *state, double *unused, double *out_524656203599624517);
void live_H_33(double *state, double *unused, double *out_9048315612766480023);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}