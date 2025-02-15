#include "nagi_foc_motor.h"

#include <stdbool.h>

#ifndef DEG_TO_RAD
#define DEG_TO_RAD (0.017453292519943295)
#endif

#ifndef RAD_TO_DEG
#define RAD_TO_DEG (57.29577951308232)
#endif

#ifndef MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif

#ifndef MAX
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#endif

#ifndef CLAMP
#define CLAMP(x, min, max) (MIN(MAX(x, min), max))
#endif

#if !defined(UNUSED)
#define UNUSED(X) (void)X      /* To avoid gcc/g++ warnings */
#endif

#ifndef M_3PI_2
#define M_3PI_2 (4.7123889803846898576939650749193)
#endif

#ifndef M_SQRT3
#define M_SQRT3 (1.7320508075688772935274463415059)
#endif

static float _normalize_angle_0_2pi(float angle) {
  // Normalize the angle to be within the range [0, 2*PI].
  angle = fmodf(angle, 2.0f * M_PI);
  if (angle < 0.0f) {
    angle += 2.0f * M_PI;
  }
  return angle;
}

/// @brief Calculate the voltage space vector pulse width modulation.
/// @param phi The angle of the rotor.
/// @param d The D axis voltage.
/// @param q The Q axis voltage.
/// @param pu The U phase voltage.
/// @param pv The V phase voltage.
/// @param pw The W phase voltage.
__attribute__((weak)) void nagi_foc_motor_calc_svpwm(float phi, float d, float q, float *pu, float *pv, float *pw) {
  // Normalize the angle phi to [0, 2pi) and compute sine and cosine in degrees.
  phi = _normalize_angle_0_2pi(phi);
  float sin_phi, cos_phi;
  arm_sin_cos_f32(phi * RAD_TO_DEG, &sin_phi, &cos_phi);

  // Inverse Park transform: transform (d,q) to alpha-beta coordinates.
  float u_alpha, u_beta;
  arm_inv_park_f32(d, q, &u_alpha, &u_beta, sin_phi, cos_phi);

  // Calculate the reference voltage magnitude and its angle (in radians).
  float Vref;
  arm_sqrt_f32(u_alpha * u_alpha + u_beta * u_beta, &Vref);
  float angle_ref;
  arm_atan2_f32(u_beta, u_alpha, &angle_ref);
  if (angle_ref < 0.0f) {
    angle_ref += 2 * M_PI;
  }

  // Determine the SVPWM sector (each sector spans pi/3 radians).
  // Use zero-based index for the lookup table.
  const int sector_index = (int)(angle_ref / (M_PI / 3.0f));
  const float angle_sector = angle_ref - sector_index * (M_PI / 3.0f);

  // Compute T1 and T2 durations for the two active vectors.
  // Here Ts is assumed to be 1.0 (i.e. T1+T2+T0 = 1).
  const float T1 = M_SQRT3 * Vref * arm_sin_f32((M_PI / 3.0f) - angle_sector);
  const float T2 = M_SQRT3 * Vref * arm_sin_f32(angle_sector);

  // Compute the zero vector time and ensure it is non-negative.
  float T0 = 1.0f - T1 - T2;
  if (T0 < 0.0f) {
    T0 = 0.0f; // simple saturation if Vref is too high.
  }
  const float t0_half = T0 * 0.5f;

  // Define a lookup table for the active vector contributions.
  // For each sector (0 to 5) the three phases (A, B, C) receive a weighted sum of T1 and T2:
  // Each LUT entry is a pair: {coefficient_for_T1, coefficient_for_T2}.
  // The computed phase duty, after subtracting the t0 baseline, is:
  // Phase_X = (coef_T1 * T1 + coef_T2 * T2) + t0_half.
  //
  // LUT breakdown (zero-based index corresponds to sector):
  // Sector 0 (SVPWM sector 1): A = T1+T2, B = T2, C = 0
  // Sector 1 (SVPWM sector 2): A = T1, B = T1+T2, C = 0
  // Sector 2 (SVPWM sector 3): A = 0, B = T1+T2, C = T2
  // Sector 3 (SVPWM sector 4): A = 0, B = T1, C = T1+T2
  // Sector 4 (SVPWM sector 5): A = T2, B = 0, C = T1+T2
  // Sector 5 (SVPWM sector 6): A = T1+T2, B = 0, C = T1
  static const float svpwm_lut[6][3][2] = {
    { {1.0f, 1.0f}, {0.0f, 1.0f}, {0.0f, 0.0f} }, // Sector 0: Phase A, B, C
    { {1.0f, 0.0f}, {1.0f, 1.0f}, {0.0f, 0.0f} }, // Sector 1
    { {0.0f, 0.0f}, {1.0f, 1.0f}, {0.0f, 1.0f} }, // Sector 2
    { {0.0f, 0.0f}, {1.0f, 0.0f}, {1.0f, 1.0f} }, // Sector 3
    { {0.0f, 1.0f}, {0.0f, 0.0f}, {1.0f, 1.0f} }, // Sector 4
    { {1.0f, 1.0f}, {0.0f, 0.0f}, {1.0f, 0.0f} }  // Sector 5
  };

  // Compute the phase duty cycles using the lookup table.
  const float t_a = svpwm_lut[sector_index][0][0] * T1 + svpwm_lut[sector_index][0][1] * T2 + t0_half;
  const float t_b = svpwm_lut[sector_index][1][0] * T1 + svpwm_lut[sector_index][1][1] * T2 + t0_half;
  const float t_c = svpwm_lut[sector_index][2][0] * T1 + svpwm_lut[sector_index][2][1] * T2 + t0_half;

  // Return the calculated duty cycles (normalized to 0.0 ~ 1.0) via the pointers.
  *pu = t_a;
  *pv = t_b;
  *pw = t_c;
}

float nagi_foc_motor_normalize_angle_diff(float angle_diff, float full_cycle) {
  const float half_cycle = full_cycle * 0.5f;
  if (angle_diff > half_cycle) {
    angle_diff -= full_cycle;
  } else if (angle_diff < -half_cycle) {
    angle_diff += full_cycle;
  }

  return angle_diff;
}

nagi_foc_error_t nagi_foc_motor_init(nagi_foc_motor_t *pmotor, const nagi_foc_motor_config_t *pconfig) {
  if (pmotor == NULL) {
    return MAGI_FOC_MOTOR_HANLDLE_NULL;
  }

  if (pconfig == NULL) {
    return NAGI_FOC_MOTOR_ERROR;
  }

  if (pconfig->set_pwm_duty_fn == NULL || pconfig->delay_fn == NULL)
  {
    return NAGI_FOC_MOTOR_POINTER_NULL;
  }

  if (pconfig->pole_pairs == 0 || pconfig->voltage_limit <= FLT_EPSILON || pconfig->position_cycle <= FLT_EPSILON) {
    return NAGI_FOC_MOTOR_INVALID_ARGUMENT;
  }

  pmotor->set_pwm_duty_fn = pconfig->set_pwm_duty_fn;
  pmotor->delay_fn = pconfig->delay_fn;
  pmotor->pole_pairs = pconfig->pole_pairs;
  pmotor->voltage_limit = pconfig->voltage_limit;
  pmotor->position_cycle = pconfig->position_cycle;
  pmotor->is_calibrated = pmotor->is_logical_angle_ready = pmotor->is_speed_ready = false;
  pmotor->i_d = pmotor->i_q = 0.0f;
  pmotor->speed = 0.0f;
  pmotor->encoder_angle = 0.0f;
  pmotor->zero_angle = 0.0f;
  pmotor->logical_angle = 0.0f;

  pmotor->pid_position.Kp = pconfig->kp_position;
  pmotor->pid_position.Ki = pconfig->ki_position;
  pmotor->pid_position.Kd = pconfig->kd_position;
  arm_pid_init_f32(&pmotor->pid_position, 1);
  pmotor->pid_speed.Kp = pconfig->kp_speed;
  pmotor->pid_speed.Ki = pconfig->ki_speed;
  pmotor->pid_speed.Kd = pconfig->kd_speed;
  arm_pid_init_f32(&pmotor->pid_speed, 1);
  pmotor->pid_current_d.Kp = pconfig->kp_current_d;
  pmotor->pid_current_d.Ki = pconfig->ki_current_d;
  pmotor->pid_current_d.Kd = pconfig->kd_current_d;
  arm_pid_init_f32(&pmotor->pid_current_d, 1);
  pmotor->pid_current_q.Kp = pconfig->kp_current_q;
  pmotor->pid_current_q.Ki = pconfig->ki_current_q;
  pmotor->pid_current_q.Kd = pconfig->kd_current_q;
  arm_pid_init_f32(&pmotor->pid_current_q, 1);

  return NAGI_FOC_MOTOR_OK;
}

nagi_foc_error_t nagi_foc_motor_set_torque(nagi_foc_motor_t *pmotor, float electrical_angle, float torque_d, float torque_q) {
  if (pmotor == NULL) {
    return MAGI_FOC_MOTOR_HANLDLE_NULL;
  }

  if (pmotor->set_pwm_duty_fn == NULL) {
    return NAGI_FOC_MOTOR_POINTER_NULL;
  }

  torque_d = CLAMP(torque_d, -1.0f, 1.0f);
  torque_q = CLAMP(torque_q, -1.0f, 1.0f);

  float u, v, w;
  nagi_foc_motor_calc_svpwm(electrical_angle, torque_d, torque_q, &u, &v, &w);
  pmotor->set_pwm_duty_fn(u, v, w);

  return NAGI_FOC_MOTOR_OK;
}

nagi_foc_error_t nagi_foc_motor_calibrate(nagi_foc_motor_t *pmotor) {
  if (pmotor == NULL) {
    return MAGI_FOC_MOTOR_HANLDLE_NULL;
  }

  if (pmotor->delay_fn == NULL) {
    return NAGI_FOC_MOTOR_POINTER_NULL;
  }

  nagi_foc_error_t err = nagi_foc_motor_set_torque(pmotor, M_3PI_2, 0.0f, 0.5f);
  if (err != NAGI_FOC_MOTOR_OK) {
    return err;
  }
  pmotor->delay_fn(1000);
  pmotor->zero_angle = pmotor->encoder_angle;
  err = nagi_foc_motor_set_torque(pmotor, M_3PI_2, 0.0f, 0.0f);
  if (err != NAGI_FOC_MOTOR_OK) {
    return err;
  }
  pmotor->delay_fn(100);
  pmotor->set_pwm_duty_fn(0.0f, 0.0f, 0.0f);
  pmotor->logical_angle = 0.0f;
  pmotor->is_calibrated = true;
  return err;
}

nagi_foc_error_t nagi_foc_motor_calculate_current(nagi_foc_motor_t *pmotor, float i_a, float i_b, float *pi_d, float *pi_q) {
  if (pmotor == NULL) {
    return MAGI_FOC_MOTOR_HANLDLE_NULL;
  }

  if (pi_d == NULL || pi_q == NULL) {
    return NAGI_FOC_MOTOR_POINTER_NULL;
  }

  // Clarke transform.
  float i_alpha, i_beta;
  arm_clarke_f32(i_a, i_b, &i_alpha, &i_beta);

  // Park transform.
  float sin_phi, cos_phi;
  arm_sin_cos_f32(nagi_foc_motor_get_elec_angle(pmotor) * RAD_TO_DEG, &sin_phi, &cos_phi);
  arm_park_f32(i_alpha, i_beta, pi_d, pi_q, sin_phi, cos_phi);

  return NAGI_FOC_MOTOR_OK;
}

static float _position_loop(nagi_foc_motor_t *pmotor, float target_angle) {
  const float angle_diff = target_angle - pmotor->logical_angle;
  const float torque_q = arm_pid_f32(&pmotor->pid_position, angle_diff);
  return torque_q;
}

nagi_foc_error_t nagi_foc_motor_position_control(nagi_foc_motor_t *pmotor, float target_angle) {
  if (pmotor == NULL) {
    return MAGI_FOC_MOTOR_HANLDLE_NULL;
  }

  const float q = _position_loop(pmotor, target_angle);
  return nagi_foc_motor_set_torque(pmotor, nagi_foc_motor_get_elec_angle(pmotor), 0.0f, q);
}

static float _speed_loop(nagi_foc_motor_t *pmotor, float target_speed) {
  const float speed_diff = target_speed - pmotor->speed;
  const float torque_q = arm_pid_f32(&pmotor->pid_speed, speed_diff);
  return torque_q;
}

nagi_foc_error_t nagi_foc_motor_speed_control(nagi_foc_motor_t *pmotor, float target_speed) {
  if (pmotor == NULL) {
    return MAGI_FOC_MOTOR_HANLDLE_NULL;
  }

  const float q = _speed_loop(pmotor, target_speed);
  return nagi_foc_motor_set_torque(pmotor, nagi_foc_motor_get_elec_angle(pmotor), 0.0f, q);
}

static float _current_d_loop(nagi_foc_motor_t *pmotor, float target_i_d) {
  const float output_limit_max = 1.0f;
  const float output_limit_min = -1.0f;
  // Anti-windup Gain.
  const float Kt = 0.7f;

  const float d_diff = target_i_d - pmotor->i_d;
  const float out_unlimited = arm_pid_f32(&pmotor->pid_current_d, d_diff);
  const float out_limited = CLAMP(out_unlimited, output_limit_min, output_limit_max);
  // Integral Anti-windup.
  const float error_integral_windup = out_limited - out_unlimited;
  pmotor->pid_current_d.state[0] -= (pmotor->pid_current_d.Ki * Kt * error_integral_windup);

  return out_limited;
}

static float _current_q_loop(nagi_foc_motor_t *pmotor, float target_i_q) {
  const float output_limit_max = 1.0f;
  const float output_limit_min = -1.0f;
  // Anti-windup Gain.
  const float Kt = 0.7f;

  const float q_diff = target_i_q - pmotor->i_q;
  const float out_unlimited = arm_pid_f32(&pmotor->pid_current_q, q_diff);
  const float out_limited = CLAMP(out_unlimited, output_limit_min, output_limit_max);
  // Integral Anti-windup.
  const float error_integral_windup = out_limited - out_unlimited;
  pmotor->pid_current_q.state[0] -= (pmotor->pid_current_q.Ki * Kt * error_integral_windup);

  return out_limited;
}

nagi_foc_error_t nagi_foc_motor_current_control(nagi_foc_motor_t *pmotor, float target_i_d, float target_i_q) {
  if (pmotor == NULL) {
    return MAGI_FOC_MOTOR_HANLDLE_NULL;
  }

  if (pmotor->set_pwm_duty_fn == NULL) {
    return NAGI_FOC_MOTOR_POINTER_NULL;
  }

  const float torque_d = _current_d_loop(pmotor, target_i_d);
  const float torque_q = _current_q_loop(pmotor, target_i_q);
  return nagi_foc_motor_set_torque(pmotor, nagi_foc_motor_get_elec_angle(pmotor), torque_d, torque_q);
}

nagi_foc_error_t nagi_foc_motor_speed_current_control(nagi_foc_motor_t *pmotor, float target_speed, float max_current) {
  if (pmotor == NULL) {
    return MAGI_FOC_MOTOR_HANLDLE_NULL;
  }

  if (pmotor->set_pwm_duty_fn == NULL) {
    return NAGI_FOC_MOTOR_POINTER_NULL;
  }

  float target_i_q = _speed_loop(pmotor, target_speed);
  if (target_i_q > max_current) {
    target_i_q = max_current;
  } else if (target_i_q < -max_current) {
    target_i_q = -max_current;
  }

  return nagi_foc_motor_current_control(pmotor, 0.0f, target_i_q);
}

nagi_foc_error_t nagi_foc_motor_position_speed_current_control(nagi_foc_motor_t *pmotor, float target_angle, float max_speed, float max_current) {
  if (pmotor == NULL) {
    return MAGI_FOC_MOTOR_HANLDLE_NULL;
  }

  if (pmotor->set_pwm_duty_fn == NULL) {
    return NAGI_FOC_MOTOR_POINTER_NULL;
  }

  float target_speed = _position_loop(pmotor, target_angle);
  if (target_speed > max_speed) {
    target_speed = max_speed;
  } else if (target_speed < -max_speed) {
    target_speed = -max_speed;
  }

  return nagi_foc_motor_speed_current_control(pmotor, target_speed, max_current);
}