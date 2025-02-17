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

/// @brief Calculate the voltage space vector pulse width modulation.
/// @param phi The angle of the rotor.
/// @param d The D axis voltage.
/// @param q The Q axis voltage.
/// @param pu The U phase voltage.
/// @param pv The V phase voltage.
/// @param pw The W phase voltage.
__attribute__((weak)) void nagi_foc_motor_calc_svpwm(float phi, float d, float q, float *pu, float *pv, float *pw) {
  // Calculate sin and cos values of the electrical angle.
  float sin_phi, cos_phi;
  arm_sin_cos_f32(phi * RAD_TO_DEG, &sin_phi, &cos_phi);

  // Calculate the alpha and beta axis.
  float i_alpha, i_beta;
  arm_inv_park_f32(d, q, &i_alpha, &i_beta, sin_phi, cos_phi);

  // Calculate the a, b, c axis.
  float i_a, i_b, i_c;
  arm_inv_clarke_f32(i_alpha, i_beta, &i_a, &i_b);
  i_c = -i_a - i_b;

  // Calculate the maximum and minimum.
  float max = MAX(MAX(i_a, i_b), i_c);
  float min = MIN(MIN(i_a, i_b), i_c);

  // Calculate the offset.
  float offset = -(max + min) * 0.5f;

  // Return the calculated duty cycles.
  *pu = 0.5f + 0.5f * (i_a + offset);
  *pv = 0.5f + 0.5f * (i_b + offset);
  *pw = 0.5f + 0.5f * (i_c + offset);
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
    return MAGI_FOC_MOTOR_HANDLE_NULL;
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
  pmotor->is_calibrated = pmotor->is_logical_angle_ready = pmotor->is_speed_ready = pmotor->is_current_ready = false;
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
    return MAGI_FOC_MOTOR_HANDLE_NULL;
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
    return MAGI_FOC_MOTOR_HANDLE_NULL;
  }

  if (pmotor->delay_fn == NULL) {
    return NAGI_FOC_MOTOR_POINTER_NULL;
  }

  nagi_foc_error_t err = nagi_foc_motor_set_torque(pmotor, 0.0f, 1.0f, 0.0f);
  if (err != NAGI_FOC_MOTOR_OK) {
    return err;
  }
  pmotor->delay_fn(1000);
  pmotor->zero_angle = pmotor->encoder_angle;
  err = nagi_foc_motor_set_torque(pmotor, 0.0f, 0.0f, 0.0f);
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
    return MAGI_FOC_MOTOR_HANDLE_NULL;
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
    return MAGI_FOC_MOTOR_HANDLE_NULL;
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
    return MAGI_FOC_MOTOR_HANDLE_NULL;
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
    return MAGI_FOC_MOTOR_HANDLE_NULL;
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
    return MAGI_FOC_MOTOR_HANDLE_NULL;
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
    return MAGI_FOC_MOTOR_HANDLE_NULL;
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
