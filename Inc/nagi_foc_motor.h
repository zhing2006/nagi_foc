#ifndef __NAGI_FOC_MOTOR_H__
#define __NAGI_FOC_MOTOR_H__

#include <stddef.h>
#include <stdbool.h>
#include <stdint.h>

#include <math.h>
#include <arm_math.h>

/// @brief FOC motor error codes.
typedef enum {
  NAGI_FOC_MOTOR_OK = 0, ///< No error.
  NAGI_FOC_MOTOR_ERROR,  ///< General error.
  MAGI_FOC_MOTOR_HANLDLE_NULL, ///< Handle is NULL.
  NAGI_FOC_MOTOR_POINTER_NULL, ///< Pointer is NULL.
  NAGI_FOC_MOTOR_INVALID_ARGUMENT, ///< Invalid argument.
} nagi_foc_error_t;

/// @brief Set PWM duty cycle function typedef.
typedef void (*nagi_foc_set_pwm_duty_fn_t)(float, float, float);
/// @brief Delay function typedef.
typedef void (*nagi_foc_delay_fn_t)(uint32_t);

/// @brief FOC motor structure configuration.
typedef struct nagi_foc_motor_config_t {
  /// @brief Set PWM duty cycle function pointer.
  nagi_foc_set_pwm_duty_fn_t set_pwm_duty_fn;
  /// @brief Delay function pointer.
  nagi_foc_delay_fn_t delay_fn;
  /// @brief Pole pairs.
  uint32_t pole_pairs;
  /// @brief Voltage limit.
  float voltage_limit;
  /// @brief Position cycle.
  float position_cycle;

  /// @brief P gain for position control.
  float kp_position;
  /// @brief I gain for position control.
  float ki_position;
  /// @brief D gain for position control.
  float kd_position;
  /// @brief P gain for speed control.
  float kp_speed;
  /// @brief I gain for speed control.
  float ki_speed;
  /// @brief D gain for speed control.
  float kd_speed;
  /// @brief P gain for d axis current control.
  float kp_current_d;
  /// @brief I gain for d axis current control.
  float ki_current_d;
  /// @brief D gain for d axis current control.
  float kd_current_d;
  /// @brief P gain for q axis current control.
  float kp_current_q;
  /// @brief I gain for q axis current control.
  float ki_current_q;
  /// @brief D gain for q axis current control.
  float kd_current_q;
} nagi_foc_motor_config_t;

/// @brief FOC motor structure.
typedef struct nagi_foc_motor_t {
  /// @brief Set PWM duty cycle function pointer.
  nagi_foc_set_pwm_duty_fn_t set_pwm_duty_fn;
  /// @brief Delay function pointer.
  nagi_foc_delay_fn_t delay_fn;
  /// @brief Pole pairs.
  uint32_t pole_pairs;
  /// @brief Voltage limit.
  float voltage_limit;
  /// @brief Position cycle.
  float position_cycle;

  /// @brief Is motor calibrated.
  bool is_calibrated;
  /// @brief Is logical angle ready.
  bool is_logical_angle_ready;
  /// @brief Is speed ready.
  bool is_speed_ready;

  /// @brief Current in d axis.
  float i_d;
  /// @brief Current in q axis.
  float i_q;
  /// @brief Speed.
  float speed;
  /// @brief Encoder angle.
  float encoder_angle;
  /// @brief Zero angle.
  float zero_angle;
  /// @brief Logical angle.
  float logical_angle;

  /// @brief PID controller for position control.
  arm_pid_instance_f32 pid_position;
  /// @brief PID controller for speed control.
  arm_pid_instance_f32 pid_speed;
  /// @brief PID controller for d axis current control.
  arm_pid_instance_f32 pid_current_d;
  /// @brief PID controller for q axis current control.
  arm_pid_instance_f32 pid_current_q;
} nagi_foc_motor_t;

/// @brief Normalize an angular difference to be within the range [-full_cycle/2, full_cycle/2].
/// @param angle_diff The input angular difference (e.g., the difference between two angles).
/// @param full_cycle The full cycle of the angle (e.g., 2*PI in radians or 360.0 in degrees).
/// @return The normalized angular difference.
float nagi_foc_motor_normalize_angle_diff(float angle_diff, float full_cycle);

/// @brief Normalize an angle to be within the range [0, 2*PI].
/// @param angle_diff The input angular difference (e.g., the difference between two angles).
/// @return The normalized angular difference.
static inline float nagi_foc_motor_normalize_angle_2pi(float angle_diff) {
  if (angle_diff > M_PI) {
    angle_diff -= (2.0f * M_PI);
  } else if (angle_diff < -M_PI) {
    angle_diff += (2.0f * M_PI);
  }

  return angle_diff;
}

/// @brief Initialize FOC motor.
/// @param[in] motor FOC motor structure.
/// @param[in] config FOC motor configuration.
/// @return FOC motor error code.
nagi_foc_error_t nagi_foc_motor_init(nagi_foc_motor_t *pmotor, const nagi_foc_motor_config_t *pconfig);

/// @brief Set the motor's torque.
/// @param[in] motor FOC motor structure.
/// @param[in] electrical_angle Electrical angle of the motor.
/// @param[in] torque_d Torque in d axis.
/// @param[in] torque_q Torque in q axis.
/// @return FOC motor error code.
nagi_foc_error_t nagi_foc_motor_set_torque(nagi_foc_motor_t *pmotor, float electrical_angle, float torque_d, float torque_q);

/// @brief Calibrate the motor's zero angle.
/// @param[in] motor FOC motor structure.
/// @return FOC motor error code.
nagi_foc_error_t nagi_foc_motor_calibrate(nagi_foc_motor_t *pmotor);

/// @brief Calculate the current in d and q axis.
/// @param[in] motor FOC motor structure.
/// @param[in] i_a Current in A phase.
/// @param[in] i_b Current in B phase.
/// @param[out] pi_d Current in d axis.
/// @param[out] pi_q Current in q axis.
nagi_foc_error_t nagi_foc_motor_calculate_current(nagi_foc_motor_t *pmotor, float i_a, float i_b, float *pi_d, float *pi_q);

/// @brief Position control.
/// @param[in] motor FOC motor structure.
/// @param[in] target_angle Target angle.
/// @return FOC motor error code.
nagi_foc_error_t nagi_foc_motor_position_control(nagi_foc_motor_t *pmotor, float target_angle);

/// @brief Speed control.
/// @param[in] motor FOC motor structure.
/// @param[in] target_speed Target speed.
/// @return FOC motor error code.
nagi_foc_error_t nagi_foc_motor_speed_control(nagi_foc_motor_t *pmotor, float target_speed);

/// @brief Motor current control.
/// @param[in] motor FOC motor structure.
/// @param[in] target_i_d Target current in d axis.
/// @param[in] target_i_q Target current in q axis.
/// @return FOC motor error code.
nagi_foc_error_t nagi_foc_motor_current_control(nagi_foc_motor_t *pmotor, float target_i_d, float target_i_q);

/// @brief Motor speed with maximum current control.
/// @param[in] motor FOC motor structure.
/// @param[in] target_speed Target speed.
/// @param[in] max_current Maximum current.
/// @return FOC motor error code.
nagi_foc_error_t nagi_foc_motor_speed_current_control(nagi_foc_motor_t *pmotor, float target_speed, float max_current);

/// @brief Motor position, speed, and current control.
/// @param[in] motor FOC motor structure.
/// @param[in] target_angle Target angle.
/// @param[in] max_speed Maximum speed.
/// @param[in] max_current Maximum current.
/// @return FOC motor error code.
nagi_foc_error_t nagi_foc_motor_position_speed_current_control(nagi_foc_motor_t *pmotor, float target_angle, float max_speed, float max_current);

/// @brief Get the mechanical angle of the motor.
/// @param[in] motor FOC motor structure.
/// @return Mechanical angle of the motor.
inline static float nagi_foc_motor_get_mech_angle(nagi_foc_motor_t *pmotor) {
  return pmotor->encoder_angle - pmotor->zero_angle;
}

/// @brief Get the electrical angle of the motor.
/// @param[in] motor FOC motor structure.
/// @return Electrical angle of the motor.
inline static float nagi_foc_motor_get_elec_angle(nagi_foc_motor_t *pmotor) {
  return nagi_foc_motor_get_mech_angle(pmotor) * pmotor->pole_pairs;
}

#endif // __NAGI_FOC_MOTOR_H__