# Code for handling the kinematics of tetra robots
#
# Copyright (C) 2016,2017  Kevin O'Connor <kevin@koconnor.net>
# Tetra version added 2018 Jonas Forssell <jonasforssell@yahoo.se>
#
# This file may be distributed under the terms of the GNU GPLv3 license.
#
# The stepper position is the same as the length of the arm. This means that
# when homing, the stepper position becomes the defined arm length
#
import math, logging
import stepper, homing

StepList = (0, 1, 2)

# Slow moves once the ratio of anchor to XY movement exceeds SLOW_RATIO
SLOW_RATIO = 3.

class TetraKinematics:
    def __init__(self, toolhead, printer, config):
        stepper_configs = [config.getsection('stepper_' + n)
                           for n in ['a', 'b', 'c']]
        
        # Set default endstop position, equal to the distance of the nozzle from the bed
        # When the nozzle is in centre (0,0) position.
        stepper_a = stepper.PrinterHomingStepper(printer, stepper_configs[0])
        stepper_b = stepper.PrinterHomingStepper(
            printer, stepper_configs[1],
            default_position=stepper_a.position_endstop)
        stepper_c = stepper.PrinterHomingStepper(
            printer, stepper_configs[2],
            default_position=stepper_a.position_endstop)
        
        # Radius of the circle where the anchors are placed
        self.steppers = [stepper_a, stepper_b, stepper_c]
        self.need_motor_enable = self.need_home = True
        self.radius = radius = config.getfloat('tetra_radius', above=0.)
        
        # Arm length (length of wire when in endstop position
        arm_length_a = stepper_configs[0].getfloat('arm_length', above=radius)
        self.arm_lengths = arm_lengths = [
            sconfig.getfloat('arm_length', arm_length_a, above=radius)
            for sconfig in stepper_configs]
        self.arm2 = [arm**2 for arm in arm_lengths]
       
        # Set endstop position of steppers in local coordinates
        self.endstops = arm_lengths
        
        # Set up building area limits (in cartesian coordinates)
        self.limit_xy2 = -1.
        self.max_z = min([s.position_endstop for s in self.steppers])
        self.min_z = config.getfloat('minimum_z_position', 0, maxval=self.max_z)
        
        # The z-coordinate when one of the arms are straight down.
        # From this point the build cylinder must become a cone
        self.limit_z = min([s.position_endstop - (arm - sqrt(arm**2 - radius**2))
                            for s, arm in zip(self.steppers, arm_lengths)])
                         
        # Logging into about printer setup
        logging.info(
            "Tetra max build height %.2fmm (radius tapered above %.2fmm)" % (
                self.max_z, self.limit_z))
                         
        # Setup stepper velocities
        self.max_velocity, self.max_accel = toolhead.get_max_velocity()
        self.max_z_velocity = config.getfloat(
            'max_z_velocity', self.max_velocity,
            above=0., maxval=self.max_velocity)
        max_halt_velocity = toolhead.get_max_axis_halt()
        for s in self.steppers:
            s.set_max_jerk(max_halt_velocity, self.max_accel)
                         
        # Determine anchor xyz locations in cartesian space
        self.angles = [sconfig.getfloat('angle', angle)
                       for sconfig, angle in zip(stepper_configs,
                                                 [210., 330., 90.])]
        
        self.anchors = [(math.cos(math.radians(angle)) * radius,
                        math.sin(math.radians(angle)) * radius,
                        es.position_endstop + sqrt(arm**2 - radius**2))
                       for angle, es, arm in zip(self.angles, self.steppers, arm_lengths)]
                         
        # Find the point where an XY move could result in excessive
        # stepper movement
        half_min_step_dist = min([s.step_dist for s in self.steppers]) * .5
        min_arm_length = min(arm_lengths)
        
        def ratio_to_dist(ratio):
            return (ratio * math.sqrt(min_arm_length**2 / (ratio**2 + 1.)
                    - half_min_step_dist**2) + half_min_step_dist)
        
        self.slow_xy2 = (ratio_to_dist(SLOW_RATIO) - radius)**2
        self.very_slow_xy2 = (ratio_to_dist(2. * SLOW_RATIO) - radius)**2
        self.max_xy2 = min(radius, min_arm_length - radius,
                           ratio_to_dist(4. * SLOW_RATIO) - radius)**2
                         
        # Info regarding printer speed configuration                 
        logging.info(
            "Tetra max build radius %.2fmm (moves slowed past %.2fmm and %.2fmm)"
            % (math.sqrt(self.max_xy2), math.sqrt(self.slow_xy2),
               math.sqrt(self.very_slow_xy2)))
        
        # Set start position
        self.set_position([0., 0., 0.], ())
                         
    def get_steppers(self, flags=""):
        return list(self.steppers)
                         
    # return length of the three arms as the actuator position by calculating the hypotenusa
    # The arm length is directly related to the stepper motor position so this
    # is the local coordinate position for the stepper
    def _cartesian_to_actuator(self, coord):
        return [math.sqrt((self.anchors[i][0] - coord[0])**2
                        + (self.anchors[i][1] - coord[1])**2 
                        + (self.anchors[i][2] - coord[2])**2)
                for i in StepList]

    # Derive the cartesian postion using triateration (end of file)
    def _actuator_to_cartesian(self, pos):
        return actuator_to_cartesian(self.anchors, pos)
    
    # Returns the current cartesian position of the nozzle
    def get_position(self):
        spos = [s.mcu_stepper.get_commanded_position() for s in self.steppers]
        return self._actuator_to_cartesian(spos)
                         
    # Sets the cartesian position of the nozzle
    def set_position(self, newpos, homing_axes):
        pos = self._cartesian_to_actuator(newpos)
        for i in StepList:
            self.steppers[i].set_position(pos[i])
                         
        self.limit_xy2 = -1.
        if tuple(homing_axes) == StepList:
            self.need_home = False
    
    # Homes the printer                         
    def home(self, homing_state):
                         
        # All axes are homed simultaneously
        homing_state.set_axes([0, 1, 2])
        endstops = [es for s in self.steppers for es in s.get_endstops()]
        s = self.steppers[0] # Assume homing speed same for all steppers
                         
        # Initial homing
        homing_speed = min(s.homing_speed, self.max_z_velocity)
        homepos = [0., 0., self.max_z, None]
        coord = list(homepos)
        coord[2] = -1.5 * math.sqrt(max(self.arm2)-self.max_xy2)
        homing_state.home(coord, homepos, endstops, homing_speed)
                         
        # Retract
        coord[2] = homepos[2] - s.homing_retract_dist
        homing_state.retract(coord, homing_speed)
                         
        # Home again
        coord[2] -= s.homing_retract_dist
        homing_state.home(coord, homepos, endstops,
                          homing_speed/2.0, second_home=True)
                         
        # Set final homed position
        spos = [ep + s.get_homed_offset()
                for ep, s in zip(self.endstops, self.steppers)]
        homing_state.set_homed_position(self._actuator_to_cartesian(spos))
                         
    def motor_off(self, print_time):
        self.limit_xy2 = -1.
        for stepper in self.steppers:
            stepper.motor_enable(print_time, 0)
        self.need_motor_enable = self.need_home = True
        
    def _check_motor_enable(self, print_time):
        for i in StepList:
            self.steppers[i].motor_enable(print_time, 1)
        self.need_motor_enable = False
        
    def check_move(self, move):
        end_pos = move.end_pos
        xy2 = end_pos[0]**2 + end_pos[1]**2
        
        if xy2 <= self.limit_xy2 and not move.axes_d[2]:
            # Normal XY move
            return
        
        if self.need_home:
            raise homing.EndstopMoveError(end_pos, "Must home first")
            
        limit_xy2 = self.max_xy2
        
        # Limit movement to a cone when above limit_z
        if end_pos[2] > self.limit_z:
            limit_xy2 = min(limit_xy2, (self.max_z - end_pos[2])**2)
        
        # Outside operating volume -> throw error
        if xy2 > limit_xy2 or end_pos[2] < self.min_z or end_pos[2] > self.max_z:
            raise homing.EndstopMoveError(end_pos)
            
        if move.axes_d[2]:
            move.limit_speed(self.max_z_velocity, move.accel)
            limit_xy2 = -1.
            
        # Limit the speed/accel of this move if is is at the extreme
        # end of the build envelope
        extreme_xy2 = max(xy2, move.start_pos[0]**2 + move.start_pos[1]**2)
        
        if extreme_xy2 > self.slow_xy2:
            r = 0.5
            if extreme_xy2 > self.very_slow_xy2:
                r = 0.25
            max_velocity = self.max_velocity
            if move.axes_d[2]:
                max_velocity = self.max_z_velocity
            move.limit_speed(max_velocity * r, self.max_accel * r)
            limit_xy2 = -1.
            
        self.limit_xy2 = min(limit_xy2, self.slow_xy2)
        
    def move(self, print_time, move):
        # Useful appendices from move class
        # _d distance (in mm)
        # _r ratio (scalar between 0.0 and 1.0)
        # _v velocity (in mm/second)
        # _v2 is velocity squared (mm^2/s^2)
        # _t time (in seconds)
        
        if self.need_motor_enable:
            self._check_motor_enable(print_time)
        # Distance between end and start position (x, y and z vectors)    
        axes_d = move.axes_d
        # Total combined length of the movements in all three directions
        move_d = move.move_d
        # xy movement ratio (from 0 to 1)
        movexy_r = 1.
        # z move ratio
        movez_r = 0.
        # Inverse of combined length
        inv_movexy_d = 1. / move_d
        
        if not axes_d[0] and not axes_d[1]:
            # Z only move
            movez_r = axes_d[2] * inv_movexy_d
            movexy_r = inv_movexy_d = 0.
            
        elif axes_d[2]:
            # XY+Z move
            movexy_d = math.sqrt(axes_d[0]**2 + axes_d[1]**2)
            movexy_r = movexy_d * inv_movexy_d
            movez_r = axes_d[2] * inv_movexy_d
            inv_movexy_d = 1. / movexy_d

        origx, origy, origz = move.start_pos[:3]

        accel = move.accel
        cruise_v = move.cruise_v
        accel_d = move.accel_r * move_d
        cruise_d = move.cruise_r * move_d
        decel_d = move.decel_r * move_d
         
        for i in StepList:
            # Calculate a virtual tower along the line of movement at
            # the point closest to this stepper's tower.
            towerx_d = self.towers[i][0] - origx
            towery_d = self.towers[i][1] - origy
            vt_startxy_d = (towerx_d*axes_d[0] + towery_d*axes_d[1])*inv_movexy_d
            tangentxy_d2 = towerx_d**2 + towery_d**2 - vt_startxy_d**2
            vt_arm_d = math.sqrt(self.arm2[i] - tangentxy_d2)
            vt_startz = origz

            # Generate steps
            step_delta = self.steppers[i].step_delta
            move_time = print_time
            if accel_d:
                step_delta(move_time, accel_d, move.start_v, accel,
                           vt_startz, vt_startxy_d, vt_arm_d, movez_r)
                vt_startz += accel_d * movez_r
                vt_startxy_d -= accel_d * movexy_r
                move_time += move.accel_t
            if cruise_d:
                step_delta(move_time, cruise_d, cruise_v, 0.,
                           vt_startz, vt_startxy_d, vt_arm_d, movez_r)
                vt_startz += cruise_d * movez_r
                vt_startxy_d -= cruise_d * movexy_r
                move_time += move.cruise_t
            if decel_d:
                step_delta(move_time, decel_d, cruise_v, -accel,
                           vt_startz, vt_startxy_d, vt_arm_d, movez_r)
  
#==== Pasted from V0.3 ===========
    def move(self, move_time, move):
        axes_d = move.axes_d
        move_d = movexy_d = move.move_d
        movexy_r = 1.
        movez_r = 0.
        inv_movexy_d = 1. / movexy_d
        if not axes_d[0] and not axes_d[1]:
            if not axes_d[2]:
                return
            movez_r = axes_d[2] * inv_movexy_d
            movexy_d = movexy_r = inv_movexy_d = 0.
        elif axes_d[2]:
            movexy_d = math.sqrt(axes_d[0]**2 + axes_d[1]**2)
            movexy_r = movexy_d * inv_movexy_d
            movez_r = axes_d[2] * inv_movexy_d
            inv_movexy_d = 1. / movexy_d

        if self.need_motor_enable:
            self._check_motor_enable(move_time)

        origx, origy, origz = move.start_pos[:3]

        accel_t = move.accel_t
        cruise_end_t = accel_t + move.cruise_t
        accel_d = move.accel_r * move_d
        cruise_end_d = accel_d + move.cruise_r * move_d

        inv_cruise_v = 1. / move.cruise_v
        inv_accel = 1. / move.accel
        accel_time_offset = move.start_v * inv_accel
        accel_multiplier = 2.0 * inv_accel
        accel_offset = move.start_v**2 * 0.5 * inv_accel
        decel_time_offset = move.cruise_v * inv_accel + cruise_end_t
        decel_offset = move.cruise_v**2 * 0.5 * inv_accel + cruise_end_d

        for i in StepList:
            # Find point on line of movement closest to tower
            towerx_d = self.towers[i][0] - origx
            towery_d = self.towers[i][1] - origy
            closestxy_d = (towerx_d*axes_d[0] + towery_d*axes_d[1])*inv_movexy_d
            tangentxy_d2 = towerx_d**2 + towery_d**2 - closestxy_d**2
            closest_height2 = self.arm_length2 - tangentxy_d2

            # Calculate accel/cruise/decel portions of move
            reversexy_d = closestxy_d + math.sqrt(closest_height2)*movez_r
            accel_up_d = cruise_up_d = decel_up_d = 0.
            accel_down_d = cruise_down_d = decel_down_d = 0.
            if reversexy_d <= 0.:
                accel_down_d = accel_d
                cruise_down_d = cruise_end_d
                decel_down_d = move_d
            elif reversexy_d >= movexy_d:
                accel_up_d = accel_d
                cruise_up_d = cruise_end_d
                decel_up_d = move_d
            elif reversexy_d < accel_d * movexy_r:
                accel_up_d = reversexy_d * move_d * inv_movexy_d
                accel_down_d = accel_d
                cruise_down_d = cruise_end_d
                decel_down_d = move_d
            elif reversexy_d < cruise_end_d * movexy_r:
                accel_up_d = accel_d
                cruise_up_d = reversexy_d * move_d * inv_movexy_d
                cruise_down_d = cruise_end_d
                decel_down_d = move_d
            else:
                accel_up_d = accel_d
                cruise_up_d = cruise_end_d
                decel_up_d = reversexy_d * move_d * inv_movexy_d
                decel_down_d = move_d

            # Generate steps
            mcu_stepper = self.steppers[i].mcu_stepper
            mcu_time = mcu_stepper.print_to_mcu_time(move_time)
            step_pos = mcu_stepper.commanded_position
            step_dist = self.steppers[i].step_dist
            height = step_pos*step_dist - origz
            if accel_up_d > 0.:
                count = mcu_stepper.step_delta_accel(
                    mcu_time - accel_time_offset, accel_up_d,
                    accel_offset, accel_multiplier, step_dist,
                    height, closestxy_d, closest_height2, movez_r)
                height += count * step_dist
            if cruise_up_d > 0.:
                count = mcu_stepper.step_delta_const(
                    mcu_time + accel_t, cruise_up_d,
                    -accel_d, inv_cruise_v, step_dist,
                    height, closestxy_d, closest_height2, movez_r)
                height += count * step_dist
            if decel_up_d > 0.:
                count = mcu_stepper.step_delta_accel(
                    mcu_time + decel_time_offset, decel_up_d,
                    -decel_offset, -accel_multiplier, step_dist,
                    height, closestxy_d, closest_height2, movez_r)
                height += count * step_dist
            if accel_down_d > 0.:
                count = mcu_stepper.step_delta_accel(
                    mcu_time - accel_time_offset, accel_down_d,
                    accel_offset, accel_multiplier, -step_dist,
                    height, closestxy_d, closest_height2, movez_r)
                height += count * step_dist
            if cruise_down_d > 0.:
                count = mcu_stepper.step_delta_const(
                    mcu_time + accel_t, cruise_down_d,
                    -accel_d, inv_cruise_v, -step_dist,
                    height, closestxy_d, closest_height2, movez_r)
                height += count * step_dist
            if decel_down_d > 0.:
                count = mcu_stepper.step_delta_accel(
                    mcu_time + decel_time_offset, decel_down_d,
                    -decel_offset, -accel_multiplier, -step_dist,
                    height, closestxy_d, closest_height2, movez_r)





    # Helper functions for DELTA_CALIBRATE script
    def get_stable_position(self):
        return [int((ep - s.mcu_stepper.get_commanded_position())
                    /  s.mcu_stepper.get_step_dist() + .5)
                * s.mcu_stepper.get_step_dist()
                for ep, s in zip(self.endstops, self.steppers)]
    
    def get_calibrate_params(self):
        return {
            'endstop_a': self.steppers[0].position_endstop,
            'endstop_b': self.steppers[1].position_endstop,
            'endstop_c': self.steppers[2].position_endstop,
            'angle_a': self.angles[0], 'angle_b': self.angles[1],
            'angle_c': self.angles[2], 'radius': self.radius,
            'arm_a': self.arm_lengths[0], 'arm_b': self.arm_lengths[1],
            'arm_c': self.arm_lengths[2] }


######################################################################
# Matrix helper functions for 3x1 matrices
######################################################################

def matrix_cross(m1, m2):
    return [m1[1] * m2[2] - m1[2] * m2[1],
            m1[2] * m2[0] - m1[0] * m2[2],
            m1[0] * m2[1] - m1[1] * m2[0]]

def matrix_dot(m1, m2):
    return m1[0] * m2[0] + m1[1] * m2[1] + m1[2] * m2[2]

def matrix_magsq(m1):
    return m1[0]**2 + m1[1]**2 + m1[2]**2

def matrix_add(m1, m2):
    return [m1[0] + m2[0], m1[1] + m2[1], m1[2] + m2[2]]

def matrix_sub(m1, m2):
    return [m1[0] - m2[0], m1[1] - m2[1], m1[2] - m2[2]]

def matrix_mul(m1, s):
    return [m1[0]*s, m1[1]*s, m1[2]*s]

def actuator_to_cartesian(anchors, pos):
    # Find nozzle position using trilateration (see wikipedia)
    s21 = matrix_sub(anchors[1], anchors[0])
    s31 = matrix_sub(anchors[2], anchors[0])

    d = math.sqrt(matrix_magsq(s21))
    ex = matrix_mul(s21, 1. / d)
    i = matrix_dot(ex, s31)
    vect_ey = matrix_sub(s31, matrix_mul(ex, i))
    ey = matrix_mul(vect_ey, 1. / math.sqrt(matrix_magsq(vect_ey)))
    ez = matrix_cross(ex, ey)
    j = matrix_dot(ey, s31)

    x = (pos[0]**2 - pos[1]**2 + d**2) / (2. * d)
    y = (pos[0]**2 - pos[2]**2 - x**2 + (x-i)**2 + j**2) / (2. * j)
    z = -math.sqrt(pos[0]**2 - x**2 - y**2)

    ex_x = matrix_mul(ex, x)
    ey_y = matrix_mul(ey, y)
    ez_z = matrix_mul(ez, z)
    return matrix_add(carriage1, matrix_add(ex_x, matrix_add(ey_y, ez_z)))

