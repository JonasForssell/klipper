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
        self.limit_z = min([s.position_endstop - (arm - math.sqrt(arm**2 - radius**2))
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
                        es.position_endstop + math.sqrt(arm**2 - radius**2))
                       for angle, es, arm in zip(self.angles, self.steppers, arm_lengths)]
        
        
        logging.info(
            "Anchor 0: X:%.2fmm Y:%.2fmm Z: %.2fmm)"
            % (self.anchors[0][0], self.anchors[0][1],
               self.anchors[0][2]))
        logging.info(
            "Anchor 1: X:%.2fmm Y:%.2fmm Z: %.2fmm)"
            % (self.anchors[1][0], self.anchors[1][1],
               self.anchors[1][2]))
        logging.info(
            "Anchor 1: X:%.2fmm Y:%.2fmm Z: %.2fmm)"
            % (self.anchors[2][0], self.anchors[2][1],
               self.anchors[2][2]))
        
        
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
    
    # Derive the cartesian postion using triateration (see wikipedia)
    def _actuator_to_cartesian(self, pos):
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
        return matrix_add(anchors[0], matrix_add(ex_x, matrix_add(ey_y, ez_z)))
        
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
        
        # Starting position (in local coordinate system)
        anchors_start = self._cartesian_to_actuator(move.start_pos)
        # Ending position 
        anchors_end = self._cartesian_to_actuator(move.end_pos)
        
        
        # Set up the movement profile consisting of three phases
        #
        #
        #       ______cruise_d____
        #      /                  \
        #     /                    \
        #    /                      \
        #  accel_d              decel_d
        #
        accel = move.accel
        cruise_v = move.cruise_v
        decel = -accel
        accel_d = move.accel_r * move_d
        cruise_d = move.cruise_r * move_d
        decel_d = move.decel_r * move_d
        
        # Now, go though each anchor and add required step times in order.
        for i in StepList:
            # Set up the step method
            step = self.steppers[i].step
            # Reset move time
            move_time = print_time
            # Reset move position
            current_pos_r = 0
            previous_pos_r = 0
            current_stepper_pos = 0
            
            beyond_reversal_point = False
            
            # Determine the vector from start point to anchor
            anchor_d = matrix_sub(self.anchors[i], move.start_pos)
            
            # Determine stepping direction for stepper position
            stepper_start_pos = self._cartesian_to_actuator(move.start_pos)
            # Calculate stepper position at the end of the move
            stepper_end_pos = self._cartesian_to_actuator(move.end_pos)
            # Check if it is moving in the right direction
            if stepper_end_pos > stepper_start_pos:
                stepper_step_distance = self.steppers[i].step_dist
            else:
                stepper_step_distance = -self.steppers[i].step_dist
            
            # Calculate reversal point if the effector passes it
            # This is achieved by orthogonal projection of the anchor point onto the line of movement.
            # https://en.wikibooks.org/wiki/Linear_Algebra/Orthogonal_Projection_Onto_a_Line
            reversal_point = matrix_dot(anchor_d, axes_d) / matrix_dot(axes_d, axes_d)
            
            # Express this in stepper coordinates
            stepper_reversal_point = math.sqrt(matrix_magsq(matrix_sub(anchor_d, matrix_mul(axes_d, reversal_point))))
            
            # Now walk along the line one step at a time and plot the time as we go along
            # Phase 1: Movement with acceleration
            while (current_pos_r < move_d):
                # Take one step on the stepper
                current_stepper_pos += stepper_step_distance                                             
                
                # Reverse step direction if we have gone past reversal point. Note, there are only min reversal points                                             
                if current_stepper_pos < stepper_reversal_point:
                    stepper_step_distance = -stepper_step_distance
                    # remember we are on the other side
                    beyond_reversal_point = True                                         
                    # and go the other way instead                                         
                    current_stepper_pos += current_stepper_pos                                         
                    current_stepper_pos += current_stepper_pos                                         
                
                # Calculate the effector position
                previous_pos_r = current_pos_r
                current_pos_r = _movement_position_from_stepper_pos(beyond_reversal_point, move.axes_d, anchor_d, current_stepper_pos)
                
                # Calculate corresponding time depending on which phase we are in
                if current_pos_r > (accel_d + cruise_d):                                             
                    move_time += math.sqrt((current_pos_r-previous_pos_r)*move_d/decel)                                         
                elif current_pos_r > accel_d:
                    move_time += (current_pos_r-previous_pos_r)/cruise_v                                         
                else:                                             
                    move_time += math.sqrt((current_pos_r-previous_pos_r)*move_d/accel)
                
                # Push time on stack
                step(move_time)
            
            # Now, repeat this for all steppers
        # All times have been pushed
        
        
    # Find the current position along line of movement which fulfils the stepper position
    # Mathematical problem is
    # P is starting point
    # A is anchor point
    # V is the (unit)vector of the line which the effector movement follows                       
    # Q is the new position when the effector has moved one step
    #
    #    A                    
    #    .                                  
    #    .  .  LS                         
    #  L .    .             ....> V          
    #    .      .    .......  
    #    .     ...Q..
    #    ......                    
    #  P      M                
    #                        
    # Next, we know the following
    # L is the line length between anchor and starting point, i.e. length of vector PA
    # LS is the line length between anchor and next point, i.e. length of vector QA 
    # M is the length between starting point and next point i.e. length of vector PQ
    #
    # We also know that Q should be along the line of movement so Q = P + M*V                       
    #
    # Let APx = Ax - Px and so on                        
    # If we set up the length of vector QA, and we substitue coordinates of Q with the equation above we get
    #
    # LS = SQRT( (APx - M*Vx)^2 + (APy - M*Vy)^2 + (APz - M*Vz)^2 )
    # LS is the same as the current stepper position (in local coordinates)
    #                        
    # Solving this equation for M gives a long expression                        
    #                        
    # M = (0.5*SQRT( (-2*APx*Vx - 2*APy*Vy - 2*APz*Vz)^2 -4*(Vx^2+Vy^2+Vz^2)*(APx^2+APy^2+APz^2-LS^2)) 
    #               + APx*Vx + APy*Vy + APz*Vz ) / (Vx^2 + Vy^2 + Vz^2)                       
    #
    # We now know the position of the effector for each step on the stepper
    # This can be used to calculate at what time we should take each step on the stepper.
    #
    def _movement_position_from_stepper_pos(beyond_reversal_point, V, PA, current_stepper_pos):
        PAx = PA[0]       
        PAy = PA[1]       
        PAz = PA[2]       
        
        Vx = V[0]
        Vy = V[1]
        Vz = V[2]        
        
        if beyond_reversal_point:
            return (0.5*math.sqrt( (-2*PAx*Vx - 2*PAy*Vy - 2*PAz*Vz)**2 -4*(Vx**2+Vy**2+Vz**2)*(PAx**2+PAy**2+PAz**2 - current_stepper_pos**2)) 
                   + PAx*Vx + PAy*Vy + PAz*Vz ) / (Vx**2 + Vy**2 + Vz**2)
        else:
            return (-0.5*math.sqrt( (-2*PAx*Vx - 2*PAy*Vy - 2*PAz*Vz)**2 -4*(Vx**2+Vy**2+Vz**2)*(PAx**2+PAy**2+PAz**2 - current_stepper_pos**2)) 
                   + PAx*Vx + PAy*Vy + PAz*Vz ) / (Vx**2 + Vy**2 + Vz**2)
    
        
    
    
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

