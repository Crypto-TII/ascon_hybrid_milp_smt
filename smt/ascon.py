# Copyright 2022 Rusydi H. Makarim and Raghvendra Rohit
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import timeit
import stp


class Trail:
    def __init__(self, trail):
        """
        Construct an instance of Trail for Ascon

        INPUT:

        - ``trail`` -- a list of Ascon state
        """
        if any(len(state) != AsconBase.nwords for state in trail):
            raise ValueError("the length of each element of trail must be equal to %d" % AsconBase.nwords)

        self._trail = trail

    def trail(self):
        """
        Return the trail
        """
        return self._trail

    def nrounds(self):
        """
        Return the no. of rounds that this trail represents
        """
        return len(self.trail()) // 2

    def input_values(self):
        """
        Return the input values of the trail
        """
        return self.trail()[0]

    def output_values(self):
        """
        Return the output values of the trail
        """
        return self.trail()[-1]

    def __repr__(self):
        return "%d-round differential trail of Ascon" % self.nrounds()

    def __str__(self):
        def format_string(s):
            return bin(s)[2:].rjust(AsconBase.word_size, '.').replace('0', '.')

        ret = []
        trail = self.trail()
        for r in range(len(trail)-1):
            x = trail[r]
            s = tuple(format_string(x[i]) for i in range(AsconBase.nwords))
            ret.append("%s" % '\n'.join(s))
            ret.append('-'*64)

            if (r % 2) == 0:
                ret[-1] += ' pS'
            else:
                ret[-1] += ' pL'

        x = self.output_values()
        s = tuple(format_string(x[i]) for i in range(AsconBase.nwords))
        ret.append("%s" % '\n'.join(s))

        return '\n'.join(ret)


class AsconBase(object):
    word_size = 64
    nwords = 5
    size = word_size * nwords

    _input_nonlinear_var_name = 'x'
    _input_linear_var_name = 'y'

    rotation_constants = ((19, 28), (61, 39), (1, 6), (10, 17), (7, 41))

    def __init__(self, nrounds, remove_outer_layer=True):
        """
        Construct an instance of AsconBase

        INPUTS::

        - ``nrounds`` -- no. of rounds
        - ``remove_outer_layer`` -- whether to remove the outer layer in the model

        NOTES::

            The outer layer refers to the first round nonlinear function and the last round function. This improves
            the runtime of the STP solver.
        """
        self._nrounds = nrounds
        self._solver = stp.Solver()
        self._additional_constraints = []
        self._solving_time = None
        self._is_outer_layer_removed = remove_outer_layer

    @property
    def nrounds(self):
        """
        Return the no. of rounds
        """
        return self._nrounds

    @property
    def solver(self):
        """
        Return the STP solver
        """
        return self._solver

    @property
    def additional_constraints(self):
        """
        Return a list of additional constraints
        """
        return self._additional_constraints

    @staticmethod
    def rotate_right(w, c):
        """
        Perform bitwise right rotation on the word `w`

        INPUT:

        - ``w`` -- input word
        - ``c`` -- rotation constant
        """
        if not 0 <= c < AsconBase.word_size:
            raise ValueError("c must be in the range 0 <= c < %d" % AsconBase.word_size)
        return (w >> c) | (w << (AsconBase.word_size - c)) & ((1 << AsconBase.word_size) - 1)

    @staticmethod
    def rotate_left(w, c):
        """
        Perform bitwise left rotation on the word 'w'

        INPUT:

        - ``w`` -- input word
        - ``c`` -- rotation constant
        """
        return AsconBase.rotate_right(w, AsconBase.word_size - c) & ((1 << AsconBase.word_size) - 1)

    def is_outer_layer_removed(self):
        """
        Return `True` if the outer layer are removed in the model and `False` otherwise.
        """
        return self._is_outer_layer_removed

    def vars(self, var_name, r):
        """
        Return the a variable

        - ``var_name`` -- name of variables
        - ``r`` -- round index
        """
        if not 0 <= r <= self.nrounds:
            raise ValueError("r must be in the range 0 <= r < %d" % self.nrounds)
        return [self.solver.bitvec("%s%02d%02d" % (var_name, r, i), width=AsconBase.word_size) for i in range(AsconBase.nwords)]

    def set_nactive_sboxes_at_round(self, r, w):
        """
        Set the no. of active SBox at round `r` equal to `w`

        INPUT:

        - ``r`` -- round index
        - ``w`` -- a positive integer
        """
        self.additional_constraints.append(self.nactive_sboxes_var(r) == w)

    def set_min_nactive_sboxes_at_round(self, r, w):
        """
        Set the minimum no. of active SBox at around `r` to be greater or equal to `w`

        INPUT:

        - ``r`` -- round index
        - ``w`` -- a positive integer
        """
        self.additional_constraints.append(self.nactive_sboxes_var(r) >= w)

    def set_max_nactive_sboxes_at_round(self, r, w):
        """
        Set the maximum no. of active SBox at round `r` to be less than or equal to `w`

        INPUT:

        - ``r`` -- round index
        - ``w`` -- a positive integer
        """
        self.additional_constraints.append(self.nactive_sboxes_var(r) <= w)

    def set_active_sboxes_at_round(self, r, active_sboxes):
        """
        Set the active S-Boxes at round `r`

        INPUT:

        - ``r`` -- round index
        - ``active_sboxes`` -- a list/tuple of indices where S-Boxes are active
        """
        if any(not 0 <= index < AsconBase.word_size for index in active_sboxes):
            raise ValueError("index of active sboxes must be in the range 0 <= index < %d" % AsconBase.word_size)

        x = self.input_linear_vars(0) if r == 0 and self.is_outer_layer_removed() else self.input_nonlinear_vars(r)
        w = 0
        for index in active_sboxes:
            w |= (1 << index)
        self.additional_constraints.append((x[0] | x[1] | x[2] | x[3] | x[4]) == w)

    def set_input_nonlinear_state_at_round(self, r, x0, x1, x2, x3, x4):
        """
        Set the input of the nonlinear function at round `r`

        INPUT:

        - ``r`` -- round index
        - ``x0, x1, x2, x3, x4`` -- input words
        """
        X = self.input_nonlinear_vars(r)
        x = (x0, x1, x2, x3, x4)
        self._additional_constraints += [X[i] == x[i] for i in range(AsconBase.nwords)]

    def set_output_nonlinear_state_at_round(self, r, x0, x1, x2, x3, x4):
        """
        Set the output of the nonlinear function at round `r`

        INPUT:

        - ``r`` -- round index
        - ``x0, x1, x2, x3, x4`` -- input words
        """
        self.set_input_linear_state_at_round(r, x0, x1, x2, x3, x4)

    def set_input_linear_state_at_round(self, r, x0, x1, x2, x3, x4):
        """
        Set the input of the linear function at round `r`

        INPUT:

        - ``r`` -- round index
        - ``x0, x1, x2, x3, x4`` -- input words
        """
        X = self.input_linear_vars(r)
        x = (x0, x1, x2, x3, x4)
        self._additional_constraints += [X[i] == x[i] for i in range(AsconBase.nwords)]

    def set_output_linear_state_at_round(self, r, x0, x1, x2, x3, x4):
        """
        Set the output of the linear function at round `r`

        INPUT:

        - ``r`` -- round index
        - ``x0, x1, x2, x3, x4`` -- input words
        """
        self.set_input_nonlinear_state_at_round(r + 1, x0, x1, x2, x3, x4)

    def set_input_round_state(self, r, x0, x1, x2, x3, x4):
        """
        Set the input of the `r`-th round function

        INPUT:

        - ``r`` -- round index
        - ``x0, x1, x2, x3, x4`` -- input words
        """
        self.set_input_nonlinear_state_at_round(r, x0, x1, x2, x3, x4)

    def set_output_round_state(self, r, x0, x1, x2, x3, x4):
        """
        Set the output of the `r`-th round function

        INPUT:

        - ``r`` -- round index
        - ``x0, x1, x2, x3, x4`` -- input words
        """
        self.set_input_round_state(r + 1, x0, x1, x2, x3, x4)

    def solving_time(self):
        """
        Return the time to run the STP solver in seconds
        """
        if self._solving_time is None:
            raise RuntimeError("the method has_trail() needs to be executed first")

        return self._solving_time

    def active_sboxes_at_round(self, r):
        """
        Return the active SBoxes at round `r`

        INPUT:

        - ``r`` -- round index
        """
        x = self.input_linear_values(0) if r == 0 and self.is_outer_layer_removed() else self.input_nonlinear_values(r)
        w = x[0] | x[1] | x[2] | x[3] | x[4]
        return tuple(i for i in range(AsconBase.word_size) if w & (1 << i))

    def nactive_sboxes_at_round(self, r):
        """
        Return the number of active S-Boxes at round `r`

        INPUT:

        - ``r`` -- round index
        """
        return len(self.active_sboxes_at_round(r))

    def nactive_sboxes(self):
        """
        Return the total number of active S-Boxes
        """
        return sum(self.nactive_sboxes_at_round(r) for r in range(self.nrounds))

    def nactive_sboxes_var(self, r):
        """
        Return the variable representing the no. of active S-Box at round `r`

        INPUT:

        - ``r`` -- round index
        """
        return self.solver.bitvec("w%02d" % r, width=AsconBase.word_size)

    def nactive_sboxes_vars(self):
        """
        Return a list of variables representing the no. of active S-Box for each round
        """
        return [self.nactive_sboxes_var(r) for r in range(self.nrounds)]

    def nactive_sboxes_constraint_using_input_nonlinear_vars(self, r):
        """
        Return a constraint for the no. of active S-Boxes using the input variables of the nonlinear layer at round `r`

        INPUT:

        - ``r`` -- round index
        """
        return self._nactive_sboxes_constraint_using(self.input_nonlinear_vars(r), self.nactive_sboxes_var(r))

    def nactive_sboxes_constraint_using_output_nonlinear_vars(self, r):
        """
        Return a constraint for the no. of active S-Boxes using the output variables of the nonlinear layer at round `r`

        INPUT:

        - ``r`` -- round index
        """
        return self._nactive_sboxes_constraint_using(self.output_nonlinear_vars(r), self.nactive_sboxes_var(r))

    def nactive_sboxes_constraint_using_input_linear_vars(self, r):
        """
        Return a constraint for the no. of active S-Boxes using the input variables of the linear layer at round `r`

        INPUT:

        - ``r`` -- round index
        """
        return self._nactive_sboxes_constraint_using(self.input_linear_vars(r), self.nactive_sboxes_var(r))

    def nactive_sboxes_constraint_using_output_linear_vars(self, r):
        """
        Return a constraint for the no. of active S-Boxes using the output variables of the linear layer at round `r`

        INPUT:

        - ``r`` -- round index
        """
        return self._nactive_sboxes_constraint_using(self.output_linear_vars(r), self.nactive_sboxes_var(r))

    def _nactive_sboxes_constraint_using(self, variables, nactive_sboxes_var):
        """
        Return a constraint representing the no. of active S-Boxes

        INPUT:

        - ``variables`` -- a list of variables of an Ascon state
        - ``nactive_sboxes_var`` -- a variable for the number of active S-Boxes
        """
        x = variables
        w = nactive_sboxes_var
        return self.hamming_weight_constraint(x[0] | x[1] | x[2] | x[3] | x[4], w)

    def input_nonlinear_vars(self, r):
        """
        Return a list of input variables for the nonlinear layer at round `r`

        INPUT:

        - ``r`` -- round index
        """
        return self.vars(AsconBase._input_nonlinear_var_name, r)

    def input_linear_vars(self, r):
        """
        Return a list of input variables for the linear layer at round `r`

        INPUT:

        - ``r`` -- round index
        """
        return self.vars(AsconBase._input_linear_var_name, r)

    def output_linear_vars(self, r):
        """
        Return a list of the output variables for the linear layer at round `r`

        INPUT:

        - ``r`` -- round index
        """
        return self.input_nonlinear_vars(r + 1)

    def output_nonlinear_vars(self, r):
        """
        Return a list of the output variables for the nonlinear layer at round `r`

        INPUT:

        - ``r`` -- round index
        """
        return self.input_linear_vars(r)

    def nonlinear_layer_constraints(self, r):
        """
        Return a list of constraints for the nonlinear layer at round `r`

        INPUT:

        - ``r`` -- round index
        """
        x = self.input_nonlinear_vars(r)
        y = self.output_nonlinear_vars(r)

        def column(i, state):
            ret = [((row & (1 << i)) >> i) for row in state]
            ret.reverse()  # the bottom word is the least-significant bit to the SBox
            return ret

        constraints = []
        for i in range(AsconBase.word_size):
            xi = column(i, x)
            yi = column(i, y)
            constraints += self.sbox_constraints(xi, yi)

        return constraints

    def round_constraints(self, r):
        """
        Return a list of constraints for the `r`-th round

        INPUT:

        - ``r`` -- round index
        """
        return self.nonlinear_layer_constraints(r) + self.linear_layer_constraints(r)

    def permutation_constraints(self):
        """
        Return a list of constraints the Ascon permutation
        """
        if self.is_outer_layer_removed():
            constraints = self.linear_layer_constraints(0)
            constraints += sum([self.round_constraints(r) for r in range(1, self.nrounds - 1)], [])
        else:
            constraints = sum([self.round_constraints(r) for r in range(self.nrounds)], [])

        return constraints

    def total_nactive_sboxes_constraints(self):
        """
        Return a list of constraints for the total no. of active S-Box
        """
        if self.is_outer_layer_removed():
            constraints = [self.nactive_sboxes_constraint_using_input_linear_vars(0)]
            constraints += [self.nactive_sboxes_constraint_using_input_nonlinear_vars(r) for r in range(1, self.nrounds)]
        else:
            constraints = [self.nactive_sboxes_constraint_using_input_nonlinear_vars(r) for r in range(self.nrounds)]

        return constraints

    def has_trail(self, nactive_sboxes=None, timeout=None):
        """
        Return `True` if Ascon has a trail

        INPUT:

        - ``nactive_sboxes`` -- no. of active S-Boxes (default: None)
        - ``timeout`` -- specified timeout to run STP in seconds (default: None)
        """
        constraints = self.permutation_constraints() + self.total_nactive_sboxes_constraints()
        constraints += self.additional_constraints

        if nactive_sboxes is not None:
            if not isinstance(nactive_sboxes, int):
                raise TypeError("nactive_sboxes must be an int")

            w = self.nactive_sboxes_vars()
            constraints += [sum(w) == nactive_sboxes]

        for constraint in constraints:
            self.solver.add(constraint)

        max_time = -1 if timeout is None else timeout

        start_time = timeit.default_timer()
        output_int = self.solver.check_with_timeout(max_time=max_time)
        self._solving_time = timeit.default_timer() - start_time

        # See the definition of the result in
        # https://github.com/stp/stp/blob/9a59a72e82d67cefeb88d8baa34965f70acb5d1c/lib/Interface/c_interface.cpp#L520
        result = None
        if output_int == 0:
            result = True
        elif output_int == 1:
            result = False
        elif output_int == 2:
            result = "error"
        elif output_int == 3:
            result = "timeout"

        return result

    def input_vars(self):
        """
        Return a list of input variables
        """
        return self.input_nonlinear_vars(0)

    def output_vars(self):
        """
        Return a list of ouput variables
        """
        return self.input_nonlinear_vars(self.nrounds)

    def output_values(self):
        """
        Return a list of output values of the permutation
        """
        return [x.value for x in self.output_vars()]

    def input_nonlinear_values(self, r):
        """
        Return a list of input values to the nonlinear layer at round `r`

        INPUT:

        - ``r`` -- round index
        """
        return [x.value for x in self.input_nonlinear_vars(r)]

    def output_nonlinear_values(self, r):
        """
        Return a list of output values to the nonlinear layer at round `r`

        INPUT:

        - ``r`` -- round index
        """
        return self.input_linear_values(r)

    def input_linear_values(self, r):
        """
        Return a list of input values to the linear layer at round `r`

        INPUT:

        - ``r`` -- round index
        """
        return [x.value for x in self.input_linear_vars(r)]

    def output_linear_values(self, r):
        """
        Return a list of output values to the linear layer at round `r`

        INPUT:

        - ``r`` -- round index
        """
        return self.input_nonlinear_values(r + 1)

    def trail_at_round(self, r):
        """
        Return the trail at round `r`

        INPUT:

        - ``r`` -- round index
        """
        return [self.input_nonlinear_values(r), self.input_linear_values(r)]

    def trail(self):
        """
        Return the trail
        """
        trail = sum([self.trail_at_round(r) for r in range(self.nrounds)], [])
        trail += [self.output_values()]
        return Trail(trail)

    def hamming_weight_constraint(self, x, weight):
        """
        Return a constraint represent the Hamming weight of `x`

        INPUT:

        - ``x`` -- a word
        - ``weight`` -- the specified Hamming weight
        """
        return sum((x & (1 << i)) >> i for i in range(AsconBase.word_size)) == weight

    def sbox_constraints(self, x, y):
        """
        Return a list of constraints that describe propagation of Ascon's SBox

        INPUT:

        - ``x`` -- input variables
        - ``y`` -- output variables
        """
        raise NotImplementedError

    def linear_layer_constraints(self, r):
        """
        Return a list of constraints that describe propagation of Ascon's linear layer

        INPUT:

        - ``r`` -- round index
        """
        raise NotImplementedError

    def __repr__(self):
        return "%d-rounds of %s" % (self.nrounds, self.__class__.__name__)


class AsconDifferential(AsconBase):
    def __init__(self, nrounds, remove_outer_layer=True):
        """
        Construct an instance of AsconDifferential

        INPUT:

        - ``nrounds`` -- no. of rounds
        - ``remove_outer_layer`` -- whether to remove the outer layer in the model (default: True)
        """
        super(AsconDifferential, self).__init__(nrounds, remove_outer_layer)

    def linear_layer_constraints(self, r):
        """
        Return a list of constraints for the linear layer at round `r`

        INPUT:

        - ``r`` -- round index
        """
        x = self.input_linear_vars(r)
        y = self.output_linear_vars(r)
        rc = AsconBase.rotation_constants
        rotate_right = AsconBase.rotate_right

        constraints = [y[i] == x[i] ^ rotate_right(x[i], rc[i][0]) ^ rotate_right(x[i], rc[i][1])
                       for i in range(AsconBase.nwords)]

        return constraints

    def sbox_constraints(self, x, y):
        """
        Return a list of constraints that describe differential propagation of Ascon's SBox

        INPUT:

        - ``x`` -- a list of input variables
        - ``y`` -- a list of output variables
        """
        if len(x) != AsconBase.nwords:
            raise ValueError("the length of x must be equal to %d" % AsconBase.nwords)

        if len(y) != AsconBase.nwords:
            raise ValueError("the length of y must be equal to %d" % AsconBase.nwords)

        constraints = []
        constraints += [((~x[0] & ~x[1] & ~x[2] & ~x[3] & ~x[4]) & (y[0] | y[1] | y[2] | y[3] | y[4])) == 0]
        constraints += [((x[0] & ~x[1] & ~x[2] & ~x[3] & ~x[4]) & (~y[3] | ~(y[0] ^ y[4]))) == 0]
        constraints += [((~x[0] & x[1] & ~x[2] & ~x[3] & ~x[4]) & (~y[0] | ~y[4])) == 0]
        constraints += [((x[0] & x[1] & ~x[2] & ~x[3] & ~x[4]) & (y[1] | ~(y[0] ^ y[4]))) == 0]
        constraints += [((~x[0] & ~x[1] & x[2] & ~x[3] & ~x[4]) & (y[0] | ~y[1] | ~y[2])) == 0]
        constraints += [((x[0] & ~x[1] & x[2] & ~x[3] & ~x[4]) & (~y[4] | ~(y[0] ^ y[2] ^ y[3]))) == 0]
        constraints += [((~x[0] & x[1] & x[2] & ~x[3] & ~x[4]) & ~y[0]) == 0]
        constraints += [((x[0] & x[1] & x[2] & ~x[3] & ~x[4]) & (~y[1] | y[4])) == 0]
        constraints += [((~x[0] & ~x[1] & ~x[2] & x[3] & ~x[4]) & (~y[1] | ~y[2])) == 0]
        constraints += [((x[0] & ~x[1] & ~x[2] & x[3] & ~x[4]) & ~(y[0] ^ y[2] ^ y[3] ^ y[4])) == 0]
        constraints += [((~x[0] & x[1] & ~x[2] & x[3] & ~x[4]) & ~(y[0] ^ y[1] ^ y[2])) == 0]
        constraints += [((x[0] & x[1] & ~x[2] & x[3] & ~x[4]) & ~y[1]) == 0]
        constraints += [((~x[0] & ~x[1] & x[2] & x[3] & ~x[4]) & (y[1] | y[2] | ~(y[0] ^ y[3] ^ y[4]))) == 0]
        constraints += [((x[0] & ~x[1] & x[2] & x[3] & ~x[4]) & ~(y[0] ^ y[3] ^ y[4])) == 0]
        constraints += [((~x[0] & x[1] & x[2] & x[3] & ~x[4]) & (y[3] | ~(y[0] ^ y[1] ^ y[2]))) == 0]
        constraints += [((x[0] & x[1] & x[2] & x[3] & ~x[4]) & (y[1] | ~y[3])) == 0]
        constraints += [((~x[0] & ~x[1] & ~x[2] & ~x[3] & x[4]) & (~y[3] | y[2] | ~(y[0] ^ y[4]))) == 0]
        constraints += [((x[0] & ~x[1] & ~x[2] & ~x[3] & x[4]) & (~y[0] | y[3] | ~y[4])) == 0]
        constraints += [((~x[0] & x[1] & ~x[2] & ~x[3] & x[4]) & ~(y[0] ^ y[4])) == 0]
        constraints += [((x[0] & x[1] & ~x[2] & ~x[3] & x[4]) & (y[0] | y[4] | ~(y[1] ^ y[2]))) == 0]
        constraints += [((~x[0] & ~x[1] & x[2] & ~x[3] & x[4]) & (~y[2] | y[4])) == 0]
        constraints += [((x[0] & ~x[1] & x[2] & ~x[3] & x[4]) & (~y[0] | ~(y[2] ^ y[3] ^ y[4]))) == 0]
        constraints += [((~x[0] & x[1] & x[2] & ~x[3] & x[4]) & ~y[4]) == 0]
        constraints += [((x[0] & x[1] & x[2] & ~x[3] & x[4]) & (y[0] | ~(y[1] ^ y[2]))) == 0]
        constraints += [((~x[0] & ~x[1] & ~x[2] & x[3] & x[4]) & ~y[2]) == 0]
        constraints += [((x[0] & ~x[1] & ~x[2] & x[3] & x[4]) &
                         ~((y[1] & (y[0] ^ y[2])) ^ ((y[3] ^ y[4]) & ~(y[0] ^ y[1] ^ y[2])))) == 0]
        constraints += [((~x[0] & x[1] & ~x[2] & x[3] & x[4]) & ~(y[0] ^ y[1] ^ y[3])) == 0]
        constraints += [((x[0] & x[1] & ~x[2] & x[3] & x[4]) & ~(y[1] ^ y[2])) == 0]
        constraints += [((~x[0] & ~x[1] & x[2] & x[3] & x[4]) & (y[2] | ~(y[0] ^ y[3] ^ y[4]))) == 0]
        constraints += [((x[0] & ~x[1] & x[2] & x[3] & x[4]) &
                         ~((y[0] & (y[1] ^ y[2])) ^ ((y[3] ^ y[4]) & ~(y[0] ^ y[1] ^ y[2])))) == 0]
        constraints += [((~x[0] & x[1] & x[2] & x[3] & x[4]) & ~y[3]) == 0]
        constraints += [((x[0] & x[1] & x[2] & x[3] & x[4]) & (y[3] | ~(y[2] ^ y[1]))) == 0]

        return constraints


class AsconLinear(AsconBase):
    def __init__(self, nrounds, remove_outer_layer=True):
        """
        Construct an instance of AsconLinear

        INPUT:

        - ``nrounds`` -- no. of rounds
        - ``remove_outer_layer`` -- whether to remove the outer layer in the model (default: True)
        """
        super(AsconLinear, self).__init__(nrounds, remove_outer_layer)

    def linear_layer_constraints(self, r):
        """
        Return a list of constraints for the linear layer at round `r`

        INPUT:

        - ``r`` -- round index
        """
        x = self.input_linear_vars(r)
        y = self.output_linear_vars(r)
        rc = AsconBase.rotation_constants
        rotate_left = AsconBase.rotate_left

        constraints = [x[i] == y[i] ^ rotate_left(y[i], rc[i][0]) ^ rotate_left(y[i], rc[i][1])
                       for i in range(AsconBase.nwords)]

        return constraints

    def sbox_constraints(self, x, y):
        """
        Return a list of constraints that describe linear propagation of Ascon's SBox

        INPUT:

        - ``x`` -- a list of input variables
        - ``y`` -- a list of output variables
        """
        if len(x) != AsconBase.nwords:
            raise ValueError("the length of x must be equal to %d" % AsconBase.nwords)

        if len(y) != AsconBase.nwords:
            raise ValueError("the length of y must be equal to %d" % AsconBase.nwords)

        constraints = []
        constraints += [((~x[0] & ~x[1] & ~x[2] & ~x[3] & ~x[4]) & (y[0] | y[1] | y[2] | y[3] | y[4])) == 0]
        constraints += [((x[0] & ~x[1] & ~x[2] & ~x[3] & ~x[4]) & ~((y[0] & y[1] & y[2] & y[3] & y[4]) ^ (y[0] & y[1] & y[2] & y[3]) ^ (y[0] & y[1] & y[2] & y[4]) ^ (y[1] & y[2] & y[3] & y[4]) ^ (y[0] & y[1] & y[2]) ^ (y[1] & y[2] & y[3]) ^ (y[1] & y[2] & y[4]) ^ (y[1] & y[2]) ^ (y[0] & y[3]) ^ (y[1] & y[3]) ^ (y[1] & y[4]) ^ (y[3] & y[4]))) == 0]
        constraints += [((~x[0] & x[1] & ~x[2] & ~x[3] & ~x[4]) & ~(y[1] & ((y[2] & y[3] & y[4]) ^ (y[2] & y[3]) ^ (y[2] & y[4]) ^ y[2] ^ y[3] ^ y[4]))) == 0]
        constraints += [((x[0] & x[1] & ~x[2] & ~x[3] & ~x[4]) & ~((y[0] ^ y[4]) & ~((y[1] & y[2] & y[3]) ^ (y[1] & y[2]) ^ (y[1] & y[3]) ^ (y[2] & y[3]) ^ y[1] ^ y[2]))) == 0]
        constraints += [((~x[0] & ~x[1] & x[2] & ~x[3] & ~x[4]) & ~((y[0] & y[1] & y[2] & y[3] & y[4]) ^ (y[0] & y[1] & y[3] & y[4]) ^ (y[0] & y[2] & y[3] & y[4]) ^ (y[0] & y[3] & y[4]) ^ (y[0] & y[1]) ^ (y[0] & y[2]) ^ (y[1] & y[3]) ^ (y[1] & y[4]) ^ (y[2] & y[4]))) == 0]
        constraints += [((x[0] & ~x[1] & x[2] & ~x[3] & ~x[4]) & ~((y[0] & y[1] & y[2] & y[3] & y[4]) ^ (y[0] & y[2] & y[3] & y[4]) ^ (y[0] & y[1]) ^ (y[0] & y[2]) ^ (y[0] & y[3]) ^ (y[2] & y[4]) ^ (y[3] & y[4]))) == 0]
        constraints += [((~x[0] & x[1] & x[2] & ~x[3] & ~x[4]) & ~((y[0] & y[1] & y[2] & y[3] & y[4]) ^ (y[0] & y[1] & y[3] & y[4]) ^ (y[0] & y[2] & y[3] & y[4]) ^ (y[0] & y[3] & y[4]) ^ (y[0] & y[1]) ^ (y[0] & y[2]) ^ (y[2] & y[4]))) == 0]
        constraints += [((x[0] & x[1] & x[2] & ~x[3] & ~x[4]) & ~((y[0] & y[1] & y[2] & y[3] & y[4]) ^ (y[0] & y[2] & y[3] & y[4]) ^ (y[0] & y[1]) ^ (y[0] & y[2]) ^ (y[0] & y[3]) ^ (y[1] & y[3]) ^ (y[1] & y[4]) ^ (y[2] & y[4]) ^ (y[3] & y[4]))) == 0]
        constraints += [((~x[0] & ~x[1] & ~x[2] & x[3] & ~x[4]) & ~(y[3]  &  ((y[0] & y[1] & y[4]) ^ (y[0] & y[4]) ^ y[1]))) == 0]
        constraints += [((x[0] & ~x[1] & ~x[2] & x[3] & ~x[4]) & ~((y[0] & y[1] & y[2] & y[3] & y[4]) ^ (y[0] & y[1] & y[2] & y[3]) ^ (y[0] & y[1] & y[2] & y[4]) ^ (y[0] & y[1] & y[2]) ^ (y[0] & y[3]) ^ (y[1] & y[4]) ^ (y[3] & y[4]))) == 0]
        constraints += [((~x[0] & x[1] & ~x[2] & x[3] & ~x[4]) & ~(y[4]  &  ((y[0] & y[1] & y[3]) ^ (y[0] & y[3]) ^ y[1]))) == 0]
        constraints += [((x[0] & x[1] & ~x[2] & x[3] & ~x[4]) & ~((y[0] & y[1] & y[2] & y[3]) ^ (y[1] & y[2] & y[3] & y[4]) ^ (y[0] & y[1] & y[2]) ^ (y[0] & y[1] & y[3]) ^ (y[0] & y[2] & y[3]) ^ (y[1] & y[2] & y[4]) ^ (y[1] & y[3] & y[4]) ^ (y[2] & y[3] & y[4]) ^ (y[0] & y[1]) ^ (y[0] & y[2]) ^ (y[1] & y[3]) ^ (y[1] & y[4]) ^ (y[2] & y[4]) ^ y[0] ^ y[4]) )== 0]
        constraints += [((~x[0] & ~x[1] & x[2] & x[3] & ~x[4]) & ~((y[1] ^ y[2])  &  ~((y[0] & y[3] & y[4]) ^ (y[0] & y[3]) ^ (y[0] & y[4]) ^ (y[3] & y[4]) ^ y[3]))) == 0]
        constraints += [((x[0] & ~x[1] & x[2] & x[3] & ~x[4]) & ~((y[0] & y[1] & y[2] & y[3] & y[4]) ^ (y[0] & y[1] & y[2] & y[3]) ^ (y[0] & y[1] & y[2] & y[4]) ^ (y[0] & y[2] & y[3] & y[4]) ^ (y[1] & y[2] & y[3] & y[4]) ^ (y[0] & y[1] & y[2]) ^ (y[0] & y[2] & y[3]) ^ (y[1] & y[2] & y[3]) ^ (y[0] & y[2] & y[4]) ^ (y[1] & y[2] & y[4]) ^ (y[2] & y[3] & y[4]) ^ (y[0] & y[1]) ^ (y[1] & y[2]) ^ (y[0] & y[3]) ^ (y[1] & y[3]) ^ (y[2] & y[3]) ^ (y[3] & y[4]) ^ y[2])) == 0]
        constraints += [((~x[0] & x[1] & x[2] & x[3] & ~x[4]) & ~((y[0] & y[1] & y[2] & y[3] & y[4]) ^ (y[0] & y[1] & y[2] & y[3]) ^ (y[0] & y[1] & y[2] & y[4]) ^ (y[0] & y[2] & y[3] & y[4]) ^ (y[1] & y[2] & y[3] & y[4]) ^ (y[0] & y[1] & y[2]) ^ (y[0] & y[2] & y[3]) ^ (y[1] & y[2] & y[3]) ^ (y[0] & y[2] & y[4]) ^ (y[1] & y[2] & y[4]) ^ (y[2] & y[3] & y[4]) ^ (y[0] & y[1]) ^ (y[1] & y[2]) ^ (y[1] & y[3]) ^ (y[2] & y[3]) ^ y[2])) == 0]
        constraints += [((x[0] & x[1] & x[2] & x[3] & ~x[4]) & ~((y[0] & y[1] & y[3] & y[4]) ^ (y[0] & y[2] & y[3] & y[4]) ^ (y[0] & y[1] & y[3]) ^ (y[0] & y[2] & y[3]) ^ (y[0] & y[1] & y[4]) ^ (y[0] & y[2] & y[4]) ^ (y[1] & y[3] & y[4]) ^ (y[2] & y[3] & y[4]) ^ (y[0] & y[3]) ^ (y[1] & y[3]) ^ (y[2] & y[3]) ^ (y[3] & y[4]) ^ y[1] ^ y[2])) == 0]
        constraints += [((~x[0] & ~x[1] & ~x[2] & ~x[3] & x[4]) & ~((y[0] & y[1] & y[2] & y[3] & y[4]) ^ (y[0] & y[1] & y[2] & y[3]) ^ (y[0] & y[1] & y[2] & y[4]) ^ (y[1] & y[2] & y[3] & y[4]) ^ (y[0] & y[1] & y[2]) ^ (y[1] & y[2] & y[3]) ^ (y[1] & y[2] & y[4]) ^ (y[1] & y[2]) ^ (y[0] & y[3]) ^ (y[2] & y[3]) ^ (y[2] & y[4]) ^ (y[3] & y[4]))) == 0]
        constraints += [((x[0] & ~x[1] & ~x[2] & ~x[3] & x[4]) & ~((y[0] & y[1] & y[2] & y[3]) ^ (y[0] & y[1] & y[2] & y[4]) ^ (y[1] & y[2] & y[3] & y[4]) ^ (y[0] & y[1] & y[3]) ^ (y[0] & y[2] & y[3]) ^ (y[1] & y[2] & y[3]) ^ (y[0] & y[1] & y[4]) ^ (y[0] & y[2] & y[4]) ^ (y[1] & y[3] & y[4]) ^ (y[2] & y[3] & y[4]) ^ (y[0] & y[3]) ^ (y[0] & y[4]) ^ (y[1] & y[4]) ^ (y[2] & y[4]) ^ (y[3] & y[4]) ^ y[3])) == 0]
        constraints += [((~x[0] & x[1] & ~x[2] & ~x[3] & x[4]) & ~((y[0] & y[1] & y[2] & y[3] & y[4]) ^ (y[0] & y[1] & y[2] & y[3]) ^ (y[0] & y[1] & y[2] & y[4]) ^ (y[0] & y[1] & y[3] & y[4]) ^ (y[0] & y[2] & y[3] & y[4]) ^ (y[0] & y[1] & y[2]) ^ (y[0] & y[1] & y[3]) ^ (y[0] & y[2] & y[3]) ^ (y[0] & y[1] & y[4]) ^ (y[0] & y[2] & y[4]) ^ (y[0] & y[3] & y[4]) ^ (y[0] & y[1]) ^ (y[0] & y[2]) ^ (y[1] & y[3]) ^ (y[2] & y[3]) ^ (y[0] & y[4]) ^ (y[1] & y[4]) ^ (y[2] & y[4]) ^ (y[3] & y[4]) ^ y[0])) == 0]
        constraints += [((x[0] & x[1] & ~x[2] & ~x[3] & x[4]) & ~(y[2]  &  ((y[1] & y[3] & y[4]) ^ (y[1] & y[3]) ^ (y[1] & y[4]) ^ y[1] ^ y[3] ^ y[4]))) == 0]
        constraints += [((~x[0] & ~x[1] & x[2] & ~x[3] & x[4]) & ~((y[0] & y[1]) ^ (y[0] & y[2]) ^ (y[0] & y[3]) ^ (y[1] & y[3]) ^ (y[2] & y[3]) ^ (y[1] & y[4]) ^ (y[3] & y[4]))) == 0]
        constraints += [((x[0] & ~x[1] & x[2] & ~x[3] & x[4]) & ~((y[0] & y[1] & y[2] & y[3] & y[4]) ^ (y[0] & y[1] & y[2] & y[4]) ^ (y[0] & y[1] & y[3] & y[4]) ^ (y[0] & y[2] & y[3] & y[4]) ^ (y[0] & y[1] & y[4]) ^ (y[0] & y[2] & y[4]) ^ (y[0] & y[3] & y[4]) ^ (y[0] & y[1]) ^ (y[0] & y[2]) ^ (y[2] & y[3]) ^ (y[0] & y[4]))) == 0]
        constraints += [((~x[0] & x[1] & x[2] & ~x[3] & x[4]) & ~((y[0] & y[1] & y[2] & y[3] & y[4]) ^ (y[0] & y[1] & y[2] & y[4]) ^ (y[0] & y[1] & y[3] & y[4]) ^ (y[0] & y[2] & y[3] & y[4]) ^ (y[1] & y[2] & y[3] & y[4]) ^ (y[0] & y[1] & y[4]) ^ (y[0] & y[2] & y[4]) ^ (y[1] & y[2] & y[4]) ^ (y[0] & y[3] & y[4]) ^ (y[1] & y[3] & y[4]) ^ (y[2] & y[3] & y[4]) ^ (y[0] & y[1]) ^ (y[0] & y[2]) ^ (y[0] & y[3]) ^ (y[2] & y[3]) ^ (y[0] & y[4]) ^ (y[1] & y[4]) ^ (y[2] & y[4]) ^ y[4])) == 0]
        constraints += [((x[0] & x[1] & x[2] & ~x[3] & x[4]) & ~((y[0] & y[1] & y[2] & y[3] & y[4]) ^ (y[0] & y[1] & y[2] & y[3]) ^ (y[0] & y[1] & y[3] & y[4]) ^ (y[0] & y[2] & y[3] & y[4]) ^ (y[1] & y[2] & y[3] & y[4]) ^ (y[0] & y[1] & y[3]) ^ (y[0] & y[2] & y[3]) ^ (y[1] & y[2] & y[3]) ^ (y[0] & y[3] & y[4]) ^ (y[1] & y[3] & y[4]) ^ (y[2] & y[3] & y[4]) ^ (y[0] & y[1]) ^ (y[0] & y[2]) ^ (y[0] & y[3]) ^ (y[1] & y[4]) ^ (y[3] & y[4]) ^ y[3])) == 0]
        constraints += [((~x[0] & ~x[1] & ~x[2] & x[3] & x[4]) & ~((y[0] & y[1] & y[2] & y[3] & y[4]) ^ (y[0] & y[1] & y[2] & y[3]) ^ (y[0] & y[1] & y[2] & y[4]) ^ (y[0] & y[1] & y[2]) ^ (y[0] & y[3]) ^ (y[1] & y[3]) ^ (y[2] & y[3]) ^ (y[2] & y[4]) ^ (y[3] & y[4]))) == 0]
        constraints += [((x[0] & ~x[1] & ~x[2] & x[3] & x[4]) & ~((y[0] & y[1] & y[2] & y[3] & y[4]) ^ (y[0] & y[1] & y[2] & y[4]) ^ (y[0] & y[1] & y[3] & y[4]) ^ (y[0] & y[2] & y[3] & y[4]) ^ (y[0] & y[1] & y[4]) ^ (y[0] & y[2] & y[4]) ^ (y[0] & y[3] & y[4]) ^ (y[2] & y[3]) ^ (y[0] & y[4]) ^ (y[1] & y[4]) ^ (y[2] & y[4]))) == 0]
        constraints += [((~x[0] & x[1] & ~x[2] & x[3] & x[4]) & ~((y[0] & y[1] & y[2] & y[3] & y[4]) ^ (y[0] & y[1] & y[2] & y[3]) ^ (y[0] & y[1] & y[2] & y[4]) ^ (y[0] & y[1] & y[3] & y[4]) ^ (y[0] & y[2] & y[3] & y[4]) ^ (y[0] & y[1] & y[2]) ^ (y[0] & y[1] & y[3]) ^ (y[0] & y[2] & y[3]) ^ (y[0] & y[1] & y[4]) ^ (y[0] & y[2] & y[4]) ^ (y[0] & y[3] & y[4]) ^ (y[0] & y[1]) ^ (y[0] & y[2]) ^ (y[2] & y[3]) ^ (y[0] & y[4]) ^ (y[1] & y[4]) ^ (y[2] & y[4]) ^ (y[3] & y[4]) ^ y[0])) == 0]
        constraints += [((x[0] & x[1] & ~x[2] & x[3] & x[4]) & ~((y[0] & y[1] & y[2] & y[3] & y[4]) ^ (y[0] & y[1] & y[2] & y[3]) ^ (y[0] & y[1] & y[3] & y[4]) ^ (y[0] & y[2] & y[3] & y[4]) ^ (y[1] & y[2] & y[3] & y[4]) ^ (y[0] & y[1] & y[3]) ^ (y[0] & y[2] & y[3]) ^ (y[1] & y[2] & y[3]) ^ (y[0] & y[3] & y[4]) ^ (y[1] & y[3] & y[4]) ^ (y[2] & y[3] & y[4]) ^ (y[0] & y[3]) ^ (y[2] & y[4]) ^ (y[3] & y[4]) ^ y[3])) == 0]
        constraints += [((~x[0] & ~x[1] & x[2] & x[3] & x[4]) & ~((y[0] & y[1] & y[2] & y[3] & y[4]) ^ (y[0] & y[1] & y[2] & y[3]) ^ (y[0] & y[1] & y[2] & y[4]) ^ (y[0] & y[1] & y[3] & y[4]) ^ (y[1] & y[2] & y[3] & y[4]) ^ (y[0] & y[1] & y[2]) ^ (y[0] & y[1] & y[3]) ^ (y[1] & y[2] & y[3]) ^ (y[0] & y[1] & y[4]) ^ (y[1] & y[2] & y[4]) ^ (y[1] & y[3] & y[4]) ^ (y[0] & y[2]) ^ (y[1] & y[2]) ^ (y[0] & y[3]) ^ (y[1] & y[3]) ^ (y[2] & y[3]) ^ (y[3] & y[4]) ^ y[1])) == 0]
        constraints += [((x[0] & ~x[1] & x[2] & x[3] & x[4]) & ~((y[0] & y[1] & y[2] & y[3]) ^ (y[0] & y[1] & y[2] & y[4]) ^ (y[1] & y[2] & y[3] & y[4]) ^ (y[0] & y[1] & y[3]) ^ (y[0] & y[2] & y[3]) ^ (y[1] & y[2] & y[3]) ^ (y[0] & y[1] & y[4]) ^ (y[0] & y[2] & y[4]) ^ (y[1] & y[3] & y[4]) ^ (y[2] & y[3] & y[4]) ^ (y[0] & y[1]) ^ (y[0] & y[2]) ^ (y[0] & y[3]) ^ (y[0] & y[4]) ^ (y[3] & y[4]) ^ y[3])) == 0]
        constraints += [((~x[0] & x[1] & x[2] & x[3] & x[4]) & ~((y[0] & y[1] & y[2] & y[3] & y[4]) ^ (y[0] & y[1] & y[2] & y[4]) ^ (y[0] & y[1] & y[3] & y[4]) ^ (y[0] & y[2] & y[3] & y[4]) ^ (y[1] & y[2] & y[3] & y[4]) ^ (y[0] & y[1] & y[4]) ^ (y[0] & y[2] & y[4]) ^ (y[1] & y[2] & y[4]) ^ (y[0] & y[3] & y[4]) ^ (y[1] & y[3] & y[4]) ^ (y[2] & y[3] & y[4]) ^ (y[0] & y[1]) ^ (y[0] & y[2]) ^ (y[0] & y[3]) ^ (y[1] & y[3]) ^ (y[2] & y[3]) ^ (y[0] & y[4]) ^ (y[1] & y[4]) ^ (y[2] & y[4]) ^ y[4])) == 0]
        constraints += [((x[0] & x[1] & x[2] & x[3] & x[4]) & ~((y[0] & y[1] & y[2] & y[3] & y[4]) ^ (y[0] & y[1] & y[2] & y[3]) ^ (y[0] & y[1] & y[2] & y[4]) ^ (y[0] & y[1] & y[3] & y[4]) ^ (y[1] & y[2] & y[3] & y[4]) ^ (y[0] & y[1] & y[2]) ^ (y[0] & y[1] & y[3]) ^ (y[1] & y[2] & y[3]) ^ (y[0] & y[1] & y[4]) ^ (y[1] & y[2] & y[4]) ^ (y[1] & y[3] & y[4]) ^ (y[0] & y[2]) ^ (y[1] & y[2]) ^ (y[1] & y[3]) ^ (y[2] & y[3]) ^ y[1])) == 0]

        return constraints
