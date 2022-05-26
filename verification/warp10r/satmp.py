"""
Applying monomial prediction technique to find integral
distinguishers of WARP in the single-key setting
Copyright (C) Jan 2, 2022 Hosein Hadipour

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import time
from itertools import combinations
from math import comb
from pysat import solvers
from pysat import formula
import random

"""
x_roundNumber_nibbleNumber_bitNumber
x_roundNumber_nibbleNumber_0: msb
x_roundNumber_nibbleNumber_3: lsb
Variable mapping:

x_r_0                       k_r  x_r_1  ...
|                            |     |
|--- y_r---> | S | ---z_r----+---->+    ...
|                                  |    ...

"""

class Warp:
    count = 0
    def __init__(self, nrounds=10, solver_name="cadical"):
        Warp.count += 1
        self.nrounds = nrounds
        self.sat_solver_name = solver_name
        self.supported_sat_solvers = list(solvers.SolverNames.cadical) + \
             list(solvers.SolverNames.glucose4) + \
                 list(solvers.SolverNames.glucose3) + \
                     list(solvers.SolverNames.lingeling) + \
                         list(solvers.SolverNames.maplesat) + \
                         list(solvers.SolverNames.maplechrono) + \
                             list(solvers.SolverNames.maplecm) + \
                                 list(solvers.SolverNames.minicard) + \
                                     list(solvers.SolverNames.minisat22) + \
                                         list(solvers.SolverNames.minisatgh)
        assert(self.sat_solver_name in self.supported_sat_solvers)
        if self.sat_solver_name in solvers.SolverNames.cadical:
            self.sat_solver = solvers.Cadical()
        elif self.sat_solver_name in solvers.SolverNames.glucose4:
            self.sat_solver = solvers.Glucose4()
        elif self.sat_solver_name in solvers.SolverNames.glucose3:
            self.sat_solver = solvers.Glucose3()
        elif self.sat_solver_name in solvers.SolverNames.lingeling:
            self.sat_solver = solvers.Lingeling()
        elif self.sat_solver_name in solvers.SolverNames.maplesat:
            self.sat_solver = solvers.Maplesat()
        elif self.sat_solver_name in solvers.SolverNames.maplechrono:
            self.sat_solver = solvers.MapleChrono()
        elif self.sat_solver_name in solvers.SolverNames.maplecm:
            self.sat_solver = solvers.MapleCM()
        elif self.sat_solver_name in solvers.SolverNames.minicard:
            self.sat_solver = solvers.Minicard()
        elif self.sat_solver_name in solvers.SolverNames.minisat22:
            self.sat_solver = solvers.Minisat22()
        elif self.sat_solver_name in solvers.SolverNames.minisatgh:
            self.sat_solver = solvers.MinisatGH()

        self.cnf_formula = formula.CNF()
        self.solver = self.solver = solvers.Solver(name=self.sat_solver_name)
        self.variables_dictionary = dict()
        self.top_variable_identifier_so_far = 0
        self.cnf_file_name = f"warp_nr_{self.nrounds}.cnf"
        self.result_file_name = f"result_nr_{self.nrounds}.txt"
        self.permute_nibbles = [31, 6, 29, 14, 1, 12, 21, 8, 27, 2, 3, 0, 25, 4, 23, 10,
                                15, 22, 13, 30, 17, 28, 5, 24, 11, 18, 19, 16, 9, 20, 7, 26]
        self.RC0 = [0x0, 0x0, 0x1, 0x3, 0x7, 0xf, 0xf, 0xf, 0xe, 0xd, 0xa, 0x5, 0xa, 0x5, 0xb, 0x6, 0xc, 0x9, 0x3, 0x6, 0xd, 0xb, 0x7, 0xe, 0xd, 0xb, 0x6, 0xd, 0xa, 0x4, 0x9, 0x2, 0x4, 0x9, 0x3, 0x7, 0xe, 0xc, 0x8, 0x1, 0x2]
        self.RC1 = [0x4, 0xc, 0xc, 0xc, 0xc, 0xc, 0x8, 0x4, 0x8, 0x4, 0x8, 0x4, 0xc, 0x8, 0x0, 0x4, 0xc, 0x8, 0x4, 0xc, 0xc, 0x8, 0x4, 0xc, 0x8, 0x4, 0x8, 0x0, 0x4, 0x8, 0x0, 0x4, 0xc, 0xc, 0x8, 0x0, 0x0, 0x4, 0x8, 0x4, 0xc]                                
        # 1: msb of input
        # 5: msb of output
        # input bits: 1, 2, 3, 4
        # output bits: 5, 6, 7, 8
        # convert monomial prediction table (MPT) to a boolean function
        # and then minimize its CNF (POS) representation via Quine-McCluskey algorithm
        self.sbox_cnf_pattern = [[-1, -4, 5, -6, -8],
                                [-1, 2, 4, 5, 6],
                                [2, 3, -7, -8],
                                [3, -6, -8],
                                [3, 4, -8],
                                [-1, -3, -4, -5, 6, -8],
                                [1, 2, -4, 5, 6, 7],
                                [2, -6, -7],
                                [2, 4, -7],
                                [1, -3, 4, 8],
                                [-3, -5, -6, 8],
                                [-1, -2, -4, 5, 6, 8],
                                [-2, -5, -6, 7],
                                [1, -2, 4, 7],
                                [-3, 5, 6, 8],
                                [-1, -2, -4, 7],
                                [-1, -2, -3, -7, 8],
                                [-2, 3, 5, 7, 8],
                                [-2, -3, -4, 6, -7],
                                [-1, 2, 3, 6, 7, 8],
                                [-2, -3, -4, 6, 8],
                                [3, -4, 6, 7, 8],
                                [1, 2, -3, -6, 8],
                                [-2, 4, -5, 7, -8],
                                [-4, 5, -6, -7, -8],
                                [2, -5, -7, -8],
                                [-2, 4, -6, 7, -8],
                                [1, 2, 3, -8],
                                [2, -3, 4, -6, 8]]
        # xor inputs: 1, 2
        # xor output: 3
        self.xor_cnf_template = [[-1, -2], [-2, 3], [-1, 3], [1, 2, -3]]

        # xor3 inputs: 1, 2, 3
        # xor3 output: 4
        self.xor3_cnf_template = [[-2, -3], [-1, -3], [-1, -2], [-3, 4], [-2, 4], [-1, 4], [1, 2, 3, -4]]
        self.xor3_with_constant_cnf_template = [[-2, -3], [-1, -3], [-1, -2], [-3, 4], [-2, 4], [-1, 4]]

    def generate_sbox_constraint(self, input_bits, output_bits):
        substitution_list = [self.variables_dictionary[x] for x in input_bits + output_bits]
        for sl in self.sbox_cnf_pattern:
            temp = []
            for index in sl:
                if index > 0:
                    temp.append(substitution_list[index - 1])
                else:
                    temp.append(-substitution_list[abs(index) - 1])
            self.cnf_formula.append(temp)

    def update_variables_dictionary(self, new_vars):
        """
        This method is used to update variables' dictionary
        """
        for nv in new_vars:
            if nv not in self.variables_dictionary.keys():
                self.top_variable_identifier_so_far += 1
                self.variables_dictionary[nv] = self.top_variable_identifier_so_far

    @staticmethod
    def flatten_state(s):
        state_bits = [s[i][j] for i in range(len(s)) for j in range(len(s[0]))]
        return state_bits

    def inv_permute_nibbles(self, state):
        temp = [0]*32
        for i in range(32):
            temp[i] = state[self.permute_nibbles[i]]
        return temp    
    
    def generate_round_x_variables(self, rn):
        """
        Generate the input variables of rn'th round
        """

        x = [[f"x_{rn}_{nibble}_{bit}" for bit in range(4)] for nibble in range(32)]
        self.update_variables_dictionary(self.flatten_state(x))
        return x
    
    def generate_round_y_z_k_variables(self, rn):
        """
        Generate the intermediate variables in rn'th round
        """

        y = [[f"y_{rn}_{nibble}_{bit}" for bit in range(4)] for nibble in range(16)]
        z = [[f"z_{rn}_{nibble}_{bit}" for bit in range(4)] for nibble in range(16)]
        k = [[f"k_{rn}_{nibble}, {bit}" for bit in range(4)] for nibble in range(16)]
        self.update_variables_dictionary(self.flatten_state(y))
        self.update_variables_dictionary(self.flatten_state(z))
        self.update_variables_dictionary(self.flatten_state(k))
        return y, z, k

    def constraints_by_fork(self, a, b1, b2):
        """
        a ---fork---> (b1, b2)
        """

        self.cnf_formula.append([-self.variables_dictionary[a],\
                                self.variables_dictionary[b1],\
                                self.variables_dictionary[b2]])
        self.cnf_formula.append([-self.variables_dictionary[b1], self.variables_dictionary[a]])
        self.cnf_formula.append([-self.variables_dictionary[b2], self.variables_dictionary[a]])
    
    def nibble_fork(self, s, s1, s2):
        constraints = ""
        for bit in range(4):
            self.constraints_by_fork(s[bit], s1[bit], s2[bit])
        return constraints

    def constraints_by_xor(self, a1, a2, b):
        """
        a1, a2 ---> b = a1 + a2
        """

        substitution_list = [self.variables_dictionary[x] for x in  [a1, a2, b]]
        for sl in self.xor_cnf_template:
            temp = []
            for index in sl:
                if index > 0:
                    temp.append(substitution_list[index - 1])
                else:
                    temp.append(-substitution_list[abs(index) - 1])
            self.cnf_formula.append(temp)
    
    def nibble_xor(self, s1, s2, s):
        for bit in range(4):
            self.constraints_by_xor(s1[bit], s2[bit], s[bit])
    
    def constraints_by_3xor(self, a1, a2, a3, b, constant=0):
        """
        a1, a2, a3 ---> b = a1 + a2 + a3
        or
        a1, a2, a3 ---> b = a1 + a2 + a3 + 1
        """
        
        assert(constant in [0, 1])
        substitution_list = [self.variables_dictionary[x] for x in  [a1, a2, a3, b]]        
        if constant == 0:
            template = self.xor3_cnf_template
        elif constant == 1:
            template = self.xor3_with_constant_cnf_template        
        for sl in template:
            temp = []
            for index in sl:
                if index > 0:
                    temp.append(substitution_list[index - 1])
                else:
                    temp.append(-substitution_list[abs(index) - 1])
            self.cnf_formula.append(temp)
    
    def nibble_3xor(self, s1, s2, s3, s):
        for bit in range(4):
            self.constraints_by_3xor(s1[bit], s2[bit], s3[bit], s[bit])
    
    def generate_constraints(self):
        """
        Generate the constraints of MILP model
        """

        for rn in range(self.nrounds):
            x_in = self.generate_round_x_variables(rn)
            y, z, k = self.generate_round_y_z_k_variables(rn)
            x_out = self.generate_round_x_variables(rn + 1)
            x_middle = self.inv_permute_nibbles(x_out)
            for nibble in range(16):                
                self.nibble_fork(x_in[2*nibble], y[nibble], x_middle[2*nibble])
                self.generate_sbox_constraint(y[nibble], z[nibble])
                self.nibble_3xor(z[nibble], k[nibble], x_in[2*nibble + 1], x_middle[2*nibble + 1])
    
    def exclude_key_independent_term(self):
        """
        Limit the key variables to not be all zero
        """

        key_vars = [self.variables_dictionary[var] \
            for var in self.variables_dictionary.keys() if "k" in var]
        self.cnf_formula.append(key_vars)
    
    def generate_sat_model(self):
        self.generate_constraints()
        # self.exclude_key_independent_term()
        # self.cnf_formula.to_file(self.cnf_file_name)
        self.sat_solver.append_formula(self.cnf_formula)

    def check_cube(self, active_indices=list(range(127)), target_output_bits = range(128)):
        # Fix the input choice vector
        input_vars = self.flatten_state(self.generate_round_x_variables(0))
        output_vars = self.flatten_state(self.generate_round_x_variables(self.nrounds))

        input_active_pattern = []
        for i in range(128):
            if i in active_indices:
                input_active_pattern.append(self.variables_dictionary[input_vars[i]])
            else:
                # input_active_pattern.append(-self.variables_dictionary[input_vars[i]])
                pass
        
        balanced_bits = []
        not_checked_bits = []
        for output_bit in target_output_bits:
            output_active_pattern = []
            for i in range(128):
                if i != output_bit:
                    output_active_pattern.append(-self.variables_dictionary[output_vars[i]])
                else:
                    output_active_pattern.append(self.variables_dictionary[output_vars[i]])
            assumptions = input_active_pattern + output_active_pattern
            ##########################
            ##########################
            result = self.sat_solver.solve(assumptions=assumptions)
            ##########################
            ##########################
            if result == True:
                # print("Output bit number {:03d} may NOT be key-independent :-(".format(output_bit))
                pass
            elif result == False:
                balanced_bits.append(output_bit)
                # print("Output bit number {:03d} is key-independent ;-)".format(output_bit))
            else:
                not_checked_bits.append(output_bit)
                # print("Output bit number {:03d} was not checked!".format(output_bit))
        return balanced_bits
    
if __name__ == '__main__':
    NUMBER_OF_ROUNDS = 10
    BRUTE_FORCE = True
    separator_line = "#"*100 + "\n"

    warp = Warp(nrounds=NUMBER_OF_ROUNDS,
                solver_name="minisat22")
    with open(warp.result_file_name, "w") as outputfile:
        outputfile.write(f"Results of applying monomial prediction method on {warp.nrounds} rounds of WARP\n")
    mp_dict = dict()
    warp.generate_sat_model()
    
    number_of_active_nibbles = 1
    distinguisher_flag = False

    start_time = time.time()
    if BRUTE_FORCE == False:
        while not distinguisher_flag:
            num_of_possible_cases = comb(32, number_of_active_nibbles)
            num_of_tries = int(num_of_possible_cases**0.5)
            print(f"Looking for a distinguisher with {number_of_active_nibbles} active nibbles ...")
            for tr in range(num_of_tries):
                active_nibbles = random.sample(range(32), number_of_active_nibbles)
                active_bits = []
                for nibble in active_nibbles:
                    active_bits.extend(list(range(4*nibble, 4*(nibble + 1))))
                balanced_bits = warp.check_cube(active_indices=active_bits)
                if balanced_bits != []:
                    distinguisher_flag = True            
                    break
            number_of_active_nibbles += 1
    else:
        while not distinguisher_flag:
            print(f"Looking for a distinguisher with {number_of_active_nibbles} active nibbles ...")
            for active_nibbles in combinations(list(range(32)), number_of_active_nibbles):                  
                active_bits = []
                for nibble in active_nibbles:
                    active_bits.extend(list(range(4*nibble, 4*(nibble + 1))))       
                balanced_bits = warp.check_cube(active_indices=active_bits)
                if balanced_bits != []:
                    distinguisher_flag = True            
                    break
            number_of_active_nibbles += 1
    elapsed_time = time.time() - start_time

    if distinguisher_flag == True:
        number_of_active_nibbles -= 1
        print("Integral distinguisher found :-)")
        time_line = "Total time to find an intergal distinguisher: {:0.02f}\n".format(elapsed_time)
        print(time_line)
        print(f"Looking for the best integral distinguisher with {len(active_nibbles)} active nibbles ...")
        best_balanced_bits = balanced_bits
        best_active_nibbles = active_nibbles
        for active_nibbles in combinations(list(range(32)), number_of_active_nibbles):
            active_bits = []
            for nibble in active_nibbles:
                active_bits.extend(list(range(4*nibble, 4*(nibble + 1))))  
            balanced_bits = warp.check_cube(active_indices=active_bits)
            if len(balanced_bits) > len(best_balanced_bits):
                best_balanced_bits = balanced_bits
                best_active_nibbles = active_nibbles        
        with open(warp.result_file_name, "a") as outputfile:
            na = len(best_active_nibbles)        
            active_nibbles = ", ".join(map(str, best_active_nibbles))
            nb = len(best_balanced_bits)
            balanced_bits = ", ".join(map(str, best_balanced_bits))
            outputfile.write(f"int NUMBER_OF_ROUNDS = {NUMBER_OF_ROUNDS};\n")
            outputfile.write(f"int na = {na};\n")
            outputfile.write(f"int nb = {nb};\n")
            outputfile.write("int active_nibbles[{}] = {{{}}};\n".format(na, active_nibbles))
            outputfile.write("int balanced_positions[{}] = {{{}}};\n".format(nb, balanced_bits))
    else:
        print("Try with more active nibbles :-(")
    
