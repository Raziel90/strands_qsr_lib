# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 17:40:39 2015

@author: ccoppola
"""


from __future__ import print_function, division
from qsrlib_qsrs.qsr_qtc_simplified_abstractclass import QSR_QTC3D_Simplified_Abstractclass


class qtc3DException(Exception):
    pass

class QSR_QTC3D_B1_Simplified(QSR_QTC3D_Simplified_Abstractclass):
    def __init__(self):
        super(QSR_QTC3D_B1_Simplified, self).__init__()
        self._unique_id = "QTC3D"
        self.qtc3Dtype= "b1"
        self.all_possible_relations = tuple(self.return_all_possible_state_combinations)
        self._dtype = "points"

    def _compute_qsr(self, data1, data2, qsr_params, **kwargs):
        return {
            data1.x < data2.x: "left",
            data1.x > data2.x: "right"
        }.get(True, "together")
        
    def return_all_possible_state_combinations(self):
        """Method that returns all possible state combinations for the qtc_type b1
        defined for this calss instance.

        :return:
            - String representation as a list of possible tuples
            - Integer representation as a list of lists of possible tuples
        """
        ret_str = []
        ret_int = []
        for i in xrange(1, 4):
            for j in xrange(1, 4):
                for k in xrange(1, 4):
                    for l in xrange(1, 4):
                        for n in xrange(1, 4):
                            ret_int.append([i-2, j-2, k-2, l-2, n-2])
                            ret_str.append(str(i-2) + "," + str(j-2) + "," + str(k-2) + "," + str(l-2) + "," + str(n-2))
        
        return [s.replace('-1','-').replace('1','+') for s in ret_str], ret_int