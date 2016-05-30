# -*- coding: utf-8 -*-
"""
Created on Thu May 26 15:34:51 2016

@author: claudio
"""
from qsrrep_hmms.qtc_hmm_abstractclass import QTCHMMAbstractclass
from qtc3d_converter.point_to_QTC3D import ind2qtc3d,qtc3d2ind
import pickle
import json
import numpy as np


class QTC3DHMM(QTCHMMAbstractclass):
    
    
    def __init__(self):
        super(QTC3DHMM, self).__init__()
        
                
    def _save_to_file(self,filename):
        with open(filename,'w') as f:
                pickle.dump(self.f)
                f.close()
    def _create_transition_matrix(self, size, **kwargs):
        """Creates a Conditional Neighbourhood Diagram for QTCB as a basis for the HMM.

        :param kwargs:
            * input_to_state: dict mapping the input symbol to the state of the hmm. Different for all 3 qtc versions.

        :return: The transition matrix only allowing transitions according to the CND
        """
        input_to_state = dict(kwargs["input_to_state"])
        num_states=max(input_to_state.values())
        tran=np.zeros((num_states,num_states))

        # Calling parent to generate actual matrix
        return super(QTC3DHMM, self)._create_transition_matrix(size=size, qtc=qtc)

    def _symbol_to_qsr(self, symbols):
        """Transforming alphabet symbols to QTCB states.

        :param symbols: A list of symbols

        :return: The list of corresponding qtc symbols
        """

        ret = []
        for s in symbols:
            

            ret.append(ind2qtc3d(s))

        return ret

    def symbol_to_qsr(self, symbol):
        print "TODO"

        return ind2qtc3d(symbol)
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
if __name__ == '__main__':
    a=QTC3DHMM()
    print a
    