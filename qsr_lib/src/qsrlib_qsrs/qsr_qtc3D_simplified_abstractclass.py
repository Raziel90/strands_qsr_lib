# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 15:32:45 2015

@author: ccoppola
"""

# -*- coding: utf-8 -*-
from abc import abstractmethod, ABCMeta
from __future__ import print_function, division
from qsrlib_qsrs.qsr_dyadic_abstractclass import QSR_Dyadic_1t_Abstractclass
import numpy as np

class qtc3DException(Exception):
    pass

class QSR_QTC3D_Simplified_Abstractclass(QSR_Dyadic_1t_Abstractclass):
    def __init__(self):
        super(QSR_QTC3D_Simplified_Abstractclass, self).__init__()
        self._unique_id = "QTC3D"
        self._dtype = "points"
        
        self.__qsr_params_defaults= {
            "distance_threshold": 0.01,
            "angle_threshold": 50
        }
        
    def _process_qsr_parameters_from_request_parameters(self, req_params, **kwargs):
        "THIS IS COPIED FROM CHRISTIAN'S QTC TO ADAPT!!!!"        
        
        qsr_params = self.__qsr_params_defaults.copy()

        try: # global namespace
            if req_params["dynamic_args"]["for_all_qsrs"]:
                for k, v in req_params["dynamic_args"]["for_all_qsrs"].items():
                    qsr_params[k] = v
        except KeyError:
            pass

        try: # General case
            if req_params["dynamic_args"][self.__global_unique_id]:
                for k, v in req_params["dynamic_args"][self.__global_unique_id].items():
                    qsr_params[k] = v
        except KeyError:
            pass

        try: # Parameters for a specific variant
            if req_params["dynamic_args"][self._unique_id]:
                for k, v in req_params["dynamic_args"][self._unique_id].items():
                    qsr_params[k] = v
        except KeyError:
            pass

        if not isinstance(qsr_params["no_collapse"], bool) or not isinstance(qsr_params["validate"], bool):
            raise TypeError("'no_collapse' and 'validate' have to be boolean values.")

        for param in qsr_params:
            if param not in self.__qsr_params_defaults and param not in self._allowed_parameters:
                raise KeyError("%s is an unknown parameter" % str(param))

        return qsr_params        
        
    def _create_qtc_representation(self, pos_k, pos_l, quantisation_factor=0):
        """Creating the QTCC representation for the given data. Uses the
        double cross to determine to which side of the lines the points are
        moving.

        :param pos_k: An array of positions for agent k, exactly 2 entries of x,y positions
        :param pos_l: An array of positions for agent l, exactly 2 entries of x,y positions
        :param quantisation_factor: The minimum distance the points have to diverge from either line to be regrded a non-0-state

        :return: The QTCC 4-tuple (q1,q2,q4,q5) for the movement of the two agents: [k[0],l[0],k[1],l[1]]

        """
        #print "######################################################"
        pos_k = np.array(pos_k).reshape(-1, 3)
        pos_l = np.array(pos_l).reshape(-1, 3)
        
        quantisation_factor=0
        distKprevLnow=np.linalg.norm(pos_k[0]-pos_l[1])
        distKnowLnow=np.linalg.norm(pos_k[1]-pos_l[1])
        distKnextLnow=np.linalg.norm(pos_k[2]-pos_l[1])
        if (distKprevLnow > distKnowLnow + quantisation_factor) & (distKnowLnow > distKnextLnow + quantisation_factor):
            #-
        elif (distKprevLnow < distKnowLnow + quantisation_factor) & (distKnowLnow < distKnextLnow + quantisation_factor)
            #+
        else #0
         
        distKnowLprev=np.linalg.norm(pos_k[1]-pos_l[0])
        distKnowtLnext=np.linalg.norm(pos_k[1]-pos_l[2])
        if (distKnowLprev > distKnowLnow + quantisation_factor) & (distKnowLnow > distKnowtLnext + quantisation_factor):
            #-
        elif (distKnowLprev < distKnowLnow + quantisation_factor) & (distKnowLnow < distKnowtLnext + quantisation_factor)
            #+
        else #0
        
        
        
        
        tk=(pos_k[2]-pos_k[1])/np.linalg.norm(pos_k[2]-pos_k[1])
        tl=(pos_l[2]-pos_l[1])/np.linalg.norm(pos_l[2]-pos_l[1])
        tkprev=(pos_k[1]-pos_k[0])/np.linalg.norm(pos_k[1]-pos_k[0])
        tlprev=(pos_l[1]-pos_l[0])/np.linalg.norm(pos_l[1]-pos_l[0])        
        
        
        
        # Creating double cross, RL_ext being the connecting line, trans_RL_k
        # and l being the orthogonal lines going through k and l respectively.
        RL_ext = np.append(
            self._translate(pos_k[-2], (pos_k[-2]-pos_l[-2])/2),
            self._translate(pos_l[-2], (pos_l[-2]-pos_k[-2])/2)
        ).reshape(-1,2)
        #print "RL_ext", RL_ext
        rot_RL = self._orthogonal_line(
            pos_k[-2],
            np.append(pos_k[-2], (pos_l[-2]-pos_k[-2]))
        ).reshape(-1,2)
        #print "rot_RL", rot_RL
        trans_RL_k = self._translate(
            [rot_RL[0], rot_RL[1]],
            (rot_RL[0]-rot_RL[1])/2
        )
        #print "transk", trans_RL_k
        trans_RL_l = self._translate(
            trans_RL_k[0:2],
            (pos_l[-2]-pos_k[-2])
        )
        #print "transl", trans_RL_l

        # Test constraints for k
        k = np.append(
            self._test_constraint(
                pos_k,
                trans_RL_k,
                quantisation_factor=quantisation_factor),
            self._test_constraint(
                pos_k,
                RL_ext,
                quantisation_factor=quantisation_factor,
                constraint="side")
        )
        #print "k", k

        # Test constraints for l
        l = np.append(
            self._test_constraint(
                pos_l,
                np.array([ # Needs to be turned around to determine correct side
                    [trans_RL_l[1,0],trans_RL_l[1,1]],
                    [trans_RL_l[0,0],trans_RL_l[0,1]]
                ]),
                quantisation_factor=quantisation_factor),
            self._test_constraint(
                pos_l,
                np.array([ # Needs to be turned around to determine correct side
                    [RL_ext[1,0],RL_ext[1,1]],
                    [RL_ext[0,0],RL_ext[0,1]]
                ]),
                quantisation_factor=quantisation_factor,
                constraint="side")
        )
        #print "l", l

        return np.array([k[0],l[0],k[1],l[1]])
        
    @abstractmethod          
    def _compute_qsr(self, data1, data2, qsr_params, **kwargs):
        return
    
    @abstractmethod    
    def return_all_possible_state_combinations(self):
        """Method that returns all possible state combinations for the qtc_type
        defined for this class instance.

        :return:
            - String representation as a list of possible tuples
            - Integer representation as a list of lists of possible tuples
        """
        return