# -*- coding: utf-8 -*-
"""
Created on Thu May 26 15:34:51 2016

@author: claudio
"""
from qsrrep_hmms.qtc_hmm_abstractclass import QTCHMMAbstractclass
from qsrrep_hmms.hmm_abstractclass import HMMAbstractclass
from qtc3d_converter.point_to_QTC3D import ind2qtc3d,qtc3d2ind
import pickle
import ghmm as gh
import numpy as np
NUM_SYMBOLS=3**5+2*3**2

class QTC3DHMM(QTCHMMAbstractclass):
    
    
    def __init__(self):
        super(QTC3DHMM, self).__init__()
        self.num_possible_states = 0
        
                
    def _save_to_file(self,filename):
        with open(filename,'w') as f:
                pickle.dump(self.f)
                f.close()
                
    def _create_emission_matrix(self, size, **kwargs):
        """Creates a Conditional Neighbourhood Diagram for QTCB as a basis for the HMM.

        :param kwargs:
            * size: unused
            * input_to_state: dict mapping the input symbol to the state of the hmm. Different for all 3 qtc versions.

        :return: The transition matrix only allowing transitions according to the CND
        """
        input_to_state = dict(kwargs["input_to_state"])
        sequences_of_activity = kwargs["qsr_seq"]
        num_states=max(input_to_state.values())
        num_emiss=NUM_SYMBOLS
        emiss=np.ones((num_states+2,num_emiss))
        for sequence in sequences_of_activity:
            for sample in sequence:
                emiss[input_to_state[sample],sample]+=1
        return emiss/emiss.sum(axis=1).reshape(-1, 1)
        
        
    def _create_transition_matrix(self, size, **kwargs):
        """Creates a Conditional Neighbourhood Diagram for QTCB as a basis for the HMM.

        :param kwargs:
            * size: unused
            * input_to_state: dict mapping the input symbol to the state of the hmm. Different for all 3 qtc versions.
            * sequences_of_activity: the sequences of QSRs. This should be a list of state chains, i.e. a list of lists
            
        :return: The transition matrix only allowing transitions according to the CND
        """
        input_to_state = dict(kwargs["input_to_state"])
        sequences_of_activity = kwargs["qsr_seq"]
        
        
        
        num_states=max(input_to_state.values())
        tran=np.ones((num_states+2,num_states+2))
        
        for sequence in sequences_of_activity:
            tran[0,input_to_state[sequence[0]]]+=1
            tran[input_to_state[sequence[-1]],-1]+=1
            for ind,sample in enumerate(sequence[:-1]):
                #print sample,sequence[ind+1]
                tran[input_to_state[sample],input_to_state[sequence[ind+1]]]+=1
        
        
        
        tran[tran != 1] = 0
        tran[tran == 0] = 0.00001
        tran[0] = 1
        tran[:, 0] = 0
        tran[:, -1] = 1
        tran[0, -1] = 0
        tran[-1] = 0
        tran += np.dot(np.eye(len(tran)), 0.00001)
        tran[0, 0] = 0    
        
        tran=tran/tran.sum(axis=1).reshape(-1, 1)
        # Calling parent to generate actual matrix
        return tran
    def _create(self, **kwargs):
        """Creates and trains (using '_train') a HMM to represent the given qtc sequences.
        Main function to create and train the hmm. Please override with special
        behaviour if necessary.

        This function is called by the library to create the hmm.

        :param **kwargs:
            * qsr_seq: the sequences of QSRs. This should be a list of state chains, i.e. a list of lists
            * input_to_state: dict mapping the input symbol to the state of the hmm. Different for all 3 qtc versions.

        :return: The trained HMM

        """
        input_to_state = kwargs["input_to_state"]
        state_seq = kwargs["qsr_seq"]
        
        num_states=max(input_to_state.values())
        #state_seq = self._qsr_to_symbol(kwargs["qsr_seq"])
        trans = self._create_transition_matrix(size=0, **kwargs)
        emi = self._create_emission_matrix(size=0, **kwargs)
        print len(trans),num_states+2
        hmm = self._train(state_seq, trans, emi, num_states+2)
        print '...done'
        return hmm        
        


    def _train(self, seq, trans, emi, num_possible_states):
        """Uses the given parameters to train a multinominal HMM to represent
        the given seqences of observations. Uses Baum-Welch training.
        Please override if special training is necessary for your QSR.

        :param seq: the sequence of observations represented by alphabet symbols
        :param trans: the transition matrix as a numpy array
        :param emi: the emission matrix as a numpy array
        :param num_possible_states: the total number of possible states

        :return: the via baum-welch training generated hmm
        """

        print 'Generating HMM:'
        print '\tCreating symbols...'
        symbols = self.generate_alphabet(NUM_SYMBOLS)
        startprob = np.zeros(num_possible_states)
        startprob[0] = 1
        print '\t\t', symbols
        print '\tCreating HMM...'
        print emi[:,0]
        hmm = gh.HMMFromMatrices(
            symbols,
            gh.DiscreteDistribution(symbols),
            trans.tolist(),
            emi.tolist(),
            startprob.tolist()
        )
        print '\tTraining...'
        hmm.baumWelch(self._create_sequence_set(seq, symbols))

        return hmm
        
        
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
        

        return ind2qtc3d(symbol)
        
        
    def _log_likelihood(self, **kwargs):
        """Computeed the loglikelihood for the given sample(s) to be produced by
        the HMM.

        :param kwargs:
            * qsr_seq: A list of lists of qsr sequences to check against the HMM
            * hmm: The to generate the loglikelihood for

        :return: The accumulated loglikelihood for all the given samples
        """

        return kwargs["hmm"].loglikelihood(self._create_sequence_set(
            qsr_seq=self._qsr_to_symbol(kwargs["qsr_seq"]),
            symbols=self.generate_alphabet(num_symbols=num_symbols)
        ))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
if __name__ == '__main__':
    HMM=QTC3DHMM()
    print HMM
    