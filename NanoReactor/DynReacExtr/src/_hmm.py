import json
import operator

import numpy as np
import hmmlearn.hmm as hmm

from tqdm import tqdm

from _setting_share import SharedSetting

class FilterNoise(SharedSetting):
    
    def __init__(self):
        SharedSetting.__init__(self)
    
    def hmmViterbi(self, obv):
        """
        Default set for hmm
        """
        states = ['break', 'form'] ; n_states = len(states)
        observations = ['break', 'form'] ; n_observations = len(observations)

        start_prob = np.array([0.5, 0.5])
        emiss_prob = np.array(self.obvmatr).reshape(2, 2)
        trans_prob = np.array(self.hidmatr).reshape(2, 2)

        model = hmm.MultinomialHMM(n_components=n_states)
        model.startprob_ = start_prob
        model.emissionprob_ = emiss_prob
        model.transmat_ = trans_prob

        logprob, box = model.decode(obv, algorithm='viterbi')
        
        return box

    def filterNoise(self, smi_node_t_list, step):
        """
        Transform smi_node's time list to observation list.
        Decode observation list using Viterbi algorithm.
        """
        obList = np.zeros((step), dtype=int)
        for t in smi_node_t_list:
            obList[t] = int(1)
        obList_filter = self.hmmViterbi(obList[np.newaxis, :].T)
        return obList_filter
    
    def getFilteredSig(self, smi_nodes, step):
        """
        Filter node dectected signal.

        :param smi_nodes -> dict
            key is str_smi_node_list and value is detetcted_time_list
        :param step -> int
            number of iterative readings of trajectory
        :return r_list_t_node, p_list_t_node -> list
            Returns reactant's and product's list composed of time and list of nodes,
            like [(time, nodes list), (time, nodes list), ...].
        """

        r_list_t_node, p_list_t_node= [], []
        for smi_node in tqdm(smi_nodes, desc='Filter sigal transition'):
            filt_sig = self.filterNoise(smi_nodes[smi_node], step)
            if len(set(filt_sig)) > 1:
                for i in range(0, step - 2):
                    node = smi_node.split('_')[1]
                    if filt_sig[i + 1] - filt_sig[i] == -1:
                        r_list_t_node.append((i + 1, json.loads(node)))
                    if filt_sig[i + 1] - filt_sig[i] == 1:
                        p_list_t_node.append((i + 1, json.loads(node)))
            
        r_list_t_node = sorted(r_list_t_node, key=operator.itemgetter(0))   
        p_list_t_node = sorted(p_list_t_node, key=operator.itemgetter(0))
        
        return r_list_t_node, p_list_t_node

    @staticmethod
    def uniteSig(list_t_node):
        """
        Unite signal at the same moment.

        :param list_t_node -> list
            list composed of time and node list
        :return dict_t_node -> dict
            Returns a dict which key is str_time and value is a list of nodes,
            like {time: [nodes list], ...}.
        """

        dict_t_node = {}
        for t_node in list_t_node:
            if str(t_node[0]) not in dict_t_node.keys():
                dict_t_node[str(t_node[0])] = t_node[1]
            else:
                dict_t_node[str(t_node[0])].extend(t_node[1])
        return dict_t_node