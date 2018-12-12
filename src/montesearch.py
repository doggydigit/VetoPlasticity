#!/usr/bin/env python

"""
    File name: montesearch.py
    Author: Matthias Tsai
    Email: matthias.chinyen.tsai@gmail.com
    Date created: 06/12/2018
    Date last modified: ...
    Python Version: 2.7
"""

import os
import dataset
from weight_class import *

param_names = ['Theta_high', 'Theta_low', 'A_LTP', 'A_LTD', 'tau_lowpass1', 'tau_lowpass2', 'tau_x']

grid_params = {'Theta_high': 11,
               'Theta_low': 11,
               'A_LTP': 9,
               'A_LTD': 9,
               'tau_lowpass1': 7,
               'tau_lowpass2': 7,
               'tau_x': 7}


def set_param(pname, index=None):
    """
    This function will either transform the given index to the corresponding parameter value or sample an index in the
    right parameter range to produce a randomly sampled value for the desired parameter.
    :param pname: Name of the parameter to set
    :param index: Index of the parameter value
    :return: return the value for desired parameter to set
    """
    if index is None:
        index = rnd.sample(xrange(grid_params[pname]), 1)[0] + 1

    if pname in ['Theta_high', 'Theta_low']:
        return (-15 - 5 * index) * mV
    elif pname in ['A_LTP', 'A_LTD']:
        return 0.001 * 3 ** index
    elif pname in ['tau_lowpass1', 'tau_lowpass2']:
        return 2 ** index * 2 * ms
    elif pname is 'tau_x':
        return 2 ** index * 0.5 * ms
    else:
        raise ValueError(pname)


def main(plasticity, neuron, veto, homeo=False, debug=False):
    """
    Simulate a stimulation protocol with the desired neuron model and plasticity rule to monitor the resulting synaptic
    weight dynamics. This can also plot or save some results of the simulation.
    :param plasticity: plasticity rule to use; can be either of 'Claire', 'Clopath' or 'Triplet'
    :param neuron: neuron model; can be either of default, LIF, exIF, Adex or Ziegler
    :param veto: bool whether or not to use a veto mechanism between LTP and LTD
    :param homeo: bool whether or not to use a homeostasis term in the plasticity
    :param debug: bool whether or not to have a more verbose output and simplified simulation
    """

    if veto:
        raise NotImplementedError

    # Stimulation protocol parameters
    protocols = [sTET_protocol_parameters, wTET_protocol_parameters,
                 sLFS_protocol_parameters, wLFS_protocol_parameters]

    # Presynaptic neuron parameters
    pre_neuron_parameters = ziegler_pre_parameters

    # Refractory Behaviour
    adaptation = {'LIF': False, 'exIF': False, 'Adex': True}[neuron]
    if neuron == 'Adex':
        neuron = 'exIF'

    # Postsynaptic neuron parameters
    post_neuron_parameters = {'default': default_post_neuron_parameters,
                              'LIF': LIF_post_parameters,
                              'exIF': exIF_post_parameters,
                              'Ziegler': ziegler_post_parameters}[neuron]

    # Plasticity parameters
    if plasticity is 'Clopath':
        raise NotImplementedError('give param names')
        # if post_neuron_parameters['model'] == 'exIF':
        #     parameters = Clopath_exIF_parameters
        # else:
        #     parameters = Clopath_LIF_parameters
    elif plasticity is 'Claire':
        if post_neuron_parameters['model'] == 'exIF':
            parameters = Claire_exIF_parameters
        else:
            parameters = Claire_LIF_parameters
    else:
        raise NotImplementedError('give param names')
        # parameters = Hippo_plasticity_parameters

    # Randomly initialize parameters
    indexes = {}
    for param_name in param_names:
        indexes[param_name] = rnd.sample(xrange(grid_params[param_name]), 1)[0] + 1
        parameters[param_name] = set_param(param_name, indexes[param_name])
    print('Parameters initialized')

    # Connect to database (Sqlite database corresponding to the plasticity model used)
    db_name = '../data/monteresults'
    db = dataset.connect('sqlite:///' + db_name)
    table_name = plasticity + '_veto' if veto else plasticity + '_noveto'
    the_table = db.create_table(table_name)

    # Monte-Carlo iterations
    current_score = 0
    nr_iterations = 100
    print('\nStarting Monte-Carlo optimization:')
    for i in range(nr_iterations):

        print('Iteration: {}'.format(i))

        # ##############################################################################################################
        #                                              Modify one parameter
        # ##############################################################################################################

        # Select parameter to modify and make copy of all parameters to work with
        new_parameters = dict(parameters)
        new_indexes = dict(indexes)
        param_name = rnd.sample(param_names, 1)[0]

        # Shift index of selected parameter and check whether it remains in accepted bounds. Else skip iteration.
        direction = bool(rnd.getrandbits(1))
        if direction:
            new_indexes[param_name] += 1
            if new_indexes[param_name] > grid_params[param_name]:
                continue

        else:
            new_indexes[param_name] -= 1
            if new_indexes[param_name] < 1:
                continue

        # Index shift is accepted, thus update parameter value
        new_parameters[param_name] = set_param(param_name)

        # Check whether this parameter configuration was already simulated.
        query = the_table.find_one(th=new_indexes['Theta_high'], tl=new_indexes['Theta_low'], ap=new_indexes['A_LTP'],
                                   ad=new_indexes['A_LTD'], t1=new_indexes['tau_lowpass1'],
                                   t2=new_indexes['tau_lowpass2'],
                                   tx=new_indexes['tau_x'])
        if query is None:

            # ##############################################################################################################
            #            Run Simulations of all 4 protocols with new parameters and save weight changes
            # ##############################################################################################################

            # Object containing the weight changes that occured in each protocol
            end_weights = {}
            new_score = 0
            broken = 0  # Actually unnecessary, but editor likes it

            # Iterate through all 4 protocols
            for protocol_parameters in protocols:

                # Silence print output of the PlasticityProtocol class
                if not debug:
                    sys.stdout = open(os.devnull, 'w')

                # Initialize main class
                ex = PlasticityProtocol(pre_neuron_parameters=pre_neuron_parameters,
                                        post_neuron_parameters=post_neuron_parameters,
                                        protocol_parameters=protocol_parameters,
                                        plasticity_parameters=new_parameters,
                                        same_connectivity_matrix=True,  # Repeatedly use same random seeds if True
                                        adaptation=adaptation,
                                        veto=veto,
                                        homeo=homeo,
                                        debug=debug)

                # Calibration of initial synaptic weights in order to produce a specific EPSP amplitude from a single
                # presynaptic spike. This value is defined in plasticity_parameters['v_increase_init'].
                ex.calibrate_w_init(std_cal=protocol_parameters['std'])
                initial_weights = ex.plasticity_parameters['init_weight']

                # Calibration of the fraction of presynaptic neurons triggered to spike by an extracellular stimulation
                #  pulse in order to lead to a specific percentage of postsynaptic neurons to fire
                ex.calibrate_amplitude(std_cal=protocol_parameters['std'])

                # Define which values will be monitored/saved during the simulation
                record_pre = None
                record_syn = 'w_ampa'
                record_post = 'v'

                # Run simulation
                m = ex.run(syn_parameters=record_syn, pre_parameters=record_pre, post_parameters=record_post,
                           pre_spikes=True, post_spikes=True)

                # Reenable printing of outputs
                if not debug:
                    sys.stdout = sys.__stdout__

                protocol = protocol_parameters['protocol_type']
                end_weights[protocol] = [np.mean(m['syn_monitor'].w_ampa[:, -1], 0),
                                         np.std(m['syn_monitor'].w_ampa[:, -1], 0)]

                if protocol in ['sTET', 'wTET']:
                    if end_weights[protocol][0] < initial_weights:
                        broken = True
                        break
                    else:
                        broken = False
                        protoscore = end_weights[protocol][0] - initial_weights
                elif protocol in ['sLFS', 'wLFS']:
                    if end_weights[protocol][0] < initial_weights:
                        broken = True
                        break
                    else:
                        broken = False
                        protoscore = initial_weights - end_weights[protocol][0]
                else:
                    raise ValueError(protocol)

                if protoscore > new_score:
                    new_score = protoscore

            if broken:
                new_score = 0

            the_table.insert(dict(th=new_indexes['Theta_high'], tl=new_indexes['Theta_low'], ap=new_indexes['A_LTP'],
                                  ad=new_indexes['A_LTD'], t1=new_indexes['tau_lowpass1'],
                                  t2=new_indexes['tau_lowpass2'],
                                  tx=new_indexes['tau_x'], score=new_score))
            db.commit()
        else:
            # Get score that was already computed
            new_score = query['score']

        # Given appropriate probability, update current state with the new state
        if new_score > current_score:
            parameters = new_parameters
            indexes = new_indexes
            current_score = new_score
        else:
            try:
                accept_prob = new_score / current_score
            except ZeroDivisionError:
                accept_prob = 1

            if rnd.uniform(0, 1) < accept_prob:
                parameters = new_parameters
                indexes = new_indexes
                current_score = new_score

        print('Score: {}'.format(current_score))


if __name__ == "__main__":
    # Simulation choices
    rule_name = 'Claire'  # can be either of 'Claire', 'Clopath' or 'Triplet'
    neuron_types = 'Adex'  # can be either of default, LIF, exIF, Ziegler
    vetoing = False  # whether or not to use a veto mechanism between LTP and LTD

    # Run
    main(rule_name, neuron=neuron_types, veto=vetoing, debug=False)