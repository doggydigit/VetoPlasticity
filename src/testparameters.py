#!/usr/bin/env python
"""
    File name: testparameters.py
    Author: Matthias Tsai
    Email: matthias.chinyen.tsai@gmail.com
    Date created: 17/12/2018
    Date last modified: 20/12/2018
    Python Version: 2.7
"""

from simulation import *


def set_param(pname, index, granu=0):
    """
    This function will either transform the given index to the corresponding parameter value or sample an index in the
    right parameter range to produce a randomly sampled value for the desired parameter.
    :param pname: Name of the parameter to set
    :param index: Index of the parameter value
    :param granu: Granularity of the parameter search, which determines what parameter search subspace
    :return: return the value for desired parameter to set
    """

    if granu == 0:
        if pname in ['Theta_high', 'Theta_low']:
            return (-15 - 7 * index) * mV
        elif pname in ['A_LTP', 'A_LTD']:
            return 0.001 * 10 ** index
        elif pname in ['tau_lowpass1', 'tau_lowpass2']:
            return 2 ** index * ms
        elif pname is 'tau_x':
            return 2 ** (index - 2) * ms
        elif pname is 'b_theta':
            return 0.4 * 5 ** index * ms
        elif pname is 'tau_theta':
            return 0.2 * 5 ** index * ms
        else:
            raise ValueError(pname)
    else:
        raise NotImplementedError


def main(plasticity, neuron, veto, homeo=False, debug=False):
    """
    Simulate a stimulation protocol with the desired neuron model and plasticity rule to monitor the resulting synaptic
    weight dynamics. This can also plot or save some results of the simulation.
    :param plasticity: plasticity rule to use; can be either of 'Claire', 'Clopath' or 'Triplet'
    :param neuron: neuron model; can be either of default, LIF, exIF, Adex or Ziegler
    :param veto: bool whether or not to use a veto mechanism between LTP and LTD
    :param homeo: bool whether or not to use a homeostasis term in the plasticity
    :param debug: bool whether or not to have a more verbose output and simplified simulation
    :param granularity: int describing how fine the param search will be (Note: the range is also decreased for that)
    """

    # Stimulation protocol parameters
    protocols = [wLFS_protocol_parameters, wTET_protocol_parameters,
                 sLFS_protocol_parameters, sTET_protocol_parameters]

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
    elif plasticity is 'Claire':
        if veto:
            param_names = ['Theta_high', 'Theta_low', 'A_LTP', 'A_LTD', 'tau_lowpass1', 'tau_lowpass2', 'tau_x',
                           'b_theta', 'tau_theta']
        else:
            param_names = ['Theta_high', 'Theta_low', 'A_LTP', 'A_LTD', 'tau_lowpass1', 'tau_lowpass2', 'tau_x']
        if post_neuron_parameters['model'] == 'exIF':
            parameters = Claire_exIF_parameters
        else:
            parameters = Claire_LIF_parameters
    else:
        raise NotImplementedError('give param names')

    # Randomly initialize parameters
    indexes = {'A_LTD': 8,
               'tau_x': 5,
               'tau_theta': 1,
               'tau_lowpass1': 4,
               'tau_lowpass2': 4,
               'b_theta': 2,
               'A_LTP': 7,
               'Theta_low': 7,
               'Theta_high': 4}
    for param_name in param_names:
        parameters[param_name] = set_param(param_name, indexes[param_name])

    print('Parameters initialized')

    # Iterate through all 4 protocols
    for protocol_parameters in protocols:

        # Initialize main class
        ex = PlasticityProtocol(pre_neuron_parameters=pre_neuron_parameters,
                                post_neuron_parameters=post_neuron_parameters,
                                protocol_parameters=protocol_parameters,
                                plasticity_parameters=parameters,
                                same_connectivity_matrix=True,  # Repeatedly use same random seeds if True
                                adaptation=adaptation,
                                veto=veto,
                                homeo=homeo,
                                debug=debug)

        # Calibration of initial synaptic weights in order to produce a specific EPSP amplitude from a
        # single presynaptic spike. This value is defined in plasticity_parameters['v_increase_init'].
        ex.calibrate_w_init(std_cal=protocol_parameters['std'])

        # Save initial weights for later computation of the protocol plasticity and score
        initial_weights = ex.plasticity_parameters['init_weight']

        # Calibration of the fraction of presynaptic neurons triggered to spike by an extracellular
        # stimulation pulse in order to lead to a specific percentage of postsynaptic neurons to fire
        ex.calibrate_amplitude(std_cal=protocol_parameters['std'])

        ex.run()

        # If TET protocol produced LTD or LFS protocol produced LTP, exit and assign score of zero
        # Else compute score for protocol
        protocol = protocol_parameters['protocol_type']
        endw = np.mean(np.array(ex.syn.w_ampa))
        endstd = np.std(ex.syn.w_ampa)

        if protocol in ['sTET', 'wTET']:
            protoscore = endw - initial_weights
        elif protocol in ['sLFS', 'wLFS']:
            protoscore = initial_weights - endw
        else:
            raise ValueError(protocol)

        print('For protocol "{}", we had a score of {} and std of {}'.format(protocol, protoscore, endstd))


if __name__ == "__main__":
    # Simulation choices
    rule_name = 'Claire'  # can be either of 'Claire', 'Clopath' or 'Triplet'
    neuron_types = 'Adex'  # can be either of default, LIF, exIF, Ziegler
    vetoing = True  # whether or not to use a veto mechanism between LTP and LTD

    # Run
    main(rule_name, neuron=neuron_types, veto=vetoing, debug=False)

    print('Done.')
