#!/usr/bin/env python

"""
    File name: montesearch.py
    Author: Matthias Tsai
    Email: matthias.chinyen.tsai@gmail.com
    Date created: 06/12/2018
    Date last modified: 14/12/2018
    Python Version: 2.7
"""

from brian2.utils.logger import catch_logs
from simulation import *
import os
import csv
import dataset
import warnings

warnings.filterwarnings("error")


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


def main(plasticity,  new_indexes, neuron, veto, homeo=False, debug=False):
    """
    Simulate a stimulation protocol with the desired neuron model and plasticity rule to monitor the resulting synaptic
    weight dynamics. This can also plot or save some results of the simulation.
    :param plasticity: plasticity rule to use; can be either of 'Claire', 'Clopath' or 'Triplet'
    :param new_indexes: indices of parameters to recompute
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
        # if post_neuron_parameters['model'] == 'exIF':
        #     parameters = Clopath_exIF_parameters
        # else:
        #     parameters = Clopath_LIF_parameters
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
        # parameters = Hippo_plasticity_parameters

    # Initialize parameters
    for param_name in param_names:
        parameters[param_name] = set_param(param_name, new_indexes[param_name])

    print('Parameters initialized')

    # Connect to database (Sqlite database corresponding to the plasticity model used)
    db_name = '../data/monteresults.db'
    db = dataset.connect('sqlite:///' + db_name)
    table_name = plasticity + '_veto' if veto else plasticity + '_noveto'
    the_table = db.create_table(table_name)

    # Randomized seed (because it might not be due to Plasticity class
    rnd.seed()

    # Check whether this parameter configuration was already simulated.
    if veto:
        query = the_table.find_one(th=new_indexes['Theta_high'], tl=new_indexes['Theta_low'],
                                   ap=new_indexes['A_LTP'], ad=new_indexes['A_LTD'], t1=new_indexes['tau_lowpass1'],
                                   t2=new_indexes['tau_lowpass2'], tx=new_indexes['tau_x'],
                                   bt=new_indexes['b_theta'], tt=new_indexes['tau_theta'])
    else:
        query = the_table.find_one(th=new_indexes['Theta_high'], tl=new_indexes['Theta_low'],
                                   ap=new_indexes['A_LTP'], ad=new_indexes['A_LTD'], t1=new_indexes['tau_lowpass1'],
                                   t2=new_indexes['tau_lowpass2'], tx=new_indexes['tau_x'])

    if query is None:
        print('DIDNT EXIST!!!')
        print(new_indexes)
        return 5

    else:

        # ##########################################################################################################
        #            Run Simulations of all 4 protocols with new parameters and save weight changes
        # ##########################################################################################################

        # Object containing the weight changes that occured in each protocol
        end_weights = {}
        broken = True  # Actually unnecessary, but editor likes it
        new_score = 666
        nan_bool = False

        # Iterate through all 4 protocols
        for protocol_parameters in protocols:

            # Silence print output (of the PlasticityProtocol class)
            if not debug:
                sys.stdout = open(os.devnull, 'w')

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

            # Monitor Brian warnings especially due to NaN numerical integration errors
            with catch_logs() as brian_warnings:

                # Make sure to catch any error due to initial weights calibration
                try:
                    # Calibration of initial synaptic weights in order to produce a specific EPSP amplitude from a
                    # single presynaptic spike. This value is defined in plasticity_parameters['v_increase_init'].
                    ex.calibrate_w_init(std_cal=protocol_parameters['std'])
                except SystemError:
                    if not debug:
                        sys.stdout = sys.__stdout__
                    print('#####################################################################################\n'
                          '                         Weights initialized too small\n'
                          '#####################################################################################\n')
                    print(new_indexes)
                    print(parameters)
                    return 1

                # Save initial weights for later computation of the protocol plasticity and score
                initial_weights = ex.plasticity_parameters['init_weight']

                # Make sure to catch any error due to extracellular stimulation calibration
                try:
                    # Calibration of the fraction of presynaptic neurons triggered to spike by an extracellular
                    # stimulation pulse in order to lead to a specific percentage of postsynaptic neurons to fire
                    ex.calibrate_amplitude(std_cal=protocol_parameters['std'])
                except SystemError:
                    if not debug:
                        sys.stdout = sys.__stdout__
                    print('#####################################################################################\n'
                          '                   Initial extracellular stimulation too small\n'
                          '#####################################################################################\n')
                    print(new_indexes)
                    print(parameters)
                    return 1

                # Run simulation
                try:
                    ex.run()
                except Warning:
                    nan_bool = True

                if len(brian_warnings) > 0:
                    nan_bool = True

            # Reenable printing of outputs
            if not debug:
                sys.stdout = sys.__stdout__

            # If TET protocol produced LTD or LFS protocol produced LTP, exit and assign score of zero
            # Else compute score for protocol
            protocol = protocol_parameters['protocol_type']
            end_weights[protocol] = np.mean(np.array(ex.syn.w_ampa))

            # If plasticity parameters lead to extreme dynamics producing NaNs, exit and assign score of zero
            if nan_bool or end_weights[protocol] != end_weights[protocol]:
                print('Nan solution\n')
                broken = True
                break

            if protocol in ['sTET', 'wTET']:
                if end_weights[protocol] < initial_weights:
                    broken = True
                    break
                else:
                    broken = False
                    protoscore = end_weights[protocol] - initial_weights
            elif protocol in ['sLFS', 'wLFS']:
                if end_weights[protocol] > initial_weights:
                    broken = True
                    break
                else:
                    broken = False
                    protoscore = initial_weights - end_weights[protocol]
            else:
                raise ValueError(protocol)

            # Only keep the worst score of all protocols
            if protoscore < new_score:
                if np.isnan(protoscore):
                    print('Escaped NaN caught')
                    broken = True
                    break
                else:
                    new_score = protoscore
            elif protoscore > 1:
                print('Weird protoscore = {}'.format(protoscore))
                print('The protocol is ' + protocol)
                print('\nParameters:')
                print(parameters)
                print('\nIndexes:')
                print(new_indexes)
                return 3

        if broken:
            new_score = 0

        if new_score == 666:
            print('#############################################################################################\n'
                  '                               666 score escaped\n'
                  '#############################################################################################\n')
            the_table.update(dict(id=query['id'], score=0), ['id'])
            db.commit()
            return 4
        else:
            the_table.update(dict(id=query['id'], score=new_score), ['id'])
            db.commit()

    print('    Score = {}'.format(new_score))

    return 0


if __name__ == "__main__":

    line = int(float(sys.argv[1]))

    inds = []
    i = 0
    with open('recompute.csv', 'rb') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in spamreader:
            i = i + 1
            if i == line:
                inds = row
                break

    indexes = {'A_LTD': int(inds[0]),
               'tau_x': int(inds[1]),
               'tau_theta': int(inds[2]),
               'tau_lowpass1': int(inds[4]),
               'tau_lowpass2': int(inds[3]),
               'b_theta': int(inds[5]),
               'A_LTP': int(inds[6]),
               'Theta_low': int(inds[7]),
               'Theta_high': int(inds[8])}

    print(indexes)

    # Simulation choices
    rule_name = 'Claire'  # can be either of 'Claire', 'Clopath' or 'Triplet'
    neuron_types = 'Adex'  # can be either of default, LIF, exIF, Ziegler
    vetoing = True  # whether or not to use a veto mechanism between LTP and LTD

    # Run
    exi = main(rule_name, indexes, neuron=neuron_types, veto=vetoing, debug=False)

    if exi is 0:
        print('\nMonte-Carlo search finished successfully!')
    else:
        print('\nAn error occured...')
