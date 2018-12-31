#!/usr/bin/env python

"""
    File name: montesearchcal.py
    Author: Matthias Tsai
    Email: matthias.chinyen.tsai@gmail.com
    Date created: 14/12/2018
    Date last modified: 20/12/2018
    Python Version: 2.7
"""

from brian2.utils.logger import catch_logs
from simulation import *
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


def main(plasticity, neuron, veto, homeo=False, debug=False, granularity=0, first_id=None):
    """
    Simulate a stimulation protocol with the desired neuron model and plasticity rule to monitor the resulting synaptic
    weight dynamics. This can also plot or save some results of the simulation.
    :param plasticity: plasticity rule to use; can be either of 'Claire', 'Clopath' or 'Triplet'
    :param neuron: neuron model; can be either of default, LIF, exIF, Adex or Ziegler
    :param veto: bool whether or not to use a veto mechanism between LTP and LTD
    :param homeo: bool whether or not to use a homeostasis term in the plasticity
    :param debug: bool whether or not to have a more verbose output and simplified simulation
    :param granularity: int describing how fine the param search will be (Note: the range is also decreased for that)
    :param first_id: ID of the parameter configuration to start with. If None, a random configuration is used.
    """

    # Connect to database (Sqlite database corresponding to the plasticity model used)
    if granularity == 0:
        db_name = '../data/monteresults.db'
    elif granularity == 1:
        db_name = '../data/monteresults_g1.db'
    else:
        raise NotImplementedError
    db = dataset.connect('sqlite:///' + db_name)
    table_name = plasticity + '_veto' if veto else plasticity + '_noveto'
    the_table = db.create_table(table_name)

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
            if granularity == 0:
                grid_params = {'Theta_high': 8, 'Theta_low': 8, 'A_LTP': 8, 'A_LTD': 8, 'tau_lowpass1': 7,
                               'tau_lowpass2': 7, 'tau_x': 7, 'b_theta': 5, 'tau_theta': 5}
            elif granularity == 1:
                grid_params = {'Theta_high': 8, 'Theta_low': 8, 'A_LTP': 8, 'A_LTD': 8, 'tau_lowpass1': 7,
                               'tau_lowpass2': 7, 'tau_x': 7, 'b_theta': 5, 'tau_theta': 5}
            else:
                raise NotImplementedError
        else:
            param_names = ['Theta_high', 'Theta_low', 'A_LTP', 'A_LTD', 'tau_lowpass1', 'tau_lowpass2', 'tau_x']
            if granularity == 0:
                grid_params = {'Theta_high': 11, 'Theta_low': 11, 'A_LTP': 7, 'A_LTD': 7, 'tau_lowpass1': 7,
                               'tau_lowpass2': 7, 'tau_x': 7}
            else:
                raise NotImplementedError

        if post_neuron_parameters['model'] == 'exIF':
            parameters = Claire_exIF_parameters
        else:
            parameters = Claire_LIF_parameters
    else:
        raise NotImplementedError('give param names')
        # parameters = Hippo_plasticity_parameters

    # Initialize parameter indices
    indexes = {}
    if first_id is None:
        for param_name in param_names:
            indexes[param_name] = rnd.sample(xrange(grid_params[param_name]), 1)[0] + 1
    else:
        translator = {'Theta_high': 'th', 'Theta_low': 'tl', 'A_LTP': 'ap', 'A_LTD': 'ad', 'tau_lowpass1': 't1',
                      'tau_lowpass2': 't2', 'tau_x': 'tx', 'b_theta': 'bt', 'tau_theta': 'tt'}
        first_indices = the_table.find_one(id=first_id)
        for param_name in param_names:
            indexes[param_name] = first_indices[translator[param_name]]

        print('First indices are:\n{}'.format(first_indices))

    # Initialize parameter values from indices according to desired grid design and specfici granularity
    for param_name in param_names:
        parameters[param_name] = set_param(param_name, indexes[param_name], granularity)

    print('Parameters initialized')

    # Monte-Carlo iterations
    current_score = 0
    nr_iterations = 300
    print('\nStarting Monte-Carlo optimization:')
    for i in range(nr_iterations):

        print('Iteration: {}'.format(i))
        sys.stdout.flush()

        # ##############################################################################################################
        #                                              Modify one parameter
        # ##############################################################################################################

        # Randomized seed (because it might not be due to Plasticity class
        rnd.seed()

        # Select parameter to modify and make copy of all parameters to work with
        new_parameters = dict(parameters)
        new_indexes = dict(indexes)
        param_name = rnd.sample(param_names, 1)[0]

        # Shift index of selected parameter and check whether it remains in accepted bounds. Else skip iteration.
        direction = bool(rnd.getrandbits(1))
        if direction:
            new_indexes[param_name] += 1
            if new_indexes[param_name] > grid_params[param_name]:
                if current_score > 0:
                    print('Wall reached with parameter {} for index {}'.format(param_name, new_indexes[param_name]))
                continue

        else:
            new_indexes[param_name] -= 1
            if new_indexes[param_name] < 1:
                if current_score > 0:
                    print('Wall reached with parameter {} for index {}'.format(param_name, new_indexes[param_name]))
                continue

        # Index shift is accepted, thus update parameter value
        new_parameters[param_name] = set_param(param_name, new_indexes[param_name])

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

            # Create that row and temporarily put a score of zero to prevent other processors to compute it again
            if veto:
                query_id = the_table.insert(dict(th=new_indexes['Theta_high'], tl=new_indexes['Theta_low'],
                                                 ap=new_indexes['A_LTP'], ad=new_indexes['A_LTD'],
                                                 t1=new_indexes['tau_lowpass1'], t2=new_indexes['tau_lowpass2'],
                                                 tx=new_indexes['tau_x'], bt=new_indexes['b_theta'],
                                                 tt=new_indexes['tau_theta'], score=-1))

            else:
                query_id = the_table.insert(dict(th=new_indexes['Theta_high'], tl=new_indexes['Theta_low'],
                                                 ap=new_indexes['A_LTP'], ad=new_indexes['A_LTD'],
                                                 t1=new_indexes['tau_lowpass1'], t2=new_indexes['tau_lowpass2'],
                                                 tx=new_indexes['tau_x'], score=-1))
            db.commit()

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

                # Monitor Brian warnings especially due to NaN numerical integration errors
                with catch_logs() as brian_warnings:

                    # Make sure to catch any error due to initial weights calibration
                    try:
                        # Calibration of initial synaptic weights in order to produce a specific EPSP amplitude from a
                        # single presynaptic spike. This value is defined in plasticity_parameters['v_increase_init'].
                        ex.calibrate_w_init(std_cal=protocol_parameters['std'])
                    except SystemError:
                        the_table.delete(id=query_id)
                        db.commit()
                        print('#####################################################################################\n'
                              '                         Weights initialized too small\n'
                              '#####################################################################################\n')
                        print(new_indexes)
                        print(new_parameters)
                        return 1

                    # Save initial weights for later computation of the protocol plasticity and score
                    initial_weights = ex.plasticity_parameters['init_weight']

                    # Make sure to catch any error due to extracellular stimulation calibration
                    try:
                        # Calibration of the fraction of presynaptic neurons triggered to spike by an extracellular
                        # stimulation pulse in order to lead to a specific percentage of postsynaptic neurons to fire
                        ex.calibrate_amplitude(std_cal=protocol_parameters['std'])
                    except SystemError:
                        the_table.delete(id=query_id)
                        db.commit()
                        print('#####################################################################################\n'
                              '                   Initial extracellular stimulation too small\n'
                              '#####################################################################################\n')
                        print(new_indexes)
                        print(new_parameters)
                        return 1

                    # Run simulation
                    try:
                        ex.run()
                    except Warning:
                        nan_bool = True

                    if len(brian_warnings) > 0:
                        nan_bool = True

                # If plasticity parameters lead to extreme dynamics producing NaNs, exit and assign score of zero
                if nan_bool:
                    print('Nan solution\n')
                    broken = True
                    break

                # If TET protocol produced LTD or LFS protocol produced LTP, exit and assign score of zero
                # Else compute score for protocol
                protocol = protocol_parameters['protocol_type']
                end_weights[protocol] = np.mean(np.array(ex.syn.w_ampa))
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
                    the_table.delete(id=query_id)
                    db.commit()
                    print('Weird protoscore = {}'.format(protoscore))
                    print('The protocol is ' + protocol)
                    print('\nParameters:')
                    print(new_parameters)
                    print('\nIndexes:')
                    print(new_indexes)
                    return 3

            if broken:
                new_score = 0

            the_table.update(dict(id=query_id, score=new_score), ['id'])
            db.commit()

        else:
            # Get score that was already computed
            new_score = query['score']
            print('    Was already simulated')

        # Given appropriate probability, update current state with the new state
        if new_score > current_score:
            parameters = new_parameters
            indexes = new_indexes
            current_score = new_score
        else:
            try:
                accept_prob = new_score / current_score
            except ZeroDivisionError or RuntimeWarning:
                accept_prob = 1

            if rnd.uniform(0, 1) < accept_prob:
                parameters = new_parameters
                indexes = new_indexes
                current_score = new_score

        print('    Score = {}'.format(current_score))
        if current_score == 666:
            print('#################################################################################################\n'
                  '                               666 score escaped\n'
                  '#################################################################################################\n')

    return 0


if __name__ == "__main__":

    g = int(sys.argv[1])  # Resolution of the grid search

    if len(sys.argv) == 3:
        fid = int(sys.argv[2])
    else:
        fid = None

    # Simulation choices
    rule_name = 'Claire'  # can be either of 'Claire', 'Clopath' or 'Triplet'
    neuron_types = 'Adex'  # can be either of default, LIF, exIF, Ziegler
    vetoing = True  # whether or not to use a veto mechanism between LTP and LTD

    # Run
    exi = main(rule_name, neuron=neuron_types, veto=vetoing, debug=False, granularity=g, first_id=fid)

    if exi is 0:
        print('\nMonte-Carlo search finished successfully!')
    else:
        print('\nAn error occured...')
