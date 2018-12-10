#!/usr/bin/env python

"""
    File name: ziegler.py
    Author: Matthias Tsai
    Email: matthias.chinyen.tsai@gmail.com
    Date created: 27/11/2018
    Date last modified: ...
    Python Version: 2.7
"""

import sys
import numpy as np
import weight_class as wc


def main(plasticity, stimulation, neuron, plotting, saving, veto, homeo, debug):
    """
    Simulate a stimulation protocol with the desired neuron model and plasticity rule to monitor the resulting synaptic
    weight dynamics. This can also plot or save some results of the simulation.
    :param plasticity: plasticity rule to use; can be either of 'Claire', 'Clopath' or 'Triplet'
    :param stimulation: stimulation protocol; can be either of test, sTET, wTET, sLFS, wLFS or STDP
    :param neuron: neuron model; can be either of default, LIF, exIF, Adex or Ziegler
    :param plotting: bool whether or not to plot some results of the simulation
    :param saving: bool whether or not to save the weight dynamics
    :param veto: bool whether or not to use a veto mechanism between LTP and LTD
    :param homeo: bool whether or not to use a homeostasis term in the plasticity
    :param debug: bool whether or not to have a more verbose output and simplified simulation
    """

    # Stimulation protocol parameters
    if stimulation:
        protocol_parameters = {'sTET': wc.sTET_protocol_parameters,
                               'wTET': wc.wTET_protocol_parameters,
                               'sLFS': wc.sLFS_protocol_parameters,
                               'wLFS': wc.wLFS_protocol_parameters,
                               'STDP': wc.default_protocol_parameters,
                               'test': wc.testing_protocol_parameters}[stimulation]
    else:
        raise ValueError

    # Presynaptic neuron parameters
    pre_neuron_parameters = wc.ziegler_pre_parameters

    # Refractory Behaviour
    adaptation = {'LIF': False, 'exIF': False, 'Adex': True}[neuron]
    if neuron == 'Adex':
        neuron = 'exIF'

    # Postsynaptic neuron parameters
    post_neuron_parameters = {'default': wc.default_post_neuron_parameters,
                              'LIF': wc.LIF_post_parameters,
                              'exIF': wc.exIF_post_parameters,
                              'Ziegler': wc.ziegler_post_parameters}[neuron]

    # Plasticity parameters
    if plasticity is 'Clopath':
        if post_neuron_parameters['model'] == 'exIF':
            plasticity_parameters = wc.Clopath_exIF_parameters
        else:
            plasticity_parameters = wc.Clopath_LIF_parameters
    elif plasticity is 'Claire':
        if post_neuron_parameters['model'] == 'exIF':
            plasticity_parameters = wc.Claire_exIF_parameters
        else:
            plasticity_parameters = wc.Claire_LIF_parameters
    else:
        plasticity_parameters = wc.Hippo_plasticity_parameters

    # Initialize main class
    ex = wc.PlasticityProtocol(pre_neuron_parameters=pre_neuron_parameters,
                               post_neuron_parameters=post_neuron_parameters,
                               protocol_parameters=protocol_parameters,
                               plasticity_parameters=plasticity_parameters,
                               same_connectivity_matrix=True,  # Repeatedly use same random seeds if True
                               adaptation=adaptation,
                               veto=veto,
                               homeo=homeo,
                               debug=debug)

    # Calibration of initial synaptic weights in order to produce a specific EPSP amplitude from a single presynaptic
    # spike. This value is defined in plasticity_parameters['v_increase_init'].
    ex.calibrate_w_init(std_cal=protocol_parameters['std'])

    # Calibration of the fraction of presynaptic neurons triggered to spike by an extracellular stimulation pulse in
    # order to lead to a specific percentage of postsynaptic neurons to fire
    ex.calibrate_amplitude(std_cal=protocol_parameters['std'])

    # Define which values will be monitored/saved during the simulation
    record_pre = None
    record_syn = 'w_ampa'
    record_post = 'v'

    # Run simulation
    m = ex.run(syn_parameters=record_syn, pre_parameters=record_pre, post_parameters=record_post, pre_spikes=True,
               post_spikes=True)
    print('Simulation terminated successully!\n')

    # Output results
    w_single_end = m['syn_monitor'].w_ampa[:, -1]  # these are the weights at the end of LTP
    print("Average weight values after protocol is " + str(np.mean(w_single_end)))

    # Plot Results
    if plotting:
        
        from plot_triplet_parameters import *

        # Plot final weight distribution
        plt.hist(w_single_end, bins=50, alpha=0.3)
        plt.ylabel('Occurrency', fontsize=14)
        plt.xlabel('Weight value after ' + stimulation + ' plasticity protocol', fontsize=14)
        plt.show()
        plt.savefig("../fig/weight_hist_after_" + stimulation + "_plasticity.pdf", format="PDF")

        # Plot spikes in time
        plot_spikes(m, title='Spikes,Nb={},d={},NB={},D={}'.format(protocol_parameters['nr_pulses'],
                                                                   protocol_parameters['hro'],
                                                                   protocol_parameters['nr_blocks'],
                                                                   ex.timestep))
        plt.show()

        # Plot weight dynamics
        nr = m['syn_monitor'].w_ampa.shape[1]
        time_axis = list(range(nr))
        w_avg = np.mean(m['syn_monitor'].w_ampa, 0)
        w_std = np.std(m['syn_monitor'].w_ampa, 0)
        kwargs = {'color': 'y'}
        plt.figure()
        plt.plot(w_avg, 'k')
        plt.fill_between(time_axis, np.clip(w_avg - w_std, 0, None), w_avg + w_std, **kwargs)
        plt.show()

    if debug:
        
        from plot_triplet_parameters import *

        # Plot voltage trace of a postsynaptic neuron
        for i in range(100):
            plt.plot(m['post_monitor'].v[i][9000:])
            plt.show()

    # Save monitored weights
    if saving:
        import pickle
        v = '_veto' if veto else ''
        with open('../data/' + protocol + v + '.pickle', 'wb') as handle:
            pickle.dump(m['syn_monitor'].w_ampa, handle, protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == "__main__":

    # Simulation choices
    if len(sys.argv) == 1:
        protocol = 'test'  # can be either of test, sTET, wTET, sLFS, wLFS or STDP
    else:
        protocol = sys.argv[1]
    rule_name = 'Claire'  # can be either of 'Claire', 'Clopath' or 'Triplet'
    neuron_types = 'Adex'  # can be either of default, LIF, exIF, Ziegler
    plot = False  # whether or not to plot some results of the simulation
    save = True  # whether or not to save the weight dynamics

    # Run
    main(rule_name, protocol, neuron=neuron_types, plotting=plot, saving=save, veto=False, homeo=False, debug=False)
