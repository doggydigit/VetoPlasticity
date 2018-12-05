#!/usr/bin/env python

"""
    File name: ziegler.py
    Author: Matthias Tsai
    Email: matthias.chinyen.tsai@gmail.com
    Date created: 27/11/2018
    Date last modified: 4/25/2013
    Python Version: 2.7
"""

import weight_class as wc
from plot_triplet_parameters import *


def main(plasticity, stimulation, neuron, path, veto, homeo, debug):

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
            plasticity_parameters = wc.Claire_parameters
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

    # Output results

    # w_single_end = m['syn_monitor'].w_ampa[:, -1]  # these are the weights at the end of LTP
    #
    # print(w_single_end)
    #
    # plt.rcParams['figure.figsize'] = (22, 5)
    #
    # # plot final weight distribution
    # plt.hist(w_single_end, bins=50, alpha=0.3)
    # plt.ylabel('Occurrency', fontsize=14)
    # plt.xlabel('Weight value at the end of LTP', fontsize=14)
    # plt.show()
    # if path is not None:
    #     plt.savefig(path + "weight_hist_after_LTP.pdf", format="PDF")
    #
    # # Plot spikes in time
    # plt.figure()
    # plot_spikes(m, title='Spikes,Nb={},d={},NB={},D={}'.format(wc.testing_protocol_parameters['nr_pulses'],
    #                                                            wc.testing_protocol_parameters['hro'],
    #                                                            wc.testing_protocol_parameters['nr_blocks'],
    #                                                            wc.testing_protocol_parameters['dt']))
    # plt.show()
    #
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

    # Plot voltage trace of a postsynaptic neuron
    for i in range(100):
        plt.figure()
        plt.plot(m['post_monitor'].v[i])
        plt.show()


if __name__ == "__main__":
    # Simulation choices
    rule_name = 'Clopath'  # can be either of 'Claire', 'Clopath' or 'Triplet'
    protocol = 'test'  # can be either of test, sTET, wTET, sLFS, wLFS or STDP
    neuron_types = 'exIF'  # can be either of default, LIF, exIF, Ziegler
    save_path = '/home/chi/Documents/Study/LCN/project/fig/'  # directory where to save resulting figures

    # Run
    main(rule_name, protocol, neuron_types, save_path, veto=True, homeo=False, debug=True)
