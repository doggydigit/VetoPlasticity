#!/usr/bin/env python

"""
    File name: plot_triplet_parameters.py
    Author: Halla Sigurthorsdottir, Christoph Blattgerste, Giorgia Dellaferrera
    Date created: 06/10/2016
    Date last modified: 18/09/2018
    Python Version: 2.7
"""

import random
import numpy as np
import matplotlib.pyplot as plt


def plot_triplet_p(m, title=None):
    plt.subplot(221)
    plt.plot(m['syn_monitor'].t, m['syn_monitor'].o_1[0])
    plt.ylabel('o_1')

    plt.subplot(222)
    plt.plot(m['syn_monitor'].t, m['syn_monitor'].o_2[0])
    plt.ylabel('o_2')

    plt.subplot(223)
    plt.plot(m['syn_monitor'].t, m['syn_monitor'].r_1[0])
    plt.ylabel('r_1')
    plt.xlabel('t (s)')

    plt.subplot(224)
    plt.plot(m['syn_monitor'].t, m['syn_monitor'].r_2[0])
    plt.ylabel('r_2')
    plt.xlabel('t (s)')

    plt.title(title)

    plt.show()


def plot_prepost_spikes(m, voltage=True, title=None):
    '''
    This function plots the pre- and postsynaptic spikes as well as the voltage of the postsynaptic neuron (if voltage is True).
    '''
    plt.plot(m['pre_spikes'].t, -0.07 * np.ones(len(m['pre_spikes'].t)), '*')
    leg = ['Presynaptic spikes']

    if len(m['post_spikes'].t) > 0:
        plt.plot(m['post_spikes'].t, -0.07 * np.ones(len(m['post_spikes'].t)), '*')
        leg += ['Postsynaptic_spikes']

    if voltage:
        for i in range(m['post_monitor'].v.shape[0]):
            plt.plot(m['post_monitor'].t, m['post_monitor'].v[i])
            leg += ['Voltage of neuron' + str(i)]
            plt.ylabel('Voltage (V)')

    plt.legend(leg, loc='upper right')
    plt.xlabel('Time (s)')
    if not title is None:
        plt.title(title)

    plt.show()


def plot_w(m, no_syn=None, title=None):
    '''
    This function plots the weight change over time.
    '''
    leg = [0] * m['syn_monitor'].w_ampa.shape[0]
    if no_syn == None:
        no_syn = m['syn_monitor'].w_ampa.shape[0]
    else:
        no_syn = no_syn

    for i in range(no_syn):
        plt.figure()
        plt.plot(m['syn_monitor'].t, m['syn_monitor'].w_ampa[i])
        leg[i] = 'Synapse ' + str(i)

        plt.xlabel('Time (s)')
        plt.ylabel('Weight')
        plt.legend(leg)
        if not title is None:
            plt.title(title)

        plt.show()


def plot_rs_os(m, no_syn=None, title=None):
    '''
    This function plots the change of r_1, r_2, o_1, o_2 over time.
    '''
    leg = [0] * m['syn_monitor'].w_ampa.shape[0]
    if no_syn == None:
        no_syn = m['syn_monitor'].w_ampa.shape[0]
    else:
        no_syn = no_syn

    rs_and_os = True

    if rs_and_os == True:
        for i in range(no_syn):
            plt.figure()
            plt.plot(m['syn_monitor'].t, m['syn_monitor'].r_1[i], label='r_1')
            plt.plot(m['syn_monitor'].t, m['syn_monitor'].r_2[i], label='r_2')
            plt.plot(m['syn_monitor'].t, m['syn_monitor'].o_1[i], label='o_1')
            plt.plot(m['syn_monitor'].t, m['syn_monitor'].o_2[i], label='o_2')
            # leg[i] = 'Synapse ' + str(i)

            plt.xlabel('Time (s)')
            plt.ylabel('Weight')
            plt.legend()
            if not title is None:
                plt.title(title)

            plt.show()


def plot_spikes(m, save_in_path=False, path=None, show=False, title=None, no_poisson_spikes=15, poisson_pop=False):
    '''
    This function plots 20 (or input) pre spike arrivals and 10 (or input) post spike times.
    '''

    if len(m['post_spikes'].spike_trains()) < 10:
        no_post_spikes = len(m['post_spikes'].spike_trains())
    else:
        no_post_spikes = 10
    if len(m['pre_spikes'].spike_trains()) < 10:
        no_pre_spikes = len(m['pre_spikes'].spike_trains())
    else:
        no_pre_spikes = 10

    no_spikes = no_pre_spikes + no_post_spikes
    print('number of pre spike trains: {}'.format(no_pre_spikes))
    print('number of post spike trains: {}'.format(no_post_spikes))

    if poisson_pop == True:
        no_spikes += no_poisson_spikes

    # Prepare pre spike list
    pre_ids = random.sample(range(len(m['pre_spikes'].spike_trains())), no_pre_spikes)
    pre_spikes = [0] * no_pre_spikes
    for i in range(no_pre_spikes):
        pre_spikes[i] = m['pre_spikes'].spike_trains()[pre_ids[i]]

    # Prepare post spike list
    print(len(m['post_spikes'].spike_trains()), no_post_spikes)

    post_ids = random.sample(range(len(m['post_spikes'].spike_trains())), no_post_spikes)
    post_spikes = [0] * no_post_spikes
    for i in range(no_post_spikes):
        post_spikes[i] = m['post_spikes'].spike_trains()[post_ids[i]]

    if poisson_pop:
        # Prepare poisson spike list
        poisson_ids = random.sample(range(len(m['poisson_spikes'].spike_trains())), no_poisson_spikes)
        poisson_spikes = [0] * no_poisson_spikes
        for i in range(no_poisson_spikes):
            poisson_spikes[i] = m['poisson_spikes'].spike_trains()[poisson_ids[i]]

    # For the legend
    plt.plot(pre_spikes[0], no_spikes * np.ones(len(pre_spikes[0])), '|k')
    plt.plot(post_spikes[0], (no_spikes - no_pre_spikes) * np.ones(len(post_spikes[0])), '|r')
    if poisson_pop:
        plt.plot(poisson_spikes[0], (no_spikes - no_pre_spikes - no_poisson_spikes) * np.ones(len(poisson_spikes[0])),
                 '*b')

    for i in range(len(pre_spikes)):
        plt.plot(pre_spikes[i], (no_spikes - i) * np.ones(len(pre_spikes[i])), '|k', linewidth=3)

    for i in range(len(post_spikes)):
        plt.plot(post_spikes[i], (no_spikes - i - no_pre_spikes) * np.ones(len(post_spikes[i])), '|r', linewidth=3)

    if poisson_pop:
        for i in range(len(poisson_spikes)):
            plt.plot(poisson_spikes[i],
                     (no_spikes - i - no_pre_spikes - no_poisson_spikes) * np.ones(len(poisson_spikes[i])), '*b',
                     linewidth=3)

        # added by GD to monitor spike times
        # print('one train of pre spikes',pre_spikes[0])
        # print('one train of pre spikes',pre_spikes[1])
        # print('one train of pre spikes',pre_spikes[2])
        # print('one train of pre spikes',pre_spikes[3])
        # print('one train of post spikes',post_spikes[0])
        # print('one train of post spikes',post_spikes[1])
        # print('one train of post spikes',post_spikes[2])
        # print('one train of post spikes',post_spikes[3])
        # if poisson_pop:
        # print('one train of poisson spikes',poisson_spikes[0])
        # print('one train of poisson spikes',poisson_spikes[1])
        # print('one train of poisson spikes',poisson_spikes[2])

    # save the result to a text file
    # path = '../GD/fill_the_grid/sLFS/varyingNb/hippoplasticity/added_poisson_input/Clopath_long_tau_refract/'
    # path = '../GD/fill_the_grid/sLFS/varyingNb/visualplast/withpoisson_long_tau_refract/'
    # path = '../GD/puppet/100pulses_08Hz/'

    path = '../data/'
    name = 'Spike trains'
    out_name = path + name
    # save pre spikes
    np.savetxt(out_name + 'pre0.txt', pre_spikes[0])
    np.savetxt(out_name + 'pre1.txt', pre_spikes[1])
    np.savetxt(out_name + 'pre2.txt', pre_spikes[2])
    np.savetxt(out_name + 'pre3.txt', pre_spikes[3])
    # save post spikes
    np.savetxt(out_name + 'post0.txt', post_spikes[0])
    np.savetxt(out_name + 'post1.txt', post_spikes[1])
    np.savetxt(out_name + 'post2.txt', post_spikes[2])
    np.savetxt(out_name + 'post3.txt', post_spikes[3])
    # save poisson spikes
    if poisson_pop:
        np.savetxt(out_name + 'poisson0.txt', poisson_spikes[0])
        np.savetxt(out_name + 'poisson1.txt', poisson_spikes[1])
        np.savetxt(out_name + 'poisson2.txt', poisson_spikes[2])
        np.savetxt(out_name + 'poisson3.txt', poisson_spikes[3])

    plt.legend(['Presynaptic spikes', 'Postsynaptic spikes'], loc='upper right', fontsize=14)
    plt.xlabel('Time (s)', fontsize=20)
    plt.ylabel('Neuron index', fontsize=20)  # added by GD
    plt.savefig('../fig/spike_train.pdf', format='PDF')

    if not title is None:
        plt.title(title, fontsize=20)
    plt.tight_layout()
    if show:
        plt.show()
    if save_in_path and path is not None:
        plt.savefig(path + "spike_plot.pdf", format="PDF")


def plot_mean_weight(m, filename=None, path=None, title=None, show=False):
    """
    This function plots the mean synaptic weight over time.
    :param m:
    :param filename:
    :param path:
    :param title:
    :param show:
    :return:
    """

    plt.plot(m['syn_monitor'].t, np.mean(m['syn_monitor'].w_ampa, 0))
    plt.ylabel('mean weight', fontsize=18)
    plt.xlabel('time $[s]$', fontsize=18)

    if title is not None:
        plt.title(title, fontsize=18)
    plt.tight_layout()
    if show:
        plt.show()

    if filename != None and path != None:
        try:
            plt.savefig(path + filename + '.pdf', format='PDF')
        except:
            print("Wrong combination of path and filename - image not saved")


def visualize_connectivity(S, save_in_path=False):  # S is a Synapses object
    '''
    Plot a two different graps to visualize connectivity between each pre- and postsynaptic neuron.
    '''

    Ns = len(S.source)  # number of source neurons
    Nt = len(S.target)  # number of target neurons
    plt.figure(figsize=(10, 4))
    plt.subplot(121)
    plt.plot(np.zeros(Ns), np.arange(Ns), 'ok', ms=10)  # draw source neurons
    plt.plot(np.ones(Nt), np.arange(Nt), 'ok', ms=10)  # draw target neurons
    for i, j in zip(S.i, S.j):
        plt.plot([0, 1], [i, j], '-k')  # draw lines between them
        plt.xticks([0, 1], ['Source', 'Target'])
    plt.ylabel('Neuron index')
    plt.xlim(-0.1, 1.1)
    plt.ylim(-1, max(Ns, Nt))
    plt.title('Connectivity: Source to Target')

    print(S.i)
    for n in range(max(S.j)):
        indizes = np.where(S.j == n)
        print('Synapses at post neuron {} - {}'.format(n, len(indizes[0])))
    plt.subplot(122)
    plt.plot(S.i, S.j, 'ok')
    plt.xlim(-1, Ns)
    plt.ylim(-1, Nt)
    plt.xlabel('Source neuron index')
    plt.ylabel('Target neuron index')
    plt.title('Connectivity Matrix')
    if save_in_path:
        plt.savefig('../fig/pre:{}_post:{}_connect.pdf'.format(max(S.i), max(S.j)), format='PDF')
    plt.show()
