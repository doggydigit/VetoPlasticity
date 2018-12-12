#!/usr/bin/env python

"""
    File name: plasticity.py
    Author: Halla Sigurthorsdottir, Christoph Blattgerste, Giorgia Dellaferrera, Matthias Tsai
    Email: matthias.chinyen.tsai@gmail.com
    Date created: 27/11/2018
    Date last modified: ...
    Python Version: 2.7
"""

from brian2 import *

# Visual cortex parameters from Pfister and Gerstner 2006 (Triplet rule)
Visual_plasticity_parameters = {'PlasticityRule': 'Triplet',
                                'v_increase_init': 2. * mV,  # The initial increase in postsynaptic voltage induced by
                                # one presynaptic spike (initial).
                                'tau_plus': 16.8 * ms,
                                'tau_minus': 33.7 * ms,
                                'tau_x': 101. * ms,
                                'tau_y': 125. * ms,
                                'A_LTP_2': 5e-10,
                                'A_LTD_2': 7e-3,
                                'A_LTP_3': 6.2e-3,
                                'A_LTD_3': 2.3e-4,
                                'init_weight': 0.5,
                                'ampa_max_cond': 5.e-8 * siemens,
                                'w_puppet': 0.125
                                }

# Hippocampal parameters from Pfister and Gerstner 2006 (Triplet rule)
Hippo_plasticity_parameters = {'PlasticityRule': 'Triplet',
                               'v_increase_init': 2. * mV,  # The initial increase in postsynaptic voltage induced by
                               # one presynaptic spike (initial).
                               'tau_plus': 16.8 * ms,
                               'tau_minus': 33.7 * ms,
                               'tau_x': 946. * ms,
                               'tau_y': 27. * ms,
                               'A_LTP_2': 6.1e-3,
                               'A_LTD_2': 1.6e-3,
                               'A_LTP_3': 6.7e-3,
                               'A_LTD_3': 1.4e-3,
                               'init_weight': 0.5,
                               'ampa_max_cond': 5.e-8 * siemens,
                               'w_puppet': 0.125
                               }

# From example code/Clopath 2010:
# REMARK:   The following parameters are not the ones reported in Clopath's article, rather they have been modified
#           in order to exploits Clopath's rule with a LIF neuron (instead of using AdEX neuron as in the paper)
Clopath_LIF_parameters = {'PlasticityRule': 'Clopath',
                          'v_increase_init': 2. * mV,
                          'Theta_low': -70.6 * mV,  # depolarization threshold for plasticity
                          'Theta_high': -62.7 * mV,
                          'x_reset': 1. * ms,  # spike trace reset value'
                          'A_LTD': 1.5e-4,  # depression amplitude
                          'A_LTP': 1.5e-2,  # potentiation amplitude
                          'tau_lowpass1': 40 * ms,  # timeconstant for low-pass filtered voltage
                          'tau_lowpass2': 30 * ms,  # timeconstant for low-pass filtered voltage
                          'tau_m': 5. * ms,
                          'tau_homeo': 1000 * ms,  # homeostatic timeconstant
                          'ampa_max_cond': 5.e-8 * siemens,  # Ampa maximal conductance
                          'w_max': 1.,
                          'init_weight': 0.5,  # initial synaptic weight
                          }

Clopath_exIF_parameters = {'PlasticityRule': 'Clopath',
                           'v_increase_init': 2. * mV,
                           'Theta_low': -55. * mV,  # depolarization threshold for plasticity
                           'Theta_high': -45. * mV,
                           'x_reset': 1. * ms,  # spike trace reset value'
                           'A_LTD': 2.7e-4,  # depression amplitude
                           'A_LTP': 1.2e-4,  # potentiation amplitude
                           'tau_m': 5. * ms,
                           'tau_lowpass1': 10.5 * ms,  # = tau_minus
                           'tau_lowpass2': 200 * ms,  # = tau_plus
                           'tau_homeo': 1000 * ms,  # homeostatic time constant
                           'ampa_max_cond': 5.e-8 * siemens,  # Ampa maximal conductance
                           'w_max': 1.,  # (MT) Maximal hard bounded value of the synaptic weights
                           'init_weight': 0.35,  # 0.0043,  # (MT) Initial synaptic weight values; before this fraction
                           # is set by the function calibrate_w_init(). For 2 mV set to 0.0043
                           'init_stimulation_fraction': 0.255  # (MT) Initial fraction of presynaptic neurons triggered
                           # to spike by a stimulation pulse; before this fraction is set by calibrate_amplitude().
                           }

Claire_LIF_parameters = {'PlasticityRule': 'Claire',
                         'v_increase_init': 2. * mV,
                         'b_theta': 1077 * ms,
                         'Theta_low': -55 * mV,  # depolarization threshold for plasticity
                         'Theta_high': -45 * mV,
                         'x_reset': 1.,  # spike trace reset value'
                         'A_LTD': 0.1,  # depression amplitude
                         'A_LTP': 0.01795,  # potentiation amplitude
                         'tau_lowpass1': 52.63 * ms,  # = tau_minus
                         'tau_lowpass2': 3.04 * ms,  # = tau_plus
                         'tau_theta': 114.9 * ms,
                         'tau_x': 4.877 * ms,
                         'ampa_max_cond': 5.e-8 * siemens,  # Ampa maximal conductance
                         'w_max': 1.,
                         'init_weight': 0.5,  # initial synaptic weight
                         }

Claire_exIF_parameters = {'PlasticityRule': 'Claire',
                          'v_increase_init': 2. * mV,
                          'b_theta': 1077 * ms,
                          'Theta_low': -55 * mV,  # depolarization threshold for plasticity
                          'Theta_high': -45 * mV,
                          'x_reset': 1.,  # spike trace reset value'
                          'A_LTD': 0.1,  # depression amplitude
                          'A_LTP': 0.01795,  # potentiation amplitude
                          'tau_lowpass1': 52.63 * ms,  # = tau_minus
                          'tau_lowpass2': 3.04 * ms,  # = tau_plus
                          'tau_theta': 114.9 * ms,
                          'tau_x': 4.877 * ms,
                          'ampa_max_cond': 5.e-8 * siemens,  # Ampa maximal conductance
                          'w_max': 1.,
                          'init_weight': 0.7,  # initial synaptic weight (good value is 0.348)
                          'init_stimulation_fraction': 0.223
                          }

default_consolidationParams = {'tau_w_ampa': 3600 * second,
                               'k_w': 1.,
                               'w0': 1.,
                               'c_w': 0.5,
                               'tau_z': 36000 * second,
                               'k_z': 1.,
                               'z0': 1.,
                               'c_z': 0.5
                               }


def get_triplet(plasticity_parameters):
    params = {'tau_plus': plasticity_parameters['tau_plus'],
              'tau_minus': plasticity_parameters['tau_minus'],
              'tau_x': plasticity_parameters['tau_x'],
              'tau_y': plasticity_parameters['tau_y'],
              'A_LTP_2': plasticity_parameters['A_LTP_2'],
              'A_LTD_2': plasticity_parameters['A_LTD_2'],
              'A_LTP_3': plasticity_parameters['A_LTP_3'],
              'A_LTD_3': plasticity_parameters['A_LTD_3'],
              'ampa_max_cond': plasticity_parameters['ampa_max_cond']}

    # equations executed at every timestep
    syn_eqs = '''dr_1/dt = -r_1/tau_plus  : 1	(clock-driven)\n'''
    syn_eqs += '''dr_2/dt = -r_2/tau_x 	 : 1	(clock-driven)\n'''
    syn_eqs += '''do_1/dt = -o_1/tau_minus : 1	(clock-driven)\n'''
    syn_eqs += '''do_2/dt = -o_2/tau_y 	 : 1	(clock-driven)\n'''
    syn_eqs += '''w_ampa : 1  # synaptic weight (ampa synapse)'''

    # equations executed only when a presynaptic spike occurs
    # g_ampa += w_ampa*ampa_max_cond
    pre_eqs = '''r_1 += 1\n'''
    pre_eqs += '''w_ampa = w_ampa - o_1*(A_LTD_2 + A_LTD_3*r_2)  # weight update\n'''
    pre_eqs += '''r_2 += 1  # updated after the weight because it is taken in the weight change before the update\n'''

    # equations executed only when a postsynaptic spike occurs
    post_eqs = '''o_1 += 1\n'''
    post_eqs += '''w_ampa = w_ampa + r_1*(A_LTP_2 + A_LTP_3*o_2)  # weight update\n'''
    post_eqs += '''o_2 += 1  # updated after the weight because it is taken in the weight change before the update'''

    return params, pre_eqs, syn_eqs, post_eqs


def get_clopath(plasticity_parameters):
    params = {'ampa_max_cond': plasticity_parameters['ampa_max_cond'],
              'A_LTP': plasticity_parameters['A_LTP'],
              'A_LTD': plasticity_parameters['A_LTD'],
              'Theta_low': plasticity_parameters['Theta_low'],
              'Theta_high': plasticity_parameters['Theta_high'],
              'tau_x': plasticity_parameters['tau_m'],
              'w_max': plasticity_parameters['w_max'],
              'x_reset': plasticity_parameters['x_reset'],
              'step': defaultclock.dt}

    # equations executed at every time step
    syn_eqs = '''dpre_x_trace/dt = -pre_x_trace/tau_x : 1 (clock-driven) # presynaptic spike\n'''  # (MT)THIS IS INEFFICIENT BECAUSE IT COMPUTES PREXTRACE FOR EACH SYNAPSE AND NOT EACH PRESYNAPTIC NEURONX
    syn_eqs += '''dw_ampa/dt = A_LTP * pre_x_trace * (v/mV - Theta_high/mV) * (v_lowpass2_post/mV - Theta_low/mV)'''
    syn_eqs += ''' * int(v/mV - Theta_high/mV > 0) * int(v_lowpass2_post/mV - Theta_low/mV > 0) /ms'''
    syn_eqs += ''' * int(w_max - w_ampa > 0) + int(w_max - w_ampa < 0) * (w_max - w_ampa) / step: 1 (clock-driven)\n'''

    # equations executed only when a pre-synaptic spike occurs
    pre_eqs = '''g_ampa += w_ampa * ampa_max_cond  # increment synaptic conductance\n'''
    pre_eqs += '''w_ampa = clip(w_ampa - A_LTD * (v_lowpass1/mV - Theta_low/mV) * '''
    pre_eqs += '''int(v_lowpass1_post/mV - Theta_low/mV > 0), 0, w_max)\n'''
    pre_eqs += '''pre_x_trace += x_reset / tau_x  # spike trace\n'''

    # equations executed only when a post-synaptic spike occurs
    post_eqs = ''''''

    return params, pre_eqs, syn_eqs, post_eqs


def get_homeoclopath(plasticity_parameters):
    params = {'ampa_max_cond': plasticity_parameters['ampa_max_cond'],
              'A_LTP': plasticity_parameters['A_LTP'],
              'A_LTD': plasticity_parameters['A_LTD'],
              'Theta_low': plasticity_parameters['Theta_low'],
              'Theta_high': plasticity_parameters['Theta_high'],
              'tau_homeo': plasticity_parameters['tau_homeo'],
              'tau_x': plasticity_parameters['tau_m'],
              'w_max': plasticity_parameters['w_max'],
              'v_target': plasticity_parameters['v_target'],
              'x_reset': plasticity_parameters['x_reset']}

    # equations executed at every time step
    syn_eqs = '''dpre_x_trace/dt = -pre_x_trace/tau_x : 1 (clock-driven) # presynaptic spike\n'''
    syn_eqs += '''dw_ampa/dt = A_LTP * pre_x_trace * (v/mV - Theta_high/mV) * (v_lowpass2/mV - Theta_low/mV) * '''
    syn_eqs += '''int(v/mV - Theta_high/mV > 0) * int(v_lowpass2/mV - Theta_low/mV > 0) /ms : 1 (clock-driven)\n'''

    # equations executed only when a pre-synaptic spike occurs
    pre_eqs = '''g_ampa += w_ampa * ampa_max_cond  # increment synaptic conductance\n'''
    pre_eqs += '''w_ampa = clip(w_ampa - A_LTD * (v_homeo**2 / v_target) * (v_lowpass1/mV - Theta_low/mV) * '''
    pre_eqs += '''int(v_lowpass1/mV - Theta_low/mV > 0), 0, w_max)\n'''
    pre_eqs += '''w_ampa = clip(w_ampa - w_minus, 0, w_max)\n'''
    pre_eqs += '''pre_x_trace += x_reset / (tau_x/ms)  # spike trace\n'''

    # equations executed only when a post-synaptic spike occurs
    post_eqs = '''v_homeo += 0.1*mV'''

    return params, pre_eqs, syn_eqs, post_eqs


def get_claire(plasticity_parameters, veto=True):
    params = {'ampa_max_cond': plasticity_parameters['ampa_max_cond'],
              'A_LTP': plasticity_parameters['A_LTP'],
              'A_LTD': plasticity_parameters['A_LTD'],
              'b_theta': plasticity_parameters['b_theta'],
              'Theta_high': plasticity_parameters['Theta_high'],
              'tau_x': plasticity_parameters['tau_x'],
              'tau_lowpass1': plasticity_parameters['tau_lowpass1'],
              'tau_lowpass2': plasticity_parameters['tau_lowpass2'],
              'tau_theta': plasticity_parameters['tau_theta'],
              'w_max': plasticity_parameters['w_max'],
              'x_reset': plasticity_parameters['x_reset']}

    # equations executed at every time step
    syn_eqs = '''dpre_x_trace/dt = -pre_x_trace/tau_x : 1 (clock-driven) # presynaptic spike\n'''
    syn_eqs += '''wLTD = A_LTD * pre_x_trace * (v_lowpass1 - Theta_low)'''
    syn_eqs += ''' * int(v_lowpass1/mV - Theta_low/mV > 0) : volt\n'''
    syn_eqs += '''wLTP = A_LTP * pre_x_trace * (v_lowpass2 - Theta_high)'''
    syn_eqs += ''' * int(v_lowpass2/mV - Theta_high/mV > 0) : volt\n'''
    syn_eqs += '''dtheta/dt = (wLTP - theta) / tau_lowpass1 : volt (clock-driven)\n'''
    syn_eqs += '''dw_ampa/dt = (wLTP - wLTD)/(mV*ms) : 1 (clock-driven)\n'''

    # Add veto equations if required
    if veto:
        params['Theta_low_zero'] = plasticity_parameters['Theta_low']
        syn_eqs += '''Theta_low = Theta_low_zero + b_theta * (wLTP - theta) / tau_theta: volt\n'''
    else:
        params['Theta_low'] = plasticity_parameters['Theta_low']

    # equations executed only when a pre-synaptic spike occurs
    pre_eqs = '''g_ampa += w_ampa * ampa_max_cond  # increment synaptic conductance\n'''
    pre_eqs += '''pre_x_trace += x_reset / (tau_x/ms)  # spike trace\n'''

    # equations executed only when a post-synaptic spike occurs
    post_eqs = ''''''

    return params, pre_eqs, syn_eqs, post_eqs
