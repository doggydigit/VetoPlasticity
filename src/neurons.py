#!/usr/bin/env python

"""
    File name: neurons.py
    Author: Halla Sigurthorsdottir, Christoph Blattgerste, Giorgia Dellaferrera, Matthias Tsai
    Date created: 27/11/2018
    Date last modified: ...
    Python Version: 2.7
"""

from brian2 import *

# (CB) Defaults thought for STDP, triplet, quadruplet
default_pre_neuron_parameters = {'model': 'spikegen', 'nr': 1, 'connect_prob': 1}
default_post_neuron_parameters = {'model': 'spikegen', 'nr': 1}

# (CB) Default Ziegler parameters
ziegler_pre_parameters = {'model': 'spikegen',
                          'nr': 500,
                          'connect_prob': 0.1,
                          }

ziegler_post_parameters = {'model': 'LIF',
                           'nr': 10,
                           'V_rest': -70. * mV,
                           'V_thresh': -55. * mV,
                           'V_reset': -80. * mV,
                           'taux': 15. * ms,
                           'gleak': 30. * nS,
                           'C': 300. * pF,
                           'tau_AMPA': 2. * ms,
                           'E_AMPA': 0. * mV,
                           'tau_refract': 50. * ms,  # added by GD   50.ms
                           'refract_0': 5. * ms,  # added by GD  5ms
                           'tau_thr': 200. * ms,  # added by GD   200 ms
                           'V_thr_0': -55. * mV  # added by GD
                           }

# (MT) Standard LIF model with standard set of parameters
LIF_post_parameters = {'model': 'LIF',
                       'nr': 10,
                       'V_rest': -70. * mV,
                       'V_reset': -65. * mV,
                       'V_thresh': -50. * mV,
                       # 'tau_m': 55. * ms,
                       'gleak': 30. * nS,
                       'C': 300. * pF,
                       'tau_AMPA': 2. * ms,
                       'tau_refract': 50. * ms,  # necessary to account for variable in ziegler parameter set
                       'refract_0': 5. * ms,  # necessary to account for variable in ziegler parameter set
                       'tau_thr': 200. * ms,  # added by GD   200 ms
                       'E_AMPA': 0. * mV,
                       }

# # (MT) Standard exIF model (with Adex) with set of parameters from Gerstner
# exIF_post_parameters = {'model': 'exIF',
#                         'nr': 100,
#                         'dt': 0.001 * ms,  # (MT) Integration time step
#                         'V_rest': -70. * mV,  # (MT) Resting potential
#                         'V_reset_Adex': -51. * mV,  # (MT) Voltage to which potential is reset after spike for Adex
#                         'V_reset_exIF': -65. * mV,  # (MT) Voltage to which potential is reset after spike for no Adex
#                         'V_thresh': 40. * mV,  # (MT) Spiking threshold potential at which voltage is reset
#                         'V_rh': -50. * mV,  # (MT) Rheobase (threshold voltage from which exponential term takes over)
#                         # 'tau_m': 5. * ms,  # (MT) Membrane time scale
#                         'gleak': 1. * nS,  # (MT) Leak conductance
#                         'C': 10. * pF,  # (MT) Membrane capacitance
#                         'delta_T': 2. * mV,  # (MT) Slope factor (also called sharpness parameter)
#                         'tau_AMPA': 2. * ms,
#                         'E_AMPA': 0. * mV,
#
#                         # (MT) Adex parameters:
#                         'ad_w': 0.5 * nS,  # (MT) weight adaption voltage coupling
#                         'b_w': 7. * pamp,  # (MT) current nodge to adaptation parameter at each postsynaptic spike
#                         'tau_w': 100. * ms,  # (MT) weight adaption time constant
#                         }

# (MT) Standard eIF model (with Adex) with set of parameters from Brian2
exIF_post_parameters = {'model': 'exIF',
                        'nr': 100,
                        'dt': 0.01 * ms,  # (MT) Integration time step
                        'V_rest': -70. * mV,  # (MT) Resting potential
                        'V_reset_Adex': -65. * mV,  # (MT) Voltage to which potential is reset after spike for Adex
                        'V_reset_exIF': -65. * mV,  # (MT) Voltage to which potential is reset after spike for no Adex
                        'V_thresh': 40. * mV,  # (MT) Spiking threshold potential at which voltage is reset
                        'V_rh': -50. * mV,  # (MT) Rheobase (threshold voltage from which exponential term takes over)
                        # 'tau_m': 5. * ms,  # (MT) Membrane time scale
                        'gleak': 30. * nS,  # (MT) Leak conductance
                        'C': 1. * nF,  # (MT) Membrane capacitance
                        'delta_T': 2. * mV,  # (MT) Slope factor (also called sharpness parameter)
                        'tau_AMPA': 2. * ms,
                        'E_AMPA': 0. * mV,

                        # (MT) Adex parameters that I chose myself to have a visually pleasing spike trace:
                        'ad_w': 10 * nS,  # (MT) weight adaption voltage coupling
                        'b_w': 5. * namp,  # (MT) current nodge to adaptation parameter at each postsynaptic spike
                        'tau_w': 2 * ms,  # (MT) weight adaption time constant
                        }

# # (MT) HH model
# start_scope()
# # Parameters
# area = 20000*umetre**2
# Cm = 1*ufarad*cm**-2 * area
# gl = 5e-5*siemens*cm**-2 * area
# El = -65*mV
# EK = -90*mV
# ENa = 50*mV
# g_na = 100*msiemens*cm**-2 * area
# g_kd = 30*msiemens*cm**-2 * area
# VT = -63*mV
# # The model
# eqs_HH = '''
# dv/dt = (gl*(El-v) - g_na*(m*m*m)*h*(v-ENa) - g_kd*(n*n*n*n)*(v-EK) + I)/Cm : volt
# dm/dt = 0.32*(mV**-1)*(13.*mV-v+VT)/
#     (exp((13.*mV-v+VT)/(4.*mV))-1.)/ms*(1-m)-0.28*(mV**-1)*(v-VT-40.*mV)/
#     (exp((v-VT-40.*mV)/(5.*mV))-1.)/ms*m : 1
# dn/dt = 0.032*(mV**-1)*(15.*mV-v+VT)/
#     (exp((15.*mV-v+VT)/(5.*mV))-1.)/ms*(1.-n)-.5*exp((10.*mV-v+VT)/(40.*mV))/ms*n : 1
# dh/dt = 0.128*exp((17.*mV-v+VT)/(18.*mV))/ms*(1.-h)-4./(1+exp((40.*mV-v+VT)/(5.*mV)))/ms*h : 1
# I : amp
# '''
# group = NeuronGroup(1, eqs_HH,
#                     threshold='v > -40*mV',
#                     refractory='v > -40*mV',
#                     method='exponential_euler')
# group.v = El
# statemon = StateMonitor(group, 'v', record=True)
# spikemon = SpikeMonitor(group, variables='v')
# figure(figsize=(9, 4))
# for l in range(5):
#     group.I = rand()*50*nA
#     run(10*ms)
#     axvline(l*10, ls='--', c='k')
# axhline(El/mV, ls='-', c='lightgray', lw=3)
# plot(statemon.t/ms, statemon.v[0]/mV, '-b')
# plot(spikemon.t/ms, spikemon.v/mV, 'ob')
# xlabel('Time (ms)')
# ylabel('v (mV)');


# (CB) Parameters from Clopath example
Clopath_LIF_post_parameters = {'model': 'LIF',
                               'nr': 10,
                               'V_rest': -70. * mV,
                               'V_reset': -80. * mV,
                               'V_thresh': -55. * mV,
                               'taux': 15. * ms,
                               'gleak': 30. * nS,
                               'C': 300. * pF,
                               'tau_AMPA': 2. * ms,
                               'tau_refract': 50. * ms,  # necessary to account for variable in ziegler parameter set
                               'refract_0': 5. * ms,  # necessary to account for variable in ziegler parameter set
                               'tau_thr': 200. * ms,  # added by GD   200 ms
                               'V_thr_0': -55. * mV,  # added by GD
                               'E_AMPA': 0. * mV,
                               'I_ampl': 1. * mA
                               # The input amplitude required for a single current spike to cause a single neuron spike.
                               # This is only required if the neuron is stimulated directly by protocol.
                               }

# (CB) Parameters from Clopath example code
Clopath_LIF_pre_parameters = {'model': 'LIF',
                              'connect_prob': 1,
                              'nr': 1,
                              'V_rest': -70. * mV,
                              'V_reset': -80. * mV,
                              'V_thresh': -55. * mV,
                              'taux': 15. * ms,
                              'gleak': 30. * nS,
                              'C': 300. * pF,
                              'tau_AMPA': 2. * ms,
                              'E_AMPA': 0. * mV,
                              'I_ampl': 1. * mA
                              # The input amplitude required for a single current spike to cause a single neuron spike.
                              # This is only required if the neuron is stimulated directly by protocol.
                              }

# (CB) parameters from Clopath et al. 2010 Table 1+2 (Fig. 5)
Clopath_exIF_post_parameters = {'model': 'exIF',
                                'nr': 10,
                                'V_rest': -70.6 * mV,
                                'V_reset': -30.4 * mV,  # threshold potential after a spike V_T_max
                                # 'V_thresh': -50.4*mV, # = V_thr_0 (initial condition for variable)
                                'taux': 9.6 * ms,
                                'gleak': 30. * nS,
                                'C': 281. * pF,
                                'delta_T': 2 * mV,
                                'tau_AMPA': 2. * ms,
                                'E_AMPA': 0. * mV,
                                'I_ampl': 1. * mA,

                                # 'theta_reset': 40*mV,   given in Brian2 AdEx Tutorial; not necessary
                                'ad_w': 4 * nS,  # weight adaption
                                'b_w': 0.805 * pamp,  # weight adaption
                                'tau_w': 144 * ms,  # weight adaption
                                'V_thr_0': -50.4 * mV,  # threshold potential at rest = V_T_rest
                                'tau_thr': 50 * ms,  # threshold adaption
                                }

# (CB) parameters from Clopath et al. 2010 Table 1+2 (Fig. 5)
Clopath_exIF_pre_parameters = {'model': 'exIF',
                               'connect_prob': 1,
                               'nr': 1,
                               'V_rest': -70.6 * mV,
                               'V_reset': -50.4 * mV,
                               'V_thresh': 29.4 * mV,  # V_rest + 100
                               'taux': 9.6 * ms,
                               'gleak': 30. * nS,
                               'C': 281. * pF,
                               'delta_T': 2 * mV,
                               'tau_AMPA': 2. * ms,
                               'E_AMPA': 0. * mV,
                               'I_ampl': 1. * mA,

                               # no need for adaption variables for pre neurons
                               }

# (CB) Parameters from Clopath example code
Clopath_Poisson_pre_parameters = {'model': 'LIF',
                                  'nr': 250,
                                  'connect_prob': 0.1,  # (GD) before it was 1
                                  'V_rest': -70. * mV,
                                  'V_thresh': -55. * mV,
                                  'taux': 15. * ms,
                                  'gleak': 30. * nS,
                                  'C': 300. * pF,
                                  'input_rates': 1.8 * Hz  # maybe this has to be changed to 2 Hz
                                  }
