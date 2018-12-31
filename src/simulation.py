#!/usr/bin/env python

"""
    File name: simulation.py
    Author: Halla Sigurthorsdottir, Christoph Blattgerste, Giorgia Dellaferrera, Matthias Tsai
    Email: matthias.chinyen.tsai@gmail.com
    Date created: 06/10/2016
    Date last modified: 20/12/2018
    Python Version: 2.7
"""

import copy as cp
import random as rnd
from neurons import *
from plasticity import *


# Change output of run to dictionary
default_protocol_parameters = {'stimulation_type': 'intracellular',
                               'protocol_type': 'STDP',
                               'hro': 1. * Hz,
                               'window': 5. * ms,  # Positive if pre spikes first
                               'nr_rep': 60,
                               }

# From Ziegler 2015
sLFS_protocol_parameters = {'stimulation_type': 'extracellular',
                            'fraction': 0.25,  # of postsynaptic that spike
                            'stimulation_amplitude': 1,
                            'protocol_type': 'sLFS',
                            'hro': 20. * Hz,  # stimulation frequency within blocks
                            'nr_pulses': 3,  # nr of pulses per block
                            'window': 1. * second,  # orginal: 1
                            'nr_blocks': 45,  # original: 900
                            'std': 1. * ms,  # standard deviation of the presynaptic jitter
                            }

# From Ziegler 2015
wLFS_protocol_parameters = {'stimulation_type': 'extracellular',
                            'fraction': 0.25,
                            'stimulation_amplitude': 1,
                            'protocol_type': 'wLFS',
                            'hro': 1. * Hz,  # Arbitrary for this stimulation type (nr_pulses = 1)
                            'nr_pulses': 45,  # original: 900
                            'window': 1. * second,  # Time between blocks. original: 1
                            'nr_blocks': 1,  # original: 1
                            'std': 1. * ms,  # standard deviation of the presynaptic jitter
                            }

# From Ziegler 2015
wTET_protocol_parameters = {'stimulation_type': 'extracellular',
                            'fraction': 0.25,
                            'stimulation_amplitude': 1,
                            'protocol_type': 'wTET',
                            'hro': 100. * Hz,  # stimulation frequency within blocks
                            'nr_pulses': 21,  # original: 21
                            'window': 0.01 * second,  # Arbitrary for this stimulation type (nr_blocks = 1)
                            'nr_blocks': 1,  # original: 1
                            'std': 1. * ms,  # standard deviation of the presynaptic jitter
                            }

# From Ziegler 2015
sTET_protocol_parameters = {'stimulation_type': 'extracellular',
                            'fraction': 0.25,
                            'stimulation_amplitude': 1,
                            'protocol_type': 'sTET',
                            'hro': 100. * Hz,  # stimulation frequency within blocks
                            'nr_pulses': 100,
                            'window': 3. * second,  # original: 10 min
                            'nr_blocks': 3,  # original: 3
                            'std': 1. * ms,  # standard deviation of the presynaptic jitter
                            # Set the clock in case spikes appear more often than one every 5ms .
                            # The default clock has a default time step of 0.1 ms.
                            }

# GD
testing_protocol_parameters = {'stimulation_type': 'extracellular',
                               'fraction': 0.25,
                               'stimulation_amplitude': 1,
                               'protocol_type': 'testing_protocol',
                               'hro': 100. * Hz,  # stimulation frequency within blocks
                               'nr_pulses': 10,
                               'window': 0.1 * second,
                               'nr_blocks': 1,
                               'std': 1. * ms,  # standard deviation of the presynaptic jitter
                               # Set the clock in case spikes appear more often than one every 5ms.
                               # The default clock has a default time step of 0.1 ms.
                               }


class PlasticityProtocol:
    def __init__(self, pre_neuron_parameters=default_pre_neuron_parameters, post_neuron_parameters=LIF_post_parameters,
                 protocol_parameters=default_protocol_parameters, plasticity_parameters=Hippo_plasticity_parameters,
                 time_around=0.1*second, same_connectivity_matrix=False, adaptation=False, veto=False, homeo=False,
                 integration_method='euler', debug=False):
        """
        A synaptic weight change simulation class, which calculates the weight change, when the stimulation protocol
        defined by protocol_parameters is executed on neurons/populations defined by *_neuron_parameters using the
        plasticity rule defined by plasticity_parameters.

        :param pre_neuron_parameters: parameters for presynaptic neurons (see neurons.py for examples)
        :param post_neuron_parameters: parameters for postsynaptic neurons (see neurons.py for examples)
        :param protocol_parameters: parameters of the stimulation protocol (see above for examples)
        :param plasticity_parameters: parameters of the pasticity rule (see plasticity.py for examples)
        :param time_around: time after the end of the protocol
        :param same_connectivity_matrix: should be set to True if for every run the same synapses are to be created,
        to False otherwise. Note that if True a rnd.seed is set and the weights have always the same shape
        :param adaptation: True in case adaptation parameter with corresponding differential equations desired.
        :param veto: whether or not to use veto mechanism in Clopath model
        :param homeo: whether or not to use homeostatic term in Clopath model
        :param integration_method: method to use to integrate differential equations with brian
        :param debug: Set to true for verbose output and shorter simulation
        """

        if debug:
            print('\n\n\nInitialization:')

            import os
            import psutil
            self.process = psutil.Process(os.getpid())
            print(self.process.memory_info().rss)

        if veto and homeo:
            raise NotImplementedError('homeo and veto cannot be set to True at the same time')

        if same_connectivity_matrix:
            seed(1)
            rnd.seed(321)
            np.random.seed(1002)

        if 'dt' in post_neuron_parameters:
            defaultclock.dt = post_neuron_parameters['dt']
            self.timestep = post_neuron_parameters['dt']
        else:
            self.timestep = defaultclock.dt

        self.protocol = protocol_parameters['protocol_type']
        self.plasticity_rule = plasticity_parameters['PlasticityRule']
        self.pre_neuron_parameters = pre_neuron_parameters
        self.post_neuron_parameters = post_neuron_parameters
        self.protocol_parameters = protocol_parameters
        self.plasticity_parameters = plasticity_parameters
        self.time_around = time_around
        self.init_weight = plasticity_parameters['init_weight']
        self.adaptation = adaptation
        self.homeo = homeo
        self.veto = veto
        self.method = integration_method
        self.debug = debug

        # The input to the spike generator group was misunderstood -> rethink whole protocol pipeline
        self.t_pre, self.pre_neuron_parameters['ids'], final_t_pre = self.make_t('pre')
        self.final_t = self.time_around + final_t_pre * second

        self.neuron_pre = SpikeGeneratorGroup(pre_neuron_parameters['nr'],
                                              pre_neuron_parameters['ids'],
                                              self.t_pre * second)
        self.neuron_post = self.make_neuron()
        self.syn = self.make_synapse()

        print('Initialized Simulation with:')
        print('- Simulation time: {}'.format(self.final_t))
        print("- Neuron model: {}".format(post_neuron_parameters["model"]))
        print("- Stimulation protocol: {}".format(protocol_parameters['protocol_type']))
        print("- Parameters {}".format(protocol_parameters))
        print("- Plasticity model: {}\n".format(self.plasticity_rule))

    def run(self, syn_parameters=False, pre_parameters=False, post_parameters=False, pre_spikes=False,
            post_spikes=False):
        """
        :param syn_parameters: a list of the parameters that should be recorded from the synapse, e.g. ['w_ampa', 'r_1']
        :param pre_parameters: a list of the parameters that should be recorded from the presynaptic neuron.
        :param post_parameters: a list of the parameters that should be recorded from the postsynaptic neuron.
        :param pre_spikes: True or False, whether the spike times of the presynaptic neuron are desired as an output.
        :param post_spikes: True or False, whether the spike times of the postsynaptic neuron are desired as an output.
        :return: parameters_out, a dictionary containing all the parameters that were recorded, with keys syn_monitor,
        pre_monitor, post_monitor, pre_spikes, post_spikes (the ones that were asked).
        """

        print('\nRunning Simulation:')

        init_weight = self.plasticity_parameters['init_weight']
        neuron_pre = self.neuron_pre
        neuron_post = self.neuron_post
        syn = self.syn
        syn.w_ampa = init_weight
        if self.debug:
            print('before initialization, w.shape = ' + str(syn.w_ampa.shape))
            print(self.process.memory_info().rss)

        # Define monitors
        parameters_out = {}

        if syn_parameters:
            # (GD) StateMonitor and SpikeMonitor are Brian2 functions that record events and spikes respectively from
            # a NeuronGroup or another event source
            syn_monitor = StateMonitor(syn, syn_parameters, record=True, dt=1*ms)
            parameters_out['syn_monitor'] = syn_monitor

        if pre_spikes:
            pre_spikes = SpikeMonitor(neuron_pre)
            parameters_out['pre_spikes'] = pre_spikes

        if post_spikes:
            post_spikes = SpikeMonitor(neuron_post)
            parameters_out['post_spikes'] = post_spikes

        if (not self.pre_neuron_parameters['model'] == 'spikegen') and pre_parameters:
            pre_monitor = StateMonitor(neuron_pre, pre_parameters, record=True, dt=0.01*ms)
            parameters_out['pre_monitor'] = pre_monitor

        if (not self.post_neuron_parameters['model'] == 'spikegen') and post_parameters:
            post_monitor = StateMonitor(neuron_post, post_parameters, record=True, dt=0.01*ms)
            parameters_out['post_monitor'] = post_monitor

        # Run
        # (GD) 'run' function runs a simulation with all the Brian objects for the given duration (self.final_t)
        # (GD) 'report' specifies how to report the progress of the simulation
        run(self.final_t, report='text')

        return parameters_out

    def make_neuron(self):
        """
        This function makes a postsynaptic neuron group based on the postsynaptic neuron parameters of the class object.
        :return: the postsynaptic neuron group
        """

        # (MT) Define that neuron group is built from the postsynaptic neuron parameters of the class object
        neuron_parameters = self.post_neuron_parameters

        # (MT) Name of neuron group
        neurongroupname = 'postynaptic_neuron_group'

        # (MT) Define parameters that will be used in the neuron equations in this dictionary
        neuronparams = {'V_rest': neuron_parameters['V_rest'],
                        'V_thresh': neuron_parameters['V_thresh'],
                        'gleak': neuron_parameters['gleak'],
                        'C': neuron_parameters['C'],
                        'tau_AMPA': neuron_parameters['tau_AMPA'],
                        'E_AMPA': neuron_parameters['E_AMPA']}

        # (MT) Conditions for spiking
        threshold = 'v>V_thresh'

        # (MT) Equations defining what happens when spiking threshold is reached
        reset = 'v = V_reset'

        # (MT) Equations describing synaptic channel conductances
        eqs = '''I_syn = g_ampa*(E_AMPA-v): amp                       # synaptic current\n'''
        eqs += '''dg_ampa/dt = -g_ampa/tau_AMPA : siemens              # synaptic conductance\n'''

        if neuron_parameters['model'] == 'LIF':
            # (MT) Leaky integrate and fire voltage equations
            neuronparams['V_reset'] = neuron_parameters['V_reset']
            eqs += '''dv/dt = (gleak*(V_rest-v) + I_syn)/C: volt   # voltage\n'''

        elif neuron_parameters['model'] == 'exIF':
            # (MT) In case of exIF or Adex
            neuronparams['V_rh'] = neuron_parameters['V_rh']
            neuronparams['delta_T'] = neuron_parameters['delta_T']

            if self.adaptation:
                # (MT) Adex case
                neuronparams['ad_w'] = neuron_parameters['ad_w']
                neuronparams['b_w'] = neuron_parameters['b_w']
                neuronparams['tau_w'] = neuron_parameters['tau_w']
                neuronparams['V_reset'] = neuron_parameters['V_reset_Adex']

                reset += '; w += b_w'
                eqs += '''dv/dt = (gleak*(V_rest-v) + gleak * delta_T * exp((v-V_rh) / delta_T) + '''
                eqs += '''I_syn - w)/C : volt # voltage\n'''
                eqs += '''dw/dt = (ad_w * (v-V_rest) - w) / tau_w : amp # synaptic weight: 1 \n'''
            else:
                # (MT) Classical exIF case
                neuronparams['V_reset'] = neuron_parameters['V_reset_exIF']
                eqs += '''dv/dt = (gleak*(V_rest-v) + gleak * delta_T * exp((v-V_rh) / delta_T) + '''
                eqs += '''I_syn)/C : volt # voltage\n'''
        else:
            raise ValueError

        # (MT) Plasticity equations
        if self.plasticity_parameters['PlasticityRule'] in ['Clopath', 'Claire']:
            neuronparams['tau_lowpass1'] = self.plasticity_parameters['tau_lowpass1']
            neuronparams['tau_lowpass2'] = self.plasticity_parameters['tau_lowpass2']
            eqs += '''dv_lowpass1/dt = (v-v_lowpass1)/tau_lowpass1 : volt  # low-pass filter\n'''
            eqs += '''dv_lowpass2/dt = (v-v_lowpass2)/tau_lowpass2 : volt  # low-pass filter\n'''

            # (MT) In case it is desired to use the homeostatic term of Clopath et al. 2010:
            # "Connectivity reflects coding: a model of voltage-based STDP with homeostasis."
            if self.homeo:
                neuronparams['tau_homeo'] = self.plasticity_parameters['tau_homeo']
                eqs += '''dv_homeo/dt = (v-V_rest-v_homeo)/tau_homeo : volt    # low-pass filter\n'''
        else:
            raise ValueError

        # (MT) Create neuron group object according to the upper defined equations
        neuron_out = NeuronGroup(N=neuron_parameters['nr'], model=eqs, threshold=threshold, reset=reset,
                                 namespace=neuronparams, name=neurongroupname, method=self.method)

        # (MT) Initialize the values of the variables defined by differential equations
        neuron_out.v = neuronparams['V_rest']
        neuron_out.g_ampa = 0
        if self.plasticity_parameters['PlasticityRule'] in ['Clopath', 'Claire']:
            neuron_out.v_lowpass1 = neuronparams['V_rest']
            neuron_out.v_lowpass2 = neuronparams['V_rest']
            if self.homeo:
                neuron_out.v_homeo = 0

        return neuron_out

    def make_t(self, prepost='pre'):
        """
        This function makes the spiking vector corresponding to the stimulation protocol

        :param prepost: 'pre' if the presynaptic neuron is stimulate
        :return:
        """

        if prepost == 'pre':
            neuron_parameters = self.pre_neuron_parameters
        elif prepost == 'post':
            neuron_parameters = self.post_neuron_parameters
        else:
            raise ValueError

        # (GD) in the cases: sTET, wTET, sLFS, wLFS, testing
        if self.protocol_parameters['stimulation_type'] == 'extracellular':
            if prepost == 'pre':
                hro = self.protocol_parameters['hro']
                nr_pulses = self.protocol_parameters['nr_pulses']
                window = self.protocol_parameters['window']
                nr_blocks = self.protocol_parameters['nr_blocks']
                st = self.protocol_parameters['std']

                # (GD) QUESTION time_around is the time after the end of the protocol,
                # passed as input in the init of ConsolidationExperiment:
                # WHY DO WE CREATE AN ARRAY THAT STARTS FROM THE END OF THE PROTOCOL AND CONTINUES FROM THERE?
                t_main = arange(self.time_around / second,
                                self.time_around / second + nr_pulses * 1 / (hro / Hz),
                                1 / (hro / Hz))
                t_orig = t_main

                # (GD) create a long array (one for the whole protocol) in which for every block you repeat the array
                # t_main but shifted by the time lapse at which it occurs wrt the start
                for i in range(nr_blocks - 1):
                    t_main = concatenate((t_main, t_orig + (i + 1) * window / second))

                # To be calibrated:
                fraction_stimulated = self.protocol_parameters['stimulation_amplitude']

                no_firing_neurons = int(self.pre_neuron_parameters['nr'] * fraction_stimulated)
                if self.debug:
                    print('nr presynaptic neurons ' + str(self.pre_neuron_parameters['nr']))
                    print('fraction presynaptic stimulated = ' + str(fraction_stimulated))
                    print('nr presynaptic neurons firing per stimulation pulse = ' + str(no_firing_neurons))

                if no_firing_neurons < 1:
                    print('Warning: No neurons are firing, firing amplitude too low or number of neurons too low')

                # Produce spking times of presynaptic neurons following first extracellular pulse
                t = t_main[0] + np.clip(absolute(np.random.normal(0, st / second, no_firing_neurons)),
                                        a_min=None, a_max=1 / (hro / Hz) - 0.002)

                # Sample which presynaptic neurons to activate from the extracellular pulse
                ids_orig = range(neuron_parameters['nr'])
                ids = rnd.sample(ids_orig, no_firing_neurons)

                # Repeat for all other extracellular pulses
                for i in range(1, len(t_main)):
                    t_next = t_main[i] + np.clip(absolute(np.random.normal(0, st / second, no_firing_neurons)),
                                                 a_min=None, a_max=1 / (hro / Hz) - 0.002)

                    t = concatenate((t, t_next))
                    # (GD) at the end ids contains the sequence of the neurons firing for the whole protocol
                    ids = concatenate((ids, rnd.sample(ids_orig, no_firing_neurons)))

                final_t = t[-1]
            elif prepost == 'post':  # (GD) why if prepost is not 'pre' we leave the list empty?
                t = []
                final_t = 0
                ids = []

            elif prepost == 'poisson':
                # TO BE CHECKED
                t = []
                final_t = 0
                ids = []

            else:
                raise ValueError
        else:
            raise ValueError

        return t, ids, final_t

    def make_synapse(self):
        """
        This function makes a synapse between the two neuron groups based on the plastictiy parameters

        :return:
        """

        # (MT) Get the plasiticity equations depending on the chosen plasticity rule
        if self.plasticity_parameters['PlasticityRule'] == 'Triplet':
            params, pre_eqs, syn_eqs, post_eqs = get_triplet(self.plasticity_parameters)

        elif self.plasticity_parameters['PlasticityRule'] == 'Clopath':
            if self.homeo:
                params, pre_eqs, syn_eqs, post_eqs = get_homeoclopath(self.plasticity_parameters)
            else:
                params, pre_eqs, syn_eqs, post_eqs = get_clopath(self.plasticity_parameters)

        elif self.plasticity_parameters['PlasticityRule'] == 'Claire':
            params, pre_eqs, syn_eqs, post_eqs = get_claire(self.plasticity_parameters, self.veto)

        else:
            raise ValueError

        # (MT) Construct the synapses according to the equations
        syn = Synapses(source=self.neuron_pre, target=self.neuron_post, model=syn_eqs, on_pre=pre_eqs, on_post=post_eqs,
                       multisynaptic_index='synapse_number', namespace=params, method=self.method)

        # (GD) this connects all neuron pairs with given prob
        syn.connect(p=self.pre_neuron_parameters['connect_prob'])
        if self.veto:
            syn.theta = 0

        if self.debug:
            print('synapse made according to ', self.plasticity_parameters['PlasticityRule'])

        return syn

    def calibrate_amplitude(self, iter_max=20, error_margin=0.2, std_cal=0. * ms):
        """
         This function sets the amplitude (if amplitude is True) of stimulation (percentage of presynaptic neurons
         firing) such that the given percentage of postsynaptic neurons ('fraction' in post_neuron_parameters) spikes in
         response to a pulse.
        :param iter_max:
        :param error_margin:
        :param std_cal:
        :return:
        """

        print('\nCalibrating presynaptic stimulation amplitude:')

        if 'init_stimulation_fraction' in self.plasticity_parameters:
            init_frac = self.plasticity_parameters['init_stimulation_fraction']
        else:
            init_frac = 1

        cal_protocol_parameters = {'stimulation_type': 'extracellular',
                                   'protocol_type': 'calibration',
                                   'hro': 20. * Hz,  # does not matter because nr_pulses == 1
                                   'nr_pulses': 1,
                                   'window': 1. * second,  # does not matter because nr_blocks == 1
                                   'nr_blocks': 1,
                                   'stimulation_amplitude': init_frac,
                                   'std': std_cal,
                                   }

        calibration_ex = cp.deepcopy(self)
        # calibration_ex.debug = False
        calibration_ex.protocol_parameters = cal_protocol_parameters
        calibration_ex = redefine_experiment(calibration_ex)

        print('Run with calibration_amplitude = ' + str(calibration_ex.protocol_parameters['stimulation_amplitude']))

        # (GD) with the protocol just set perform the function 'run' of this class.
        # 'm' is the output ('parameters_out') of 'run'
        m = calibration_ex.run(post_spikes=True, pre_spikes=True)

        # Calculate the fraction of spiking postsynaptic neurons vs all postsynaptic neurons
        fraction = sum(m['post_spikes'].count > 0) / (self.post_neuron_parameters['nr'] * 1.0)
        if self.debug:
            print('fraction = ' + str(sum(m['post_spikes'].count > 0)) + '/' + str(self.post_neuron_parameters['nr']))
            print('presynaptic spike fraction = ' + str(sum(m['pre_spikes'].count > 0) * 1.0) + '/' +
                  str(self.pre_neuron_parameters['nr']) + ' = ' +
                  str(sum(m['pre_spikes'].count > 0)*1.0/self.pre_neuron_parameters['nr']))
        count = 1  # (GD) counter to keep trace of how many times the stimulation amplitude has been recalibrated

        # (GD) Check whether the computed fraction is too big or too small (too_big,too_small = True or False)
        too_big = fraction > self.protocol_parameters['fraction'] * (1 + error_margin)
        too_small = fraction < self.protocol_parameters['fraction'] * (1 - error_margin)

        if too_small:
            raise SystemError('Calibration failed, synaptic weights are too low to reach a postsynaptic fraction of ',
                              self.protocol_parameters['fraction'] * 100, '%')
        else:
            min_ampl = 0
            max_ampl = calibration_ex.protocol_parameters['stimulation_amplitude']
            while (too_big or too_small) and count < iter_max:
                if too_big:
                    max_ampl = calibration_ex.protocol_parameters['stimulation_amplitude']
                else:
                    min_ampl = calibration_ex.protocol_parameters['stimulation_amplitude']

                # (GD) Here we change the amplitude since the one previously tried was too big or too small
                calibration_ex.protocol_parameters['stimulation_amplitude'] = (max_ampl + min_ampl) / 2.
                calibration_ex = redefine_experiment(calibration_ex)

                if self.debug:
                    print(max_ampl, min_ampl, calibration_ex.protocol_parameters['stimulation_amplitude'])
                    print('Fraction is ' + str(fraction))
                    print('Running with calibration_amplitude = ',
                          calibration_ex.protocol_parameters['stimulation_amplitude'])

                # (GD) re'run' and check again if the fraction is too big of too small
                m = calibration_ex.run(post_spikes=True, pre_spikes=True)
                fraction = sum(m['post_spikes'].count > 0) / (1.0 * self.post_neuron_parameters['nr'])
                if self.debug:
                    print('fraction = ' + str(sum(m['post_spikes'].count > 0)) + '/' +
                          str(self.post_neuron_parameters['nr']) + ' = ' + str(fraction))
                    print('presynaptic spike fraction = ' + str(sum(m['pre_spikes'].count > 0) * 1.0) + '/' +
                          str(self.pre_neuron_parameters['nr']) + ' = ' +
                          str(sum(m['pre_spikes'].count > 0) * 1.0 / self.pre_neuron_parameters['nr']))

                # Update
                too_big = fraction > self.protocol_parameters['fraction'] * (1 + error_margin)
                too_small = fraction < self.protocol_parameters['fraction'] * (1 - error_margin)
                count += 1

        if count == iter_max:
            raise SystemError('Calibration failed, maximum number of iterations: ' + str(iter_max))

        # (GD) With this new amplitude create the spiking vector, the neuron and the synapses
        self.protocol_parameters['stimulation_amplitude'] = calibration_ex.protocol_parameters['stimulation_amplitude']
        self.t_pre, self.pre_neuron_parameters['ids'], final_t_pre = self.make_t('pre')
        self.neuron_pre = SpikeGeneratorGroup(self.pre_neuron_parameters['nr'],
                                              self.pre_neuron_parameters['ids'],
                                              self.t_pre * second)
        if self.debug:
            print('First spike at {}'.format(min(self.t_pre)))
        self.syn = self.make_synapse()

        print('*******************************************************************************************************')
        print('Stimulation amplitude (fraction) set to ' + str(self.protocol_parameters['stimulation_amplitude']) +
              ', which gives a postsynaptic firing fraction of ' + str(fraction))
        print('*******************************************************************************************************')

    def calibrate_w_init(self, error_margin=0.05, iter_max=20, std_cal=0. * ms):
        """
        This function sets the initial weight (if weight is True) such that each presynaptic spike (if connected) causes
        a postsynaptic voltage increase of plasticity parameters' v_increase_init.

        :param error_margin:
        :param iter_max:
        :param std_cal:
        :return:
        """

        print('\nCalibrating weights:')

        # Set the experiment up as a single connection and presynaptic will spike once
        cal_protocol_parameters = {'stimulation_type': 'extracellular',
                                   'protocol_type': 'calibration',
                                   'hro': 1. * Hz,  # does not matter because nr_pulses == 1
                                   'nr_pulses': 1,
                                   'window': 1. * second,  # does not matter because nr_blocks == 1
                                   'nr_blocks': 1,
                                   'stimulation_amplitude': 1,
                                   'std': std_cal,
                                   }

        calibration_ex = cp.deepcopy(self)
        # calibration_ex.debug = False
        calibration_ex.protocol_parameters = cal_protocol_parameters
        calibration_ex.pre_neuron_parameters['nr'] = 1
        calibration_ex.post_neuron_parameters['nr'] = 1
        calibration_ex.pre_neuron_parameters['connect_prob'] = 1
        calibration_ex = redefine_experiment(calibration_ex)

        print('Initial weights = ' + str(calibration_ex.plasticity_parameters['init_weight']))

        m = calibration_ex.run(post_parameters='v', pre_spikes=True)

        # Calculate the voltage increase
        # print(m['post_monitor'].v[0][m['post_monitor'].t == m['pre_spikes'].t], m['post_monitor'].v[0][0])
        v_increase = max(m['post_monitor'].v[0]) - m['post_monitor'].v[0][0]

        count = 1  # (GD) Counter to keep trace of how many times the weights are recalibrated

        # (GD) Check if the voltage increase is smaller or bigger than the one required
        too_big = v_increase > self.plasticity_parameters['v_increase_init'] * (1 + error_margin)
        too_small = v_increase < self.plasticity_parameters['v_increase_init'] * (1 - error_margin)
        if too_small:
            raise SystemError('Initial weight is too small to reach desired postsynaptic voltage nodge '
                              '-> Increase it.\n The achieved voltage increase was: ' + str(v_increase))

        max_w = calibration_ex.plasticity_parameters['init_weight']  # 4 is handpicked but should work
        min_w = 0
        while (too_big or too_small) and count < iter_max:

            if too_big:
                max_w = calibration_ex.plasticity_parameters['init_weight']
            elif too_small:
                min_w = calibration_ex.plasticity_parameters['init_weight']

            # (GD) QUESTION Why is the update done in such a way?
            calibration_ex.plasticity_parameters['init_weight'] = (max_w + min_w) / 2.

            calibration_ex = redefine_experiment(calibration_ex)

            # (GD) Re'run' and check if now the requirements are met, otherwise repeat
            m = calibration_ex.run(post_parameters='v')
            v_increase = max(m['post_monitor'].v[0]) - m['post_monitor'].v[0][0]

            # Update
            too_big = v_increase > self.plasticity_parameters['v_increase_init'] * (1 + error_margin)
            too_small = v_increase < self.plasticity_parameters['v_increase_init'] * (1 - error_margin)
            count += 1

        if count == iter_max:
            raise SystemError('Weight calibration failed with voltage nodge of: ' + str(v_increase) +
                              ' instead of ' + str(self.plasticity_parameters['v_increase_init']))
        self.plasticity_parameters['init_weight'] = calibration_ex.plasticity_parameters['init_weight']
        print('*******************************************************************************************************')
        print('Initial weights set to ' + str(self.plasticity_parameters['init_weight']) +
              ', which gives a postsynaptic response of ' + str(v_increase))
        print('*******************************************************************************************************')


def redefine_experiment(calibration_ex):
    """
    This is a helper function for calibrate_amplitude and ccalibrate_w_init. It basically reinitializes some of the
    contents of the copied class object.
    :param calibration_ex: copied class object to reinitialize
    :return: reinitialized class object
    """

    calibration_ex.t_pre, calibration_ex.pre_neuron_parameters['ids'], final_t_pre = calibration_ex.make_t('pre')

    calibration_ex.final_t = calibration_ex.time_around + final_t_pre * second

    calibration_ex.neuron_pre = SpikeGeneratorGroup(calibration_ex.pre_neuron_parameters['nr'],
                                                    calibration_ex.pre_neuron_parameters['ids'],
                                                    calibration_ex.t_pre * second)

    calibration_ex.neuron_post = calibration_ex.make_neuron()

    calibration_ex.syn = calibration_ex.make_synapse()

    return calibration_ex
