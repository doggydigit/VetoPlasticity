from brian2 import *
from matplotlib import *
import random as rnd
from numpy import *
import copy
import parameters as p

# Change output of run to dictionary

# Defaults thought for STDP, triplet, quadruplet
default_pre_neuron_parameters = {'model': 'spikegen',
					'nr': 1,
					'connect_prob': 1
}

default_post_neuron_parameters = {'model': 'spikegen',
					'nr': 1
}

ziegler_pre_parameters = {'model': 'spikegen',
					'nr': 100,
					'connect_prob': 0.1  # change back to 0.1
}

ziegler_post_parameters = {'model': 'LIF',
			'nr': 10,
			'V_rest': -70.*mV,
			'V_thresh': -55.*mV,
			'V_reset': -80.*mV,
			'taux': 15.*ms,
			'taux_LTP': 15.*ms,
			'taux_LTD': 15.*ms,
			'gleak': 30.*nS,
			'C': 300.*pF,
			'tau_AMPA': 2.*ms,
			'E_AMPA': 0.*mV,
			'tau_refract' : 50.*ms, #added by GD   50.ms
			'refract_0': 5.*ms, #added by GD  5ms
			'tau_thr': 200.*ms, # added by GD   200 ms
			'V_thr_0': -55.*mV # added by GD
}

# Parameters from Clopath example
LIF_post_parameters = {'model': 'LIF',
			'nr': 1,
			'V_rest': -70.*mV,
			'V_reset': -80.*mV,
			'V_thresh': -55.*mV,
			'taux': 15.*ms,
			'taux_LTP': 15.*ms,
			'taux_LTD': 15.*ms,
			'gleak': 30.*nS,
			'C': 300.*pF,
			'tau_AMPA': 2.*ms,
			'E_AMPA': 0.*mV,
			'I_ampl': 1.*mA, # The input amplitude required for a single current spike to cause a single neuron spike. This is only required if the neuron is stimulated directly by protocol.
			'tau_refract' : 50.*ms, #added by GD
			'refract_0': 5.*ms, # added by GD
			'tau_thr': 200.*ms, # added by GD
			'V_thr_0': -55.*mV # added by GD
}

# Parameters from Clopath example code
LIF_pre_parameters = {'model': 'LIF',
			'connect_prob': 1,
			'nr': 1,
			'V_rest': -70.*mV,
			'V_reset': -80.*mV,
			'V_thresh': -55.*mV,
			'taux': 15.*ms,
			'gleak': 30.*nS,
			'C': 300.*pF,
			'tau_AMPA': 2.*ms,
			'E_AMPA': 0.*mV,
			'I_ampl': 1.*mA # The input amplitude required for a single current spike to cause a single neuron spike. This is only required if the neuron is stimulated directly by protocol.
}

# parameters from Clopath et al. 2010 Table 1+2 (Fig. 5)
exIF_post_parameters = {'model': 'exIF',
			'nr': 10,
			'V_rest': -70.6*mV,
			'V_reset': -30.4*mV, # threshold potential after a spike V_T_max
			# 'V_thresh': -50.4*mV, # = V_thr_0 (initial condition for variable)
			'taux': 9.6*ms,
			'taux_LTP': 15.*ms,
			'taux_LTD': 15.*ms,
			'gleak': 30.*nS,
			'C': 281.*pF,
			'delta_T': 2*mV,
			'tau_AMPA': 2.*ms,
			'E_AMPA': 0.*mV,
			'I_ampl': 1.*mA,

			# 'theta_reset': 40*mV,   given in Brian2 AdEx Tutorial; not necessary
			'ad_w': 4*nS, # weight adaption
			'b_w': 0.805*pamp, # weight adaption
			'tau_w': 144*ms, # weight adaption
			'V_thr_0': -50.4*mV, # threshold potential at rest = V_T_rest
			'tau_thr': 50*ms, # threshold adaption
}

# parameters from Clopath et al. 2010 Table 1+2 (Fig. 5)
exIF_pre_parameters = {'model': 'exIF',
			'connect_prob': 1,
			'nr': 1,
			'V_rest': -70.6*mV,
			'V_reset': -50.4*mV,
			'V_thresh': 29.4*mV, # V_rest + 100
			'taux': 9.6*ms,
			'gleak': 30.*nS,
			'C': 281.*pF,
			'delta_T': 2*mV,
			'tau_AMPA': 2.*ms,
			'E_AMPA': 0.*mV,
			'I_ampl': 1.*mA,
			# no need for adaption variables for pre neurons
}

# Parameters from Clopath example code
poisson_pre_parameters = {'model': 'LIF',
			'nr': 250,
			'connect_prob': 0.1,  #(GD) before it was 1
			'V_rest': -70.*mV,
			'V_thresh': -55.*mV,
			'taux': 15.*ms,
			'gleak': 30.*nS,
			'C': 300.*pF,
			'input_rates': 1.8*Hz # maybe this has to be changed to 2 Hz
}

default_protocol_parameters = {'stimulation_type': 'intracellular',
			'protocol_type': 'STDP',
			'hro': 1.*Hz,
			'dt': 10.*ms, #Positive if pre spikes first
			'nr_rep': 60,
			'clock': defaultclock.dt
}

# From Ziegler 2015
sLFS_protocol_parameters = {'stimulation_type': 'extracellular',
				'fraction': 0.25, # of postsynaptic that spike
				'stimulation_amplitude': 1,
				'protocol_type': 'sLFS',
				'hro': 20.*Hz, # stimulation frequency within blocks
				'nr_pulses': 3,
				'dt': 1.*second,
				'nr_blocks': 100,
				'std': 3.*ms, #standard deviation of the presynaptic jitter
				'clock': defaultclock.dt

}

# From Ziegler 2015

wLFS_protocol_parameters = {'stimulation_type': 'extracellular',
				'fraction': 0.25,
				'stimulation_amplitude': 1,
				'protocol_type': 'wLFS',
				'hro': 1.*Hz, # Arbitrary for this stimulation type (nr_pulses = 1)
				'nr_pulses': 1,
				'dt': 1.*second, # Time between blocks
				'nr_blocks': 100, # original: 900
				'std': 3.*ms, #standard deviation of the presynaptic jitter
				'clock': defaultclock.dt

}

# From Ziegler 2015
wTET_protocol_parameters = {'stimulation_type': 'extracellular',
				'fraction': 0.25,
				'stimulation_amplitude': 1,
				'protocol_type': 'wTET',
				'hro': 100.*Hz, # stimulation frequency within blocks
				'nr_pulses': 21,
				'dt': 60.*60.*second, # Arbitrary for this stimulation type (nr_blocks = 1)
				'nr_blocks': 2,
				'std': 3.*ms, #standard deviation of the presynaptic jitter
				'clock': defaultclock.dt

}

# From Ziegler 2015
sTET_protocol_parameters = {'stimulation_type': 'extracellular',
				'fraction': 0.25,
				'stimulation_amplitude': 1,
				'protocol_type': 'sTET',
				'hro': 100.*Hz, # stimulation frequency within blocks
				'nr_pulses':100,
				'dt': 60*second, # 10 min
				'nr_blocks': 2,
				'std': 3.*ms, #standard deviation of the presynaptic jitter
				'clock': defaultclock.dt # Set the clock in case spikes appear more often than one every 5ms . The default clock has a default time step of 0.1 ms.

}

# GD
testing_protocol_parameters = {'stimulation_type': 'extracellular',
				'fraction': 0.25,
				'stimulation_amplitude': 1,
				'protocol_type': 'testing_protocol',
				'hro': 100.*Hz, # stimulation frequency within blocks
				'nr_pulses':100,
				'dt': 1.*second,
				'nr_blocks': 3,
				'std': 3.*ms, #standard deviation of the presynaptic jitter
				'clock': defaultclock.dt # Set the clock in case spikes appear more often than one every 5ms . The default clock has a default time step of 0.1 ms.

}

# Visual cortex parameters from Pfister and Gerstner 2006 (Triplet rule)
Visual_plasticity_parameters = {'PlasticityRule': 'Triplet',
					'v_increase_init': 2.*mV, # The increase of postsynaptic voltage induced by one spike (initial), this is only used by calibrate if weight is true.
					'tau_plus': 16.8*ms,
					'tau_minus': 33.7*ms,
					'tau_x': 101.*ms,
					'tau_y': 125.*ms,
					'A_LTP_2': 5e-10,
					'A_LTD_2': 7e-3,
					'A_LTP_3': 6.2e-3,
					'A_LTD_3': 2.3e-4,
					'init_weight': 0.5,
					'w_max': 1.,
					'w_min': 0.,
					'ampa_max_cond': 5.e-8*siemens,
					'w_puppet': 0.125
}

# Hippocampal parameters from Pfister and Gerstner 2006 (Triplet rule)
Hippo_plasticity_parameters = {'PlasticityRule': 'Triplet',
					'v_increase_init': 2.*mV, # The increase of postsynaptic voltage induced by one spike (initial), this is only used by calibrate if weight is true.
					'tau_plus': 16.8*ms,
					'tau_minus': 33.7*ms,
					'tau_x': 946.*ms,
					'tau_y': 27.*ms,
					'A_LTP_2': 6.1e-3,
					'A_LTD_2': 1.6e-3,
					'A_LTP_3': 6.7e-3,
					'A_LTD_3': 1.4e-3,
					'init_weight': 0.5,
					'w_max': 1.,
					'w_min': 0.,
					'ampa_max_cond': 5.e-8*siemens,
					'w_puppet': 0.125
}

# From example code/Clopath 2010:
# REMARK: The following parameters are not the ones reported in Clopath's article, rather they have been modified in order to exploits
#		  Clopath's rule with a LIF neuron (instead of using AdEX neuron as in the paper)
Clopath_LIF_parameters = {'PlasticityRule': 'Clopath',
				'v_increase_init': 2.*mV,
				'Theta_low': -70.6*mV,    			# depolarization threshold for plasticity (this is the value given in the paper, -70.*mV in example code)
				'Theta_high': -50.3*mV,			# threshold for potentiation
				'x_reset': 1.,       				# spike trace reset value'
				'A_LTD': 1.5e-4,     				# depression amplitude
				'A_LTP': 1.5e-2,    				# potentiation amplitude
				'tau_lowpass1': 40*ms,   			# timeconstant for low-pass filtered voltage
				'tau_lowpass2': 30*ms,   			# timeconstant for low-pass filtered voltage
				'tau_homeo': 1000*ms,  				# homeostatic timeconstant
				'v_target': 12*mV**2,
				'ampa_max_cond': 5.e-8*siemens,  	# Ampa maximal conductance
				'w_max': 1.,
				'w_min': 0.,
				'init_weight': 0.5,       			# initial synaptic weight
				'a_prot': 1				  			# coefficient made up since NA in the article
}

# Parameters of Sjoestroem case, as in Meissner article
Clopath_exIF_parameters = {'PlasticityRule': 'Clopath',
				'v_increase_init': 2.*mV,
				'Theta_low': -70.6*mV,      	# depolarization threshold for plasticity (this is the value given in the paper, -70.*mV in example code)
				'Theta_high': -52.*mV,			# threshold for potentiation
				'x_reset': 1,       			# spike trace reset value'
				'A_LTD': 6.e-4,     			# depression amplitude
				'A_LTP': 8.e-5,    			# potentiation amplitude
				'tau_lowpass1': 10*ms,   		# = tau_minus
				'tau_lowpass2': 7*ms,   		# = tau_plus
				'tau_homeo': 1000*ms,  			# homeostatic timeconstant
				'v_target': 12*mV**2,
				'ampa_max_cond': 5.e-8*siemens, # Ampa maximal conductance
				'w_max': 1.,
				'w_min': 0.,
				'init_weight': 0.5,             # initial synaptic weight
				'a_prot': 1				  		# coefficient made up since NA in the article (equal to Clopath with a=1)
}

# Parameters of Letzkus case, as in Meissner article
Veto_parameters = {'PlasticityRule': 'Clopath',
                               'v_increase_init': 2.*mV,
                               'Theta_low': 1.6*mV,      		# depolarization threshold for plasticity (this is the value given in the paper, -70.*mV in example code)
							   'Theta_high': 70.6*mV,			# threshold for potentiation
							   'x_reset': 1.,       			# spike trace reset value'
                               'A_LTD': 2.e-4,     				# depression amplitude
                               'A_LTP': 75e-4,    				# potentiation amplitude     # WARNING: Here the value should be much smaller for AdEx neuron (10.e-5) but since LIF is used the same value of clopath is kept
                               'tau_lowpass1': 68*ms,   		# timeconstant for low-pass filtered voltage
                               'tau_lowpass2': 210*ms,   		# timeconstant for low-pass filtered voltage
                               'tau_homeo': 1000*ms,  			# homeostatic timeconstant
                               'v_target': 12*mV**2,    		#(GD) be careful that the reference value is already squared here
                               'ampa_max_cond': 5.e-8*siemens,  # Ampa maximal conductance
                               'w_max': 1.,
							   'w_min': 0.,
                               'init_weight': 0.5,             	# initial synaptic weight
                               'v_lowpass1': -70.*mV,
                               'v_lowpass2': -70.*mV,
                               'v_homeo': 0,
                               'x_trace_LTP':0.,
                               'x_trace_LTD':0.,         		# trace for LTD, introduced with the veto rule
                               'a_prot': 0.18                 	# coefficient introduced with veto rule
}

default_consolidationParams = { 'tau_w_ampa' : 3600*second,
		'k_w' : 1.,
		'w0' : 1.,
		'c_w' : 0.5,
		'tau_z' : 36000*second,
		'k_z' : 1.,
		'z0' : 1.,
		'c_z' : 0.5
}

# parameters taken from Ziegler et al. 2015
triplet_consolidation = { 'tau_weight': 200*second,
		'tau_T': 200*second,
		'tau_z': 200*second,
		'k_w': 3,
		'sigma': 1e-4,
		'a_wT': 3.5,
		'a_Tw': 1.3,
		'a_zT': 0.95,
		'tau_gamma': 600*second,
		'Theta_gamma': 0.37,
		'A_LTP': 5*e-4,
		'A_LTD': 2*e-4,
		'tau_plus': 16.8*ms,
		'tau_minus': 33.7*ms,
		'tau_triplett': 40*ms,
		'k_up': 1*Hz,
		'k_down': 1/7200*Hz
}
