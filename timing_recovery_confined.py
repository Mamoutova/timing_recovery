import numpy as np
import matplotlib.pyplot as plot

from math import modf

RETARD = "retard"
ADVANCE = "advance"
HOLD = "hold"

class FIR_filter(object):
    """
    Basic FIR filter
    """
    def __init__(self, weights):
        self.weights = np.copy(weights)
        self.Areg = np.zeros(len(weights), dtype=np.float32)

    def update(self, x_k):
        self.Areg = np.roll(self.Areg, 1)
        self.Areg[0] = x_k
        return np.inner(self.Areg, self.weights)

class TED(object):
    def __init__(self, loop_filter_length, threshold):
        self.x_k_prev = 0
        self.a_k_prev = 0
        self.loop_filter_length = loop_filter_length
        self.loop_filter = np.zeros(self.loop_filter_length, dtype=np.float32)
        self.error = 0
        self.threshold = threshold

    def calc(self, x_k, a_k):
        # approximation of timing function        
        loop_filter_in = ((a_k**2.0 - 5.0)*(a_k*self.x_k_prev - self.a_k_prev*x_k)/16.0 )
        # Loop filter 
        self.loop_filter = np.roll(self.loop_filter, 1)
        self.loop_filter[0] = loop_filter_in        
        self.error =  sum(self.loop_filter) / self.loop_filter_length
        # Save values for next calc
        self.x_k_prev = x_k
        self.a_k_prev = a_k
        # Make a decision of phase shift
        if self.error > self.threshold: return RETARD
        elif self.error < -self.threshold: return ADVANCE
        else: return HOLD

class counter_divider(object):
    def __init__(self, data_f, symbol_f, phase_step, loop_filter_length):
        self.loop_filter_counter = loop_filter_length
        self.correction_counter = 0
        self.L = loop_filter_length
        self.phase_step = phase_step
        # Calculate decimation factor
        self.fract, self.decimation_factor = modf(data_f/symbol_f)
        if self.fract!=0:  self.decimation_factor_correction = int(1/self.fract)
        else: self.decimation_factor_correction = 0        


    def calc_next_sample(self, phase_correction):
        # calculate next sampling point for a decimation factor
        if (self.correction_counter == self.decimation_factor_correction-1) and (self.decimation_factor_correction != 0):
            self.correction_counter = 0
            next_s_point = int (self.decimation_factor + 1)
        else:
            self.correction_counter += 1
            next_s_point = int(self.decimation_factor)
        # shift next sampling point based on TED decision once every L samples
        if self.loop_filter_counter == 0:
            self.loop_filter_counter = self.L-1
            if phase_correction==ADVANCE:   next_s_point += self.phase_step
            elif phase_correction==RETARD:  next_s_point -= self.phase_step
            elif phase_correction==HOLD:    next_s_point = next_s_point                
        else:
            self.loop_filter_counter -= 1
        return next_s_point

if __name__ == "__main__":

    # Uncomment to select input sample signal
    filename = "data_sample_short_line.npy"
    # filename = "data_sample_long_line.npy"
    
    # Select digital filter parameters, uncomment to select weights
    ffe_ena = 0                     # filter on/off
    ffe_weights = [1.0, -0.75]
    # ffe_weights = [1.0, -1.325]     # filter weights

    # Signal parameters
    data_f = 12.5 * 10**6           # sampling frequency
    symbol_f = 320 * 10**3          # baud frequency

    # Initial phase position, in samples
    initial_phase_offset = 0

    # Load data
    data = np.load(filename)
    graph_length = data.shape[0]

    # Create dsp instances: digital filter, TED calculator and counter-divider
    if ffe_ena: ffe = FIR_filter(weights=ffe_weights)
    ted = TED(loop_filter_length=120, threshold=0)
    counter_divider = counter_divider(data_f=data_f, symbol_f=symbol_f, phase_step=1, loop_filter_length=120)    

    # Process all samples
    i = initial_phase_offset
    sampling_points = []
    array_t = []
    while True:
        if i >= graph_length: break
        # mark sampling point
        sampling_points.append(data[i])
        # read sample
        x_k = data[i]
        # preequalizer
        if ffe_ena:
            x_k = ffe.update(x_k)
            x_k = int(x_k)
        # slicer 
        if x_k >= 0: a_k = 1
        else: a_k = -1
        # TED
        phase_correction = ted.calc(x_k, a_k)
        # phase correction
        next_s_point = counter_divider.calc_next_sample(phase_correction)
        # set next sampling point
        array_t.append(i) 
        i += next_s_point     

    fig, axs = plot.subplots(1)
    axs.grid(b=True, which='major', axis='y')
    axs.plot(data)
    axs.stem(array_t, sampling_points, linefmt='r-', markerfmt='ro')
    plot.show()
