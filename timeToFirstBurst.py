import h5py
import numpy as np
import matplotlib.pyplot as plt
import bluepyopt as bpop
import bluepyopt.ephys as ephys
import efel

def computeTimeToFirstBurst(V, t):
    V = np.array(V)
    t = np.array(t)

    trace = {
        'T': t,
        'V': V,
        'stim_start': [t[0]],
        'stim_end': [t[-1]]
    }

    traces = [trace]


    # Define feature list
    features = ['peak_time']
    #efel.set_setting('peak_threshold', 30.0)
    # Extract features
    results = efel.get_feature_values(traces, features)

    #############
    peakTimes = np.array(results[0]['peak_time'])
    ISI = np.diff(peakTimes)

    # Find last ISI that is too long
    long_isi_indices = np.where(ISI >= 3000)[0]

    if len(long_isi_indices) == 0:
        # all ISIs are < 1000
        regular_burst_start_time = peakTimes[0]
    else:
        last_long_isi_index = long_isi_indices[-1]
        regular_burst_start_time = peakTimes[last_long_isi_index + 1]

    #print(f"Regular bursting begins at: {regular_burst_start_time:.2f} ms")

    return regular_burst_start_time

print("timeToFirstBurst.py loaded")
