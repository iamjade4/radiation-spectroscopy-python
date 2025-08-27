import numpy as np
import matplotlib.pyplot as plt
from particles.photon import photon
from detectors.scintillators.scintillator import Scintillator
from detectors.scintillators.scintillators import *
from detectors.solidstate.si import Si
from multiprocessing import Pool

#def get_detectors():
    # return [
    #     NaITl([600, 100, 0], [6, 1, 0]/(np.sqrt(37)), 60, 120),  #cylinder with a radius of 6cm, height of 12cm, axis facing the origin
    #     Si(100, 0, 100, 60, 60, 2)       # 6x6x0.2cm Si
    # ]
    


def simulate_batch(batch_size: int, E: float, detectors, gain):
    theta = np.random.uniform(0, 2 * np.pi, batch_size)
    phi = np.random.uniform(-np.pi/2, np.pi/2, batch_size)
    directions = photon.batch_calculate_direction(theta, phi)
    origins = np.zeros((batch_size, 3))
    total = [0] * len(detectors)
    detected_counts = [0] * len(detectors)
    energies = [[] for _ in detectors]
    out_energies = [[] for _ in detectors]
    positions = [[] for _ in detectors]
    for idx, det in enumerate(detectors):
        mask, batch_energies, total_dtd, out_energies_s, positions_s = det.detects_batch(
            origins, directions, theta, phi, E, batch_size, gain
        )
        total[idx] += np.sum(total_dtd)
        detected_counts[idx] += mask #mask is no longer a list of True/False
        energies[idx].extend(batch_energies)# this adds the cumulative data, previously you were only plotting the last batch (unless thats what you wanted to do in which case, sorry)
        out_energies[idx].extend(out_energies_s)
        positions[idx].extend(positions_s)

    return detected_counts, energies, total, positions, out_energies #energies are no longer actual energies, but are the number of scintillation photons produced/100 *gain

def combine_results(results):
    n_detectors = len(results[0][0])
    total_detected = [0] * n_detectors
    total_energies = [[] for _ in range(n_detectors)]
    total_positions = [[] for _ in range(n_detectors)]
    total_out_energies = [[] for _ in range(n_detectors)]
    total_totals = [0] * n_detectors

    for detected_counts, energies, total, positions, out_energies in results:
        for i in range(n_detectors):
            total_detected[i] += detected_counts[i]
            total_totals[i] += total[i]
            total_energies[i].extend(energies[i])
            total_positions[i].extend(positions[i])
            total_out_energies[i].extend(out_energies[i])

    return total_detected, total_energies, total_totals, total_positions, total_out_energies

def simulate(n_photons: int, batch_size: int, E: float, detectors: list, gain, progress_value=None):
    if E == 0:
        print("WARNING. 0 PASSED INTO E.")

    n_batches = n_photons // batch_size
    args = [(batch_size, E, detectors, gain) for _ in range(n_batches)]
    results = []

    with Pool() as pool:
        for i, result in enumerate(pool.imap(simulate_batch_wrapper, args)):
            results.append(result)
            if progress_value is not None:
                progress_value.value = int((i + 1) / n_batches * 100)

    return combine_results(results)

def simulate_batch_wrapper(args):
    return simulate_batch(*args)

def plot_spectra(energies: list, bins: int = 1024, energy_range=(0, 1000)):
    num_det = len(energies)
    fig, axs = plt.subplots(num_det, 1, figsize=(10, 4 * num_det), squeeze=False)
    for i, (ax, detector_energies) in enumerate(zip(axs.flat, energies)):
        ax.hist(detector_energies, bins=bins, range=energy_range, histtype='step')
    plt.tight_layout()
    return fig
