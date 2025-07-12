import numpy as np
import matplotlib.pyplot as plt
from particles.photon import photon
from detectors.naitl import NaITl
from detectors.si import Si
from multiprocessing import Pool, Manager
from functools import partial

def get_detectors():
    return [
        NaITl(100, 100, 0, 60, 60, 60),  # 6x6x6cm NaITl
        Si(100, 0, 100, 60, 60, 2)       # 6x6x0.2cm Si
    ]

def simulate_batch(batch_size: int, E: float, detectors):
    theta = np.random.uniform(0, 2 * np.pi, batch_size)
    phi = np.random.uniform(-np.pi/2, np.pi/2, batch_size)
    directions = photon.batch_calculate_direction(theta, phi)
    origins = np.zeros((batch_size, 3))
    total = [0] * len(detectors)
    detected_counts = [0] * len(detectors)
    energies = [[] for _ in detectors]

    for idx, det in enumerate(detectors):
        mask, batch_energies, total_dtd = det.detects_batch(origins, directions, theta, phi, E, batch_size)
        total[idx] += np.sum(total_dtd)
        detected_counts[idx] += np.sum(mask)
        energies[idx].extend(batch_energies)

    return detected_counts, energies, total

def combine_results(results):
    n_detectors = len(results[0][0])
    total_detected = [0]*n_detectors
    total_energies = [[] for _ in range(n_detectors)]
    total_totals = [0]*n_detectors

    for detected_counts, energies, total in results:
        for i in range(n_detectors):
            total_detected[i] += detected_counts[i]
            total_totals[i] += total[i]
            total_energies[i].extend(energies[i])

    return total_detected, total_energies, total_totals

def simulate(n_photons: int, batch_size: int, E: float, detectors: list, progress_callback=None):
    if E == 0:
        print("WARNING. 0 PASSED INTO E.")

    n_batches = n_photons // batch_size
    args = [(batch_size, E, detectors) for _ in range(n_batches)]

    with Pool() as pool:
        results = []
        for i, result in enumerate(pool.starmap(simulate_batch, args)):
            results.append(result)
            if progress_callback is not None:
                progress = int((i + 1) / n_batches * 100)
                progress_callback(progress)

    return combine_results(results)

def plot_spectra(energies: list, bins: int = 1024, energy_range=(0, 1000)):
    num_det = len(energies)
    fig, axs = plt.subplots(num_det, 1, figsize=(10, 4 * num_det), squeeze=False)
    for i, (ax, detector_energies) in enumerate(zip(axs.flat, energies)):
        ax.hist(detector_energies, bins=bins, range=energy_range, histtype='step')
    plt.tight_layout()
    return fig
