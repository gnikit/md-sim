import os
import sys
import numpy as np
from scipy.interpolate import interp1d
from scipy.special import erf
from scipy import stats


def mb_speed_distribution():
    """
    Maxwell-Boltzmann's speed distribution
    Takes in speeds range, temperature and mass.
    :return: numpy array of M-B speed distribution
    """
    term1 = (MASS / (2 * np.pi * K_BOLTZMANN * TEMPERATURE)) ** 1.5
    term2 = 4 * np.pi * RANGE_OF_SPEEDS_MONTE_CARLO ** 2
    term3 = np.exp(-MASS * RANGE_OF_SPEEDS_MONTE_CARLO ** 2 / (2 * K_BOLTZMANN * TEMPERATURE))
    return term1 * term2 * term3


def mb_cdf():
    """
    Maxwell-Boltzmann's speed distribution Cumulative Distribution Function (CDF)
    Takes in speeds range, temperature and mass.
    :return: numpy array of CDF
    """
    term0 = np.sqrt(K_BOLTZMANN * TEMPERATURE / MASS)
    term1 = erf(RANGE_OF_SPEEDS_MONTE_CARLO / (np.sqrt(2) * term0))
    term2 = np.sqrt(2 / np.pi) * RANGE_OF_SPEEDS_MONTE_CARLO
    term3 = np.exp(-RANGE_OF_SPEEDS_MONTE_CARLO ** 2 / (2 * term0 ** 2)) / term0
    return term1 - term2 * term3


class VelGen:
    def __init__(self):
        # Generates random speeds according to
        # Maxwell-Boltzmann's speed distribution
        if TEMPERATURE == 0:
            self.speeds = np.zeros(NUMBER_OF_PARTICLES)
        else:
            # Create interpolation function to CDF
            inv_cdf = interp1d(mb_cdf(), RANGE_OF_SPEEDS_MONTE_CARLO)
            # Generate random speeds
            self.speeds = inv_cdf(np.random.random(NUMBER_OF_PARTICLES))

    def get_velocities(self):
        """
        Converts random speeds into their x, y, z components.
        Takes in random speeds
        :return: numpy arrays of speeds and their x,y,z components
        """
        theta = np.arccos(np.random.uniform(-1, 1, NUMBER_OF_PARTICLES))
        phi = np.random.uniform(0, 2 * np.pi, NUMBER_OF_PARTICLES)
        _vx = self.speeds * np.cos(phi)
        _vy = self.speeds * np.sin(phi)
        _vz = self.speeds * np.cos(theta)
        if NUMBER_OF_DIMENSIONS == 3:
            _vx *= np.sin(theta)
            _vy *= np.sin(theta)
            return np.array([_vx, _vy, _vz], dtype=np.double)
        return np.array([_vx, _vy], dtype=np.double)


if __name__ == '__main__':
    RANGE_OF_SPEEDS_MONTE_CARLO = np.arange(0, 1000, 0.001)
    # high enough upper range for speeds to make interpolation of high
    # temperatures possible, otherwise add ", fill_value="extrapolate"" to interp1d
    # at the cost of potential accuracy

    K_BOLTZMANN = 1.
    MASS = 1.
    NUMBER_OF_DIMENSIONS = 3
    print(str(sys.argv))
    try:
        if len(sys.argv) != 3:
            NUMBER_OF_PARTICLES = 10 ** 3
            TEMPERATURE = 1.0
            raise ValueError('Invalid number or type of arguments\nArguments initialised with default values')
        NUMBER_OF_PARTICLES = int(sys.argv[1])  # 1st argument
        TEMPERATURE = float(sys.argv[2])        # 2nd argument
    except ValueError:
        NUMBER_OF_PARTICLES = 10 ** 3
        TEMPERATURE = 1
    file_id = '_particles_' + str(NUMBER_OF_PARTICLES) + '_T_' + "{:.4f}".format(TEMPERATURE) + '.txt'
    # TODO: add os handler
    print(os.getcwd())
    os.chdir('../data')
    a = VelGen()
    vx, vy, vz = a.get_velocities()
    np.savetxt('vx' + file_id, vx, delimiter='\n')
    np.savetxt('vy' + file_id, vy, delimiter='\n')
    np.savetxt('vz' + file_id, vz, delimiter='\n')
else:
    print("MD module loaded")
