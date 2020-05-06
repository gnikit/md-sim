import unittest
import numpy as np
import os
from subprocess import call


try:
    # Change to the file directory
    os.chdir(os.path.dirname(os.path.realpath(__file__)))
except OSError:
    print('Cannot resolve directory where md-sim is installed')


def run_file(test_dir, xml_input):

    # See if the test directory exists
    try:
        os.chdir(f'./{test_dir}')

    except OSError:
        print('-' * 80)
        print(f'Test directory: {test_dir}, does not exist.')
        return False

    # See if the executable exists
    if not os.path.isfile('../../bin/md'):
        print('-' * 80)
        print('Binary md not found in bin, calling make')
        call(['make', '-j8', '-C', '../..'])

    # See if the input test file exists
    try:
        call(['make'])
        call(['../../bin/md', f'{xml_input}'])
        os.chdir('..')
        return True

    except FileNotFoundError:
        print('-'*80)
        print(f'Input file: {xml_input} could not be found.')
        os.chdir('..')
        return False

    # Catch all
    except:
        print('-'*80)
        print(f'Error running test: {test_dir}/{xml_input}')
        return False


def get_test_result(test_results):
    if False in test_results:
        print('Failure: ', test_results)
        return False
    else:
        return True


class TestPotentials(unittest.TestCase):

    def setUp(self):
        print('-'*80)

    def tearDown(self):
        call(['make', 'clean'])
        os.chdir('..')
        print('-'*80)

    def test_bounded_inverse_power(self):
        f_dir = 'bip_short'
        print(f'Entering: {f_dir}')

        setup_ok = run_file(f_dir, 'bip_potential_test.xml')
        self.assertTrue(setup_ok)

        os.chdir(f_dir)
        test_results = []

        test_string = '_step_2000_particles_512_rho_0.5000_T_1.0000' \
            + '_n_8.00_A_0.50000.log'
        sim_name = 'bip_potential_test_'

        # Test variables
        rdf_test = np.loadtxt(
            f'{sim_name}RDF{test_string}', delimiter=',', unpack=True,
            skiprows=1)
        data_test = np.loadtxt(
            f'{sim_name}Data{test_string}', delimiter=',', unpack=True,
            skiprows=1)
        pos_test = np.loadtxt(
            f'{sim_name}Positions_Velocities{test_string}', delimiter=',',
            unpack=True)

        # Reference variables
        rdf_ref = np.loadtxt('test_rdf.log', delimiter=',', unpack=True)
        data_ref = np.loadtxt('test_data.log', delimiter=',', unpack=True)
        pos_ref = np.loadtxt('test_positions.log', delimiter=',', unpack=True)

        # Test data for consistency
        result = np.allclose(rdf_ref, rdf_test)
        print(f'Testing RDF for consistency with reference: {result}')
        test_results.append(result)

        result = np.allclose(data_ref, data_test)
        print(f'Testing all Data for consistency with reference: {result}')
        test_results.append(result)

        result = np.allclose(pos_ref, pos_test)
        print(
            f'Testing Positions & Velocities for consistency with reference:',
            f'{result}')
        test_results.append(result)

        self.assertTrue(get_test_result(test_results))

        print(f'Exiting: {f_dir}')

    def test_lennard_jones(self):
        f_dir = 'lj_short'
        print(f'Entering: {f_dir}')

        setup_ok = run_file(f_dir, 'lj_potential_test.xml')
        self.assertTrue(setup_ok)

        os.chdir(f_dir)
        test_results = []

        test_string = '_step_2000_particles_512_rho_0.5000_T_1.0000.log'
        sim_name = 'lj_potential_test_'

        # Test variables
        rdf_test = np.loadtxt(
            f'{sim_name}RDF{test_string}', delimiter=',', unpack=True)
        data_test = np.loadtxt(
            f'{sim_name}Data{test_string}', delimiter=',', unpack=True)
        pos_test = np.loadtxt(
            f'{sim_name}Positions_Velocities{test_string}', delimiter=',', unpack=True)

        # Reference variables
        rdf_ref = np.loadtxt('test_rdf.log', delimiter=',', unpack=True)
        data_ref = np.loadtxt('test_data.log', delimiter=',', unpack=True)
        pos_ref = np.loadtxt('test_positions.log', delimiter=',', unpack=True)

        # Test data for consistency
        result = np.allclose(rdf_ref, rdf_test)
        print(f'Testing RDF for consistency with reference: {result}')
        test_results.append(result)

        result = np.allclose(data_ref, data_test)
        print(f'Testing all Data for consistency with reference: {result}')
        test_results.append(result)

        result = np.allclose(pos_ref, pos_test)
        print(
            f'Testing Positions & Velocities for consistency with reference:',
            f'{result}')

        test_results.append(result)

        self.assertTrue(get_test_result(test_results))

        print(f'Exiting: {f_dir}')

    def test_exponential(self):
        f_dir = 'exp_short'
        print(f'Entering: {f_dir}')

        setup_ok = run_file(f_dir, 'exp_potential_test.xml')
        self.assertTrue(setup_ok)

        os.chdir(f_dir)
        test_results = []

        test_string = '_step_2000_particles_512_rho_0.5000_T_1.0000' \
            + '_n_8.00_A_0.50000.log'
        sim_name = 'exp_potential_test_'

        # Test variables
        rdf_test = np.loadtxt(
            f'{sim_name}RDF{test_string}', delimiter=',', unpack=True)
        data_test = np.loadtxt(
            f'{sim_name}Data{test_string}', delimiter=',', unpack=True)
        pos_test = np.loadtxt(
            f'{sim_name}Positions_Velocities{test_string}', delimiter=',', unpack=True)

        # Reference variables
        rdf_ref = np.loadtxt('test_rdf.log', delimiter=',', unpack=True)
        data_ref = np.loadtxt('test_data.log', delimiter=',', unpack=True)
        pos_ref = np.loadtxt('test_positions.log', delimiter=',', unpack=True)

        # Test data for consistency
        result = np.allclose(rdf_ref, rdf_test)
        print(f'Testing RDF for consistency with reference: {result}')
        test_results.append(result)

        result = np.allclose(data_ref, data_test)
        print(f'Testing all Data for consistency with reference: {result}')
        test_results.append(result)

        result = np.allclose(pos_ref, pos_test)
        print(
            f'Testing Positions & Velocities for consistency with reference:',
            f'{result}')

        test_results.append(result)

        self.assertTrue(get_test_result(test_results))

        print(f'Exiting: {f_dir}')

    def test_gaussian_core_model(self):
        f_dir = 'gcm_short'
        print(f'Entering: {f_dir}')

        setup_ok = run_file(f_dir, 'gcm_potential_test.xml')
        self.assertTrue(setup_ok)

        os.chdir(f_dir)
        test_results = []

        test_string = '_step_2000_particles_512_rho_0.5000_T_1.0000.log'
        sim_name = 'gcm_potential_test_'

        # Test variables
        rdf_test = np.loadtxt(
            f'{sim_name}RDF{test_string}', delimiter=',', unpack=True)
        data_test = np.loadtxt(
            f'{sim_name}Data{test_string}', delimiter=',', unpack=True)
        pos_test = np.loadtxt(
            f'{sim_name}Positions_Velocities{test_string}', delimiter=',', unpack=True)

        # Reference variables
        rdf_ref = np.loadtxt('test_rdf.log', delimiter=',', unpack=True)
        data_ref = np.loadtxt('test_data.log', delimiter=',', unpack=True)
        pos_ref = np.loadtxt('test_positions.log', delimiter=',', unpack=True)

        # Test data for consistency
        result = np.allclose(rdf_ref, rdf_test)
        print(f'Testing RDF for consistency with reference: {result}')
        test_results.append(result)

        result = np.allclose(data_ref, data_test)
        print(f'Testing all Data for consistency with reference: {result}')
        test_results.append(result)

        result = np.allclose(pos_ref, pos_test)
        print(
            f'Testing Positions & Velocities for consistency with reference:',
            f'{result}')
        test_results.append(result)

        self.assertTrue(get_test_result(test_results))

        print(f'Exiting: {f_dir}')


class TestLatticeStructures(unittest.TestCase):

    def setUp(self):
        print('-'*80)

    def tearDown(self):
        call(['make', 'clean'])
        os.chdir('..')
        print('-'*80)

    def test_simple_cubic_lattice(self):
        f_dir = 'simple_cubic_lattice_short'
        print(f'Entering: {f_dir}')

        setup_ok = run_file(f_dir, 'sc_test.xml')
        self.assertTrue(setup_ok)

        os.chdir(f_dir)
        test_results = []

        test_string = '_step_1_particles_512_rho_0.5000_T_1.0000' \
            + '_n_8.00_A_0.50000.log'
        sim_name = 'sc_test_'

        # Test variables
        pos_test = np.loadtxt(
            f'{sim_name}Positions_Velocities{test_string}', usecols=(1, 2, 3),
            delimiter=',', unpack=True)

        # Reference variables
        pos_ref = np.loadtxt('test_positions.log',
                             usecols=(1, 2, 3), delimiter=',', unpack=True)

        # Test data for consistency
        result = np.allclose(pos_ref, pos_test)
        print(
            f'Testing Positions & Velocities for consistency with reference:',
            f'{result}')
        test_results.append(result)

        self.assertTrue(get_test_result(test_results))

        print(f'Exiting: {f_dir}')

    def test_fcc_lattice(self):
        f_dir = 'fcc_lattice_short'
        print(f'Entering: {f_dir}')

        setup_ok = run_file('fcc_lattice_short', 'fcc_test.xml')
        self.assertTrue(setup_ok)

        os.chdir('fcc_lattice_short')
        test_results = []

        test_string = '_step_1_particles_2048_rho_0.5000_T_1.0000' \
            + '_n_8.00_A_0.50000.log'
        sim_name = 'fcc_test_'

        # Test variables
        pos_test = np.loadtxt(
            f'{sim_name}Positions_Velocities{test_string}', usecols=(1, 2, 3),
            delimiter=',', unpack=True)

        # Reference variables
        pos_ref = np.loadtxt('test_positions.log',
                             usecols=(1, 2, 3), delimiter=',', unpack=True)

        # Test data for consistency
        result = np.allclose(pos_ref, pos_test)
        print(
            f'Testing Positions & Velocities for consistency with reference:',
            f'{result}')
        test_results.append(result)

        self.assertTrue(get_test_result(test_results))

        print(f'Exiting: {f_dir}')

    def test_bcc_lattice(self):
        f_dir = 'bcc_lattice_short'
        print(f'Entering: {f_dir}')

        setup_ok = run_file('bcc_lattice_short', 'bcc_test.xml')
        self.assertTrue(setup_ok)

        os.chdir('bcc_lattice_short')
        test_results = []

        test_string = '_step_1_particles_1024_rho_0.5000_T_1.0000' \
            + '_n_8.00_A_0.50000.log'
        sim_name = 'bcc_test_'

        # Test variables
        pos_test = np.loadtxt(
            f'{sim_name}Positions_Velocities{test_string}', usecols=(1, 2, 3),
            delimiter=',', unpack=True)

        # Reference variables
        pos_ref = np.loadtxt('test_positions.log',
                             usecols=(1, 2, 3), delimiter=',', unpack=True)

        # Test data for consistency
        result = np.allclose(pos_ref, pos_test)
        print(
            f'Testing Positions & Velocities for consistency with reference:',
            f'{result}')
        test_results.append(result)

        self.assertTrue(get_test_result(test_results))

        print(f'Exiting: {f_dir}')


class Test3DVisualisation(unittest.TestCase):

    def setUp(self):
        print('-'*80)

    def tearDown(self):
        call(['make', 'clean'])
        os.chdir('..')
        print('-'*80)

    def test_visualisation(self):
        f_dir = '3D_visualisation'
        print(f'Entering: {f_dir}')

        setup_ok = run_file(f_dir, '3D_visualisation.xml')
        self.assertTrue(setup_ok)

        os.chdir(f_dir)
        test_results = []

        steps = 10
        sim_name = 'visual_'

        # Test variables
        xyz_tracks = []
        xyz_tracks_ref = []
        for i in range(steps):

            # Load test data
            xyz_tracks.append(np.loadtxt(
                f'{sim_name}xyz_data_{i}.csv', delimiter=',', skiprows=1))

            # Reference variables
            xyz_tracks_ref.append(np.loadtxt(
                f'test_xyz_{i}.csv', delimiter=',', skiprows=1))

        # Test for consistency
        result = np.allclose(xyz_tracks, xyz_tracks_ref, rtol=1.0e-4)
        print(f'Testing x-axis tracking consistency with reference: {result}')
        test_results.append(result)

        self.assertTrue(get_test_result(test_results))

        print(f'Exiting: {f_dir}')


class TestConstructors(unittest.TestCase):

    def setUp(self):
        print('-'*80)

    def tearDown(self):
        call(['make', 'clean'])
        os.chdir('..')
        print('-'*80)

    def test_schema_missing_options(self):
        f_dir = 'schema_missing_options'
        print(f'Entering: {f_dir}')

        setup_ok = run_file(f_dir, 'missing_options.xml')
        self.assertTrue(setup_ok)

        os.chdir(f_dir)
        test_results = []

        test_string = 'Data_step_100_particles_512_rho_0.5000_T_0.5000.log'
        sim_name = 'missing_options_'

        # Test variables
        data = np.loadtxt(f'{sim_name}{test_string}',
                          delimiter=',')

        # Reference variables
        data_ref = np.loadtxt('test_data.log', delimiter=',')

        # Test for consistency
        result = np.allclose(data, data_ref)
        print(f'Testing Data file consistency with reference: {result}')

        test_results.append(result)

        self.assertTrue(get_test_result(test_results))

        print(f'Exiting: {f_dir}')


class TestIOOptions(unittest.TestCase):

    def setUp(self):
        print('-'*80)

    def tearDown(self):
        call(['make', 'clean'])
        os.chdir('..')
        print('-'*80)

    def test_io_disable(self):
        f_dir = 'io_option_flags'
        print(f'Entering: {f_dir}')

        setup_ok = run_file(f_dir, 'disable_all_io.xml')
        self.assertTrue(setup_ok)

        os.chdir(f_dir)
        test_results = []

        test_string = 'Data_step_100_particles_125_rho_1.0000_T_1.0000.log'
        sim_name = 'disable_all_io_'

        # Test variables
        data = np.loadtxt(f'{sim_name}{test_string}',
                          delimiter=',')

        # Reference variables
        data_ref = np.loadtxt('test_data.log', delimiter=',')

        # Test for consistency
        result = np.allclose(data, data_ref)
        print(f'Testing Data file consistency with reference: {result}')

        test_results.append(result)

        # Test if an RDF file has been generated in the output dir
        files = os.listdir('./')
        result = not 'disable_all_io_RDF_step_100_particles_125_rho_1.0000_T_1.0000.log' in files
        test_results.append(result)

        self.assertTrue(get_test_result(test_results))

        print(f'Exiting: {f_dir}')


class TestCompression(unittest.TestCase):

    def setUp(self):
        print('-'*80)

    def tearDown(self):
        call(['make', 'clean'])
        os.chdir('..')
        print('-'*80)

    def test_compression_summary_stats(self):
        f_dir = 'compress_summary_stats'
        print(f'Entering: {f_dir}')

        setup_ok = run_file(f_dir, 'compress_stats_test.xml')
        self.assertTrue(setup_ok)

        os.chdir(f_dir)
        test_results = []

        test_string = ('Compression_statistics_step_10_particles_64'
                       '_rho_0.5000_T_0.5000_n_0.00_A_0.00000.log')
        sim_name = 'compress_stats_test_'

        # Test variables
        data = np.loadtxt(f'{sim_name}{test_string}', delimiter=',')

        # Reference variables
        data_ref = np.loadtxt('test_data.log', delimiter=',')

        # Test for consistency
        result = np.allclose(data, data_ref)
        print(f'Testing Data file consistency with reference: {result}')

        test_results.append(result)

        self.assertTrue(get_test_result(test_results))

        print(f'Exiting: {f_dir}')


class TestIterativeAlgorithms(unittest.TestCase):

    def setUp(self):
        print('-'*80)

    def tearDown(self):
        call(['make', 'clean'])
        os.chdir('..')
        print('-'*80)

    def test_velocity_verlet(self):
        f_dir = 'velocity_verlet'
        print(f'Entering: {f_dir}')

        setup_ok = run_file(f_dir, 'verlet.xml')
        self.assertTrue(setup_ok)

        setup_ok = run_file(f_dir, 'velocity_verlet.xml')
        self.assertTrue(setup_ok)

        os.chdir(f_dir)
        test_results = []

        test_string = ('Data_step_5000_particles_512_rho_0.5000_T_1.0000.log')
        sim_name = 'velocity_verlet_'

        # Test variables
        data = np.loadtxt(f'{sim_name}{test_string}', delimiter=',')

        # Reference variables
        data_ref = np.loadtxt(f'verlet_{test_string}', delimiter=',')

        # BUG: potential bug I think these should also agree
        skip_columns = [8, 9, 10, 11]  # skip VAF and SF
        headers = ['step', 'rho', 'T', 'U', 'K', 'Pc',
                   'Pk', 'MSD', 'VAF', 'SFx', 'SFy', 'SFz']
        # Test for consistency
        for j in range(len(data[0])):

            if j in skip_columns:
                print(f'  Skipping column {j}: {headers[j]}')
                pass
            else:
                result = np.allclose(data[-1][j], data_ref[-1][j], rtol=1e-2)
                print(f'  Testing velocity verlet data file column {j}: {headers[j]}'
                      f' against results with explicit Verlet: {result}')
                test_results.append(result)

        self.assertTrue(get_test_result(test_results))

        print(f'Exiting: {f_dir}')


# class TestBoundaryConditions(unittest.TestCase):

#     def setUp(self):
#         print('-'*80)

#     def tearDown(self):
#         call(['make', 'clean'])
#         os.chdir('..')
#         print('-'*80)

#     def test_reflective_boundary(self):
#         f_dir = 'reflective_boundary'
#         print(f'Entering: {f_dir}')

#         setup_ok = run_file(f_dir, 'reflective_bcs_test.xml')
#         self.assertTrue(setup_ok)

#         os.chdir(f_dir)
#         test_results = []

#         steps = 100
#         sim_name = 'reflective_bcs_'

#         # Test variables
#         xyz_tracks = []
#         xyz_tracks_ref = []
#         for i in range(steps):

#             # Load test data
#             xyz_tracks.append(np.loadtxt(
#                 f'{sim_name}xyz_data_{i}.csv', delimiter=',', skiprows=1))

#             # Reference variables
#             xyz_tracks_ref.append(np.loadtxt(
#                 f'test_xyz_{i}.csv', delimiter=',', skiprows=1))

#         # Test for consistency
#         result = np.allclose(xyz_tracks, xyz_tracks_ref, rtol=1.0e-4)
#         print(f'Testing x-axis tracking consistency with reference: {result}')
#         test_results.append(result)

#         self.assertTrue(get_test_result(test_results))

#         print(f'Exiting: {f_dir}')


if __name__ == '__main__':
    unittest.main()
