#include "MD.h"



MD::MD()
{
	power = 6;
	A = 1;

	Nx = Ny = Nz = 8;
	N = Nx*Ny*Nz;
	cut_off = L / 2;

	rg = cut_off;
	Nhist = 100;							         // might change # of part in bin
	std::vector <double> gr(Nhist + 1);
	std::fill(gr.begin(), gr.end(), 0);              // sets all gr elements to 0  // IMPORTANT GOES OUT OF SCOPE SOMETIMES
	fx.resize(N, 0);
	fy.resize(N, 0);
	fz.resize(N, 0);
	
}
MD::MD(int POWER, double A_cst, size_t run_number)
{
	power = POWER;
	A = A_cst;
	N_max = run_number; 

	Nx = Ny = Nz = 8;
	N = Nx*Ny*Nz;
	scale = pow((N / rho), (1.0 / 3.0)) / Nx; // scalling factor for length of box
	L = pow((N / rho), 1.0 / 3.0);			  // L depends on rho
	Vol = N / rho;

	cut_off = L / 2;

	rg = cut_off;
	Nhist = 100;		    // might change # of part in bin
	dr = rg / Nhist;

	gr.resize(Nhist + 1, 0);             // sets all gr elements to 0  // IMPORTANT GOES OUT OF SCOPE SOMETIMES
	
	fx.resize(N, 0);
	fy.resize(N, 0);
	fz.resize(N, 0);

	/*  NAMES USED FOR PYTHON PROCESSING	        	*/
	run = std::to_string(power);	
	separator = "~";
	std::stringstream stream;
	stream << std::fixed << std::setprecision(2) << A;
	A_par = stream.str();

	path = "../../Archives of Data/Density 0.5/Isothermal~step 5000/";
	file_type = ".txt";
	KIN = "KinEn"; POT = "PotEn"; TOT = "TotEn";
	PRESSUREC = "PressureC"; PRESSUREK = "PressureK"; PCKTOT = "PCK"; TEMPERATURE = "Temperature";
	HIST = "Hist"; _VAF = "VAF"; _MSD = "MSD"; data = "data";
	Ldrx = "Loadrx"; Ldry = "Loadry"; Ldrz = "Loadrz"; Ldvx = "Loadvx"; Ldvy = "Loadvy"; Ldvz = "Loadvz";

	// Path addition
	KIN = path + KIN; POT = path + POT; TOT = path + TOT; PRESSUREC = path + PRESSUREC; PRESSUREK = path + PRESSUREK; PCKTOT = path + PCKTOT;
	TEMPERATURE = path + TEMPERATURE; HIST = path + HIST; _VAF = path + _VAF; _MSD = path + _MSD; data = path + data;
	Ldrx = path + Ldrx; Ldry = path + Ldry; Ldrz = path + Ldrz; Ldvx = path + Ldvx; Ldvy = path + Ldvy; Ldvz = path + Ldvz;

	KIN = KIN + run + separator + A_par + file_type;
	POT += run +separator + A_par + file_type;
	TOT += run + separator + A_par + file_type;
	PRESSUREC += run + separator + A_par + file_type;
	PRESSUREK += run + separator + A_par + file_type;
	PCKTOT += run + separator + A_par + file_type;
	TEMPERATURE += run + separator + A_par + file_type;
	HIST += run + separator + A_par + file_type;
	_VAF += run + separator + A_par + file_type;
	_MSD += run + separator + A_par + file_type;
	//data += run + separator + A_par + file_type;		// for future use
	Ldrx += run + separator + A_par + file_type;
	Ldry += run + separator + A_par + file_type;
	Ldrz += run + separator + A_par + file_type;
	Ldvx += run + separator + A_par + file_type;
	Ldvy += run + separator + A_par + file_type;
	Ldvz += run + separator + A_par + file_type;


}
MD::~MD()
{
}

// Methods for MD Analysis
void MD::Initialise(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, 
					std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz)
/*
	Initialises the:
	+ Position Arrays
	+ Velocity Arrays (assign random velocities)
	+ Conserves/ Scales momentum == 0
	+ Temperature
	+ Velocity Autocorrelaion Function

*/
{
	// Initialise position matrix and velocity matrix
	size_t n = 0;
	for (size_t i = 0; i < Nx; i++)
	{
		for (size_t j = 0; j < Ny; j++)
		{
			for (size_t k = 0; k < Nz; k++)
			{
				x.push_back((i + 0.5)*scale);
				y.push_back((j + 0.5)*scale);
				z.push_back((k + 0.5)*scale);

				rrx.push_back((i + 0.5)*scale);
				rry.push_back((j + 0.5)*scale);
				rrz.push_back((k + 0.5)*scale);

				// Add random initial velocities 
				vx.push_back(((double)rand() / (RAND_MAX)) + 1);
				vy.push_back(((double)rand() / (RAND_MAX)) + 1);
				vz.push_back(((double)rand() / (RAND_MAX)) + 1);

				n++;
			}
		}
	}
	// scale of x, y, z
	double mean_vx = 0;	double mean_vy = 0;	double mean_vz = 0;

	// Momentum conservation array
	for (size_t i = 0; i < N; i++)
	{
		mean_vx += vx[i] / N; // Calculating Average velocity for each dimension
		mean_vy += vy[i] / N;
		mean_vz += vz[i] / N;
	}
	for (size_t i = 0; i < N; i++)
	{
		vx[i] = vx[i] - mean_vx; // Subtracting Av. velocities from each particle
		vy[i] = vy[i] - mean_vy;
		vz[i] = vz[i] - mean_vz;
	}
	// T Calc
	KE = 0;
	for (size_t i = 0; i < N; i++)
	{
		KE += 0.5*(vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
	}
	T = KE / (1.5*N);
	scale_v = sqrt(T0 / T); // scalling factor

	// Velocity scaling
	for (size_t i = 0; i < N; i++)
	{
		vx[i] *= scale_v;
		vy[i] *= scale_v;
		vz[i] *= scale_v;
	}
	// MSD initialasation
	MSDx = x; MSDy = y; MSDz = z;

	// VAF initialasation
	Cvx = vx; Cvy = vy; Cvz = vz;
	double first_val = 0;
	for (int i = 0; i < N; i++)
	{
		first_val += (Cvx[i] * Cvx[i] + Cvy[i] * Cvy[i] + Cvz[i] * Cvz[i]) / N;
	}
	first_val /= N;
	Cr.push_back(first_val); // entry for first value
	VAF << first_val << std::endl;
}

void MD::VerletAlgorithm(std::vector<double> &rx, std::vector<double> &ry, std::vector<double> &rz,
						 std::vector<double> &vx, std::vector<double> &vy, std::vector<double> &vz,
						 std::vector<double> &rrx, std::vector<double> &rry, std::vector<double> &rrz)
{
	using namespace std;

	for (size_t i = 0; i < N; i++)
	{
		vx[i] = vx[i]* scale_v + fx[i] * dt;
		vy[i] = vy[i]* scale_v + fy[i] * dt;
		vz[i] = vz[i]* scale_v + fz[i] * dt;
		rx[i] = rx[i] + vx[i] * dt;
		ry[i] = ry[i] + vy[i] * dt;
		rz[i] = rz[i] + vz[i] * dt;

		rrx[i] = rrx[i] + vx[i] * dt;
		rry[i] = rry[i] + vy[i] * dt;
		rrz[i] = rrz[i] + vz[i] * dt;

		// Kinetic Energy Calculation
		KE += 0.5*(vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);

		// Boundary conditions Updated
		if (rx[i] > L) { rx[i] = rx[i] - L; }
		if (ry[i] > L) { ry[i] = ry[i] - L; }
		if (rz[i] > L) { rz[i] = rz[i] - L; }
		if (rx[i] < 0.0) { rx[i] = rx[i] + L; }
		if (ry[i] < 0.0) { ry[i] = ry[i] + L; }
		if (rz[i] < 0.0) { rz[i] = rz[i] + L; }
	}
}

void MD::VelocityAutocorrelationFunction(std::vector<double>& Cvx, std::vector<double>& Cvy, std::vector<double>& Cvz)
{
	double temp_=0; // resets every time step
	for (int i = 0; i < N; i++)
	{
		temp_ += (Cvx[i] * vx[i] + Cvy[i] * vy[i] + Cvz[i] * vz[i]);
	}
	temp_ /= N;
	Cr.push_back(temp_);
	VAF << temp_ << std::endl; // writes to file
}

void MD::RadialDistributionFunction()
{
	double R = 0;
	double dem;
	for (int i = 0; i < Nhist + 1; i++)
	{
		R = rg * i / Nhist; // dr = rg*Nhist
		dem = (rho * 2 * pi*R*R*N*N_step*dr);
		gr[i] = gr[i] / dem;
		Hist << gr[i] << std::endl;
	}
	//std::cout << "RDF saved" << std::endl;
}
void MD::MeanSquareDisplacement(std::vector<double>& MSDx, std::vector<double>& MSDy, std::vector<double>& MSDz)
{
	double msd_temp = 0;
	for (size_t i = 0; i < N; ++i)
	{
		msd_temp += (pow((rrx[i] - MSDx[i]), 2) + pow((rry[i] - MSDy[i]), 2) + pow((rrz[i] - MSDz[i]), 2));
	}
	msd_temp /= N;
	msd.push_back(msd_temp);
	MSD << msd_temp << std::endl;
}

// MD Simulation
void MD::Simulation()
{
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	OpenFiles();
	Initialise(rx, ry, rz, vx, vy, vz);


	std::cout << "MD Simulation running..." << std::endl;
	//std::cout << "step:\tT:\tKE:\tU:\tU+K:\tPC:\tPK:\t(PK+PC):" << std::endl;

	double xx, yy, zz, rr;
	for (N_step = 1; N_step < N_max; N_step++)
	{
		// Forces loop

		// Re-setting forces
		std::fill(fx.begin(), fx.end(), 0);
		std::fill(fy.begin(), fy.end(), 0);
		std::fill(fz.begin(), fz.end(), 0);


		U = 0; // seting Potential U to 0
		PC = 0;

		for (size_t i = 0; i < N - 2; i++)
		{
			for (size_t j = i + 1; j < N - 1; j++)
			{
				x = rx[i] - rx[j]; // Separation distance
				y = ry[i] - ry[j]; // between particles i and j
				z = rz[i] - rz[j]; // in Cartesian

				xx = rrx[i] - rrx[j];
				yy = rry[i] - rry[j];
				zz = rrz[i] - rrz[j];

				// Transposing elements with Periodic BC
				if (x > (0.5*L)) { x = x - L; }
				if (y > (0.5*L)) { y = y - L; }
				if (z > (0.5*L)) { z = z - L; }
				if (x < (-0.5*L)) { x = x + L; }
				if (y < (-0.5*L)) { y = y + L; }
				if (z < (-0.5*L)) { z = z + L; }

				/////////// for MSD  /////////////////
				if (xx >(0.5*L)) { xx = xx - L; }
				if (yy >(0.5*L)) { yy = yy - L; }
				if (zz >(0.5*L)) { zz = zz - L; }
				if (xx < (-0.5*L)) { xx = xx + L; }
				if (yy < (-0.5*L)) { yy = yy + L; }
				if (zz < (-0.5*L)) { zz = zz + L; }

				r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
				double q = sqrt(r*r + A);

				// Force loop
				if (r <= cut_off)													 // for particles within the cut off range
				{
					double ff = (power)* r * pow(q, (-power - 2));					 // Force for particles // change to dv

					fx[i] += x*ff / r;
					fx[j] -= x*ff / r;												// Canceling the ij and ji pairs
					fy[i] += y*ff / r;												// Taking the lower triangular matrix
					fy[j] -= y*ff / r;
					fz[i] += z*ff / r;
					fz[j] -= z*ff / r;

					PC += r*ff;
					U += pow(q, -power);											// Potential Calculation

																					// Radial Distribution 
					igr = Nhist* r / rg;
					gr[igr] = gr[igr] + 1;
					//rn = (igr - 0.5)*dr;
				}
			}
		}

		U /= N;																		 // Average Potential Energy per particle
		PC = PC / (3 * Vol);

		// Isothermal Calibration
		scale_v = sqrt(T0 / T);														 // using T & KE from prev timestep
		//scale_v = 1;
		KE = 0;																		 // set 0 for each step

		// Verlet Algorithm
		VerletAlgorithm(rx, ry, rz, vx, vy, vz, rrx, rry, rrz);
		///////////////////// (M)EAN (S)QUARE (D)ISPLACEMENT ///////////////////////////
		MeanSquareDisplacement(MSDx, MSDy, MSDz);

		/////////////////// (V)ELOCITY (A)UTOCORRELATION (F)UNCTION /////////////////////
		VelocityAutocorrelationFunction(Cvx, Cvy, Cvz);

		T = KE / (1.5*N);																// Average T
		PK = rho*T;																        // Kinetic part of pressure
		KE /= N;

		WriteToFiles();  // writes data to open streams
		//ShowRun(500);  // shows every 500 steps 

	}
	// Saving Last Position
	for (int el = 0; el < rx.size(); el++)
	{
		Loadrx << rx[el] << std::endl;
		Loadry << ry[el] << std::endl;
		Loadrz << rz[el] << std::endl;
		Loadvx << vx[el] << std::endl;
		Loadvy << vy[el] << std::endl;
		Loadvz << vz[el] << std::endl;
	}
	//std::cout << "Last Positions saved" << std::endl;

	RadialDistributionFunction();

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::cout << "CPU run time = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() / 60 << " min " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() % 60 << "s" << std::endl;
}

// File Handling
void MD::OpenFiles()
{
	Temperature.open(TEMPERATURE, std::ios::out | std::ios::trunc); // opens for output and deletes prev content
	KinEn.open(KIN, std::ios::out | std::ios::trunc);
	PotEn.open(POT, std::ios::out | std::ios::trunc);
	TotEn.open(TOT, std::ios::out | std::ios::trunc);
	PressureC.open(PRESSUREC, std::ios::out | std::ios::trunc);
	PressureK.open(PRESSUREK, std::ios::out | std::ios::trunc);
	PCK.open(PCKTOT, std::ios::out | std::ios::trunc);
	Hist.open(HIST, std::ios::out | std::ios::trunc);
	VAF.open(_VAF, std::ios::out | std::ios::trunc);
	MSD.open(_MSD, std::ios::out | std::ios::trunc);
	//DATA.open(data, std::ios::out | std::ios::trunc);

	Loadrx.open(Ldrx, std::ios::out | std::ios::trunc);
	Loadry.open(Ldry, std::ios::out | std::ios::trunc);
	Loadrz.open(Ldrz, std::ios::out | std::ios::trunc);
	Loadvx.open(Ldvx, std::ios::out | std::ios::trunc);
	Loadvy.open(Ldvy, std::ios::out | std::ios::trunc);
	Loadvz.open(Ldvz, std::ios::out | std::ios::trunc);

}
void MD::WriteToFiles()
{
	Temperature << T << std::endl;
	KinEn << KE << std::endl;
	PotEn << U << std::endl;
	TotEn << (U + KE) << std::endl;
	PressureC << PC << std::endl;
	PressureK << PK << std::endl;
	PCK << (PC + PK) << std::endl;
}
void MD::ShowRun(size_t step_size_show)
{
	if (N_step % step_size_show == 0 || N_step == 1) // print every 500 results
	{
		std::cout.precision(5);
		std::cout << N_step << "\t" << T << "\t" << KE << "\t" << U << "\t" << (U + KE) << "\t" << PC << "\t" << PK << "\t" << (PK + PC) << std::endl;
	}
}
