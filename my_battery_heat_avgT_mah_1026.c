/*my_battery_heat_source_trial.c*/
#include "udf.h"
#include "stdio.h"
#include "math.h"
#include "mem.h"
#include "time.h"
#include "sg_pollut.h"

// Inputs for First Submodel
#define V_bat 0.0000165  // Battery volume (m^3)
#define H_chem 11740     // Total heat release of chemical reaction (J) (?)
#define H_sc 17560       // Total heat release of short circuit (J) (?)
#define rho 2600         // Density of battery (kg/m^3)
#define Cp 1720          // Specific heat capacity (J/kg.K)
#define A 0.92           // Pre-exponential factor (unitless)
#define b 28.5           // Reaction rate constant (unitless)
#define T_ref 533.15     // Normalized temperature coefficient (K) (?)
#define T_chem 413.15    // Start temperature of chemical reactions (K)
#define T_tr 533.15      // Triggering temperature of thermal runaway (K)
#define C_chem_max 12    
#define C_elec_max 12    
#define NUM_UDM 25 // (?)

// Inputs for third submodel

#define MW_H2 2.016       // Molar mass of H2 (g/mol)
#define MW_C2H4 28.05     // Molar mass of C2H4(g/mol)
#define MW_CH4 16.04      // Molar mass of CH4(g/mol)
#define MW_CO 28.01      // Molar mass of CO(g/mol)
#define MW_CO2 44.01      // Molar mass of CO2(g/mol)
#define MW_O2 32      // Molar mass of O2(g/mol)
#define P_atm 101325      // Aymospheric pressure (Pa)
#define V_void_init 0.00000152      // Initial void volume (m3)
#define d_void 0.000007594      // volume increase (m3)
#define Area_valve 0.00003848 // Safety valve area (dia 7mm) (m^2)

// Global variable of generation rates and temperature
real T_avg;  // Initial assumption (K)
real G_x_1; // H2 Gas Generation Rate (mol/s m3)
real G_x_2; // C2H4 Gas Generation Rate
real G_x_3; // CH4 Gas Generation Rate
real G_x_4; // CO Gas Generation Rate
real G_x_5; // CO2 Gas Generation Rate
real G_x_6; // O2 Gas Generation Rate
real G_total; // Total Gas Generation Rate (mol/s m3)

// number of generated moles (instantaneous)
real n_x_1_gen; // H2 number of moles
real n_x_2_gen; // C2H4 number of moles
real n_x_3_gen; // CH4 number of moles
real n_x_4_gen; // CO number of moles
real n_x_5_gen; // CO2 number of moles
real n_x_6_gen; // O2 number of moles
real n_total_gen; // Total number of generated moles

// number of generated moles (Accumulated)
real n_x_1_gen_accu; // H2 number of moles
real n_x_2_gen_accu; // C2H4 number of moles
real n_x_3_gen_accu; // CH4 number of moles
real n_x_4_gen_accu; // CO number of moles
real n_x_5_gen_accu; // CO2 number of moles
real n_x_6_gen_accu; // O2 number of moles
real n_total_gen_accu; // Total number of generated moles

// number of ejected moles (accumulated)
real n_x_1_exit_accu; // H2 number of moles
real n_x_2_exit_accu; // C2H4 number of moles
real n_x_3_exit_accu; // CH4 number of moles
real n_x_4_exit_accu; // CO number of moles
real n_x_5_exit_accu; // CO2 number of moles
real n_x_6_exit_accu; // O2 number of moles
real n_exit_total_accu; // total number of moles exit the battery (moles)

// Rate of ejected moles
real n_x_1_exit_t; // H2 rate of exit moles (mole/sec)
real n_x_2_exit_t; // C2H4 rate of exit moles (mole/sec)
real n_x_3_exit_t; // CH4 rate of exit moles (mole/sec)
real n_x_4_exit_t; // CO rate of exit moles (mole/sec)
real n_x_5_exit_t; // CO2 rate of exit moles (mole/sec)
real n_x_6_exit_t; // O2 rate of exit moles (mole/sec)
real n_exit_t; // rate of exit (moles/sec)

//number of moles in the battery (generated - ejected)
real n_x_1; // H2 number of moles
real n_x_2; // C2H4 number of moles
real n_x_3; // CH4 number of moles
real n_x_4; // CO number of moles
real n_x_5; // CO2 number of moles
real n_x_6; // O2 number of moles
real n_total_inside;     // total number of moles inside the battery (accumulated - exit)
real n_exit_inis; // instaneous gas vent
real R_gas = 8.314;  // Universal gas constant (J/mol.K) (?)

//Mass Flux rates
real m_x_1; // H2 mass flux(Kg/m2 sec)
real m_x_2; // C2H4 mass flux(Kg/m2 sec)
real m_x_3; // CH4 mass flux(Kg/m2 sec)
real m_x_4; // CO mass flux(Kg/m2 sec)
real m_x_5; // CO2 mass flux(Kg/m2 sec)
real m_x_6; // O2 mass flux(Kg/m2 sec)
real m_total;     // total mass flux(Kg/m2 sec)

// Variables for the third submodel
real P_avg;  // Average internal pressure
real V_void;  // void volume
real rho_mix;  // void volume
real gama_H2 = 1.4;
real gama_C2H4 = 1.19;
real gama_CH4 = 1.3;
real gama_CO = 1.4;
real gama_CO2 = 1.29;
real gama_O2 = 1.39;
real gama;
real Ma;
real P2;
real T2;
real P3;
real u2;
real u3;
real TINY = 1e-30;
// Function to define initial values
DEFINE_INIT(init_concentrations, domain)
{
    Thread* t;
    cell_t c;

    thread_loop_c(t, domain)
    {
        if (THREAD_ID(t) == 6)
        {
            begin_c_loop(c, t)
            {
                // Initialize UDM slots with starting concentrations
                C_UDMI(c, t, 0) = 1;        // c_chem_old initial
                C_UDMI(c, t, 1) = 1;        // c_elec_old initial
                C_UDMI(c, t, 2) = 2.2e5;    // c_x_1 initial
                C_UDMI(c, t, 3) = 1.02e3;   // c_x_2 initial
                C_UDMI(c, t, 4) = 5.6e3;    // c_x_3 initial
                C_UDMI(c, t, 5) = 8e3;      // c_x_4 initial
                C_UDMI(c, t, 6) = 1.02e3;   // c_x_5 initial
                C_UDMI(c, t, 7) = 7.6e2;    // c_x_6 initial
                C_UDMI(c, t, 8) = 5.6e3;    // c_x_7 initial
                C_UDMI(c, t, 9) = 1e3;      // c_x_8 initial
                C_UDMI(c, t, 10) = 2.4e3;   // c_x_9 initial
                C_UDMI(c, t, 11) = 0.12;    // c_x_10 initial
                C_UDMI(c, t, 12) = 0.12;    // c_x_11 initial
                P_avg = P_atm;
                V_void = V_void_init;
                n_exit_t = 0.0;
                n_exit_total_accu = 0.0;
                n_x_1_exit_accu = 0.0;
                n_x_2_exit_accu = 0.0;
                n_x_3_exit_accu = 0.0;
                n_x_4_exit_accu = 0.0;
                n_x_5_exit_accu = 0.0;
                n_x_6_exit_accu = 0.0;
                n_x_1 = 0.0;
                n_x_2 = 0.0;
                n_x_3 = 0.0;
                n_x_4 = 0.0;
                n_x_5 = 0.0;
                n_x_6 = 0.0;
            } end_c_loop(c, t)
        }
    }
}

// Function to Compute Average Temperature in Battery Zone (2D)
DEFINE_EXECUTE_AT_END(average_temp) 
{
	Domain* d;
	cell_t c;
	real area = 0.0;
	real area_tot = 0.0;
	real T;
	int ID = 6;  // Set the battery zone ID
	Thread* t;
	int zone_ID;
	d = Get_Domain(1);
	// Get thread corresponding to battery zone
	t = Lookup_Thread(d, ID);
	T_avg = 0.0;  // Reset T_avg for the new calculation
	begin_c_loop(c, t)
	{
		area = C_VOLUME(c, t);
		T = C_T(c, t);   // Get the temperature from the cell
		area_tot += area;
		T_avg += T * area;
	}
	end_c_loop(c, t)
	T_avg /= area_tot; // Compute area-weighted avg temp
	printf("1.T_avg = %g area_tot = %g\n", T_avg, area_tot);
	begin_c_loop(c, t)
	{
		C_T(c, t) = T_avg;   // Return Average temp to the cell
        C_UDMI(c, t, 18) = T_avg; // Saving the average temperature in memory
	}
	end_c_loop(c, t)
}

// UDF to define heat generation source term
DEFINE_SOURCE(my_battery_heat_source,c,t,ds,eqn)
{
    // Declaring Local variable for chemical and electrical energy concentration
    real C_chem_old = C_UDMI(c,t,0);  // UDM[0] for C_chem
    real C_elec_old = C_UDMI(c,t,1);  // UDM[1] for C_elec
    real C_chem;
    real C_elec;
    real dC_chem_dt = 0.0;
    real dC_elec_dt = 0.0;
    real q_chem;
    real q_elec;
    real qv;
    real A_x, E_a_x, n_x1, n_x2;
    real g_x;  // chemical reaction correction term
    real c_x_1_stoich = 2.2e5;
    real c_x_2_stoich = 1.02e3;
    real c_x_3_stoich = 5.6e3;
    real c_x_4_stoich = 8e3;
    // declaring reaction rates
    real R_x_1 = 0.0; // chemical reaction one (mol/s m3)
    real R_x_2 = 0.0; // chemical reaction two
    real R_x_3 = 0.0; // chemical reaction three
    real R_x_4 = 0.0; // chemical reaction four
    real R_x_5 = 0.0; // chemical reaction five
    real R_x_6 = 0.0; // chemical reaction six
    real R_x_7 = 0.0; // chemical reaction seven
    real R_x_8 = 0.0; // chemical reaction eight
    real R_x_9 = 0.0; // chemical reaction nine
    real R_x_10 = 0.0; // chemical reaction ten
    real R_x_11 = 0.0; // chemical reaction eleven
    // declaring UDMs for concentrations
    real c_x_1_old = C_UDMI(c, t, 2);   // UDM[2] for concentration of equation 1
    real c_x_2_old = C_UDMI(c, t, 3);   // UDM[3] for concentration of equation 2
    real c_x_3_old = C_UDMI(c, t, 4);   // UDM[4] for concentration of equation 3
    real c_x_4_old = C_UDMI(c, t, 5);   // UDM[5] for concentration of equation 4
    real c_x_5_old = C_UDMI(c, t, 6);   // UDM[6] for concentration of equation 5
    real c_x_6_old = C_UDMI(c, t, 7);   // UDM[7] for concentration of equation 6
    real c_x_7_old = C_UDMI(c, t, 8);   // UDM[8] for concentration of equation 7
    real c_x_8_old = C_UDMI(c, t, 9);   // UDM[9] for concentration of equation 8
    real c_x_9_old = C_UDMI(c, t, 10);  // UDM[10] for concentration of equation 9
    real c_x_10_old = C_UDMI(c, t, 11); // UDM[11] for concentration of equation 10
    real c_x_11_old = C_UDMI(c, t, 12); // UDM[12] for concentration of equation 11

    real c_x_1; // (mol/m3)
    real c_x_2;
    real c_x_3;
    real c_x_4;
    real c_x_5;
    real c_x_6;
    real c_x_7;
    real c_x_8;
    real c_x_9;
    real c_x_10;
    real c_x_11;
	real T_avg = C_T(c, t);   // Return Average temp from the cell
    real dt = CURRENT_TIMESTEP;  // Get the simulation time step
    // Function to compute the gas generation rates
    if (T_avg > 413.15)
    {
        // H2 and CH4 Gas Generation Rate = reaction one and three
        R_x_1 = (1.917e25 * pow(c_x_1_old, 1) * pow(1 - c_x_1_old, 0) * exp(-2.86e5 / (R_gas * T_avg)) * exp(-c_x_1_old / c_x_1_stoich));
        c_x_1 = c_x_1_old - (R_x_1 * dt);
        if (c_x_1 < .00001)
        {
            c_x_1 = 0.0;
        }
        R_x_3 = (4e28 * pow(c_x_3_old, 1) * pow(1 - c_x_3_old, 0) * exp(-4.7e5 / (R_gas * T_avg)) * exp(-c_x_3_old / c_x_3_stoich));
        c_x_3 = c_x_3_old - (R_x_3 * dt);
        if (c_x_3 < .00001)
        {
            c_x_3 = 0.0;
        }
        // C2H4 Gas Generation Rate = reation two and nine
        R_x_2 = (1.95e20 * pow(c_x_2_old, 1) * pow(1 - c_x_2_old, 0) * exp(-2e5 / (R_gas * T_avg)) * exp(-c_x_2_old / c_x_2_stoich));
        c_x_2 = c_x_2_old - (R_x_2 * dt);
        if (c_x_2 < .00001)
        {
            c_x_2 = 0.0;
        }
        // CO Gas Generation Rate
        R_x_4 = (4e28 * pow(c_x_4_old, 1) * pow(1 - c_x_4_old, 0) * exp(-4.7e5 / (R_gas * T_avg)) * exp(-c_x_4_old / c_x_4_stoich));
        c_x_4 = c_x_4_old - (2 * R_x_4 * dt);
        if (c_x_4 < .00001)
        {
            c_x_4 = 0.0;
        }
        // CO2 Gas Generation Rate
        R_x_5 = (5.14e25 * pow(c_x_5_old, 1) * pow(1 - c_x_5_old, 0) * exp(-2.74e5 / (R_gas * T_avg)) * 1);
        c_x_5 = c_x_5_old - (R_x_5 * dt);
        if (c_x_5 < .00001)
        {
            c_x_5 = 0.0;
        }
        R_x_6 = (5.14e25 * pow(c_x_6_old, 1) * pow(1 - c_x_6_old, 0) * exp(-2.74e5 / (R_gas * T_avg)) * 1);
        c_x_6 = c_x_6_old - (R_x_6 * dt);
        if (c_x_6 < .00001)
        {
            c_x_6 = 0.0;
        }
        R_x_7 = (5.14e25 * pow(c_x_7_old, 1) * pow(1 - c_x_7_old, 0) * exp(-2.74e5 / (R_gas * T_avg)) * 1);
        c_x_7 = c_x_7_old - (R_x_7 * dt);
        if (c_x_7 < .00001)
        {
            c_x_7 = 0.0;
        }
        R_x_8 = (5.14e25 * pow(c_x_8_old, 1) * pow(1 - c_x_8_old, 0) * exp(-2.74e5 / (R_gas * T_avg)) * 1);
        c_x_8 = c_x_8_old - (R_x_8 * dt);
        if (c_x_8 < .00001)
        {
            c_x_8 = 0.0;
        }
        // O2 Gas Generation Rate
        R_x_10 = (1.75e9 * pow(c_x_10_old, 1) * pow(1 - c_x_10_old, 0) * exp(-1.1495e5 / (R_gas * T_avg)) * 1);
        c_x_10 = c_x_10_old - (R_x_10 * dt);
        if (c_x_10 < .00001)
        {
            c_x_10 = 0.0;
        }
        R_x_11 = (1.077e12 * pow(c_x_11_old, 1) * pow(1 - c_x_11_old, 0) * exp(-1.5888e5 / (R_gas * T_avg)) * 1);
        c_x_11 = c_x_11_old - (R_x_11 * dt);
        if (c_x_11 < .00001)
        {
            c_x_11 = 0.0;
        }
        C_UDMI(c, t, 2) = c_x_1;
        C_UDMI(c, t, 3) = c_x_2;
        C_UDMI(c, t, 4) = c_x_3;
        C_UDMI(c, t, 5) = c_x_4;
        C_UDMI(c, t, 6) = c_x_5;
        C_UDMI(c, t, 7) = c_x_6;
        C_UDMI(c, t, 8) = c_x_7;
        C_UDMI(c, t, 9) = c_x_8;
        C_UDMI(c, t, 11) = c_x_10;
        C_UDMI(c, t, 12) = c_x_11;
    }
    if (T_avg > 343.15 && T_avg <413.15)
    {
        R_x_9 = (7.88e36 * pow(c_x_9_old, 1) * pow(1 - c_x_9_old, 0) * exp(-2.81e5 / (R_gas * T_avg)) * 1);
        c_x_9 = c_x_9_old - (R_x_9 * dt);
        if (c_x_9 < .00001)
        {
            c_x_9 = 0.0;
        }
        C_UDMI(c, t, 10) = c_x_9;
    }
 // Store the new concentration
   G_x_1 = (0.5 * R_x_1) - (1 * R_x_3); // Generation rate for H2
   if (G_x_1 < 0.0)
    {
      G_x_1 = 0.0;
    }
   G_x_2 = R_x_2 + R_x_9; // Generation rate for C2H4
    if (G_x_2 < 0.0)
    {
      G_x_2 = 0.0;
    }
   G_x_3 = 2 * R_x_3; // Generation rate for CH4
    if (G_x_3 < 0.0)
    {
      G_x_3 = 0.0;
    }
   G_x_4 = R_x_4; // Generation rate for CO
    if (G_x_4 < 0.0)
    {
      G_x_4 = 0.0;
    }
   G_x_5 = (3 * R_x_5) + (4 * R_x_6) + (3 * R_x_7) + R_x_8 + R_x_9 - (2 * R_x_4);
    if (G_x_5 < 0.0)
    {
      G_x_5 = 0.0;
    }
   G_x_6 = R_x_10 + R_x_11 + (0.5 * R_x_9) - (2.5 * R_x_5) - (4 * R_x_6) - (3 * R_x_7);
    if (G_x_6 < 0.0)
    {
      G_x_6 = 0.0;
    }
 // Calculate the total gas generation rate and total number of moles
   G_total = G_x_1 + G_x_2 + G_x_3 + G_x_4 + G_x_5 + G_x_6; // Total generated species(mol/m3 sec)
   n_x_1_gen_accu += G_x_1 * dt * V_bat; // Accumulated number of moles for H2 (moles)
   n_x_2_gen_accu += G_x_2 * dt * V_bat; // Accumulated number of moles for C2H4 (moles)
   n_x_3_gen_accu += G_x_3 * dt * V_bat; // Accumulated number of moles for CH4 (moles)
   n_x_4_gen_accu += G_x_4 * dt * V_bat; // Accumulated number of moles for CO (moles)
   n_x_5_gen_accu += G_x_5 * dt * V_bat; // Accumulated number of moles for CO2 (moles)
   n_x_6_gen_accu += G_x_6 * dt * V_bat; // Accumulated number of moles for O2 (moles)
   n_total_gen_accu += G_total * dt * V_bat; // Accumulated number of moles for all species (moles)

 // Instantaneous
   n_x_1_gen = G_x_1 * dt * V_bat; // number of moles for H2 (moles)
   n_x_2_gen = G_x_2 * dt * V_bat; // number of moles for C2H4 (moles)
   n_x_3_gen = G_x_3 * dt * V_bat; // number of moles for CH4 (moles)
   n_x_4_gen = G_x_4 * dt * V_bat; // number of moles for CO (moles)
   n_x_5_gen = G_x_5 * dt * V_bat; // number of moles for CO2 (moles)
   n_x_6_gen = G_x_6 * dt * V_bat; // number of moles for O2 (moles)
   n_total_gen = G_total * dt * V_bat; // number of moles for all species (moles)

 // Calculations for the third submodel
   /*V_void = V_void_init + ((n_exit_total_accu / n_total_gen_accu) * d_void);*/
   if (n_total_gen_accu > TINY) 
   {
       V_void = V_void_init + ((n_exit_total_accu / n_total_gen_accu) * d_void);
   }
   else {
       V_void = V_void_init;
   }
   P_avg = ((n_total_inside * R_gas * T_avg) / V_void) + P_atm;

   //Calculating number of moles inside the battery (generated - vent)
   n_total_inside = n_total_gen_accu - n_exit_total_accu; //(moles)
   if (n_total_inside < .0000000001)
   {
       n_total_inside = 0.0;
   }
   n_x_1 = n_x_1_gen_accu - n_x_1_exit_accu;
   if (n_x_1 < .0000000001)
   {
       n_x_1 = 0.0;
   }
   n_x_2 = n_x_2_gen_accu - n_x_2_exit_accu;
   if (n_x_2 < .0000000001)
   {
       n_x_2 = 0.0;
   }
   n_x_3 = n_x_3_gen_accu - n_x_3_exit_accu;
   if (n_x_3 < .0000000001)
   {
       n_x_3 = 0.0;
   }
   n_x_4 = n_x_4_gen_accu - n_x_4_exit_accu;
   if (n_x_4 < .0000000001)
   {
       n_x_4 = 0.0;
   }
   n_x_5 = n_x_5_gen_accu - n_x_5_exit_accu;
   if (n_x_5 < .0000000001)
   {
       n_x_5 = 0.0;
   }
   n_x_6 = n_x_6_gen_accu - n_x_6_exit_accu;
   if (n_x_6 < .0000000001)
   {
       n_x_6 = 0.0;
   }
   if ((P_atm / P_avg) > pow(2 / (gama + 1), gama / (gama - 1)))
   {
        Ma = sqrt((2 / (gama - 1)) * (pow(P_avg / P_atm, (gama - 1) / gama) - 1));
   }
   else
   {
       Ma = 1.0;
   }
   C_UDMI(c, t, 27) = Ma;
   if (P_avg > 1.5e6)
   {
       n_exit_t = (u3 * Area_valve * P2) / (R_gas * T2); // rate of ejected moles (moles /sec)
       n_exit_total_accu += n_exit_t * dt; // number of ejected moles (moles)(accumualted vent)
       n_exit_inis = n_exit_t * dt; // number of ejected moles (moles) (inistaneously)
       //Calculating number of moles inside the battery (generated - vent)
       n_x_1_exit_accu = (n_x_1_gen / n_total_gen) * n_exit_total_accu;
       n_x_2_exit_accu = (n_x_2_gen / n_total_gen) * n_exit_total_accu;
       n_x_3_exit_accu = (n_x_3_gen / n_total_gen) * n_exit_total_accu;
       n_x_4_exit_accu = (n_x_4_gen / n_total_gen) * n_exit_total_accu;
       n_x_5_exit_accu = (n_x_5_gen / n_total_gen) * n_exit_total_accu;
       n_x_6_exit_accu = (n_x_6_gen / n_total_gen) * n_exit_total_accu;
       n_x_1_exit_t = (n_x_1_gen / n_total_gen) * n_exit_t;  //(moles/sec)
       n_x_2_exit_t = (n_x_2_gen / n_total_gen) * n_exit_t;
       n_x_3_exit_t = (n_x_3_gen / n_total_gen) * n_exit_t;
       n_x_4_exit_t = (n_x_4_gen / n_total_gen) * n_exit_t;
       n_x_5_exit_t = (n_x_5_gen / n_total_gen) * n_exit_t;
       n_x_6_exit_t = (n_x_6_gen / n_total_gen) * n_exit_t;
   }
   else
   {
       n_exit_t = 0.0;
       n_exit_total_accu = 0.0;
       n_x_1_exit_accu = 0.0;
       n_x_2_exit_accu = 0.0;
       n_x_3_exit_accu = 0.0;
       n_x_4_exit_accu = 0.0;
       n_x_5_exit_accu = 0.0;
       n_x_6_exit_accu = 0.0;
   }

   //UDMS for measured quantities (Stores)
   C_UDMI(c, t, 13) = P_avg; // UDM[13] for internal average pressure within the battery
   C_UDMI(c, t, 14) = n_total_inside; // UDM[14] for total gas inside the battery (moles)
   C_UDMI(c, t, 15) = n_total_gen_accu; // UDM[15] for accumulated generated gas (moles)
   C_UDMI(c, t, 16) = n_exit_total_accu; // UDM[16] for accumualted exit gas (moles
   C_UDMI(c, t, 17) = V_void;
   C_UDMI(c, t, 19) = G_x_1; //UDM for gas generation rate H2 (mol/m3)
   C_UDMI(c, t, 20) = G_x_2; //UDM for gas generation rate C2H4 (mol/m3)
   C_UDMI(c, t, 21) = G_x_3; //UDM for gas generation rate CH4 (mol/m3)
   C_UDMI(c, t, 22) = G_x_4; //UDM for gas generation rate CO (mol/m3)
   C_UDMI(c, t, 23) = G_x_5; //UDM for gas generation rate CO2 (mol/m3)
   C_UDMI(c, t, 24) = G_x_6; //UDM for gas generation rate O2 (mol/m3)

   if (((n_x_1 * MW_H2) + (n_x_2 * MW_C2H4) + (n_x_3 * MW_CH4) + (n_x_4 * MW_CO) + (n_x_5 * MW_CO2) + (n_x_6 * MW_O2)) > 0)
   {
       gama = ((gama_H2 * n_x_1 * MW_H2) + (gama_C2H4 * n_x_2 * MW_C2H4) + (gama_CH4 * n_x_3 * MW_CH4) + (gama_CO * n_x_4 * MW_CO) + (gama_CO2 * n_x_5 * MW_CO2) + (gama_O2 * n_x_6 * MW_O2)) / ((n_x_1 * MW_H2) + (n_x_2 * MW_C2H4) + (n_x_3 * MW_CH4) + (n_x_4 * MW_CO) + (n_x_5 * MW_CO2) + (n_x_6 * MW_O2));
   }
   else
   {
       gama = 0.0;
   }
       C_UDMI(c, t, 25) = gama;
   rho_mix = (((n_x_1 * MW_H2) + (n_x_2 * MW_C2H4) + (n_x_3 * MW_CH4) + (n_x_4 * MW_CO) + (n_x_5 * MW_CO2) + (n_x_6 * MW_O2)) / V_void);
   C_UDMI(c, t, 26) = rho_mix;
   // Function to compute the temperature and velocity at the orifice
   T2 = T_avg / (1.0 + ((gama - 1.0) / 2.0) * Ma * Ma);
   u2 = Ma * sqrt(gama * R_gas * T2);
   if (Ma <= 1)
   {
       P2 = P_atm;
       u3 = u2 * 0.7;
   }
   else
   {
       P2 = P_avg / pow((1.0 + (gama - 1.0) / 2.0), (gama / (gama - 1)));
       P3 = P_atm;
       u3 = 0.7 * (u2 - ((P3 - P2) / rho_mix * u2));    
   }
   
    //Calculating Mass Flux
    m_x_1 = (n_x_1_exit_t * MW_H2) / (1000 * Area_valve); //(Kg/m2 sec)
    m_x_2 = (n_x_2_exit_t * MW_C2H4) / (1000 * Area_valve); //(Kg/m2 sec)
    m_x_3 = (n_x_3_exit_t * MW_CH4) / (1000 * Area_valve); //(Kg/m2 sec)
    m_x_4 = (n_x_4_exit_t * MW_CO) / (1000 * Area_valve); //(Kg/m2 sec)
    m_x_5 = (n_x_5_exit_t * MW_CO2) / (1000 * Area_valve); //(Kg/m2 sec)
    m_x_6 = (n_x_6_exit_t * MW_O2) / (1000 * Area_valve); //(Kg/m2 sec)


    
    
    
 // Function to compute the rate of change of C_chem (dC_chem/dt)
    if (T_avg > T_chem && T_avg <= T_tr && C_chem_old >= 0 && C_chem_old <= 1) {
        dC_chem_dt = (Cp / H_chem) * A * rho * V_bat * pow((T_avg / T_ref), b);
    }
    else if (T_avg > T_tr && C_chem_old >= 0 && C_chem_old <= 1) {
        dC_chem_dt = C_chem_max;
    }
    else {
        dC_chem_dt = 0.0;
    }
    /*printf("dC_chem_dt = %e", dC_chem_dt);*/
    // Update the chemical energy concentration
    C_chem = C_chem_old - (dC_chem_dt * dt);

    // Function to compute volumetric chemical heat generation rate (q_chem)
    q_chem = (1.0 / V_bat) * H_chem * dC_chem_dt;
    C_UDMI(c, t, 0) = C_chem;
    if (C_chem < 0.0)
    {
        C_chem_old = 0.0;
        dC_chem_dt = 0.0;
    }
    /*printf("C_chem = %e", C_chem);*/
    // Function to compute short-circuit reaction rate (dC_elec/dt)
    if (T_avg > T_tr && C_elec_old >= 0 && C_elec_old <= 1) {
        dC_elec_dt = C_elec_max;
    }
    else {
        dC_elec_dt = 0.0;
    }
    // Update the electrical energy concentration
    C_elec = C_elec_old - (dC_elec_dt * dt);
    // Function to compute volumetric electric heat generation rate (q_elec)
    q_elec = (1.0 / V_bat) * H_sc * dC_elec_dt;
    C_UDMI(c, t, 1) = C_elec;
    if (C_elec < 0.0)
    {
        C_elec_old = 0.0;
        dC_elec_dt = 0.0;
    }
    /*printf("C_elec = %e", C_elec);*/
    // Compute heat generation rates
    qv = q_chem + q_elec; // Total heat generation rate (W/m³)
    ds[eqn] = 0;
    /*printf("q_chem = %g q_elec = %g\n", q_chem, q_elec);*/
    return qv; // Return total heat source term for Fluent
}

// Define profiles for the boundry condition
//Define H2 Profile
DEFINE_PROFILE(mass_flux_H2, t, i)
{
    real x[ND_ND]; /* this will hold the position vector */
    real y;
    face_t f;
    begin_f_loop(f, t)
    {
        if (G_total > 0)
        {
            F_CENTROID(x, f, t);
            y = x[1];
            F_PROFILE(f, t, i) = m_x_1;
        }
        else
        {
            F_PROFILE(f, t, i) = 0.0;
        }
     }
    end_f_loop(f, t)
}
//Define C2H4 Profile
DEFINE_PROFILE(mass_flux_C2H4, t, i)
{
    real x[ND_ND]; /* this will hold the position vector */
    real y;
    face_t f;
    begin_f_loop(f, t)
    {
        if (G_total > 0)
        {
            F_CENTROID(x, f, t);
            y = x[1];
            F_PROFILE(f, t, i) = m_x_2;
        }
        else
        {
            F_PROFILE(f, t, i) = 0.0;
        }
     }
    end_f_loop(f, t)
}
//Define CH4 Profile
DEFINE_PROFILE(mass_flux_CH4, t, i)
{
    real x[ND_ND]; /* this will hold the position vector */
    real y;
    face_t f;
    begin_f_loop(f, t)
    {
        if (G_total > 0)
        {
            F_CENTROID(x, f, t);
            y = x[1];
            F_PROFILE(f, t, i) = m_x_3;
        }
        else
        {
            F_PROFILE(f, t, i) = 0.0;
        }
     }
    end_f_loop(f, t)
}

//Define CO Profile
DEFINE_PROFILE(mass_flux_CO, t, i)
{
    real x[ND_ND]; /* this will hold the position vector */
    real y;
    face_t f;
    begin_f_loop(f, t)
    {
        if (G_total > 0)
        {
            F_CENTROID(x, f, t);
            y = x[1];
            F_PROFILE(f, t, i) = m_x_4;
        }
        else
        {
            F_PROFILE(f, t, i) = 0.0;
        }
     }
    end_f_loop(f, t)
}

//Define CO2 Profile
DEFINE_PROFILE(mass_flux_CO2, t, i)
{
    real x[ND_ND]; /* this will hold the position vector */
    real y;
    face_t f;
    begin_f_loop(f, t)
    {
        if (G_total > 0)
        {
            F_CENTROID(x, f, t);
            y = x[1];
            F_PROFILE(f, t, i) = m_x_5;
        }
        else
        {
            F_PROFILE(f, t, i) = 0.0;
        }
     }
    end_f_loop(f, t)
}

//Define O2 Profile
DEFINE_PROFILE(mass_flux_O2, t, i)
{
    real x[ND_ND]; /* this will hold the position vector */
    real y;
    face_t f;
    begin_f_loop(f, t)
    {
        if (G_total > 0)
        {
            F_CENTROID(x, f, t);
            y = x[1];
            F_PROFILE(f, t, i) = m_x_6;
        }
        else
        {
            F_PROFILE(f, t, i) = 0.0;
        }
     }
    end_f_loop(f, t)
}

//Define T2 Profile
DEFINE_PROFILE(inlet_temperature, t, i)
{
    real x[ND_ND]; /* this will hold the position vector */
    real y;
    face_t f;
    begin_f_loop(f, t)
    {
        if (G_total > 0)
        {
            F_CENTROID(x, f, t);
            y = x[1];
            F_PROFILE(f, t, i) = T2;
        }
        else
        {
            F_PROFILE(f, t, i) = 0.0;
        }
     }
    end_f_loop(f, t)
}

//Define U3 Profile
DEFINE_PROFILE(inlet_velocity, t, i)
{
    real x[ND_ND]; /* this will hold the position vector */
    real y;
    face_t f;
    begin_f_loop(f, t)
    {
        if (G_total > 0)
        {
            F_CENTROID(x, f, t);
            y = x[1];
            F_PROFILE(f, t, i) = u3;
        }
        else
        {
            F_PROFILE(f, t, i) = 0.0;
        }
     }
    end_f_loop(f, t)
}