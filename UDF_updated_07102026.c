/*my_battery_heat_source_trial.c*/
#include "udf.h"
#include "stdio.h"
#include "math.h"
#include "mem.h"
#include "time.h"
#include "sg_pollut.h"

#define BATTERY_ZONE_ID        6   // battery cell zone ID


// Inputs for First Submodel
#define V_bat 0.0000165  // Battery volume (m^3)
#define H_chem 11740.0     // Total heat release of chemical reaction (J) (?)
#define H_sc 17560.0       // Total heat release of short circuit (J) (?)
#define rho 2600.0         // Density of battery (kg/m^3)
#define Cp 1720.0          // Specific heat capacity (J/kg.K)
#define A_chem 0.92           // Pre-exponential factor (unitless)
#define b_chem 28.5           // Reaction rate constant (unitless)
#define T_ref 533.15     // Normalized temperature coefficient (K) (?)
#define T_chem 413.15    // Start temperature of chemical reactions (K)
#define T_tr 533.15      // Triggering temperature of thermal runaway (K)
#define C_chem_max 12    
#define C_elec_max 12    

// Inputs for third submodel

#define MW_H2 2.016       // Molar mass of H2 (g/mol)
#define MW_C2H4 28.05     // Molar mass of C2H4(g/mol)
#define MW_CH4 16.04      // Molar mass of CH4(g/mol)
#define MW_CO 28.01      // Molar mass of CO(g/mol)
#define MW_CO2 44.01      // Molar mass of CO2(g/mol)
#define MW_O2 32.0      // Molar mass of O2(g/mol)

#define c_x_1_stoich 2.2e5
#define c_x_2_stoich 1.02e3
#define c_x_3_stoich 5.6e3
#define c_x_4_stoich 8e3

#define P_atm 101325.0      // Aymospheric pressure (Pa)
#define P_vent 1.5e6       // vent opening pressure (Pa) 
#define V_void_init 0.00000152      // Initial void volume (m3)
#define d_void 0.000007594      // volume increase (m3)
#define Area_valve 0.00003848 // Safety valve area (dia 7mm) (m^2)
#define R_gas        8.314          // universal gas constant (J/mol K)
#define TINY_VAL     1e-30
#define m_total 9.0  // 9 gram total mass decomposed

int printed = 0;
real T_avg = 300.0;     // average battery temperature (K)
real q_stored = 0.0;    // heat source used by DEFINE_SOURCE (W/m3)

real C_chem, C_elec;    // remaining chemical / electrical energy fraction
real c_x[11];           // concentrations c_x_1 ... c_x_11 (mol/m3) -> c_x[0] ... c_x[10]

/* Species index used in all arrays below:
   0 = H2, 1 = C2H4, 2 = CH4, 3 = CO, 4 = CO2, 5 = O2 */
real MW[6]   = {MW_H2, MW_C2H4, MW_CH4, MW_CO, MW_CO2, MW_O2};
real gama_sp[6] = {1.4, 1.19, 1.3, 1.4, 1.29, 1.39};

real G_x[6], G_total;                  // gas generation rate (mol/m3 s)
real n_gen_accu[6], n_total_gen_accu;  // generated moles, accumulated
real n_exit_accu[6], n_exit_total_accu;// vented moles, accumulated
real n_exit_rate[6], n_exit_t;         // vent rate (mol/s)
real n_inside[6], n_total_inside;      // moles inside the battery
real m_x[6];                           // vent mass flux per species (kg/m2 s)
real m_tot_gen_accum;                      // accumulated mass generated
real m_x_total = 0;
real valve_open = 0;
real P_avg, V_void, rho_mix, gama, Ma, T2, u2, u3, P2;



/* ======================================================================
   1) INITIALIZATION: runs once when Initialize.
      Resets ALL global variables.
   ====================================================================== */
DEFINE_INIT(init_concentrations, domain)
{
    int k;
    Thread *t;
    cell_t c;

    // starting values of the battery state
    C_chem = 1.0;  // c_chem_old initial
    C_elec = 1.0;  // c_elec_old initial
    c_x[0] = 2.2e5;  // c_x_1 initial
    c_x[1] = 1.02e3;  // c_x_2 initial
    c_x[2] = 5.6e3;  // c_x_3 initial
    c_x[3] = 8e3;   // c_x_4 initial
    c_x[4] = 1.02e3;  // c_x_5 initial
    c_x[5] = 7.6e2; // c_x_6 initial
    c_x[6] = 5.6e3;  // c_x_7 initial
    c_x[7] = 1e3;  // c_x_8 initial
    c_x[8] = 2.4e3;  // c_x_9 initial
    c_x[9] = 0.12;    // c_x_10 initial
    c_x[10] = 0.12;  // c_x_11 initial

    for (k = 0; k < 6; k++)
    {
        G_x[k] = 0.0;
        n_gen_accu[k] = 0.0;
        n_exit_accu[k] = 0.0;
        n_exit_rate[k] = 0.0;
        n_inside[k] = 0.0;
        m_x[k] = 0.0;
    }
    G_total = 0.0;
    n_total_gen_accu = 0.0;
    n_exit_total_accu = 0.0;
    n_exit_t = 0.0;
    n_total_inside = 0.0;
    m_tot_gen_accum = 0.0;
    P_avg = P_atm;
    V_void = V_void_init;
    rho_mix = 0.0;  
    gama = 0.0;  
    Ma = 0.0;
    T2 = 300.0;  
    u2 = 0.0;  
    u3 = 0.0;  
    P2 = 0.0;
    q_stored = 0.0;
    valve_open = 0.0;
    m_x_total = 0;

#if !RP_HOST
    t = Lookup_Thread(domain, BATTERY_ZONE_ID);
    begin_c_loop(c, t)
    {
        C_UDMI(c, t, 0) = C_chem;
        C_UDMI(c, t, 1) = C_elec;
        for (k = 0; k < 11; k++) 
        {
            C_UDMI(c, t, 2 + k) = c_x[k];
        }
    }
    end_c_loop(c, t)
#endif

}

/* ======================================================================
   2) HEAT SOURCE: Fluent calls this for every battery cell in every
      iteration. It ONLY returns the stored value.
   ====================================================================== */
DEFINE_SOURCE(battery_heat_source, c, t, ds, eqn)
{
    ds[eqn] = 0.0;
    return q_stored;    // W/m3
}


/* ======================================================================
   3) END OF TIME STEP: runs ONCE per time step.
      a) average battery temperature
      b) battery model advanced by one time step
      c) every battery cell set to T_avg
      d) results saved in UDMs for plotting
   ====================================================================== */
DEFINE_EXECUTE_AT_END(average_temp)
{

#if !RP_HOST

    Domain *d = Get_Domain(1);
    Thread *t = Lookup_Thread(d, BATTERY_ZONE_ID);
    cell_t c;
    int k;
    real dt = CURRENT_TIMESTEP;
    real T;
    real vol_sum = 0.0, TV_sum = 0.0;
    real R[11];                     // reaction rates R_x_1 ... R_x_11
    real n_gen[6], n_total_gen;     // moles generated in this time step
    real mass, gama_sum, P3;
    real dC_chem_dt, dC_elec_dt, q_chem, q_elec;

    // ---------- a) average temperature ----------
    begin_c_loop_int(c, t)
    {
        TV_sum  += C_T(c, t) * C_VOLUME(c, t);
        vol_sum += C_VOLUME(c, t);
    }
    end_c_loop_int(c, t)
    TV_sum  = PRF_GRSUM1(TV_sum);  
    vol_sum = PRF_GRSUM1(vol_sum);
    T_avg = TV_sum / vol_sum;
    T = T_avg;

    // ---------- b) battery model per time step ----------

    /* gas-generating reactions */
    for (k = 0; k < 11; k++) R[k] = 0.0;

    if (T > 413.15 && m_tot_gen_accum < m_total)
    {
         // H2 and CH4 Gas Generation Rate = reaction one and three
        R[0] = 1.917e25 * c_x[0] * exp(-2.86e5 / (R_gas * T)) * exp(-c_x[0] / c_x_1_stoich);
        if (R[0] * dt > c_x[0]) R[0] = c_x[0] / (dt);
        c_x[0] = c_x[0] - R[0] * dt;
        if (c_x[0] < .00001) c_x[0] = 0.0;

        R[2] = 4e28 * c_x[2] * exp(-4.7e5 / (R_gas * T)) * exp(-c_x[2] / c_x_3_stoich);
        if (R[2] * dt > c_x[2]) R[2] = c_x[2] / (dt);
        c_x[2] = c_x[2] - R[2] * dt;
        if (c_x[2] < .00001) c_x[2] = 0.0;

        // C2H4 Gas Generation Rate = reation two and nine
        R[1] = 1.95e20 * c_x[1] * exp(-2e5 / (R_gas * T)) * exp(-c_x[1] / c_x_2_stoich);
        if (R[1] * dt > c_x[1]) R[1] = c_x[1] / (dt);
        c_x[1] = c_x[1] - R[1] * dt;
        if (c_x[1] < .00001) c_x[1] = 0.0;

        // CO Gas Generation Rate
        R[3] = 4e28 * c_x[3] * exp(-4.7e5 / (R_gas * T)) * exp(-c_x[3] / c_x_4_stoich);
        if (R[3] * dt > c_x[3]) R[3] = c_x[3] / (dt);
        c_x[3] = c_x[3] - R[3] * dt;
        if (c_x[3] < .00001) c_x[3] = 0.0;

        // CO2 Gas Generation Rate
        for (k = 4; k <= 7; k++)    // reactions 5 to 8 (same constants)
        {
            R[k] = 5.14e25 * c_x[k] * exp(-2.74e5 / (R_gas * T));
            if (R[k] * dt > c_x[k]) R[k] = c_x[k] / (dt);
            c_x[k] = c_x[k] - R[k] * dt;
            if (c_x[k] < .00001) c_x[k] = 0.0;
        }

        // O2 Gas Generation Rate
        R[9] = 1.75e9 * c_x[9] * exp(-1.1495e5 / (R_gas * T));
        if (R[9] * dt > c_x[9]) R[9] = c_x[9] / (dt);
        c_x[9] = c_x[9] - R[9] * dt;
        if (c_x[9] < .00001) c_x[9] = 0.0;

        R[10] = 1.077e12 * c_x[10] * exp(-1.5888e5 / (R_gas * T));
        if (R[10] * dt > c_x[10]) R[10] = c_x[10] / (dt);
        c_x[10] = c_x[10] - R[10] * dt;
        if (c_x[10] < .00001) c_x[10] = 0.0;
        
        R[8] = 7.88e36 * c_x[8] * exp(-2.81e5 / (R_gas * T));
        if (R[8] * dt > c_x[8]) R[8] = c_x[8] / dt;      // NEW: limit to what is left
        c_x[8] = c_x[8] - R[8] * dt;
        if (c_x[8] < .00001) c_x[8] = 0.0;
    }


    /* generation rate of each gas species */
    if (m_tot_gen_accum < m_total)
    {
        G_x[0] = 0.5 * R[0] - R[2];                                           // H2
        G_x[1] = R[1] + R[8];                                                 // C2H4
        G_x[2] = 2 * R[2];                                                    // CH4
        G_x[3] = R[3];                                                        // CO
        G_x[4] = 3 * R[4] + 4 * R[5] + 3 * R[6] + R[7] + R[8] - 2 * R[3];     // CO2
        G_x[5] = R[9] + R[10] + 0.5 * R[8] - 2.5 * R[4] - 4 * R[5] - 3 * R[6]; // O2
    }
    else {
        G_x[0] = 0.0;                 // H2
        G_x[1] = 0.0;              // C2H4
        G_x[2] = 0.0;          // CH4
        G_x[3] = 0.0;       // CO
        G_x[4] = 0.0;     // CO2
        G_x[5] = 0.0;   // O2
    }


    G_total = 0.0;
    for (k = 0; k < 6; k++)
    {
        if (G_x[k] < 0.0) G_x[k] = 0.0;
        G_total += G_x[k];
    }

    /* moles generated: this step, and accumulated */
    for (k = 0; k < 6; k++)
    {
        n_gen[k] = G_x[k] * dt * V_bat;
        n_gen_accu[k] += n_gen[k];
    }
    n_total_gen = G_total * dt * V_bat;
    n_total_gen_accu += n_total_gen;
    for (k = 0; k < 6; k++)
    { 
        m_tot_gen_accum += n_gen[k] * MW[k];
    
    }

    // void volume and internal pressure 
    if (m_x_total <= 0)
    {V_void = V_void_init;}
    else {
        if (n_total_gen_accu > TINY_VAL)
            V_void = V_void_init + (m_x_total / m_total) * d_void;
        else
            V_void = V_void_init;
    }
    

    // moles inside the battery = generated - vented
    n_total_inside = n_total_gen_accu - n_exit_total_accu;
    if (n_total_inside < 1e-10) n_total_inside = 0.0;
    for (k = 0; k < 6; k++)
    {
        n_inside[k] = n_gen_accu[k] - n_exit_accu[k];
        if (n_inside[k] < 1e-10) n_inside[k] = 0.0;
    }

    if (n_total_inside > 0 && V_void > 0)
    {P_avg = (n_total_inside * R_gas * T) / V_void + P_atm;}
    else {
        P_avg = P_atm;
    }

    if (P_avg >= 1.5e6)
    {valve_open = 1;}

    // Mach number at the valve
    if (valve_open == 1 && P_avg > P_atm)
    {
        if ((P_atm / P_avg) > pow(2 / (gama + 1), gama / (gama - 1)))
            Ma = sqrt((2 / (gama - 1)) * (pow(P_avg / P_atm, (gama - 1) / gama) - 1));
        else
            Ma = 1.0;
    }


    // mixture gamma and density
    mass = 0.0; // current mass inside the battery (gram)
    gama_sum = 0.0;
    //n_inside_total = 0.0;
    for (k = 0; k < 6; k++)
    {
        mass     += n_inside[k] * MW[k];
        gama_sum += gama_sp[k] * n_inside[k] * MW[k];
        /* n_inside_total += n_inside[k] */
    }
    if (mass > 0.0) gama = gama_sum / mass;
    else          gama = 1.4;
    rho_mix = (mass/1000) / V_void;

    // temperature and velocity at the valve
    T2 = T / (1.0 + ((gama - 1.0) / 2.0) * Ma * Ma);
    if (mass>0 && n_total_inside>0)
    u2 = Ma * sqrt(gama * ((R_gas*1000)/(mass/n_total_inside)) * T2);
    else u2 = 0;

    if (Ma < 1)
    {
        P2 = P_atm;
        u3 = u2 * 0.7;
    }
    else
    {
        P2 = P_avg / pow((1.0 + (gama - 1.0) / 2.0), (gama / (gama - 1)));
        P3 = P_atm;
        u3 = 0.7 * (u2 - ((P3 - P2) / (rho_mix * u2)));
    }

       // venting (opens above 1.5 MPa)
       // n_exit_accu[k] = (n_gen[k] / n_total_gen) * n_exit_total_accu;
    if (valve_open == 1 && T2 > 0 && u2 > 0 && n_total_inside > 0)
    {
        n_exit_t = (u3 * Area_valve * P2) / (R_gas * T2);
        if (n_exit_t * dt > n_total_inside) n_exit_t = n_total_inside / dt;  // can't vent more than inside
        n_exit_total_accu += n_exit_t * dt;
        
        for (k = 0; k < 6; k++)
        {
            n_exit_rate[k] = (n_gen[k] / n_total_gen) * n_exit_t;
            n_exit_accu[k] += n_exit_rate[k] * dt; // for mixture gamma and density calculation
        }
    }
    else
    {
        n_exit_t = 0.0;
        for (k = 0; k < 6; k++) n_exit_rate[k] = 0.0;    // stop the flow, but KEEP the totals
    }


    // vent mass flux of each species (kg/m2 s)
    if (valve_open == 1 && m_x_total <= m_total)
    {
        for (k = 0; k < 6; k++)
        {   m_x[k] = (n_exit_rate[k] * MW[k]) / (1000 * Area_valve);
            m_x_total += n_exit_rate[k] * MW[k] * dt; // gram
        }
    }

    
    // chemical heat
    if (T > T_chem && T <= T_tr && C_chem > 0 && C_chem <= 1)
        dC_chem_dt = (Cp / H_chem) * A_chem * rho * V_bat * pow(T / T_ref, b_chem);
    else if (T > T_tr && C_chem > 0 && C_chem <= 1)
        dC_chem_dt = C_chem_max;
    else
        dC_chem_dt = 0.0;
    if (C_chem < 0) C_chem = 0;
    if (dC_chem_dt * dt > C_chem) dC_chem_dt = C_chem / dt;    // NEW: can't release more than is left
    C_chem = C_chem - dC_chem_dt * dt;
    q_chem = (1.0 / V_bat) * H_chem * dC_chem_dt;

    // electrical (short-circuit heat)
    if (T > T_tr && C_elec >= 0 && C_elec <= 1)
        dC_elec_dt = C_elec_max;
    else
        dC_elec_dt = 0.0;
    if (dC_elec_dt * dt > C_elec) dC_elec_dt = C_elec / dt;    // NEW
    C_elec = C_elec - dC_elec_dt * dt;
    q_elec = (1.0 / V_bat) * H_sc * dC_elec_dt;

    // heat used by DEFINE_SOURCE during the NEXT time step
    q_stored = q_chem + q_elec;

    /* ---------- c and d set cells to T_avg, save results in UDM ---------- */
    /* TEST: what is in UDM 13 of one battery cell before we write it? */
    begin_c_loop_int(c, t)
    {
        if (!printed)
        {
            Message("check: UDM13 before write = %g   P_avg now = %g\n", C_UDMI(c, t, 13), P_avg);
            printed = 1;
        }
    }
    end_c_loop_int(c, t)

    begin_c_loop(c, t)
    {
        C_T(c, t) = T_avg;                 // lumped battery temperature

        C_UDMI(c, t, 0) = C_chem;
        C_UDMI(c, t, 1) = C_elec;
        for (k = 0; k < 11; k++) C_UDMI(c, t, 2 + k) = c_x[k];
        C_UDMI(c, t, 13) = P_avg;
        C_UDMI(c, t, 14) = n_total_inside;
        C_UDMI(c, t, 15) = n_total_gen_accu;
        C_UDMI(c, t, 16) = n_exit_total_accu;
        C_UDMI(c, t, 17) = V_void;
        C_UDMI(c, t, 18) = T_avg;
        for (k = 0; k < 6; k++) C_UDMI(c, t, 19 + k) = G_x[k];
        C_UDMI(c, t, 25) = gama;
        C_UDMI(c, t, 26) = rho_mix;
        C_UDMI(c, t, 27) = Ma;
        C_UDMI(c, t, 28) = m_tot_gen_accum;
    }
    end_c_loop(c, t)

    if (m_tot_gen_accum >= m_total)
    {
        Message("9 gram limit reached");
    }

   Message("t = %g s  T = %g K  P = %g Pa  n_in = %g mol  m_gen = %g g  q = %g W/m3\n",
           CURRENT_TIME, T_avg, P_avg, n_total_inside, m_tot_gen_accum, q_stored);
#endif

}



/* ======================================================================
   4) BOUNDARY PROFILES
   ====================================================================== */
DEFINE_PROFILE(mass_flux_H2, t, i)
{
    face_t f;
    begin_f_loop(f, t)
    {
        if (n_exit_t > 0) F_PROFILE(f, t, i) = m_x[0];
        else             F_PROFILE(f, t, i) = 0.0;
    }
    end_f_loop(f, t)
}

DEFINE_PROFILE(mass_flux_C2H4, t, i)
{
    face_t f;
    begin_f_loop(f, t)
    {
        if (n_exit_t > 0) F_PROFILE(f, t, i) = m_x[1];
        else             F_PROFILE(f, t, i) = 0.0;
    }
    end_f_loop(f, t)
}

DEFINE_PROFILE(mass_flux_CH4, t, i)
{
    face_t f;
    begin_f_loop(f, t)
    {
        if (n_exit_t > 0) F_PROFILE(f, t, i) = m_x[2];
        else             F_PROFILE(f, t, i) = 0.0;
    }
    end_f_loop(f, t)
}

DEFINE_PROFILE(mass_flux_CO, t, i)
{
    face_t f;
    begin_f_loop(f, t)
    {
        if (n_exit_t > 0) F_PROFILE(f, t, i) = m_x[3];
        else             F_PROFILE(f, t, i) = 0.0;
    }
    end_f_loop(f, t)
}

DEFINE_PROFILE(mass_flux_CO2, t, i)
{
    face_t f;
    begin_f_loop(f, t)
    {
        if (n_exit_t > 0) F_PROFILE(f, t, i) = m_x[4];
        else             F_PROFILE(f, t, i) = 0.0;
    }
    end_f_loop(f, t)
}

DEFINE_PROFILE(mass_flux_O2, t, i)
{
    face_t f;
    begin_f_loop(f, t)
    {
        if (n_exit_t > 0) F_PROFILE(f, t, i) = m_x[5];
        else             F_PROFILE(f, t, i) = 0.0;
    }
    end_f_loop(f, t)
}

DEFINE_PROFILE(inlet_temperature, t, i)
{
    face_t f;
    begin_f_loop(f, t)
    {
        if (n_exit_t > 0) F_PROFILE(f, t, i) = T2;
        else             F_PROFILE(f, t, i) = 0.0;   // 0 K when closed: original, fix in step 2
    }
    end_f_loop(f, t)
}

DEFINE_PROFILE(inlet_velocity, t, i)
{
    face_t f;
    begin_f_loop(f, t)
    {
        if (n_exit_t > 0) F_PROFILE(f, t, i) = u3;
        else             F_PROFILE(f, t, i) = 0.0;
    }
    end_f_loop(f, t)
}

  DEFINE_PROFILE(heater_flux, t, i)
  {
      face_t f;
      begin_f_loop(f, t)
      {
          if (T_avg < 533.15) F_PROFILE(f, t, i) = 1e5;   // W/m2
          else                F_PROFILE(f, t, i) = 0.0;
      }
      end_f_loop(f, t)
  }
