#define h 1.0e-4          //time step
#define c_p 500e-12       //capacitance of P cell
#define c_Ia 200e-12      //            of Ia cell
#define c_Ib 600e-12      //            of Ib cell
#define c_gl 45e-12       //            of glial cell
#define g_p 25e-9         //conductance of P cells
#define g_Ia 20e-9        //            of Ia cell
#define g_Ib 15e-9        //            of Ib cell
#define g_gl 9e-9         //            of glial cell
#define gh_AMPA 0.5e-9    //maximal conductance for AMPA receptor
#define gh_GABA 0.7e-9    //                    for GABA receptor
#define a_p 290e-12       //sensory input current
#define a_AMPA 1.1e+6     //channel opening rate for AMPA receptor
#define a_GABA 5.0e+6     //                     for GABA receptor
#define b_AMPA 190        //channel closing rate for AMPA receptor
#define b_GABA 180        //                     for GABA receptor
#define n_p 270           //steepness of sigmoid function for P cell
#define n_Ia 240          //                              for Ia cell
#define n_Ib 360          //                              for Ib cell
#define n_p_DMN 230       //                              for P (DMN) cell
#define n_Ib_DMN 320      //                              for Ib (DMN) cell
#define s_p -0.0335       //threshold of sigmoid function for P cell
#define s_Ia -0.0390      //                              for Ia cell
#define s_Ib -0.040       //                              for Ib cell
#define s_p_DMN -0.0330   //                              for P (DMN) cell
#define s_Ib_DMN -0.040   //                              for Ib (DMN) cell
#define gam 2.5           //decay constant for ambient GABA concentration 
#define T_GL 140e+7       //GABA transfer coefficient 
#define GABA_0 1.0e-6     //basal ambinet GABA concentration
#define GABA_max 1.5e-6   //maximal ambient GABA concentration
#define GABA_min 0.0      //minimal ambient GABA concentration
#define u_AMPA_rev 0.0    //reversal potential of AMPA receptor
#define u_GABA_rev -80e-3 //                   of GABA receptor
#define u_gl_rev -70e-3   //                   of glial transpoter
#define o_p 60.0e+1       //amount of extrasynaptic GABAa receptor
#define u_p_rest -65e-3   //resting potential of P cell
#define u_Ia_rest -70e-3  //                  of Ia cell
#define u_Ib_rest -70e-3  //                  of Ib cell
#define u_gl_rest -70e-3  //                  of glial cell
#define inp 3.0           //stimulus relevant cell assembly 
#define t_p 0.01          //broadness of input 
#define w_pp 6.5          //synaptic weight from P to P
#define w_pp_DMN 8.0      //synaptic weight from P (DMN) to P (DMN)
#define w_pIb 4.0         //synaptic weight from Ib to P
#define w_pIb_DMN 14.0    //synaptic weight from Ib (DMN) to P
#define w_Ibp 60.0        //synaptic weight from P to Ib
#define w_Ibp_DMN 4.50    //synaptic weight from P (DMN) to Ib
#define w_Ia 30.0         //synaptic weight from P to Ia
#define	w_gl_Ia 20.0      //synaptic weight from Ia to glia
#define w_Ia_DMN 0.0      //synaptic weight from P (DMN) to Ia/////strength of tonic excitation (#: 4.0)
#define	w_p_DMN_Nsen 0.0   //synaptic weight from P (DMN) to P//////strength of phasic excitation(*: 0.535)
#define inp_t 2.0         //input time
#define inp_t_len 0.5     //input time length
#define t_end 3.0         //simulation time
#define seed 5            //seed of random number
#define m 10              //cell unit number of output data (data.csv)
