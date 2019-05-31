#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include <string.h>
#include"MT19937.h"
#include"parameter.h"

double rung(double x,double y,double z,int n);
double equation(double x,double y,double z,int n);
struct cell{
		double mem;//////////////////////////////////////////////membrane potential
		double curs;/////////////////////////////////////////////current which cell recieve
		double r;////////////////////////////////////////////////fraction of receptors in open state
		double conc;/////////////////////////////////////////////concentration of transmitters
		int fire;////////////////////////////////////////////////state of cell (firing stete, resting stete or not)
};

void firing(struct cell *str,int n);
struct conection{
		double cur;//////////////////////////////////////////////current which flows via each connection
		double sum;//////////////////////////////////////////////summation of fraction of open state receptors
};
int main(void)
{
	FILE *fp;
	FILE *fp_ras;
	FILE *fp_p;
	int i,j,v,n,z,rast_plot;
	double t;
	struct conection I_pp[8][20], I_pp_DMN[20], I_pIb[8][20], I_pIb_DMN[20], I_pext[8][20], I_inp[8][20], I_Iap[8][20], I_Iap_DMN[8][20], I_Ibp[8][20], I_Ibp_DMN[20], I_gl_p[8][20], I_gl_Ia[8][20],I_P_sen_P_DMN[8][20];
	struct cell P_sen[8][20],Ib_sen[8][20],Ia[8][20],glia[8][20],P_DMN[20],Ib_DMN[20];
	fp=fopen("data.csv","w");
	fp_ras=fopen("raster_plot.dat","w");
	fp_p=fopen("up.csv","w");
	init_genrand(seed);
//////////////////////////////////////////////////////initalize all value 
	for(i=0;i<=7;i++)
	{
		for(j=0;j<=19;j++)
		{
			P_sen[i][j].mem=u_p_rest;
			Ib_sen[i][j].mem=u_Ib_rest;
			Ia[i][j].mem=u_Ia_rest;
			glia[i][j].mem=u_gl_rest;
			P_DMN[j].mem=u_p_rest;
			Ib_DMN[j].mem=u_Ib_rest;
			P_sen[i][j].r=0.0;
			Ib_sen[i][j].r=0.0;
			Ia[i][j].r=0.0;
			glia[i][j].r=0.0;
			P_DMN[j].r=0.0;
			Ib_DMN[j].r=0.0;
			P_sen[i][j].conc=0.0;
			Ib_sen[i][j].conc=0.0;
			Ia[i][j].conc=0.0;
			glia[i][j].conc=GABA_0;
			P_DMN[j].conc=0.0;
			Ib_DMN[j].conc=0.0;
			P_sen[i][j].fire=0;
			Ib_sen[i][j].fire=0;
			Ia[i][j].fire=0;
			glia[i][j].fire=0;
			P_DMN[j].fire=0;
			Ib_DMN[j].fire=0;
		}
	}
//////////////////////////////////////////////////////set label of each date
	for(i=0;i<=179;i++)
	{
		fprintf(fp_ras,"ras_wave_%d_%d,",i/20,i%20);
	}
	for(i=0;i<=159;i++)
	{
		fprintf(fp_p,"up_wave_%d_%d,",i/20,i%20);
	}
	fprintf(fp_ras,"\n");
	fprintf(fp_p,"\n");
	fprintf(fp,"up_DMN, ,uIb_DMN, ,up_sen0,up_sen1,up_sen2,up_sen3,up_sen4,up_sen5,up_sen6,up_sen7, ,uIa0,uIa1,uIa2,uIa3,uIa4,uIa5,uIa6,uIa7, ,uIb_sen0,uIb_sen1,uIb_sen2,uIb_sen3,uIb_sen4,uIb_sen5,uIb_sen6,uIb_sen7, ,ugl0,ugl1,ugl2,ugl3,ugl4,ugl5,ugl6,ugl7, ,GABA_ext0,GABA_ext1,GABA_ext2,GABA_ext3,GABA_ext4,GABA_ext5,GABA_ext6,GABA_ext7\n");
	for(t=0.0;t<=t_end;t+=h)
	{
//////////////////////////////////////////////////////output membrane potentials of P cell(DMN), Ib (DMN), P (Nsen), Ia, Ib (Nsen) and glia, and concentration of ambient-GABA
		fprintf(fp,"%.12f, ,%.12f, ,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f, ,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f, ,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f, ,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f, ,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f\n",
		P_DMN[m].mem,
		Ib_DMN[m].mem,
		P_sen[0][m].mem,P_sen[1][m].mem,P_sen[2][m].mem,P_sen[3][m].mem,P_sen[4][m].mem,P_sen[5][m].mem,P_sen[6][m].mem,P_sen[7][m].mem,
		Ia[0][m].mem,Ia[1][m].mem,Ia[2][m].mem,Ia[3][m].mem,Ia[4][m].mem,Ia[5][m].mem,Ia[6][m].mem,Ia[7][m].mem,
		Ib_sen[0][m].mem,Ib_sen[1][m].mem,Ib_sen[2][m].mem,Ib_sen[3][m].mem,Ib_sen[4][m].mem,Ib_sen[5][m].mem,Ib_sen[6][m].mem,Ib_sen[7][m].mem,
		glia[0][m].mem,glia[1][m].mem,glia[2][10].mem,glia[3][m].mem,glia[4][m].mem,glia[5][m].mem,glia[6][m].mem,glia[7][m].mem,
		glia[0][m].conc,glia[1][m].conc,glia[2][m].conc,glia[3][m].conc,glia[4][m].conc,glia[5][m].conc,glia[6][m].conc,glia[7][m].conc);
//////////////////////////////////////////////////////output raster plot and all membrane potentials of P cells
		rast_plot=0;
		for(i=0;i<=7;i++)
		{
			
			for(j=0;j<=19;j++)
			{
		    	rast_plot++;
				if(P_sen[i][j].fire==1)
				{
					fprintf(fp_ras,"%d,",rast_plot);
				}
				else
				{
					fprintf(fp_ras,"0,");
				}
				
				if(t>=0.5 && t<=(inp_t-h))
				{
					if(P_sen[i][j].fire>=1 && P_sen[i][j].fire<=10)
					{
						fprintf(fp_p,",");
					}
					else
					{
						fprintf(fp_p,"%.12f,",P_sen[i][j].mem);
					}
				}
			}
		}
		for(j=0;j<=19;j++)
		{
			rast_plot++;
			if(P_DMN[j].fire==1)
			{
				fprintf(fp_ras,"%d,",rast_plot);
			}
			else
			{
				fprintf(fp_ras,"0,");
			}
		}
		fprintf(fp_ras,"\n");
		if(t>=0.5 && t<=(inp_t-h))
		{
			fprintf(fp_p,"\n");
		}
//////////////////////////////////////////////////////calucuate sigma
		for(i=0;i<=7;i++)
		{
			for(j=0;j<=19;j++)	
			{
				I_pp[i][j].sum=0.0;
				I_pp_DMN[j].sum=0.0;
				I_pIb[i][j].sum=0.0;
				I_pIb_DMN[j].sum=0.0;
				I_Ibp[i][j].sum=0.0;
				I_Ibp_DMN[j].sum=0.0;
				I_Iap_DMN[i][j].sum=0.0;
				I_P_sen_P_DMN[i][j].sum=0.0;
			}
		}
		for(i=0;i<=7;i++)
		{
			for(j=0;j<=19;j++)
			{
				for(v=0;v<=19;v++)
				{
					if(v != j)
					{
						I_pp[i][j].sum += w_pp * P_sen[i][v].r;
						if(i==0)
						{
							I_pp_DMN[j].sum += w_pp_DMN * P_DMN[v].r;
						}	
					}
					I_Ibp_DMN[j].sum += w_Ibp_DMN * P_sen[i][v].r;
					I_pIb[i][j].sum += w_pIb * Ib_sen[i][v].r;
					if(i==0)
					{
						I_pIb_DMN[j].sum += w_pIb_DMN * Ib_DMN[j].r;
					}
				}
				for(v=0;v<=7;v++)
				{
					if(v != i)
					{
						I_Ibp[i][j].sum += w_Ibp * P_sen[v][j].r;
					}
				}
				if(i==0)
				{
					I_Iap_DMN[0][0].sum += w_Ia_DMN * P_DMN[j].r;
					I_P_sen_P_DMN[0][0].sum += w_p_DMN_Nsen * P_DMN[j].r;
				}
			}
		}
//////////////////////////////////////////////////////solve current and membrane potential equaton
		for(i=0;i<=7;i++)
		{
			for(j=0;j<=19;j++)
			{
				I_pp[i][j].cur=-gh_AMPA*(P_sen[i][j].mem-u_AMPA_rev)*I_pp[i][j].sum;
				I_pIb[i][j].cur=-gh_GABA*(P_sen[i][j].mem-u_GABA_rev)*I_pIb[i][j].sum;
				I_pext[i][j].cur=-gh_GABA*(P_sen[i][j].mem-u_GABA_rev)*o_p*glia[i][j].r;
				if(t<=inp_t || t>=inp_t+inp_t_len )
				{
					I_inp[i][j].cur = 0.0;
				}
				else
				{
					I_inp[i][j].cur = a_p*exp(-1*((i-inp)/t_p)*((i-inp)/t_p));
				}
				I_Iap[i][j].cur = -gh_AMPA*(Ia[i][j].mem-u_AMPA_rev)*w_Ia * P_sen[i][j].r;
				I_Ibp[i][j].cur = -gh_AMPA*(Ib_sen[i][j].mem-u_AMPA_rev)*I_Ibp[i][j].sum;
				I_gl_Ia[i][j].cur = -gh_GABA*(glia[i][j].mem-u_GABA_rev)*w_gl_Ia*Ia[i][j].r;
				I_Iap_DMN[i][j].cur = -gh_AMPA*(Ia[i][j].mem-u_AMPA_rev)*I_Iap_DMN[0][0].sum;
				I_P_sen_P_DMN[i][j].cur = -gh_AMPA*(P_sen[i][j].mem-u_AMPA_rev)*I_P_sen_P_DMN[0][0].sum;
				
				P_sen[i][j].curs = I_pp[i][j].cur + I_pIb[i][j].cur + I_pext[i][j].cur + I_inp[i][j].cur + I_P_sen_P_DMN[i][j].cur;
				Ia[i][j].curs = I_Iap[i][j].cur + I_Iap_DMN[i][j].cur;
				Ib_sen[i][j].curs = I_Ibp[i][j].cur;
				glia[i][j].curs = I_gl_Ia[i][j].cur;
				
	    		P_sen[i][j].mem += rung(t,P_sen[i][j].mem,P_sen[i][j].curs,0);
			    Ia[i][j].mem += rung(t,Ia[i][j].mem,Ia[i][j].curs,1);
			    Ib_sen[i][j].mem += rung(t,Ib_sen[i][j].mem,Ib_sen[i][j].curs,2);
			    glia[i][j].mem += rung(t,glia[i][j].mem,glia[i][j].curs,3);
				
				firing(&P_sen[i][j],0);
				firing(&Ia[i][j],1);
				firing(&Ib_sen[i][j],2);
				if(i==0)
				{
					I_pp_DMN[j].cur=-gh_AMPA*(P_DMN[j].mem-u_AMPA_rev)*I_pp_DMN[j].sum;
					I_pIb_DMN[j].cur=-gh_GABA*(P_DMN[j].mem-u_GABA_rev)*I_pIb_DMN[j].sum;
					I_Ibp_DMN[j].cur=-gh_AMPA*(Ib_DMN[j].mem-u_AMPA_rev)*I_Ibp_DMN[j].sum;
			
					P_DMN[j].curs = I_pp_DMN[j].cur + I_pIb_DMN[j].cur;
					Ib_DMN[j].curs= I_Ibp_DMN[j].cur;
			
					P_DMN[j].mem += rung(t,P_DMN[j].mem,P_DMN[j].curs,0);
		   			Ib_DMN[j].mem += rung(t,Ib_DMN[j].mem,Ib_DMN[j].curs,2);
			
					firing(&P_DMN[j],3);
					firing(&Ib_DMN[j],4);
				}
			}		
		}
//////////////////////////////////////////////////////solve equation of open state receptor's fraction and ambient-GABA connentration
		for(i=0;i<=7;i++)
		{
			for(j=0;j<=19;j++)
			{
		    	P_sen[i][j].r += rung(t,P_sen[i][j].r,P_sen[i][j].conc,4);
		    	Ia[i][j].r += rung(t,Ia[i][j].r,Ia[i][j].conc,5);
				Ib_sen[i][j].r += rung(t,Ib_sen[i][j].r,Ib_sen[i][j].conc,5);
				glia[i][j].conc += rung(t,glia[i][j].conc,glia[i][j].mem,6);
			    glia[i][j].r+= rung(t,glia[i][j].r,glia[i][j].conc,7);
				if(i==0)
				{
				    P_DMN[j].r += rung(t,P_DMN[j].r,P_DMN[j].conc,4);
					Ib_DMN[j].r+= rung(t,Ib_DMN[j].r,Ib_DMN[j].conc,5);
				}
			}
		}
	}
	fclose(fp);		
	fclose(fp_ras);
	fclose(fp_p);
}

double rung(double x,double y,double z, int n)///////////////////Runge-Kutta algorithm
{
	double dy;
	double k[4];
	k[0]=equation(x,y,z,n)*h;
	k[1]=equation(x+h/2,y+k[0]/2,z,n)*h;
	k[2]=equation(x+h/2,y+k[1]/2,z,n)*h;
	k[3]=equation(x+h,y+k[2],z,n)*h;
	dy=(k[0]+2*k[1]+2*k[2]+k[3])/6;
	
	return dy;
}

double equation(double x,double y,double z,int n)////////////////differential equations
{
	switch(n)
	{
		case 0:
			return (-g_p*(y-u_p_rest)+z)/c_p;
			break;
		case 1:
			return (-g_Ia*(y-u_Ia_rest)+z)/c_Ia;
			break;
		case 2:
			return (-g_Ib*(y-u_Ib_rest)+z)/c_Ib;
			break;
		case 3:
			return (-g_gl*(y-u_gl_rest)+z)/c_gl;
			break;
		case 4:
			return a_AMPA*z*(1-y)-b_AMPA*y;
			break;
		case 5:
			return a_GABA*z*(1-y)-b_GABA*y;
			break;
		case 6:
			return -gam*(y-GABA_0)+T_GL*(GABA_max-y)*(y-GABA_min)*(z-u_gl_rev);
			break;
		case 7:
			return a_GABA*z*(1-y)-b_GABA*y;
			break;
	}
}

void firing(struct cell *str,int n)/////////////////////////////generator of action potentials
{
	double ran,Prob,rest;
	switch(n)
	{
		case 0:
			rest=u_p_rest;
			Prob=1.0/(1.0+exp(-1.0*n_p*(str->mem-s_p)));
			break;
		case 1:
			rest=u_Ia_rest;
			Prob=1.0/(1.0+exp(-1.0*n_Ia*(str->mem-s_Ia)));
			break;
		case 2:
			rest=u_Ib_rest;
			Prob=1.0/(1.0+exp(-1.0*n_Ib*(str->mem-s_Ib)));
			break;
		case 3:
			rest=u_p_rest;
			Prob=1.0/(1.0+exp(-1.0*n_p_DMN*(str->mem-s_p_DMN)));
			break;
		case 4:
			rest=u_Ib_rest;
			Prob=1.0/(1.0+exp(-1.0*n_Ib_DMN*(str->mem-s_Ib_DMN)));
			break;
	}
	if(str->fire==0)
	{
		ran=genrand_real2();
		if(Prob >= ran)
		{
				str->fire++;
				str->conc=1.0e-3;
				str->mem=-10e-3;
		}
	}
	else
	{
		if(str->fire==20)
		{
			str->fire=0;
			str->conc=0.0;
		}
		else
		{
			if(str->fire<10)
			{
				str->conc=1e-3;
				str->mem=-10e-3;
				str->fire++;
			}
			else
			{
				str->conc=0.0;
				str->mem=rest;
				str->fire++;
			}
		}
	}
}
