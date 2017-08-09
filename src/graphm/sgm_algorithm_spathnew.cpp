/***************************************************************************
 *   Copyright (C) 2008 by Mikhail Zaslavskiy   *
 *   mikhail.zaslavskiy@ensmp.fr   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "algorithm_path.h"
#include "sgm_algorithm_spathnew.h"

sgm_algorithm_spathnew::sgm_algorithm_spathnew(int num_seeds) :sgm_algorithm(num_seeds)
 {

}
match_result sgm_algorithm_spathnew::match(graph &g,graph &h,gsl_matrix* gm_P_i, gsl_matrix* gm_ldh,double dalpha_ldh){
 return(this->match_with_seeds(g, h ,gm_P_i, gm_ldh, dalpha_ldh, 0))	;
}

match_result sgm_algorithm_spathnew::match_with_seeds(graph& g, graph& h, gsl_matrix* gm_P_i, gsl_matrix* gm_ldh,double dalpha_ldh, unsigned int m_seeds)
{
	//int m  = get_param_i("num_seeds");
	//this->m = m_seeds;

    bool bblast_match=(get_param_i("blast_match")==1);
    bool bblast_match_end=(get_param_i("blast_match_proj")==1);
	bool bbest_path_proj=(get_param_i("best_path_proj_sol")==1);
	bool bbest_path_blast_proj=(get_param_i("best_path_blast_proj_sol")==1);
	bool bbest_path_greedy=(get_param_i("best_path_greedy_sol")==1);
	bool bbest_path_blast_greedy=(get_param_i("best_path_blast_greedy_sol")==1);
	bool bbest_path=bbest_path_proj or bbest_path_blast_proj or bbest_path_greedy or bbest_path_blast_greedy;
	std::ofstream fverbose;
	double dfw_xeps=get_param_d("algo_fw_xeps");
	double dfw_feps=get_param_d("algo_fw_feps");
	double dlambda_M=get_param_d("qcvqcc_lambda_M");
	double dlambda_min=get_param_d("qcvqcc_lambda_min");
	bool bverbose=(get_param_i("verbose_mode")==1);
	std::string verbose_fname = get_param_s("verbose_file");
	double dhung_max=get_param_d("hungarian_max");
  double bgreedy=(get_param_i("hungarian_greedy")==1);

	if (bverbose)
		*gout<<"Seeded Path matching"<<std::endl;
	if (bverbose) {

		*gout<<"this_m: "<< this->m <<std::endl; //

	}
	//some duplicate variables
	gsl_matrix* gm_Ag_d=g.get_descmatrix(cdesc_matrix);
	gsl_matrix* gm_Ah_d=h.get_descmatrix(cdesc_matrix);
	if (pdebug.ivalue) gsl_matrix_printout(gm_Ag_d,"Ag",pdebug.strvalue);
	if (pdebug.ivalue) gsl_matrix_printout(gm_Ah_d,"Ah",pdebug.strvalue);
	//laplacian construction
	gsl_matrix* gm_Lg_d=gsl_matrix_alloc(N,N);
	gsl_matrix* gm_Lh_d=gsl_matrix_alloc(N,N);
	gsl_matrix_memcpy(gm_Lg_d,gm_Ag_d);
	gsl_matrix_memcpy(gm_Lh_d,gm_Ah_d);
	gsl_vector* gv_ones = gsl_vector_alloc(N);
	gsl_vector* gv_res = gsl_vector_alloc(N);
	gsl_vector_set_all(gv_ones,1);
	gsl_blas_dgemv(CblasNoTrans,1,gm_Lg_d,gv_ones,0,gv_res);
	gsl_matrix_scale(gm_Lg_d,-1);
	for (int i=0;i<N;i++)
		gsl_matrix_set(gm_Lg_d,i,i,gsl_vector_get(gv_res,i));

	gsl_blas_dgemv(CblasNoTrans,1,gm_Lh_d,gv_ones,0,gv_res);
	gsl_matrix_scale(gm_Lh_d,-1);
	for (int i=0;i<N;i++)
		gsl_matrix_set(gm_Lh_d,i,i,gsl_vector_get(gv_res,i));

	gsl_vector_free(gv_res);
	gsl_vector_free(gv_ones);
	//Delta matrix
	gsl_matrix* gm_Delta=gsl_matrix_alloc(N,N);
	gsl_matrix_set_zero(gm_Delta);
	for (int i=m;i<N;i++)
		for (int j=m;j<N;j++){
			gsl_matrix_set(gm_Delta,i,j,pow(gsl_matrix_get(gm_Lg_d,i,i)-gsl_matrix_get(gm_Lh_d,j,j),2));
		};

	gsl_matrix_transpose(gm_Delta);
	gsl_matrix_scale(gm_Delta,1-dalpha_ldh);

	gsl_matrix* gm_Ag_dc = gsl_matrix_alloc(N,N);
	gsl_matrix_memcpy(gm_Ag_dc, gm_Ag_d);
	gsl_matrix* gm_Ah_dc = gsl_matrix_alloc(N,N);
        gsl_matrix_memcpy(gm_Ah_dc, gm_Ah_d);

        gsl_matrix* gm_Lg_dc=gsl_matrix_alloc(N,N);
        gsl_matrix* gm_Lh_dc=gsl_matrix_alloc(N,N);
        gsl_matrix_memcpy(gm_Lg_dc,gm_Lg_d);
        gsl_matrix_memcpy(gm_Lh_dc,gm_Lh_d);

        gsl_matrix* gm_Deltac = gsl_matrix_alloc(N,N);
        gsl_matrix_memcpy(gm_Deltac, gm_Delta);

	gsl_matrix_view gvi_Ag_d_sub11 = gsl_matrix_submatrix(gm_Ag_dc,0,0,m,m);
	gsl_matrix* gm_Ag_d_sub11 = &(gvi_Ag_d_sub11.matrix);

        gsl_matrix_view gvi_Ag_d_sub21t = gsl_matrix_submatrix(gm_Ag_dc,0,m,m,N-m);
        gsl_matrix* gm_Ag_d_sub21t = &(gvi_Ag_d_sub21t.matrix);

	gsl_matrix_view gvi_Ag_d_sub21 = gsl_matrix_submatrix(gm_Ag_dc,m,0,N-m,m);
	gsl_matrix* gm_Ag_d_sub21=&(gvi_Ag_d_sub21.matrix);

	gsl_matrix_view gvi_Ag_d_sub22=gsl_matrix_submatrix(gm_Ag_dc,m,m,N-m,N-m);
	gsl_matrix* gm_Ag_d_sub22=&(gvi_Ag_d_sub22.matrix);

        gsl_matrix_view gvi_Ah_d_sub11 = gsl_matrix_submatrix(gm_Ah_dc,0,0,m,m);
        gsl_matrix* gm_Ah_d_sub11 = &(gvi_Ah_d_sub11.matrix);

        gsl_matrix_view gvi_Ah_d_sub21t = gsl_matrix_submatrix(gm_Ah_dc,0,m,m,N-m);
        gsl_matrix* gm_Ah_d_sub21t = &(gvi_Ah_d_sub21t.matrix);

        gsl_matrix_view gvi_Ah_d_sub21 = gsl_matrix_submatrix(gm_Ah_dc,m,0,N-m,m);
        gsl_matrix* gm_Ah_d_sub21=&(gvi_Ah_d_sub21.matrix);

        gsl_matrix_view gvi_Ah_d_sub22=gsl_matrix_submatrix(gm_Ah_dc,m,m,N-m,N-m);
        gsl_matrix* gm_Ah_d_sub22=&(gvi_Ah_d_sub22.matrix);

        gsl_matrix_view gvi_Lg_d_sub11 = gsl_matrix_submatrix(gm_Lg_dc,0,0,m,m);
        gsl_matrix* gm_Lg_d_sub11 = &(gvi_Lg_d_sub11.matrix);

        gsl_matrix_view gvi_Lg_d_sub21t = gsl_matrix_submatrix(gm_Lg_dc,0,m,m,N-m);
        gsl_matrix* gm_Lg_d_sub21t = &(gvi_Lg_d_sub21t.matrix);

        gsl_matrix_view gvi_Lg_d_sub21 = gsl_matrix_submatrix(gm_Lg_dc,m,0,N-m,m);
        gsl_matrix* gm_Lg_d_sub21=&(gvi_Lg_d_sub21.matrix);

        gsl_matrix_view gvi_Lg_d_sub22=gsl_matrix_submatrix(gm_Lg_dc,m,m,N-m,N-m);
        gsl_matrix* gm_Lg_d_sub22=&(gvi_Lg_d_sub22.matrix);

	gsl_matrix_view gvi_Lh_d_sub11 = gsl_matrix_submatrix(gm_Lh_dc,0,0,m,m);
        gsl_matrix* gm_Lh_d_sub11 = &(gvi_Lh_d_sub11.matrix);

        gsl_matrix_view gvi_Lh_d_sub21t = gsl_matrix_submatrix(gm_Lh_dc,0,m,m,N-m);
        gsl_matrix* gm_Lh_d_sub21t = &(gvi_Lh_d_sub21t.matrix);

        gsl_matrix_view gvi_Lh_d_sub21 = gsl_matrix_submatrix(gm_Lh_dc,m,0,N-m,m);
        gsl_matrix* gm_Lh_d_sub21=&(gvi_Lh_d_sub21.matrix);

        gsl_matrix_view gvi_Lh_d_sub22=gsl_matrix_submatrix(gm_Lh_dc,m,m,N-m,N-m);
        gsl_matrix* gm_Lh_d_sub22=&(gvi_Lh_d_sub22.matrix);

        gsl_matrix_view gvi_Delta_sub11 = gsl_matrix_submatrix(gm_Deltac,0,0,m,m);
        gsl_matrix* gm_Delta_sub11 = &(gvi_Delta_sub11.matrix);

        gsl_matrix_view gvi_Delta_sub21t = gsl_matrix_submatrix(gm_Deltac,0,m,m,N-m);
        gsl_matrix* gm_Delta_sub21t = &(gvi_Delta_sub21t.matrix);

        gsl_matrix_view gvi_Delta_sub21 = gsl_matrix_submatrix(gm_Deltac,m,0,N-m,m);
        gsl_matrix* gm_Delta_sub21=&(gvi_Delta_sub21.matrix);

        gsl_matrix_view gvi_Delta_sub22=gsl_matrix_submatrix(gm_Deltac,m,m,N-m,N-m);
        gsl_matrix* gm_Delta_sub22=&(gvi_Delta_sub22.matrix);


	gm_Ag_d=gm_Ag_d_sub22;
	gm_Ah_d=gm_Ah_d_sub22;
	gm_Lg_d=gm_Lg_d_sub22;
	gm_Lh_d=gm_Lh_d_sub22;
	gm_Delta=gm_Delta_sub22;
	gsl_vector_view gvv_Delta=gsl_vector_view_array(gm_Delta->data,(N-m)*(N-m));

df_norm=pow(gsl_matrix_norm(gm_Ag_d,2),2)+pow(gsl_matrix_norm(gm_Ah_d,2),2);
df_norm=(df_norm>EPSILON)?df_norm:EPSILON;



	//memory allocation
	bool bstop_algo=false;
	gsl_vector* gv_C=gsl_vector_alloc((N-m)*(N-m));
	gsl_vector* gv_temp=gsl_vector_alloc((N-m)*(N-m));
	gsl_vector* gv_temp2=gsl_vector_alloc((N-m)*(N-m));
	gsl_matrix_view gmv_temp2=gsl_matrix_view_array(gv_temp2->data,(N-m),(N-m));
	gsl_matrix* gm_temp2=&gmv_temp2.matrix;
	gsl_matrix_view gmv_temp=gsl_matrix_view_array(gv_temp->data,(N-m),(N-m));
	gsl_matrix* gm_temp=&gmv_temp.matrix;
	gsl_matrix* C;
	gsl_matrix_view gmv_C;
	gsl_vector_view gvv_P,gvv_P_prev,gvv_dP,gvv_P_lambda;
	gsl_matrix* gm_P=gsl_matrix_alloc((N-m),(N-m));
	gsl_matrix* gm_P_prev=gsl_matrix_alloc((N-m),(N-m));
	gsl_matrix* gm_P_lambda=gsl_matrix_alloc((N-m),(N-m));
	gsl_matrix* gm_dP=gsl_matrix_alloc((N-m),(N-m));
	gsl_matrix_set_zero(gm_P_prev);

	if (gm_P_i==NULL){
		gsl_matrix_set_all(gm_P,1.0/(N-m));
		if (m > 0) {
			nonseededPtoseededP(gm_P,this->m);
		}
	}
	else {
		gsl_matrix_memcpy(gm_P,gm_P_i);
	    if (m > 0) {
			nonseededPtoseededP(gm_P,this->m);
		}
	}

        //perm matrix transformation into vector
	gvv_P=gsl_vector_view_array(gm_P->data,(N-m)*(N-m));
	gvv_P_prev=gsl_vector_view_array(gm_P_prev->data,(N-m)*(N-m));
	gvv_P_lambda=gsl_vector_view_array(gm_P_lambda->data,(N-m)*(N-m));
	gvv_dP=gsl_vector_view_array(gm_dP->data,(N-m)*(N-m));
	//and in opposite direction for gradient
	gmv_C=gsl_matrix_view_vector(gv_C,(N-m),(N-m));
	C=&gmv_C.matrix;
	gsl_vector*gv_debug_trace=gsl_vector_alloc(3);//debug trace information
	gsl_vector_set_zero(gv_debug_trace);
	//extern cycle over dlambda_cvcc
	//some temp variables
	gsl_matrix* gm_Delta_tmp=gsl_matrix_alloc((N-m),(N-m));
	gsl_matrix_memcpy(gm_Delta_tmp,gm_Delta);
	//temp variables to trace the best solution along the algorithm path
	gsl_matrix* gm_P_bp_temp=NULL;
	gsl_matrix* gm_P_bp=NULL;
	double fbest_path=1e+300;

	if (bbest_path)
	{
		gm_P_bp_temp=gsl_matrix_alloc((N-m),(N-m));
		gm_P_bp=gsl_matrix_alloc((N-m),(N-m));
	};
        //*************ALGORITHM***********************
	bool bpath_continue=true;
	double dlambda_add=dlambda_min;
	double dlambda_fix;
	double df_value_old,df_value,dP_norm,ddP_norm,dtemp;
	double dlambda=1;
	double dlambda_M_c=1;
	bool binc_lambda=false;
	bool bmax_lambda;
	gsl_matrix_memcpy(gm_P_prev,gm_P);
	if (pdebug.ivalue) gsl_matrix_printout(gm_Lg_d,"Lg",pdebug.strvalue);
	if (pdebug.ivalue) gsl_matrix_printout(gm_Lh_d,"Lh",pdebug.strvalue);
	if (pdebug.ivalue) gsl_matrix_printout(gm_Ag_d,"Ag",pdebug.strvalue);
	if (pdebug.ivalue) gsl_matrix_printout(gm_Ah_d,"Ah",pdebug.strvalue);
	if (pdebug.ivalue) gsl_matrix_printout(gm_P,"gm_P",pdebug.strvalue);
	double dt0,dt1,dt2,dt3;
	while (bpath_continue){
	//LINEAR COMBINATION
	//new step initialisation

	gsl_matrix_memcpy(gm_Delta,gm_Delta_tmp);

	gsl_matrix_scale(gm_Delta,1-dlambda);
	//remove df_norm
	gsl_matrix_scale(gm_Delta,1/df_norm);


	//linear term due to label similarities
	if (dalpha_ldh>0){
		gsl_matrix_transpose(gm_Delta);
		gsl_matrix_scale(gm_ldh,dalpha_ldh);
		gsl_matrix_add(gm_Delta,gm_ldh);
		gsl_matrix_scale(gm_ldh,(1.0/dalpha_ldh));
		gsl_matrix_transpose(gm_Delta);
	};
	//main optimization cycle
	int icounter=0;
	ddP_norm=1;
	df_value_old=f_qcvqcc(gm_Ag_d,gm_Ah_d,gm_Lg_d,gm_Lh_d,gm_Lg_d_sub21,gm_Lh_d_sub21t, gm_Delta,gm_P_prev,dlambda, gm_temp,gm_temp2, gm_Ag_d_sub21, gm_Ag_d_sub21t, gm_Ah_d_sub21,  gm_Ah_d_sub21t);
	while(!bstop_algo)
	{
	dt0=clock();
	//the default gsl representation is made by rows
	//gradient estimation: Agh*P, Here instead of P we have to use gvv_P_prev.
	if (pdebug.ivalue) gsl_matrix_printout(&gvv_P_prev.vector,"gvv_P_prev",pdebug.strvalue);
	sgm_algorithm_spathnew::qcvqcc_gradient(gm_Ag_d,gm_Ah_d,gm_Lg_d,gm_Lh_d,&gvv_P_prev.vector,gv_C,dlambda,gv_temp);
	if (pdebug.ivalue) gsl_matrix_printout(gv_C,"gv_C",pdebug.strvalue);
	//scaling without minus because before we have changed the Lghsign and Delta substraction
	//also we do not need to transpose Delta
	gsl_matrix_scale(C,2);
	gsl_matrix_sub(C,gm_Delta);// C is the gradient

	// Update the gradient in the seeded version.
	gsl_matrix* gm_temp5=gsl_matrix_alloc(N-m,N-m);
	gsl_matrix* gm_temp6=gsl_matrix_alloc(N-m,N-m);
	gsl_matrix* gm_temp7=gsl_matrix_alloc(m,N-m);
        gsl_matrix* gm_temp8=gsl_matrix_alloc(N-m,m);


 	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_Ag_d_sub21t,gm_P_prev,0,gm_temp7);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_Ag_d_sub21,gm_temp7,0,gm_temp6);
	gsl_matrix_scale(gm_temp6,2*dlambda/df_norm);
	gsl_matrix_add(C,gm_temp6);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_P_prev,gm_Ah_d_sub21,0,gm_temp8);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_temp8,gm_Ah_d_sub21t,0,gm_temp6);
	gsl_matrix_scale(gm_temp6,2*dlambda/df_norm);
        gsl_matrix_add(C,gm_temp6);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_Ag_d_sub21,gm_Ah_d_sub21t,0,gm_temp5);
	gsl_matrix_scale(gm_temp5,4*dlambda/df_norm);
        gsl_matrix_sub(C,gm_temp5);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_Lg_d_sub21,gm_Lh_d_sub21t,0,gm_temp5);
	gsl_matrix_scale(gm_temp5,2*(1-dlambda)/df_norm);
        gsl_matrix_sub(C,gm_temp5);

	//debug trace
	if (pdebug.ivalue) gsl_vector_set(gv_debug_trace,0,gsl_vector_get(gv_debug_trace,0)+1);
	if (pdebug.ivalue) gsl_vector_set(gv_debug_trace,1,gsl_matrix_norm(C,1));
	if (pdebug.ivalue) gsl_vector_set(gv_debug_trace,2,graph_dist(g,h,gm_P,cscore_matrix));
	if (pdebug.ivalue) gsl_matrix_printout(gv_debug_trace,"debug_trace",pdebug.strvalue);

	//result save
	gsl_matrix_transpose(C);
	if (pdebug.ivalue) gsl_matrix_printout(C,"C=gradient",pdebug.strvalue);
        //update_C_hungarian(C);
	double dscale_factor =gsl_matrix_max_abs(C);
	dscale_factor=(dscale_factor>EPSILON)?dscale_factor:EPSILON;
	dscale_factor=dhung_max/dscale_factor;
	gsl_matrix_scale(C,dscale_factor);

	if (pdebug.ivalue) gsl_matrix_printout(C,"scale(C)",pdebug.strvalue);
	//hungarian, before  C matrix must be transposed
	gsl_matrix_transpose(C);
	dt1=clock();
	// constrain gm_P to be seeded by solving the hungarian problem with the nonseeded portion of gradient.


	nonseededPtoseededP(gm_P, m);
	gsl_matrix_view gm_P_sub = gsl_matrix_submatrix(gm_P,m,m,N-m,N-m);



	gsl_matrix_hungarian(C, &gm_P_sub.matrix,NULL,NULL,false,(bblast_match?gm_ldh:NULL),bgreedy);
	dt2=clock();
	if (pdebug.ivalue) gsl_matrix_printout(&gm_P_sub.matrix,"gm_P_sub",pdebug.strvalue);
	if (pdebug.ivalue) gsl_matrix_printout(C,"C_sub",pdebug.strvalue);
	gsl_matrix_transpose(C);
	gsl_matrix_scale(C,1/dscale_factor);

	if (pdebug.ivalue) gsl_matrix_printout(gm_P,"gm_P",pdebug.strvalue);
	if (pdebug.ivalue) gsl_matrix_printout(gm_P_prev,"gm_P_prev",pdebug.strvalue);

	//line search
	gsl_matrix_memcpy(gm_dP,gm_P);
	gsl_matrix_sub(gm_dP,gm_P_prev);

	if (pdebug.ivalue) gsl_matrix_printout(gm_dP,"gm_dP",pdebug.strvalue);

	double a,b1,b2,b3,b4;
	//transpose all matrices but not Delta
	nonseededGradtoseededGrad(gm_dP,m);
	gsl_matrix_transpose(gm_dP);
	gsl_matrix_transpose(gm_P_prev);
	gsl_matrix_transpose(C);

	gsl_matrix_add(C,gm_Delta);//return to the original C for line search


	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,4*dlambda ,gm_Ag_d_sub21,gm_Ah_d_sub21t,0,gm_temp5);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,2*(1-dlambda),gm_Lg_d_sub21,gm_Lh_d_sub21t,0,gm_temp6);
	gsl_matrix_add(gm_temp5, gm_temp6);
	gsl_matrix_scale(gm_temp5, 1/df_norm);
	gsl_matrix_add(C, gm_temp5);



	gsl_vector_scale(gv_C,0.5);
	gsl_blas_ddot(&gvv_dP.vector,gv_C,&b1);

        gsl_matrix* gm_cb4=gsl_matrix_alloc((N-m),(N-m));
	gsl_matrix_memcpy(gm_cb4,gm_temp5);
        gsl_vector_view gvv_cb4=gsl_vector_view_array(gm_cb4->data,(N-m)*(N-m));
	gsl_blas_ddot(&gvv_dP.vector,&gvv_cb4.vector,&b4);
	gsl_blas_ddot(&gvv_dP.vector,&gvv_Delta.vector,&b3);
	b3 = b3 + b4;
	gsl_matrix_transpose(gm_dP);

	qcvqcc_gradient(gm_Ag_d,gm_Ah_d,gm_Lg_d,gm_Lh_d,&gvv_dP.vector,gv_temp,dlambda,gv_temp2);

        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_Ag_d_sub21t,gm_dP,0,gm_temp7);
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_Ag_d_sub21,gm_temp7,0,gm_temp6);
        gsl_matrix_scale(gm_temp6,dlambda/df_norm);
	gsl_matrix_add(gm_temp,gm_temp6);
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_dP,gm_Ah_d_sub21,0,gm_temp8);
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_temp8,gm_Ah_d_sub21t,0,gm_temp6);
        gsl_matrix_scale(gm_temp6,dlambda/df_norm);
        gsl_matrix_add(gm_temp,gm_temp6);


	if (pdebug.ivalue) gsl_matrix_printout(gm_temp,"gm_temp",pdebug.strvalue);
	gsl_matrix_transpose(gm_dP);
	//qcvqcc_gradient_sparse(spm_A,&gvv_dP.vector, gv_temp,v_x,v_res);
	gsl_blas_ddot(&gvv_P_prev.vector,gv_temp,&b2);
	b1 += b2-b3;
	gsl_blas_ddot(&gvv_dP.vector,gv_temp,&a);
	//transpose to the initial state
	gsl_matrix_transpose(gm_dP);
	gsl_matrix_transpose(gm_P_prev);
	gsl_matrix_transpose(C);
	df_value=sgm_algorithm_spathnew::f_qcvqcc(gm_Ag_d,gm_Ah_d,gm_Lg_d,gm_Lh_d,gm_Lg_d_sub21,gm_Lh_d_sub21t, gm_Delta,gm_P,dlambda,gm_temp,gm_temp2, gm_Ag_d_sub21, gm_Ag_d_sub21t, gm_Ah_d_sub21,  gm_Ah_d_sub21t);

	gsl_matrix *a_mat = gsl_matrix_alloc(1,1);
	gsl_matrix_set_all(a_mat,a);

	if (pdebug.ivalue) gsl_matrix_printout(a_mat,"a_quadratic",pdebug.strvalue);
	double alpha = -4;
	if (a>0)//projection is convex
	{
		alpha = -b1/(2*a);
		gsl_matrix *alpha_mat = gsl_matrix_alloc(1,1);
		gsl_matrix_set_all(alpha_mat,alpha);
		if (pdebug.ivalue) gsl_matrix_printout(alpha_mat,"alpha",pdebug.strvalue);
		gsl_matrix *b1_mat = gsl_matrix_alloc(1,1);
		gsl_matrix_set_all(b1_mat,b1);
		if (pdebug.ivalue) gsl_matrix_printout(b1_mat,"b_1",pdebug.strvalue);
		if (alpha >1) {
			alpha = 1.0 ;
		}
		if ((alpha<1) and (alpha>0))
			{
				gsl_matrix_scale(gm_dP,(1-alpha));
				gsl_matrix_sub(gm_P,gm_dP);
			};
		if (!(alpha>0))
			gsl_matrix_memcpy(gm_P,gm_P_prev);
	}
	//if (abs(a)<1e-50)
	//	if (!(b1>0))
	//		gsl_matrix_memcpy(gm_P,gm_P_prev);
	if ((a<0) or (abs(a)<1e-50)) //projection is concave
	{
		if (df_value>df_value_old)
			{gsl_matrix_memcpy(gm_P,gm_P_prev);
			 df_value=df_value_old;};
	};


	if (pdebug.ivalue) gsl_matrix_printout(gm_P,"gm_P_step_finish",pdebug.strvalue);


	//stop criterion
	dP_norm=gsl_matrix_norm(gm_P_prev,1);
	df_value=f_qcvqcc(gm_Ag_d,gm_Ah_d,gm_Lg_d,gm_Lh_d,gm_Lg_d_sub21,gm_Lh_d_sub21t, gm_Delta,gm_P,dlambda,gm_temp,gm_temp2, gm_Ag_d_sub21, gm_Ag_d_sub21t, gm_Ah_d_sub21,  gm_Ah_d_sub21t);

	gsl_matrix_memcpy(gm_temp,gm_P_prev);
	gsl_matrix_sub(gm_P_prev,gm_P);
	if (df_value>df_value_old)
		int dbg=1;
	ddP_norm=gsl_matrix_norm(gm_P_prev,1);
	bool bstop_algo_1=((ddP_norm<dfw_xeps*(N-m)*dlambda_M_c) and ((abs(df_value-df_value_old)<dfw_feps*abs(df_value_old)*dlambda_M_c) or (ddP_norm==0)));
	bool bstop_algo_2 = ((icounter>pow(N-m,4)));
	bool bstop_algo_3=( (abs(df_value-df_value_old)<1e-30));
	bool bstop_algo_4=(((icounter>0) and binc_lambda));//if we are on the increment step do not repreat many times
	bstop_algo = bstop_algo_1 or bstop_algo_2 or bstop_algo_3 or bstop_algo_4;
	icounter++;
	dt3=clock();
	long lnum_constraints=0;
	for (long li=0;li<(N-m)*(N-m);li++)
		 if (gm_P->data[li]<1e-30)
			{
			lnum_constraints++;
			};
	if (bverbose){
		*gout<<"iter="<<icounter<<", x="<<dP_norm<<", dx="<<ddP_norm<<", f="<<df_value<<", df="<<df_value_old-df_value<<", "<<"grad="<<(df_value_old-df_value)/ddP_norm<<". #Act.Constr="<<lnum_constraints<<". Timing="<<(dt1-dt0)/CLOCKS_PER_SEC<<" "<<(dt2-dt1)/CLOCKS_PER_SEC<<" "<<(dt3-dt2)/CLOCKS_PER_SEC<<std::endl;
		*gout<<"a="<<a<<";"<<"b1="<<b1<<";alpha="<<alpha<<std::endl;
	   if (bstop_algo_1) {
	   	*gout<<"stop because ((ddP_norm<dfw_xeps*N*dlambda_M_c) and ((abs(df_value-df_value_old)<dfw_feps*abs(df_value_old)*dlambda_M_c) or (ddP_norm==0)));"<<std::endl;
	   }
	   if (bstop_algo_2) {
	   	*gout<<"stop because ((icounter>pow(N,4)));"<<std::endl;
	   }
	   	if (bstop_algo_3) {

	   		*gout<<"stop because  (abs(df_value-df_value_old)<1e-30))"<<std::endl;
	   	}
	   	if (bstop_algo_4) {
	   		*gout<<"stop because ((icounter>0) and binc_lambda));"<<std::endl;
	   }
	   }
	df_value_old=df_value;

	//now we test different projection to estimate the best permutation
	if (bbest_path_proj)
	{
		//permuation projection
		gsl_matrix_transpose_memcpy(gm_P_bp_temp,gm_P);
		gsl_matrix_scale(gm_P_bp_temp,-10000);


		gsl_matrix_hungarian(gm_P_bp_temp,gm_P_prev,NULL,NULL,false,NULL,bgreedy);
		double df_bp_new=f_qcv1(gm_Ag_d,gm_Ah_d,gm_P_prev,gm_temp2,true);
		if (df_bp_new<fbest_path)
		{
			fbest_path=df_bp_new;
			gsl_matrix_memcpy(gm_P_bp,gm_P_prev);
		};
	};

	if (bbest_path_blast_proj)
	{
		//permuation projection
		gsl_matrix_transpose_memcpy(gm_P_bp_temp,gm_P);
		gsl_matrix_scale(gm_P_bp_temp,-10000);
		gsl_matrix_hungarian(gm_P_bp_temp,gm_P_prev,NULL,NULL,false,gm_ldh,bgreedy);
		double df_bp_new=f_qcv1(gm_Ag_d,gm_Ah_d,gm_P_prev,gm_temp2,true);
		if (df_bp_new<fbest_path)
		{
			fbest_path=df_bp_new;
			gsl_matrix_memcpy(gm_P_bp,gm_P_prev);
		};
	};

	if (bbest_path_greedy)
	{
	};

	if (bbest_path_blast_greedy)
	{
	};
	gsl_matrix_memcpy(gm_P_prev,gm_P);
	};//end Frank-Wolfe cycle
	bpath_continue=(dlambda>0);//stop if it is a concave function
	//we decrement lambda and repeat the Frank-Wolfe step
	//we are on the lambda-increment step
	if (!(dlambda>0) and (icounter>1) and (dlambda_add>2*dlambda_min))//if dlambda_add is very small then we have to stop
		{
			binc_lambda=true;
			bpath_continue=true;
			dlambda_M_c=dlambda_M;//for the next times
		};
	if (binc_lambda)
	{
		if ((bmax_lambda) and (icounter>1)) //there was two and more steps
		{
			dlambda_add/=2;//return to the previos step
			bmax_lambda=false;
		}
		else{
			if (icounter>1)	dlambda_add/=2;
			if (bmax_lambda) dlambda_add*=2;//continue lambda step
			//we change the current point if it is a limit step or if step is too small

		    };
		if (((!bmax_lambda) and (icounter==1)) or (dlambda_add<dlambda_min)){ binc_lambda=false;dlambda_M_c=1;//the first time we use the original stop criterion
			}; //just  change the lambda_fix
 		dlambda=dlambda_fix-dlambda_add;
		dlambda=(dlambda<0)?0:dlambda;
		if (dlambda==0){ //concave function
			binc_lambda=false;dlambda_M_c=1;
				};

		if (bverbose) *gout<<"lambda_fix="<<dlambda<<", dl_add="<<dlambda_add<<std::endl;
		gsl_matrix_memcpy(gm_P_prev,gm_P_lambda);
	}
	else
	{
		bmax_lambda=true;
		dlambda_fix=dlambda;
		gsl_matrix_memcpy(gm_P_lambda,gm_temp);
		binc_lambda=true;
		dlambda_M_c=dlambda_M;//for the next times
		bool beigen_trace=false;double dmin_eval=0;
		if (beigen_trace)
		{
		gsl_matrix* gm_Hessian=gsl_matrix_alloc((N-m)*(N-m),(N-m)*(N-m));
		gsl_matrix* gm_Hessian_2=gsl_matrix_alloc((N-m)*(N-m),(N-m)*(N-m));
		gsl_matrix* gm_Hessian_3=gsl_matrix_alloc((N-m)*(N-m),(N-m)*(N-m));

		for (int i1=0;i1<(N-m);i1++)
			for (int i2=0;i2<(N-m);i2++)
			   for (int j1=0;j1<(N-m);j1++)
				for (int j2=0;j2<(N-m);j2++)
			  {
				double dvalue=0;
				if (i2==j2)
					dvalue+=gm_Ag_d->data[i1*(N-m)+j1];
				if (i1==j1)
					dvalue-=gm_Ah_d->data[i2*(N-m)+j2];
				gsl_matrix_set(gm_Hessian,i1+i2*(N-m),j1+j2*(N-m),dvalue);
			  };
		gsl_matrix_transpose_memcpy(gm_Hessian_2,gm_Hessian);
		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_Hessian_2,gm_Hessian,0,gm_Hessian_3);
		gsl_matrix_transpose_memcpy(gm_Hessian,gm_Hessian_3);
		gsl_matrix_scale(gm_Hessian,dlambda);
		for (int i1=0;i1<(N-m);i1++)
			for (int i2=0;i2<(N-m);i2++)
			   for (int j1=0;j1<(N-m);j1++)
				for (int j2=0;j2<(N-m);j2++)
			  {
				double dvalue=0;
				dvalue+=-(1-dlambda)*gm_Lg_d->data[i1*(N-m)+j1]*gm_Lh_d->data[i2*(N-m)+j2];
				gsl_matrix_set(gm_Hessian,i1+i2*(N-m),j1+j2*(N-m),gsl_matrix_get(gm_Hessian,i1+i2*(N-m),j1+j2*(N-m))+dvalue);
			   };
		gsl_eigen_symm_workspace * gesw= gsl_eigen_symm_alloc ((N-m)*(N-m));
		gsl_vector* eval_h=gsl_vector_alloc((N-m)*(N-m));
		gsl_eigen_symm (gm_Hessian, eval_h,gesw);
		gsl_matrix_free(gm_Hessian);
		gsl_matrix_free(gm_Hessian_2);
		gsl_matrix_free(gm_Hessian_3);
		dmin_eval=gsl_vector_min(eval_h);
		gsl_eigen_symm_free (gesw);
		gsl_vector_free(eval_h);
		};
		if (bverbose) *gout<<"lambda_fix="<<dlambda_fix<<", dl_add="<<dlambda_add<<", min eval="<<dmin_eval<<std::endl;
	};
	if (dlambda_add<dlambda_min)
		dlambda_add=dlambda_min;
	if (bverbose) *gout<<"lambda="<<dlambda<<", dl_add="<<dlambda_add<<std::endl;
	bstop_algo=false;
	};//end pathway cycle
	match_result mres;
	mres.gm_P_exact=gsl_matrix_alloc((N-m),(N-m));
	gsl_matrix_memcpy(mres.gm_P_exact,gm_P);

	//permutation projection
	gsl_matrix_transpose_memcpy(gm_P_prev,gm_P);
	//nonseededPtoseededP(gm_P_prev,m);
	gsl_matrix_scale(gm_P_prev,-10000);

	//gsl_matrix_view gm_P_prev_sub = gsl_matrix_submatrix(gm_P_prev,m,m,N-m,N-m);
	//nonseededPtoseededP(gm_P,m);
	//gsl_matrix_view gm_P_sub = gsl_matrix_submatrix(gm_P,m,m,N-m,N-m);

	//gsl_matrix_hungarian(&gm_P_prev_sub.matrix,&gm_P_sub.matrix,NULL,NULL,false,(bblast_match_end?gm_ldh:NULL),bgreedy);
	gsl_matrix_hungarian(gm_P_prev,gm_P,NULL,NULL,false,(bblast_match_end?gm_ldh:NULL),bgreedy);
	if (bbest_path){
		double df_bp_new=f_qcv1(gm_Ag_d,gm_Ah_d,gm_P,gm_temp,true);
		if (df_bp_new>fbest_path)
		{
			*gout<<"Best path solution is used"<<std::endl;
			gsl_matrix_memcpy(gm_P,gm_P_bp);
		};
	};

	if (pdebug.ivalue) gsl_matrix_printout(gm_P,"gm_P_projected",pdebug.strvalue);
	//memory free
	gsl_matrix_free(gm_P_lambda);
	gsl_matrix_free(gm_P_prev);
	gsl_matrix_free(gm_dP);
	gsl_matrix_free(gm_Lg_d);
	gsl_matrix_free(gm_Lh_d);
	gsl_matrix_free(gm_Ag_d);
	gsl_matrix_free(gm_Ah_d);
	gsl_matrix_free(gm_Delta_tmp);
	gsl_matrix_free(gm_Delta);
	gsl_vector_free(gv_C);
	gsl_vector_free(gv_temp);
	gsl_vector_free(gv_temp2);
	gsl_vector_free(gv_debug_trace);
	if (bbest_path)
	{
		gsl_matrix_free(gm_P_bp_temp);
		gsl_matrix_free(gm_P_bp);
	};

	mres.gm_P=gm_P;

	//initial score
	mres.vd_trace.push_back(graph_dist(g,h,cscore_matrix));
	//final score
	mres.vd_trace.push_back(graph_dist(g,h,gm_P,cscore_matrix));
	//other output parameters
	mres.dres=mres.vd_trace.at(1);
	//transpose matrix save
	mres.gm_P=gm_P;
	return mres;
}

//convex-concave function value
double sgm_algorithm_spathnew::f_qcvqcc(gsl_matrix *gm_Ag_d,gsl_matrix *gm_Ah_d,gsl_matrix *gm_Lg_d,gsl_matrix *gm_Lh_d, gsl_matrix* gm_Lg_d_sub21, gsl_matrix* gm_Lh_d_sub21t, gsl_matrix* gm_Delta,gsl_matrix* gm_P,double dlambda,gsl_matrix * gm_temp,gsl_matrix *gm_temp2, gsl_matrix* gm_Ag_d_sub21, gsl_matrix* gm_Ag_d_sub21t, gsl_matrix* gm_Ah_d_sub21, gsl_matrix* gm_Ah_d_sub21t)
{
	double f1=f_qcv1(  gm_Ag_d,gm_Ah_d,gm_P,gm_temp,false);
	double f2=f_qcc(  gm_Lg_d,gm_Lh_d,gm_Delta,gm_P,gm_temp,gm_temp2);
	double dres=dlambda*f1+(1-dlambda)*f2;
	//QCV update
	gsl_matrix* gm_temp10 = gsl_matrix_alloc((N-m),m);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_P,gm_Ah_d_sub21,0,gm_temp10);
	gsl_matrix* gm_temp12 = gsl_matrix_alloc(m,(N-m));
	gsl_matrix_transpose_memcpy(gm_temp12, gm_temp10);
	gsl_matrix_mul_elements(gm_temp12, gm_Ag_d_sub21t);
	dres=dres-4*gsl_matrix_sum(gm_temp12)*dlambda/df_norm;

	gsl_matrix* gm_temp11 = gsl_matrix_alloc(m,(N-m));
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_Ag_d_sub21t,gm_P,0,gm_temp11);
	gsl_matrix_mul_elements(gm_temp11, gm_temp11);
	dres=dres+gsl_matrix_sum(gm_temp11)*dlambda/df_norm;

	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_P,gm_Ah_d_sub21,0,gm_temp10);
	gsl_matrix_mul_elements(gm_temp10, gm_temp10);
	dres=dres+gsl_matrix_sum(gm_temp10)*dlambda/df_norm;
	//QCC update
	gsl_matrix_transpose_memcpy(gm_temp,gm_Delta);
	gsl_matrix_mul_elements(gm_temp,gm_P);
	dres=dres-gsl_matrix_sum(gm_temp)/df_norm;
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_Lg_d_sub21,gm_Lh_d_sub21t,0,gm_temp2);
	gsl_matrix_transpose(gm_P);
	gsl_matrix_mul_elements(gm_temp2,gm_P);
	dres=dres-gsl_matrix_sum(gm_temp2)*(1-dlambda)/df_norm;
	gsl_matrix_transpose(gm_P);
	return dres;
}

//gradient for qcvqcc function
void sgm_algorithm_spathnew::qcvqcc_gradient(gsl_matrix *gm_Ag_d,gsl_matrix *gm_Ah_d,gsl_matrix *gm_Lg_d,gsl_matrix *gm_Lh_d,gsl_vector* gv_P, gsl_vector* gv_grad,double dlambda,gsl_vector * gv_temp)
{
	gsl_matrix_view gmv_temp=gsl_matrix_view_array(gv_temp->data,(N-m),(N-m));
	qcv_gradient_opt(gm_Ag_d,gm_Ah_d,gv_P,gv_grad,&gmv_temp.matrix);
	gsl_vector_scale(gv_grad,dlambda);
	qcc_gradient_opt(gm_Lg_d,gm_Lh_d,gv_P,gv_grad,-(1-dlambda),&gmv_temp.matrix);

}

//gradient for qcc function, here gm_A*_d are the laplacian matrices,update current gradient value
void sgm_algorithm_spathnew::qcc_gradient_opt(gsl_matrix *gm_Ag_d,gsl_matrix *gm_Ah_d,gsl_vector* gv_P, gsl_vector* gv_grad,double dmult,gsl_matrix * gm_temp)
{
	gsl_matrix_view gmv_grad=gsl_matrix_view_array(gv_grad->data,(N-m),(N-m));
	gsl_matrix_transpose(&gmv_grad.matrix);
	gsl_matrix_view gmv_P=gsl_matrix_view_array(gv_P->data,(N-m),(N-m));
	gsl_matrix_transpose(&gmv_P.matrix);
	//if (bnosymm) gsl_matrix_transpose(gm_Ah_d);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_Ah_d,&gmv_P.matrix,0,gm_temp);

	//if (bnosymm) gsl_matrix_transpose(gm_Ah_d);
	gsl_matrix_transpose(&gmv_P.matrix);
	gsl_matrix_transpose(gm_temp);

	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,2*dmult*(1-dalpha_ldh)/df_norm,gm_Ag_d,gm_temp,1,&gmv_grad.matrix);
	gsl_matrix_transpose(&gmv_grad.matrix);
}

//optimized version of gradient calculations: tensor product tricks
// Note: returns  the transpose of gradient
void sgm_algorithm_spathnew::qcv_gradient_opt(gsl_matrix *gm_Ag_d,gsl_matrix *gm_Ah_d,gsl_vector* gv_P,  gsl_vector* gv_grad, gsl_matrix * gm_temp)
{


	// This block computes 1st term of convex component  in 3.5.2  of TPAMI paper
	if (bnosymm) gsl_matrix_transpose(gm_Ah_d);
	gsl_matrix *gm_1=gsl_matrix_alloc(N,N);
	gsl_matrix_view gmv_grad=gsl_matrix_view_array(gv_grad->data,N,N);
	gsl_matrix_view gmv_P=gsl_matrix_view_array(gv_P->data,N,N);
	//nonseededPtoseededP(&gmv_P.matrix, m);
	// gm_temp= Ag*P
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_Ag_d,&gmv_P.matrix,0,gm_temp);
	if (bnosymm) gsl_matrix_transpose(gm_Ag_d);
	// \Grad <- Ag*g_temp = Ag*(Ag*P)
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_Ag_d,gm_temp,0,&gmv_grad.matrix);
	//I gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,gm_Ah_d,gm_Ah_d,0,gm_temp);
	//I gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,&gmv_P.matrix,gm_temp,0,gm_1);


	// This block computes 4th term of convex component  in 3.5.2  of TPAMI paper
	// gmv_P <- P'
	gsl_matrix_transpose(&gmv_P.matrix);
	// gm_temp <- Ah*P'
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_Ah_d,&gmv_P.matrix,0,gm_temp);
	if (bnosymm) gsl_matrix_transpose(gm_Ah_d);
	// gm_1 <- Ah * gm_temp = Ah*Ah*P'
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_Ah_d,gm_temp,0,gm_1);
	// gm_1 <- gm_1'
	// gm_1  = P*Ah*Ah
	gsl_matrix_transpose(gm_1);
	// \Grad = Ag*(Ag*P) +P*Ah*Ah
	gsl_matrix_add(&gmv_grad.matrix,gm_1);

	// third term
	if (bnosymm) gsl_matrix_transpose(gm_Ah_d);
	//gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,&gmv_P.matrix,gm_Ah_d,0,gm_temp);
	// gm_temp <- Ah*P'
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_Ah_d,&gmv_P.matrix,0,gm_temp);
	// gm_temp = P*Ah
	gsl_matrix_transpose(gm_temp);
	// \Grad  +=  -Ag*P*Ah
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,-1,gm_Ag_d,gm_temp,1,&gmv_grad.matrix);
	//second term
	if (bnosymm) gsl_matrix_transpose(gm_Ag_d);
	if (bnosymm) gsl_matrix_transpose(gm_Ah_d);
	// gm_temp = Ah*P'
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_Ah_d,&gmv_P.matrix,0,gm_temp);
	// gm_temp = P*Ah'
	gsl_matrix_transpose(gm_temp);

	// -Ag*P*Ah'
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,-1,gm_Ag_d,gm_temp,1,&gmv_grad.matrix);
  // gmv_P <-P
  //
	gsl_matrix_transpose(&gmv_P.matrix);

	if (bnosymm) gsl_matrix_transpose(gm_Ah_d);
	if (bnosymm) gsl_matrix_transpose(gm_Ah_d);  //Sancar: should this be Ag?
	//III gsl_blas_dgemm_sym(CblasLeft,CblasUpper,CblasNoTrans,CblasNoTrans,-2,gm_Ag_d,gm_temp,1,&gmv_grad.matrix);
	gsl_matrix_free(gm_1);
  //  gmv_grad <- \Grad_transpose
  //  gmv_grad <- \Grad'
	gsl_matrix_transpose(&gmv_grad.matrix);
	gsl_matrix_scale(&gmv_grad.matrix,1-dalpha_ldh);//label cost function scaling
	//remove df_norm
	gsl_vector_scale(gv_grad,1/df_norm);//normalization
	nonseededGradtoseededGrad(&gmv_grad.matrix,m);

}



// //concave function value

double sgm_algorithm_spathnew::f_qcc( gsl_matrix *gm_Ag_d,gsl_matrix *gm_Ah_d,gsl_matrix* gm_Delta,gsl_matrix* gm_P,gsl_matrix * gm_temp,gsl_matrix *gm_temp2)
{
	bool bnosymm = false;
	gsl_matrix_transpose(gm_P);
	if (bnosymm) gsl_matrix_transpose(gm_Ah_d);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_Ah_d,gm_P,0,gm_temp);
	if (bnosymm) gsl_matrix_transpose(gm_Ah_d);
	gsl_matrix_transpose(gm_P);
	gsl_matrix_transpose(gm_temp);
	if (bnosymm) gsl_matrix_transpose(gm_Ag_d);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_Ag_d,gm_temp,0,gm_temp2);
	if (bnosymm) gsl_matrix_transpose(gm_Ag_d);
	gsl_matrix_memcpy(gm_temp,gm_P);
	gsl_matrix_mul_elements(gm_temp,gm_temp2);
	double dres=-2*gsl_matrix_sum(gm_temp);
	gsl_matrix_transpose_memcpy(gm_temp,gm_Ah_d);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_temp,gm_Ah_d,0,gm_temp2);
	double dconst_add=0;
	for (int i=0;i<gm_temp2->size1;i++)
		dconst_add+=gsl_matrix_get(gm_temp2,i,i);
	gsl_matrix_transpose_memcpy(gm_temp,gm_Ag_d);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_temp,gm_Ag_d,0,gm_temp2);
	for (int i=0;i<gm_temp2->size1;i++)
		dconst_add+=gsl_matrix_get(gm_temp2,i,i);
	dres+=dconst_add;
	dres=dres*(1-dalpha_ldh);
	dres/=df_norm;
	return dres;
}

//convex-concave function value
double sgm_algorithm_spathnew::f_qcvqcc(gsl_matrix *gm_Ag_d,gsl_matrix *gm_Ah_d,gsl_matrix *gm_Lg_d,gsl_matrix *gm_Lh_d, gsl_matrix* gm_Lg_d_sub21, gsl_matrix* gm_Lh_d_sub21t, gsl_matrix* gm_Delta,gsl_matrix* gm_P,double dlambda,gsl_matrix * gm_temp,gsl_matrix *gm_temp2, gsl_matrix* gm_Ag_d_sub21, gsl_matrix* gm_Ag_d_sub21t, gsl_matrix* gm_Ah_d_sub21, gsl_matrix* gm_Ah_d_sub21t)
{
	double f1=f_qcv(gm_Ag_d,gm_Ah_d,gm_P,gm_temp,false);
	gsl_matrix *A_g_temp=gsl_matrix_alloc(N,N);
	gsl_matrix *A_h_temp=gsl_matrix_alloc(N,N);
	double f2=f_qcc(gm_Lg_d,gm_Lh_d,gm_Delta,gm_P,gm_temp,gm_temp2);//,A_g_temp,A_h_temp);
	double dres=dlambda*f1+(1-dlambda)*f2;
	gsl_matrix_transpose_memcpy(gm_temp,gm_Delta);
	gsl_matrix_mul_elements(gm_temp,gm_P);
	dres=dres-gsl_matrix_sum(gm_temp);
	gsl_matrix_free(A_g_temp);
	gsl_matrix_free(A_h_temp);
	return dres;
}



// //convex representation of the objective function ||AP-PA||^2_F
// double sgm_algorithm_spathnew::f_qcv(gsl_matrix *gm_Ag_d,gsl_matrix *gm_Ah_d,gsl_matrix* gm_P,gsl_matrix * gm_temp,bool bqcv)
// {

// 	//nonseededPtoseededP(gm_P,m);
// 	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_Ag_d,gm_P,0,gm_temp);
// 	//gsl_blas_dgemm_sym(CblasLeft,CblasUpper,CblasNoTrans,CblasNoTrans,1,gm_Ag_d,gm_P,0,gm_temp);
// 	gsl_matrix_transpose(gm_P);
// 	gsl_matrix_transpose(gm_temp);
//         if (bnosymm) gsl_matrix_transpose(gm_Ah_d);
// 	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,-1,gm_Ah_d,gm_P,1,gm_temp);
// 	if (bnosymm) gsl_matrix_transpose(gm_Ah_d);
// 	gsl_matrix_transpose(gm_P);
// 	//gsl_matrix_view gm_temp_upper_left = gsl_matrix_submatrix(gm_temp,0,0,m,m);
// 	//gsl_matrix_set_zero(&gm_temp_upper_left.matrix);
// 	//gsl_blas_dgemm_sym(CblasLeft,CblasUpper,CblasTrans,CblasTrans,-1,gm_Ah_d,gm_P,1,gm_temp);
// 	//double dres= pow(gsl_matrix_norm(gm_temp_upper_left,2),2);
// 	double dres= pow(gsl_matrix_norm(gm_temp,2),2);
// 	dres=dres*(1-dalpha_ldh);
// 	dres/=df_norm;
// 	if ((bqcv) and (dalpha_ldh>0))//if just QCV optimization then alpha_ldh must be integrated, in QCVQCC and P cases we use Delta for these purposes
// 	{gsl_vector_view gvv_ldh=gsl_vector_view_array(gm_ldh->data,N*N);
// 	gsl_vector_view gvv_P=gsl_vector_view_array(gm_P->data,N*N);
// 	double dtr;
// 	gsl_blas_ddot(&gvv_ldh.vector,&gvv_P.vector,&dtr);
// 	dtr*=dalpha_ldh;
// 	dres-=dtr;};
// 	//normalization
// 	return dres;
// }


//convex representation of the objective function ||AP-PA||^2_F
double sgm_algorithm_spathnew::f_qcv(gsl_matrix *gm_Ag_d,gsl_matrix *gm_Ah_d,gsl_matrix* gm_P,gsl_matrix * gm_temp,bool bqcv)
{
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_Ag_d,gm_P,0,gm_temp);
	//gsl_blas_dgemm_sym(CblasLeft,CblasUpper,CblasNoTrans,CblasNoTrans,1,gm_Ag_d,gm_P,0,gm_temp);
	gsl_matrix_transpose(gm_P);
	gsl_matrix_transpose(gm_temp);
        if (bnosymm) gsl_matrix_transpose(gm_Ah_d);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,-1,gm_Ah_d,gm_P,1,gm_temp);
	if (bnosymm) gsl_matrix_transpose(gm_Ah_d);
	gsl_matrix_transpose(gm_P);
	//gsl_blas_dgemm_sym(CblasLeft,CblasUpper,CblasTrans,CblasTrans,-1,gm_Ah_d,gm_P,1,gm_temp);
	double dres= pow(gsl_matrix_norm(gm_temp,2),2);
	dres=dres*(1-dalpha_ldh);
	//remove df_norm
	dres/=df_norm;
	if ((bqcv) and (dalpha_ldh>0))//if just QCV optimization then alpha_ldh must be integrated, in QCVQCC and P cases we use Delta for these purposes
	{gsl_vector_view gvv_ldh=gsl_vector_view_array(gm_ldh->data,N*N);
	gsl_vector_view gvv_P=gsl_vector_view_array(gm_P->data,N*N);
	double dtr;
	gsl_blas_ddot(&gvv_ldh.vector,&gvv_P.vector,&dtr);
	dtr*=dalpha_ldh;
	dres-=dtr;};
	//normalization
	return dres;
}



//concave function value
double sgm_algorithm_spathnew::f_qcc(gsl_matrix *gm_Ag_d,gsl_matrix *gm_Ah_d,gsl_matrix* gm_Delta,gsl_matrix* gm_P,gsl_matrix * gm_temp,gsl_matrix *gm_temp2)
{
	gsl_matrix_transpose(gm_P);
	if (bnosymm) gsl_matrix_transpose(gm_Ah_d);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_Ah_d,gm_P,0,gm_temp);
	if (bnosymm) gsl_matrix_transpose(gm_Ah_d);
	gsl_matrix_transpose(gm_P);
	gsl_matrix_transpose(gm_temp);
	if (bnosymm) gsl_matrix_transpose(gm_Ag_d);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_Ag_d,gm_temp,0,gm_temp2);
	if (bnosymm) gsl_matrix_transpose(gm_Ag_d);
	gsl_matrix_memcpy(gm_temp,gm_P);
	gsl_matrix_mul_elements(gm_temp,gm_temp2);
	double dres=-2*gsl_matrix_sum(gm_temp);
	gsl_matrix_transpose_memcpy(gm_temp,gm_Ah_d);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_temp,gm_Ah_d,0,gm_temp2);
	double dconst_add=0;
	for (int i=0;i<gm_temp2->size1;i++)
		dconst_add+=gsl_matrix_get(gm_temp2,i,i);
	gsl_matrix_transpose_memcpy(gm_temp,gm_Ag_d);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_temp,gm_Ag_d,0,gm_temp2);
	for (int i=0;i<gm_temp2->size1;i++)
		dconst_add+=gsl_matrix_get(gm_temp2,i,i);
	dres+=dconst_add;
	dres=dres*(1-dalpha_ldh);
	//remove df_norm
	dres/=df_norm;
	return dres;
}

void sgm_algorithm_spathnew::nonseededGradtoseededGrad(gsl_matrix *Grad, int m_seeds){
	gsl_matrix_view sub_grad = gsl_matrix_submatrix(Grad,0,0,m_seeds,m_seeds);
	gsl_matrix_set_zero(&sub_grad.matrix);
	sub_grad = gsl_matrix_submatrix(Grad,0,m_seeds,m_seeds,N-m_seeds);
	gsl_matrix_set_zero(&sub_grad.matrix);
	sub_grad = gsl_matrix_submatrix(Grad,m_seeds,0,N-m_seeds,m_seeds);
	gsl_matrix_set_zero(&sub_grad.matrix);
}
