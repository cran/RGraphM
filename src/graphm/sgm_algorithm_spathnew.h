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
#ifndef SGM_ALGORITHM_sPATHNEW_H
#define SGM_ALGORITHM_sPATHNEW_H
#include "sgm_algorithm.h"
#include "algorithm_qcv.h"

/**
A path following graph matching algorithm

	@author Sancar Adali <sadali@gmail.com>;Mikhail Zaslavskiy <mikhail.zaslavskiy@ensmp.fr>
*/
class sgm_algorithm_spathnew : public sgm_algorithm
{

public:
	//virtual match_result match(graph &g,graph &h,int m, gsl_matrix* gm_P_i=NULL, gsl_matrix* gm_ldh=NULL,double dalpha_ldh=-1);
	//virtual void qcc_gradient_opt(int m, gsl_matrix *gm_Ag_d,gsl_matrix *gm_Ah_d,gsl_vector* gv_P, gsl_vector* gv_grad,double dmult,gsl_matrix * gm_temp);
	//void qcvqcc_gradient(int m, gsl_matrix *gm_Ag_d,gsl_matrix *gm_Ah_d,gsl_matrix *gm_Lg_d,gsl_matrix *gm_Lh_d,gsl_vector* gv_P, gsl_vector* gv_grad,double dlambda,gsl_vector * gv_temp);
 //double f_qcvqcc(int m, gsl_matrix *gm_Ag_d,gsl_matrix *gm_Ah_d,gsl_matrix *gm_Lg_d,gsl_matrix *gm_Lh_d, gsl_matrix* gm_Lg_d_sub21, gsl_matrix* gm_Lh_d_sub21t, gsl_matrix* gm_Delta,gsl_matrix* gm_P,double dlambda,gsl_matrix * gm_temp,gsl_matrix *gm_temp2, gsl_matrix* gm_Ag_d_sub21, gsl_matrix* gm_Ag_d_sub21t, gsl_matrix* gm_Ah_d_sub21, gsl_matrix* gm_Ah_d_sub21t);

	// double f_qcvqcc(gsl_matrix *gm_Ag_d,gsl_matrix *gm_Ah_d,gsl_matrix *gm_Lg_d,gsl_matrix *gm_Lh_d,gsl_matrix* gm_Delta,gsl_matrix* gm_P,double dlambda,gsl_matrix * gm_temp,gsl_matrix *gm_temp2);
	//concave function value
 //virtual double f_qcc(gsl_matrix *gm_Ag_d,gsl_matrix *gm_Ah_d,gsl_matrix* gm_Delta,gsl_matrix* gm_P,gsl_matrix * gm_temp,gsl_matrix *gm_temp2);
	// gsl_matrix*  submatrix(gsl_matrix* mainMatrix, int i, int j, int h, int w);




	  sgm_algorithm_spathnew(int);
    virtual match_result match(graph &g,graph &h,gsl_matrix* gm_P_i=NULL, gsl_matrix* gm_ldh=NULL,double dalpha_ldh=-1);

		match_result match_with_seeds(graph& g, graph& h, gsl_matrix* gm_P_i=NULL, gsl_matrix* gm_ldh=NULL,double dalpha_ldh=-1, unsigned int m_seeds=0);

    void qcv_gradient_opt(gsl_matrix *gm_Ag_d,gsl_matrix *gm_Ah_d,gsl_vector* gv_P, gsl_vector* gv_grad,gsl_matrix * gm_temp);
    void qcc_gradient_opt(gsl_matrix *gm_Ag_d,gsl_matrix *gm_Ah_d,gsl_vector* gv_P, gsl_vector* gv_grad,double dmult,gsl_matrix * gm_temp);
    void qcvqcc_gradient(gsl_matrix *gm_Ag_d,gsl_matrix *gm_Ah_d,gsl_matrix *gm_Lg_d,gsl_matrix *gm_Lh_d,gsl_vector* gv_P, gsl_vector* gv_grad,double dlambda,gsl_vector * gv_temp);
    double f_qcvqcc(gsl_matrix *gm_Ag_d,gsl_matrix *gm_Ah_d,gsl_matrix *gm_Lg_d,gsl_matrix *gm_Lh_d, gsl_matrix* gm_Lg_d_sub21, gsl_matrix* gm_Lh_d_sub21t, gsl_matrix* gm_Delta,gsl_matrix* gm_P,double dlambda,gsl_matrix * gm_temp,gsl_matrix *gm_temp2, gsl_matrix* gm_Ag_d_sub21, gsl_matrix* gm_Ag_d_sub21t, gsl_matrix* gm_Ah_d_sub21, gsl_matrix* gm_Ah_d_sub21t);

    //concave function value
//    double f_qcc(gsl_matrix *gm_Ag_d,gsl_matrix *gm_Ah_d,gsl_matrix* gm_Delta,gsl_matrix* gm_P,gsl_matrix * gm_temp,gsl_matrix *gm_temp2, gsl_matrix *gm_tempAg, gsl_matrix *gm_tempAh);
    double f_qcv1(gsl_matrix *gm_Ag_d,gsl_matrix *gm_Ah_d,gsl_matrix* gm_P,gsl_matrix * gm_temp,bool bqcv=false);
    double f_qcc(gsl_matrix *gm_Ag_d,gsl_matrix *gm_Ah_d,gsl_matrix* gm_Delta,gsl_matrix* gm_P,gsl_matrix * gm_temp,gsl_matrix *gm_temp2);
     void nonseededGradtoseededGrad(gsl_matrix *Grad, int m_seeds);
};

#endif
