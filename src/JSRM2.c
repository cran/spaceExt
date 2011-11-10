#include <stdlib.h>
void JSRM2( int* NN, int * PP, float * L1, float * L2, float *Y_data, float *sigma_sr, int * N_iter, float *haction, float *beta_new)
{
 /////////////////////////  In put data Y.m /////////////////////////////
 

 int n, p, n_iter;
 float lambda1, lambda2;
 float *Y_m;

 float *meanx;
 float *normx;
 float *B_s;
 float *h_a;
 float *B;
 float *Yhat_m;
 float *E_m;
 int *pick;


 float temp, temp1;
 float  Aij, Aji;
 float  beta_next;
 float beta_change;
 float maxdif;
 float tempChange;
 int i,j,k;
 int iter, nbeta;
 int nrow_pick, jrow;
 int change_i, change_j;
 int cur_i,cur_j;

 float eps1;       //   tolerant value;


 iter=0;
 eps1 = 1e-6;
 maxdif = -100.0;


 n=*NN;
 p=*PP;
 n_iter=*N_iter;
 nbeta = p*(p-1)/2;

 lambda1=*L1;
 lambda2=*L2;

 meanx=(float *) malloc(p*sizeof(float));
 normx=(float *) malloc(p*sizeof(float));  //kesi
 Y_m=(float *) malloc(n*p*sizeof(float));

 for(i=0; i<n; i++)
   for(j=0; j<p; j++)
     Y_m[i*p+j]=Y_data[i*p+j];


 ////// normalize each column of Y_m into mean=0 and norm=1

    for(j=0;j<p;j++)
        meanx[j] = 0 ;
    for(j=0;j<p;j++)
        for(i=0;i<n;i++)
           meanx[j] = meanx[j]+ Y_m[i*p+j] ;
    for(j=0;j<p;j++)
        meanx[j] = meanx[j]/n;
    for(j=0;j<p;j++)
        for(i=0;i<n;i++)
            Y_m[i*p+j] = Y_m[i*p+j] - meanx[j];



		for(k=0;k<p;k++)
		{
			normx[k]=0;  
			for(i=0;i<n;i++)
				normx[k] = normx[k]+ Y_m[i*p+k]*Y_m[i*p+k];
		}
 //////////////////////////////     Step 0      //////////////////////////////
 ////////////////////////////// Get initial value (Equation 3)



 B_s=(float *) malloc(p*p*sizeof(float));
 h_a=(float *) malloc(p*p*sizeof(float));
 B=(float *) malloc(p*p*sizeof(float));


  for(i = 0;i<p;i++)
    {
      for(j=0;j<p;j++)  //ratio of sig_i/sig_j
        B[i*p+j] = (sigma_sr[i]/sigma_sr[j]);
   }

 for(i = 0;i<p-1;i++)  //denominator s-11 with haction
     for(j=i;j<p;j++)
      {
      B_s[i*p+j]=B[i*p+j]*B[i*p+j]*normx[i] +  B[j*p+i]*B[j*p+i]*normx[j];
      B_s[j*p+i]=B_s[i*p+j];
	  h_a[i*p+j]=haction[j*p+j]/(sigma_sr[i]*sigma_sr[i])+haction[i*p+i]/(sigma_sr[j]*sigma_sr[j]);
	  h_a[j*p+i]=h_a[i*p+j];
      }


 ///////////////////////// End of Step 0 ///////////////////////////////////

 E_m=(float *) malloc(n*p*sizeof(float));
 pick=(int *) malloc(nbeta*2*sizeof(int));

/////////////////////// Step 1:   Get Initial E
/////////////////////// Equation (4)-(6)




 Yhat_m=(float *) malloc(n*p*sizeof(float));

 for(k=0; k<n; k++)
   for(j=0; j<p; j++)
     {
       Yhat_m[k*p+j]=0;
       for(i=0; i<p; i++)
         Yhat_m[k*p+j]=Yhat_m[k*p+j]+Y_m[k*p+i]*beta_new[i*p+j]*B[i*p+j];
       E_m[k*p+j]=Y_m[k*p+j]-Yhat_m[k*p+j];
     }

 free(Yhat_m);

///////////////////////// Step 2: update one beta  //////////////////////////////



  for(i = 0;i<nbeta;i++)
    {
     for(j=0;j<2;j++)   // pick[q] the q(th) beta_ij i=pick[q*2+0] j=pick[q*2+1]
        pick[i*2+j] = 0 ;
    }

  ///////// Get active Set;
     k = 0;
     for(j = p-1; j>=1; j--)
      {
       for(i = j-1; i>=0; i--)
       {
         if( beta_new[i*p+j]>  eps1  || beta_new[i*p+j] < -eps1 ){
			 pick[k*2+0] =i; pick[k*2+1] = j;
			 k = k + 1;
			 break;
       	 }
       }
       if (k>0)
          break;
      }

 

	  
 if(k>0) //otherwise, converge to 0 at the initial step.
 {


  ///////// change one beta
     cur_i = pick[0];
     cur_j = pick[1];

  //////// Equation (s-16)
     Aij=0;
     Aji=0;
	 float numer_ij,numer_ji;
	 numer_ij=-haction[cur_i*p+cur_j];
	 numer_ji=numer_ij;
     for(k=0; k<n; k++){
		 Aij=Aij+E_m[k*p+cur_j]*Y_m[k*p+cur_i];
		 Aji=Aji+E_m[k*p+cur_i]*Y_m[k*p+cur_j];
     }
	 for(k=0;k<p;k++){
		 if(k!=cur_i && k!=cur_j){
			 numer_ij=numer_ij+beta_new[k*p+cur_i]*haction[k*p+cur_j];
			 numer_ji=numer_ji+beta_new[k*p+cur_j]*haction[k*p+cur_i];
		}
	 }
     Aij=Aij*B[cur_i*p+cur_j]-numer_ij/(sigma_sr[cur_i]*sigma_sr[cur_i]);
     Aji=Aji*B[cur_j*p+cur_i]-numer_ji/(sigma_sr[cur_j]*sigma_sr[cur_j]);

  //////// Equation (10)
     beta_next=Aij+Aji+beta_new[cur_i*p+cur_j]*B_s[cur_i*p+cur_j];  //beta_new not updated yet

      ///shrink beta_next
      temp1 = beta_next;
      if ( beta_next > 0 )
            temp = beta_next - lambda1;
      else
            temp = - beta_next - lambda1;
      if(temp < 0 )
            temp = 0;
      else
        {
            temp =  temp /(B_s[cur_i*p+cur_j]+h_a[cur_i*p+cur_j]+lambda2);   // question???
			if(temp1 < 0 )
				temp = (-1) * temp;
        }

		  beta_change=beta_new[cur_i*p+cur_j]-temp;
		  tempChange=beta_change;
		  if(tempChange<0) tempChange=-tempChange;
		  if(tempChange>maxdif) maxdif=tempChange;
          beta_new[cur_i*p+cur_j] = temp;
          beta_new[cur_j*p+cur_i] = temp;
          change_i=cur_i;
          change_j=cur_j;

 //////////////////////////////////////////////////////////////////////////////////
 ///////////////////////// Step 3: begin to iterate  //////////////////////////////

 
 for ( iter = 0; iter < n_iter; iter++ ) { // iteration for all loop
 

	maxdif=-10;
    for(i = 0;i<nbeta;i++)
    {
      for(j=0;j<2;j++)
        pick[i*2+j] = 0 ;
    }

    
	// Get active Set;
	
   k = 0;
   for(j = p-1; j>=1; j--)
     for(i = j-1; i>=0; i--)
     {
       if( beta_new[i*p+j]>  eps1  || beta_new[i*p+j] < -eps1 ){
		   pick[k*2] =i; pick[k*2+1] = j; k = k + 1; 
	   }
     }


  nrow_pick = k;
   	
 if(nrow_pick>0)  // otherwise, go to all loop directly.
  {



   /////////////////////// loop for active set ///////////////////
   for(jrow = 0; jrow<=nrow_pick - 1; jrow++)
   {

      cur_i = pick[jrow*2];
      cur_j = pick[jrow*2+1];

     /////////// Update Residue   /////////////////
     /////////// Equation (11)

    for(k=0; k<n; k++)
    {
       E_m[k*p+change_i]=E_m[k*p+change_i]+ Y_m[k*p+change_j]*beta_change*B[change_j*p+change_i];
       E_m[k*p+change_j]=E_m[k*p+change_j]+ Y_m[k*p+change_i]*beta_change*B[change_i*p+change_j];
    }

    ///////////// update beta
    //////////// Equation (12)
     Aij=0;
     Aji=0;
	 numer_ij=-haction[cur_i*p+cur_j];
	 numer_ji=numer_ij;
     for(k=0; k<n; k++){
       Aij=Aij+E_m[k*p+cur_j]*Y_m[k*p+cur_i];
       Aji=Aji+E_m[k*p+cur_i]*Y_m[k*p+cur_j];
     }
	 for(k=0;k<p;k++){
		 if(k!=cur_i && k!=cur_j){
			 numer_ij=numer_ij+beta_new[k*p+cur_i]*haction[k*p+cur_j];
			 numer_ji=numer_ji+beta_new[k*p+cur_j]*haction[k*p+cur_i];
		}
	 }
     Aij=Aij*B[cur_i*p+cur_j]-numer_ij/(sigma_sr[cur_i]*sigma_sr[cur_i]);
     Aji=Aji*B[cur_j*p+cur_i]-numer_ji/(sigma_sr[cur_j]*sigma_sr[cur_j]);

     beta_next=Aij+Aji+beta_new[cur_i*p+cur_j]*B_s[cur_i*p+cur_j];

    ///shrink beta_next
      temp1 = beta_next;
      if ( beta_next > 0 )
            temp = beta_next - lambda1;
      else
            temp = - beta_next - lambda1;
      if(temp < 0 )
            temp = 0;
      else
        {
            temp =  temp /(B_s[cur_i*p+cur_j]+h_a[cur_i*p+cur_j]+lambda2);   // question???
         if(temp1 < 0 )
            temp = (-1) * temp;
        }

			  
		  beta_change=beta_new[cur_i*p+cur_j]-temp;
		  tempChange=beta_change;
		  if(tempChange<0) tempChange=-tempChange;
		  if(tempChange>maxdif) maxdif=tempChange;
          beta_new[cur_i*p+cur_j] = temp;
          beta_new[cur_j*p+cur_i] = temp;
          change_i=cur_i;
          change_j=cur_j;


   }///////////////////////  End of loop for active set ///////////////////
  
   }// end of if(nrow_pick>0)

   // if convergent on active set, move to full set
//////////////$$$$$$$$$$$  stop full set:: 
   if ( maxdif < 1e-6 || nrow_pick<1)
     {
	     maxdif=-10;
       /////////////////////// loop for all set //////////////////////
	   


	for( cur_i = 0; cur_i < p-1 ; cur_i++ )
	    for ( cur_j = cur_i + 1; cur_j < p ; cur_j++ )
	    {


	   if(beta_change< -eps1 || beta_change> eps1)
	   {
	  
	     ///////////    Update E_m   /////////////////
	     //////////// Equation (11)
	       for(k=0; k<n; k++)
	         {
	           E_m[k*p+change_i]=E_m[k*p+change_i]+ Y_m[k*p+change_j]*beta_change*B[change_j*p+change_i];
	           E_m[k*p+change_j]=E_m[k*p+change_j]+ Y_m[k*p+change_i]*beta_change*B[change_i*p+change_j];
             }

	   }

	 ///////////// update beta
    
	Aij=0;
     Aji=0;
	 numer_ij=-haction[cur_i*p+cur_j];
	 numer_ji=numer_ij;
    for(k=0; k<n; k++){
		Aij=Aij+E_m[k*p+cur_j]*Y_m[k*p+cur_i];
		Aji=Aji+E_m[k*p+cur_i]*Y_m[k*p+cur_j];
    }
	for(k=0;k<p;k++){
		 if(k!=cur_i && k!=cur_j){
			 numer_ij=numer_ij+beta_new[k*p+cur_i]*haction[k*p+cur_j];
			 numer_ji=numer_ji+beta_new[k*p+cur_j]*haction[k*p+cur_i];
		}
	 }
     Aij=Aij*B[cur_i*p+cur_j]-numer_ij/(sigma_sr[cur_i]*sigma_sr[cur_i]);
     Aji=Aji*B[cur_j*p+cur_i]-numer_ji/(sigma_sr[cur_j]*sigma_sr[cur_j]);
     
  //////// Equation (10)
     beta_next=Aij+Aji+beta_new[cur_i*p+cur_j]*B_s[cur_i*p+cur_j];

      ///shrink beta_next

      if ( beta_next > 0 )
            temp = beta_next - lambda1;
      else
            temp = - beta_next - lambda1;
      if(temp < 0 )
            temp = 0;
      else
        {
            temp =  temp /(B_s[cur_i*p+cur_j]+h_a[cur_i*p+cur_j]+lambda2);   // question???
         if(beta_next < 0 )
            temp = (-1) * temp;
        }
		  beta_change=beta_new[cur_i*p+cur_j]-temp;
		  tempChange=beta_change;
		  if(tempChange<0) tempChange=-tempChange;
		  if(tempChange>maxdif) maxdif=tempChange;
          beta_new[cur_i*p+cur_j] = temp;
          beta_new[cur_j*p+cur_i] = temp;
          change_i=cur_i;
          change_j=cur_j;

	 }///////////////////////  End of loop for all set ///////////////////

      
 

     }   // end of if(maxdif<1e-6)

	if( maxdif <1e-06)
          break ;
  }// end of iter
 
 
 
 
}// end of if(k>0)
 /////////////////////   End of Step 1,2,...........///////////////////////////






 free(E_m);
 free(B_s);
 free(B);
 free(Y_m);
 free(pick);
 free(meanx);
 free(normx);
 }  // End of jsrm.shoot();

