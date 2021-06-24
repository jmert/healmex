

#define NMT_MAX(a,b)  (((a)>(b)) ? (a) : (b)) // maximum
#define NMT_MIN(a,b)  (((a)<(b)) ? (a) : (b)) // minimum

using namespace std;

//Returns all non-zero wigner-3j symbols
// il2 (in) : l2
// il3 (in) : l3
// im2 (in) : m2
// im3 (in) : m3
// l1min_out (out) : min value for l1
// l1max_out (out) : max value for l1
// thrcof (out) : array with the values of the wigner-3j
// size (in) : size allocated for thrcof
int drc3jj(int il2,int il3,int im2, int im3,int *l1min_out, int *l1max_out,std::vector<double> &thrcof,int size)
{
  int sign1,sign2,nfin,im1,l1max,l1min,ii,lstep;
  int converging,nstep2,nfinp1,index,nlim;
  double newfac,c1,c2,sum1,sum2,a1,a2,a1s,a2s,dv,denom,c1old,oldfac,l1,l2,l3,m1,m2,m3;
  double x,x1,x2,x3,y,y1,y2,y3,sumfor,sumbac,sumuni,cnorm,thresh,ratio;
  double huge=sqrt(1.79E308/20.0);
  double srhuge=sqrt(huge);
  double tiny=1./huge;
  double srtiny=1./srhuge;

  im1=-im2-im3;
  l2=(double)il2; l3=(double)il3;
  m1=(double)im1; m2=(double)im2; m3=(double)im3;
  
  if((abs(il2+im2-il3+im3))%2==0)
    sign2=1;
  else
    sign2=-1;
  
  //l1 bounds
  l1max=il2+il3;
  l1min=NMT_MAX((abs(il2-il3)),(abs(im1)));
  *l1max_out=l1max;
  *l1min_out=l1min;

  if((il2-abs(im2)<0)||(il3-abs(im3)<0)) {
    for(ii=0;ii<=l1max-l1min;ii++)
      thrcof[ii]=0;
    return 0;
  }
  
  if(l1max-l1min<0) //Check for meaningful values
    throw("WTF?\n");
  
  if(l1max==l1min) { //If it's only one value:
    thrcof[0]=sign2/sqrt(l1min+l2+l3+1);
    return 0;
  }
  else {
    nfin=l1max-l1min+1;
    if(nfin>size) //Check there's enough space
      throw("Output array is too small %d\n",nfin);
    else {
      l1=l1min;
      newfac=0.;
      c1=0.;
      sum1=(l1+l1+1)*tiny;
      thrcof[0]=srtiny;
      
      lstep=0;
      converging=1;
      while((lstep<nfin-1)&&(converging)) { //Forward series
	lstep++;
	l1++; //order
	
	oldfac=newfac;
	a1=(l1+l2+l3+1)*(l1-l2+l3)*(l1+l2-l3)*(-l1+l2+l3+1);
	a2=(l1+m1)*(l1-m1);
	newfac=sqrt(a1*a2);
	
	if(l1>1) {
	  dv=-l2*(l2+1)*m1+l3*(l3+1)*m1+l1*(l1-1)*(m3-m2);
	  denom=(l1-1)*newfac;
	  if(lstep>1)
	    c1old=fabs(c1);
	  c1=-(l1+l1-1)*dv/denom;
	}
	else {
	  c1=-(l1+l1-1)*l1*(m3-m2)/newfac;
	}
	
	if(lstep<=1) {
	  x=srtiny*c1;
	  thrcof[1]=x;
	  sum1+=tiny*(l1+l1+1)*c1*c1;
	}
	else {
	  c2=-l1*oldfac/denom;
	  x=c1*thrcof[lstep-1]+c2*thrcof[lstep-2];
	  thrcof[lstep]=x;
	  sumfor=sum1;
	  sum1+=(l1+l1+1)*x*x;
	  if(lstep<nfin-1) {
	    if(fabs(x)>=srhuge) {
	      for(ii=0;ii<=lstep;ii++) {
		if(fabs(thrcof[ii])<srtiny)
		  thrcof[ii]=0;
		thrcof[ii]/=srhuge;
	      }
	      sum1/=huge;
	      sumfor/=huge;
	      x/=srhuge;
	    }
	    
	    if(c1old<=fabs(c1))
	      converging=0;
	  }
	}
      }
      
      if(nfin>2) {
	x1=x;
	x2=thrcof[lstep-1];
	x3=thrcof[lstep-2];
	nstep2=nfin-lstep-1+3;
	
	nfinp1=nfin+1;
	l1=l1max;
	thrcof[nfin-1]=srtiny;
	sum2=tiny*(l1+l1+1);
	
	l1+=2;
	lstep=0;
	while(lstep<nstep2-1) { //Backward series
	  lstep++;
	  l1--;
	  
	  oldfac=newfac;
	  a1s=(l1+l2+l3)*(l1-l2+l3-1)*(l1+l2-l3-1)*(-l1+l2+l3+2);
	  a2s=(l1+m1-1)*(l1-m1-1);
	  newfac=sqrt(a1s*a2s);
	  
	  dv=-l2*(l2+1)*m1+l3*(l3+1)*m1+l1*(l1-1)*(m3-m2);
	  denom=l1*newfac;
	  c1=-(l1+l1-1)*dv/denom;
	  if(lstep<=1) {
	    y=srtiny*c1;
	    thrcof[nfin-2]=y;
	    sumbac=sum2;
	    sum2+=tiny*(l1+l1-3)*c1*c1;
	  }
	  else {
	    c2=-(l1-1)*oldfac/denom;
	    y=c1*thrcof[nfin-lstep]+c2*thrcof[nfinp1-lstep]; //is the index ok??
	    if(lstep!=nstep2-1) {
	      thrcof[nfin-lstep-1]=y; //is the index ok??
	      sumbac=sum2;
	      sum2+=(l1+l1-3)*y*y;
	      if(fabs(y)>=srhuge) {
		for(ii=0;ii<=lstep;ii++) {
		  index=nfin-ii-1; //is the index ok??
		  if(fabs(thrcof[index])<srtiny)
		    thrcof[index]=0;
		  thrcof[index]=thrcof[index]/srhuge;
		}
		sum2/=huge;
		sumbac/=huge;
	      }
	    }
	  }
	}
	
	y3=y;
	y2=thrcof[nfin-lstep]; //is the index ok??
	y1=thrcof[nfinp1-lstep]; //is the index ok??
	
	ratio=(x1*y1+x2*y2+x3*y3)/(x1*x1+x2*x2+x3*x3);
	nlim=nfin-nstep2+1;
	
	if(fabs(ratio)<1) {
	  nlim++;
	  ratio=1./ratio;
	  for(ii=nlim-1;ii<nfin;ii++) //is the index ok??
	    thrcof[ii]*=ratio;
	  sumuni=ratio*ratio*sumbac+sumfor;
	}
	else {
	  for(ii=0;ii<nlim;ii++)
	    thrcof[ii]*=ratio;
	  sumuni=ratio*ratio*sumfor+sumbac;
	}
      }
      else
	sumuni=sum1;
      
      cnorm=1./sqrt(sumuni);
      if(thrcof[nfin-1]<0) sign1=-1;
      else sign1=1;
      
      if(sign1*sign2<=0)
	cnorm=-cnorm;
      if(fabs(cnorm)>=1) {
	for(ii=0;ii<nfin;ii++)
	  thrcof[ii]*=cnorm;
	return 0;
      }
      else {
	thresh=tiny/fabs(cnorm);
	for(ii=0;ii<nfin;ii++) {
	  if(fabs(thrcof[ii])<thresh)
	    thrcof[ii]=0;
	  thrcof[ii]*=cnorm;
	}
	return 0;
      }
    } //Size is good
  } //Doing for many l1s
  
  return 2;
}

static int lend_toeplitz(int l2, int l_toeplitz, int l_exact, int dl_band, int lmax)
{
  int l_end;

  if(l_toeplitz > 0) {
    if(l2<=l_exact)
      l_end = lmax;
    else if(l2<=l_toeplitz)
      l_end = l2+dl_band;
    else
      l_end = l2;
  }
  else
    l_end = lmax;

  return fmin(l_end, lmax);
}

/*
static void populate_toeplitz(mcm_calculator *c, int lmax, int lmax_mask, int npcl, int s1, int s2, std::vector<double> pcl_masks, int lt)
{
  int ic,l2;
  double ***tplz_00, ***tplz_0s, ***tplz_pp, ***tplz_mm;
  if(c->has_00) {
    tplz_00=my_malloc(npcl*sizeof(double **));
    for(ic=0;ic<npcl;ic++) {
      tplz_00[ic]=my_malloc(2*sizeof(double *));
      tplz_00[ic][0]=my_calloc((lmax+1),sizeof(double));
      tplz_00[ic][1]=my_calloc((lmax+1),sizeof(double));
    }
  }
  if(c->has_0s) {
    tplz_0s=my_malloc(npcl*sizeof(double **));
    for(ic=0;ic<npcl;ic++) {
      tplz_0s[ic]=my_malloc(2*sizeof(double *));
      tplz_0s[ic][0]=my_calloc((lmax+1),sizeof(double));
      tplz_0s[ic][1]=my_calloc((lmax+1),sizeof(double));
    }
  }
  if(c->has_ss) {
    tplz_pp=my_malloc(npcl*sizeof(double **));
    tplz_mm=my_malloc(npcl*sizeof(double **));
    for(ic=0;ic<npcl;ic++) {
      tplz_pp[ic]=my_malloc(2*sizeof(double *));
      tplz_pp[ic][0]=my_calloc((lmax+1),sizeof(double));
      tplz_pp[ic][1]=my_calloc((lmax+1),sizeof(double));
      tplz_mm[ic]=my_malloc(2*sizeof(double *));
      tplz_mm[ic][0]=my_calloc((lmax+1),sizeof(double));
      tplz_mm[ic][1]=my_calloc((lmax+1),sizeof(double));
    }
  }

  int lstart=0;
  int max_spin=NMT_MAX(s1, s2);
  int c->has_ss2=(s1!=0) && (s2!=0) && (s1!=s2);
  if(!(c->has_00))
    lstart=max_spin;

  #pragma omp parallel default(none) shared(c, s1, s2, lmax, lmax_mask, npcl, lstart, pcl_masks, lt) shared(tplz_00,tplz_0s,tplz_pp,tplz_mm)
  {
    int il3,ll2,icc;
    int l3_list[2];
    double *wigner_00=NULL,*wigner_ss1=NULL,*wigner_ss2=NULL;
    if(c->has_00 || c->has_0s)
      wigner_00=my_malloc(2*(lmax_mask+1)*sizeof(double));
    if(c->has_0s || c->has_ss)
      wigner_ss1=my_malloc(2*(lmax_mask+1)*sizeof(double));
    if(c->has_ss2)
      wigner_ss2=my_malloc(2*(lmax_mask+1)*sizeof(double));
    else
      wigner_ss2=wigner_ss1;

    l3_list[0]=0;
    l3_list[1]=lt;

    #pragma omp for schedule(dynamic)
    for(ll2=lstart;ll2<=lmax;ll2++) {
      l3_list[0]=ll2;  //Diagonal first, then column
      for(il3=0;il3<2;il3++) {
        int ll3=l3_list[il3];
        int jj,l1,lmin_here,lmax_here;
	    int lmin_00=0,lmax_00=2*(lmax_mask+1)+1;
	    int lmin_ss1=0,lmax_ss1=2*(lmax_mask+1)+1;
	    int lmin_ss2=0,lmax_ss2=2*(lmax_mask+1)+1;
	    int lmin_12=0,lmax_12=2*(lmax_mask+1)+1;
	    int lmin_02=0,lmax_02=2*(lmax_mask+1)+1;
        lmin_here=abs(ll2-ll3);
        lmax_here=ll2+ll3;

	    if(c->has_00 || c->has_0s)
	      drc3jj(ll2,ll3,0,0,&lmin_00,&lmax_00,wigner_00,2*(lmax_mask+1));
	    if(c->has_0s || c->has_ss)
	      drc3jj(ll2,ll3,s1,-s1,&lmin_ss1,&lmax_ss1,wigner_ss1,2*(lmax_mask+1));
          if(has_ss2)
            drc3jj(ll2,ll3,s2,-s2,&lmin_ss2,&lmax_ss2,wigner_ss2,2*(lmax_mask+1));
          else {
            lmin_ss2=lmin_ss1;
            lmax_ss2=lmax_ss1;
          }
	    
	    for(l1=lmin_here;l1<=lmax_here;l1++) {
          int ipp;
          if(l1<=lmax_mask) {
            double wfac;
            double w00=0,wss1=0,wss2=0,w12=0,w02=0;
            int j00=l1-lmin_00;
            int jss1=l1-lmin_ss1;
            int jss2=l1-lmin_ss2;
            if(c->has_00 || c->has_0s)
              w00=j00 < 0 ? 0 : wigner_00[j00];
            if(c->has_ss || c->has_0s) {
              wss1=jss1 < 0 ? 0 : wigner_ss1[jss1];
              wss2=jss2 < 0 ? 0 : wigner_ss2[jss2];
            }
	      
            for(icc=0;icc<npcl;icc++) {
              if(c->has_00) {
                wfac=pcl_masks[l1]*w00*w00;
                tplz_00[icc][il3][ll2]+=wfac;
              }
              if(c->has_0s) {
                wfac=pcl_masks[l1]*wss1*w00;
                tplz_0s[icc][il3][ll2]+=wfac;
              }
              if(c->has_ss) {
                int suml=l1+ll2+ll3;
                wfac=pcl_masks[l1]*wss1*wss2;
	      
                if(suml & 1) //Odd sum
                  tplz_mm[icc][il3][ll2]+=wfac;
                else
                  tplz_pp[icc][il3][ll2]+=wfac;
              }
            }
          }
        }
      }
    } //end omp for
    free(wigner_00);
    free(wigner_ss1);
    if(has_ss2)
      free(wigner_ss2);
  } //end omp parallel

  for(ic=0;ic<npcl;ic++) {
    // Take absolute value of the diagonal to avoid sqrt(-1) later
    for(l2=0;l2<=lmax;l2++) {
      if(c->has_00)
        tplz_00[ic][0][l2]=fabs(tplz_00[ic][0][l2]);
      if(c->has_0s)
        tplz_0s[ic][0][l2]=fabs(tplz_0s[ic][0][l2]);
      if(c->has_ss) {
        tplz_pp[ic][0][l2]=fabs(tplz_pp[ic][0][l2]);
        tplz_mm[ic][0][l2]=fabs(tplz_mm[ic][0][l2]);
      }
    }

    // Compute column correlation coefficient
    for(l2=0;l2<=lmax;l2++) {
      double d1,d2;
      if(c->has_00) {
        d1=tplz_00[ic][0][l2];
        d2=tplz_00[ic][0][lt];
        if((d1>0) && (d2>0))
          tplz_00[ic][1][l2]=tplz_00[ic][1][l2]/sqrt(d1*d2);
        else
          tplz_00[ic][1][l2]=0;
      }
      if(c->has_0s) {
        d1=tplz_0s[ic][0][l2];
        d2=tplz_0s[ic][0][lt];
        if((d1>0) && (d2>0))
          tplz_0s[ic][1][l2]=tplz_0s[ic][1][l2]/sqrt(d1*d2);
        else
          tplz_0s[ic][1][l2]=0;
      }
      if(c->has_ss) {
        d1=tplz_pp[ic][0][l2];
        d2=tplz_pp[ic][0][lt];
        if((d1>0) && (d2>0))
          tplz_pp[ic][1][l2]=tplz_pp[ic][1][l2]/sqrt(d1*d2);
        else
          tplz_pp[ic][1][l2]=0;
        d1=tplz_mm[ic][0][l2];
        d2=tplz_mm[ic][0][lt];
        if((d1>0) && (d2>0))
          tplz_mm[ic][1][l2]=tplz_mm[ic][1][l2]/sqrt(d1*d2);
        else
          tplz_mm[ic][1][l2]=0;
      }
    }

    // Populate matrices
    #pragma omp parallel default(none) shared(c, lmax, ic, lt, tplz_00, tplz_0s, tplz_pp, tplz_mm)
    {
      int ll2, ll3;

      #pragma omp for schedule(dynamic)
      for(ll2=0;ll2<=lmax;ll2++) {
        for(ll3=0;ll3<=ll2;ll3++) {
          int il=toeplitz_wrap(ll2+lt-ll3,lmax+1);
          if(c->has_00)
            xi_00[ic][ll2][ll3]=tplz_00[ic][1][il]*sqrt(tplz_00[ic][0][ll2]*tplz_00[ic][0][ll3]);
          if(c->has_0s)
            c->xi_0s[ic][0][ll2][ll3]=tplz_0s[ic][1][il]*sqrt(tplz_0s[ic][0][ll2]*tplz_0s[ic][0][ll3]);
          if(c->has_ss) {
            c->xi_pp[ic][0][ll2][ll3]=tplz_pp[ic][1][il]*sqrt(tplz_pp[ic][0][ll2]*tplz_pp[ic][0][ll3]);
            c->xi_mm[ic][0][ll2][ll3]=tplz_mm[ic][1][il]*sqrt(tplz_mm[ic][0][ll2]*tplz_mm[ic][0][ll3]);
          }
          if(ll3!=ll2) {
            if(c->has_00)
              c->xi_00[ic][ll3][ll2]=c->xi_00[ic][ll2][ll3];
            if(c->has_0s)
              c->xi_0s[ic][0][ll3][ll2]=c->xi_0s[ic][0][ll2][ll3];
            if(c->has_ss) {
              c->xi_pp[ic][0][ll3][ll2]=c->xi_pp[ic][0][ll2][ll3];
              c->xi_mm[ic][0][ll3][ll2]=c->xi_mm[ic][0][ll2][ll3];
            }
          }
        }
      } //end omp for
    } //end omp parallel
  }

  if(c->has_ss) {
    for(ic=0;ic<npcl;ic++) {
      free(tplz_pp[ic][0]);
      free(tplz_pp[ic][1]);
      free(tplz_pp[ic]);
      free(tplz_mm[ic][0]);
      free(tplz_mm[ic][1]);
      free(tplz_mm[ic]);
    }
    free(tplz_pp);
    free(tplz_mm);
  }
  if(c->has_0s) {
    for(ic=0;ic<npcl;ic++) {
      free(tplz_0s[ic][0]);
      free(tplz_0s[ic][1]);
      free(tplz_0s[ic]);
    }
    free(tplz_0s);
  }
  if(c->has_00) {
    for(ic=0;ic<npcl;ic++) {
      free(tplz_00[ic][0]);
      free(tplz_00[ic][1]);
      free(tplz_00[ic]);
    }
    free(tplz_00);
  }
}
*/

void compute_mcm_namaster(int lmax, int lmax_mask, int npcl, std::vector<double>  pcl_masks, int s1, int s2, int pure_e1, int pure_b1, int pure_e2, int pure_b2, int do_teb, int l_toeplitz, int l_exact, int dl_band, std::vector<std::vector<std::vector<double>>> &xi_00, std::vector<std::vector<std::vector<std::vector<double>>>> &xi_0s, std::vector<std::vector<std::vector<std::vector<double>>>> &xi_pp, std::vector<std::vector<std::vector<std::vector<double>>>> &xi_mm) {
  int ic, ip, ii;
  int pure_any;
  int npure_0s;
  int npure_ss;

  pure_any=pure_e1 || pure_b1 || pure_e2 || pure_b2;

  if(pure_any) {
    npure_0s=2;
    npure_ss=3;
  }
  else {
    npure_0s=1;
    npure_ss=1;
  }

  if(s1==0) {
    if(s2==0) {
      s1=0; s2=0;
    }
    else {
      s1=s2; s2=0;
    }
  }
  else {
    s1=s1;
    s2=s2;
  }


  // if(l_toeplitz>0)
    // populate_toeplitz(c, lmax, lmax_mask, npcl, s1, s2, pcl_masks, l_toeplitz);

  int lstart=0;
  int max_spin=NMT_MAX(s1, s2);
  int has_ss2=(s1!=0) && (s2!=0) && (!do_teb) && (s1!=s2);
  
  lstart=max_spin;

  // #pragma omp parallel default(none) shared(xi_00, xi_0s, xi_pp, xi_mm, s1, s2, lmax, lmax_mask, npcl, pure_any, pure_e1, pure_e2, pure_b1, pure_b2, npure_0s, npure_ss, lstart, do_teb, pcl_masks, has_ss2) shared(l_toeplitz, l_exact, dl_band)
  {
    int ll2,ll3,icc;
    std::vector<double> wigner_00,wigner_ss1,wigner_12,wigner_02,wigner_ss2;
    int pe1=pure_e1,pe2=pure_e2,pb1=pure_b1,pb2=pure_b2;
    wigner_00=std::vector<double>(2*(lmax_mask+1),0.0);
    wigner_ss1=std::vector<double>(2*(lmax_mask+1),0.0);
    wigner_ss2=std::vector<double>(2*(lmax_mask+1),0.0);

    if(pure_any) {
      wigner_12=std::vector<double>(2*(lmax_mask+1),0.0);
      wigner_02=std::vector<double>(2*(lmax_mask+1),0.0);
    }

    // #pragma omp for schedule(dynamic)
    for(ll2=lstart;ll2<=lmax;ll2++) {
      int l3_end=lend_toeplitz(ll2, l_toeplitz, l_exact, dl_band, lmax);
      int l3_start=lstart;
      if(!(pure_any)) //We can use symmetry
        l3_start=ll2;
      for(ll3=l3_start;ll3<=l3_end;ll3++) {
        int jj,l1,lmin_here,lmax_here;
        int lmin_00=0,lmax_00=2*(lmax_mask+1)+1;
        int lmin_ss1=0,lmax_ss1=2*(lmax_mask+1)+1;
        int lmin_ss2=0,lmax_ss2=2*(lmax_mask+1)+1;
        int lmin_12=0,lmax_12=2*(lmax_mask+1)+1;
        int lmin_02=0,lmax_02=2*(lmax_mask+1)+1;
        lmin_here=abs(ll2-ll3);
        lmax_here=ll2+ll3;

        if(l_toeplitz > 0) {
          // Set all elements that will be recomputed to zero
          for(icc=0;icc<npcl;icc++) {
            xi_00[icc][ll2][ll3]=0;
            xi_0s[icc][0][ll2][ll3]=0;
            xi_pp[icc][0][ll2][ll3]=0;
            xi_mm[icc][0][ll2][ll3]=0;
          }
        }

        drc3jj(ll2,ll3,0,0,&lmin_00,&lmax_00,wigner_00,2*(lmax_mask+1));
        drc3jj(ll2,ll3,s1,-s1,&lmin_ss1,&lmax_ss1,wigner_ss1,2*(lmax_mask+1));
        drc3jj(ll2,ll3,s2,-s2,&lmin_ss2,&lmax_ss2,wigner_ss2,2*(lmax_mask+1));

	    if(pure_any) {
	      drc3jj(ll2,ll3,1,-2,&lmin_12,&lmax_12,wigner_12,2*(lmax_mask+1));
	      drc3jj(ll2,ll3,0,-2,&lmin_02,&lmax_02,wigner_02,2*(lmax_mask+1));
	    }
	    
        for(l1=lmin_here;l1<=lmax_here;l1++) {
          int ipp;
          if(l1<=lmax_mask) {
            double wfac,fac_12=0,fac_02=0;
            double w00=0,wss1=0,wss2=0,w12=0,w02=0;
            int j02,j12;
            int j00=l1-lmin_00;
            int jss1=l1-lmin_ss1;
            int jss2=l1-lmin_ss2;
            w00=j00 < 0 ? 0 : wigner_00[j00];
            wss1=jss1 < 0 ? 0 : wigner_ss1[jss1];
            wss2=jss2 < 0 ? 0 : wigner_ss2[jss2];
				    
            if(pure_any) {
	          j12=l1-lmin_12;
	          j02=l1-lmin_02;
	          if(ll2>1.) {
	            fac_12=2*sqrt((l1+1.)*(l1+0.)/((ll2+2)*(ll2-1.)));
	            if(l1>1.) fac_02=sqrt((l1+2.)*(l1+1.)*(l1+0.)*(l1-1.)/((ll2+2.)*(ll2+1.)*(ll2+0.)*(ll2-1.)));
	            else fac_02=0;
	          }
	          else {
	            fac_12=0;
	            fac_02=0;
	          }
	          if(j12<0) { //If out of range, w12 is just 0
	            fac_12=0;
	            j12=0;
	          }
	          if(j02<0) { //if out of range, w02 is just 0
	            fac_02=0;
	            j02=0;
	          }
              w12=j12 < 0 ? 0 : wigner_12[j12];
              w02=j02 < 0 ? 0 : wigner_02[j02];
	        }
	    
            for(icc=0;icc<npcl;icc++) {
              wfac=pcl_masks[l1]*w00*w00;
              xi_00[icc][ll2][ll3]+=wfac;
              double wfac_ispure[2];
              wfac_ispure[0]=wss1;
              wfac_ispure[0]*=pcl_masks[l1]*w00;
              if(pure_any) {
                wfac_ispure[1]=wss1+fac_12*w12+fac_02*w02;
                wfac_ispure[1]*=pcl_masks[l1]*w00;
              }
              for(ipp=0;ipp<npure_0s;ipp++) xi_0s[icc][ipp][ll2][ll3]+=wfac_ispure[ipp];
              double wfac_ispure2[3];
              int suml=l1+ll2+ll3;
              wfac_ispure2[0]=wss1;
              wfac_ispure2[0]*=wss2*pcl_masks[l1];
              if(pure_any) {
                wfac_ispure2[1]=wss1+fac_12*w12+fac_02*w02;
                wfac_ispure2[2]=wfac_ispure2[1]*wfac_ispure2[1]*pcl_masks[l1];
                wfac_ispure2[1]*=wss2*pcl_masks[l1];
              }
	    
              if(suml & 1) { //Odd sum
                for(ipp=0;ipp<npure_ss;ipp++)
                  xi_mm[icc][ipp][ll2][ll3]+=wfac_ispure2[ipp];
              }
              else {
                for(ipp=0;ipp<npure_ss;ipp++)
                  xi_pp[icc][ipp][ll2][ll3]+=wfac_ispure2[ipp];
              }
	        }
	      }
	    }

        if((!(pure_any)) && (ll2 != ll3)) { //Can use symmetry
          for(icc=0;icc<npcl;icc++) {
            xi_00[icc][ll3][ll2]=xi_00[icc][ll2][ll3];
            xi_0s[icc][0][ll3][ll2]=xi_0s[icc][0][ll2][ll3];
            xi_pp[icc][0][ll3][ll2]=xi_pp[icc][0][ll2][ll3];
            xi_mm[icc][0][ll3][ll2]=xi_mm[icc][0][ll2][ll3];
          }
        }
      }
    } //end omp for
  } //end omp parallel

  // Fill out lower triangle
  if(l_toeplitz > 0) {
    int l2, l3;
    for(l2=lmax+l_exact-l_toeplitz;l2<=lmax;l2++) {
      for(l3=l_exact;l3<=l2+l_toeplitz-lmax;l3++){
        std::vector<std::vector<double>> mat;
        double m;
        int lx=l_exact+l2-l3;
        for(ic=0;ic<npcl;ic++) {
          mat = xi_00[ic];
          m=mat[lx][l_exact]*sqrt(fabs(mat[l2][l2]*mat[l3][l3]/(mat[lx][lx]*mat[l_exact][l_exact])));
          mat[l2][l3]=m;
          mat[l3][l2]=m;
          mat = xi_0s[ic][0];
          m=mat[lx][l_exact]*sqrt(fabs(mat[l2][l2]*mat[l3][l3]/(mat[lx][lx]*mat[l_exact][l_exact])));
          mat[l2][l3]=m;
          mat[l3][l2]=m;
          mat = xi_pp[ic][0];
          m=mat[lx][l_exact]*sqrt(fabs(mat[l2][l2]*mat[l3][l3]/(mat[lx][lx]*mat[l_exact][l_exact])));
          mat[l2][l3]=m;
          mat[l3][l2]=m;
          mat = xi_mm[ic][0];
          m=mat[lx][l_exact]*sqrt(fabs(mat[l2][l2]*mat[l3][l3]/(mat[lx][lx]*mat[l_exact][l_exact])));
          mat[l2][l3]=m;
          mat[l3][l2]=m;
        }
      }
    }
  }
}