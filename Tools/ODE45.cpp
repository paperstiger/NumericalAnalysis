#include"ODE45.h"
#define abs(x) (((x) > 0) ? (x) : (-(x)))
double eps(double t)
{
	if(t==0)
		return pow(2.0, -1074.0);
	else
		return pow(2.0, -52.0+floor(log(fabs(t))/log(2.0)));
}
//Çó1,2,inf-·¶Êý
double norm(const double* vec, int n, int type)
{
	int i;
	double res=0.0;
	if(type==1)
	{
		for(i=0;i<n;i++)
		{
			if(vec[i]>=0)
				res+=vec[i];
			else
				res-=vec[i];
		}
	}
	else if(type==2)
	{
		for(i=0;i<n;i++)
			res+=vec[i]*vec[i];
		res=sqrt(res);
	}
	else
	{
		for(i=0;i<n;i++) 
		{
			if(vec[i]>=0.0)
			{
				if(vec[i]>res)
					res=vec[i];
			}
			else
			{
				if(-vec[i]>res)
					res=-vec[i];
			}
		}
	}
	return res;
}

int ode45nw(int (*fun)(double t, const double* x, double* dx, const double* para), double* x, const double* para,
		  double t0, double tf, int n, int& NumPoint, const double* AbsTol, double RelTol,
		  int NormControl, double MaxStep, double InitialStep, FILE* fid)
{
	// Stats
	int flag=0, odeflag=0;
	int nsteps, nfailed, nfevals, tdir, nofailed;
	double rtol, normx, t, powth, hmax, hmin, absh, rh, h, tnew, normxnew, errwt, err, temp;
	int i, j, k;
	double* work[10];
	for(i=0;i<10;i++)
		work[i]=new double[n];
	for(i=0;i<10;i++)
		for(j=0;j<n;j++)
			work[i][j]=0.0;
	double* threshold=work[7];
	double* xnew=work[9];
	double hB[7][6];
	nsteps  = 0;
	nfailed = 0;
	nfevals = 0; 
	static const double A[6]={1.0/5.0, 3.0/10.0, 4.0/5.0, 8.0/9.0, 1.0, 1.0};
	static const double B[7][6]=
	{	{1.0/5.0,	3.0/40.0,	44.0/45.0,	19372.0/6561.0,		9017.0/3168.0,		35.0/384.0		},
		{0.0,		9.0/40.0,	-56.0/15.0,	-25360.0/2187.0,	-355.0/33.0,		0.0				},
		{0.0,		0.0,		32.0/9.0,	64448.0/6561.0,		46732.0/5247.0,		500.0/1113.0	},
		{0.0,		0.0,		0.0,		-212.0/729.0,		49.0/176.0,			125.0/192.0		},
		{0.0,		0.0,		0.0,		0.0,				-5103.0/18656.0,	-2187.0/6784.0	},
		{0.0,		0.0,		0.0,		0.0,				0.0,				11.0/84.0		},
		{0.0,		0.0,		0.0,		0.0,				0.0,				0.0				}};
	static const double E[7]={71.0/57600.0, 0.0, -71.0/16695.0, 71.0/1920.0, -17253.0/339200.0, 22.0/525.0, -1.0/40.0};
	// Handle solver arguments
	tdir=sign(tf-t0);
	flag=fun(t0,x,work[0],para);	
	if(flag<1)
	{ odeflag=0;goto Final;}
	rtol=RelTol;
	if(rtol<100.0*eps(1.0))
		rtol=100.0*eps(1.0);
	normx=0.0;
	threshold[0]=AbsTol[0]/rtol;
	if(NormControl)
		normx=norm(x,n,2);
	else
		for(i=1;i<n;i++)
			threshold[i]=AbsTol[i]/rtol;
	if(MaxStep>0.0)
		hmax=min(fabs(tf-t0),MaxStep);
	else
		hmax=min(fabs(tf-t0), fabs(0.1*(tf-t0)));
	
	nfevals=nfevals+1;
	
	t=t0;
	// Allocate memory if we're generating output.
	// alloc in chunks
	NumPoint = 1;
	if(fid)
	{
		fprintf(fid,"%22.14e",t);
		for(i=0;i<n;i++)
			fprintf(fid,"%24.14e",x[i]);
		fflush(fid);
	}

	// Initialize method parameters.
	powth=1.0/5.0;
	
	hmin=16.0*eps(t);
	if(InitialStep<=0.0)
	{
		// Compute an initial step size h using y'(t).
		absh=min(hmax, fabs(tf-t0));
		if(NormControl)
			rh=(norm(work[0],n,2)/max(normx,threshold[0]))/(0.8*pow(rtol, powth));		
		else
		{
			for(i=0;i<n;i++)
				xnew[i]=work[0][i]/max(fabs(x[i]), threshold[i]);
			rh=norm(xnew, n, 3)/(0.8*pow(rtol, powth));			
		}
		if(absh*rh>1.0)
			absh=1.0/rh;
		absh=max(absh, hmin);
	}
	else
		absh=min(hmax, max(hmin, InitialStep));
	
	// THE MAIN LOOP
	int done;
	done=0;
	while(!done)
	{
		// By default, hmin is a small number such that t+hmin is only slightly
		// different than t.  It might be 0 if t is 0.
		hmin = 16.0*eps(t);
		absh = min(hmax, max(hmin, absh));    // couldn't limit absh until new hmin
		h = tdir * absh;

		// Stretch the step if within 10// of tf-t.
		if(1.1*absh >= fabs(tf - t))
		{
			h = tf - t;
			absh = fabs(h);
			done = 1;
		}
 
		// LOOP FOR ADVANCING ONE STEP.
		nofailed = 1;                      // no failed attempts
		while(1)
		{
			for(i=0;i<7;i++)
				for(j=i;j<6;j++)
					hB[i][j]=h*B[i][j];
			for(k=0;k<5;k++)
			{
				for(i=0;i<n;i++)
				{
					work[8][i]=x[i];
					for(j=0;j<=k;j++) work[8][i]+=work[j][i]*hB[j][k];					
				}
				flag=fun(t+h*A[k],work[8],work[k+1],para);
				if(flag<1){ odeflag=0;goto Final;}				
			}			
			//f(:,2) = feval(fun,t+h*A(1),y+f*hB(:,1),para);
			//f(:,3) = feval(fun,t+h*A(2),y+f*hB(:,2),para);
			//f(:,4) = feval(fun,t+h*A(3),y+f*hB(:,3),para);
			//f(:,5) = feval(fun,t+h*A(4),y+f*hB(:,4),para);
			//f(:,6) = feval(fun,t+h*A(5),y+f*hB(:,5),para);
			
			tnew=t+h*A[5];
			if(done)
				tnew=tf;   // Hit end point exactly.
			//h = tnew - t;      // Purify h.
			
			for(i=0;i<n;i++)
			{
				xnew[i]=x[i];
				for(j=0;j<=5;j++) xnew[i]+=work[j][i]*hB[j][5];				
			}
			flag=fun(tnew,xnew,work[6],para);
			if(flag<1) { odeflag=0;goto Final;}			
			
			nfevals = nfevals + 6;              

			// Estimate the error.
			for(i=0;i<n;i++)
			{
				work[8][i]=0.0;
				for(j=0;j<7;j++) work[8][i]+=work[j][i]*E[j];					
			}
			if(NormControl)
			{
				normxnew = norm(xnew, n, 2);
				errwt = max(max(normx,normxnew),threshold[0]);				
				err = absh * (norm(work[8], n, 2) / errwt);
			}
			else
			{
				for(i=0;i<n;i++)
					work[8][i]/=max(max(fabs(x[i]),fabs(xnew[i])),threshold[i]);
				err = absh * norm(work[8], n, 3);              
			}

			// Accept the solution only if the weighted error is no more than the
			// tolerance rtol.  Estimate an h that will yield an error of rtol on
			// the next step or the next try at taking this step, as the case may be,
			// and use 0.8 of this value to avoid failures.
			if(err>rtol)                      // Failed step
			{
				nfailed = nfailed + 1;            
				if(absh <= hmin)
				{
					printf("Failure at t=%e.,...Unable to meet integration tolerances without reducing the step size below the smallest value allowed (%e) at time t.",t,hmin);
					odeflag=0;
					goto Final;
				}

				if(nofailed)
				{
					nofailed = false;
					absh = max(hmin, absh * max(0.1, 0.8*pow(rtol/err, powth))); 
				}
				else
					absh = max(hmin, 0.5 * absh);
				h = tdir * absh;
				done = false;
			}
			else        // Successful step
				break;      
      
		}
		
		nsteps = nsteps + 1;
		
		NumPoint = NumPoint + 1; 
		if(fid)
		{
			fprintf(fid,"\n");
			fprintf(fid,"%22.14e",tnew);
			for(i=0;i<n;i++)
				fprintf(fid,"%24.14e",xnew[i]);			
			fflush(fid);
		}
		
		if(done)
			break;
		

		// If there were no failures compute a new h.
		if(nofailed)
		{
			// Note that absh may shrink by 0.8, and that err may be 0.
			temp = 1.25*pow(err/rtol,powth);
			if(temp > 0.2)
				absh = absh / temp;
			else
				absh = 5.0*absh;			
		}

		// Advance the integration one step.
		t = tnew;
		for(i=0;i<n;i++) x[i] = xnew[i];
		if(NormControl)
			normx = normxnew;
		for(i=0;i<n;i++) work[0][i]= work[6][i];  // Already have f(tnew,ynew)
	}
	for(i=0;i<n;i++) x[i] = xnew[i];	
	odeflag=1;
Final:
	for(i=0;i<10;i++)
		delete[] work[i];	
	return odeflag;
}


int ode45(int (*fun)(double t, const double* x, double* dx, const double* para), double* x, const double* para,
		  double t0, double tf, int n, int& NumPoint, double* work, const double* AbsTol, double RelTol,
		  int MaxIter, int NormControl, double MaxStep, double InitialStep, FILE* fid)
{
	// Stats
	//printf("Enter ode45\n");
	int flag=0, odeflag=0;
	int nsteps, nfailed, nfevals, tdir, nofailed;
	double rtol, normx, t, powth, hmax, hmin, absh, rh, h, tnew, normxnew, errwt, err, temp;
	int i, j, k;		
	double* xnew=&work[9*n];
	double hB[7][6]={0.0};
	nsteps  = 0;
	nfailed = 0;
	nfevals = 0; 
	static const double A[6]={1.0/5.0, 3.0/10.0, 4.0/5.0, 8.0/9.0, 1.0, 1.0};
	static const double B[7][6]=
	{	{1.0/5.0,	3.0/40.0,	44.0/45.0,	19372.0/6561.0,		9017.0/3168.0,		35.0/384.0		},
		{0.0,		9.0/40.0,	-56.0/15.0,	-25360.0/2187.0,	-355.0/33.0,		0.0				},
		{0.0,		0.0,		32.0/9.0,	64448.0/6561.0,		46732.0/5247.0,		500.0/1113.0	},
		{0.0,		0.0,		0.0,		-212.0/729.0,		49.0/176.0,			125.0/192.0		},
		{0.0,		0.0,		0.0,		0.0,				-5103.0/18656.0,	-2187.0/6784.0	},
		{0.0,		0.0,		0.0,		0.0,				0.0,				11.0/84.0		},
		{0.0,		0.0,		0.0,		0.0,				0.0,				0.0				}};
	static const double E[7]={71.0/57600.0, 0.0, -71.0/16695.0, 71.0/1920.0, -17253.0/339200.0, 22.0/525.0, -1.0/40.0};
	// Handle solver arguments
	tdir=sign(tf-t0);
	flag=fun(t0,x, work,para);	
	if(flag<1)
	{ odeflag=0;goto Final;}
	rtol=RelTol;
	if(rtol<100.0*eps(1.0))
		rtol=100.0*eps(1.0);
	normx=0.0;
	work[7*n]=AbsTol[0]/rtol;
	if(NormControl)
		normx=norm(x,n,2);
	else
		for(i=1;i<n;i++)
			work[7*n+i]=AbsTol[i]/rtol;
	if(MaxStep>0.0)
		hmax=min(fabs(tf-t0),MaxStep);
	else
		hmax=min(fabs(tf-t0), fabs(0.1*(tf-t0)));
	
	nfevals=nfevals+1;
	
	t=t0;
	// Allocate memory if we're generating output.
	// alloc in chunks
	NumPoint = 1;
	if(fid)
	{
		fprintf(fid,"%22.14e",t);
		for(i=0;i<n;i++)
			fprintf(fid,"%24.14e",x[i]);
		fflush(fid);
	}

	// Initialize method parameters.
	powth=1.0/5.0;
	hmin=16.0*eps(t);
	if(InitialStep<=0.0)
	{
		// Compute an initial step size h using y'(t).
		absh=min(hmax, fabs(tf-t0));
		if(NormControl)
			rh=(norm(work,n,2)/max(normx,work[7*n]))/(0.8*pow(rtol, powth));		
		else
		{
			for(i=0;i<n;i++)
				xnew[i]=work[i]/max(fabs(x[i]), work[7*n+i]);
			rh=norm(xnew, n, 3)/(0.8*pow(rtol, powth));
		}
			/*rh=norm(work,n,3)/max(fabs(x[i]),work[7*n+i])/(0.8*pow(rtol, powth));*/
		if(absh*rh>1.0)
			absh=1.0/rh;
		absh=max(absh, hmin);
	}
	else
		absh=min(hmax, max(hmin, InitialStep));
	
	// THE MAIN LOOP
	int done;
	done = 0;
	while(!done)
	{
		// By default, hmin is a small number such that t+hmin is only slightly
		// different than t.  It might be 0 if t is 0.
		hmin = 16.0*eps(t);
		absh = min(hmax, max(hmin, absh));    // couldn't limit absh until new hmin
		h = tdir * absh;

		// Stretch the step if within 10// of tf-t.
		if(1.1*absh >= fabs(tf - t))
		{
			h = tf - t;
			absh = fabs(h);
			done = 1;
		}
  
		// LOOP FOR ADVANCING ONE STEP.
		nofailed = 1;                      // no failed attempts
		while(1)
		{
			for(i=0;i<7;i++)
				for(j=i;j<6;j++)
					hB[i][j]=h*B[i][j];
			for(k=0;k<5;k++)
			{
				for(i=0;i<n;i++)
				{
					work[8*n+i]=x[i];
					for(j=0;j<=k;j++) work[8*n+i]+=work[j*n+i]*hB[j][k];					
				}
				flag=fun(t+h*A[k],&work[8*n],&work[n*(k+1)],para);
				if(flag<1){ odeflag=0;goto Final;}				
			}			
			//f(:,2) = feval(fun,t+h*A(1),y+f*hB(:,1),para);
			//f(:,3) = feval(fun,t+h*A(2),y+f*hB(:,2),para);
			//f(:,4) = feval(fun,t+h*A(3),y+f*hB(:,3),para);
			//f(:,5) = feval(fun,t+h*A(4),y+f*hB(:,4),para);
			//f(:,6) = feval(fun,t+h*A(5),y+f*hB(:,5),para);
			
			tnew=t+h*A[5];
			if(done)
				tnew=tf;   // Hit end point exactly.
			//h = tnew - t;      // Purify h.
			
			for(i=0;i<n;i++)
			{
				xnew[i]=x[i];
				for(j=0;j<=5;j++) xnew[i]+=work[j*n+i]*hB[j][5];				
			}
			flag=fun(tnew,xnew,&work[6*n],para);
			if(flag<1) { odeflag=0;goto Final;}			
			
			nfevals = nfevals + 6;              

			// Estimate the error.
			for(i=0;i<n;i++)
			{
				work[8*n+i]=0.0;
				for(j=0;j<7;j++) work[8*n+i]+=work[j*n+i]*E[j];					
			}
			if(NormControl)
			{
				normxnew = norm(xnew, n, 2);
				errwt = max(max(normx,normxnew),work[7*n]);				
				err = absh * (norm(&work[8*n], n, 2) / errwt);
			}
			else
			{
				for(i=0;i<n;i++)
					work[8*n+i]/=max(max(fabs(x[i]),fabs(xnew[i])),work[7*n+i]);
				err = absh * norm(&work[8*n], n, 3);              
			}

			// Accept the solution only if the weighted error is no more than the
			// tolerance rtol.  Estimate an h that will yield an error of rtol on
			// the next step or the next try at taking this step, as the case may be,
			// and use 0.8 of this value to avoid failures.
			if(err>rtol)                      // Failed step
			{
				nfailed = nfailed + 1;            
				if(absh <= hmin)
				{
					printf("Failure at t=%e.,...Unable to meet integration tolerances without reducing the step size below the smallest value allowed (%e) at time t.",t,hmin);
					odeflag=0;
					goto Final;
				}

				if(nofailed)
				{
					nofailed = false;
					absh = max(hmin, absh * max(0.1, 0.8*pow(rtol/err, powth))); 
				}
				else
					absh = max(hmin, 0.5 * absh);
				h = tdir * absh;
				done = false;
			}
			else        // Successful step
				break;      
       
		}
		
		nsteps = nsteps + 1;
		if(nsteps>MaxIter)
		{odeflag=0;goto Final;}
		
		NumPoint = NumPoint + 1; 
		if(fid)
		{
			fprintf(fid,"\n");
			fprintf(fid,"%22.14e",tnew);
			for(i=0;i<n;i++)
				fprintf(fid,"%24.14e",xnew[i]);			
			fflush(fid);
		}
		
		if(done)
			break;
		

		// If there were no failures compute a new h.
		if(nofailed)
		{
			// Note that absh may shrink by 0.8, and that err may be 0.
			temp = 1.25*pow(err/rtol,powth);
			if(temp > 0.2)
				absh = absh / temp;
			else
				absh = 5.0*absh;			
		}

		// Advance the integration one step.
		t = tnew;
		for(i=0;i<n;i++) x[i] = xnew[i];
		if(NormControl)
			normx = normxnew;
		for(i=0;i<n;i++) work[i]= work[6*n+i];  // Already have f(tnew,ynew)
	}
	for(i=0;i<n;i++) x[i] = xnew[i];	
	odeflag=1;
Final:	
	//printf("Leave ode45\n");
	return odeflag;
}
//One subroutine that is often used and make it a subroutine
//Input h, x , calc xnew
void onestep(double t0, double tf, int n, double *work, double *x, int (*fun)(double t, const double* x, double* dx, const double* para), const double *A, const double *para, double *xnew)
{
	double t = t0, h = tf - t0;
	double hB[7][6]={0.0};
	static const double B[7][6]=
	{	{1.0/5.0,	3.0/40.0,	44.0/45.0,	19372.0/6561.0,		9017.0/3168.0,		35.0/384.0		},
		{0.0,		9.0/40.0,	-56.0/15.0,	-25360.0/2187.0,	-355.0/33.0,		0.0				},
		{0.0,		0.0,		32.0/9.0,	64448.0/6561.0,		46732.0/5247.0,		500.0/1113.0	},
		{0.0,		0.0,		0.0,		-212.0/729.0,		49.0/176.0,			125.0/192.0		},
		{0.0,		0.0,		0.0,		0.0,				-5103.0/18656.0,	-2187.0/6784.0	},
		{0.0,		0.0,		0.0,		0.0,				0.0,				11.0/84.0		},
		{0.0,		0.0,		0.0,		0.0,				0.0,				0.0				}};
	for(int i=0;i<7;i++)
		for(int j=i;j<6;j++)
			hB[i][j]=h*B[i][j];
	for(int k=0;k<5;k++)
	{
		for(int i=0;i<n;i++)
		{
			work[8*n+i]=x[i];
			for(int j=0;j<=k;j++) work[8*n+i]+=work[j*n+i]*hB[j][k];					
		}
		int flag=fun(t+h*A[k],&work[8*n],&work[n*(k+1)],para);			
	}			
	for(int i=0;i<n;i++)
	{
		xnew[i]=x[i];
		for(int j=0;j<=5;j++) xnew[i]+=work[j*n+i]*hB[j][5];				
	}
}
//ODE45 by myself which can handle events
//Another fun to determine the satisfactory of condition is needed and named condfun
//int(*eventfun)(t, x, para),this fun return 0 before the state changed(when t=t0,ignore it)
//Note that event usually equal to 0, while ode process may jump over it, find new time
//when condfun return 1, and use False Method to predict the time, tol is satisfied of course
//eventflag == 1 if event satisfied

//I think that USING NEWTON iter is faster because analytic grad is easy to calc
int ode45(int (*fun)(double t, const double* x, double* dx, const double* para), int (*eventfun)(double t, const double *x, const double *eventPara, double &fval, double &fgrad), double* x, const double* para,
 const double *eventPara, const double eventTol, int &eventFlag, 
		  double t0, double &tf, int n, int& NumPoint, double* work, const double* AbsTol, double RelTol,
		  int MaxIter, int NormControl, double MaxStep, double InitialStep, FILE* fid)
{
	// Stats
	//printf("Enter ode45\n");
	int flag=0, odeflag=0, tmpEventFlag = 0, bisIter = 0;
	int nsteps, nfailed, nfevals, tdir, nofailed;
	double rtol, normx, t, powth, hmax, hmin, absh, rh, h, tnew, normxnew, errwt, err, temp;
	double eventFval, eventGrad;
	int i, j, k;		
	double* xnew=&work[9*n];
	double hB[7][6]={0.0};
	nsteps  = 0;
	nfailed = 0;
	nfevals = 0; 
	static const double A[6]={1.0/5.0, 3.0/10.0, 4.0/5.0, 8.0/9.0, 1.0, 1.0};
	static const double B[7][6]=
	{	{1.0/5.0,	3.0/40.0,	44.0/45.0,	19372.0/6561.0,		9017.0/3168.0,		35.0/384.0		},
		{0.0,		9.0/40.0,	-56.0/15.0,	-25360.0/2187.0,	-355.0/33.0,		0.0				},
		{0.0,		0.0,		32.0/9.0,	64448.0/6561.0,		46732.0/5247.0,		500.0/1113.0	},
		{0.0,		0.0,		0.0,		-212.0/729.0,		49.0/176.0,			125.0/192.0		},
		{0.0,		0.0,		0.0,		0.0,				-5103.0/18656.0,	-2187.0/6784.0	},
		{0.0,		0.0,		0.0,		0.0,				0.0,				11.0/84.0		},
		{0.0,		0.0,		0.0,		0.0,				0.0,				0.0				}};
	static const double E[7]={71.0/57600.0, 0.0, -71.0/16695.0, 71.0/1920.0, -17253.0/339200.0, 22.0/525.0, -1.0/40.0};
	// Handle solver arguments
	tdir=sign(tf-t0);
	//eval fun at t0, what's the point??No eventfun eval at t0
	flag=fun(t0,x, work,para);	
	if(flag<1)
	{ odeflag=0;goto Final;}
	rtol=RelTol;
	if(rtol<100.0*eps(1.0))
		rtol=100.0*eps(1.0);
	normx=0.0;
	work[7*n]=AbsTol[0]/rtol;
	if(NormControl)
		normx=norm(x,n,2);
	else
		for(i=1;i<n;i++)
			work[7*n+i]=AbsTol[i]/rtol;
	if(MaxStep>0.0)
		hmax=min(fabs(tf-t0),MaxStep);
	else
		hmax=min(fabs(tf-t0), fabs(0.1*(tf-t0)));
	
	nfevals=nfevals+1;
	
	t=t0;
	// Allocate memory if we're generating output.
	// alloc in chunks
	NumPoint = 1;
	if(fid)
	{
		fprintf(fid,"%22.14e",t);
		for(i=0;i<n;i++)
			fprintf(fid,"%24.14e",x[i]);
		fflush(fid);
	}

	// Initialize method parameters.
	powth=1.0/5.0;
	hmin=16.0*eps(t);
	if(InitialStep<=0.0)
	{
		// Compute an initial step size h using y'(t).
		absh=min(hmax, fabs(tf-t0));
		if(NormControl)
			rh=(norm(work,n,2)/max(normx,work[7*n]))/(0.8*pow(rtol, powth));		
		else
		{
			for(i=0;i<n;i++)
				xnew[i]=work[i]/max(fabs(x[i]), work[7*n+i]);
			rh=norm(xnew, n, 3)/(0.8*pow(rtol, powth));
		}
			/*rh=norm(work,n,3)/max(fabs(x[i]),work[7*n+i])/(0.8*pow(rtol, powth));*/
		if(absh*rh>1.0)
			absh=1.0/rh;
		absh=max(absh, hmin);
	}
	else
		absh=min(hmax, max(hmin, InitialStep));
	
	// THE MAIN LOOP
	int done;
	done = 0;
	while(!done)
	{
		// By default, hmin is a small number such that t+hmin is only slightly
		// different than t.  It might be 0 if t is 0.
		hmin = 16.0*eps(t);
		absh = min(hmax, max(hmin, absh));    // couldn't limit absh until new hmin
		h = tdir * absh;

		// Stretch the step if within 10// of tf-t.
		if(1.1*absh >= fabs(tf - t))
		{
			h = tf - t;
			absh = fabs(h);
			done = 1;
		}
  
		// LOOP FOR ADVANCING ONE STEP.
		nofailed = 1;                      // no failed attempts
		while(1)
		{
			for(i=0;i<7;i++)
				for(j=i;j<6;j++)
					hB[i][j]=h*B[i][j];
			for(k=0;k<5;k++)
			{
				for(i=0;i<n;i++)
				{
					work[8*n+i]=x[i];
					for(j=0;j<=k;j++) work[8*n+i]+=work[j*n+i]*hB[j][k];					
				}
				flag=fun(t+h*A[k],&work[8*n],&work[n*(k+1)],para);
				if(flag<1){ odeflag=0;goto Final;}				
			}			
			//f(:,2) = feval(fun,t+h*A(1),y+f*hB(:,1),para);
			//f(:,3) = feval(fun,t+h*A(2),y+f*hB(:,2),para);
			//f(:,4) = feval(fun,t+h*A(3),y+f*hB(:,3),para);
			//f(:,5) = feval(fun,t+h*A(4),y+f*hB(:,4),para);
			//f(:,6) = feval(fun,t+h*A(5),y+f*hB(:,5),para);
			
			tnew=t+h*A[5];
			if(done)
				tnew=tf;   // Hit end point exactly.
			//h = tnew - t;      // Purify h.
			
			for(i=0;i<n;i++)
			{
				xnew[i]=x[i];
				for(j=0;j<=5;j++) xnew[i]+=work[j*n+i]*hB[j][5];				
			}
			flag=fun(tnew,xnew,&work[6*n],para);
			if(flag<1) { odeflag=0;goto Final;}			
			
			nfevals = nfevals + 6;              

			// Estimate the error.
			for(i=0;i<n;i++)
			{
				work[8*n+i]=0.0;
				for(j=0;j<7;j++) work[8*n+i]+=work[j*n+i]*E[j];					
			}
			if(NormControl)
			{
				normxnew = norm(xnew, n, 2);
				errwt = max(max(normx,normxnew),work[7*n]);				
				err = absh * (norm(&work[8*n], n, 2) / errwt);
			}
			else
			{
				for(i=0;i<n;i++)
					work[8*n+i]/=max(max(fabs(x[i]),fabs(xnew[i])),work[7*n+i]);
				err = absh * norm(&work[8*n], n, 3);              
			}

			// Accept the solution only if the weighted error is no more than the
			// tolerance rtol.  Estimate an h that will yield an error of rtol on
			// the next step or the next try at taking this step, as the case may be,
			// and use 0.8 of this value to avoid failures.
			if(err>rtol)                      // Failed step
			{
				nfailed = nfailed + 1;            
				if(absh <= hmin)
				{
					printf("Failure at t=%e.,...Unable to meet integration tolerances without reducing the step size below the smallest value allowed (%e) at time t.",t,hmin);
					odeflag=0;
					goto Final;
				}

				if(nofailed)
				{
					nofailed = false;
					absh = max(hmin, absh * max(0.1, 0.8*pow(rtol/err, powth))); 
				}
				else
					absh = max(hmin, 0.5 * absh);
				h = tdir * absh;
				done = false;
			}
			else        // Successful step
				break;      
       
		}
		
		nsteps = nsteps + 1;
		if(nsteps>MaxIter)
		{odeflag=0;goto Final;}
		
		NumPoint = NumPoint + 1;
		//Now we are sure this step is ok, calc eventfun at tnew
		tmpEventFlag = eventfun(tnew, xnew, eventPara, eventFval, eventGrad);
		if(abs(eventFval) < eventTol)
		{
			done = 1;
			eventFlag = 1;
			tf = tnew;
			break;
		}
		//next step should be calced carefully
		if(tmpEventFlag)
		{
			double middletime;
			double *xmiddle = new double[n];
			//First calc value at t and x
			//tmpEventFlag = eventfun(t, x, eventPara, eventFval, eventGrad);
			double t0 = t;
			middletime = tnew;
			//First try Using Newton method from tnew because i don't want to handle tstep 
			for(int newtonIter = 0;newtonIter < 20;newtonIter++)
			{
				middletime -= eventFval/eventGrad;
				onestep(t0, middletime, n, work, x, fun, A, para, xmiddle);
				tmpEventFlag = eventfun(t, x, eventPara, eventFval, eventGrad);
				if(abs(eventFval) < eventTol)
				{
					done = 1;
					eventFlag = 1;
					tnew = middletime;
					tf = tnew;
					for(int ii = 0;ii < n;ii++) xnew[ii] = xmiddle[ii];
					delete[] xmiddle;
					goto tofile;
				}
			}
			
			
			for(bisIter = 0;bisIter <= 50;bisIter++)
			{
				middletime = (t + tnew)/2.0;
				//First integrate to middletime so we can get accurate xmiddle
				onestep(t0, middletime, n, work, x, fun, A, para, xmiddle);
				tmpEventFlag = eventfun(middletime, xmiddle, eventPara, eventFval, eventGrad);
				if(abs(eventFval) < eventTol)
				{
					done = 1;
					eventFlag = 1;
					tnew = middletime;
					tf = tnew;
					for(int ii = 0;ii < n;ii++) xnew[ii] = xmiddle[ii];
					break;
				}
				else
				{
					//determine which time be changed,Eventflag=1 means closer with tnew
					if(tmpEventFlag == 1)
					{
						tnew = middletime;
					}
					else
					{
						t = middletime;
					}
				}
			}
			delete[] xmiddle;
		}
		if(bisIter == 21)
		{
			done = 1;
			printf("Failed to find event point\n");
		}
tofile:
		if(fid)
		{
			fprintf(fid,"\n");
			fprintf(fid,"%22.14e",tnew);
			for(i=0;i<n;i++)
				fprintf(fid,"%24.14e",xnew[i]);			
			fflush(fid);
		}
		
		if(done)
			break;
		

		// If there were no failures compute a new h.
		if(nofailed)
		{
			// Note that absh may shrink by 0.8, and that err may be 0.
			temp = 1.25*pow(err/rtol,powth);
			if(temp > 0.2)
				absh = absh / temp;
			else
				absh = 5.0*absh;			
		}

		// Advance the integration one step.
		t = tnew;
		for(i=0;i<n;i++) x[i] = xnew[i];
		if(NormControl)
			normx = normxnew;
		for(i=0;i<n;i++) work[i]= work[6*n+i];  // Already have f(tnew,ynew)
	}
	for(i=0;i<n;i++) x[i] = xnew[i];	
	odeflag=1;
	if(eventFlag)
		odeflag = 100;
Final:	
	//printf("Leave ode45\n");
	return odeflag;
}
