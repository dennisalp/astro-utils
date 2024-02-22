#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <string.h>
#include <complex.h>

/*
2018-05-19, Dennis Alp, dalp@kth.se
Detection of pulsation in photon counting data. 
1st argument: file name
2nd argument: number of events? Recommended to launch with PATH2FILE
$(wc -l < PATH2FILE) as arguments. The ./run script should take care
of this.
*/

#define MEMORY 1000 // GB of ram
#define OVERLAP 10
#define ROOT 0
#define DBG_LEVEL 1
#define ND 10000
#define DIS_RES 100
#define TOP 0
#define NTOP 10
#define SPARSE 10.
#define PRINT_RATE 100000
#define PMIN 0.03
#define PMAX 1.e2
#define TREF 0.
#define SPIN_DOWN_LIM 6.042267201371935e-117
#define AGE_LIM 5e-119
#define MMAX 1

#define INFO(STR, pp, LL) \
  if (DBG_LEVEL >= LL) printf("ID=%d: %s\n", pp, STR);
#define GET_NP(n_tot, PP, pp) \
  np = n_tot/PP;				\
  np += pp < n_tot%PP ? 1 : 0;

void send_mail(char *ff, const double ptmp, const double pdot, const double HH) {
  char cmd[100];  // to hold the command.
  char to[] = "example@gmail.com"; // email id of the recepient.
  char body[999];
  sprintf(body, "Subject: ALERT cas.c\nTo: example@gmail.com\nFrom:87a.c\n%s\nPP=%.17e\nPd=%.17e\nHH=%.17e", ff, ptmp, pdot, HH);    // email body.
  printf("%lu", sizeof(body)/sizeof(char));
  char tmp[] = "mail.tmp";
 
  FILE *fp = fopen(tmp,"w"); // open it for writing.
  fprintf(fp,"%s\n",body);        // write body to it.
  fclose(fp);             // close it.
  
  sprintf(cmd,"sendmail %s < %s",to, tmp); // prepare command.
  system(cmd);     // execute it.
}

int read_file(double *tt, double *ww, const int pp, const char *ff, const int nn) {
    int ii = 0;
    FILE *fp = fopen(ff, "r");
    if (fp == 0) {
      printf("ERROR: Failed to open file: %s\n", ff);
      return 1;
    }
    
    while ((ii != nn) && (2 == fscanf(fp, "%lf %lf", tt+ii, ww+ii))) ii++;

    if (((ii != nn) || 0 == fscanf(fp, "%*lf %*lf")) && (pp == ROOT)) {
      printf("ERROR: Provided number of events does not match number of events in file! Recommended to launch with PATH2FILE $(wc -l < PATH2FILE) as arguments. The ./run script should take care of this.\n");
      return 1;
    }
    
    return 0;
}

int parse_input(char *argv[], const int pp, char **ff, int *nn) {
  (*ff) = argv[1];
  (*nn) = (int) strtoumax(argv[2], NULL, 10);

  if (pp == ROOT) printf("ID=%d: Input data: %s\nID=%d: Number of lines to be read: %d\n", pp, (*ff), pp, (*nn));
  
  return 0;
}

double get_pdot_lim(const double ptmp) {
  if (SPIN_DOWN_LIM*pow(ptmp, 3) < AGE_LIM*ptmp) {
    return SPIN_DOWN_LIM*pow(ptmp, 3);
  } else {
    return AGE_LIM*ptmp;
  }
}

int get_plim(const int PP, const int pp, long *nt, double *pmin, double *pmax, double TT) {
  long np, n_tot = 0, ii, jj;
  double ptmp = PMIN, pdot, dpdot, pdot_lim;

  // One loop takes toughly 3 minutes, this makes it easier to compare to Python, so it's worth it.
  while (ptmp < PMAX) {
    pdot = 0.;
    pdot_lim = get_pdot_lim(ptmp);
    dpdot = 0.1*(ptmp*ptmp)/(TT*TT);

    while (pdot < pdot_lim) {
      if (pdot == 0.) n_tot++;
      else n_tot += 2;
      pdot += dpdot;
    }
    
    ptmp += SPARSE*ptmp*ptmp/(OVERLAP*TT);
  }

  // Step past periods of earlier processes
  ptmp = PMIN;
  jj=0;
  for (ii = 0; ii < pp; ii++) {
    GET_NP(n_tot, PP, ii);
    
    while (jj < np) {
      pdot_lim = get_pdot_lim(ptmp);
      dpdot = 0.1*(ptmp*ptmp)/(TT*TT);
      pdot = 0.;
      
      while (pdot < pdot_lim) {
	if (pdot == 0.) jj++;
	else jj += 2;
	pdot += dpdot;
      }
      // DBG      printf("ASDASDSADAS ID=%d: %ld %ld %ld\n", pp, jj, np, ii);
      ptmp += SPARSE*ptmp*ptmp/(OVERLAP*TT);
    }
    jj -= np;
  }
  
  // Assign local pmin and step to find local pmax
  (*pmin) = ptmp;
  GET_NP(n_tot, PP, pp);
  ii = 0;
  while (ii < np-jj) {
    pdot_lim = get_pdot_lim(ptmp);
    dpdot = 0.1*(ptmp*ptmp)/(TT*TT);
    pdot = 0.;
      
    while (pdot < pdot_lim) {
      if (pdot == 0.) ii++;
      else ii += 2;
      pdot += dpdot;
    }
    
    ptmp += SPARSE*ptmp*ptmp/(OVERLAP*TT);
  }

  (*nt) = ii;
  (*pmax) = ptmp;
  // DBG printf("ID=%d: %.17f %.17f %ld\n", pp, (*pmin), (*pmax), (*nt));
  return 0;
}

int get_par(const int PP, const int pp, const double *tt, const double *ww, const int nn, long *nt, double *pmin, double *pmax, double *k2, double *TT) {
  int ii;
  (*TT) = tt[nn-1] - tt[0];

  if (get_plim(PP, pp, nt, pmin, pmax, (*TT))) return 1;
  
  for (ii = 0; ii < nn; ii++) {
    (*k2) += ww[ii]*ww[ii];
  }
  (*k2) *= 0.5;

  //DBG  printf("ID=%d: %.17f %.17f %f %f\n", pp, (*pmin), (*pmax), (*k2), (*TT));
  //DBG  (*pmin) = 9.99999;
  //DBG  (*pmax) = 10.;
  return 0;
}

double get_phase(const double tt, const double ptmp, const double pdtmp) {
  double tp=tt/ptmp;
  return tp-0.5*pdtmp*tp*tp;
}

double get_Pn(const double *phase, const double *ww, const double k2, const double nn, const double mm) {
  int ii;
  double complex walk = 0.0;
  double complex coef = -2*M_PI*I*mm;

  for (ii = 0; ii < nn; ii++) walk += ww[ii]*cexp(coef*phase[ii]);
  walk = creal(walk)*creal(walk)+cimag(walk)*cimag(walk);
  return walk/k2;
}

int my_alloc(double **topP, double **topD, double **topH, const int PP, const int pp, const long nt) {
  long tmp, hlp;

  if (TOP) tmp = NTOP;
  else tmp = nt;

  // Allocates a buffer when writing to file.
  // char *buffer = malloc((nt * 73 + 1) * sizeof(char));
  hlp = PP*tmp*3*sizeof(double) + PP*tmp*73*sizeof(char);
  hlp = PP*tmp*3*sizeof(double);
  if (hlp > (long) MEMORY*1000000000) {
    printf("ERROR: Trying to allocate more than %d GB of RAM.\n", MEMORY);
    return 1;
  } else if (pp == 0) printf("ID=%d: Allocating %.3f GB of RAM.\n", pp, (double) hlp/1e9);
  
  (*topP) = calloc(tmp, sizeof(double));
  (*topD) = calloc(tmp, sizeof(double));
  (*topH) = calloc(tmp, sizeof(double));
  return 0;
}

double h_stat(const double *phase, const double *ww, const double k2, const double nn) {
  int ii, jj;
  double tmp, HH = 0.;
  double Z2[MMAX];
  
  for (ii = 0; ii < MMAX; ii++) {
    Z2[ii] = get_Pn(phase, ww, k2, nn, ii+1);
  }

  for (ii = 0; ii < MMAX; ii++) {
    tmp = 0.;
    for (jj = 0; jj <= ii; jj++) tmp += Z2[jj];
    tmp += - 4*(ii+1) + 4;
    if (tmp > HH) HH = tmp;
  }
  return HH;
}

void save2top(long *top_i, double *top_min, double *topP, double *topD, double *topH, const double ptmp, const double pdot, const double HH){
  int ii;
  
  if (HH > (*top_min)) {
    topP[(*top_i)] = ptmp;
    topD[(*top_i)] = pdot;
    topH[(*top_i)] = HH;
    (*top_min) = HH;
    
    for (ii = 0; ii < NTOP; ii++) {  
      if (topH[ii] < (*top_min)) {
	(*top_min) = topH[ii];
	(*top_i) = ii;
      }
    }
  }
}

int search_periods(const int pp, double *topP, double *topD, double *topH, int *dis, double *tt, const double *ww, const long nt, const int nn, const double pmin, const double k2, const double TT, char *ff) {
  long ii, jj, kk, npd, top_i = 0;
  double HH, start, *phase, pdot_lim, dpdot, ptmp = pmin, pdot, top_min = 0., top_max = 20.;
  phase = calloc(nn, sizeof(double));
  for (ii = 0; ii < nn; ii++) tt[ii] = tt[ii]-TREF;
  start = MPI_Wtime();

  ii = 0;
  while (ii < nt) {
    pdot_lim = get_pdot_lim(ptmp);
    pdot = 0.;
    npd = 0;
    dpdot = 0.1*(ptmp*ptmp)/(TT*TT);
    
    while (pdot < pdot_lim) {
      pdot += dpdot;
      npd++;
    }
    npd--;

    for (jj = -npd; jj <= npd; jj++) {
      pdot = jj*dpdot;
      for (kk = 0; kk < nn; kk++) phase[kk] = get_phase(tt[kk], ptmp, pdot);
      HH = h_stat(phase, ww, k2, nn);
      // DBG      HH = k2+ww[0];
      dis[(int) fmin(DIS_RES*HH, ND-1.)] += 1;

      if (TOP) save2top(&top_i, &top_min, topP, topD, topH, ptmp, pdot, HH);
      else {
	if (HH > top_max) {
	  //	  send_mail(ff, ptmp, pdot, HH);
	  (void)(ff); // suppress warning
	  //	  top_max = HH;
	  printf("%.17f %.17f %.17f\n", ptmp, pdot, HH);
	}
	topP[ii] = ptmp;
	topD[ii] = pdot;
	topH[ii] = HH;
      }
	
      // DBG printf("%d %.17f %.6e %.6e %.6e %.6f\n", pp, ptmp, pdot, dpdot, pdot_lim, HH);
      ii++;
      if ((pp == ROOT) && (ii%PRINT_RATE == 0) && (ii > 0)) printf("ID=%d: %8ld/%ld, %10.0f s, %10.0f s, %10.0f s\n", pp, ii, nt, MPI_Wtime()-start, (MPI_Wtime()-start)*(1/((double) ii/nt)-1), (MPI_Wtime()-start)/((double) ii/nt));
    }
    
    ptmp += SPARSE*ptmp*ptmp/(OVERLAP*TT);
  }

  free(phase);
  return 0;
}

int write2file(double *topP, double *topD, double *topH, const int *dis, const int pp, const long nt) {
  int ii;
  FILE *fp;
  char buf[71];

  MPI_File output;
  MPI_Status status;
  int sum[ND] = {0.0};
  MPI_Reduce(dis, sum, ND, MPI_INT, MPI_SUM, ROOT, MPI_COMM_WORLD);
  if (pp == ROOT) {
    sprintf(buf, "out/dis.txt");
    fp = fopen(buf, "w");
    for(ii = 0; ii < ND; ii++) fprintf(fp, "%d\n", sum[ii]);
    fclose(fp);
  }

  MPI_File_open(MPI_COMM_WORLD, "out/top.txt", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL,
                          &output);
  // delete any existing content
  MPI_File_set_size(output, 0);

  // fill buffers with formated data
  char *buffer = malloc((nt * 73 + 1) * sizeof(char));
  int count = 0;
  
  for (ii = 0; ii < nt; ii++) {
    sprintf(buffer + count, "%.17e %24.17e %.17e\n", topP[ii], topD[ii], topH[ii]);
    count += 73;
  }
  MPI_File_write_ordered(output, buffer, count, MPI_CHAR, &status);
  free(buffer);

  return 0;
}

int main(int argc, char *argv[]) {
  int PP, pp, nn, ii;
  long nt;
  char *ff = NULL;
  double pmin, pmax, k2 = 0.0, TT;
  double *topP, *topD, *topH, *tt, *ww;
  int dis[ND] = {0.0};
  

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &PP);
  MPI_Comm_rank(MPI_COMM_WORLD, &pp);
  INFO("Started", pp, 1);

  if (parse_input(argv, pp, &ff, &nn)) return 1;
  
  tt = calloc(nn, sizeof(double));
  ww = calloc(nn, sizeof(double));
  if (read_file(tt, ww, pp, ff, nn)) return 1;

  if (get_par(PP, pp, tt, ww, nn, &nt, &pmin, &pmax, &k2, &TT)) return 1;
  printf("ID=%d: pmin=%g, pmax=%g, nt=%ld, TT=%g\n", pp, pmin, pmax, nt, TT);
  
  if (my_alloc(&topP, &topD, &topH, PP, pp, nt)) return 1;

  //for (ii = 0; ii < nn; ii++) tt[ii] = 3328.44*rand()/RAND_MAX;
  (void)ii;
  
  if (search_periods(pp, topP, topD, topH, dis, tt, ww, nt, nn, pmin, k2, TT, ff)) return 1;

  if (write2file(topP, topD, topH, dis, pp, nt)) return 1;
  
  INFO("Cleaning", pp, 1);
  free(topP);
  free(topD);
  free(topH);
  if (pp == ROOT) INFO("All processes done. Program ending", pp, 1);
  return MPI_Finalize();
}
