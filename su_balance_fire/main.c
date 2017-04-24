#include <stdio.h>
#include <stdlib.h>
#include <su.h>
#include <segy.h>
#include <float.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

#include "expfit.c"

#define MAXLEN 80

int parseargs(int argc, char* argv[], char* sustack, char* sumig, char* suoutfile);
int verbose;

char **sdoc;

int main(int argc, char** argv) {
    char sustack[MAXLEN], sumig[MAXLEN], suoutfile[MAXLEN];
    FILE *fpstack, *fpmig, *fpout;
    FILE *gnuplotPipe;
    FILE *fptextdump;

    segy tr;
    int nt;
    float dt;
    float *st_sum, *mig_sum;
    double trener;
    int i, nfit;
    int nooftraces;

    float t_min, t_max;
    char c;

    int query;

    /* solver things, from GSL manual */
    double *y, *t;
    const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
    gsl_multifit_fdfsolver *s;
    gsl_vector *res_f;
    double chi, chi0;

    int status, info;
    gsl_multifit_function_fdf f;

    double x_init[] = { 0.61, -0.04 };
    gsl_vector_view x = gsl_vector_view_array(x_init, 2);

    parseargs(argc, argv, sustack, sumig, suoutfile);
    gnuplotPipe = (FILE*) popen("gnuplot -persistent", "w");

    if ((fpstack=fopen(sustack,"r"))==NULL) {
        fprintf(stderr, "Cannot open stack file.\n");
        exit(1);
    }
    if (verbose) printf("File %s opened.\n",sustack);


    if (!fgettr(fpstack,&tr)) {
        fprintf(stderr,"Can't get first trace.\n");
        exit(1);
    }

    nt = (int) tr.ns;
    dt = ((double) tr.dt)/1000000.0;

    if (verbose) printf("nt = %d, dt = %f\n",nt,dt);

    if ((st_sum=malloc(sizeof(float)*nt))==NULL) {
        fprintf(stderr,"Cannot allocate memory for summed stack.\n");
        exit(1);
    }

    /* Read in stack. */
    nooftraces=0;
    do {
        trener=0.0;
        for (i=0; i<nt; i++) {
            trener+=pow(tr.data[i],2.0);
        }
        trener/=nt;
        if (trener<0.5) {
            if (verbose) printf("Excluding cdp = %d, average energy: %f\n",tr.cdp,trener);
        } else {
            nooftraces++;
            for (i=0; i<nt; i++) {
                st_sum[i]+=fabs(tr.data[i]);
                /* IF RMS: st_sum[i]+=pow(tr.data[i],2.0); */
            }
        }

    } while (fgettr(fpstack, &tr));

    if (verbose) printf("Stack has %d traces.\n",nooftraces);
    /* Find average energy function for stack, to be mimicked in the output file, and plot. */

    for (i=0; i<nt; i++) {
        st_sum[i]/=nooftraces;
        printf("st_sum[%d] = %f\n",i,st_sum[i]);
    }

    fclose(fpstack);

    /* QUERY WHAT THE USER WANTS TO DO */
    /* REMOVE FROM HERE if you want to automate the action */
    printf("What would you like to do?\n");
    printf("   1 = view amplitude decay, 2 = balance section\n");
    if (scanf("%d",&query)<0) {
        fprintf(stderr,"Error.\n");
        exit(1);
    }
    /* REMOVE UNTIL HERE, but remember to set query ! */

    if (query==2) {

    /* Read in a migrated SU file. */
    if ((fpmig=fopen(sumig,"r"))==NULL) {
        fprintf(stderr, "Cannot open migrated file.\n");
    }
    if (verbose) printf("File %s opened.\n",sumig);


    if (!fgettr(fpmig,&tr)) {
        fprintf(stderr,"Can't get first trace.\n");
        exit(1);
    }

    nt = tr.ns < nt ? tr.ns : nt;
    dt = ((double) tr.dt)/1000000.0;

    if (verbose) printf("nt = %d, dt = %f\n",nt,dt);

    if ((mig_sum=malloc(sizeof(float)*nt))==NULL) {
        fprintf(stderr,"Cannot allocate memory for migrated stack.\n");
        exit(1);
    }

    nooftraces=0;
    do {
        trener=0.0;
        for (i=0; i<nt; i++) {
            trener+=pow(tr.data[i],2.0);
        }
        trener/=nt;
        if (isnan(trener)) {
            if (verbose) printf("Excluding cdp = %d, average energy: %f\n",tr.cdp,trener);
        } else {
            nooftraces++;
            for (i=0; i<nt; i++) {
                mig_sum[i]+=fabs(tr.data[i]);
                /* if RMS: mig_sum[i]+=pow(tr.data[i],2.0); */
            }
        }

    } while (fgettr(fpmig, &tr));

    if (verbose) printf("Migrated stack has %d traces.\n",nooftraces);
    for (i=0; i<nt; i++) {
        printf("mig_sum[i] = %f\n",mig_sum[i]);
        mig_sum[i]/=nooftraces;
    }

    fclose(fpmig);

    /* end of if (query==2) */
    }

    /* collect amplitude decay data in this datatype */
        if ((y=malloc(sizeof(double)*nt))==NULL) {
        fprintf(stderr,"Cannot allocate space for y.\n");
    }
        if ((t=malloc(sizeof(double)*nt))==NULL) {
        fprintf(stderr,"Cannot allocate space for t.\n");
    }

    nfit=0;
    /* Plot amplitude decay function and initialise data for fitting. */

    /* IF THE USER ONLY WANTS TO VIEW AN AMPLITUDE DECAY CURVE */
    if (query==1) {
    if ((fptextdump=fopen("dump.txt","w"))==NULL) {
        fprintf(stderr, "Cannot open dump file. Exit.\n");
        exit(1);
    }

    fprintf(gnuplotPipe, "set multiplot\n");
    fprintf(gnuplotPipe, "plot '-' \n");
        for (i=0; i<nt; i++) {
            if ((float)i*dt<29.0)
            fprintf(gnuplotPipe, "%f %f\n", (float)i*dt, st_sum[i]);
            fprintf(fptextdump, "%f %f\n", (float)i*dt, st_sum[i]);
        }
        fprintf(gnuplotPipe, "e\n");

    fclose(fptextdump);
    }

    else if (query==2) {
    fprintf(gnuplotPipe, "set multiplot\n");
    fprintf(gnuplotPipe, "plot '-','-' with lines \n");


    printf("Enter t_min t_max.\n");
    if (scanf("%f %f",&t_min,&t_max)<0) {
        fprintf(stderr,"Error!\n");
        exit(1);
    }

    /* Suitable for OpenFIRE NMO stacks:
    t_min=4.0;
    t_max=22.0;
    */

    /* Suitable for OpenFIRE DMO stacks:
    t_min=2.0;
    t_max=12.0;
    */

    for (i=4; i<0.8*nt; i++) {
        if (mig_sum[i]!=0.0) {
            printf("st_sum = %f, mig_sum = %f\n",st_sum[i],mig_sum[i]);
            fprintf(gnuplotPipe, "%f %f\n", (float)i*dt, st_sum[i]/mig_sum[i]);
            if ((float)i*dt <= t_max && (float)i*dt >= t_min) {
                t[nfit]=(double)i*dt;
                y[nfit]=(double)(st_sum[i]/mig_sum[i]);
                nfit+=1;
            }
        }
    }
    fprintf(gnuplotPipe, "e\n");



    /* initialize solver function f and Jacobian */
    /* closely follows the GSL manual */
    struct data d = { nfit, t, y };
    f.f = &expb_f;
    f.df = &expb_df;
    f.n = nfit;
    f.p = 2;
    f.params = &d;

    /* initialize solver */
    s = gsl_multifit_fdfsolver_alloc(T, nfit, 2);
    gsl_multifit_fdfsolver_set(s, &f, &x.vector);

    /* compute initial residual norm */
    res_f = gsl_multifit_fdfsolver_residual(s);
    chi0=gsl_blas_dnrm2(res_f);

        /* solve the system */
        status = gsl_multifit_fdfsolver_driver(s, 20, 1.0e-8, 1.0e-8, 0.0, &info);

    /* error limits */
    gsl_matrix *J = gsl_matrix_alloc(nfit, 2);
    gsl_matrix *covar = gsl_matrix_alloc(2,2);
    gsl_multifit_fdfsolver_jac(s, J);
    gsl_multifit_covar (J, 0.0, covar);

    chi=gsl_blas_dnrm2(res_f);

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

    fprintf(stderr, "Number of iterations: %zu\n", gsl_multifit_fdfsolver_niter(s));
    fprintf(stderr, "Function evaluations: %zu\n", f.nevalf);
    fprintf(stderr, "Jacobian evaluations: %zu\n", f.nevaldf);
    fprintf(stderr, "Reason for stopping: %s\n", (info == 1) ? "small step size" : "small gradient");
    fprintf(stderr, "initial |f(x)| = %g\n", chi0);
    fprintf(stderr, "final   |f(x)| = %g\n", chi);
    fprintf(stderr,"chisq/dof = %g\n",pow(chi,2.0)/(nfit-3));
    fprintf(stderr," model: f(x) = a*exp(b*x)\n");
        fprintf(stderr," a = %.5f +/- %.5f,\n b = %.5f +/- %.5f\n",FIT(0),ERR(0),FIT(1),ERR(1));
    fprintf(stderr, "Exit status = %s\n",gsl_strerror(status));


    for (i=4; i<0.8*nt; i++) {
        fprintf(gnuplotPipe, "%f %f\n", (float)i*dt, FIT(0)*exp(FIT(1)*i*dt));
    }
    fprintf(gnuplotPipe,"e\n");
    free(y);
    free(t);
    free(st_sum);
    free(mig_sum);

    /* ready to start writing new migrated file */

    /* REMOVE FROM HERE if you want to skip confirmation */
    printf("Happy with fit? Continue? (y/N)\n");
    while ((c=getchar()) != '\n')
    if (c != 'y') {
        printf("You typed no.\n");
        exit(1);
    }
    /* REMOVE UNTIL HERE */

    /* Read in the migrated SU once again, this time for output too. */
    if ((fpmig=fopen(sumig,"r"))==NULL) {
        fprintf(stderr, "Cannot open migrated file.\n");
    }
    if (verbose) printf("File %s opened.\n",sumig);

    if ((fpout=fopen(suoutfile,"w"))==NULL) {
        fprintf(stderr, "Cannot open outfile.\n");
    }
    if (verbose) printf("File %s opened.\n",suoutfile);

    if (!fgettr(fpmig, &tr)) {
        fprintf(stderr, "Cannot get trace from migrated file. Exit.\n");
        exit(1);
    }

    nt = (int) tr.ns;
    dt = ((double) tr.dt)/1000000.0;

    do {
        for (i=0; i<nt; i++) {
            tr.data[i]=tr.data[i]*(FIT(0)*exp(FIT(1)*i*dt));
        }
        fputtr(fpout, &tr);
    } while (fgettr(fpmig, &tr));

    gsl_multifit_fdfsolver_free(s);
    gsl_matrix_free(covar);
    gsl_matrix_free(J);

    fclose(fpmig);
    fclose(fpout);

    /* end of query==2 */
    }
    fprintf(gnuplotPipe, "set xlabel \"TWT (seconds)\"\n");
    fprintf(gnuplotPipe, "set ylabel \"mean amplitude of stack\"\n");
    fprintf(gnuplotPipe, "set title \"Stack - amplitude decay\"\n");
    fflush(gnuplotPipe);
    return 0;
}
