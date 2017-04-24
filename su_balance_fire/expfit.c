/* expfit.c -- model functions for exponential + background */
/* mostly from GSL manual */

struct data {
    size_t n;
    double *t;
    double *y;
};

int expb_f(const gsl_vector *x, void *data, gsl_vector *f) {
    size_t n = ((struct data *)data)->n;
    double *t = ((struct data *)data)->t;
    double *y = ((struct data *)data)->y;

    double a = gsl_vector_get(x, 0);
    double b = gsl_vector_get(x, 1);

    size_t i;

    for (i=0; i<n; i++) {
        /* model: Yi = a * exp(b * t) */
        double Yi = a*exp(b*t[i]);
        gsl_vector_set(f, i, Yi-y[i]);
    }
    return GSL_SUCCESS;
}

int expb_df (const gsl_vector *x, void *data, gsl_matrix *J) {
    size_t n = ((struct data*)data)->n;
    double *t = ((struct data *)data)->t;

    double a = gsl_vector_get(x, 0);
    double b = gsl_vector_get(x, 1);

    size_t i;

    for (i=0; i<n; i++) {
        /* Jacobian matrix. */
        gsl_matrix_set(J, i, 0, exp(b*t[i]));
        gsl_matrix_set(J, i, 1, a*t[i]*exp(b*t[i]));
    }
    return GSL_SUCCESS;
}
