#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define Abs(a) ((a)>=0?(a):-(a))

static int total_point(FILE *fp) {
    int n=0;
    char line[100];
    while (fgets(line,100,fp)!=NULL) {
        n++;
    }
    rewind(fp);
    return n;
}

static double sum_x_pow(const double *x,int n,int p) {
    double s=0;
    for (int i=0; i<n; i++) {
        s+=pow(x[i],p*1.0);
    }
    return s;
}

static double sum_y_x_pow(const double *x,const double *y,const int n,const int p) {
    double s=0;
    for (int i=0; i<n; i++) {
        s+=pow(x[i],p*1.0)*y[i];
    }
    return s;
}

static void lu_decomposition_without_file(double **matrix, double *soln, const int n) {
    for (int step=0; step<n-1; step++) {
        double col_max=Abs(matrix[step][step]);

        for (int i=step+1; i<n; i++) {
            if (Abs(matrix[i][step])>Abs(col_max)) {
                col_max=Abs(matrix[i][step]);
            }
        }

        if (col_max<1e-10) {
            printf("Singular Matrix. No unique Solution!\n");
            return;
        }

        if (Abs(matrix[step][step])!=col_max) {
            for (int i=step+1; i<n; i++) {
                if (col_max==Abs(matrix[i][step])) {
                    double *temp= malloc(sizeof(double)*(n+1));
                    if (!temp) {
                        return;
                    }

                    for (int it=0; it<n+1; it++) {
                        temp[it]=matrix[i][it];
                    }
                    for (int it=0; it<n+1; it++) {
                        matrix[i][it]=matrix[step][it];
                    }
                    for (int it=0; it<n+1; it++) {
                        matrix[step][it]=temp[it];
                    }
                    // pivoting done. now gaussian kora need
                    free(temp);
                    break;
                }
            }
        }

        // gaussian elimination here
        for (int i=step+1; i<n; i++) {
            double factor= matrix[i][step]/matrix[step][step];
            for(int j=step; j<n; j++) {
                if (j==step) {
                    matrix[i][j]= factor;
                }
                else {
                    matrix[i][j]+= -factor*matrix[step][j];
                }
            }
        }
    }

    if (Abs(matrix[n-1][n-1])<1e-10) {
        printf("Singular Matrix. No unique Solution!\n");
        return;
    }

    double *y= malloc(sizeof(double)*n);
    if (!y) {
        return;
    }
    for (int i=0; i<n; i++) {
        y[i]=0;
    }

    for (int i=0; i<n; i++) {
        double s=0;
        for (int j=0; j<i; j++) {
            s+=matrix[i][j]*y[j];
        }
        y[i]=matrix[i][n]-s;
    }

    for (int i=0; i<n; i++) {
        soln[i]=0;
    }

    for (int i=n-1; i>=0; i--) {
        double s=0;
        for (int j=n-1; j>=i+1; j--) {
            s+=matrix[i][j]*soln[j];
        }
        soln[i]=(y[i]-s)/matrix[i][i];
    }

    free(y);
}

static void polyfit(FILE *fp,double *soln, const int fit_rank) {
    const int n= total_point(fp);
    double **matrix= malloc(sizeof(double*)*(fit_rank+1));
    if (!matrix) {
        return;
    }
    for (int i=0; i<=fit_rank; i++) {
        matrix[i]=malloc(sizeof(double)*(fit_rank+2));
        if (!matrix[i]) {
            free(matrix);
            return;
        }
    }

    double *x,*y;
    x=malloc(sizeof(double)*n);
    if (!x) {
        for (int i=0; i<=fit_rank; i++) {
            free(matrix[i]);
        }
        free(matrix);
        return;
    }
    y=malloc(sizeof(double)*n);
    if (!y) {
        for (int i=0; i<=fit_rank; i++) {
            free(matrix[i]);
        }
        free(matrix);
        free(x);
        return;
    }

    for (int i=0; i<n; i++) {
        fscanf(fp,"%lf\t%lf",&x[i],&y[i]);
    }
    rewind(fp);

    for (int i=0; i<=fit_rank; i++) {
        for (int j=0; j<=fit_rank; j++) {
            matrix[i][j]=sum_x_pow(x,n,i+j);
        }
    }

    for (int i=0; i<=fit_rank; i++) {
        matrix[i][fit_rank+1]=sum_y_x_pow(x,y,n,i);
    }

    lu_decomposition_without_file(matrix,soln,fit_rank+1);

    // freeing here
    for (int i=0; i<=fit_rank; i++) {
        free(matrix[i]);
    }
    free(matrix);
    free(x);
    free(y);
}

static void polyfit_without_file(const double *x, const double *y, double *soln,
    int n, int fit_rank) {
    double **matrix= malloc(sizeof(double*)*(fit_rank+1));
    if (!matrix) {
        return;
    }
    for (int i=0; i<=fit_rank; i++) {
        matrix[i]=malloc(sizeof(double)*(fit_rank+2));
        if (!matrix[i]) {
            free(matrix);
            return;
        }
    }

    for (int i=0; i<=fit_rank; i++) {
        for (int j=0; j<=fit_rank; j++) {
            matrix[i][j]=sum_x_pow(x,n,i+j);
        }
    }

    for (int i=0; i<=fit_rank; i++) {
        matrix[i][fit_rank+1]=sum_y_x_pow(x,y,n,i);
    }

    lu_decomposition_without_file(matrix,soln,fit_rank+1);

    // freeing here
    for (int i=0; i<=fit_rank; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

static double goodness(FILE *fp, const double *co_eff, const int fit_rank) {
    int n= total_point(fp);

    double*x=malloc(sizeof(double)*n);
    if (!x) {
        return NAN;
    }
    double *y=malloc(sizeof(double)*n);
    if (!y) {
        free(x);
        return NAN;
    }

    for (int i=0; i<n; i++) {
        fscanf(fp,"%lf\t%lf",&x[i],&y[i]);
    }
    rewind(fp);

    double avg_y= sum_y_x_pow(x,y,n,0)/n;
    // printf("%.4lf\n",avg_y);

    double st=0;
    for (int i=0; i<n; i++) {
        st+=(avg_y-y[i])*(avg_y-y[i]);
    }

    double *fitted_y= malloc(sizeof(double)*n);
    if (!fitted_y) {
        free(x);
        free(y);
        return NAN;
    }

    for (int i=0; i<n; i++) {
        double s=0;
        for (int j=0; j<fit_rank+1; j++) {
            s+=co_eff[j]*pow(x[i],j*1.0);
        }
        fitted_y[i]=s;
    }

    double sr=0;
    for (int i=0; i<n; i++) {
        sr+=(fitted_y[i]-y[i])*(fitted_y[i]-y[i]);
    }

    free(x);
    free(y);
    free(fitted_y);
    return (1-sr/st);
}

static double goodness_without_file(const double *x, const double *y,const double *co_eff,
    const int n, const int fit_rank) {
    double avg_y= sum_y_x_pow(x,y,n,0)/n;
    double st=0;
    for (int i=0; i<n; i++) {
        st+=(avg_y-y[i])*(avg_y-y[i]);
    }

    double *fitted_y= malloc(sizeof(double)*n);
    if (!fitted_y) {
        return NAN;
    }

    for (int i=0; i<n; i++) {
        double s=0;
        for (int j=0; j<fit_rank+1; j++) {
            s+=co_eff[j]*pow(x[i],j*1.0);
        }
        fitted_y[i]=s;
    }

    double sr=0;
    for (int i=0; i<n; i++) {
        sr+=(fitted_y[i]-y[i])*(fitted_y[i]-y[i]);
    }

    free(fitted_y);
    return (1-sr/st);
}

static int best_fit(FILE *fp, const int high) {
    double good=0;
    int rank=0;
    for (int i=1; i<=high; i++) {
        double *co_eff= malloc(sizeof(double)*(i+1));
        if (!co_eff) {
            return NAN;
        }
        polyfit(fp,co_eff,i);
        double gd=goodness(fp,co_eff,i);

        if (gd>good) {
            good=gd;
            rank=i;
        }
        free(co_eff);
    }
    return rank;
}

static int best_fit_without_file(const double *x, const double *y,
    const int n, const int high) {
    double good=0;
    int rank=0;
    for (int i=1; i<=high; i++) {
        double *co_eff= malloc(sizeof(double)*(i+1));
        if (!co_eff) {
            return NAN;
        }
        polyfit_without_file(x,y,co_eff,n,i);
        const double gd=goodness_without_file(x,y,co_eff,n,i);

        if (gd>good) {
            good=gd;
            rank=i;
        }
        free(co_eff);
    }
    return rank;
}


int main() {
    FILE *fp=fopen("curve_fitting_data.txt","r");
    if (!fp) {
        return 1;
    }
    int total= total_point(fp);

    // int n=5;
    // double *co_eff= malloc(sizeof(double)*(n+1));
    // if (!co_eff) {
    //     fclose(fp);
    //     return 1;
    // }
    //
    // polyfit(fp,co_eff,total,n);
    // for (int i=0; i<=n; i++) {
    //     printf("%.10lf\n",co_eff[i]);
    // }
    // printf("%.6lf\n",goodness(fp,co_eff,n));
    // printf("%d\n%.6lf\n",rank,good);

    printf("%d\n",best_fit(fp,5));

    fclose(fp);
    // free(co_eff);
}
