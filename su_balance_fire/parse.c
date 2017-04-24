#include <stdlib.h>
#include <stdio.h>


extern int verbose;

int parseargs(int argc, char* argv[], char* sustack, char* sumig, char* suoutfile) {
    int er1;
    char opti[82];
    int usage();

    argc--;
    er1=0;
    if (argc==0) {
        usage();
        fprintf(stderr,"\nNo arguments! Exit.\n");
        exit(1);
    }

    while (argc) {
        argv++;
        sscanf(*argv,"%s",opti);
        if (opti[0]=='-') {
            switch(opti[1]) {
                case 's':
                    if (argc) {
                        if ((sscanf(*++argv,"%80s",sustack))!=1) {
                            fprintf(stderr,"Error! Stack file invalid.\n");
                            usage();
                            exit(1);
                        }
                        argc--;
                        er1++;
                    }
                    break;
                case 'm':
                    if (argc) {
                        if ((sscanf(*++argv,"%80s",sumig))!=1) {
                            fprintf(stderr,"Error! Migrated file invalid.\n");
                            usage();
                            exit(1);
                        }
                        argc--;
                    }
                    break;
                case 'o':
                    if (argc) {
                        if ((sscanf(*++argv,"%80s",suoutfile))!=1) {
                            fprintf(stderr,"Error! Outfile name invalid.\n");
                            usage();
                            exit(1);
                        }
                        argc--;
                    }
                    break;
                case 'v':
                    verbose=1;
                    break;
            }
        }
        argc--;
    }
    if (er1!=1) {
        usage();
        exit(1);
    }
    return 1;
}

int usage() {
    printf("\nsu_balance_fire by S. Vakeva (2016-2017)\n\nUsage: su_balance_fire -s [STACK.su] -m [MIGRATED.su] -o [OUTFILE.su] [-v]\n");
    printf("\nScales the migrated SU file so that the section retains the same energy distribution in depth as the original stack.\n");
    printf("\nAverages all traces, determines a function of form a*exp(b*t) by which the migrated section is then scaled.\n");
    return 1;
}
