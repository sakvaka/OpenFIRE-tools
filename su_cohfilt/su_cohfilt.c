#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <su.h>
#include <segy.h>
#include <cwp.h>

#define LENMAX 16384
#define PANMAX 120
#define PTMAX 500
/* 
 *     the command is run with
 *
 *     su_cohfilt <input.su >output.su 
 * 
 *     run with unlimited stack size!
 */
int smooth(int iwin, int nleng, float* ros);
int semb(float* ros, float* rnim, float* rcoh, float rnorm, int nleng);

char **sdoc;

int main(void) {
    int i,j;
    segy tr, tr2;
    float ros[LENMAX],rnim[LENMAX];
    float ramppi,rcoh[LENMAX];
    float trace2[LENMAX],trace[LENMAX],aputra[LENMAX];
    float migtra[LENMAX][PANMAX];
    int cdpind[PANMAX],trsind[PANMAX];
    float distar[PANMAX];
    float bell[PANMAX],ptbell[PTMAX];
    int final;
    float srate;

    float staspace,slomin,slomax,slodel;
    int nwid,ibell,iwindow,iptbell;
    float kwin,apu,pii;
    float rnorm;
    int nofpt,ipanel,iover,itrace,nmid;
    int nleng,ndatlen,iii,kloop,iviive,jmin,jmax,jk;
    int itra;
    float paino,weiwei,dista,xtee,amp;
    float slow;

    final=0;

    /* set: staspace, slomin, slomax, slodel, nwid, ibell, iwindow, iptbell */
    /* TODO: lue tiedostosta */
    staspace=25.0;
    slomin=-0.8;
    slomax=0.8;
    slodel=0.025;
    nwid=31;
    ibell=5;
    iwindow=25;
    iptbell=5;

    nmid=(nwid-1.0)/2.0;
    kwin=(iwindow-1.0)/2.0;

    if(slomax<slomin) {
        apu=slomax;
        slomax=slomin;
        slomin=apu;
    }

    nofpt=(int)((slomax-slomin)/slodel)+1;

    fprintf(stderr,"slomin=%f\n",slomin);
    fprintf(stderr,"slomax=%f\n",slomax);
    fprintf(stderr,"nofpt=%d\n",nofpt);
    pii=4.0*atan(1.0);

    for (i=0;i<PANMAX;i++) {
        for (j=0;j<LENMAX;j++) {
            migtra[j][i]=0.0;
        }
    }

    for (j=0;j<PANMAX;j++) {
        bell[j]=1.0;
    }

    for (j=0;j<PANMAX;j++) { 
        bell[j]=1.0;
    }

    for (j=0;j<ibell;j++) {
        ramppi=(1.0-cos(pii*(float)(j+1.0)/(float)ibell))/2.0;
        bell[j]=ramppi;
        bell[nwid-1-j]=ramppi;
    }

    if(nwid<PANMAX) {
        for (j=nwid;j<PANMAX;j++) {
            bell[j]=0.0;
        }
    }

    for (j=0;j<nwid;j++) {
        fprintf(stderr,"j %d %f\n",j,bell[j]);
    }

    for (j=0;j<PTMAX;j++) {
        ptbell[j]=1.0;
    }

    for (j=0;j<iptbell;j++) {
        ramppi=(1.0-cos(pii*(float)(j+1.0)/(float)iptbell))/2.0;
        ptbell[j]=ramppi;
        ptbell[nofpt-1-j]=ramppi;
    }

    if(nofpt<PTMAX) {
        for (j=nofpt;j<PTMAX;j++) {
            ptbell[j]=0.0;
        }
    }

    fprintf(stderr,"-------------------------\n");
    for (j=0;j<nofpt;j++) {
        fprintf(stderr,"j %d %f\n",j,ptbell[j]);
    }

    fprintf(stderr,"-------------------------\n");
    for (ipanel=0;ipanel<nwid;ipanel++) {
        distar[ipanel]=(float)(ipanel-nmid)*staspace;
    }
    rnorm=1.0/(float)nwid; 

    iover=0;

    for (itrace=0;itrace<nmid;itrace++) {

        /* read trace */
        if (!gettr(&tr)) {
            final=1;
            goto loop1;
        }

        ipanel=nmid+itrace;    

        for (j=0;j<tr.ns;j++) {
            migtra[j][ipanel]=tr.data[j];
        }

        srate=(float)tr.dt/1000.0;
        cdpind[ipanel]=tr.cdp;
        trsind[ipanel]=tr.tracl;

    }

loop2: 

    /* read trace */
    if (!gettr(&tr)) {
        final=1;
        goto loop1;
    }

    ipanel=nwid-1;

    for (j=0;j<tr.ns;j++) {
        migtra[j][ipanel]=tr.data[j];
    }

    srate=(float)tr.dt/1000.0;
    cdpind[ipanel]=tr.cdp;
    trsind[ipanel]=tr.tracl;

loop1:

    /*----  start pt */                      

    for (j=0;j<LENMAX;j++) {
        trace2[j]=0.0;
    }

    nleng=tr.ns;
    if (nleng>LENMAX) {
        nleng=LENMAX;
    } 
    ndatlen=tr.ns;

    for (iii=0;iii<nofpt;iii++) {

        paino=ptbell[iii];
        slow=slomin+(float)(iii)*slodel;

        for (j=0;j<LENMAX;j++) {
            aputra[j]=0.0;
            trace[j]=0.0;
        }

        for (kloop=0;kloop<nwid;kloop++) {
            weiwei=bell[kloop]*rnorm;
            dista=distar[kloop];
            xtee=(dista*slow)/srate;   
            iviive=(int)xtee;

            jmin = (iviive<0) ? (-iviive) : 0;
            jmax = ((ndatlen-iviive)>nleng) ? (nleng-1) : (ndatlen-1-iviive);

            for (j=jmin;j<=jmax;j++) {
                amp=migtra[j+iviive][kloop]*weiwei;
                trace[j]=amp+trace[j];
                aputra[j]=aputra[j]+amp*amp;
            }
        }

        for (jk=0;jk<nleng;jk++) {
            ros[jk]=pow(trace[jk],2.0);
            rnim[jk]=aputra[jk];
            rcoh[jk]=1.0;
        }

        /*  calculate semblance of the traces ! */

        if(iwindow>0) {
            smooth(kwin,nleng,ros);
            smooth(kwin,nleng,rnim);
            semb(ros,rnim,rcoh,rnorm,nleng);  
        }
        /*
           add the calculated pt trace to earlier ones
           */
        for (jk=0;jk<nleng;jk++) {
            trace2[jk]=trace2[jk]+trace[jk]*rcoh[jk]*paino;
        }

    }
    /*
       all p values are now gone through and sum is calculated
       save mid-trace
       */
    tr2.cdp=cdpind[nmid];
    tr2.tracl=trsind[nmid];
    tr2.ns=tr.ns;
    tr2.dt=tr.dt;
    if((tr2.cdp%100)==0) fprintf(stderr,"cdpnum=%d\n",tr2.cdp);
    if ((tr2.cdp==0)) {
        fprintf(stderr, "Error! CDP = 0\n");
        for (i=nmid-5;i<=nmid+5;i++) {
            fprintf(stderr,"i=%d cdp=%d\n",i,cdpind[i]);
        }
    }

    for (i=0;i<tr.ns;i++) {
        tr2.data[i]=trace2[i];
    }
    puttr(&tr2);

    /* update tables */

    for (itra=0;itra<nwid-1;itra++) {
        cdpind[itra]=cdpind[itra+1];
        trsind[itra]=trsind[itra+1];
        for (jk=0;jk<tr.ns;jk++) {
            migtra[jk][itra]=migtra[jk][itra+1];
        }
    }
    for (jk=0;jk<LENMAX;jk++) {
        migtra[jk][nwid]=0.0;
    }

    /* go read a new trace */

    if(final) {
        iover=iover+1;
        if(iover<(nmid+1)) goto loop1;
    } else { 
        goto loop2;
    }

    return 0;

} 
int smooth(int iwin, int nleng, float* ros) {
    float aputra[24000];
    int j,jk,kk;
    float smo;
    int ialku,ilopp;
    for (j=0;j<nleng;j++) { 
        smo=0.0;
        ialku = (j>iwin) ? (j-iwin) : 0;
        ilopp = ((nleng-1)<(j+iwin)) ? nleng-1 : j+iwin;

        for (kk=ialku;kk<=ilopp;kk++) {
            smo=smo+ros[kk];
        }
        aputra[j]=smo;
    }

    for (jk=0;jk<nleng;jk++) {
        ros[jk]=aputra[jk];
    }
    return 0;
}

int semb(float* ros, float* rnim, float* rcoh, float rnorm, int nleng) {
    float semmes;
    int jk;
    for (jk=0;jk<nleng;jk++) {
        if (rnim[jk]>0.0) {
            semmes=ros[jk]/rnim[jk];
        } else {
            semmes=0.0;
        }
        rcoh[jk]=semmes*rnorm;
    }
    return 0;
}
