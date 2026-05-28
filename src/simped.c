/****************************************************************************************
SimPed: A simulation program to generate haplotype and genotype data for pedigree structures.

The SimPed program generates haplotype and/or genotype data for pedigrees as follows.  All of
the founders within the pedigree are assigned haplotypes and/or alleles conditional on the
user specified frequencies. Once assignment is completed each founder has two haplotypes.
Starting at the top of the pedigree structure the first offspring of the founder is randomly
assigned one of the founder's haplotypes. The allele at the first marker from this haplotype
is assigned to the offspring.  It is then determined based upon the genetic map whether a
recombination event has occurred between the first and second marker loci. If with probability,
a recombination event has occurred then at the second marker locus the allele from the founder's
other haplotype is assigned to the offspring.  If recombination event has not occurred with
probability (1- ) then the allele from the founder's same haplotype is assigned at the offspring's
second marker locus.  This procedure is repeated until alleles for all markers loci have been
assigned from one founder to their offspring.  The process is then repeated this time assigning
alleles to the offspring from their other parent.  In this manner the haplotypes flow down the
pedigree tree as all nonfounders are assigned haplotypes conditional on parental haplotypes.
Once all individuals within the pedigree have been assigned haplotypes, for those individuals/
marker loci for which it was specified that they are unavailable the genotypes are made unknown
(i.e. 0 0).

Original authors: Xiaoming Leal & Suzanne M. Leal (2005).
This file is a modernized build of the original simped.c so it compiles cleanly with
contemporary C compilers (clang / gcc -std=c99+).  Algorithmic behavior is unchanged.
********************************************************************************************/

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#define MAX_PED         100     /* maximum number of pedigrees                 */
#define MAX_IND         3000    /* maximum number of individuals               */
#define MAX_LOCUS       10000   /* maximum number of loci                      */
#define MAX_HAPLOTYPES  1000    /* maximum number of haplotypes                */
#define MAX_ALLELE      20      /* maximum number of alleles at a locus        */

struct genotype {
    int geno_chars[MAX_LOCUS][2];
};

typedef struct person {
    char pedId[20];
    char id[20];
    char father[20];
    char mother[20];
    int  sex;
    int  filled;
    struct genotype *geno;
} Person;

static Person *indi[MAX_IND];
static int    membersOfFamily[MAX_PED];
static int    seeds[3];
static float  recomb_fractions[MAX_LOCUS];
static int    totalInd, totalPed, totalLoci;
static int    extra_col, replicates, inter_dist_type;
static int    haplotype[MAX_HAPLOTYPES][MAX_LOCUS];
static float  possibilities[MAX_HAPLOTYPES], allele_freq[MAX_ALLELE];
static long   offset;            /* current file position in the parameter file */

/* Renamed from `index` to avoid a clash with POSIX's index() declared in
   <strings.h> (transitively pulled in via <string.h> on macOS / BSD). */
static int    founder_marker_idx; /* reset in main() each replicate; advanced by
                                     buildFounderGeno() to track how many marker
                                     columns have already been filled for founders. */
static int    familyNo;           /* re-numbers pedigrees across replicates       */
static int    iteration;          /* 2nd value on line 5 of the parameter file    */

/* ---- forward declarations ---- */
static float randomnum(void);
static void  buildFounderGeno(int type, int loci, int haplo_or_allele);
static void  assignGenotype(int thisPerson, int parent, char c);
static void  readPed(const char *pedfile);
static void  readParam(const char *paramfile);
static void  readMarkerInfo(const char *paramfile);
static void  simulate(void);
static void  writeResult(const char *pedfile, FILE *fp2);


/************************************************/
/*                                              */
/*              randomnum                       */
/*                                              */
/************************************************/
static float randomnum(void)  /* Wichmann-Hill style RNG */
{
    /* Global variables used: seeds[3] */
    float r;
    seeds[0] = seeds[0] % 177 * 171 - seeds[0] / 177 * 2;
    seeds[1] = seeds[1] % 176 * 172 - seeds[1] / 176 * 35;
    seeds[2] = seeds[2] % 178 * 170 - seeds[2] / 178 * 63;
    if (seeds[0] < 0) seeds[0] += 30269;
    if (seeds[1] < 0) seeds[1] += 30307;
    if (seeds[2] < 0) seeds[2] += 30323;
    r = seeds[0] / 30269.00f + seeds[1] / 30307.00f + seeds[2] / 30323.00f;

    return (float)(r - (long)r);
}

/************************************************/
/*                                              */
/*              buildFounderGeno                */
/*                                              */
/************************************************/
/* called by readMarkerInfo(); fills founder genotypes using the RNG. */
static void buildFounderGeno(int type, int loci, int haplo_or_allele)
{
    float randnum;
    int   i, j, t, tmp;
    int   count;

    tmp = 0;
    for (i = 0; i < totalInd; i++)
    {
        if (strcmp(indi[i]->father, "0") == 0 &&
            strcmp(indi[i]->mother, "0") == 0)
        {
            count = 0;
            while (count < 2)
            {
                if (type == 1)   /* haplotype information */
                {
                    randnum = randomnum();
                    for (j = 0; j < haplo_or_allele; j++)
                        if (randnum < possibilities[j]) { tmp = j; break; }
                    for (t = 0; t < loci; t++)
                        indi[i]->geno->geno_chars[founder_marker_idx + t][count] = haplotype[tmp][t];
                }
                if (type == 2)   /* genotype information */
                {
                    randnum = randomnum();
                    for (j = 0; j < haplo_or_allele; j++)
                        if (randnum < allele_freq[j]) { tmp = j + 1; break; }
                    indi[i]->geno->geno_chars[founder_marker_idx][count] = tmp;
                }
                count++;
            }
            indi[i]->filled = 1;
        }
    } /* end of assigning genotypes to founders */
    founder_marker_idx += loci;
}

/************************************************/
/*                                              */
/*              assignGenotype                  */
/*                                              */
/************************************************/
static void assignGenotype(int thisPerson, int parent, char c)
{
    int j, k = 0, side;
    float randnum;

    if (c == 'f') k = 0;          /* from father */
    if (c == 'm') k = 1;          /* from mother */

    /* 1st locus: pick which parental haplotype to start from */
    randnum = randomnum();
    if (randnum < 0.5f) {
        indi[thisPerson]->geno->geno_chars[0][k] = indi[parent]->geno->geno_chars[0][0];
        side = 0;
    } else {
        indi[thisPerson]->geno->geno_chars[0][k] = indi[parent]->geno->geno_chars[0][1];
        side = 1;
    }

    /* remaining loci: switch sides on recombination */
    for (j = 1; j < totalLoci; j++)
    {
        randnum = randomnum();
        if (randnum <= recomb_fractions[j - 1])
            side = 1 - side;
        indi[thisPerson]->geno->geno_chars[j][k] = indi[parent]->geno->geno_chars[j][side];
    }
}

/************************************************/
/*                                              */
/*              readPed                         */
/*                                              */
/************************************************/
static void readPed(const char *pedfile)
{
    int  i;
    char currentFamily[20];
    FILE *fp;

    i = 0;
    totalPed = 0;
    indi[i] = (Person *) malloc(sizeof(Person));

    if ((fp = fopen(pedfile, "r")) == NULL)
    {
        printf("ERROR: cannot open file \'%s\'\n\n", pedfile);
        exit(1);
    }

    if (fscanf(fp, "%19s", indi[i]->pedId) != 1) {
        printf("ERROR: pedigree file \'%s\' is empty or unreadable.\n\n", pedfile);
        exit(1);
    }
    strcpy(currentFamily, indi[i]->pedId);
    while (!feof(fp))
    {
        if (strcmp(currentFamily, indi[i]->pedId) == 0)
        {
            membersOfFamily[totalPed]++;
        } else {
            totalPed++;
            membersOfFamily[totalPed] = 1;
            strcpy(currentFamily, indi[i]->pedId);
        }
        if (fscanf(fp, "%19s%19s%19s%d%*[^\n]",
                   indi[i]->id, indi[i]->father, indi[i]->mother, &indi[i]->sex) != 4) {
            break;
        }

        i++;
        if (i >= MAX_IND) {
            printf("\nERROR: maximum number %d of individuals exceeded.\n\n", MAX_IND);
            exit(1);
        }
        indi[i] = (Person *) malloc(sizeof(Person));
        if (fscanf(fp, "%19s", indi[i]->pedId) != 1) break;
    }
    fclose(fp);

    totalInd = i;
    if (totalInd > MAX_IND) {
        printf("\nERROR: maximum number %d of individuals exceeded.\n\n", MAX_IND);
        exit(1);
    }
    totalPed += 1;
    if (totalPed > MAX_PED) {
        printf("\nERROR: maximum number %d of pedigrees exceeded.\n\n", MAX_PED);
        exit(1);
    }
}

/************************************************/
/*                                              */
/*              readParam                       */
/*                                              */
/************************************************/
static void readParam(const char *paramfile)
{
    int   i = 0, j;
    int   pattern, lp, mod;
    float f;
    FILE  *fp;
    char  file1[50], file2[50];
    float tmp_array[MAX_LOCUS];

    fp = fopen(paramfile, "r");
    if (fp == NULL) {
        printf("ERROR: cannot open file \'%s\'\n\n", paramfile);
        exit(1);
    }

    /* line 1: input pedigree filename, output filename */
    fscanf(fp, "%49s%49s%*[^\n]", file1, file2);

    /* line 2: three random seeds */
    fscanf(fp, "%d%d%d%*[^\n]", &seeds[0], &seeds[1], &seeds[2]);

    /* line 3: number of extra columns between sex and the 1st marker column */
    fscanf(fp, "%d%*[^\n]", &extra_col);

    /* line 4: number of replicates */
    fscanf(fp, "%d%*[^\n]", &replicates);

    /* line 5: total number of loci, number of times pattern is repeated */
    fscanf(fp, "%d%d%*[^\n]", &totalLoci, &iteration);
    if (totalLoci > MAX_LOCUS) {
        printf("\nERROR: maximum number %d of markers exceeded.\n\n", MAX_LOCUS);
        exit(1);
    }
    if ((totalLoci % iteration) != 0) {
        printf("ERROR in parameter file: the total number of loci divided by");
        printf(" the number of \ntimes the pattern should be repeated");
        printf(" must be a whole number.\n\n");
        exit(1);
    }

    /* line 6: interval distance type (1=theta, 2=Kosambi, 3=Haldane) */
    fscanf(fp, "%d%*[^\n]", &inter_dist_type);

    /* line 7: pattern of recombination fractions or map distances */
    fscanf(fp, "%d", &pattern);
    for (i = 0; i < pattern; i++)
    {
        fscanf(fp, "%f", &f);
        if (inter_dist_type == 1) {
            tmp_array[i] = f;
        } else if (inter_dist_type == 2) {        /* Kosambi -> theta */
            tmp_array[i] = 0.5f * (expf(4 * f) - 1) / (expf(4 * f) + 1);
        } else if (inter_dist_type == 3) {        /* Haldane -> theta */
            tmp_array[i] = 0.5f * (1 - expf(-2 * f));
        }
    }
    fscanf(fp, "%*[^\n]");

    /* expand pattern across all (totalLoci - 1) inter-marker intervals */
    mod = (totalLoci - 1) % pattern;
    lp  = (totalLoci - 1) / pattern;
    for (i = 0; i < lp; i++)
        for (j = 0; j < pattern; j++)
            recomb_fractions[i * pattern + j] = tmp_array[j];
    if (mod != 0)
        for (j = 0; j < mod; j++)
            recomb_fractions[i * pattern + j] = tmp_array[j];

    offset = ftell(fp);
    fclose(fp);
}

/************************************************/
/*                                              */
/*              readMarkerInfo                  */
/*                                              */
/************************************************/
static void readMarkerInfo(const char *paramfile)
{
    int   t1, t2;
    int   repeat_h, repeat_g;
    FILE  *fp;
    int   ii, jj = 0, type;
    int   no_of_haplo, no_of_loci;
    float sum, freq;
    float tmp_matrix[MAX_LOCUS][MAX_ALLELE];
    int   ary_no_of_allele[MAX_LOCUS];

    int   tmp_loci_of_set = 0;

    fp = fopen(paramfile, "r");
    if (fp == NULL) {
        printf("ERROR: cannot open file \'%s\'\n\n", paramfile);
        exit(1);
    }

    /* skip the first seven header lines via the offset recorded by readParam */
    if (fseek(fp, offset, SEEK_SET) != 0) {
        printf("ERROR: cannot move file pointer in parameter file.\n\n");
        exit(1);
    }

    /* line 8 onward: alternating blocks of haplotype / genotype info */
    if (fscanf(fp, "%d", &type) != 1) { fclose(fp); return; }

    while (!feof(fp))
    {
        if (type == 1) /* haplotype block */
        {
            fscanf(fp, "%d%d%d%*[^\n]", &no_of_haplo, &no_of_loci, &repeat_h);
            tmp_loci_of_set += no_of_loci * repeat_h;

            if (no_of_haplo > MAX_HAPLOTYPES) {
                printf("\nERROR: maximum number %d of haplotypes exceeded.\n\n", MAX_HAPLOTYPES);
                exit(1);
            }

            for (ii = 0; ii < no_of_haplo; ii++)
                fscanf(fp, "%f", &possibilities[ii]);
            fscanf(fp, "%*[^\n]");

            for (ii = 0; ii < no_of_haplo; ii++)
            {
                if (ii > 0)
                    possibilities[ii] = possibilities[ii - 1] + possibilities[ii];
                for (jj = 0; jj < no_of_loci; jj++)
                    fscanf(fp, "%d", &haplotype[ii][jj]);
                fscanf(fp, "%*[^\n]");
            }

            for (ii = 0; ii < repeat_h; ii++)
                buildFounderGeno(type, no_of_loci, no_of_haplo);
        }

        if (type == 2)  /* genotype block (per-marker allele frequencies) */
        {
            fscanf(fp, "%d%d%*[^\n]", &no_of_loci, &repeat_g);
            tmp_loci_of_set += no_of_loci * repeat_g;

            for (ii = 0; ii < no_of_loci; ii++)
            {
                fscanf(fp, "%d", &t1);
                ary_no_of_allele[ii] = t1;
                if (ary_no_of_allele[ii] > MAX_ALLELE) {
                    printf("\nERROR: maximum number %d of alleles exceeded.\n\n", MAX_ALLELE);
                    exit(1);
                }
                sum = 0;
                for (jj = 0; jj < t1 - 1; jj++)
                {
                    fscanf(fp, "%f", &freq);
                    sum += freq;
                    tmp_matrix[ii][jj] = sum;
                }
                tmp_matrix[ii][jj] = 1;
                fscanf(fp, "%*[^\n]");
            }
            for (ii = 0; ii < repeat_g; ii++)
                for (jj = 0; jj < no_of_loci; jj++)
                {
                    for (t2 = 0; t2 < ary_no_of_allele[jj]; t2++)
                        allele_freq[t2] = tmp_matrix[jj][t2];
                    buildFounderGeno(type, 1, ary_no_of_allele[jj]);
                }
        }

        if (fscanf(fp, "%d", &type) != 1) break;
    }

    if (tmp_loci_of_set != totalLoci / iteration) {
        printf("ERROR in parameter file: the total number of loci on line 5 does not match\n");
        printf("the sum of the number of loci in proceeding lines.\n\n");
        exit(1);
    }
    fclose(fp);
}

/************************************************/
/*                                              */
/*              simulate                        */
/*                                              */
/************************************************/
static void simulate(void)
{
    int ii, jj, k;
    int tmp, counter;
    int currentOne, currentPa = 0, currentMa = 0;

    tmp = 0;
    for (jj = 0; jj < totalPed; jj++)
    {
        while (1)
        {
            counter = 0;
            for (ii = 0; ii < membersOfFamily[jj]; ii++)
            {
                currentOne = tmp + ii;
                if (indi[currentOne]->filled == 1) {
                    counter++;
                } else {
                    for (k = 0; k < membersOfFamily[jj]; k++)
                    {
                        if (strcmp(indi[tmp + k]->id, indi[currentOne]->father) == 0)
                            currentPa = tmp + k;
                        if (strcmp(indi[tmp + k]->id, indi[currentOne]->mother) == 0)
                            currentMa = tmp + k;
                    }
                    if (indi[currentPa]->filled == 1 && indi[currentMa]->filled == 1)
                    {
                        assignGenotype(currentOne, currentPa, 'f');
                        assignGenotype(currentOne, currentMa, 'm');
                        indi[currentOne]->filled = 1;
                    }
                }
            }
            if (counter == membersOfFamily[jj]) break;
        }
        tmp += membersOfFamily[jj];
    }
}

/************************************************/
/*                                              */
/*              writeResult                     */
/*                                              */
/************************************************/
static void writeResult(const char *pedfile, FILE *fp2)
{
    int   i, j, lineNo;
    int   sex;
    char  sped[20], currentFamily[20];
    char  sperson[20], fa[20], ma[20];
    char  sequence[MAX_LOCUS];
    char  str[64][32];                 /* extra columns between sex and markers */
    int   c;
    FILE  *fp1;

    if (extra_col > 64) {
        printf("ERROR: extra_col=%d exceeds compiled-in limit (64).\n\n", extra_col);
        exit(1);
    }

    fp1 = fopen(pedfile, "r");
    if (fp1 == NULL) {
        printf("ERROR: cannot open file \'%s\'\n\n", pedfile);
        exit(1);
    }

    if (fscanf(fp1, "%19s", sped) != 1) { fclose(fp1); return; }
    strcpy(currentFamily, sped);
    lineNo = 0;
    while (!feof(fp1))
    {
        if (strcmp(currentFamily, sped) != 0)
        {
            familyNo++;
            strcpy(currentFamily, sped);
        }
        fprintf(fp2, "%d  ", familyNo + 1);
        if (fscanf(fp1, "%19s%19s%19s%d", sperson, fa, ma, &sex) != 4) break;
        fprintf(fp2, "%s  %s  %s  %d", sperson, fa, ma, sex);
        for (i = 0; i < extra_col; i++)
        {
            fscanf(fp1, "%31s", str[i]);
            fprintf(fp2, "  %s", str[i]);
        }
        fprintf(fp2, "  ");

        /* read the 0/1 mask of which markers to emit for this individual */
        j = 0;
        while ((c = fgetc(fp1)) != '\n' && c != EOF)
        {
            if (c != ' ' && c != '\t' && j < MAX_LOCUS)
            {
                sequence[j] = (char)c;
                j++;
            }
        }
        if (j == 1) {
            if (sequence[0] == '1')
                for (i = 0; i < totalLoci; i++)
                    fprintf(fp2, " %d %d", indi[lineNo]->geno->geno_chars[i][0], indi[lineNo]->geno->geno_chars[i][1]);
            else
                for (i = 0; i < totalLoci; i++)
                    fprintf(fp2, " 0 0");
        } else {
            for (i = 0; i < totalLoci; i++)
                if (sequence[i] == '1')
                    fprintf(fp2, " %d %d", indi[lineNo]->geno->geno_chars[i][0], indi[lineNo]->geno->geno_chars[i][1]);
                else
                    fprintf(fp2, " 0 0");
        }

        fprintf(fp2, "\n");

        if (fscanf(fp1, "%19s", sped) != 1) break;
        lineNo++;
    }
    fclose(fp1);
}

/************************************************/
/*                                              */
/*              main                            */
/*                                              */
/************************************************/
int main(int argc, char *argv[])
{
    char parameter_file[256], pedigree_file[256];
    char output[256];
    FILE *fp;
    int  i, j;

    printf("\n  ********************************************************\n");
    printf("  *                                                      *\n");
    printf("  *               Program   SimPed                       *\n");
    printf("  *                                                      *\n");
    printf("  ********************************************************\n");
    printf("\nUsage: simped parameter_file\n");

    if (argc > 2) {
        printf("ERROR: Too many command line arguments!\n\n");
        return 1;
    }
    if (argc == 2) {
        strncpy(parameter_file, argv[1], sizeof(parameter_file) - 1);
        parameter_file[sizeof(parameter_file) - 1] = '\0';
    } else {
        printf("\nParameter file -> ");
        if (scanf("%255s", parameter_file) != 1) {
            printf("ERROR: cannot read parameter filename\n\n");
            return 1;
        }
    }

    if ((fp = fopen(parameter_file, "r")) == NULL) {
        printf("ERROR: cannot open file \'%s\'\n\n", parameter_file);
        return 1;
    }
    fscanf(fp, "%255s%255s%*[^\n]", pedigree_file, output);
    fclose(fp);

    printf("\nConstants in effect:\n");
    printf("  Maximum number of pedigrees\t\t\t%d\n", MAX_PED);
    printf("  Maximum number of individuals\t\t\t%d\n", MAX_IND);
    printf("  Maximum number of loci\t\t\t%d\n", MAX_LOCUS);
    printf("  Maximum number of haplotypes\t\t\t%d\n", MAX_HAPLOTYPES);
    printf("  Maximum number of alleles at a locus\t\t%d\n", MAX_ALLELE);
    printf("\nInput pedigree file: %s\n", pedigree_file);
    printf("Output file: %s\n\n", output);

    if ((fp = fopen(output, "w")) == NULL) {
        printf("ERROR: cannot open file \'%s\'\n\n", output);
        return 1;
    }

    readPed(pedigree_file);
    readParam(parameter_file);

    for (i = 0; i < totalInd; i++)
    {
        indi[i]->geno = (struct genotype *) malloc(sizeof(struct genotype));
        if (indi[i]->geno == NULL) {
            printf("ERROR: Cannot allocate memory\n\n");
            return 1;
        }
    }

    familyNo = 0;
    for (i = 0; i < replicates; i++)
    {
        founder_marker_idx = 0;

        for (j = 1; j <= iteration; j++)
            readMarkerInfo(parameter_file);

        simulate();

        writeResult(pedigree_file, fp);

        familyNo++;
        for (j = 0; j < totalInd; j++)
            indi[j]->filled = 0;
    }
    fclose(fp);
    return 0;
}
