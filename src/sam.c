#include "sam.h"

#include <string.h>

void print_sam(FILE *f, const char *rname, const char *read, const char *chrom, int pos)
{
    int m = strlen(read);
    fprintf(f, "%s\t%s\t%d\t%dM\t%s\n", rname, chrom, pos + 1, m, read);
}
