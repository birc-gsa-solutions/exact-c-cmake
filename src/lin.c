#include <stdio.h>
#include <stdlib.h>

#include "parsers.h"
#include "sam.h"

// FIXME: still naive; will fix it soon
int next_match(const char *x, size_t n,
               const char *p, size_t m,
               int *i)
{
    for (; *i < n; (*i)++)
    {
        for (int j = 0; j < m; j++)
        {
            if (x[*i + j] != p[j])
                break;
            if (j == m - 1)
            {
                // a match
                (*i)++; // next time, start from here
                return *i - 1;
            }
        }
    }
    return -1; // If we get here, we are done.
}

int main(int argc, char const *argv[])
{
    if (argc != 3)
    {
        fprintf(stderr, "%s genome reads", argv[0]);
    }
    const char *genome_fname = argv[1];
    const char *reads_fname = argv[2];

    char *genome_buf = load_file(genome_fname);
    struct fasta_rec *recs = parse_recs(genome_buf);

    FILE *reads = fopen(reads_fname, "r");
    struct fastq_rec read;
    while (next_fastq_rec(reads, &read))
    {
        for (struct fasta_rec *chrom = recs; chrom; chrom = chrom->next)
        {
            for (int i = 0, pos = next_match(chrom->seq, chrom->len, read.seq, read.len, &i);
                 pos >= 0; pos = next_match(chrom->seq, chrom->len, read.seq, read.len, &i))
            {
                print_sam(stdout, read.name, read.seq, chrom->name, pos);
            }
        }
    }

    fclose(reads);
    free_recs(recs);
    free(genome_buf);

    return 0;
}
