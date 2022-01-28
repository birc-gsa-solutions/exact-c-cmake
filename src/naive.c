#include <stdio.h>
#include <stdlib.h>

#include "parsers.h"
#include "sam.h"

struct match_iter
{
    const char *x;
    size_t n;
    const char *p;
    size_t m;
    int i;
};

static struct match_iter match_iter(const char *x, size_t n,
                                    const char *p, size_t m)
{
    return (struct match_iter){.x = x, .n = n, .p = p, .m = m, .i = 0};
}

static int next_match(struct match_iter *itr)
{

    for (int i = itr->i; i < itr->n; i++)
    {
        for (int j = 0; j < itr->m; j++)
        {
            if (itr->x[i + j] != itr->p[j])
                break;
            if (j == itr->m - 1)
            {
                // a match
                itr->i = i + 1; // next time, start from here
                return i;
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
            struct match_iter itr = match_iter(chrom->seq, chrom->len, read.seq, read.len);
            for (int i = next_match(&itr); i >= 0; i = next_match(&itr))
            {
                print_sam(stdout, read.name, read.seq, chrom->name, i);
            }
        }
    }

    fclose(reads);
    free_recs(recs);
    free(genome_buf);

    return 0;
}
