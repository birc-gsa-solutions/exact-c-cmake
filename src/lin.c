#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "parsers.h"
#include "sam.h"

static void compute_border_array(const char *p, size_t m, int ba[m])
{
    // Border array
    ba[0] = 0;
    for (int i = 1; i < m; ++i)
    {
        int b = ba[i - 1];
        while (b > 0 && p[i] != p[b])
            b = ba[b - 1];
        ba[i] = (p[i] == p[b]) ? b + 1 : 0;
    }

    // restricted border array
    for (int i = 0; i < m - 1; i++)
    {
        if (ba[i] > 0 && p[ba[i]] == p[i + 1])
            ba[i] = ba[ba[i] - 1];
    }
}

struct match_iter
{
    const char *x;
    size_t n;
    const char *p;
    size_t m;
    int *ba;
    int b, i;
};

static struct match_iter match_iter(const char *x, size_t n,
                                    const char *p, size_t m)
{
    int *ba = malloc(m * sizeof *ba);
    assert(ba); // unlikely, but allocations can fail.
    compute_border_array(p, m, ba);
    return (struct match_iter){
        .x = x, .n = n, .p = p, .m = m, .ba = ba, .b = 0, .i = 0};
}

// Frees the border array; assumes the actual iter is stack allocated.
static void free_iter(struct match_iter *const itr)
{
    free(itr->ba);
}

// Iterators give us ugly code; macros help, although they are still a bit
// ugly...
#define X (itr->x)
#define N (itr->n)
#define P (itr->p)
#define M (itr->m)
#define B (itr->b)
#define I (itr->i)
#define BA (itr->ba)
static int next_match(struct match_iter *itr)
{
    for (; I < N; ++I)
    {
        while (B > 0 && X[I] != P[B])
            B = BA[B - 1];
        B = (X[I] == P[B]) ? B + 1 : 0;
        if (B == M)
        {
            I++; // ready for next match = i + 1;
            return I - M;
        }
    }
    return -1;
}
#undef X
#undef N
#undef P
#undef M
#undef B
#undef I
#undef BA

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
            free_iter(&itr);
        }
    }

    fclose(reads);
    free_recs(recs);
    free(genome_buf);

    return 0;
}
