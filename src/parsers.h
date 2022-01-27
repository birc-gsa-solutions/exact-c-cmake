#ifndef PARSERS_H
#define PARSERS_H

#include <stdbool.h>
#include <stdio.h>

char *load_file(const char *fname);
char *read_fasta_header(char **buf);
char *read_fasta_sequence(char **buf);

struct fasta_rec
{
    const char *name;
    const char *seq;
    size_t len;
    struct fasta_rec *next;
};
struct fasta_rec *parse_recs(char *buf);
void free_recs(struct fasta_rec *recs);

// This is a bit ugly, but it also saves a little time to fix
// the strings so I don't have to alloc and dealloc more than
// necessary
#define MAX_STRING_LEN 1000
struct fastq_rec
{
    char name[MAX_STRING_LEN];
    char seq[MAX_STRING_LEN];
    size_t len;
};

bool next_fastq_rec(FILE *f, struct fastq_rec *rec);

#endif
