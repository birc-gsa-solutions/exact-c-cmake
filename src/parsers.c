#include <assert.h>
#include <ctype.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "parsers.h"

// This only works if we have an actual file. If we get
// the data on a stream, such as stdin, we cannot work out what
// the size of the total buffer is, and then we would need
// another solution.
char *load_file(const char *fname)
{
    FILE *f = fopen(fname, "r");
    if (!f)
    {
        return 0;
    }

    // get the size of the file so we can read
    // it all...
    fseek(f, 0, SEEK_END);
    long fsize = ftell(f);
    fseek(f, 0, SEEK_SET);

    // Read the file and zero-terminate to
    // make it into a string.
    char *buf = malloc(fsize + 1);
    assert(buf); // quick and dirty error handling
    fread(buf, fsize, 1, f);
    if (ferror(f))
    {
        free(buf);
        buf = 0;
    }
    buf[fsize] = '\0';
    fclose(f);

    return buf;
}

// Reads a header from a string in fasta format.
// Move buf forward to the beginning of the following sequence
// and return the header.
char *read_fasta_header(char **buf)
{
    char *header, *header_end;

    // get to the first non-space character
    *buf += strspn(*buf, " \t\n\r");

    // if well-formed, we should point at a header
    assert(**buf == '>');
    (*buf)++; // skip '>'

    // skip leading space
    *buf += strspn(*buf, " \t");
    header = *buf; // this is where the header starts.

    // read until we get space, the name is only up to there
    *buf += strcspn(*buf, " \t\n\r");
    header_end = *buf;

    // skip to the newline that indicates the end of the header
    *buf = strchr(*buf, '\n');

    // remove trailing spaces
    header_end = (*buf)++; // get header and move buf past newline
    while (isspace(*header_end))
        header_end--;
    *(header_end + 1) = '\0'; // stop at the space to the right

    return header;
}

// Reads a sequence from a string in fasta format. Like above, buf moves
// forward to the next header or the end of file while we compact
// the sequence and return a pointer to it.
char *read_fasta_sequence(char **buf)
{
    char *seq = *buf;

    // Read the sequence; compact spaces until we see "\n>" or
    // the end of the buffer
    char *s = seq;
    for (;;)
    {
        // compact one line
        while (**buf && **buf != '\n')
        {
            *s++ = *(*buf)++; // compact
        }

        if (**buf == '\0' || *++(*buf) == '>')
        {
            // If pointing to the terminal zero we are done.
            // The same if we are looking at '>' after we skip
            // past the newline.
            break;
        }
    }
    *s = '\0'; // terminate sequence

    return seq;
}

static struct fasta_rec *
push_rec(const char *name, const char *seq, struct fasta_rec *next)
{
    struct fasta_rec *rec = malloc(sizeof *rec);
    assert(rec);
    *rec = (struct fasta_rec){
        .name = name,
        .seq = seq,
        .len = strlen(seq),
        .next = next};
    return rec;
}

void free_recs(struct fasta_rec *rec)
{
    // Don't free the strings, they are in a separate buffer
    struct fasta_rec *next;
    for (; rec; rec = next)
    {
        next = rec->next;
        free(rec);
    }
}

static char const *get_seq(struct fasta_rec *rec, char const *name)
{
    for (; rec; rec = rec->next)
    {
        if (strcmp(name, rec->name) == 0)
            return rec->seq;
    }
    return 0; // no such record
}

struct fasta_rec *parse_recs(char *buf)
{
    struct fasta_rec *recs = NULL;
    while (*buf)
    {
        char *header = read_fasta_header(&buf);
        char *seq = read_fasta_sequence(&buf);
        recs = push_rec(header, seq, recs);
    }
    return recs;
}

bool next_fastq_rec(FILE *f, struct fastq_rec *rec)
{
    static char buf[MAX_STRING_LEN];

    if (fgets(buf, MAX_STRING_LEN, f) == NULL)
        return false;
    strcpy(rec->name, strtok(buf + 1, "\n"));

    if (fgets(buf, MAX_STRING_LEN, f) == NULL)
        return false;
    strcpy(rec->seq, strtok(buf, "\n"));

    rec->len = strlen(rec->seq);

    return true;
}
