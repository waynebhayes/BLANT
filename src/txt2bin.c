#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>

int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s k [directed]\n", argv[0]); return 1; }
    int k = atoi(argv[1]);
    int directed = (argc > 2 && strcmp(argv[2], "directed") == 0);

    char dir[256];
    if (directed) sprintf(dir, "canon_maps/directed");
    else sprintf(dir, "canon_maps");

    char mapFile[256], listFile[256], binFile[256], permBinFile[256];
    sprintf(mapFile, "%s/canon_map%d.txt", dir, k);
    sprintf(listFile, "%s/canon_list%d.txt", dir, k);
    sprintf(binFile, "%s/canon_map%d.bin", dir, k);
    sprintf(permBinFile, "%s/perm_map%d.bin", dir, k);

    FILE *fp = fopen(listFile, "r");
    if (!fp) { perror(listFile); return 1; }
    int numCanon;
    if (fscanf(fp, "%d", &numCanon) != 1) { fprintf(stderr, "Bad canon_list\n"); return 1; }
    uint64_t *canonList = malloc(numCanon * sizeof(uint64_t));
    for (int i = 0; i < numCanon; i++) {
        if (fscanf(fp, "%" SCNu64, &canonList[i]) != 1) { fprintf(stderr, "Bad entry %d\n", i); return 1; }
        int c; while ((c = fgetc(fp)) != EOF && c != '\n');
    }
    fclose(fp);

    int bits = directed ? k * (k - 1) : k * (k - 1) / 2;
    unsigned long long Bk = 1ULL << bits;

    fp = fopen(mapFile, "r");
    if (!fp) { perror(mapFile); return 1; }

    FILE *fbin = fopen(binFile, "wb");
    if (!fbin) { perror(binFile); return 1; }
    FILE *fperm = fopen(permBinFile, "wb");
    if (!fperm) { perror(permBinFile); return 1; }

    char buf[4096];
    for (unsigned long long raw = 0; raw < Bk; raw++) {
        if (!fgets(buf, sizeof(buf), fp)) { fprintf(stderr, "EOF at %llu\n", raw); return 1; }
        char *tab1 = strchr(buf, '\t');
        if (!tab1) { fprintf(stderr, "Bad line %llu\n", raw); return 1; }
        *tab1 = '\0';
        uint64_t canonical = strtoull(buf, NULL, 10);

        int lo = 0, hi = numCanon - 1, ord = -1;
        while (lo <= hi) {
            int mid = (lo + hi) / 2;
            if (canonList[mid] == canonical) { ord = mid; break; }
            if (canonList[mid] < canonical) lo = mid + 1;
            else hi = mid - 1;
        }
        if (ord < 0) { fprintf(stderr, "ord not found for %llu at %llu\n", canonical, raw); return 1; }

        uint16_t u16 = (uint16_t)ord;
        fwrite(&u16, 2, 1, fbin);

        char *permStr = tab1 + 1;
        char *tab2 = strchr(permStr, '\t');
        if (tab2) *tab2 = '\0';
        int perm[8];
        int plen = (int)strlen(permStr);
        for (int i = 0; i < k && i < plen; i++) perm[i] = permStr[i] - '0';

        uint32_t i32 = 0;
        for (int j = 0; j < k; j++) i32 |= (perm[j] << (3 * j));
        unsigned char pbytes[3] = { i32 & 0xFF, (i32 >> 8) & 0xFF, (i32 >> 16) & 0xFF };
        fwrite(pbytes, 3, 1, fperm);

        if (raw % 1000000 == 0 && raw > 0)
            fprintf(stderr, "  %llu / %llu (%.0f%%)\n", raw, Bk, 100.0 * raw / Bk);
    }
    fclose(fp); fclose(fbin); fclose(fperm);
    fprintf(stderr, "Done: %llu entries (%lluB ord + %lluB perm)\n", Bk, Bk*2ULL, Bk*3ULL);
    free(canonList);
    return 0;
}
