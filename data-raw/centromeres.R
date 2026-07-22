# ============================================================
# data-raw/centromeres.R
#
# Reproducible build script for the `centromeres` package dataset:
# real centromere positions, arm boundaries, and chromosome lengths
# for hg19 (GRCh37) and hg38 (GRCh38), chromosomes 1-22, X, Y.
#
# Run manually with `Rscript data-raw/centromeres.R` from the package
# root. Not run automatically at R CMD build/check time (network
# access is forbidden during checks) -- its only output artifact is
# the checked-in data/centromeres.rda.
#
# Sources (downloaded 2026-07-22):
#   hg19 centromeres : https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/gap.txt.gz
#                       (rows where type == "centromere"; one row per chromosome)
#   hg38 centromeres : https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz
#                       (hg38 moved centromere annotation out of gap.txt.gz and
#                       into its own table, with multiple contig fragments per
#                       chromosome -- aggregated below into one outer span each)
#   hg19 chrom sizes : https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
#   hg38 chrom sizes : https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
#
# License/provenance: these are basic genome-assembly structural
# tracks (gap locations, chromosome sizes) produced by the public
# GRC/UCSC Genome Browser project, not copyrighted gene annotation
# data; UCSC distributes them for unrestricted reuse. See
# https://genome.ucsc.edu/conditions.html.
# ============================================================

library(dplyr)

MAIN_CHRS <- c(as.character(1:22), "X", "Y")

.strip_chr_prefix <- function(x) sub("^chr", "", x)

read_gz_tsv <- function(url, col_names) {
  tmp <- tempfile(fileext = ".txt.gz")
  utils::download.file(url, tmp, quiet = TRUE, mode = "wb")
  on.exit(unlink(tmp))
  utils::read.delim(gzfile(tmp), header = FALSE, col.names = col_names,
                     colClasses = "character")
}

read_chrom_sizes <- function(url) {
  tmp <- tempfile(fileext = ".chrom.sizes")
  utils::download.file(url, tmp, quiet = TRUE, mode = "wb")
  on.exit(unlink(tmp))
  sizes <- utils::read.delim(tmp, header = FALSE,
                              col.names = c("chrom", "chr_length"),
                              colClasses = c("character", "integer"))
  sizes$chr <- .strip_chr_prefix(sizes$chrom)
  sizes[sizes$chr %in% MAIN_CHRS, c("chr", "chr_length")]
}

# ---- hg19: gap.txt.gz, type == "centromere" (already one row/chromosome) ----
gap_hg19 <- read_gz_tsv(
  "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/gap.txt.gz",
  col_names = c("bin", "chrom", "chromStart", "chromEnd", "ix", "n", "size", "type", "bridge")
)
hg19_cent <- gap_hg19 %>%
  filter(type == "centromere") %>%
  transmute(
    chr = .strip_chr_prefix(chrom),
    centromere_start = as.integer(chromStart) + 1L,  # UCSC is 0-based half-open
    centromere_end   = as.integer(chromEnd)
  ) %>%
  filter(chr %in% MAIN_CHRS)

# ---- hg38: centromeres.txt.gz, multiple contig fragments per chromosome ----
cent_hg38 <- read_gz_tsv(
  "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz",
  col_names = c("bin", "chrom", "chromStart", "chromEnd", "name")
)
hg38_cent <- cent_hg38 %>%
  mutate(chr = .strip_chr_prefix(chrom)) %>%
  filter(chr %in% MAIN_CHRS) %>%
  group_by(chr) %>%
  summarise(
    centromere_start = min(as.integer(chromStart)) + 1L,
    centromere_end   = max(as.integer(chromEnd)),
    .groups = "drop"
  )

hg19_sizes <- read_chrom_sizes("https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes")
hg38_sizes <- read_chrom_sizes("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes")

build_table <- function(cent, sizes, genome_name) {
  cent %>%
    left_join(sizes, by = "chr") %>%
    mutate(
      genome = genome_name,
      centromere_mid = (centromere_start + centromere_end) / 2,
      p_arm_start = 1L,
      p_arm_end = centromere_start - 1L,
      q_arm_start = centromere_end + 1L,
      q_arm_end = chr_length
    ) %>%
    select(genome, chr, chr_length, centromere_start, centromere_end,
           centromere_mid, p_arm_start, p_arm_end, q_arm_start, q_arm_end)
}

centromeres <- bind_rows(
  build_table(hg19_cent, hg19_sizes, "hg19"),
  build_table(hg38_cent, hg38_sizes, "hg38")
) %>%
  arrange(genome, factor(chr, levels = MAIN_CHRS)) %>%
  as.data.frame()

# Sanity checks before freezing -- fail loudly, not silently.
stopifnot(
  all(MAIN_CHRS %in% centromeres$chr[centromeres$genome == "hg19"]),
  all(MAIN_CHRS %in% centromeres$chr[centromeres$genome == "hg38"]),
  nrow(centromeres) == 2L * length(MAIN_CHRS),
  all(centromeres$centromere_start < centromeres$centromere_end),
  all(centromeres$centromere_end < centromeres$chr_length),
  !anyNA(centromeres)
)

usethis::use_data(centromeres, overwrite = TRUE)
