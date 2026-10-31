#!/usr/bin/env Rscript
# Independent reference values for the datasets written by prepare_data.pl,
# computed with PopGenome. Writes popgenome_values.tsv.
suppressPackageStartupMessages(library(PopGenome))
here <- normalizePath(".")
rows <- list()
add <- function(ds, stat, value) rows[[length(rows) + 1]] <<- data.frame(dataset = ds, statistic = stat, popgenome = as.numeric(value))

neutrality <- function(ds, ingroup, outgroup) {
    g <- readData(file.path(here, "data", ds), progress_bar_switch = FALSE)
    g <- set.populations(g, list(ingroup))
    g <- set.outgroup(g, outgroup)
    g <- neutrality.stats(g, detail = TRUE)
    g <- diversity.stats(g, pi = TRUE)
    nt <- get.neutrality(g, theta = TRUE)[[1]]
    add(ds, "n_sequences", length(ingroup))
    add(ds, "n_sites", g@n.sites)
    add(ds, "segregating_sites", nt[, "n.segregating.sites"])
    add(ds, "pi", g@Pi[1, 1])
    add(ds, "theta", g@theta_Watterson[1, 1])
    add(ds, "tajima_D", nt[, "Tajima.D"])
    add(ds, "fu_li_D", nt[, "Fu.Li.D"])
    add(ds, "fu_li_F", nt[, "Fu.Li.F"])
    # D* and F* : same statistics without an outgroup
    g2 <- readData(file.path(here, "data", ds), progress_bar_switch = FALSE)
    g2 <- set.populations(g2, list(ingroup))
    g2 <- neutrality.stats(g2, detail = TRUE)
    nt2 <- get.neutrality(g2)[[1]]
    add(ds, "fu_li_D_star", nt2[, "Fu.Li.D"])
    add(ds, "fu_li_F_star", nt2[, "Fu.Li.F"])
}

neutrality("popstats", c("A01", "A02", "A03", "A04", "A05"), "out")
neutrality("diploid", c("i1.1", "i1.2", "i2.1", "i2.2", "i3.1", "i3.2"), "o1")

# McDonald-Kreitman: ingroup sim sequences, outgroup mel, whole alignment is one CDS.
# For an alignment read with a GFF, PopGenome classifies sites as synonymous or
# not while reading; set.synnonsyn is only for SNP data.
for (ds in c("mk_CG2865", "mk_CG11099")) {
    mkt <- tryCatch({
        mk <- readData(file.path(here, "data", ds), gffpath = file.path(here, "data", paste0(ds, "_gff")),
                       progress_bar_switch = FALSE)
        ids <- get.individuals(mk)[[1]]
        mk <- set.populations(mk, list(ids[grepl("^sim", ids)], ids[grepl("^mel", ids)]))
        mk <- MKT(mk)
        r <- get.MKT(mk)
        if (is.list(r)) r <- r[[1]]
        r
    }, error = function(e) { message(ds, ": MKT step failed: ", conditionMessage(e)); NULL })
    if (!is.null(mkt)) {
        add(ds, "poly_N",  mkt[1, "P1_nonsyn"])
        add(ds, "poly_S",  mkt[1, "P1_syn"])
        add(ds, "fixed_N", mkt[1, "D_nonsyn"])
        add(ds, "fixed_S", mkt[1, "D_syn"])
    }
}

out <- do.call(rbind, rows)
write.table(out, file.path(here, "popgenome_values.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", file.path(here, "popgenome_values.tsv"), "(", nrow(out), "values )\n")
