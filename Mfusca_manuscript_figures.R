library(tidyverse)
library(ggnewscale)
library(patchwork)
library(Cairo)

format_genomic <- function(...) {
    # Format a vector of numeric values according
    # to the International System of Units.
    # http://en.wikipedia.org/wiki/SI_prefix
    #
    # Based on code by Ben Tupper
    # https://stat.ethz.ch/pipermail/r-help/2012-January/299804.html
    # Args:
    #   ...: Args passed to format()
    #
    # Returns:
    #   A function to format a vector of strings using
    #   SI prefix notation
    #
    function(x) {
        limits <- c(1e0,   1e3, 1e6)
        #prefix <- c("","Kb","Mb")
        # Vector with array indices according to position in intervals
        i <- findInterval(abs(x), limits)
        # Set prefix to " " for very small values < 1e-24
        i <- ifelse(i==0, which(limits == 1e0), i)
        paste(format(round(x/limits[i], 1),
                     trim=TRUE, scientific=FALSE, ...)
              #  ,prefix[i]
        )
    }
}

###### Fig 1 ##################
source("map_temp.R")

pic1 <- cowplot::ggdraw() + 
    cowplot::draw_image(image = "fusca_pics/pic1.png")

pic2 <- cowplot::ggdraw() + 
    cowplot::draw_image(image = "fusca_pics/pic2.png")

pic3 <- cowplot::ggdraw() + 
    cowplot::draw_image(image = "fusca_pics/pic3.png")

gs_spectra <- cowplot::ggdraw() + 
    cowplot::draw_image(image = "fusca_gs_crop.png") +
    theme(plot.margin = unit(c(-0.5, 1, -0.5, -0.5), "cm"))

layout <- "
AAABBBCCC
DDEEEEEEE
"

layout <- c(
    area(t = 1, l = 1, b = 1, r = 3),
    area(t = 1, l = 4, b = 1, r = 7),
    area(t = 1, l = 8, b = 1, r = 11),
    area(t = 2, l = 1, b = 2, r = 4),
    area(t = 2, l = 4, b = 2, r = 11)
)

pics <- cowplot::plot_grid(pic2, pic3, pic1, ncol = 3, labels = "AUTO")
plots <- cowplot::plot_grid(map, gs_spectra, 
                            ncol = 2, 
                            rel_widths = c(6, 4),
                            labels = c("D", "E"))
fig1 <- cowplot::plot_grid(pics, plots, ncol = 1, rel_heights = c(3, 4))

ggsave(plot = fig1, filename = "Fig_1_v2.pdf", device = cairo_pdf, width = 10, height = 6, units = "in")



#### gene, TE and SV distribution ####
Mfus_gff <- read_tsv("../../data/Mfusca/annotation/Mfusca_H1.ab2.maker.noseq.gff",
                     col_names = c(
                         "seqid",
                         "source",
                         "type",
                         "start",
                         "end",
                         "score",
                         "strand",
                         "phase",
                         "attributes"
                     ),
                     comment = "#"
)

Mfus_TEs <-
    read_tsv(
        "../../data/TE_anno/Mfusca_hap1.onlyChrs.chr05rc.fa.mod.EDTA.TEanno.gff3",
        col_names = c(
            "seqid",
            "source",
            "type",
            "start",
            "end",
            "score",
            "strand",
            "phase",
            "attributes"
        ),
        comment = "#"
    )


Mfus_SVs <- read_tsv("../../data/Mfusca/SV_analysis/M_fusca_hap1_hap2_50k.Assemblytics_structural_variants.bed") %>% 
    rename(seqid = reference,
           start = ref_start,
           end = ref_stop)

genesTEs <-
    bind_rows(
        Mfus_gff %>% filter(type == "gene"),
        Mfus_TEs %>% filter(!str_detect(type, "long_terminal_repeat|repeat_region|target_site_duplication")),
        Mfus_SVs
        ) %>% 
    mutate(anno = case_when(type == "gene" ~ "Genes",
                            type %in% (Mfus_SVs %>% group_by(type) %>% summarise() %>% pull()) ~ "SVs",
                            #str_detect(type, "_retrotransposon") ~ "retTE",
                            #str_detect(type, "_transposon") ~ "TEs",
                            TRUE ~ "TEs"))


scale_this <- function(x){
    (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}



binGenesTEs <- genesTEs %>% 
    filter(grepl("chr", seqid)) %>% 
    # filter(grepl("chr01", seqid)) %>% 
    mutate(facetLabel = paste0("Chr", str_extract(seqid, '[0-9][0-9]')),
           length = end - start) %>% 
    group_by(seqid, facetLabel) %>%
    mutate(bin = cut_width(start, 1e5)) %>% 
    group_by(anno, seqid, facetLabel, bin) %>% 
    count(name = "count", wt = length) %>% 
    group_by(anno, seqid, facetLabel) %>% 
    mutate(scaled_n = scale_this(count), 
           binEnd = as.numeric(str_replace(bin, "(^.*),([0-9].*)\\]|\\)", 
                                           replacement = "\\2")),
           binStart = lag(binEnd+1,default = 0),
           rel_length = count / (binEnd - binStart),
           rel_length = ifelse(rel_length > 1, 1, rel_length)
    ) %>% 
    ungroup()

# get max bins per chrom
bins <- genesTEs %>% 
    group_by(seqid) %>% 
    summarize(maxChr = max(end)) %>% 
    group_by(seqid) %>% 
    do(cutChr =seq(1, .$maxChr, by = 1e5)) %>% 
    unnest(cols = c(cutChr)) %>% 
    mutate(bin = cut_width(cutChr, 1e5))

#add empty bins to all categories 
allBins <- bins %>% 
    left_join(., binGenesTEs %>% filter(anno == "SVs"), by = c("seqid", "bin")) %>% 
    mutate(rel_length = replace_na(rel_length, 0)) %>% 
    fill(anno, facetLabel, .direction = "downup") %>% 
    bind_rows(bins %>% 
                  left_join(., binGenesTEs %>% filter(anno == "TEs"), by = c("seqid", "bin")) %>% 
                  mutate(rel_length = replace_na(rel_length, 0)) %>% 
                  fill(anno, facetLabel, .direction = "downup") ,
              bins %>% 
                  left_join(., binGenesTEs %>% filter(anno == "Genes"), by = c("seqid", "bin")) %>% 
                  mutate(rel_length = replace_na(rel_length, 0)) %>% 
                  fill(anno, facetLabel, .direction = "downup")
              ) %>% 
    group_by(anno, seqid, facetLabel) %>%
    mutate(binEnd = as.numeric(str_replace(bin, "(^.*),([0-9].*)\\]|\\)", 
                                           replacement = "\\2")),
           binStart = lag(binEnd+1,default = 0))
    
GSTdist_plot <- ggplot(allBins) +
    geom_rect(data = allBins %>% filter(anno == "Genes"),
              aes(
                  xmin = binStart ,
                  xmax = binEnd,
                  ymin = as.numeric(as.factor(anno)) + 2,
                  ymax = as.numeric(as.factor(anno)) + 2 + 0.9,
                  fill = rel_length
              )
    ) +
    scale_fill_viridis_c(aesthetics = "fill",
                         # guide = "legend",
                         name = "Relative\nlength genes",
                         limits = c(0, 1),
                         guide = guide_colorbar(order = 1),
                         rescaler = function(x, to = c(0, 1), from = NULL) {
                             ifelse(x<0.75, 
                                    scales::rescale(x,
                                                    to = to,
                                                    from = c(min(x, na.rm = TRUE), 0.75)),
                                    1)}
    ) +
    new_scale_fill() +
    geom_rect(data = allBins %>% filter(anno == "TEs"),
              aes(
                  xmin = binStart ,
                  xmax = binEnd,
                  ymin = as.numeric(as.factor(anno)) + 1,
                  ymax = as.numeric(as.factor(anno)) + 1 + 0.9,
                  fill = rel_length
              )) +
    scale_colour_viridis_c(aesthetics = "fill", 
                           # guide = "legend", 
                           name = "Relative\nlength TEs", 
                           option = "A",
                           limits = c(0, 1),
                           guide = guide_colorbar(order = 2)
    ) +
    new_scale_fill() +
    geom_rect(data = allBins %>% filter(anno == "SVs"),
              aes(
                  xmin = binStart ,
                  xmax = binEnd,
                  ymin = as.numeric(as.factor(anno)),
                  ymax = as.numeric(as.factor(anno)) + 0.9,
                  fill = rel_length
              )
    ) +
    scale_fill_viridis_c(aesthetics = "fill",
                         # guide = "legend",
                         name = "Relative\nlength SVs",
                         option = "C",
                         limits = c(0, 1),
                         guide = guide_colorbar(order = 3)
    ) +
    facet_wrap(facetLabel ~ ., strip.position = "top", scales = "free_y", nrow = 1) +
    coord_flip() +
    scale_y_reverse(breaks = 1:3 + 0.5, labels = c("SVs", "TEs", "Genes")) +
    scale_x_reverse(labels=format_genomic(), expand = c(0, 0)) +
    cowplot::theme_cowplot(font_size = 11) +
    cowplot::panel_border() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "bottom", 
          legend.direction = "horizontal",
          legend.title = element_text(vjust = 1.25),
          legend.box.margin = margin(l = 0.5, r = 0.5),
          legend.spacing.x = unit(15, "pt"),
          legend.text = element_text(angle = 45, hjust = 1))

# Sv size distribution
SV_plot <- Mfus_SVs %>% 
    filter(between(size, 500, 20000)) %>% 
    mutate(type = str_replace(type, "_", "\n")) %>% 
    ggplot() +
    geom_histogram(aes(x = size, fill = type), bins = 100) + 
    facet_grid(type ~ ., scales = "free_y") +
    scale_fill_viridis_d(option = "C") +
    labs(x = "Structural variant size (bp)", y = "Count") +
    cowplot::theme_cowplot(font_size = 11) +
    cowplot::panel_border() +
    theme(legend.position = "none")


SV_Ref_dist_plot <- Mfus_SVs %>% 
    mutate(bin = cut(size, breaks = c(0, 100, 500, 1000, 2500, 5000, 10000, 50000))) %>% 
    group_by(type) %>% 
    add_count(name = "total") %>% 
    mutate(facet_label = paste0(str_replace(type, "_", " "), " (n=", total, ")")) %>% 
    group_by(bin, facet_label) %>% 
    count() %>% 
    ggplot() +
    geom_bar(aes(x = bin, y = n, fill = facet_label), stat = "identity") + 
    facet_wrap(~ facet_label,
               ncol = 1,
               scales = "free_y") +
    cowplot::theme_cowplot(font_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    cowplot::panel_border() +
    labs(x = "Variant size (bp)", y = "Count") + 
    scale_x_discrete(labels = c("0-100", "100-500", "500-1000", "1000-2500", "2500-5000", "5000-10000", "10000-50000")) +
    guides(fill = guide_none()) + 
    scale_fill_viridis_d(option = "C")


# alignment dotplot ####
readDelta <- function(deltafile){
    lines = scan(deltafile, 'a', sep='\n', quiet=TRUE)
    lines = lines[-1]
    lines.l = strsplit(lines, ' ')
    lines.len = lapply(lines.l, length) %>% as.numeric
    lines.l = lines.l[lines.len != 1]
    lines.len = lines.len[lines.len != 1]
    head.pos = which(lines.len == 4)
    head.id = rep(head.pos, c(head.pos[-1], length(lines.l)+1)-head.pos)
    mat = matrix(as.numeric(unlist(lines.l[lines.len==7])), 7)
    res = as.data.frame(t(mat[1:5,]))
    colnames(res) = c('rs','re','qs','qe','error')
    res$qid = unlist(lapply(lines.l[head.id[lines.len==7]], '[', 2))
    res$rid = unlist(lapply(lines.l[head.id[lines.len==7]], '[', 1)) %>% gsub('^>', '', .)
    res$strand = ifelse(res$qe-res$qs > 0, '+', '-')
    res
}

filterMum <- function(df, minl=1000, flanks=1e4){
    coord = df %>% filter(abs(re-rs)>minl) %>% group_by(qid, rid) %>%
        summarize(qsL=min(qs)-flanks, qeL=max(qe)+flanks, rs=median(rs)) %>%
        ungroup %>% arrange(desc(rs)) %>%
        mutate(qid=factor(qid, levels=unique(qid))) %>% select(-rs)
    merge(df, coord) %>% filter(qs>qsL, qe<qeL) %>%
        mutate(qid=factor(qid, levels=levels(coord$qid))) %>% select(-qsL, -qeL)
}
delta <- readDelta("../../data/Mfusca/genome_alignments/fusca_vs_gddh13_mummer.delta.filter")
delta_chr <- delta
delta_chr <- filter(delta, str_detect(rid, "Chr"), re-rs > 1e4)  #%>% arrange(rid, rs, qs) %>%
#rename(A = rid, B = qid, AStart = rs, AEnd =  re, BStart = qs, BEnd =  qe)

diagMum <- function(df){
    ## Find best qid order
    rid.o = df %>% group_by(qid, rid) %>% summarize(base=sum(abs(qe-qs)),
                                                    rs=weighted.mean(rs, abs(qe-qs))) %>%
        ungroup %>% arrange(desc(base)) %>% group_by(qid) %>% do(head(., 1)) %>%
        ungroup %>% arrange(desc(rid), desc(rs)) %>%
        mutate(qid=factor(qid, levels=unique(qid)))
    ## Find best qid strand
    major.strand = df %>% group_by(qid) %>%
        summarize(major.strand=ifelse(sum(sign(qe-qs)*abs(qe-qs))>0, '+', '-'),
                  maxQ=max(c(qe, qs)))
    merge(df, major.strand) %>% mutate(qs=ifelse(major.strand=='-', maxQ-qs, qs),
                                       qe=ifelse(major.strand=='-', maxQ-qe, qe),
                                       qid=factor(qid, levels=levels(rid.o$qid)))
}

delta_chr.diag <- diagMum(delta_chr) %>%
    mutate(rlab = paste0("Chr", str_extract(rid, "[0-9][0-9]"))) %>%
    mutate(Similarity = 1 - error / abs(qe - qs)) %>%
    arrange(desc(as.numeric(qid))) %>%
    filter(rid != "Chr00")

ctg_lengths <-
    delta_chr.diag %>% group_by(qid) %>% summarise(ctg_length = max(qe)) %>%
    ungroup() %>%
    arrange(desc(as.numeric(qid))) %>%
    mutate(ypad = lag(ctg_length, default = 0),
           ycumsumpad = cumsum(ypad),
           labelypos = ycumsumpad + ctg_length/2,
           label = paste0("Chr", str_extract(qid, "[0-9][0-9]"))
    )

chrm_lengths <-
    delta_chr.diag %>% group_by(rid) %>% summarise(length = max(re)) %>%
    ungroup() %>%
    mutate(xpad = lag(length, default = 0),
           xcumsumpad = cumsum(xpad),
           labelpos = xcumsumpad + length/2,
           label = paste0("Chr", str_extract(rid, "[0-9][0-9]")))

delta_chr.diag <- delta_chr.diag %>%
    left_join(chrm_lengths, by = "rid") %>%
    left_join(ctg_lengths, by = "qid")

hapDotPlot <- delta_chr.diag %>% 
    ggplot() +
    geom_point(
        data = filter(delta_chr.diag, Similarity <= 0.94),
        aes(
            x = xcumsumpad + rs,
            y = ycumsumpad + qs,
            color = Similarity,
            size = ctg_length
        ),
        alpha = 0.5
    ) +
    geom_point(
        data = filter(delta_chr.diag, between(Similarity, 0.94, 0.98)),
        aes(
            x = xcumsumpad + rs,
            y = ycumsumpad + qs,
            color = Similarity,
            size = ctg_length
        ),
        alpha = 0.5
    ) +
    geom_point(
        data = filter(delta_chr.diag, Similarity >= 0.98),
        aes(
            x = xcumsumpad + rs,
            y = ycumsumpad + qs,
            color = Similarity,
            size = ctg_length
        ),
        alpha = 0.5
    ) +
    geom_vline(data = chrm_lengths,
               aes(xintercept = xcumsumpad),
               linetype = 2, 
               color = "grey85") +
    geom_hline(data = ctg_lengths,
               aes(yintercept = ycumsumpad),
               linetype = 2, 
               color = "grey85"
               ) +
    geom_text(data = chrm_lengths, aes(x = labelpos, y = -30e6, label = label), 
              size = 3, 
              angle = 90
              ) +
    geom_text(data = ctg_lengths, aes(y = labelypos, x = -30e6, label = label), 
              size = 3
              ) +
    scale_size_continuous(name = "Alignment\nlength") +
    labs(x = "M. domestica (GDDH13)", y = "M. fusca hap1") +
    scale_color_viridis_c(option = "A") +
    cowplot::theme_cowplot(font_size = 11) +
    theme(legend.position = c(.80, .40), 
          legend.box.background = element_rect(fill = "white", color = "white"),
          legend.margin = margin(6, 6, 6, 6))


####### deletion exon overlap ######

dels_exon <- read_tsv("../../data/Mfusca/SV_analysis/deletions50kb_intersect_genes.bed", col_names = F)

dels_exon %>% 
    filter(X20 != ".") %>% 
    mutate(igv = paste0(X1, ":", X2, "-", X3),
           SV_cmd = paste0(
               "-o ",
               X1,
               "_",
               X2,
               "_",
               X3,
               " -c ",
               str_extract(X1, "[0-9][0-9]"),
               " -s ",
               X2 , #- floor(X5 / 2),
               " -e ",
               X3, #+ floor(X5 / 2),
               " -w",
               floor((X3 - X2) / 2)
               )
           ) %>%
    group_by(X4) %>% summarise(size = max(X5), unique(SV_cmd)) %>% 
    arrange(desc(size)) 
    

#### deletions vs random #####
    
rand <- read_tsv("../../data/Mfusca/SV_analysis/all_coverage.txt", 
                 col_names = c("filename", colnames(Mfus_SVs), "Coverage")) %>% 
    mutate(filename = str_remove(filename, ".bed_coverage.txt"),
           Region = ifelse(filename == "deletions50k", "Deletions", "Random")) 

del_obs <- rand %>% filter(Region == "Deletions", Coverage < 300) %>% pull(Coverage)
rand_obs <- rand %>% filter(Region == "Random", Coverage < 300) %>% pull(Coverage)

ks <- Matching::ks.boot(del_obs, rand_obs, nboots = 1000)
summary(ks)

delplot <- rand %>% 
    filter(Coverage < 300) %>% 
    group_by(Region) %>% 
    ggplot() +
    geom_violin(aes(x = Region, y = Coverage, fill = Region)) +
    scale_fill_manual(values = c("#0D0887", "grey85"), guide = "none") +
    cowplot::theme_cowplot(font_size = 11) +
    ggpubr::geom_bracket(
        xmin = "Deletions", xmax = "Random", y.position = 325,
        label = "Bootstrap KS,\np < 2.22e-16"
    ) +
    ylim(0, 350)

####### Fig 2 Build ##########

SV_valid <- cowplot::ggdraw() + 
    cowplot::draw_image(image = "SV_validation.png")


layout <- "
AAAAAA
BCDDDD
BCDDDD
"
p <- GSTdist_plot + SV_Ref_dist_plot + delplot + SV_valid +
    plot_layout(design = layout) +
    plot_annotation(tag_levels = 'A')

ggsave(plot = p, filename = "Fig_2.pdf", device = cairo_pdf, width = 12, height = 10, units = "in")


############# fig3 Build ###########

syn_wild <- cowplot::ggdraw() + 
    cowplot::draw_image(image = magick::image_read_pdf("fusca_vs_sie_sylv.pdf", density = 600))

syn_Mfu <- cowplot::ggdraw() + 
    cowplot::draw_image(image = magick::image_read_pdf("Fig3panelD.pdf")) +
    theme(plot.margin = unit(c(-5, 0, -0.5, 0), "cm"))

tree <- cowplot::ggdraw() + 
    cowplot::draw_image(image = magick::image_read_pdf("all_rgenes_w_wild.fasta.final_tree_edit.pdf"))

layout <- "
AABBB
AABBB
AABBB
CCCCC
CCCCC
DDDDD
DDDDD
DDDDD
"

fig3 <- hapDotPlot + syn_wild + tree + syn_Mfu + 
    plot_layout(design = layout) +
    plot_annotation(tag_levels = 'A')

ggsave(plot = fig3, filename = "Fig_3_v2.pdf", device = cairo_pdf, width = 12, height = 12, units = "in")



