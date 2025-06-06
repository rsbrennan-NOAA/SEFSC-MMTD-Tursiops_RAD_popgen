library(ggplot2)
library(tidyr)
library(ggbeeswarm)
library(ggdist)


# maintenance of species boundaries


 dat <- read.csv("depths.csv")

 pops <- read.table("analysis/population_assignments_summary.txt", header=T)
 
 df <- dat %>%
   inner_join(pops, by = c("id" = "indiv")) %>%
   filter(corrected_depth != "stranding") %>%
   mutate(corrected_depth = as.numeric(corrected_depth),
          logdepth = log10(corrected_depth*-1)*-1) %>%
   mutate(logdepth = if_else(logdepth > 0, 0, logdepth)) %>%
   mutate(sixpop = str_replace_all(sixpop, "_", "\n"))  
 
 summary_stats <- df %>%
   group_by(sixpop) %>%
   summarise(
     n = n(),
     mean_depth = mean(corrected_depth, na.rm = TRUE),
     min_depth = min(corrected_depth, na.rm = TRUE),
     max_depth = max(corrected_depth, na.rm = TRUE),
     range_depth = max_depth - min_depth,
     sd_depth = sd(corrected_depth, na.rm = TRUE),
     se_depth = sd_depth / sqrt(n),
     ci_lower = mean_depth - 1.96 * se_depth,
     ci_upper = mean_depth + 1.96 * se_depth,
     .groups = 'drop'
   ) %>%
   mutate(
     ci_95 = paste0("[", round(ci_lower, 2), ", ", round(ci_upper, 2), "]")
   ) %>%
   select(sixpop, n, mean_depth,,sd_depth, se_depth, min_depth, max_depth, range_depth, ci_95)
 summary_stats
 
 summary_stats$max_depth[2:6]
 
 df%>% filter(sixpop == "Intermediate\nAtlantic") %>%
  filter(corrected_depth > -10)
 

 fill_colors <- c(
   "Coastal\nAtlantic" = "#A6DDF0",
   "Coastal\nGulf" = "#276FBF",
   "Intermediate\nAtlantic" = "#B4ED50",
   "Intermediate\nGulf" = "#2E8B57", 
   "Offshore\nAtlantic" = "#FFDD33", 
   "Offshore\nGulf" = "#C49E45"
 )
 
 y_breaks <- c(0, -log10(10), -log10(100), -log10(500), -log10(1000), -log10(2000), -log10(3000))
 y_labels <- c("1", "10", "100", "500", "1000", "2000", "3000")

 depth_plot <- ggplot(df, aes(x = sixpop, y = logdepth, fill=sixpop)) +
   geom_hline(yintercept = y_breaks, color = "grey80", linetype = "solid", linewidth = 0.3) +
   geom_violin() +
   #geom_swarm(overflow = "compress", color="black",size=0.5) +
   theme_classic(base_size = 14) +
   scale_y_continuous(
     breaks = y_breaks,
     labels = y_labels,
     limits = c(min(df$logdepth, na.rm = TRUE), max(df$logdepth, na.rm = TRUE))
   ) +
   labs(y = "Depth (m)",
          x = NULL) +
   scale_fill_manual(values=fill_colors) +
   theme(legend.position = "none",
         axis.text.x = element_text(angle = 45, hjust = 1)
   )
 
 depth_plot


# repeat for distance from shore, for supplemental:
 
 
 
 dat <- read.csv("analysis/distance.csv")
 
 pops <- read.table("analysis/population_assignments_summary.txt", header=T)
 
 df <- dat %>%
   inner_join(pops, by = c("id" = "indiv")) %>%
   filter(distance_to_shore != "stranding") %>%
   mutate(distance_to_shore = as.numeric(distance_to_shore),
          logdistance = log10(distance_to_shore)) %>%
   mutate(sixpop = str_replace_all(sixpop, "_", "\n"))  
 
 df$distance_km <- df$distance_to_shore/1000
 
 summary_stats <- df %>%
   group_by(sixpop) %>%
   summarise(
     n = n(),
     mean_depth = mean(distance_km, na.rm = TRUE),
     min_depth = min(distance_km, na.rm = TRUE),
     max_depth = max(distance_km, na.rm = TRUE),
     range_depth = max_depth - min_depth,
     sd_depth = sd(distance_to_shore, na.rm = TRUE),
     se_depth = sd_depth / sqrt(n),
     ci_lower = mean_depth - 1.96 * se_depth,
     ci_upper = mean_depth + 1.96 * se_depth,
     .groups = 'drop'
   ) %>%
   mutate(
     ci_95 = paste0("[", round(ci_lower, 2), ", ", round(ci_upper, 2), "]")
   ) %>%
   select(sixpop, n, mean_depth,,sd_depth, se_depth, min_depth, max_depth, range_depth, ci_95)
 summary_stats
 
 #df%>% filter(sixpop == "Intermediate\nAtlantic") %>%
#   filter(corrected_depth > -10)
 
 
 fill_colors <- c(
   "Coastal\nAtlantic" = "#A6DDF0",
   "Coastal\nGulf" = "#276FBF",
   "Intermediate\nAtlantic" = "#B4ED50",
   "Intermediate\nGulf" = "#2E8B57", 
   "Offshore\nAtlantic" = "#FFDD33", 
   "Offshore\nGulf" = "#C49E45"
 )
 

 y_breaks <- c(1,10,100,200,300)
 y_labels <- c("1", "10", "100", "200", "300")
 
dist_plot <- ggplot(df, aes(x = sixpop, y = distance_km, fill=sixpop)) +
   geom_hline(yintercept = y_breaks, color = "grey80", linetype = "solid", linewidth = 0.3) +
   geom_violin() +
   #geom_swarm(overflow = "compress", color="black",size=0.5) +
   theme_classic(base_size = 14) +
   labs(y = "Distance from shore (km)",
        x = NULL) +
   scale_fill_manual(values=fill_colors) +
   theme(legend.position = "right",
         legend.title = element_blank(),
         axis.text.x = element_text(angle = 45, hjust = 1)
   )
 
dist_plot

ggsave("figures/distance_from_shore.png", h=4, w=7)
ggsave("figures/distance_from_shore.pdf", h=4, w=7)

 


# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------

 
 #----------------------------------------------------
 # bgc:
 
 setwd("~/projects/Tursiops_RAD_popgen/")
 
 
 # load environment
 load(file = "analysis/bgc.RData")
 
 
 
 
 
 
 
 
 #-------------------------------------------------------------------------------
 #-------------------------------------------------------------------------------
 #-------------------------------------------------------------------------------
 #-------------------------------------------------------------------------------
 # analyze the actual results
 #-------------------------------------------------------------------------------
 #-------------------------------------------------------------------------------
 #-------------------------------------------------------------------------------
 
 #--------------
 ## cline plots
 #--------------
 
 o<-sz_out
 
 library(ggplot2)
 library(scales) # For alpha 
 
 # makehybrid index values
 h <- seq(0, 1, 0.01)
 
 # make empty df to store curves
 plot_data <- data.frame()
 
 # non-significant 
 nonsig <- (o$gradient[,2] <= 1)
 nonsig_center <- o$cent[nonsig,1]
 nonsig_v <- o$gradient[nonsig,1]
 nonsig_u <- log(nonsig_center/(1 - nonsig_center)) * nonsig_v
 
 # create non-sig curves df
 for (i in 1:length(nonsig_v)) {
   phi <- (h^nonsig_v[i])/(h^nonsig_v[i] + (1 - h)^nonsig_v[i] * exp(nonsig_u[i]))
   temp_df <- data.frame(h = h, phi = phi, 
                         group = paste0("nonsig_", i), 
                         type = "nonsig")
   plot_data <- rbind(plot_data, temp_df)
 }
 
 #---------
 # shallow curves
 shallow <- (o$gradient[,3] < 1)
 shallow_center <- o$cent[shallow,1]
 shallow_v <- o$gradient[shallow,1]
 shallow_u <- log(shallow_center/(1 - shallow_center)) * shallow_v
 
 for (i in 1:length(shallow_v)) {
   phi <- (h^shallow_v[i])/(h^shallow_v[i] + (1 - h)^shallow_v[i] * exp(shallow_u[i]))
   temp_df <- data.frame(h = h, phi = phi, 
                         group = paste0("shallow_", i), 
                         type = "shallow")
   plot_data <- rbind(plot_data, temp_df)
 }
 
 #---------------
 #  steep curves
 sig <- (o$gradient[,2] > 1)
 sig_center <- o$cent[sig,1]
 sig_v <- o$gradient[sig,1]
 sig_u <- log(sig_center/(1 - sig_center)) * sig_v
 
 for (i in 1:length(sig_v)) {
   phi <- (h^sig_v[i])/(h^sig_v[i] + (1 - h)^sig_v[i] * exp(sig_u[i]))
   temp_df <- data.frame(h = h, phi = phi, 
                         group = paste0("sig_", i), 
                         type = "sig")
   plot_data <- rbind(plot_data, temp_df)
 }
 
 # Add diagonal line data, just for plotting
 diag_line <- data.frame(h = c(0, 1), phi = c(0, 1), group = "diagonal", 
                         type = "diagonal")
 plot_data <- rbind(plot_data, diag_line)
 
 head(plot_data)
 # Create the plot
 p <- ggplot() +
   # Plot the curves by type with appropriate colors and line widths
   geom_line(data = subset(plot_data, type == "nonsig"), 
             aes(x = h, y = phi, group = group), 
             color = "grey45", linewidth = 0.5, alpha = 0.5) +
   geom_line(data = subset(plot_data, type == "shallow"), 
             aes(x = h, y = phi, group = group), 
             color = "skyblue3", linewidth = 1, alpha = 0.5) +
   geom_line(data = subset(plot_data, type == "sig"), 
             aes(x = h, y = phi, group = group), 
             color = "firebrick3", linewidth = 1, alpha = 0.5) +
   # add fake points
   geom_point(data = data.frame(h = -10, phi = -10),
              aes(x = h, y = phi, fill = "Neutral"),
              shape = 21, size = 4, color = "black", stroke = 0.5) +
   geom_point(data = data.frame(h = -10, phi = -10),
              aes(x = h, y = phi, fill = "Shallow"),
              shape = 21, size = 4, color = "black", stroke = 0.5) +
   geom_point(data = data.frame(h = -10, phi = -10),
              aes(x = h, y = phi, fill = "Steep"),
              shape = 21, size = 4, color = "black", stroke = 0.5) +
   
   # Add the diagonal line
   geom_line(data = subset(plot_data, type == "diagonal"), 
             aes(x = h, y = phi), 
             linetype = 2, linewidth = 1, color = "black") +
   scale_color_manual(name = NULL,
                      values = c("Neutral" = "grey45", 
                                 "Shallow" = "skyblue3", 
                                 "Steep" = "firebrick3")) +
   scale_fill_manual(name = NULL,
                     values = c("Neutral" = "grey45", 
                                "Shallow" = "skyblue3", 
                                "Steep" = "firebrick3")) +
   
   # Combine both legends
   guides(
     color = guide_legend(order = 1),
     fill = guide_legend(order = 2)
   ) +
   
   labs(x = "Hybrid index", 
        y = "Ancestry probability") +
   scale_x_continuous(limits = c(0, 1)) +
   scale_y_continuous(limits = c(0, 1)) +
   theme_classic(base_size = 12) +
   theme(
     
     legend.position = "bottom",
     
   )
 
 p
 #
 
 
 ggsave(file="figures/genomic_clines.png", p, 
        h=6.5, w=7)
 
 
 #---------------------------------------------------
 #---------------------------------------------------
 #---------------------------------------------------
 #---------------------------------------------------
 # plot gradients (V) CI and result:
 
 # gompert was ordering by linkage group
 #sig<-1 + (o$gradient[order(snpLgSz[,1]),2] > 1)
 
 split_cols <- strsplit(colnames(genos), "_")
 # Extract the accession numbers (first element of each split)
 chroms <- sapply(split_cols, function(x) paste0(x[1], "_", x[2]))
 SNP <- colnames(genos)
 # ggplot
 plot_data <- data.frame(
   SNP_id <- SNP,
   CHR <- chroms,
   SNP_number = 1:nrow(o$gradient),
   log_gradient = log(o$gradient[,1]),
   lower_ci = log(o$gradient[,2]),
   upper_ci = log(o$gradient[,3])
 )
 
 colnames(plot_data) <- c("SNP", "CHR", "SNP_number", "log_gradient", "lower_ci", "upper_ci")
 
 plot_data$category <- "Neutral"
 plot_data$category[o$gradient[,2] > 1] <- "Steep"  # significant
 plot_data$category[o$gradient[,3] < 1] <- "Shallow" # shallow
 
 # make a factor for CHR 
 plot_data$CHR <- factor(plot_data$CHR, levels = unique(plot_data$CHR))
 
 # colors for chromosomes
 chr_colors <- rep(c("black", "grey1"), length.out = length(levels(plot_data$CHR)))
 names(chr_colors) <- levels(plot_data$CHR)
 
 # order by chr then position:
 
 
 # Split data by category
 neutral_data <- subset(plot_data, category == "Neutral")
 shallow_data <- subset(plot_data, category == "Shallow")
 steep_data <- subset(plot_data, category == "Steep")
 
 chr_bounds <- plot_data %>%
   # We need to join with the chr_colors information
   mutate(color = chr_colors[as.character(CHR)]) %>%
   # Only keep chromosomes with black color
   filter(color == "black") %>%
   group_by(CHR) %>%
   summarize(min_snp = min(SNP_number),
             max_snp = max(SNP_number)) %>%
   mutate(ymin = -Inf, ymax = Inf)
 
 # for chr labels
 chr_levels <- levels(plot_data$CHR)
 n_chrs <- length(chr_levels)
 
 # Create new chromosome labels (1:21 plus X)
 new_chr_labels <- c(as.character(1:(n_chrs-1)), "X")
 
 # Calculate midpoints for each chromosome for label placement
 chr_midpoints <- plot_data %>%
   group_by(CHR) %>%
   summarize(mid_pos = mean(c(min(SNP_number), max(SNP_number)))) %>%
   arrange(match(CHR, chr_levels))  # Maintain original order
 
 # Assign new labels to the midpoints dataframe
 chr_midpoints$new_label <- new_chr_labels
 
 # Create the plot with 3 geom_point layers - one for each category
 p2 <- ggplot() +
   geom_rect(data = chr_bounds,
             aes(xmin = min_snp, xmax = max_snp,
                 ymin = ymin, ymax = ymax,
                 fill = CHR),
             alpha = 1, fill="grey80") +
   # Add points for each category with alternating chromosome colors
   # Add confidence interval segments colored by category
   geom_segment(data = neutral_data,
                aes(x = SNP_number, y = lower_ci, 
                    xend = SNP_number, yend = upper_ci),
                color = "grey30", linewidth = 0.3) +
   geom_segment(data = shallow_data,
                aes(x = SNP_number, y = lower_ci, 
                    xend = SNP_number, yend = upper_ci),
                color = "skyblue3", linewidth = 0.3) +
   geom_segment(data = steep_data,
                aes(x = SNP_number, y = lower_ci, 
                    xend = SNP_number, yend = upper_ci),
                color = "firebrick3", linewidth = 0.3) +
   geom_point(data = neutral_data, 
              aes(x = SNP_number, y = log_gradient, fill = "Neutral"),
              size = 2, shape = 21, color="black") +
   geom_point(data = shallow_data, 
              aes(x = SNP_number, y = log_gradient, fill = "Shallow"),
              size = 2.5, shape = 21, color="black") +
   geom_point(data = steep_data, 
              aes(x = SNP_number, y = log_gradient, , fill = "Steep"),
              size = 2.5, shape = 21, color="black") +
   #  geom_point(data = shallow_data, 
   #             aes(x = SNP_number, y = log_gradient),
   #             color = "skyblue3",
   #             size = 3, shape = 21, stroke =2) +
   #  geom_point(data = steep_data, 
   #             aes(x = SNP_number, y = log_gradient),
   #             color = "firebrick3",
   #             size = 3, shape = 21, stroke =2) +
   # Add horizontal line at y=0
   geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.8) +
   
   # Set point colors for chromosomes
   scale_color_manual(values = chr_colors, guide = "none") +  # Hide chromosome colors from legend
   
   # Add chromosome labels
   scale_x_continuous(
     breaks = chr_midpoints$mid_pos,
     labels = chr_midpoints$new_label,
     expand = c(0.01, 0.01)
   ) +
   
   # Create a dummy aesthetic for the legend
   #annotate("point", x = -Inf, y = -Inf, color = "grey30", size = 3, shape = 19) +
   #annotate("point", x = -Inf, y = -Inf, color = "skyblue3", size = 3, shape = 19) +
   #annotate("point", x = -Inf, y = -Inf, color = "firebrick3", size = 3, shape = 19) +
   
   labs(x = "Chromosome", 
        y = "Log gradient (v)") +
   theme_classic(base_size = 12) +
   theme(
     legend.position = "none",
   ) +
   scale_fill_manual(
     name = NULL,
     values = c("Neutral" = "grey50", "Shallow" = "skyblue3", "Steep" = "firebrick3"),
     guide = guide_legend(
       override.aes = list(size = 4)
     )
   ) 
   
 p2
 
 # Add manual legend
 #p2 <- p2 + guides(color = guide_legend(
 #  override.aes = list(
 #    color = c("grey30", "skyblue3", "firebrick3"),
 #    shape = c(19, 19, 19),
 #    size = c(3, 3, 3)
 #  ),
 #  labels = c("Neutral", "Shallow", "Steep")
 #))
 
 #p2
 
 
 
 
 
 
 plt1 <- ggarrange(depth_plot, p, common.legend=FALSE, nrow=1, labels="AUTO",
                   widths=c(1.7, 1))
 plt1
 final_plot <- ggarrange(plt1, p2, heights=c(1.7, 1), common.legend=FALSE, 
                         nrow=2, labels=c("", "C"))
 final_plot
 
 # Save the final combined plot
 ggsave(filename="figures/species_maintenance.pdf",
        final_plot,
        height=5, width=7)
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
  
