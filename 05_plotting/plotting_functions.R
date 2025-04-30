filter_relevel  <- function(dataframe, column, models_list){
# Utility function to filter a df based on the presence of "column" in "models_list", and relevelling the factors of the column to match the order of "models_list"
  dataframe %>% 
  filter({{column}} %in% models_list) %>% 
  mutate({{column}} := fct_relevel({{column}}, models_list))
}

metrics_to_wide  <- function(dataframe, metrics_list){
# Utility function to pivot a dataframe from long to wide format, based on the presence of "metrics" in the list "metrics_list"
  dataframe %>% 
  filter(metric %in% metrics_list) %>% 
  pivot_wider(names_from = metric, values_from = estimate)
}

add_y_break <- function(y_breaks = c(0, 0.5)) {
  list(
    scale_y_break(y_breaks), # ggbreak
    theme(
      axis.text.y.right = element_blank(),
      axis.ticks.y.right = element_blank(),
      axis.title.y.right = element_blank(),
      axis.title.y = element_text(hjust = 0.75)
    )
  )
}

#-
#   split1 = "test_iPSYCH1",
#   split2 = "iPSYCH2"
# ) 
calculate_performance_diff <- function(metrics_df, features_to_compare, cols_to_remove, metrics_cols, category_col = "Split_eval", split1, split2) {

  # metrics_df <- metrics_df %>% 
  #   metrics_to_wide(metrics_cols)

  # Validate inputs
  if (!all(metrics_cols %in% pull(metrics_df, metric))) {
    stop("Not all specified metrics columns exist in the dataframe")
  }
  
  if (!category_col %in% colnames(metrics_df)) {
    stop("Split column not found in dataframe")
  }
  
  # Extract metric name prefixes (e.g., "roc_auc" from "avg_roc_auc")
  metric_prefixes <- unique(gsub("^(avg_|std_|lower_|upper_)", "", metrics_cols))
  
  
  # Process each metric
  metrics_df <- metrics_df %>%
    filter(Features %in% features_to_compare) %>%
    metrics_to_wide(metrics_cols) %>% 
    select(-any_of(cols_to_remove)) %>% 
    # Pivot wider for the specified metrics
    pivot_wider(
      names_from = all_of(category_col),
      values_from = all_of(metrics_cols)
    )
  
  # For each metric prefix, calculate the differences and CIs
  for (prefix in metric_prefixes) {
    avg1_col <- paste0("avg_", prefix, "_", split1)
    avg2_col <- paste0("avg_", prefix, "_", split2)
    std1_col <- paste0("std_", prefix, "_", split1)
    std2_col <- paste0("std_", prefix, "_", split2)
    
    metrics_df <- metrics_df %>%
      mutate(
        diff_avg = .data[[avg1_col]] - .data[[avg2_col]],
        std_diff = sqrt(.data[[std1_col]]^2 + .data[[std2_col]]^2)
      ) %>%
      mutate(
        lower_diff = diff_avg - 1.96 * std_diff,
        upper_diff = diff_avg + 1.96 * std_diff
      )
  }
  
  return(metrics_df)
}

# df %>% 
# select(ID, all_of(contains("PGS")), all_of(contains("FGRS")), all_of(contains("bigstatsr"))) %>%
# skim()
# count 


#-
keep_runs_inboth <- function(df, column_to_join_by) {
  # Function to keep only the runs that appear in both models
  df_toplot_combs <- df %>%
    group_by({{column_to_join_by}}) %>%
    tally()

  n_diagnoses  <- length(unique(df$Diagnosis))
  
  threshold = case_when(
    n_diagnoses == 1 ~ n_diagnoses*2,
    n_diagnoses == 2 ~ n_diagnoses*2,
    n_diagnoses == 3 ~ n_diagnoses*2,
    n_diagnoses == 4 ~ n_diagnoses*2,
    n_diagnoses == 5 ~ n_diagnoses*2,
    TRUE ~ NA_real_
  )
  df_filtered <- df %>%
    inner_join(df_toplot_combs %>% filter(n == threshold), by = {{column_to_join_by}})
  
  return(df_filtered)
}

plot_bars_basic  <- function(dataframe, X, Y, FILL){
    if (is.null(Y) || length(Y) != 1L || !Y %in% c("roc_auc", "r2")) {
    stop("'Y' argument must be one of: 'roc_auc','r2'")
    }
    Y  <- glue("avg_{Y}")
    Y_lower  <- glue("lower_{Y}")
    Y_upper  <- glue("upper_{Y}")

    print(Y, Y_lower, Y_upper)
    ggplot(dataframe, aes(x = {{X}}, y = {{Y}}, fill = {{FILL}})) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = {{Y_lower}}, ymax = {{Y_upper}}), width = 0.2, position = position_dodge(0.9))
}



library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)


global_size = 20


# show_col(hue_pal(h = c(0, 9))(6))
# 

# COLORS <- scales::hue_pal(h = c(0,360))(6)
# methods.color <- setNames(
#   c(COLORS, COLORS[c(4, 2, 5)], "black"), 
#   c("PLR", "C+T-max", "C+T-stringent", "T-Trees", "PLR3", "C+T-all",
#     "C+T-max-0.05", "C+T-max-0.2", "C+T-max-0.8", "biglasso")
# ) 


coolors  <- c("#FF9505", "#A63D40", "#F7A9A8", "#90A959", "#6494AA")
nice_blue_red_colors  <- c("steelblue", "darkred")

cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


color_data <- data.frame(
  Color = cbp1,
  Name = cbp1,
  Count = 1:length(cbp1)
) %>% 
mutate(Name = factor(Name, levels = cbp1))

# Generate the bar plot
ggplot(color_data, aes(x = factor(Count), y = 1, fill = Color)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_identity() +
  geom_text(aes(label = Name),
            position = position_stack(vjust = 0.5), size = 5, color = "black") +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, face = "bold")) +
  labs(title = "Color blindCustom Color Palette Barplot")



# show_col(cbp1)
# colors_df  <- tibble(
#   Model = c("EIR", "Logreg", "bigstatsr"),
#   Color = c(coolors[1], coolors[3], coolors[2])
# )

cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#377EB8","#CC79A7")

gln_colors  <- tribble(
  ~Features, ~Color,
  "GLN+bigstatsr", "#8CCBFF",
  "GLN+PGS+FGRS", "#56B4E9",
  "GLN+bigstatsr+PGS+FGRS", "#0072B2"
)

gln_colors_set  <- setNames(gln_colors$Color, gln_colors$Features)


features_colors  <- tribble(
  ~Features, ~Color,
  "GLN",                    "#B7DDFD",
  "bigstatsr",              "#009E73",
  "GLN+bigstatsr",          "#56B4E9",
  "PGS",                    "#F0E442",
  "FGRS",                   "#E69F00",
  "PGS+FGRS",               "#D55E00",
  "GLN+PGS+FGRS",           "#377EB8",
  "bigstatsr+PGS+FGRS",     "#A63D40",
  "GLN+bigstatsr+PGS+FGRS", "#00588A"
)

features_colors_palette  <- setNames(features_colors$Color, features_colors$Features)


colors_df  <- tibble(
    Model = c("EIR", "GLN", "Logreg", "bigstatsr"),
    Color = c(cbp1[6], "#15616d", cbp1[7], "#ff7d00")
)
colors_set  <- setNames(colors_df$Color, colors_df$Model)

models_colors  <- tibble(
    Model = c("EIR", "GLN", "Logreg", "bigstatsr", "bigstatsr+PCs", "bigstatsr+sex+age+PCs"),
    Color = c(cbp1[6], "#B7DDFD", cbp1[7], "#009E73", "#015f46", "#013a2b")
)
models_colors_palette  <- setNames(models_colors$Color, models_colors$Model)



# #-
# colors_df  <- tibble(
#     Model = c("EIR", "GLN", "Logreg", "bigstatsr"),
#     Color = c(cbp1[6], "#005F8F", cbp1[7], "#F54100")
# )
# colors_set  <- setNames(colors_df$Color, colors_df$Model)


colors_df  <- tibble(
    Model_type = c("DL", "GLN", "Linear", "bigstatsr"),
    Color = c(cbp1[6], "#15616d", cbp1[7], "#ff7d00")
)
model_type_colors_set  <- setNames(colors_df$Color, colors_df$Model_type)





# colors_df  <- tibble(
#     Model = c("EIR", "GLN", "Logreg", "bigstatsr"),
#     Color = c(cbp1[6], "#00afb9", cbp1[7], "#f07167")
# )


# colors_df  <- tibble(
#     Model = c("EIR", "GLN", "Logreg", "bigstatsr"),
#     Color = c(cbp1[6], cbp1[4], cbp1[7], cbp1[8])
# )

colors_df_split  <- tibble(
    Split_eval = c("iPSYCH1-test", "iPSYCH2"),
    Color = c(cbp1[4], cbp1[8]))
colors_set_split  <- setNames(colors_df_split$Color, colors_df_split$Split_eval)

cbp1

colors_inter_df  <- tibble(
    Interaction = c("No-Inter", "Inter"),
    Color = c(cbp1[1], cbp1[2])

)

colors_set_inter  <- setNames(colors_inter_df$Color, colors_inter_df$Interaction)




plot_clevel <- function(df_toplot) {
  df_toplot  %>%
  ggplot(aes(y = avg_roc_auc, x = Run, color = Model), color = "black") +
  geom_point(size = 3, aes(color = Model), position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = lower_roc_auc, ymax = upper_roc_auc), width = 0.05, position = position_dodge(width = 0.3)) +
  labs(y = "AVG-AUROC", x = "", title = "") +
  scale_y_continuous(breaks = seq(0.5, 0.8, by = 0.05), limits = c(0.45, 0.8)) +
  theme_bw() +
  theme(
    text = element_text(size = 30),
    panel.border = element_rect(linewidth = 1.5),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.ticks.x = element_line(linewidth = 0.5)
  ) +
  scale_color_manual(values = colors_set)
}


#-
plot_clevel_better <- function(df_toplot, x_column, y_column, lower_column, upper_column, color_column, color_palette, ytitle, plot_title = "") {

  # print({{y_column}})
  
  # extracts string of y_column variable
  y_column_string  <- deparse(substitute(y_column))
  if(y_column_string == "avg_roc_auc"){
    scales_y  <- scale_y_continuous(breaks = seq(0.5, 0.72, by = 0.05), limits = c(0.45, 0.72))
  } else if(y_column_string == "avg_r2"){
    scales_y  <- scale_y_continuous(breaks = seq(0, 0.225, by = 0.025), limits = c(0, 0.225))
  }
  
  df_toplot  %>%
  ggplot(aes(y = {{y_column}}, x = {{x_column}}, color = {{color_column}}), color = "black") +
  geom_point(size = 3, aes(color = {{color_column}}), position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = {{lower_column}}, ymax = {{upper_column}}), width = 0.05, position = position_dodge(width = 0.3)) +
  labs(y = ytitle, x = "", title = plot_title) +
  scales_y +
  # scale_y_continuous(breaks = seq(0.5, 0.8, by = 0.05), limits = c(0.45, 0.8)) +
  my_theme +
  theme(
    # add gridlines
    panel.grid.major = element_line(color = "grey92", size = 0.5),
    panel.grid.minor = element_line(color = "grey92", size = 0.25),  text = element_text(size = 30),
    # panel.border = element_rect(linewidth = 1.5),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    # axis.ticks.x = element_line(linewidth = 0.5),
    axis.text = element_text(color = "black")
  ) +
  # reduce size text of legend title
  theme(legend.title = element_text(size = 28)) +
  scale_color_manual(values = color_palette)
}



#-
plot_barplots <- function(dataframe, to_plot, x_column, fill_column, title, 
                         facet = FALSE, facet_type = "grid", facet_var = NULL,
                         metrics = NULL, 
                         x_angle = 45, x_hjust = 1, 
                         fill_colors = NULL, 
                         legend_position = "none", text_size = 11, title_hjust = 0.5,
                         x_title = NULL, y_title = NULL,
                         add_text = FALSE, size_text = 12, bar_width = bar_width) {
  
  plot_data <- dataframe %>%
    filter(Label %in% to_plot) %>%
    mutate(Label = fct_relevel(Label, to_plot)) %>%
    filter(metric %in% metrics) %>%
    pivot_wider(names_from = metric, values_from = estimate)

  # plot_data  <- dataframe %>% 
  # filter(Features %in% to_plot) %>% 
  # mutate(Features = fct_relevel(Features, to_plot)) %>% 
  # filter(metric %in% metrics) %>% 
  # pivot_wider(names_from = metric, values_from = estimate)




  if(str_detect(metrics, "roc_auc")[1]) {
    metric_name  <- "roc_auc"
    y_max = 0.65
    y_breaks = seq(0, y_max, by = 0.05)
  } else if(str_detect(metrics, "r2")[1]) {
    metric_name  <- "r2"
    y_max = 0.055
    y_breaks  <- seq(0, y_max, by = 0.025)
  }

  # metric_name  <- "r2"
  # y_breaks  <- seq(0, 0.2, by = 0.025)

  
  p <- ggplot(plot_data, aes(x = {{x_column}}, y = .data[[glue("avg_{metric_name}")]], fill = {{fill_column}})) +
    geom_bar(stat = "identity", position = position_dodge(width = bar_width, preserve = "single"), color = "black", width = bar_width) +
    geom_errorbar(aes(ymin = .data[[glue("lower_{metric_name}")]], ymax = .data[[glue("upper_{metric_name}")]]), 
                  width = 0.2, position = position_dodge(width = bar_width)) +
    scale_y_continuous(breaks = y_breaks, limits = c(0, y_max)) +
    labs(y = y_title)+
    theme(
    legend.key.height = unit(1, "cm"),
    legend.key.width = unit(1, "cm")
  )
    # my_theme +
    # theme(axis.text.x = element_text(angle = x_angle, hjust = x_hjust)) +
    # theme(plot.title = element_text(hjust = title_hjust)) +
    # theme(legend.position = legend_position) +
    # theme(text=element_text(size=size_text))

  if (add_text) {
    p <- p + geom_text(
      aes(label = Label),
      position = position_stack(vjust = 0.5),
      angle = 90,
      color = "black",
      size = 6
    )
  }
  
  if (!is.null(fill_colors)) {
    p <- p + scale_fill_manual(values = fill_colors)
  }
  
  if (facet) {
    if (facet_type == "grid") {
      p <- p + facet_grid(reformulate(facet_var))
    } else if (facet_type == "wrap") {
      p <- p + facet_wrap(reformulate(facet_var))
    }
  }
  
  if (!is.null(x_title)) {
    p <- p + xlab(x_title)
  }
  
  return(p)
}


#-
plot_barplots_ggpubr <- function(dataframe, x_column, fill_column, fill_label, title, 
                         facet = FALSE, facet_type = "grid", facet_var = NULL,
                         metrics = NULL, 
                         x_angle = 45, x_hjust = 1, 
                         fill_colors = NULL, 
                         legend_position = "none", text_size = 11, title_hjust = 0.5,
                         x_title = NULL, y_title = NULL,
                         add_text = FALSE, pvalue_text_size = 7,
                         padding = 0.01, 
                         theme = NULL, bar_width = 0.9) {

  require(ggpubr)
  require(ggbreak)
                          
  
  # Process data
  plot_data <- dataframe %>%
    filter(metric %in% metrics) %>%
    pivot_wider(names_from = metric, values_from = estimate) %>% 
    distinct()

  # Determine metric type and set y-axis parameters
  if(str_detect(metrics[1], "auc")) {
    metric_name <- "roc_auc"
    # Find max value including error bars
    y_max_data <- max(plot_data[[paste0("upper_", metric_name)]], na.rm = TRUE)
    # Add padding and round up to nearest 0.05
    y_max <- ceiling((y_max_data * (1 + padding)) * 20) / 20
    y_breaks = seq(0, y_max, by = 0.05)
  } else if(str_detect(metrics[1], "r2")) {
    metric_name <- "r2"
    # Find max value including error bars
    y_max_data <- max(plot_data[[paste0("upper_", metric_name)]], na.rm = TRUE)
    # Add padding and round up to nearest 0.025
    y_max <- ceiling((y_max_data * (1 + padding)) * 40) / 40
    y_breaks = seq(0, y_max, by = 0.025)
  }
  y_column <- glue("avg_{metric_name}")
  # Create base plot using ggpubr
  p <- ggbarplot(
        plot_data,
        x = x_column,                        
        y = y_column,
        fill = fill_column,   
        facet.by = facet_var,
        color = "black",
        width = bar_width,
        position = position_dodge(width = bar_width, preserve = "single"),
        add = ""                         
      ) + theme + 
      labs(fill = fill_label)

  # Add error bars using custom values
  p <- p + geom_errorbar(
    aes(
      ymin = !!sym(paste0("lower_", metric_name)),
      ymax = !!sym(paste0("upper_", metric_name)),
      group = interaction(!!sym(x_column), !!sym(fill_column))
    ),
    width = 0.2,
    position = position_dodge(0.9)
  )



  # Apply theme and customizations
  # p <- p + theme(
  #   legend.key.height = unit(1, "cm"),
  #   legend.key.width = unit(1, "cm"),
  #   legend.position = legend_position,
  #   text = element_text(size = text_size),
  #   axis.text.x = element_text(angle = x_angle, hjust = x_hjust),
  #   plot.title = element_text(hjust = title_hjust)
  # )

  # Add text if requested
  if (add_text) {
    p <- p + geom_text(
      aes(label = Label),
      position = position_stack(vjust = 0.5),
      angle = 90,
      color = "black",
      size = size_text
    )
  }

  # Add custom fill colors if provided
  if (!is.null(fill_colors)) {
    p <- p + scale_fill_manual(values = fill_colors)
  }

  # Add faceting if requested
  if (facet) {
    if (facet_type == "grid") {
      p <- p + facet_grid(reformulate(facet_var))
    } else if (facet_type == "wrap") {
      p <- p + facet_wrap(reformulate(facet_var))
    }
  }

  # Add axis labels
  if (!is.null(x_title)) {
    p <- p + xlab(x_title)
  }
  if (!is.null(y_title)) {
    p <- p + ylab(y_title)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }


  return(p)
}





#-
plot_barplots3 <- function(dataframe, x_column, fill_column, title, 
                         facet = FALSE, facet_type = "grid", facet_var = NULL,
                         metrics = NULL, 
                         x_angle = 45, x_hjust = 1, 
                         fill_colors = NULL, 
                         legend_position = "none", text_size = 11, title_hjust = 0.5,
                         x_title = NULL, y_title = NULL,
                         add_text = FALSE, size_text = 12,
                         padding = 0.01, expand = NULL) {
  
  plot_data <- dataframe %>%
    filter(metric %in% metrics) %>%
    pivot_wider(names_from = metric, values_from = estimate)

  if(str_detect(metrics[1], "roc_auc")) {
    metric_name <- "roc_auc"
    # Find max value including error bars
    y_max_data <- max(plot_data[[paste0("upper_", metric_name)]], na.rm = TRUE)
    # Add padding and round up to nearest 0.05
    y_max <- ceiling((y_max_data * (1 + padding)) * 20) / 20
    y_breaks = seq(0, y_max, by = 0.05)
  } else if(str_detect(metrics[1], "r2")) {
    metric_name <- "r2"
    # Find max value including error bars
    y_max_data <- max(plot_data[[paste0("upper_", metric_name)]], na.rm = TRUE)
    # Add padding and round up to nearest 0.025
    y_max <- ceiling((y_max_data * (1 + padding)) * 40) / 40
    y_breaks = seq(0, y_max, by = 0.025)
  }

  if (!is.null(expand)) {
    scale_y_cont  <- scale_y_continuous(breaks = y_breaks, limits = c(0, y_max), expand = expand)
  } else {
    scale_y_cont  <- scale_y_continuous(breaks = y_breaks, limits = c(0, y_max))
  }
  
  p <- ggplot(plot_data, 
              aes(x = !!sym(x_column), 
                  y = !!sym(paste0("avg_", metric_name)), 
                  fill = !!sym(fill_column))) +
    geom_bar(stat = "identity", position = "dodge", color = "black") +
    geom_errorbar(
      aes(ymin = !!sym(paste0("lower_", metric_name)), 
          ymax = !!sym(paste0("upper_", metric_name))), 
      width = 0.2, 
      position = position_dodge(0.9)
    ) +
    scale_y_cont +
    my_theme+
    # scale_y_continuous(breaks = y_breaks, limits = c(0, y_max), expand = c(0, 0)) +
    labs(y = y_title) +
    theme(
      legend.key.height = unit(1, "cm"),
      legend.key.width = unit(1, "cm")
    )

  if (add_text) {
    p <- p + geom_text(
      aes(label = Label),
      position = position_stack(vjust = 0.5),
      angle = 90,
      color = "black",
      size = 6
    )
  }
  
  if (!is.null(fill_colors)) {
    p <- p + scale_fill_manual(values = fill_colors)
  }
  
  if (facet) {
    if (facet_type == "grid") {
      p <- p + facet_grid(reformulate(facet_var))
    } else if (facet_type == "wrap") {
      p <- p + facet_wrap(reformulate(facet_var))
    }
  }
  
  if (!is.null(x_title)) {
    p <- p + xlab(x_title)
  }
  
  return(p)
}


#######





plot_roc_auc_old <- function(dataframe, to_plot, x_column, fill_column, title) {
  dataframe %>%
    filter(Label %in% to_plot) %>%
    mutate(Label = fct_relevel(Label, to_plot)) %>%
    filter(metric %in% metrics_roc_auc) %>%
    pivot_wider(names_from = metric, values_from = estimate) %>%
    ggplot(., aes(x = {{x_column}}, y = avg_roc_auc, fill = {{fill_column}}), color = "black") +
    # ggplot(., aes(x = {{x_column}}, y = avg_roc_auc, fill = Model), color = "black") +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = avg_roc_auc - std_roc_auc, ymax = avg_roc_auc + std_roc_auc), width = 0.2, position = position_dodge(0.9)) +
    facet_grid("~Diagnosis") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_text(
      aes(label = Label),
      position = position_stack(vjust = 0.5),
      angle = 90,
      color = "black",
      size = 6
    ) +
    theme(text = element_text(size = 11)) +
    scale_y_continuous(breaks = seq(0, 0.8, by = 0.05)) +
    theme(axis.ticks.x = element_blank()) +
    ylab("ROC-AUC") +
    ggtitle(title) +
    # center title
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "none")+
    # color barplot based on column Model
    scale_fill_manual(values = c("EIR" = cbp1[2], "Logreg" = cbp1[1]))+
    theme(axis.text.x = element_blank())
}

plot_roc_auc_eir_logreg <- function(dataframe, to_plot, x_column, fill_column, title) {
  dataframe %>%
    filter(Label %in% to_plot) %>%
    mutate(Label = fct_relevel(Label, to_plot)) %>%
    filter(metric %in% metrics_roc_auc) %>%
    pivot_wider(names_from = metric, values_from = estimate) %>%
    # ggplot(., aes(x = {{x_column}}, y = avg_roc_auc, fill = {{fill_column}}), color = "black") +
    ggplot(., aes(x = {{x_column}}, y = avg_roc_auc, fill = Model), color = "black") +

    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = lower_roc_auc, ymax = upper_roc_auc), width = 0.2, position = position_dodge(0.9)) +
    facet_grid("~Diagnosis") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_text(
      aes(label = Label),
      position = position_stack(vjust = 0.5),
      angle = 90,
      color = "black",
      size = 6
    ) +
    theme(text = element_text(size = 11)) +
    scale_y_continuous(breaks = seq(0, 0.8, by = 0.05)) +
    theme(axis.ticks.x = element_blank()) +
    ylab("ROC-AUC") +
    ggtitle(title) +
    # center title
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "none")+
    # color barplot based on column Model
    scale_fill_manual(values = c("EIR" = cbp1[2], "Logreg" = cbp1[1]))
    
}
metrics_r2  <- c("avg_r2", "lower_r2", "upper_r2")
plot_r2 <- function(dataframe, to_plot, x_column, fill_column, title) {
  dataframe %>%
    filter(Label %in% to_plot) %>%
    mutate(Label = fct_relevel(Label, to_plot)) %>%
    filter(metric %in% metrics_r2) %>%
    pivot_wider(names_from = metric, values_from = estimate) %>%
    ggplot(., aes(x = {{x_column}}, y = avg_r2, fill = {{fill_column}}), color = "black") +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = lower_r2, ymax = upper_r2), width = 0.2, position = position_dodge(0.9)) +
    facet_grid("~Diagnosis") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    # geom_text(
    #   aes(label = Label),
    #   position = position_stack(vjust = 0.5),
    #   angle = 90,
    #   color = "black",
    #   size = 8
    # ) +
    theme(text = element_text(size = 11)) +
    scale_y_continuous(breaks = seq(0, 0.8, by = 0.05)) +
    theme(axis.ticks.x = element_blank()) +
    ylab("AVG R2 on Liability scale") +
    ggtitle(title) +
    theme(legend.position = "none")
}
# plot1r2  <- plot_r2(dataframe, to_plot1, x_column = Label, fill_column = Features_set, "R2 on Liability scale after ensembl with averaging the scores"); plot1r2


models_single_feature  <- c(
  "EIR: FGRS",
  "EIR: GLN",
  "EIR: PGS",
  "EIR: bigstatsr",
  "Logreg: FGRS",
  "Logreg: GLN",
  "Logreg: PGS",
  "Logreg: bigstatsr",
  "bigstatsr: bigstatsr"
)

# cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#882255")

metrics_roc_auc  <- c("avg_roc_auc", "lower_roc_auc", "upper_roc_auc", "std_roc_auc")

count_df  <- function(dataframe){
  dataframe  %>% 
  count(Diagnosis, Model, Interaction, CV, Split_trained, Split_eval) %>% 
  print(n = Inf)
}

count_df_labels  <- function(dataframe){
  dataframe  %>% 
  count(Diagnosis, Model, Interaction, CV, Split_trained, Split_eval, Label) %>% 
  print(n = Inf)
}



